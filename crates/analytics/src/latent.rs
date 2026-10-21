//! Bounded CPU latent-variable methods with explicitly scoped estimators.
//!
//! MDS is classical metric/Gower scaling, not JMP's optimized Stress1/nonmetric
//! solver. Correspondence analysis accepts a nonnegative contingency matrix.
//! Factor analysis is iterated principal-axis factoring on correlations, with
//! ONE/SMC initial communalities and optional Kaiser-normalized varimax; no ML,
//! oblique rotation, or factor significance tests are implied. CCA whitens two
//! numeric blocks. Gaussian mixtures use full covariances and EM under a minimum
//! eigenvalue constraint in standardized units. SOM uses PCA grid initialization
//! and batch Gaussian-kernel local linear smoothing; singular local regressions
//! explicitly fall back to weighted means. Both clustering models save portable
//! predictions; factor/CCA score vectors are compatible with `glm_predict`.

use std::cmp::Ordering;

const MAX_CELLS: usize = 1_000_000;
const MAX_WORK: usize = 200_000_000;
const MAX_FEATURES: usize = 64;
const MAX_SEED: f64 = 9_007_199_254_740_991.;
const MAX_PARAMETERS: usize = 16384;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 10000
        || width > 128
        || data.len() > MAX_CELLS / width
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
        || args.len() > MAX_PARAMETERS
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("latent methods require finite rectangular data, at most 10000 rows,128 columns,1000000 cells".into());
    }
    let result = match operation {
        "mds" => mds(data, args),
        "correspondence_analysis" => correspondence(data, args),
        "factor_analysis" => factor(data, args),
        "canonical_correlation" => canonical(data, args),
        "gaussian_mixture" => mixture(data, args),
        "gaussian_mixture_predict" => mixture_predict(data, args),
        "self_organizing_map" => som(data, args),
        "som_predict" => som_predict(data, args),
        _ => Err(format!("unknown latent operation: {operation}")),
    }?;
    if result
        .split([':', '[', ']', ',', '}', ' '])
        .any(|t| matches!(t, "inf" | "-inf" | "NaN"))
    {
        return Err("latent calculation exceeds finite numeric range".into());
    }
    Ok(result)
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numeric range"))
    }
}
fn integer(x: f64, lo: usize, hi: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("{name} must be an integer in [{lo},{hi}]"))
    } else {
        Ok(x as usize)
    }
}
fn unit(x: f64, name: &str) -> Result<bool, String> {
    match x {
        0. => Ok(false),
        1. => Ok(true),
        _ => Err(format!("{name} must be 0 or 1")),
    }
}
fn tolerance(x: f64) -> Result<f64, String> {
    if (1e-12..=1e-2).contains(&x) {
        Ok(x)
    } else {
        Err("tolerance must be in [1e-12,1e-2]".into())
    }
}
fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.saturating_mul(b).saturating_mul(c) > MAX_WORK {
        Err("latent computation exceeds numerical work limit".into())
    } else {
        Ok(())
    }
}
fn numbers(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(v: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}
fn norm(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(0_f64, f64::hypot)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn mean(values: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut m = 0.;
    for (i, v) in values.enumerate() {
        m = finite(
            m + finite(v - m, "mean difference")? / (i + 1) as f64,
            "mean",
        )?;
    }
    Ok(m)
}
fn identity(p: usize) -> Vec<Vec<f64>> {
    (0..p)
        .map(|i| (0..p).map(|j| f64::from(i == j)).collect())
        .collect()
}
fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j]).collect())
        .collect()
}
fn multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let bt = transpose(b);
    a.iter()
        .map(|r| {
            bt.iter()
                .map(|c| finite(dot(r, c), "matrix product"))
                .collect()
        })
        .collect()
}

struct Frame {
    x: Vec<Vec<f64>>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
}
impl Frame {
    fn new(data: &[Vec<f64>], standardize: bool) -> Result<Self, String> {
        let p = data[0].len();
        if p > MAX_FEATURES || data.len() < 2 {
            return Err("fit needs >=2 rows and <=64 features".into());
        }
        let origins = data[0].clone();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|r| r[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = if standardize {
                    finite(
                        norm(data.iter().map(|r| (r[j] - origins[j]) - offsets[j]))
                            / (data.len() as f64).sqrt(),
                        "feature scale",
                    )?
                } else {
                    1.
                };
                Ok(if s == 0. { 1. } else { s })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let x = data
            .iter()
            .map(|r| {
                (0..p)
                    .map(|j| {
                        finite(
                            ((r[j] - origins[j]) - offsets[j]) / scales[j],
                            "standardized feature",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            x,
            origins,
            offsets,
            scales,
        })
    }
    fn original(&self, row: &[f64]) -> Result<Vec<f64>, String> {
        (0..row.len())
            .map(|j| {
                finite(
                    self.origins[j] + (self.offsets[j] + self.scales[j] * row[j]),
                    "original feature",
                )
            })
            .collect()
    }
    fn transform(&self, data: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
        data.iter()
            .map(|r| {
                (0..r.len())
                    .map(|j| {
                        finite(
                            ((r[j] - self.origins[j]) - self.offsets[j]) / self.scales[j],
                            "prediction feature",
                        )
                    })
                    .collect()
            })
            .collect()
    }
    fn score_models(&self, weights: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
        transpose(weights)
            .iter()
            .map(|column| {
                let mut model = vec![0., 0., 0.];
                model.extend_from_slice(&self.origins);
                model.extend_from_slice(&self.offsets);
                for (v, s) in column.iter().zip(&self.scales) {
                    model.push(finite(v / s, "score coefficient")?);
                }
                Ok(model)
            })
            .collect()
    }
}
fn cross_covariance(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    denominator: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let p = a[0].len();
    let q = b[0].len();
    work(a.len(), p, q)?;
    (0..p)
        .map(|i| {
            (0..q)
                .map(|j| {
                    finite(
                        a.iter().zip(b).map(|(x, y)| x[i] * y[j]).sum::<f64>() / denominator as f64,
                        "covariance",
                    )
                })
                .collect()
        })
        .collect()
}

struct Eigen {
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}
#[allow(clippy::needless_range_loop)]
fn eigen(matrix: &[Vec<f64>]) -> Result<Eigen, String> {
    let p = matrix.len();
    work(p, p, 64 * p)?;
    let scale = matrix
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(0_f64, f64::max);
    if scale == 0. {
        return Ok(Eigen {
            values: vec![0.; p],
            vectors: identity(p),
        });
    }
    let mut a = matrix
        .iter()
        .map(|r| r.iter().map(|v| v / scale).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut v = identity(p);
    let mut converged = p == 1;
    for _ in 0..64 {
        let mut maximum = 0_f64;
        for i in 0..p {
            for j in i + 1..p {
                let off = a[i][j];
                maximum = maximum.max(off.abs());
                if off.abs() <= 2e-15 {
                    continue;
                }
                let tau = (a[j][j] - a[i][i]) / (2. * off);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / ((-tau) + tau.hypot(1.))
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                a[i][i] -= t * off;
                a[j][j] += t * off;
                a[i][j] = 0.;
                a[j][i] = 0.;
                for k in 0..p {
                    if k != i && k != j {
                        let x = a[k][i];
                        let y = a[k][j];
                        a[k][i] = c * x - s * y;
                        a[i][k] = a[k][i];
                        a[k][j] = s * x + c * y;
                        a[j][k] = a[k][j];
                    }
                    let x = v[k][i];
                    let y = v[k][j];
                    v[k][i] = c * x - s * y;
                    v[k][j] = s * x + c * y;
                }
            }
        }
        if maximum <= 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("symmetric eigensolver did not converge".into());
    }
    let mut order = (0..p).collect::<Vec<_>>();
    order.sort_by(|i, j| a[*j][*j].total_cmp(&a[*i][*i]).then(i.cmp(j)));
    let values = order
        .iter()
        .map(|j| finite(a[*j][*j] * scale, "eigenvalue"))
        .collect::<Result<Vec<_>, _>>()?;
    let mut vectors = (0..p)
        .map(|i| order.iter().map(|j| v[i][*j]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for j in 0..p {
        let anchor = (0..p)
            .max_by(|i, k| {
                vectors[*i][j]
                    .abs()
                    .total_cmp(&vectors[*k][j].abs())
                    .then(k.cmp(i))
            })
            .unwrap();
        if vectors[anchor][j] < 0. {
            for row in &mut vectors {
                row[j] = -row[j];
            }
        }
    }
    Ok(Eigen { values, vectors })
}
fn symmetric_power(a: &[Vec<f64>], power: f64) -> Result<Vec<Vec<f64>>, String> {
    let e = eigen(a)?;
    let p = a.len();
    let threshold = e.values[0].abs() * 1e-12;
    if e.values.iter().any(|v| *v <= threshold) {
        return Err(
            "matrix is singular or numerically unidentified; increase regularization".into(),
        );
    }
    (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    finite(
                        (0..p)
                            .map(|k| e.vectors[i][k] * e.values[k].powf(power) * e.vectors[j][k])
                            .sum(),
                        "symmetric matrix power",
                    )
                })
                .collect()
        })
        .collect()
}
struct Svd {
    values: Vec<f64>,
    u: Vec<Vec<f64>>,
    v: Vec<Vec<f64>>,
}
#[allow(clippy::needless_range_loop)] // Gram-Schmidt traverses matrix columns.
fn svd(a: &[Vec<f64>]) -> Result<Svd, String> {
    let p = a.len();
    let q = a[0].len();
    let k = p.min(q);
    work(p, q, 64 * q)?;
    let scale = a.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
    let mut b = a
        .iter()
        .map(|r| {
            r.iter()
                .map(|v| if scale == 0. { 0. } else { v / scale })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut rotations = identity(q);
    let overall = norm(b.iter().flatten().copied());
    let cutoff = overall * 1e-12;
    let mut converged = false;
    // One-sided Jacobi orthogonalizes A itself, avoiding A'A and its squared
    // condition number. Tiny numerical null columns are completed separately.
    for _ in 0..64 {
        let mut changed = false;
        for i in 0..q {
            for j in i + 1..q {
                let alpha = b.iter().map(|r| r[i] * r[i]).sum::<f64>();
                let beta = b.iter().map(|r| r[j] * r[j]).sum::<f64>();
                if alpha.sqrt() <= cutoff || beta.sqrt() <= cutoff {
                    continue;
                }
                let gamma = b.iter().map(|r| r[i] * r[j]).sum::<f64>();
                if gamma.abs() <= 2e-14 * alpha.sqrt() * beta.sqrt() {
                    continue;
                }
                let tau = (beta - alpha) / (2. * gamma);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / ((-tau) + tau.hypot(1.))
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                for row in b.iter_mut().chain(&mut rotations) {
                    let x = row[i];
                    let y = row[j];
                    row[i] = c * x - s * y;
                    row[j] = s * x + c * y;
                }
                changed = true;
            }
        }
        if !changed {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("one-sided Jacobi SVD did not converge".into());
    }
    let lengths = (0..q)
        .map(|j| norm(b.iter().map(|r| r[j])))
        .collect::<Vec<_>>();
    let mut order = (0..q).collect::<Vec<_>>();
    order.sort_by(|i, j| lengths[*j].total_cmp(&lengths[*i]).then(i.cmp(j)));
    let mut values = Vec::with_capacity(k);
    let mut u = vec![vec![0.; k]; p];
    let v = rotations
        .iter()
        .map(|r| order[..k].iter().map(|j| r[*j]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for axis in 0..k {
        let column = order[axis];
        if lengths[column] > cutoff && lengths[column] > 0. {
            values.push(finite(lengths[column] * scale, "singular value")?);
            for i in 0..p {
                u[i][axis] = b[i][column] / lengths[column];
            }
        } else {
            values.push(0.);
            let mut best = vec![0.; p];
            let mut best_length = 0.;
            for anchor in 0..p {
                let mut candidate = vec![0.; p];
                candidate[anchor] = 1.;
                for _ in 0..2 {
                    for previous in 0..axis {
                        let projection = (0..p).map(|i| candidate[i] * u[i][previous]).sum::<f64>();
                        for i in 0..p {
                            candidate[i] -= projection * u[i][previous];
                        }
                    }
                }
                let length = norm(candidate.iter().copied());
                if length > best_length {
                    best = candidate;
                    best_length = length;
                }
            }
            if best_length <= 1e-12 {
                return Err("SVD null-space completion failed".into());
            }
            for i in 0..p {
                u[i][axis] = best[i] / best_length;
            }
        }
    }
    Ok(Svd { values, u, v })
}

fn mds(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 3 || data.len() < 2 || data.len() > 128 {
        return Err(
            "mds requires 2..128 rows and [dimensions=min(2,n-1),input0features/1distances=0,standardize=0]"
                .into(),
        );
    }
    let n = data.len();
    let k = integer(
        args.first().copied().unwrap_or((n - 1).min(2) as f64),
        1,
        (n - 1).min(32),
        "MDS dimensions",
    )?;
    let distances = unit(args.get(1).copied().unwrap_or(0.), "distance input")?;
    let standardized = unit(args.get(2).copied().unwrap_or(0.), "standardize")?;
    let d = if distances {
        if standardized || data[0].len() != n {
            return Err(
                "distance input must be a square matrix without feature standardization".into(),
            );
        }
        let scale = data.iter().flatten().copied().fold(0_f64, f64::max);
        let tol = f64::EPSILON * 128. * scale;
        for (i, row) in data.iter().enumerate() {
            if row[i].abs() > tol {
                return Err("distance diagonal must be zero".into());
            }
            for (j, value) in row.iter().enumerate() {
                if *value < 0. || (*value - data[j][i]).abs() > tol {
                    return Err("distances must be nonnegative and symmetric".into());
                }
            }
        }
        data.to_vec()
    } else {
        let f = Frame::new(data, standardized)?;
        f.x.iter()
            .map(|a| {
                f.x.iter()
                    .map(|b| {
                        finite(
                            norm(a.iter().zip(b).map(|(a, b)| a - b)),
                            "pairwise distance",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<_>, _>>()?
    };
    let scale = d.iter().flatten().copied().fold(0_f64, f64::max);
    let squares = d
        .iter()
        .map(|r| {
            r.iter()
                .map(|v| if scale == 0. { 0. } else { (v / scale).powi(2) })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let means = squares
        .iter()
        .map(|r| r.iter().sum::<f64>() / n as f64)
        .collect::<Vec<_>>();
    let grand = means.iter().sum::<f64>() / n as f64;
    let gram = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| -0.5 * (squares[i][j] - means[i] - means[j] + grand))
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let e = eigen(&gram)?;
    let coordinates: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            (0..k)
                .map(|j| {
                    finite(
                        e.vectors[i][j] * e.values[j].max(0.).sqrt() * scale,
                        "MDS coordinate",
                    )
                })
                .collect()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let eigenvalues = e
        .values
        .iter()
        .map(|v| finite(v * scale * scale, "MDS eigenvalue"))
        .collect::<Result<Vec<_>, _>>()?;
    let negative = e
        .values
        .iter()
        .filter(|v| **v < 0.)
        .map(|v| -v)
        .sum::<f64>();
    let absolute = e.values.iter().map(|v| v.abs()).sum::<f64>();
    let mut numerator = 0.;
    let mut denominator = 0.;
    for i in 0..n {
        for j in 0..i {
            let fitted = norm(
                coordinates[i]
                    .iter()
                    .zip(&coordinates[j])
                    .map(|(a, b)| (a - b) / if scale == 0. { 1. } else { scale }),
            );
            let target = if scale == 0. { 0. } else { d[i][j] / scale };
            numerator += (fitted - target).powi(2);
            denominator += target * target;
        }
    }
    let stress = if denominator == 0. {
        0.
    } else {
        (numerator / denominator).sqrt()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"classical_metric_gower_scaling\",\"dimensions\":{k},\"coordinates\":{},\"eigenvalues\":{},\"normalized_raw_stress\":{stress},\"stress_denominator\":\"sum_observed_squared_distances\",\"negative_eigenvalue_mass\":{},\"negative_eigenvalue_fraction\":{},\"is_euclidean\":{}}}",matrix(&coordinates),numbers(&eigenvalues),finite(negative*scale*scale,"negative eigenvalue mass")?,if absolute==0.{0.}else{negative/absolute},negative<=1e-10*absolute.max(f64::MIN_POSITIVE)))
}

fn correspondence(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let r = data.len();
    let c = data[0].len();
    if r < 2
        || c < 2
        || r > 64
        || c > 64
        || args.len() > 2
        || data.iter().flatten().any(|v| *v < 0.)
    {
        return Err("correspondence_analysis needs nonnegative 2..64 by2..64 counts and [dimensions=min(2,r-1,c-1),scaling0principal/1symmetric=0]".into());
    }
    let max_axes = (r - 1).min(c - 1);
    let k = integer(
        args.first().copied().unwrap_or(max_axes.min(2) as f64),
        1,
        max_axes,
        "CA dimensions",
    )?;
    let symmetric = unit(args.get(1).copied().unwrap_or(0.), "symmetric scaling")?;
    let largest = data.iter().flatten().copied().fold(0_f64, f64::max);
    if largest == 0. {
        return Err("contingency total is zero".into());
    }
    let sum = data.iter().flatten().map(|v| v / largest).sum::<f64>();
    let total = finite(largest * sum, "contingency total")?;
    let p = data
        .iter()
        .map(|row| row.iter().map(|v| (v / largest) / sum).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let rm = p
        .iter()
        .map(|row| row.iter().sum::<f64>())
        .collect::<Vec<_>>();
    let cm = (0..c)
        .map(|j| p.iter().map(|row| row[j]).sum::<f64>())
        .collect::<Vec<_>>();
    if rm.iter().chain(&cm).any(|v| *v <= 0.) {
        return Err(
            "zero row/column margins are not supported; remove empty categories explicitly".into(),
        );
    }
    let s = (0..r)
        .map(|i| {
            (0..c)
                .map(|j| (p[i][j] - rm[i] * cm[j]) / rm[i].sqrt() / cm[j].sqrt())
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let decomposition = svd(&s)?;
    let singular = decomposition.values[..max_axes].to_vec();
    let inertia = singular.iter().map(|s| s * s).collect::<Vec<_>>();
    let total_inertia = inertia.iter().sum::<f64>();
    let rows = (0..r)
        .map(|i| {
            (0..k)
                .map(|j| {
                    decomposition.u[i][j]
                        * if symmetric {
                            singular[j].sqrt()
                        } else {
                            singular[j]
                        }
                        / rm[i].sqrt()
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let columns = (0..c)
        .map(|i| {
            (0..k)
                .map(|j| {
                    decomposition.v[i][j]
                        * if symmetric {
                            singular[j].sqrt()
                        } else {
                            singular[j]
                        }
                        / cm[i].sqrt()
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let portions = inertia
        .iter()
        .map(|v| {
            if total_inertia == 0. {
                0.
            } else {
                v / total_inertia
            }
        })
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"simple_correspondence_analysis_svd\",\"total_frequency\":{total},\"dimensions\":{k},\"scaling\":\"{}\",\"row_masses\":{},\"column_masses\":{},\"singular_values\":{},\"inertias\":{},\"inertia_proportions\":{},\"total_inertia\":{total_inertia},\"chi_square\":{},\"row_coordinates\":{},\"column_coordinates\":{}}}",if symmetric{"symmetric_sqrt_singular_value"}else{"principal"},numbers(&rm),numbers(&cm),numbers(&singular),numbers(&inertia),numbers(&portions),finite(total*total_inertia,"CA chi square")?,matrix(&rows),matrix(&columns)))
}

fn factor(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 6 {
        return Err("factor_analysis requires [factors,rotation0none/1varimax=1,maxiter=200,tol=1e-8,prior0ONE/1SMC=1,input0rows/1correlation=0]".into());
    }
    let p = data[0].len();
    if !(2..=64).contains(&p) {
        return Err("factor analysis needs 2..64 variables".into());
    }
    let k = integer(args[0], 1, p - 1, "factor count")?;
    let rotate = unit(args.get(1).copied().unwrap_or(1.), "varimax")?;
    let maxiter = integer(
        args.get(2).copied().unwrap_or(200.),
        1,
        2000,
        "factor iterations",
    )?;
    let tol = tolerance(args.get(3).copied().unwrap_or(1e-8))?;
    let smc = unit(args.get(4).copied().unwrap_or(1.), "SMC prior")?;
    let input = unit(args.get(5).copied().unwrap_or(0.), "correlation input")?;
    work(p, p, 64 * p * maxiter)?;
    let frame = if input {
        Frame {
            x: vec![],
            origins: vec![0.; p],
            offsets: vec![0.; p],
            scales: vec![1.; p],
        }
    } else {
        Frame::new(data, true)?
    };
    let r = if input {
        if data.len() != p {
            return Err("correlation input must be square".into());
        }
        for (i, row) in data.iter().enumerate() {
            if (row[i] - 1.).abs() > 1e-10 {
                return Err("correlation diagonal must equal one".into());
            }
            for (j, value) in row.iter().enumerate() {
                if (*value - data[j][i]).abs() > 1e-10 || value.abs() > 1. + 1e-10 {
                    return Err("correlation input must be symmetric with entries in [-1,1]".into());
                }
            }
        }
        data.to_vec()
    } else {
        cross_covariance(&frame.x, &frame.x, data.len())?
    };
    let inverse = symmetric_power(&r, -1.)?;
    let mut h = (0..p)
        .map(|j| {
            if smc {
                (1. - 1. / inverse[j][j]).clamp(0., 1.)
            } else {
                1.
            }
        })
        .collect::<Vec<_>>();
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..maxiter {
        let mut reduced = r.clone();
        for j in 0..p {
            reduced[j][j] = h[j];
        }
        let e = eigen(&reduced)?;
        if e.values[k - 1] <= 0. {
            return Err(
                "requested factors do not have positive common-variance eigenvalues".into(),
            );
        }
        let next = (0..p)
            .map(|i| {
                (0..k)
                    .map(|j| e.vectors[i][j].powi(2) * e.values[j])
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        if next.iter().any(|v| *v > 1. + 1e-8) {
            return Err("Heywood case: communality exceeds one".into());
        }
        let change = next
            .iter()
            .zip(&h)
            .map(|(a, b)| (a - b).abs())
            .fold(0_f64, f64::max);
        h = next;
        iterations = iteration + 1;
        if change <= tol {
            converged = true;
            break;
        }
    }
    let mut reduced = r.clone();
    for j in 0..p {
        reduced[j][j] = h[j];
    }
    let e = eigen(&reduced)?;
    let unrotated = (0..p)
        .map(|i| {
            (0..k)
                .map(|j| e.vectors[i][j] * e.values[j].max(0.).sqrt())
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let (loadings, rotation, rotation_converged) = if rotate && k > 1 {
        varimax(&unrotated, maxiter, tol)?
    } else {
        (unrotated.clone(), identity(k), true)
    };
    let communalities = loadings.iter().map(|row| dot(row, row)).collect::<Vec<_>>();
    if communalities.iter().any(|v| *v > 1. + 1e-8) {
        return Err("Heywood case after final factor extraction".into());
    }
    let uniqueness = communalities
        .iter()
        .map(|v| (1. - v).max(0.))
        .collect::<Vec<_>>();
    let mut reproduced = multiply(&loadings, &transpose(&loadings))?;
    for j in 0..p {
        reproduced[j][j] += uniqueness[j];
    }
    let residual = finite(
        norm(
            (0..p)
                .flat_map(|i| (0..p).map(move |j| (i, j)))
                .filter(|(i, j)| i != j)
                .map(|(i, j)| r[i][j] - reproduced[i][j]),
        ) / ((p * (p - 1)) as f64).sqrt(),
        "factor residual RMS",
    )?;
    let weights = multiply(&inverse, &loadings)?;
    let scores = if input {
        vec![]
    } else {
        multiply(&frame.x, &weights)?
    };
    let models = frame.score_models(&weights)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"iterated_principal_axis_factoring\",\"factors\":{k},\"rotation\":\"{}\",\"prior_communality\":\"{}\",\"converged\":{},\"principal_axis_converged\":{converged},\"rotation_converged\":{rotation_converged},\"iterations\":{iterations},\"loadings\":{},\"unrotated_loadings\":{},\"rotation_matrix\":{},\"communalities\":{},\"uniquenesses\":{},\"reproduced_correlation\":{},\"off_diagonal_residual_rms\":{residual},\"score_method\":\"regression_R_inverse_loadings\",\"scores\":{},\"score_predict_parameters\":{},\"score_input\":\"{}\"}}",if rotate{"kaiser_normalized_varimax"}else{"none"},if smc{"SMC"}else{"ONE"},converged&&rotation_converged,matrix(&loadings),matrix(&unrotated),matrix(&rotation),numbers(&communalities),numbers(&uniqueness),matrix(&reproduced),matrix(&scores),matrix(&models),if input{"standardized_variables"}else{"original_variables"}))
}
type RotatedLoadings = (Vec<Vec<f64>>, Vec<Vec<f64>>, bool);
fn varimax(loadings: &[Vec<f64>], maxiter: usize, tol: f64) -> Result<RotatedLoadings, String> {
    let p = loadings.len();
    let k = loadings[0].len();
    let scales = loadings
        .iter()
        .map(|r| norm(r.iter().copied()))
        .collect::<Vec<_>>();
    let mut l = loadings
        .iter()
        .zip(&scales)
        .map(|(r, s)| {
            r.iter()
                .map(|v| if *s == 0. { 0. } else { v / s })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut rotation = identity(k);
    let mut converged = false;
    for _ in 0..maxiter {
        let mut max_angle = 0_f64;
        for a in 0..k {
            for b in a + 1..k {
                let u = l
                    .iter()
                    .map(|r| r[a] * r[a] - r[b] * r[b])
                    .collect::<Vec<_>>();
                let v = l.iter().map(|r| 2. * r[a] * r[b]).collect::<Vec<_>>();
                let su = u.iter().sum::<f64>();
                let sv = v.iter().sum::<f64>();
                let aa = dot(&u, &u) - dot(&v, &v) - (su * su - sv * sv) / p as f64;
                let bb = 2. * dot(&u, &v) - 2. * su * sv / p as f64;
                let angle = 0.25 * bb.atan2(aa);
                max_angle = max_angle.max(angle.abs());
                let c = angle.cos();
                let s = angle.sin();
                for r in l.iter_mut().chain(&mut rotation) {
                    let x = r[a];
                    let y = r[b];
                    r[a] = c * x + s * y;
                    r[b] = -s * x + c * y;
                }
            }
        }
        if max_angle <= tol {
            converged = true;
            break;
        }
    }
    for (i, r) in l.iter_mut().enumerate() {
        for v in r {
            *v = finite(*v * scales[i], "rotated loading")?;
        }
    }
    Ok((l, rotation, converged))
}

fn canonical(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 3 || data.len() < 3 || data[0].len() > 64 {
        return Err("canonical_correlation requires >=3 rows, <=64 columns and [x_columns,dimensions=min(p,q),ridge=0]".into());
    }
    let p = integer(args[0], 1, data[0].len() - 1, "X columns")?;
    let q = data[0].len() - p;
    if p > 32 || q > 32 {
        return Err("CCA supports <=32 variables per block".into());
    }
    let k = integer(
        args.get(1).copied().unwrap_or(p.min(q) as f64),
        1,
        p.min(q),
        "CCA dimensions",
    )?;
    let ridge = args.get(2).copied().unwrap_or(0.);
    if ridge < 0. {
        return Err("CCA ridge must be nonnegative".into());
    }
    let xf = Frame::new(
        &data.iter().map(|r| r[..p].to_vec()).collect::<Vec<_>>(),
        true,
    )?;
    let yf = Frame::new(
        &data.iter().map(|r| r[p..].to_vec()).collect::<Vec<_>>(),
        true,
    )?;
    let mut xx = cross_covariance(&xf.x, &xf.x, data.len() - 1)?;
    let mut yy = cross_covariance(&yf.x, &yf.x, data.len() - 1)?;
    let xy = cross_covariance(&xf.x, &yf.x, data.len() - 1)?;
    for (i, r) in xx.iter_mut().enumerate() {
        r[i] = finite(r[i] + ridge, "CCA regularization")?;
    }
    for (i, r) in yy.iter_mut().enumerate() {
        r[i] = finite(r[i] + ridge, "CCA regularization")?;
    }
    let wx = symmetric_power(&xx, -0.5)?;
    let wy = symmetric_power(&yy, -0.5)?;
    let decomposition = svd(&multiply(&multiply(&wx, &xy)?, &wy)?)?;
    let correlations = decomposition.values[..k]
        .iter()
        .map(|s| {
            if *s <= 1. + 1e-9 {
                Ok(s.clamp(0., 1.))
            } else {
                Err("CCA correlation exceeds one numerically".into())
            }
        })
        .collect::<Result<Vec<_>, String>>()?;
    let a = multiply(
        &wx,
        &decomposition
            .u
            .iter()
            .map(|r| r[..k].to_vec())
            .collect::<Vec<_>>(),
    )?;
    let b = multiply(
        &wy,
        &decomposition
            .v
            .iter()
            .map(|r| r[..k].to_vec())
            .collect::<Vec<_>>(),
    )?;
    let xs = multiply(&xf.x, &a)?;
    let ys = multiply(&yf.x, &b)?;
    let actual = (0..k)
        .map(|j| {
            let x = xs.iter().map(|r| r[j]).collect::<Vec<_>>();
            let y = ys.iter().map(|r| r[j]).collect::<Vec<_>>();
            finite(
                dot(&x, &y) / norm(x.iter().copied()) / norm(y.iter().copied()),
                "canonical score correlation",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"sample_covariance_whitened_canonical_correlation\",\"x_columns\":{p},\"y_columns\":{q},\"dimensions\":{k},\"ridge\":{ridge},\"canonical_correlations\":{},\"score_correlations\":{},\"x_scores\":{},\"y_scores\":{},\"x_score_predict_parameters\":{},\"y_score_predict_parameters\":{},\"normalization\":\"unit_regularized_sample_covariance_variance\",\"inference\":\"not_supplied\"}}",numbers(&correlations),numbers(&actual),matrix(&xs),matrix(&ys),matrix(&xf.score_models(&a)?),matrix(&yf.score_models(&b)?)))
}

struct Random(u64);
impl Random {
    fn new(seed: f64) -> Result<Self, String> {
        if !(0. ..=MAX_SEED).contains(&seed) || seed.fract() != 0. {
            Err("seed must be an exact integer in [0,2^53-1]".into())
        } else {
            Ok(Self(seed as u64))
        }
    }
    fn uniform(&mut self) -> f64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        ((z ^ (z >> 31)) >> 11) as f64 / 9_007_199_254_740_992.
    }
    fn index(&mut self, n: usize) -> usize {
        (self.uniform() * n as f64) as usize
    }
}
#[derive(Clone)]
struct Gaussian {
    weights: Vec<f64>,
    means: Vec<Vec<f64>>,
    covariances: Vec<Vec<Vec<f64>>>,
}
struct Density {
    whitening: Vec<Vec<f64>>,
    logdet: f64,
}
impl Density {
    fn new(cov: &[Vec<f64>]) -> Result<Self, String> {
        let e = eigen(cov)?;
        let p = cov.len();
        if e.values
            .iter()
            .any(|v| *v <= 0. || *v <= e.values[0] * 1e-14)
        {
            return Err("mixture covariance is singular; increase covariance floor".into());
        }
        let whitening = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        finite(
                            e.vectors[j][i] / e.values[i].sqrt(),
                            "mixture whitening coefficient",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            whitening,
            logdet: finite(
                e.values.iter().map(|v| v.ln()).sum(),
                "covariance log determinant",
            )?,
        })
    }
    fn log(&self, x: &[f64], mean: &[f64]) -> Result<f64, String> {
        let d = x.iter().zip(mean).map(|(x, m)| x - m).collect::<Vec<_>>();
        let quadratic = finite(
            norm(self.whitening.iter().map(|row| dot(row, &d))).powi(2),
            "Gaussian quadratic form",
        )?;
        finite(
            -0.5 * (x.len() as f64 * (2. * std::f64::consts::PI).ln() + self.logdet + quadratic),
            "Gaussian log density",
        )
    }
}
fn mixture_evaluate(model: &Gaussian, x: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, Vec<f64>), String> {
    let densities = model
        .covariances
        .iter()
        .map(|c| Density::new(c))
        .collect::<Result<Vec<_>, _>>()?;
    let k = model.weights.len();
    let mut probabilities = Vec::with_capacity(x.len());
    let mut likelihoods = Vec::with_capacity(x.len());
    for row in x {
        let logs = (0..k)
            .map(|c| {
                densities[c]
                    .log(row, &model.means[c])
                    .map(|v| v + model.weights[c].ln())
            })
            .collect::<Result<Vec<_>, _>>()?;
        let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let terms = logs.iter().map(|v| (v - max).exp()).collect::<Vec<_>>();
        let total = terms.iter().sum::<f64>();
        likelihoods.push(finite(max + total.ln(), "mixture sample log likelihood")?);
        probabilities.push(terms.iter().map(|v| v / total).collect());
    }
    Ok((probabilities, likelihoods))
}
fn floored_covariance(cov: &[Vec<f64>], floor: f64) -> Result<Vec<Vec<f64>>, String> {
    let e = eigen(cov)?;
    let p = cov.len();
    (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    finite(
                        (0..p)
                            .map(|k| e.vectors[i][k] * e.values[k].max(floor) * e.vectors[j][k])
                            .sum(),
                        "floored covariance",
                    )
                })
                .collect()
        })
        .collect()
}
fn initialize_mixture(
    x: &[Vec<f64>],
    k: usize,
    floor: f64,
    rng: &mut Random,
) -> Result<Gaussian, String> {
    let n = x.len();
    let mut means = vec![x[rng.index(n)].clone()];
    while means.len() < k {
        let distances = x
            .iter()
            .map(|r| {
                means
                    .iter()
                    .map(|m| norm(r.iter().zip(m).map(|(a, b)| a - b)).powi(2))
                    .fold(f64::INFINITY, f64::min)
            })
            .collect::<Vec<_>>();
        let sum = distances.iter().sum::<f64>();
        let index = if sum == 0. {
            means.len() % n
        } else {
            let target = rng.uniform() * sum;
            let mut cumulative = 0.;
            distances
                .iter()
                .position(|d| {
                    cumulative += d;
                    cumulative > target
                })
                .unwrap_or(n - 1)
        };
        means.push(x[index].clone());
    }
    let covariance = floored_covariance(&cross_covariance(x, x, n)?, floor)?;
    Density::new(&covariance)?;
    Ok(Gaussian {
        weights: vec![1. / k as f64; k],
        means,
        covariances: vec![covariance; k],
    })
}
fn mixture_step(
    x: &[Vec<f64>],
    probabilities: &[Vec<f64>],
    floor: f64,
) -> Result<Gaussian, String> {
    let n = x.len();
    let p = x[0].len();
    let k = probabilities[0].len();
    let counts = (0..k)
        .map(|c| probabilities.iter().map(|r| r[c]).sum::<f64>())
        .collect::<Vec<_>>();
    if counts.iter().any(|v| *v <= 1e-12 * n as f64) {
        return Err("mixture component collapsed; reduce component count or change seed".into());
    }
    let weights = counts.iter().map(|v| v / n as f64).collect();
    let means = (0..k)
        .map(|c| {
            (0..p)
                .map(|j| {
                    finite(
                        x.iter()
                            .zip(probabilities)
                            .map(|(r, w)| w[c] * r[j] / counts[c])
                            .sum(),
                        "component mean",
                    )
                })
                .collect()
        })
        .collect::<Result<Vec<Vec<f64>>, _>>()?;
    let mut covariances = Vec::with_capacity(k);
    for c in 0..k {
        let covariance = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        finite(
                            x.iter()
                                .zip(probabilities)
                                .map(|(r, w)| {
                                    w[c] * (r[i] - means[c][i]) * (r[j] - means[c][j]) / counts[c]
                                })
                                .sum(),
                            "component covariance",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<_>, _>>()?;
        covariances.push(floored_covariance(&covariance, floor)?);
    }
    Ok(Gaussian {
        weights,
        means,
        covariances,
    })
}
fn labels(probabilities: &[Vec<f64>]) -> Vec<f64> {
    probabilities
        .iter()
        .map(|r| {
            let mut best = 0;
            for j in 1..r.len() {
                if r[j] > r[best] {
                    best = j;
                }
            }
            best as f64
        })
        .collect()
}
fn mixture(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 6 || data[0].len() > 32 {
        return Err("gaussian_mixture requires <=32 features and [components,maxiter=200,tol=1e-6,covariance_floor=1e-6,seed=0,n_init=1]".into());
    }
    let n = data.len();
    let p = data[0].len();
    let k = integer(args[0], 1, n.min(32), "mixture components")?;
    if 2 + 3 * p + k + k * p + k * p * p > MAX_PARAMETERS {
        return Err("saved Gaussian mixture would exceed the native16384-parameter limit".into());
    }
    let maxiter = integer(
        args.get(1).copied().unwrap_or(200.),
        1,
        2000,
        "EM iterations",
    )?;
    let tol = tolerance(args.get(2).copied().unwrap_or(1e-6))?;
    let floor = args.get(3).copied().unwrap_or(1e-6);
    if floor < 0. {
        return Err("covariance floor must be nonnegative".into());
    }
    let seed = args.get(4).copied().unwrap_or(0.);
    let mut rng = Random::new(seed)?;
    let starts = integer(args.get(5).copied().unwrap_or(1.), 1, 16, "mixture starts")?;
    work(n, k * p * p, maxiter * starts)?;
    let frame = Frame::new(data, true)?;
    let mut best = None;
    let mut last_error = "no mixture initialization".to_string();
    for _ in 0..starts {
        let trial = (|| -> Result<_, String> {
            let mut model = initialize_mixture(&frame.x, k, floor, &mut rng)?;
            let (mut probabilities, mut ll) = mixture_evaluate(&model, &frame.x)?;
            let mut value = ll.iter().sum::<f64>();
            let mut history = vec![value];
            let mut converged = false;
            let mut iterations = 0;
            for iteration in 0..maxiter {
                let next = mixture_step(&frame.x, &probabilities, floor)?;
                let (new_probabilities, new_ll) = mixture_evaluate(&next, &frame.x)?;
                let next_value = new_ll.iter().sum::<f64>();
                let improvement = next_value - value;
                if improvement < -1e-8 * (1. + value.abs()) {
                    return Err("EM likelihood decreased beyond numerical tolerance".into());
                }
                model = next;
                probabilities = new_probabilities;
                ll = new_ll;
                value = next_value;
                history.push(value);
                iterations = iteration + 1;
                if improvement.abs() <= tol * (1. + value.abs()) {
                    converged = true;
                    break;
                }
            }
            Ok((model, ll, history, value, converged, iterations))
        })();
        match trial {
            Ok(candidate) => {
                if best.as_ref().is_none_or(
                    |old: &(Gaussian, Vec<f64>, Vec<f64>, f64, bool, usize)| candidate.3 > old.3,
                ) {
                    best = Some(candidate);
                }
            }
            Err(error) => last_error = error,
        }
    }
    let (model, _, mut history, _, converged, iterations) = best.ok_or(last_error)?;
    let mut order = (0..k).collect::<Vec<_>>();
    order.sort_by(|a, b| {
        model.means[*a]
            .iter()
            .zip(&model.means[*b])
            .map(|(a, b)| a.total_cmp(b))
            .find(|c| *c != Ordering::Equal)
            .unwrap_or(a.cmp(b))
    });
    let model = Gaussian {
        weights: order.iter().map(|i| model.weights[*i]).collect(),
        means: order.iter().map(|i| model.means[*i].clone()).collect(),
        covariances: order
            .iter()
            .map(|i| model.covariances[*i].clone())
            .collect(),
    };
    let (probabilities, mut ll) = mixture_evaluate(&model, &frame.x)?;
    let log_scale = finite(
        frame.scales.iter().map(|s| s.ln()).sum(),
        "feature log scale",
    )?;
    for value in &mut ll {
        *value -= log_scale;
    }
    for value in &mut history {
        *value -= n as f64 * log_scale;
    }
    let likelihood = finite(ll.iter().sum(), "mixture log likelihood")?;
    let parameters = k - 1 + k * p + k * p * (p + 1) / 2;
    let means = model
        .means
        .iter()
        .map(|m| frame.original(m))
        .collect::<Result<Vec<_>, _>>()?;
    let covariance = model
        .covariances
        .iter()
        .map(|c| {
            (0..p)
                .flat_map(|i| (0..p).map(move |j| (i, j)))
                .map(|(i, j)| {
                    finite(
                        c[i][j] * frame.scales[i] * frame.scales[j],
                        "original component covariance",
                    )
                })
                .collect()
        })
        .collect::<Result<Vec<Vec<f64>>, _>>()?;
    let mut portable = vec![p as f64, k as f64];
    portable.extend_from_slice(&frame.origins);
    portable.extend_from_slice(&frame.offsets);
    portable.extend_from_slice(&frame.scales);
    portable.extend_from_slice(&model.weights);
    portable.extend(model.means.iter().flatten());
    portable.extend(model.covariances.iter().flatten().flatten());
    Ok(format!("{{\"version\":1,\"algorithm\":\"full_covariance_gaussian_mixture_constrained_em\",\"components\":{k},\"covariance_floor\":{floor},\"floor_units\":\"standardized_covariance_eigenvalues\",\"weights\":{},\"component_means\":{},\"component_covariances_flat\":{},\"probabilities\":{},\"labels\":{},\"sample_log_likelihood\":{},\"log_likelihood\":{likelihood},\"log_likelihood_history\":{},\"aic\":{},\"bic\":{},\"converged\":{converged},\"iterations\":{iterations},\"n_init\":{starts},\"seed\":{seed},\"predict_parameters\":{}}}",numbers(&model.weights),matrix(&means),matrix(&covariance),matrix(&probabilities),numbers(&labels(&probabilities)),numbers(&ll),numbers(&history),finite(2.*parameters as f64-2.*likelihood,"mixture AIC")?,finite(parameters as f64*(n as f64).ln()-2.*likelihood,"mixture BIC")?,numbers(&portable)))
}
fn portable_frame(
    data: &[Vec<f64>],
    args: &[f64],
    max_k: usize,
) -> Result<(Frame, usize, usize), String> {
    if args.len() < 2 {
        return Err("saved model is missing its dimension header".into());
    }
    let p = integer(args[0], 1, 64, "model features")?;
    let k = integer(args[1], 1, max_k, "model components")?;
    if data[0].len() != p || args.len() < 2 + 3 * p {
        return Err("model feature count or portable length is invalid".into());
    }
    let scales = args[2 + 2 * p..2 + 3 * p].to_vec();
    if scales.iter().any(|s| *s <= 0.) {
        return Err("model scales must be positive".into());
    }
    let mut frame = Frame {
        x: vec![],
        origins: args[2..2 + p].to_vec(),
        offsets: args[2 + p..2 + 2 * p].to_vec(),
        scales,
    };
    frame.x = frame.transform(data)?;
    Ok((frame, p, k))
}
fn mixture_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let (frame, p, k) = portable_frame(data, args, 32)?;
    if p > 32 || args.len() != 2 + 3 * p + k + k * p + k * p * p {
        return Err("invalid saved Gaussian mixture length".into());
    }
    work(data.len(), k, p * p)?;
    let start = 2 + 3 * p;
    let weights = args[start..start + k].to_vec();
    if weights.iter().any(|v| *v <= 0.) || (weights.iter().sum::<f64>() - 1.).abs() > 1e-10 {
        return Err("saved mixture weights must be positive and sum to one".into());
    }
    let means = args[start + k..start + k + k * p]
        .chunks_exact(p)
        .map(|r| r.to_vec())
        .collect();
    let covariances = args[start + k + k * p..]
        .chunks_exact(p * p)
        .map(|r| r.chunks_exact(p).map(|r| r.to_vec()).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for cov in &covariances {
        let scale = cov.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
        for (i, row) in cov.iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                if (*value - cov[j][i]).abs() > 1e-12 * scale {
                    return Err("saved covariance is not symmetric".into());
                }
            }
        }
    }
    let model = Gaussian {
        weights,
        means,
        covariances,
    };
    let (probabilities, mut ll) = mixture_evaluate(&model, &frame.x)?;
    let log_scale = frame.scales.iter().map(|s| s.ln()).sum::<f64>();
    for value in &mut ll {
        *value = finite(*value - log_scale, "prediction log likelihood")?;
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"gaussian_mixture_predict\",\"labels\":{},\"probabilities\":{},\"sample_log_likelihood\":{}}}",numbers(&labels(&probabilities)),matrix(&probabilities),numbers(&ll)))
}

fn nearest(x: &[f64], codebook: &[Vec<f64>]) -> (usize, usize, f64) {
    let mut distances = codebook
        .iter()
        .enumerate()
        .map(|(i, c)| (norm(x.iter().zip(c).map(|(a, b)| a - b)), i))
        .collect::<Vec<_>>();
    distances.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
    (
        distances[0].1,
        distances.get(1).unwrap_or(&distances[0]).1,
        distances[0].0,
    )
}
fn som(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 2 || args.len() > 6 {
        return Err("self_organizing_map needs [grid_rows,grid_columns,maxiter=100,bandwidth=0.5,standardize=1,tolerance=1e-6]".into());
    }
    let rows = integer(args[0], 1, 32, "grid rows")?;
    let columns = integer(args[1], 1, 32, "grid columns")?;
    let k = rows * columns;
    if k > 256 {
        return Err("SOM supports at most256 nodes".into());
    }
    let maxiter = integer(
        args.get(2).copied().unwrap_or(100.),
        1,
        2000,
        "SOM iterations",
    )?;
    let bandwidth = args.get(3).copied().unwrap_or(0.5);
    if bandwidth <= 0. {
        return Err("SOM bandwidth must be positive".into());
    }
    let standardize = unit(args.get(4).copied().unwrap_or(1.), "SOM standardize")?;
    let tol = tolerance(args.get(5).copied().unwrap_or(1e-6))?;
    let n = data.len();
    let p = data[0].len();
    if 2 + 3 * p + k * p > MAX_PARAMETERS {
        return Err("saved SOM would exceed the native16384-parameter limit".into());
    }
    work(n, k * p, maxiter)?;
    work(k, k * p, maxiter)?;
    let frame = Frame::new(data, standardize)?;
    let e = eigen(&cross_covariance(&frame.x, &frame.x, n)?)?;
    let grid = (0..rows)
        .flat_map(|r| {
            (0..columns).map(move |c| {
                vec![
                    if rows == 1 {
                        0.
                    } else {
                        2. * r as f64 / (rows - 1) as f64 - 1.
                    },
                    if columns == 1 {
                        0.
                    } else {
                        2. * c as f64 / (columns - 1) as f64 - 1.
                    },
                ]
            })
        })
        .collect::<Vec<_>>();
    let mut codebook = grid
        .iter()
        .map(|g| {
            (0..p)
                .map(|j| {
                    2.5 * e.values[0].max(0.).sqrt()
                        * e.vectors[j][0]
                        * if columns > 1 { g[1] } else { g[0] }
                        + if columns > 1 && rows > 1 && p > 1 {
                            2.5 * e.values[1].max(0.).sqrt() * e.vectors[j][1] * g[0]
                        } else {
                            0.
                        }
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut converged = false;
    let mut iterations = 0;
    let mut fallbacks = 0;
    for iteration in 0..maxiter {
        let assigned = frame
            .x
            .iter()
            .map(|x| nearest(x, &codebook).0)
            .collect::<Vec<_>>();
        let mut counts = vec![0_usize; k];
        let mut means = vec![vec![0.; p]; k];
        for (x, c) in frame.x.iter().zip(&assigned) {
            counts[*c] += 1;
            for j in 0..p {
                means[*c][j] += x[j];
            }
        }
        for c in 0..k {
            if counts[c] > 0 {
                for value in &mut means[c] {
                    *value /= counts[c] as f64;
                }
            }
        }
        let mut next = vec![vec![0.; p]; k];
        for target in 0..k {
            let logs = (0..k)
                .map(|c| {
                    if counts[c] == 0 {
                        f64::NEG_INFINITY
                    } else {
                        (counts[c] as f64).ln()
                            - 0.5
                                * ((grid[c][0] - grid[target][0]).powi(2)
                                    + (grid[c][1] - grid[target][1]).powi(2))
                                / bandwidth
                                / bandwidth
                    }
                })
                .collect::<Vec<_>>();
            let max = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
            let weights = logs.iter().map(|v| (v - max).exp()).collect::<Vec<_>>();
            let total = weights.iter().sum::<f64>();
            if !total.is_finite() || total == 0. {
                return Err("SOM bandwidth is outside representable range".into());
            }
            let mut local = Vec::with_capacity(k);
            for g in &grid {
                let mut z = vec![1.];
                if rows > 1 {
                    z.push(g[0] - grid[target][0]);
                }
                if columns > 1 {
                    z.push(g[1] - grid[target][1]);
                }
                local.push(z);
            }
            let q = local[0].len();
            let information = (0..q)
                .map(|i| {
                    (0..q)
                        .map(|j| {
                            (0..k)
                                .map(|c| weights[c] * local[c][i] * local[c][j] / total)
                                .sum()
                        })
                        .collect()
                })
                .collect::<Vec<Vec<f64>>>();
            match symmetric_power(&information, -1.) {
                Ok(inverse) => {
                    for j in 0..p {
                        let rhs = (0..q)
                            .map(|i| {
                                (0..k)
                                    .map(|c| weights[c] * local[c][i] * means[c][j] / total)
                                    .sum::<f64>()
                            })
                            .collect::<Vec<_>>();
                        next[target][j] =
                            finite(dot(&inverse[0], &rhs), "SOM local linear estimate")?;
                    }
                }
                Err(_) => {
                    fallbacks += 1;
                    for j in 0..p {
                        next[target][j] = finite(
                            (0..k).map(|c| weights[c] * means[c][j] / total).sum(),
                            "SOM local constant fallback",
                        )?;
                    }
                }
            }
        }
        let change = next
            .iter()
            .zip(&codebook)
            .map(|(a, b)| norm(a.iter().zip(b).map(|(a, b)| a - b)))
            .fold(0_f64, f64::max);
        let scale = codebook
            .iter()
            .map(|r| norm(r.iter().copied()))
            .fold(0_f64, f64::max);
        codebook = next;
        iterations = iteration + 1;
        if change <= tol * (1. + scale) {
            converged = true;
            break;
        }
    }
    let winners = frame
        .x
        .iter()
        .map(|r| nearest(r, &codebook))
        .collect::<Vec<_>>();
    let assigned = winners.iter().map(|v| v.0 as f64).collect::<Vec<_>>();
    let errors = winners.iter().map(|v| v.2).collect::<Vec<_>>();
    let topographic = winners
        .iter()
        .filter(|(a, b, _)| {
            let ar = a / columns;
            let ac = a % columns;
            let br = b / columns;
            let bc = b % columns;
            k > 1 && ar.abs_diff(br) + ac.abs_diff(bc) != 1
        })
        .count() as f64
        / n as f64;
    let original = codebook
        .iter()
        .map(|r| frame.original(r))
        .collect::<Result<Vec<_>, _>>()?;
    let quantization = finite(
        errors.iter().sum::<f64>() / n as f64,
        "SOM quantization error",
    )?;
    let mut portable = vec![p as f64, k as f64];
    portable.extend_from_slice(&frame.origins);
    portable.extend_from_slice(&frame.offsets);
    portable.extend_from_slice(&frame.scales);
    portable.extend(codebook.iter().flatten());
    Ok(format!("{{\"version\":1,\"algorithm\":\"batch_local_linear_gaussian_som\",\"grid_rows\":{rows},\"grid_columns\":{columns},\"bandwidth\":{bandwidth},\"initialization\":\"PCA_grid_plus_minus_2_5_population_standard_deviations\",\"grid_coordinates\":{},\"codebook\":{},\"labels\":{},\"quantization_error\":{quantization},\"distance_units\":\"fit_preprocessed_features\",\"topographic_error\":{topographic},\"converged\":{converged},\"iterations\":{iterations},\"local_linear_fallbacks\":{fallbacks},\"predict_parameters\":{}}}",matrix(&grid),matrix(&original),numbers(&assigned),numbers(&portable)))
}
fn som_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let (frame, p, k) = portable_frame(data, args, 256)?;
    if args.len() != 2 + 3 * p + k * p {
        return Err("invalid saved SOM length".into());
    }
    work(data.len(), k, p)?;
    let codebook = args[2 + 3 * p..]
        .chunks_exact(p)
        .map(|r| r.to_vec())
        .collect::<Vec<_>>();
    let winners = frame
        .x
        .iter()
        .map(|r| nearest(r, &codebook))
        .collect::<Vec<_>>();
    let assigned = winners.iter().map(|v| v.0 as f64).collect::<Vec<_>>();
    let distances = winners.iter().map(|v| v.2).collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"som_predict\",\"labels\":{},\"distances\":{},\"quantization_error\":{}}}",numbers(&assigned),numbers(&distances),finite(distances.iter().sum::<f64>()/data.len()as f64,"SOM prediction error")?))
}
