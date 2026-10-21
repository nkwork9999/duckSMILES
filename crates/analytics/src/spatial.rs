//! Isotropic spatial covariance, empirical variograms and known-covariance GLS.
//! Covariance parameters are supplied; their estimation uncertainty is excluded.
type Matrix = Vec<Vec<f64>>;
fn integer(x: f64, cap: usize) -> Result<usize, String> {
    if x >= 0. && x <= cap as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("invalid spatial integer parameter".into())
    }
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("spatial result exceeds binary64; rescale inputs".into())
    }
}
fn array(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn matrix(v: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|r| array(r)).collect::<Vec<_>>().join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(&x, &y)| x * y).sum()
}
fn distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).fold(0_f64, |h, (&x, &y)| h.hypot(x - y))
}
#[derive(Clone, Copy)]
struct Kernel {
    kind: usize,
    sill: f64,
    nugget: f64,
    rho: f64,
}
impl Kernel {
    fn new(a: &[f64]) -> Result<Self, String> {
        if a.len() != 4 {
            return Err("spatial kernel requires [kind,sill,nugget,rho]".into());
        }
        let kind = integer(a[0], 3)?;
        if a[1] < 1e-100
            || a[1] > 1e100
            || a[2] < 0.
            || a[2] > 1e100
            || a[3] <= 0.
            || a[3] > 1e100
            || (kind == 3 && a[3] >= 1.)
        {
            return Err("invalid positive spatial sill/range or nonnegative nugget; power rho must be below one".into());
        }
        Ok(Self {
            kind,
            sill: a[1],
            nugget: a[2],
            rho: a[3],
        })
    }
    fn corr(self, d: f64) -> f64 {
        if d == 0. {
            return 1.;
        }
        let u = d / self.rho;
        match self.kind {
            0 if u >= 1. => 0.,
            0 => (1. - u).powi(2) * (1. + 0.5 * u),
            1 => (-u).exp(),
            2 => (-(u * u)).exp(),
            _ => (d * self.rho.ln()).exp(),
        }
    }
    fn scale(self) -> f64 {
        self.sill + self.nugget
    }
    fn semivariance(self, distance: f64) -> f64 {
        let u = distance / self.rho;
        let process = match self.kind {
            0 if u >= 1. => self.sill,
            // Compute physical small-distance products before dimensionless
            // underflow; the discarded relative correction is below binary64.
            0 if u < 1e-150 => (1.5 * (self.sill / self.rho)) * distance,
            0 => self.sill * u * (1.5 - 0.5 * u * u),
            1 if u < 1e-150 => (self.sill / self.rho) * distance,
            1 => self.sill * -(-u).exp_m1(),
            2 if u < 1e-75 => (self.sill.sqrt() * u).powi(2),
            2 => self.sill * (-(u * u)).exp_m1().abs(),
            _ => {
                let log_rho = -self.rho.ln();
                let t = distance * log_rho;
                if t < 1e-150 {
                    distance * (self.sill * log_rho)
                } else {
                    self.sill * -(-t).exp_m1()
                }
            }
        };
        self.nugget + process
    }
    fn covariance(self, points: &Matrix) -> Matrix {
        points
            .iter()
            .enumerate()
            .map(|(i, a)| {
                points
                    .iter()
                    .enumerate()
                    .map(|(j, b)| {
                        self.sill * self.corr(distance(a, b))
                            + if i == j { self.nugget } else { 0. }
                    })
                    .collect()
            })
            .collect()
    }
}
fn cholesky(a: &Matrix) -> Result<Matrix, String> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        if a[i][i] <= 0. {
            return Err("spatial covariance/design must be positive definite".into());
        }
        for j in 0..=i {
            let v = a[i][j] - dot(&l[i][..j], &l[j][..j]);
            if i == j {
                if v <= 1e-12 * a[i][i] || !v.is_finite() {
                    return Err("spatial covariance/design is numerically singular".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    Ok(l)
}
fn forward(l: &Matrix, y: &[f64]) -> Vec<f64> {
    let mut x = vec![0.; y.len()];
    for i in 0..y.len() {
        x[i] = (y[i] - dot(&l[i][..i], &x[..i])) / l[i][i];
    }
    x
}
fn solve(l: &Matrix, y: &[f64]) -> Vec<f64> {
    let mut x = forward(l, y);
    for i in (0..y.len()).rev() {
        x[i] = (x[i] - (i + 1..y.len()).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
    }
    x
}
fn inverse(l: &Matrix) -> Matrix {
    let n = l.len();
    let mut a = vec![vec![0.; n]; n];
    for j in 0..n {
        let mut unit = vec![0.; n];
        unit[j] = 1.;
        for (i, &v) in solve(l, &unit).iter().enumerate() {
            a[i][j] = v;
        }
    }
    for i in 0..n {
        let (before, after) = a.split_at_mut(i);
        let row = &mut after[0];
        for (j, prior) in before.iter_mut().enumerate() {
            let v = 0.5 * (row[j] + prior[i]);
            row[j] = v;
            prior[i] = v;
        }
    }
    a
}
fn pack(a: &Matrix) -> Vec<f64> {
    a.iter()
        .enumerate()
        .flat_map(|(i, r)| r[..=i].iter().copied())
        .collect()
}
fn unpack(a: &[f64], n: usize, symmetric: bool) -> Matrix {
    let mut r = vec![vec![0.; n]; n];
    let mut k = 0;
    for (i, row) in r.iter_mut().enumerate() {
        row[..=i].copy_from_slice(&a[k..k + i + 1]);
        k += i + 1;
    }
    if symmetric {
        for i in 0..n {
            let (before, after) = r.split_at_mut(i);
            let row = &after[0];
            for (j, prior) in before.iter_mut().enumerate() {
                prior[i] = row[j];
            }
        }
    }
    r
}
fn covariance(data: &Matrix, args: &[f64]) -> Result<String, String> {
    if data.len() > 256 || data[0].len() > 3 {
        return Err("spatial covariance supports up to256 locations in1..3 dimensions".into());
    }
    let kernel = Kernel::new(args)?;
    let c = kernel.covariance(data);
    let variogram: Matrix = data
        .iter()
        .enumerate()
        .map(|(i, row)| {
            data.iter()
                .enumerate()
                .map(|(j, other)| {
                    if i == j {
                        0.
                    } else {
                        kernel.semivariance(distance(row, other))
                    }
                })
                .collect()
        })
        .collect();
    Ok(format!("{{\"version\":1,\"method\":\"spatial_covariance\",\"covariance\":{},\"semivariance\":{},\"positive_definite_at_tolerance\":{},\"nugget_scope\":\"independent_observation_records_even_when_coordinates_match\"}}",matrix(&c),matrix(&variogram),cholesky(&c).is_ok()))
}
fn variogram(data: &Matrix, args: &[f64]) -> Result<String, String> {
    let n = data.len();
    let d = data[0].len() - 1;
    if !(2..=512).contains(&n) || !(1..=3).contains(&d) || args.len() > 3 {
        return Err(
            "variogram requires2..512 rows[coordinates1..3,response]; [bins,lower,upper]".into(),
        );
    }
    let requested = integer(args.first().copied().unwrap_or(10.), 100)?;
    if requested == 0 {
        return Err("variogram bins must be positive".into());
    }
    let mut pairs = Vec::new();
    for i in 0..n {
        for j in 0..i {
            pairs.push((
                distance(&data[i][..d], &data[j][..d]),
                0.5 * (data[i][d] - data[j][d]).powi(2),
            ));
        }
    }
    let lower = args
        .get(1)
        .copied()
        .unwrap_or_else(|| pairs.iter().map(|r| r.0).fold(f64::INFINITY, f64::min));
    let upper = args
        .get(2)
        .copied()
        .unwrap_or_else(|| pairs.iter().map(|r| r.0).fold(0., f64::max));
    if lower < 0. || upper < lower {
        return Err("invalid variogram distance bounds".into());
    }
    let bins = if lower == upper { 1 } else { requested };
    let width = (upper - lower) / bins as f64;
    let edges = (0..=bins)
        .map(|i| {
            if i == bins {
                upper
            } else {
                lower + width * i as f64
            }
        })
        .collect::<Vec<_>>();
    let mut counts = vec![0.; bins];
    let mut distances = vec![0.; bins];
    let mut sums = vec![0.; bins];
    for (h, v) in pairs {
        if h < lower || h > upper {
            continue;
        }
        let i = if width == 0. || h == upper {
            bins - 1
        } else {
            ((h - lower) / width).floor() as usize
        };
        let i = i.min(bins - 1);
        counts[i] += 1.;
        distances[i] += h;
        sums[i] += v;
    }
    if counts.iter().sum::<f64>() == 0. {
        return Err("no point pairs in variogram distance bounds".into());
    }
    let optional = |v: &[f64]| {
        format!(
            "[{}]",
            v.iter()
                .zip(&counts)
                .map(|(&s, &c)| if c > 0. {
                    (s / c).to_string()
                } else {
                    "null".into()
                })
                .collect::<Vec<_>>()
                .join(",")
        )
    };
    Ok(format!("{{\"version\":1,\"method\":\"spatial_variogram\",\"edges\":{},\"pair_counts\":{},\"mean_distances\":{},\"semivariances\":{},\"estimator\":\"Matheron_half_mean_squared_pair_difference_no_trend_removal\"}}",array(&edges),array(&counts),optional(&distances),optional(&sums)))
}
fn gls(data: &Matrix, args: &[f64]) -> Result<String, String> {
    if args.len() != 6 {
        return Err("spatial GLS args [dimensions,fixed_features,kind,sill,nugget,rho]".into());
    }
    let d = integer(args[0], 3)?;
    let p = integer(args[1], 12)?;
    let q = p + 1;
    let n = data.len();
    if d == 0 || n <= q || n > 128 || data[0].len() != d + p + 1 {
        return Err(
            "spatial GLS requires n>p+1,n<=128, rows[coords,fixed_features,response]".into(),
        );
    }
    let kernel = Kernel::new(&args[2..])?;
    let scale = kernel.scale();
    let points: Matrix = data.iter().map(|r| r[..d].to_vec()).collect();
    let c = kernel.covariance(&points);
    let normalized: Matrix = c
        .iter()
        .map(|r| r.iter().map(|v| v / scale).collect())
        .collect();
    let l = cholesky(&normalized)?;
    let origins = data[0][d..d + p].to_vec();
    let mut scales = Vec::new();
    let mut offsets = Vec::new();
    for j in 0..p {
        let s = data
            .iter()
            .map(|r| (r[d + j] - origins[j]).abs())
            .fold(0., f64::max);
        if s == 0. {
            return Err("constant spatial fixed feature is aliased with intercept".into());
        }
        scales.push(s);
        offsets.push(
            data.iter()
                .map(|r| (r[d + j] - origins[j]) / s / n as f64)
                .sum::<f64>(),
        );
    }
    let x: Matrix = data
        .iter()
        .map(|r| {
            std::iter::once(1.)
                .chain((0..p).map(|j| (r[d + j] - origins[j]) / scales[j] - offsets[j]))
                .collect()
        })
        .collect();
    let yo = data[0][d + p];
    let ys = data
        .iter()
        .map(|r| (r[d + p] - yo).abs())
        .fold(0., f64::max)
        .max(1e-100);
    let yoff = data
        .iter()
        .map(|r| (r[d + p] - yo) / ys / n as f64)
        .sum::<f64>();
    let y = data
        .iter()
        .map(|r| (r[d + p] - yo) / ys - yoff)
        .collect::<Vec<_>>();
    let wy = forward(&l, &y);
    let wx = (0..q)
        .map(|j| forward(&l, &x.iter().map(|r| r[j]).collect::<Vec<_>>()))
        .collect::<Matrix>();
    let gram = wx
        .iter()
        .map(|a| wx.iter().map(|b| dot(a, b)).collect())
        .collect();
    let gl = cholesky(&gram)?;
    let gi = inverse(&gl);
    let beta = solve(&gl, &wx.iter().map(|v| dot(v, &wy)).collect::<Vec<_>>());
    let residual = y
        .iter()
        .zip(&x)
        .map(|(&v, r)| v - dot(r, &beta))
        .collect::<Vec<_>>();
    let wr = forward(&l, &residual);
    let sse = finite(dot(&wr, &wr) * ys * ys / scale)?;
    let ll = finite(
        -0.5 * (n as f64 * ((2. * std::f64::consts::PI).ln() + scale.ln())
            + 2. * l.iter().enumerate().map(|(i, r)| r[i].ln()).sum::<f64>()
            + sse),
    )?;
    let alpha = solve(&l, &residual)
        .iter()
        .map(|v| finite(v * ys))
        .collect::<Result<Vec<_>, _>>()?;
    let mut b = beta
        .iter()
        .map(|v| finite(v * ys))
        .collect::<Result<Vec<_>, _>>()?;
    b[0] = finite(yo + ys * (yoff + beta[0]))?;
    let coefficients = (0..p)
        .map(|j| finite(b[j + 1] / scales[j]))
        .collect::<Result<Vec<_>, _>>()?;
    let means = (0..p)
        .map(|j| origins[j] + offsets[j] * scales[j])
        .collect::<Vec<_>>();
    let intercept = finite(b[0] - dot(&coefficients, &means))?;
    let mut t = vec![vec![0.; q]; q];
    t[0][0] = 1.;
    for j in 0..p {
        t[0][j + 1] = -means[j] / scales[j];
        t[j + 1][j + 1] = 1. / scales[j];
    }
    // Apply the physical covariance scale on each factor before multiplying;
    // the unscaled quadratic form can overflow although the answer is finite.
    for row in &mut t {
        for value in row {
            *value *= scale.sqrt();
        }
    }
    let mut cov = vec![vec![0.; q]; q];
    for i in 0..q {
        for j in 0..q {
            cov[i][j] = finite((0..q).map(|k| t[i][k] * dot(&gi[k], &t[j])).sum::<f64>())?;
        }
    }
    let marginal = x
        .iter()
        .map(|r| finite(dot(r, &b)))
        .collect::<Result<Vec<_>, _>>()?;
    let conditional = marginal
        .iter()
        .enumerate()
        .map(|(i, &v)| {
            finite(
                v + kernel.sill / scale
                    * (0..n)
                        .map(|j| kernel.corr(distance(&points[i], &points[j])) * alpha[j])
                        .sum::<f64>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut params = vec![
        1.,
        n as f64,
        d as f64,
        p as f64,
        kernel.kind as f64,
        kernel.sill,
        kernel.nugget,
        kernel.rho,
    ];
    params.extend(origins);
    params.extend(offsets);
    params.extend(scales);
    params.extend(&b);
    params.extend(points.iter().flatten());
    params.extend(x.iter().flatten());
    params.extend(alpha);
    params.extend(pack(&l));
    params.extend(pack(&gi));
    Ok(format!("{{\"version\":1,\"method\":\"spatial_gls\",\"intercept\":{},\"coefficients\":{},\"coefficient_covariance\":{},\"marginal_fitted\":{},\"conditional_fitted\":{},\"mahalanobis_sse\":{},\"log_likelihood\":{},\"predict_parameters\":{},\"covariance_parameters_estimated\":false,\"inference_scope\":\"Gaussian_known_absolute_covariance_no_residual_scale_or_covariance_parameter_estimation\"}}",intercept,array(&coefficients),matrix(&cov),array(&marginal),array(&conditional),sse,ll,array(&params)))
}
fn predict(data: &Matrix, args: &[f64]) -> Result<String, String> {
    if args.len() < 8 || args[0] != 1. {
        return Err("invalid saved spatial model version".into());
    }
    let n = integer(args[1], 128)?;
    let d = integer(args[2], 3)?;
    let p = integer(args[3], 12)?;
    let q = p + 1;
    let len = 8 + 3 * p + q + n * (d + q + 1) + n * (n + 1) / 2 + q * (q + 1) / 2;
    if d == 0 || n <= q || args.len() != len || data[0].len() != d + p {
        return Err("invalid saved spatial model/query dimensions".into());
    }
    let kernel = Kernel::new(&args[4..8])?;
    let scale = kernel.scale();
    let origins = &args[8..8 + p];
    let offsets = &args[8 + p..8 + 2 * p];
    let scales = &args[8 + 2 * p..8 + 3 * p];
    if scales.iter().any(|&s| s <= 0.) {
        return Err("invalid saved spatial feature scale".into());
    }
    let mut pos = 8 + 3 * p;
    let b = &args[pos..pos + q];
    pos += q;
    let points = args[pos..pos + n * d]
        .chunks_exact(d)
        .map(<[f64]>::to_vec)
        .collect::<Matrix>();
    pos += n * d;
    let x = args[pos..pos + n * q]
        .chunks_exact(q)
        .map(<[f64]>::to_vec)
        .collect::<Matrix>();
    pos += n * q;
    let alpha = &args[pos..pos + n];
    pos += n;
    let l = unpack(&args[pos..pos + n * (n + 1) / 2], n, false);
    pos += n * (n + 1) / 2;
    if l.iter().enumerate().any(|(i, r)| r[i] <= 0.) {
        return Err("invalid saved spatial covariance factor".into());
    }
    let gi = unpack(&args[pos..], q, true);
    cholesky(&gi)?;
    let mut predictions = Vec::new();
    let mut latent = Vec::new();
    let mut observation = Vec::new();
    for r in data {
        let query = std::iter::once(1.)
            .chain((0..p).map(|j| (r[d + j] - origins[j]) / scales[j] - offsets[j]))
            .collect::<Vec<_>>();
        let k = points
            .iter()
            .map(|s| kernel.sill / scale * kernel.corr(distance(&r[..d], s)))
            .collect::<Vec<_>>();
        predictions.push(finite(dot(&query, b) + dot(&k, alpha))?);
        let v = forward(&l, &k);
        let a = solve(&l, &k);
        let delta = (0..q)
            .map(|j| query[j] - (0..n).map(|i| x[i][j] * a[i]).sum::<f64>())
            .collect::<Vec<_>>();
        let correction = scale * (0..q).map(|j| delta[j] * dot(&gi[j], &delta)).sum::<f64>();
        let removed = scale * dot(&v, &v);
        let variance = finite(kernel.sill - removed + correction)?;
        let comparison_scale = kernel.sill.max(removed.abs()).max(correction.abs());
        let relative_tolerance = 1e-9
            * (kernel.sill / comparison_scale
                + removed.abs() / comparison_scale
                + correction.abs() / comparison_scale);
        if variance / comparison_scale < -relative_tolerance {
            return Err("saved spatial model implies negative predictive variance".into());
        }
        // With no measurement nugget, an observed record with the same trend
        // row has known latent value. Preserve this structural zero instead of
        // magnifying a squared rounding residual by a large physical sill.
        let exact_record = kernel.nugget == 0.
            && points
                .iter()
                .zip(&x)
                .any(|(point, design)| point.as_slice() == &r[..d] && design == &query);
        let variance = if exact_record { 0. } else { variance.max(0.) };
        latent.push(variance);
        observation.push(finite(variance + kernel.nugget)?);
    }
    Ok(format!("{{\"version\":1,\"method\":\"spatial_predict\",\"predictions\":{},\"latent_variances\":{},\"observation_variances\":{},\"variance_scope\":\"universal_kriging_with_fixed_effect_uncertainty_given_known_covariance_independent_query_nugget\"}}",array(&predictions),array(&latent),array(&observation)))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 10000
        || data[0].len() > 1024
        || data.len().saturating_mul(data[0].len()) > 1_000_000
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data.iter().any(|r| {
            r.len() != data[0].len() || r.iter().any(|v| !v.is_finite() || v.abs() > 1e100)
        })
    {
        return Err("spatial operations require bounded finite rectangular inputs".into());
    }
    let data = data.to_vec();
    match op {
        "spatial_covariance" => covariance(&data, args),
        "spatial_variogram" => variogram(&data, args),
        "spatial_gls" => gls(&data, args),
        "spatial_predict" => predict(&data, args),
        _ => Err("unknown spatial operation".into()),
    }
}
