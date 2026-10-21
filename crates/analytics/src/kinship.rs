//! Gaussian genomic kernel ML/REML with explicit bounded variance ratio.
type Matrix = Vec<Vec<f64>>;
fn array(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("kinship result outside finite binary64 range".into())
    }
}
// Keep physical prediction terms outside binary64's exponent range until
// after summation. This also retains a small term between cancelling large
// terms; merely changing multiplication order would not do that.
#[derive(Clone, Copy)]
struct Scaled {
    mantissa: f64,
    exponent: i32,
}
impl Scaled {
    fn from(x: f64) -> Self {
        if x == 0. {
            return Self {
                mantissa: 0.,
                exponent: 0,
            };
        }
        let bits = x.to_bits();
        let exponent = ((bits >> 52) & 0x7ff) as i32;
        if exponent == 0 {
            let mut value = Self::from(x * 18014398509481984.);
            value.exponent -= 54;
            return value;
        }
        Self {
            mantissa: f64::from_bits((bits & !(0x7ff_u64 << 52)) | (1022_u64 << 52)),
            exponent: exponent - 1022,
        }
    }
    fn product(numerators: &[f64], denominators: &[f64]) -> Self {
        let mut value = Self::from(1.);
        for &x in numerators {
            if x == 0. {
                return Self::from(0.);
            }
            let factor = Self::from(x);
            let normalized = Self::from(value.mantissa * factor.mantissa);
            value = Self {
                mantissa: normalized.mantissa,
                exponent: value.exponent + factor.exponent + normalized.exponent,
            };
        }
        for &x in denominators {
            let factor = Self::from(x);
            let normalized = Self::from(value.mantissa / factor.mantissa);
            value = Self {
                mantissa: normalized.mantissa,
                exponent: value.exponent - factor.exponent + normalized.exponent,
            };
        }
        value
    }
    fn add(self, other: Self) -> (Self, Self) {
        if self.mantissa == 0. {
            return (other, self);
        }
        if other.mantissa == 0. {
            return (self, other);
        }
        let (a, b) = if self.exponent > other.exponent
            || (self.exponent == other.exponent && self.mantissa.abs() >= other.mantissa.abs())
        {
            (self, other)
        } else {
            (other, self)
        };
        let distance = a.exponent - b.exponent;
        if distance > 54 {
            return (a, b);
        }
        let low = b.mantissa * 2_f64.powi(-distance);
        let sum = a.mantissa + low;
        let error = low - (sum - a.mantissa);
        let normalize = |x| {
            let mut value = Self::from(x);
            if x != 0. {
                value.exponent += a.exponent;
            }
            value
        };
        (normalize(sum), normalize(error))
    }
    fn physical(self) -> Result<f64, String> {
        if self.mantissa == 0. || self.exponent < -1074 {
            return Ok(0.);
        }
        if self.exponent > 1024 {
            return Err("kinship result outside finite binary64 range".into());
        }
        let value = if self.exponent == 1024 {
            (self.mantissa * 2_f64.powi(1023)) * 2.
        } else if self.exponent < -1022 {
            (self.mantissa * 2_f64.powi(self.exponent + 1074)) * f64::from_bits(1)
        } else {
            self.mantissa * 2_f64.powi(self.exponent)
        };
        finite(value)
    }
}
fn physical_sum(terms: impl IntoIterator<Item = Scaled>) -> Result<f64, String> {
    let mut expansion: Vec<Scaled> = Vec::new();
    for term in terms {
        let mut value = term;
        let mut next = Vec::with_capacity(expansion.len() + 1);
        for component in expansion {
            let (sum, error) = value.add(component);
            if error.mantissa != 0. {
                next.push(error);
            }
            value = sum;
        }
        if value.mantissa != 0. {
            next.push(value);
        }
        expansion = next;
    }
    expansion
        .into_iter()
        .fold(Scaled::from(0.), |s, x| s.add(x).0)
        .physical()
}
#[allow(clippy::needless_range_loop)]
fn eigen(a: &[Vec<f64>]) -> Result<(Vec<f64>, Matrix), String> {
    let n = a.len();
    let mut a = a.to_vec();
    let mut u = vec![vec![0.; n]; n];
    for i in 0..n {
        u[i][i] = 1.;
    }
    for _ in 0..64 {
        let mut maximum = 0_f64;
        for i in 0..n {
            for j in i + 1..n {
                let off = a[i][j];
                maximum = maximum.max(off.abs());
                if off.abs() < 1e-15 {
                    continue;
                }
                let tau = (a[j][j] - a[i][i]) / (2. * off);
                let t = tau.signum() / (tau.abs() + tau.hypot(1.));
                let t = if tau == 0. { 1. } else { t };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                a[i][i] -= t * off;
                a[j][j] += t * off;
                a[i][j] = 0.;
                a[j][i] = 0.;
                for k in 0..n {
                    if k != i && k != j {
                        let x = a[k][i];
                        let y = a[k][j];
                        a[k][i] = c * x - s * y;
                        a[i][k] = a[k][i];
                        a[k][j] = s * x + c * y;
                        a[j][k] = a[k][j];
                    }
                    let x = u[k][i];
                    let y = u[k][j];
                    u[k][i] = c * x - s * y;
                    u[k][j] = s * x + c * y;
                }
            }
        }
        if maximum < 1e-12 {
            let d = (0..n).map(|i| a[i][i]).collect::<Vec<_>>();
            let high = d.iter().copied().fold(1., f64::max);
            if d.iter().any(|&v| v < -1e-12 * high) {
                return Err("kinship matrix is not positive semidefinite".into());
            }
            return Ok((d.iter().map(|&v| v.max(0.)).collect(), u));
        }
    }
    Err("kinship symmetric eigensolver did not converge".into())
}
#[allow(clippy::needless_range_loop)]
fn cholesky(a: &[Vec<f64>]) -> Result<(Matrix, f64), String> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    let mut logdet = 0.;
    for i in 0..n {
        for j in 0..=i {
            let value = a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if value <= 1e-12 * a[i][i] || !value.is_finite() {
                    return Err("whitened fixed design is rank deficient or ill-conditioned".into());
                }
                l[i][j] = value.sqrt();
                logdet += value.ln();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    Ok((l, logdet))
}
fn solve(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let n = l.len();
    let mut x = b.to_vec();
    for i in 0..n {
        x[i] = (x[i] - (0..i).map(|j| l[i][j] * x[j]).sum::<f64>()) / l[i][i];
    }
    for i in (0..n).rev() {
        x[i] = (x[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
    }
    x
}
struct Profile {
    objective: f64,
    beta: Vec<f64>,
    variance: f64,
    residual: Vec<f64>,
}
struct Model {
    n: usize,
    p: usize,
    reml: bool,
    eigenvalues: Vec<f64>,
    u: Matrix,
    rotated_x: Matrix,
    rotated_y: Vec<f64>,
}
impl Model {
    fn profile(&self, h: f64) -> Result<Profile, String> {
        let q = self.p + 1;
        let d = self
            .eigenvalues
            .iter()
            .map(|&v| 1. - h + h * v)
            .collect::<Vec<_>>();
        if d.iter().any(|&v| v <= 0. || !v.is_finite()) {
            return Err("nonpositive marginal covariance".into());
        }
        let mut gram = vec![vec![0.; q]; q];
        let mut b = vec![0.; q];
        for (i, row) in self.rotated_x.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                b[j] += v * self.rotated_y[i] / d[i];
                for (k, &other) in row.iter().enumerate() {
                    gram[j][k] += v * other / d[i];
                }
            }
        }
        let (l, logdet_fixed) = cholesky(&gram)?;
        let beta = solve(&l, &b);
        let residual = self
            .rotated_x
            .iter()
            .zip(&self.rotated_y)
            .map(|(row, &y)| y - row.iter().zip(&beta).map(|(&v, &b)| v * b).sum::<f64>())
            .collect::<Vec<_>>();
        let rss = residual
            .iter()
            .zip(&d)
            .map(|(&v, &d)| v * v / d)
            .sum::<f64>();
        let whitened_response_norm = self
            .rotated_y
            .iter()
            .zip(&d)
            .map(|(&v, &d)| v * v / d)
            .sum::<f64>();
        if rss <= 1e-24 * whitened_response_norm {
            return Err("response is numerically exactly explained; positive-scale Gaussian MLE does not exist within relative residual norm1e-12".into());
        }
        let df = self.n - if self.reml { q } else { 0 };
        let variance = rss / df as f64;
        if variance <= 0. || !variance.is_finite() {
            return Err("zero/nonfinite profiled residual variance".into());
        }
        let objective = 0.5
            * (df as f64 * ((2. * std::f64::consts::PI).ln() + 1. + variance.ln())
                + d.iter().map(|v| v.ln()).sum::<f64>()
                + if self.reml { logdet_fixed } else { 0. });
        Ok(Profile {
            objective: finite(objective)?,
            beta,
            variance,
            residual: residual.iter().zip(&d).map(|(&r, &v)| r / v).collect(),
        })
    }
}
fn fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 5 {
        return Err("kinship_model expects [fixed_covariates,method0ML/1REML=1,h_upper=0.99999999,grid=128,tolerance=1e-10]".into());
    }
    let n = data.len();
    let p = integer(args[0], 16)?;
    let reml = integer(*args.get(1).unwrap_or(&1.), 1)? == 1;
    let upper = *args.get(2).unwrap_or(&0.99999999);
    let grid = integer(*args.get(3).unwrap_or(&128.), 512)?;
    let tolerance = *args.get(4).unwrap_or(&1e-10);
    if !(4..=128).contains(&n)
        || n <= p + 2
        || data[0].len() != n + p + 1
        || !(0.01..=0.999999999999).contains(&upper)
        || grid < 16
        || !(1e-12..=1e-4).contains(&tolerance)
    {
        return Err("invalid kinship dimensions or bounded search parameters".into());
    }
    let kernel: Matrix = data.iter().map(|r| r[1 + p..].to_vec()).collect();
    let maxdiag = (0..n).map(|i| kernel[i][i]).fold(0., f64::max);
    if maxdiag <= 0. {
        return Err("kinship diagonal must have positive mean".into());
    }
    let kernel_scale = maxdiag
        * (0..n)
            .map(|i| kernel[i][i] / maxdiag / n as f64)
            .sum::<f64>();
    if kernel_scale <= 0. || !kernel_scale.is_finite() {
        return Err("invalid kernel scale".into());
    }
    for (i, row) in kernel.iter().enumerate() {
        for j in 0..n {
            if row[j] != kernel[j][i] {
                return Err("kinship matrix must be exactly symmetric; asymmetric roundoff correction is unsupported, symmetrize explicitly before fitting".into());
            }
            if kernel[i][i] < 0. || kernel[i][j].abs() / maxdiag > 1e4 {
                return Err(
                    "kinship matrix must be positive semidefinite within eigenvalue tolerance"
                        .into(),
                );
            }
        }
    }
    let normalized: Matrix = kernel
        .iter()
        .map(|r| r.iter().map(|v| v / kernel_scale).collect())
        .collect();
    let (eigenvalues, u) = eigen(&normalized)?;
    let y_origin = data[0][0];
    let y_scale = data
        .iter()
        .map(|r| (r[0] - y_origin).abs())
        .fold(0., f64::max);
    if y_scale == 0. {
        return Err("response variance is zero".into());
    }
    let y_offset = data
        .iter()
        .map(|r| ((r[0] - y_origin) / y_scale) / n as f64)
        .sum::<f64>();
    let y = data
        .iter()
        .map(|r| (r[0] - y_origin) / y_scale - y_offset)
        .collect::<Vec<_>>();
    let origins = data[0][1..1 + p].to_vec();
    let mut scales = Vec::new();
    let mut offsets = Vec::new();
    for j in 0..p {
        let scale = data
            .iter()
            .map(|r| (r[j + 1] - origins[j]).abs())
            .fold(0., f64::max);
        if scale == 0. {
            return Err("fixed feature is constant".into());
        }
        scales.push(scale);
        offsets.push(
            data.iter()
                .map(|r| ((r[j + 1] - origins[j]) / scale) / n as f64)
                .sum::<f64>(),
        );
    }
    let x: Matrix = data
        .iter()
        .map(|r| {
            std::iter::once(1.)
                .chain((0..p).map(|j| (r[j + 1] - origins[j]) / scales[j] - offsets[j]))
                .collect()
        })
        .collect();
    let rotated_y = (0..n)
        .map(|j| (0..n).map(|i| u[i][j] * y[i]).sum::<f64>())
        .collect();
    let rotated_x: Matrix = (0..n)
        .map(|j| {
            (0..p + 1)
                .map(|k| (0..n).map(|i| u[i][j] * x[i][k]).sum::<f64>())
                .collect()
        })
        .collect();
    let model = Model {
        n,
        p,
        reml,
        eigenvalues,
        u,
        rotated_x,
        rotated_y,
    };
    let mut sampled = Vec::new();
    for i in 0..=grid {
        let h = upper * i as f64 / grid as f64;
        sampled.push((h, model.profile(h)?.objective));
    }
    let (mut best_h, mut best) = sampled
        .iter()
        .copied()
        .min_by(|a, b| a.1.total_cmp(&b.1))
        .unwrap();
    let mut evaluations = sampled.len();
    let sampled_max = sampled
        .iter()
        .map(|r| r.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let flat = sampled_max - best <= 1e-9 * (1. + best.abs());
    // Endpoint samples can be closest to an interior optimum. Refine their
    // adjacent intervals as well, retaining the exact endpoints as candidates.
    for i in 0..=grid {
        let left = i.saturating_sub(1);
        let right = (i + 1).min(grid);
        if !flat && sampled[i].1 <= sampled[left].1 && sampled[i].1 <= sampled[right].1 {
            let (mut lo, mut hi) = (sampled[left].0, sampled[right].0);
            let ratio = (5_f64.sqrt() - 1.) / 2.;
            let mut c = hi - ratio * (hi - lo);
            let mut d = lo + ratio * (hi - lo);
            let mut fc = model.profile(c)?.objective;
            let mut fd = model.profile(d)?.objective;
            evaluations += 2;
            for _ in 0..100 {
                if hi - lo <= tolerance {
                    break;
                }
                if fc < fd {
                    hi = d;
                    d = c;
                    fd = fc;
                    c = hi - ratio * (hi - lo);
                    fc = model.profile(c)?.objective;
                } else {
                    lo = c;
                    c = d;
                    fc = fd;
                    d = lo + ratio * (hi - lo);
                    fd = model.profile(d)?.objective;
                }
                evaluations += 1;
                if evaluations > 30000 {
                    return Err("kinship profile evaluation limit reached".into());
                }
            }
            let candidate = if fc < fd { (c, fc) } else { (d, fd) };
            if candidate.1 < best {
                best_h = candidate.0;
                best = candidate.1;
            }
        }
    }
    let max_objective = sampled
        .iter()
        .map(|r| r.1)
        .fold(f64::NEG_INFINITY, f64::max);
    let identified = max_objective - best > 1e-9 * (1. + best.abs());
    if !identified {
        best_h = 0.;
    }
    let fit = model.profile(best_h)?;
    let at_upper = upper - best_h <= tolerance * 4.;
    let boundary = best_h <= tolerance * 4. || at_upper;
    let genetic_variance = finite(if best_h == 0. {
        0.
    } else {
        (fit.variance.ln() + best_h.ln() - kernel_scale.ln() + 2. * y_scale.ln()).exp()
    })?;
    let residual_variance = finite(fit.variance * (1. - best_h) * y_scale * y_scale)?;
    if residual_variance == 0. || (best_h > 0. && genetic_variance == 0.) {
        return Err("positive fitted variance underflows physical binary64 units; rescale response or kernel".into());
    }
    let kernel_weights = (0..n)
        .map(|i| {
            finite(
                best_h
                    * (0..n)
                        .filter(|&j| model.eigenvalues[j] > 0.)
                        .map(|j| model.u[i][j] * fit.residual[j])
                        .sum::<f64>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let genetic_values = model
        .u
        .iter()
        .map(|row| {
            physical_sum(row.iter().enumerate().map(|(j, &u)| {
                Scaled::product(
                    &[y_scale, best_h, u, model.eigenvalues[j], fit.residual[j]],
                    &[],
                )
            }))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let fixed = x
        .iter()
        .map(|r| r.iter().zip(&fit.beta).map(|(&v, &b)| v * b).sum::<f64>())
        .collect::<Vec<_>>();
    let fitted = fixed
        .iter()
        .map(|&v| finite(y_origin + y_scale * (y_offset + v)))
        .collect::<Result<Vec<_>, _>>()?;
    let conditional = fixed
        .iter()
        .zip(&genetic_values)
        .map(|(&v, &g)| finite(y_origin + y_scale * (y_offset + v + g / y_scale)))
        .collect::<Result<Vec<_>, _>>()?;
    let coefficients = fit.beta[1..]
        .iter()
        .zip(&scales)
        .map(|(&b, &s)| finite(b * y_scale / s))
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        y_origin + y_scale * (y_offset + fit.beta[0])
            - coefficients
                .iter()
                .enumerate()
                .map(|(j, &b)| b * (origins[j] + scales[j] * offsets[j]))
                .sum::<f64>(),
    )?;
    let df = n - if reml { p + 1 } else { 0 };
    let log_likelihood = finite(
        -fit.objective
            - df as f64 * y_scale.ln()
            - if reml {
                scales.iter().map(|v| v.ln()).sum::<f64>()
            } else {
                0.
            },
    )?;
    let mut params = vec![
        1.,
        n as f64,
        p as f64,
        y_origin,
        y_offset,
        y_scale,
        fit.beta[0],
        kernel_scale,
    ];
    params.extend(&origins);
    params.extend(&offsets);
    params.extend(&scales);
    params.extend(&fit.beta[1..]);
    params.extend(&kernel_weights);
    Ok(format!("{{\"version\":1,\"method\":\"kinship_model\",\"estimation\":\"{}\",\"intercept\":{},\"coefficients\":{},\"genetic_variance\":{},\"residual_variance\":{},\"kernel_mean_diagonal\":{},\"normalized_genetic_fraction\":{},\"log_likelihood\":{},\"genetic_values\":{},\"marginal_fitted\":{},\"conditional_fitted\":{},\"predict_parameters\":{},\"variance_components_identified\":{},\"boundary\":{},\"upper_bound_active\":{},\"converged\":{},\"likelihood_evaluations\":{},\"variance_parameter_uncertainty\":null,\"search\":\"bounded_grid_and_local_profile_refinement_no_global_optimum_certificate\"}}",if reml{"REML"}else{"ML"},intercept,array(&coefficients),genetic_variance,residual_variance,kernel_scale,best_h,log_likelihood,array(&genetic_values),array(&fitted),array(&conditional),array(&params),identified,boundary,at_upper,identified&&!at_upper,evaluations))
}
fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 8 || args[0] != 1. {
        return Err("invalid saved kinship model version".into());
    }
    let n = integer(args[1], 128)?;
    let p = integer(args[2], 16)?;
    if n < 4
        || args.len() != 8 + 4 * p + n
        || data[0].len() != p + n
        || args[5] <= 0.
        || args[7] <= 0.
        || args[8 + 2 * p..8 + 3 * p].iter().any(|&v| v <= 0.)
    {
        return Err("invalid saved kinship model dimensions/scales".into());
    }
    let mut predictions = Vec::new();
    let mut marginal = Vec::new();
    for r in data {
        let mut terms = vec![
            Scaled::from(args[3]),
            Scaled::product(&[args[5], args[4]], &[]),
            Scaled::product(&[args[5], args[6]], &[]),
        ];
        for j in 0..p {
            terms.push(Scaled::product(
                &[args[5], r[j] - args[8 + j], args[8 + 3 * p + j]],
                &[args[8 + 2 * p + j]],
            ));
            terms.push(Scaled::product(
                &[-args[5], args[8 + p + j], args[8 + 3 * p + j]],
                &[],
            ));
        }
        marginal.push(physical_sum(terms.iter().copied())?);
        terms.extend(
            r[p..]
                .iter()
                .zip(&args[8 + 4 * p..])
                .map(|(&v, &a)| Scaled::product(&[args[5], v, a], &[args[7]])),
        );
        predictions.push(physical_sum(terms)?);
    }
    Ok(format!("{{\"version\":1,\"method\":\"kinship_predict\",\"predictions\":{},\"marginal_predictions\":{},\"fit_scope\":\"saved_fixed_effects_and_training_kernel_weights_only\"}}",array(&predictions),array(&marginal)))
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
            r.len() != data[0].len() || r.iter().any(|v| !v.is_finite() || v.abs() > 1e150)
        })
    {
        return Err("kinship requires bounded finite rectangular observations".into());
    }
    match op {
        "kinship_model" => fit(data, args),
        "kinship_predict" => predict(data, args),
        _ => Err("unknown kinship operation".into()),
    }
}
