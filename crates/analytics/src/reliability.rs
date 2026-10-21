//! Right-censored parametric AFT regression, power-law NHPP and independent blocks.
//!
//! AFT models log time location as an intercept plus linear covariate effects,
//! with common SEV scale (Weibull), common normal scale (lognormal), or fixed SEV
//! scale1 (exponential). Exact failure densities and right-censor survival terms
//! form the likelihood. No left/interval censoring, delayed entry, competing
//! risks, scale covariates, frailty, repair simulation or general RBD topology.

const MAX_ROWS: usize = 100_000;
const MAX_CELLS: usize = 1_000_000;
const MAX_WORK: usize = 200_000_000;
const LOG_SQRT_TAU: f64 = 0.918_938_533_204_672_7;
type Matrix = Vec<Vec<f64>>;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > MAX_ROWS
        || p == 0
        || p > 1024
        || data.len() > MAX_CELLS / p
        || data
            .iter()
            .any(|r| r.len() != p || r.iter().any(|v| !v.is_finite()))
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("reliability requires finite rectangular nonempty data within100000rows,1024columns,1000000cells and16384arguments".into());
    }
    match operation {
        "aft_exponential" => fit(data, args, 0),
        "aft_weibull" => fit(data, args, 1),
        "aft_lognormal" => fit(data, args, 2),
        "aft_predict" => predict(data, args),
        "nhpp_powerlaw" => nhpp(data, args),
        "nhpp_predict" => nhpp_predict(data, args),
        "reliability_blocks" => blocks(data, args),
        _ => Err(format!("unknown reliability operation: {operation}")),
    }
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn integer(x: f64, min: usize, max: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < min as f64 || x > max as f64 {
        Err(format!("{name} must be integer in[{min},{max}]"))
    } else {
        Ok(x as usize)
    }
}
fn arguments(a: &[f64], min: usize, max: usize) -> Result<(), String> {
    if (min..=max).contains(&a.len()) {
        Ok(())
    } else {
        Err(format!("operation requires{min}to{max}arguments"))
    }
}
fn vector(v: &[f64]) -> Result<String, String> {
    if v.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite reliability result".into());
    }
    Ok(format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    ))
}
fn matrix(v: &[Vec<f64>]) -> Result<String, String> {
    Ok(format!(
        "[{}]",
        v.iter()
            .map(|r| vector(r))
            .collect::<Result<Vec<_>, _>>()?
            .join(",")
    ))
}
fn optional_vector(v: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|v| v.map_or_else(|| "null".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn optional_exp(x: f64) -> Option<f64> {
    let value = x.exp();
    (value.is_finite() && value > 0.).then_some(value)
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(a.iter().zip(b).map(|(a, b)| a * b).sum(), "inner product")
}
fn center(values: &[f64]) -> Result<(f64, f64, Vec<f64>), String> {
    let origin = values[0];
    let mut offset = 0.;
    for (i, &v) in values.iter().enumerate() {
        let d = finite(v - origin, "origin-relative feature")?;
        offset = finite(offset + (d - offset) / (i + 1) as f64, "feature mean")?;
    }
    Ok((
        origin,
        offset,
        values
            .iter()
            .map(|v| finite((v - origin) - offset, "centered feature"))
            .collect::<Result<Vec<_>, _>>()?,
    ))
}
fn log_ratio(a: f64, b: f64) -> f64 {
    if a > b * 0.5 && a < b * 2. {
        ((a - b) / b).ln_1p()
    } else {
        a.ln() - b.ln()
    }
}
fn identity(n: usize, scale: f64) -> Matrix {
    let mut m = vec![vec![0.; n]; n];
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = scale;
    }
    m
}
fn multiply(a: &Matrix, b: &Matrix) -> Result<Matrix, String> {
    let n = a.len();
    let mut out = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..n {
            out[i][j] = finite(
                (0..n).map(|k| a[i][k] * b[k][j]).sum(),
                "covariance transform",
            )?;
        }
    }
    Ok(out)
}
fn cholesky(a: &Matrix) -> Result<Matrix, String> {
    let n = a.len();
    let scale = (0..n).map(|j| a[j][j].abs()).fold(0., f64::max);
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let value = a[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if !value.is_finite() || value <= scale * 1e-12 {
                    return Err(
                        "observed information is singular, indefinite or numerically unidentified"
                            .into(),
                    );
                }
                l[i][j] = value.sqrt();
            } else {
                l[i][j] = finite(value / l[j][j], "information Cholesky")?;
            }
        }
    }
    Ok(l)
}
fn inverse_info(a: &Matrix) -> Result<Matrix, String> {
    let l = cholesky(a)?;
    let n = a.len();
    let cols = (0..n)
        .map(|j| {
            let mut v = vec![0.; n];
            v[j] = 1.;
            for i in 0..n {
                v[i] = (v[i] - (0..i).map(|k| l[i][k] * v[k]).sum::<f64>()) / l[i][i];
            }
            for i in (0..n).rev() {
                v[i] = finite(
                    (v[i] - (i + 1..n).map(|k| l[k][i] * v[k]).sum::<f64>()) / l[i][i],
                    "inverse observed information",
                )?;
            }
            Ok(v)
        })
        .collect::<Result<Vec<_>, String>>()?;
    Ok((0..n)
        .map(|i| cols.iter().map(|c| c[i]).collect())
        .collect())
}
struct Problem {
    design: Matrix,
    logtime: Vec<f64>,
    logtime_relative: Vec<f64>,
    events: Vec<f64>,
    deaths: usize,
    reference: f64,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    distribution: usize,
}
struct Evaluation {
    ll: f64,
    gradient: Vec<f64>,
    information: Matrix,
}
struct Budget(usize);
impl Budget {
    fn charge(&mut self, n: usize, d: usize) -> Result<(), String> {
        let cost = n
            .checked_mul(d)
            .and_then(|v| v.checked_mul(d))
            .ok_or("AFT work overflow")?;
        self.0 = self.0.checked_add(cost).ok_or("AFT work overflow")?;
        if self.0 > MAX_WORK {
            Err("AFT evaluations exceed200000000work units".into())
        } else {
            Ok(())
        }
    }
}
impl Problem {
    fn new(data: &[Vec<f64>], distribution: usize) -> Result<Self, String> {
        let p = data[0]
            .len()
            .checked_sub(2)
            .ok_or("AFT data rows are[features...,positive_time,event0or1]")?;
        if p > 16 {
            return Err("AFT supports at most16linear covariates".into());
        }
        let d = p + 1 + usize::from(distribution > 0);
        if data.len() <= d {
            return Err("AFT requires more observations than estimated parameters".into());
        }
        if data
            .iter()
            .any(|r| r[p] <= 0. || r[p + 1] != 0. && r[p + 1] != 1.)
        {
            return Err("AFT lifetimes must be positive; exact failure1/right censor0".into());
        }
        let deaths = data.iter().filter(|r| r[p + 1] == 1.).count();
        if deaths < if distribution == 0 { 1 } else { 2 } {
            return Err(
                "AFT requires at least one exponential or two shape-estimating failures".into(),
            );
        }
        let reference = data[0][p].ln();
        let logtime = data.iter().map(|r| r[p].ln()).collect::<Vec<_>>();
        let logtime_relative = data
            .iter()
            .map(|r| log_ratio(r[p], data[0][p]))
            .collect::<Vec<_>>();
        let mut design = vec![vec![1.; p + 1]; data.len()];
        let mut origins = Vec::new();
        let mut offsets = Vec::new();
        let mut scales = Vec::new();
        for j in 0..p {
            let values = data.iter().map(|r| r[j]).collect::<Vec<_>>();
            let (o, m, c) = center(&values)?;
            let largest = c.iter().map(|v| v.abs()).fold(0., f64::max);
            if largest == 0. {
                return Err(
                    "constant AFT covariates are not identifiable with an intercept".into(),
                );
            }
            let scale = finite(
                largest
                    * (c.iter().map(|v| (v / largest).powi(2)).sum::<f64>() / data.len() as f64)
                        .sqrt(),
                "feature scale",
            )?;
            if scale == 0. {
                return Err("feature scale underflow".into());
            }
            for i in 0..data.len() {
                design[i][j + 1] = c[i] / scale;
            }
            origins.push(o);
            offsets.push(m);
            scales.push(scale);
        }
        // Check the supplied fixed-effect design independently of the fitted weights.
        let mut gram = vec![vec![0.; p + 1]; p + 1];
        for row in &design {
            for j in 0..=p {
                for k in 0..=p {
                    gram[j][k] += row[j] * row[k];
                }
            }
        }
        cholesky(&gram)?;
        Ok(Self {
            design,
            logtime,
            logtime_relative,
            events: data.iter().map(|r| r[p + 1]).collect(),
            deaths,
            reference,
            origins,
            offsets,
            scales,
            distribution,
        })
    }
    fn evaluate(&self, theta: &[f64], budget: &mut Budget) -> Result<Evaluation, String> {
        let m = self.scales.len() + 1;
        let d = theta.len();
        budget.charge(self.design.len(), d)?;
        let h = if self.distribution == 0 { 0. } else { theta[m] };
        if !(-30. ..=30.).contains(&h) {
            return Err("AFT log sigma outside supported[-30,30] identifiable search range".into());
        }
        let sigma = h.exp();
        let mut ll = 0.;
        let mut gradient = vec![0.; d];
        let mut info = vec![vec![0.; d]; d];
        for i in 0..self.design.len() {
            let mu = dot(&theta[..m], &self.design[i])?;
            let z = finite(
                (self.logtime_relative[i] - mu) / sigma,
                "standardized log lifetime",
            )?;
            let event = self.events[i];
            let (value, gmu, gh, imm, imh, ihh) = if self.distribution < 2 {
                let e = finite(z.exp(), "Weibull cumulative hazard")?;
                let a = e - event;
                (
                    event * (-self.logtime[i] - h + z) - e,
                    a / sigma,
                    a * z - event,
                    e / (sigma * sigma),
                    (e * z + a) / sigma,
                    e * z * z + a * z,
                )
            } else if event == 1. {
                (
                    -self.logtime[i] - h - 0.5 * z * z - LOG_SQRT_TAU,
                    z / sigma,
                    z * z - 1.,
                    1. / (sigma * sigma),
                    2. * z / sigma,
                    2. * z * z,
                )
            } else {
                let (logsurvival, r, rprime) = normal_tail(z)?;
                (
                    logsurvival,
                    r / sigma,
                    r * z,
                    rprime / (sigma * sigma),
                    (rprime * z + r) / sigma,
                    (rprime * z + r) * z,
                )
            };
            ll += value;
            for j in 0..m {
                gradient[j] += gmu * self.design[i][j];
                for (cell, value) in info[j].iter_mut().take(m).zip(&self.design[i]) {
                    *cell += imm * self.design[i][j] * value;
                }
                if self.distribution > 0 {
                    info[j][m] += imh * self.design[i][j];
                    info[m][j] = info[j][m];
                }
            }
            if self.distribution > 0 {
                gradient[m] += gh;
                info[m][m] += ihh;
            }
        }
        finite(ll, "AFT log likelihood")?;
        if gradient
            .iter()
            .chain(info.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err("AFT derivatives exceeded finite range".into());
        }
        Ok(Evaluation {
            ll,
            gradient,
            information: info,
        })
    }
}
fn penalized(
    e: &Evaluation,
    theta: &[f64],
    p: usize,
    lambda: f64,
) -> Result<(f64, Vec<f64>), String> {
    let value = finite(
        e.ll - 0.5 * lambda * theta[1..=p].iter().map(|v| v * v).sum::<f64>(),
        "penalized AFT likelihood",
    )?;
    let mut g = e.gradient.clone();
    for j in 1..=p {
        g[j] -= lambda * theta[j];
    }
    Ok((value, g))
}
fn fit(data: &[Vec<f64>], args: &[f64], distribution: usize) -> Result<String, String> {
    arguments(args, 0, 3)?;
    let lambda = args.first().copied().unwrap_or(0.);
    if !(0. ..=1e12).contains(&lambda) {
        return Err("AFT ridge lambda must be in[0,1e12]".into());
    }
    let maxiter = integer(
        args.get(1).copied().unwrap_or(200.),
        1,
        1000,
        "AFT maximum iterations",
    )?;
    let tolerance = args.get(2).copied().unwrap_or(1e-8);
    if !(1e-12..=1e-3).contains(&tolerance) {
        return Err("AFT tolerance must be in[1e-12,1e-3]".into());
    }
    let problem = Problem::new(data, distribution)?;
    let p = problem.scales.len();
    let d = p + 1 + usize::from(distribution > 0);
    let n = data.len();
    let mut theta = vec![0.; d];
    theta[0] = problem.logtime_relative.iter().sum::<f64>() / n as f64;
    if distribution > 0 {
        let spread = (problem
            .logtime_relative
            .iter()
            .map(|v| (v - theta[0]).powi(2))
            .sum::<f64>()
            / n as f64)
            .sqrt();
        if spread == 0. {
            return Err("lifetime shape is unidentified without log-time variation".into());
        }
        theta[p + 1] = spread.max(0.1).ln();
    }
    let mut budget = Budget(0);
    let mut e = problem.evaluate(&theta, &mut budget)?;
    let mut inverse = identity(d, 1. / n as f64);
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..maxiter {
        iterations = iteration + 1;
        let (value, g) = penalized(&e, &theta, p, lambda)?;
        let mut direction = inverse
            .iter()
            .map(|r| dot(r, &g))
            .collect::<Result<Vec<_>, _>>()?;
        let mut slope = dot(&direction, &g)?;
        if slope <= 0. {
            inverse = identity(d, 1. / n as f64);
            direction = g.iter().map(|v| v / n as f64).collect();
            slope = dot(&direction, &g)?;
        }
        let gradient = g.iter().map(|v| v.abs()).fold(0., f64::max);
        let maxstep = direction.iter().map(|v| v.abs()).fold(0., f64::max);
        let maxparameter = theta.iter().map(|v| v.abs()).fold(0., f64::max);
        if gradient <= tolerance * (1. + problem.deaths as f64)
            && maxstep <= tolerance * (1. + maxparameter)
        {
            converged = true;
            break;
        }
        let mut step = if maxstep > 5. { 5. / maxstep } else { 1. };
        let mut accepted = None;
        for _ in 0..50 {
            let next = theta
                .iter()
                .zip(&direction)
                .map(|(b, v)| b + step * v)
                .collect::<Vec<_>>();
            if let Ok(next_eval) = problem.evaluate(&next, &mut budget) {
                let (next_value, next_g) = penalized(&next_eval, &next, p, lambda)?;
                if next_value >= value + 1e-4 * step * slope - 1e-12 {
                    accepted = Some((next, next_eval, next_g));
                    break;
                }
            }
            step *= 0.5;
        }
        let (next, next_eval, next_g) =
            accepted.ok_or("AFT line search failed; model may be unidentified")?;
        let s = next
            .iter()
            .zip(&theta)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let y = g
            .iter()
            .zip(&next_g)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let ys = dot(&y, &s)?;
        let norm = dot(&s, &s)?.sqrt() * dot(&y, &y)?.sqrt();
        if ys > 1e-12 * norm && ys > 0. {
            let hy = inverse
                .iter()
                .map(|r| dot(r, &y))
                .collect::<Result<Vec<_>, _>>()?;
            let yhy = dot(&y, &hy)?;
            let rho = 1. / ys;
            for i in 0..d {
                for j in 0..d {
                    inverse[i][j] = finite(
                        inverse[i][j] + (rho + rho * rho * yhy) * s[i] * s[j]
                            - rho * (hy[i] * s[j] + s[i] * hy[j]),
                        "BFGS inverse curvature",
                    )?;
                }
            }
        } else {
            inverse = identity(d, 1. / n as f64);
        }
        theta = next;
        e = next_eval;
    }
    if !converged {
        return Err("AFT fit did not converge within iteration limit".into());
    }
    let mut information = e.information.clone();
    for (j, row) in information.iter_mut().enumerate().take(p + 1).skip(1) {
        row[j] += lambda;
    }
    let standardized_covariance = inverse_info(&information)?;
    let mut transform = identity(d, 1.);
    let mut coefficients = vec![problem.reference + theta[0]];
    for j in 0..p {
        let b = finite(theta[j + 1] / problem.scales[j], "original-unit AFT slope")?;
        coefficients.push(b);
        coefficients[0] -= b * (problem.origins[j] + problem.offsets[j]);
        transform[0][j + 1] = finite(
            -(problem.origins[j] + problem.offsets[j]) / problem.scales[j],
            "intercept covariance transform",
        )?;
        transform[j + 1][j + 1] = 1. / problem.scales[j];
    }
    finite(coefficients[0], "original-unit intercept")?;
    let temp = multiply(&transform, &standardized_covariance)?;
    let transpose = (0..d)
        .map(|j| transform.iter().map(|r| r[j]).collect())
        .collect::<Matrix>();
    let covariance = multiply(&temp, &transpose)?;
    let se = (0..d)
        .map(|j| finite(covariance[j][j].sqrt(), "AFT standard error"))
        .collect::<Result<Vec<_>, _>>()?;
    let h = if distribution == 0 { 0. } else { theta[p + 1] };
    let mut parameters = coefficients.clone();
    if distribution > 0 {
        parameters.push(h);
    }
    let z = parameters
        .iter()
        .zip(&se)
        .map(|(b, s)| finite(b / s, "AFT Wald statistic"))
        .collect::<Result<Vec<_>, _>>()?;
    let regularized = lambda > 0. && p > 0;
    let (pvalues, intervals, inference) = if regularized {
        (
            format!("[{}]", vec!["null"; d].join(",")),
            "null".into(),
            "penalized_curvature_no_calibrated_inference",
        )
    } else {
        let values = z
            .iter()
            .map(|v| normal_two_sided(*v))
            .collect::<Result<Vec<_>, _>>()?;
        let intervals = parameters
            .iter()
            .zip(&se)
            .map(|(b, s)| {
                Ok(vec![
                    finite(b - 1.959_963_984_540_054 * s, "AFT lower Wald limit")?,
                    finite(b + 1.959_963_984_540_054 * s, "AFT upper Wald limit")?,
                ])
            })
            .collect::<Result<Matrix, String>>()?;
        (
            vector(&values)?,
            matrix(&intervals)?,
            "observed_information_wald_asymptotic",
        )
    };
    let mut model = vec![
        distribution as f64,
        p as f64,
        problem.reference,
        h,
        theta[0],
    ];
    model.extend(&problem.origins);
    model.extend(&problem.offsets);
    model.extend(&problem.scales);
    model.extend(&theta[1..=p]);
    let name = ["exponential", "weibull", "lognormal"][distribution];
    let aic = if regularized {
        "null".into()
    } else {
        finite(2. * d as f64 - 2. * e.ll, "AFT AIC")?.to_string()
    };
    let bic = if regularized {
        "null".into()
    } else {
        finite((n as f64).ln() * d as f64 - 2. * e.ll, "AFT BIC")?.to_string()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"right_censored_{name}_accelerated_failure_time\",\"distribution\":\"{name}\",\"n_samples\":{n},\"n_events\":{},\"n_features\":{p},\"coefficients\":{},\"parameter_order\":\"intercept_then_original_unit_slopes_then_log_sigma_if_estimated\",\"log_sigma\":{h},\"sigma\":{},\"covariance_matrix\":{},\"standard_errors\":{},\"z_statistics\":{},\"p_values\":{pvalues},\"confidence_intervals\":{intervals},\"confidence_level\":0.95,\"inference_method\":\"{inference}\",\"log_likelihood\":{},\"aic\":{aic},\"bic\":{bic},\"parameter_count\":{d},\"lambda\":{lambda},\"penalty_scale\":\"sum_loglik_minus_lambda_half_squared_standardized_slopes\",\"predict_parameters\":{},\"converged\":true,\"iterations\":{iterations},\"work_units\":{}}}",problem.deaths,vector(&coefficients)?,h.exp(),matrix(&covariance)?,vector(&se)?,vector(&z)?,e.ll,vector(&model)?,budget.0))
}
fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 6 {
        return Err("AFT prediction args are[mode,predict_parameters...]".into());
    }
    let mode = integer(args[0], 0, 1, "prediction mode")?;
    let distribution = integer(args[1], 0, 2, "saved distribution")?;
    let p = integer(args[2], 0, 16, "saved feature count")?;
    if args.len() != 6 + 4 * p || data[0].len() != p + 1 {
        return Err("AFT saved model length or feature-plus-query width mismatch".into());
    }
    let reference = args[3];
    let h = args[4];
    if !(-30. ..=30.).contains(&h) || distribution == 0 && h != 0. {
        return Err("invalid saved AFT log sigma".into());
    }
    let intercept = args[5];
    let origins = &args[6..6 + p];
    let offsets = &args[6 + p..6 + 2 * p];
    let scales = &args[6 + 2 * p..6 + 3 * p];
    let slopes = &args[6 + 3 * p..];
    if scales.iter().any(|&v| v <= 0.) {
        return Err("saved AFT feature scales must be positive".into());
    }
    let sigma = h.exp();
    let mut locations = Vec::new();
    let mut survival = Vec::new();
    let mut quantiles = Vec::new();
    let mut logquantiles = Vec::new();
    for row in data {
        let mu = finite(
            reference
                + (intercept
                    + (0..p)
                        .map(|j| slopes[j] * (((row[j] - origins[j]) - offsets[j]) / scales[j]))
                        .sum::<f64>()),
            "predicted AFT log-time location",
        )?;
        locations.push(mu);
        let query = row[p];
        if mode == 0 {
            if query < 0. {
                return Err("survival query time must be nonnegative".into());
            }
            let s = if query == 0. {
                1.
            } else {
                let z = finite((query.ln() - mu) / sigma, "predicted standardized log time")?;
                if distribution < 2 {
                    if z > 709. {
                        0.
                    } else {
                        (-z.exp()).exp()
                    }
                } else {
                    normal_tail(z)?.0.exp()
                }
            };
            survival.push(s);
        } else {
            if !(0. ..1.).contains(&query) || query == 0. {
                return Err("quantile probability must be strictly between0and1".into());
            }
            let q = if distribution < 2 {
                mu + sigma * (-(-query).ln_1p()).ln()
            } else {
                mu + sigma * normal_quantile(query)?
            };
            let q = finite(q, "predicted log quantile")?;
            quantiles.push(optional_exp(q));
            logquantiles.push(q);
        }
    }
    let extra = if mode == 0 {
        format!(",\"survival_probabilities\":{}", vector(&survival)?)
    } else {
        format!(
            ",\"time_quantiles\":{},\"log_time_quantiles\":{}",
            optional_vector(&quantiles),
            vector(&logquantiles)?
        )
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"saved_aft_distribution_prediction\",\"mode\":{mode},\"log_time_locations\":{}{extra},\"unrepresentable_positive_quantile\":\"null_with_finite_log_quantile\"}}",vector(&locations)?))
}
fn nhpp(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 1, 1)?;
    if data[0].len() != 1
        || data.len() < 2
        || args[0] <= 0.
        || data.iter().any(|r| r[0] <= 0. || r[0] > args[0])
    {
        return Err(
            "power NHPP needs>=2 positive exact recurrence times <=positive observation_end".into(),
        );
    }
    let end = args[0];
    let n = data.len() as f64;
    let sum = data.iter().map(|r| log_ratio(r[0], end)).sum::<f64>();
    if sum >= 0. {
        return Err("NHPP shape is unidentified when all event times equal observation end".into());
    }
    let shape = finite(-n / sum, "NHPP shape")?;
    let logeta = finite(end.ln() - n.ln() / shape, "NHPP log eta")?;
    let ll = finite(
        n * n.ln() + n * shape.ln() - n * end.ln() + (shape - 1.) * sum - n,
        "NHPP point-process log likelihood",
    )?;
    let se = shape / n.sqrt();
    let hse = 1. / n.sqrt();
    let logeta_se = ((1. + n.ln().powi(2)) / n).sqrt() / shape;
    let scale = optional_exp(logeta).map_or_else(|| "null".into(), |v| v.to_string());
    let interval = vector(&[
        shape * (-1.959_963_984_540_054 * hse).exp(),
        shape * (1.959_963_984_540_054 * hse).exp(),
    ])?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"time_truncated_single_system_powerlaw_nhpp_mle\",\"observation_end\":{end},\"event_count\":{n},\"shape\":{shape},\"log_eta\":{logeta},\"eta\":{scale},\"shape_standard_error\":{se},\"log_shape_standard_error\":{hse},\"log_eta_standard_error\":{logeta_se},\"shape_confidence_interval\":{interval},\"log_likelihood\":{ll},\"aic\":{},\"bic\":{},\"expected_events_at_end\":{n},\"intensity_at_end\":{},\"growth_improving\":{},\"predict_parameters\":{},\"converged\":true,\"scope\":\"single_system_exact_events_time_truncation_not_phase_reset_or_interval_counts\"}}",4.-2.*ll,2.*n.ln()-2.*ll,finite((n/end)*shape,"NHPP endpoint intensity")?,shape<1.,vector(&[shape,logeta])?))
}
fn nhpp_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 2, 2)?;
    if data[0].len() != 1 || args[0] <= 0. || data.iter().any(|r| r[0] < 0.) {
        return Err(
            "NHPP prediction requires nonnegative scalar times and[positive_shape,log_eta]".into(),
        );
    }
    let shape = args[0];
    let logeta = args[1];
    let mut cumulative = Vec::new();
    let mut intensity = Vec::new();
    let mut logc = Vec::new();
    for r in data {
        if r[0] == 0. {
            cumulative.push(Some(0.));
            logc.push(None);
            intensity.push(if shape > 1. {
                Some(0.)
            } else if shape == 1. {
                optional_exp(-logeta)
            } else {
                None
            });
        } else {
            let lc = finite(shape * (r[0].ln() - logeta), "NHPP log cumulative count")?;
            logc.push(Some(lc));
            let c = lc.exp();
            cumulative.push(c.is_finite().then_some(c));
            let li = finite(shape.ln() - r[0].ln() + lc, "NHPP log intensity")?;
            let rate = li.exp();
            intensity.push(rate.is_finite().then_some(rate));
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"powerlaw_nhpp_prediction\",\"expected_cumulative_events\":{},\"log_expected_cumulative_events\":{},\"intensities\":{},\"unrepresentable_or_infinite_values\":\"null_log_cumulative_retained\"}}",optional_vector(&cumulative),optional_vector(&logc),optional_vector(&intensity)))
}
fn blocks(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 2)?;
    let mode = integer(args.first().copied().unwrap_or(0.), 0, 2, "block structure")?;
    if data.iter().flatten().any(|&p| !(0. ..=1.).contains(&p)) {
        return Err("component reliabilities must be in[0,1]".into());
    }
    let n = data[0].len();
    let k = if mode == 2 {
        arguments(args, 2, 2)?;
        integer(args[1], 1, n, "k-of-n threshold")?
    } else {
        if args.len() > 1 {
            return Err("k is only accepted for k-of-n structure2".into());
        }
        if mode == 0 {
            n
        } else {
            1
        }
    };
    if mode == 2
        && data
            .len()
            .checked_mul(n)
            .and_then(|v| v.checked_mul(n))
            .is_none_or(|v| v > 50_000_000)
    {
        return Err("k-of-n block work exceeds50000000".into());
    }
    let mut reliability = Vec::new();
    let mut working = Vec::new();
    for row in data {
        let value = if mode == 0 {
            if row.contains(&0.) {
                0.
            } else {
                row.iter().map(|p| p.ln()).sum::<f64>().exp()
            }
        } else if mode == 1 {
            -row.iter().map(|p| (-p).ln_1p()).sum::<f64>().exp_m1()
        } else {
            let mut mass = vec![0.; n + 1];
            mass[0] = 1.;
            for (i, &p) in row.iter().enumerate() {
                for j in (0..=i + 1).rev() {
                    mass[j] = mass[j] * (1. - p) + if j > 0 { mass[j - 1] * p } else { 0. };
                }
            }
            mass[k..].iter().sum::<f64>()
        };
        reliability.push(finite(value, "system reliability")?.clamp(0., 1.));
        working.push(row.iter().sum::<f64>());
    }
    let name = ["series", "parallel", "k_of_n"][mode];
    Ok(format!("{{\"version\":1,\"algorithm\":\"independent_component_reliability_blocks\",\"structure\":\"{name}\",\"k\":{k},\"n_components\":{n},\"reliabilities\":{},\"expected_working_components\":{},\"assumption\":\"mutually_independent_component_states\",\"scope\":\"single_block_nodes_no_dependency_or_repair_simulation\"}}",vector(&reliability)?,vector(&working)?))
}
fn log_gamma(z: f64) -> f64 {
    const C: [f64; 8] = [
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if z < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * z).sin().ln()
            - log_gamma(1. - z);
    }
    let z = z - 1.;
    let mut a = 0.999_999_999_999_809_9;
    for (i, c) in C.iter().enumerate() {
        a += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    LOG_SQRT_TAU + (z + 0.5) * t.ln() - t + a.ln()
}
fn gamma_q(a: f64, x: f64) -> Result<f64, String> {
    if x == 0. {
        return Ok(1.);
    }
    if x.is_infinite() {
        return Ok(0.);
    }
    let front = (a * x.ln() - x - log_gamma(a)).exp();
    if x < a + 1. {
        let mut term = 1. / a;
        let mut sum = term;
        for i in 1..=10000 {
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() < sum.abs() * 2e-14 {
                return Ok((1. - sum * front).clamp(0., 1.));
            }
        }
    } else {
        let nz = |v: f64| {
            if v.abs() < 1e-300 {
                1e-300_f64.copysign(v)
            } else {
                v
            }
        };
        let mut b = x + 1. - a;
        let mut c = 1e300;
        let mut d = 1. / nz(b);
        let mut h = d;
        for i in 1..=10000 {
            let i = i as f64;
            let an = -i * (i - a);
            b += 2.;
            d = 1. / nz(an * d + b);
            c = nz(b + an / c);
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-14 {
                return finite(front * h, "normal tail probability");
            }
        }
    }
    Err("normal gamma tail did not converge".into())
}
// Stable log survival and inverse Mills ratio, including its derivative.
// Laplace's continued fraction avoids cancellation in r*(r-z) for large z.
fn normal_tail(z: f64) -> Result<(f64, f64, f64), String> {
    if z >= 8. {
        let mut r = z;
        let mut derivative = 1.;
        for k in (1..=100).rev() {
            let old = r;
            let old_d = derivative;
            r = z + k as f64 / old;
            derivative = 1. - k as f64 * old_d / (old * old);
        }
        let logsf = finite(-0.5 * z * z - LOG_SQRT_TAU - r.ln(), "normal log survival")?;
        return Ok((logsf, r, derivative));
    }
    let small = 0.5 * gamma_q(0.5, z * z / 2.)?;
    let logsf = if z >= 0. {
        small.ln()
    } else {
        (-small).ln_1p()
    };
    let r = (-0.5 * z * z - LOG_SQRT_TAU - logsf).exp();
    Ok((logsf, r, r * (r - z)))
}
fn normal_two_sided(z: f64) -> Result<f64, String> {
    Ok((2. * normal_tail(z.abs())?.0.exp()).min(1.))
}
fn normal_quantile(p: f64) -> Result<f64, String> {
    if p == 0.5 {
        return Ok(0.);
    }
    let target = if p < 0.5 { p } else { 1. - p };
    let logarithm = target.ln();
    let mut low = 0.;
    let mut high = 1.;
    while normal_tail(high)?.0 > logarithm {
        high *= 2.;
        if high > 128. {
            return Err("normal quantile tail outside supported range".into());
        }
    }
    for _ in 0..256 {
        let mid = (low + high) / 2.;
        if normal_tail(mid)?.0 > logarithm {
            low = mid;
        } else {
            high = mid;
        }
        if high - low < 2e-14 * high.max(1.) {
            break;
        }
    }
    let q = (low + high) / 2.;
    Ok(if p < 0.5 { -q } else { q })
}
