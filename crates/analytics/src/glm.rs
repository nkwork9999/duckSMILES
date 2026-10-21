//! Ordinary linear inference, two canonical GLMs, elastic net, and feature expansion.
//!
//! Fit rows are `[features..., response]`. All fit preprocessing uses only the
//! supplied training rows. GLM penalties minimize summed negative log likelihood
//! plus `lambda / 2 * ||standardized_slopes||^2`; the intercept is unpenalized.
//! Elastic net minimizes `SSE/(2*n) + lambda * (alpha*L1 + (1-alpha)*L2/2)`
//! on population-standardized features. OLS confidence intervals assume independent
//! homoscedastic Gaussian errors; neither weights nor categorical coding is inferred.
//!
//! Prediction accepts `[link, intercept, slopes...]` (0 identity, 1 logit, 2 log),
//! or the returned `stable_predict_parameters`: `[link, response_origin,
//! response_mean_offset, feature_origins..., feature_mean_offsets..., slopes...]`.
//! The latter retains small variation at large absolute feature/response offsets.

const MAX_CELLS: usize = 1_000_000;
const MAX_FEATURES: usize = 256;
const MAX_WORK: usize = 200_000_000;
const MAX_INTEGER: f64 = 9_007_199_254_740_991.;

/// Executes a supported operation and returns finite, version-1 JSON.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = validate_matrix(data)?;
    if args.iter().any(|x| !x.is_finite()) {
        return Err("parameters must be finite".into());
    }
    match operation {
        "linear_model" => linear_model(data, args),
        "logistic" | "poisson" => generalized_linear_model(operation, data, args),
        "elastic_net" => elastic_net(data, args),
        "glm_predict" => predictions(data, width, args),
        "polynomial_features" => polynomial_features(data, width, args),
        _ => Err(format!("unknown GLM operation: {operation}")),
    }
}

fn validate_matrix(data: &[Vec<f64>]) -> Result<usize, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0 || data.len() > MAX_CELLS / width {
        return Err("data must be nonempty with at most 1000000 cells".into());
    }
    if data
        .iter()
        .any(|row| row.len() != width || row.iter().any(|x| !x.is_finite()))
    {
        return Err("data must be rectangular and finite".into());
    }
    Ok(width)
}

fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}

fn mean(values: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut result = 0.;
    for (i, x) in values.enumerate() {
        let difference = finite(x - result, "mean difference")?;
        result = finite(result + difference / (i + 1) as f64, "mean")?;
    }
    Ok(result)
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}

fn norm(x: impl Iterator<Item = f64>) -> f64 {
    x.fold(0_f64, f64::hypot)
}

fn numbers(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|x| nullable(*x))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn nullable(x: f64) -> String {
    if x.is_finite() {
        x.to_string()
    } else {
        "null".into()
    }
}

fn matrix(values: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|row| numbers(row))
            .collect::<Vec<_>>()
            .join(",")
    )
}

struct Prepared {
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    response_origin: f64,
    response_offset: f64,
}

impl Prepared {
    fn new(data: &[Vec<f64>], center_response: bool) -> Result<Self, String> {
        let p = data[0].len().saturating_sub(1);
        if data.len() < 2 || p == 0 || p > MAX_FEATURES {
            return Err("fit requires at least two rows and 1 to 256 features".into());
        }
        let origins = data[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|row| row[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = finite(
                    norm(data.iter().map(|row| (row[j] - origins[j]) - offsets[j]))
                        / (data.len() as f64).sqrt(),
                    "feature scale",
                )?;
                Ok(if s == 0. { 1. } else { s })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let x = data
            .iter()
            .map(|row| {
                (0..p)
                    .map(|j| {
                        finite(
                            ((row[j] - origins[j]) - offsets[j]) / scales[j],
                            "scaled feature",
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        let response_origin = if center_response { data[0][p] } else { 0. };
        let response_offset = if center_response {
            mean(data.iter().map(|row| row[p] - response_origin))?
        } else {
            0.
        };
        let y = data
            .iter()
            .map(|row| {
                finite(
                    (row[p] - response_origin) - response_offset,
                    "centered response",
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            x,
            y,
            origins,
            offsets,
            scales,
            response_origin,
            response_offset,
        })
    }

    fn original_coefficients(&self, scaled: &[f64]) -> Result<Vec<f64>, String> {
        scaled
            .iter()
            .zip(&self.scales)
            .map(|(b, s)| finite(b / s, "coefficient"))
            .collect()
    }

    fn intercept(&self, coefficients: &[f64], center: f64) -> Result<f64, String> {
        finite(
            (self.response_origin - dot(&self.origins, coefficients))
                + (self.response_offset + center - dot(&self.offsets, coefficients)),
            "intercept",
        )
    }

    fn stable_parameters(
        &self,
        link: f64,
        center: f64,
        coefficients: &[f64],
    ) -> Result<Vec<f64>, String> {
        let mut values = vec![
            link,
            self.response_origin,
            finite(self.response_offset + center, "prediction center")?,
        ];
        values.extend_from_slice(&self.origins);
        values.extend_from_slice(&self.offsets);
        values.extend_from_slice(coefficients);
        Ok(values)
    }

    fn means(&self) -> Result<Vec<f64>, String> {
        self.origins
            .iter()
            .zip(&self.offsets)
            .map(|(a, b)| finite(a + b, "feature mean"))
            .collect()
    }
}

struct Qr {
    coefficients: Vec<f64>,
    r: Vec<Vec<f64>>,
}

// Householder QR keeps condition numbers unsquared, including weighted GLM steps.
#[allow(clippy::needless_range_loop)]
fn qr_solve(mut x: Vec<Vec<f64>>, mut y: Vec<f64>) -> Result<Qr, String> {
    let p = x[0].len();
    if x.len() < p {
        return Err("rank-deficient design".into());
    }
    let tolerance =
        f64::EPSILON * 64. * (x.len().max(p) as f64).sqrt() * norm(x.iter().flatten().copied());
    for j in 0..p {
        let column_norm = finite(norm(x[j..].iter().map(|row| row[j])), "QR norm")?;
        if column_norm <= tolerance {
            return Err("rank-deficient or numerically singular design".into());
        }
        let diagonal = -column_norm.copysign(x[j][j]);
        let mut reflector = x[j..]
            .iter()
            .map(|row| row[j] / column_norm)
            .collect::<Vec<_>>();
        reflector[0] -= diagonal / column_norm;
        let length = norm(reflector.iter().copied());
        for v in &mut reflector {
            *v /= length;
        }
        for k in j..p {
            let projection = 2.
                * reflector
                    .iter()
                    .enumerate()
                    .map(|(i, v)| v * x[j + i][k])
                    .sum::<f64>();
            for (i, v) in reflector.iter().enumerate() {
                x[j + i][k] -= projection * v;
            }
        }
        let projection = 2. * dot(&reflector, &y[j..]);
        for (i, v) in reflector.iter().enumerate() {
            y[j + i] -= projection * v;
        }
        x[j][j] = diagonal;
    }
    let mut coefficients = vec![0.; p];
    for j in (0..p).rev() {
        coefficients[j] = finite(
            (y[j] - dot(&x[j][j + 1..], &coefficients[j + 1..])) / x[j][j],
            "QR coefficient",
        )?;
    }
    Ok(Qr {
        coefficients,
        r: x[..p].to_vec(),
    })
}

// Triangular inverse entries depend on multiple rows and columns at each step.
#[allow(clippy::needless_range_loop)]
fn linear_model(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 1 {
        return Err("linear_model requires [confidence=0.95]".into());
    }
    let confidence = args.first().copied().unwrap_or(0.95);
    if !(0. < confidence && confidence < 1.) {
        return Err("confidence must be between 0 and 1".into());
    }
    let f = Prepared::new(data, true)?;
    let n = data.len();
    let p = f.scales.len();
    if n <= p + 1 {
        return Err("OLS inference requires positive residual degrees of freedom".into());
    }
    if n.saturating_mul(p).saturating_mul(p) > MAX_WORK {
        return Err("OLS exceeds numerical work limit".into());
    }
    let qr = qr_solve(f.x.clone(), f.y.clone())?;
    let coefficients = f.original_coefficients(&qr.coefficients)?;
    let intercept = f.intercept(&coefficients, 0.)?;
    let centered_fits =
        f.x.iter()
            .map(|row| finite(dot(row, &qr.coefficients), "fitted value"))
            .collect::<Result<Vec<_>, _>>()?;
    let fitted = centered_fits
        .iter()
        .map(|v| finite(f.response_origin + (f.response_offset + v), "fitted value"))
        .collect::<Result<Vec<_>, _>>()?;
    let residuals =
        f.y.iter()
            .zip(&centered_fits)
            .map(|(y, v)| finite(y - v, "residual"))
            .collect::<Result<Vec<_>, _>>()?;
    let sse = finite(norm(residuals.iter().copied()).powi(2), "SSE")?;
    let sst = finite(norm(f.y.iter().copied()).powi(2), "SST")?;
    let df = n - p - 1;
    let sigma_squared = sse / df as f64;
    let r2 = if sst > 0. { 1. - sse / sst } else { f64::NAN };
    let adjusted = 1. - (1. - r2) * (n - 1) as f64 / df as f64;
    let mut inverse = vec![vec![0.; p]; p];
    for k in 0..p {
        for j in (0..p).rev() {
            let rhs = if j == k { 1. } else { 0. };
            inverse[j][k] = finite(
                (rhs - dot(
                    &qr.r[j][j + 1..],
                    &(j + 1..p).map(|i| inverse[i][k]).collect::<Vec<_>>(),
                )) / qr.r[j][j],
                "coefficient covariance",
            )?;
        }
    }
    let standard_errors = (0..p)
        .map(|j| {
            finite(
                sigma_squared.sqrt() * norm(inverse[j].iter().copied()) / f.scales[j],
                "coefficient standard error",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    // Intercept uncertainty uses the original-coordinate feature mean and R^-T.
    let means = f.means()?;
    let mut intercept_weights = vec![0.; p];
    for j in 0..p {
        let rhs = finite(means[j] / f.scales[j], "intercept covariance")?;
        let projected = (0..j)
            .map(|k| qr.r[k][j] * intercept_weights[k])
            .sum::<f64>();
        intercept_weights[j] = finite((rhs - projected) / qr.r[j][j], "intercept covariance")?;
    }
    let intercept_se = finite(
        sigma_squared.sqrt()
            * (1. / n as f64 + norm(intercept_weights.iter().copied()).powi(2)).sqrt(),
        "intercept standard error",
    )?;
    let quantile = student_quantile(confidence, df as f64)?;
    let t_stats = coefficients
        .iter()
        .zip(&standard_errors)
        .map(|(b, se)| b / se)
        .collect::<Vec<_>>();
    let p_values = t_stats
        .iter()
        .map(|t| student_two_sided(*t, df as f64))
        .collect::<Result<Vec<_>, _>>()?;
    let intervals = coefficients
        .iter()
        .zip(&standard_errors)
        .map(|(b, se)| {
            Ok(vec![
                finite(b - quantile * se, "confidence interval")?,
                finite(b + quantile * se, "confidence interval")?,
            ])
        })
        .collect::<Result<Vec<_>, String>>()?;
    let intercept_t = intercept / intercept_se;
    let intercept_p = student_two_sided(intercept_t, df as f64)?;
    let intercept_ci = vec![
        finite(intercept - quantile * intercept_se, "intercept interval")?,
        finite(intercept + quantile * intercept_se, "intercept interval")?,
    ];
    let overall_f = if sst > 0. {
        ((sst - sse).max(0.) / p as f64) / sigma_squared
    } else {
        f64::NAN
    };
    let overall_p = if overall_f.is_nan() {
        f64::NAN
    } else if overall_f.is_infinite() {
        0.
    } else {
        regularized_beta(
            df as f64 / (df as f64 + p as f64 * overall_f),
            df as f64 / 2.,
            p as f64 / 2.,
        )?
    };
    let stable = f.stable_parameters(0., 0., &coefficients)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"linear_model\",\"n_samples\":{n},\"n_features\":{p},\"intercept\":{intercept},\"coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"sse\":{sse},\"residual_df\":{df},\"sigma_squared\":{sigma_squared},\"r2\":{},\"adjusted_r2\":{},\"confidence\":{confidence},\"intercept_standard_error\":{intercept_se},\"coefficient_standard_errors\":{},\"intercept_t_statistic\":{},\"coefficient_t_statistics\":{},\"intercept_p_value\":{},\"coefficient_p_values\":{},\"intercept_confidence_interval\":{},\"coefficient_confidence_intervals\":{},\"overall_f\":{},\"overall_p_value\":{},\"stable_predict_parameters\":{},\"feature_means\":{},\"feature_scales\":{},\"standardization\":\"training_population_standard_deviation\"}}",numbers(&coefficients),numbers(&fitted),numbers(&residuals),nullable(r2),nullable(adjusted),numbers(&standard_errors),nullable(intercept_t),numbers(&t_stats),nullable(intercept_p),numbers(&p_values),numbers(&intercept_ci),matrix(&intervals),nullable(overall_f),nullable(overall_p),numbers(&stable),numbers(&means),numbers(&f.scales)))
}

fn iteration_args(args: &[f64], default_iterations: usize) -> Result<(f64, usize, f64), String> {
    if args.len() > 3 {
        return Err("GLM requires [lambda=0, max_iterations=100, tolerance=1e-8]".into());
    }
    let lambda = args.first().copied().unwrap_or(0.);
    let iterations = args.get(1).copied().unwrap_or(default_iterations as f64);
    let tolerance = args.get(2).copied().unwrap_or(1e-8);
    if lambda < 0. {
        return Err("lambda must be nonnegative".into());
    }
    if !(1. ..=10_000.).contains(&iterations) || iterations.fract() != 0. {
        return Err("max_iterations must be an integer between 1 and 10000".into());
    }
    if !(0. < tolerance && tolerance < 1.) {
        return Err("tolerance must be between 0 and 1".into());
    }
    Ok((lambda, iterations as usize, tolerance))
}

fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}

fn softplus(x: f64) -> f64 {
    x.max(0.) + (-x.abs()).exp().ln_1p()
}

fn glm_objective(
    method: &str,
    z: &[Vec<f64>],
    y: &[f64],
    b: &[f64],
    lambda: f64,
) -> Result<f64, String> {
    let mut result = 0.;
    for (row, target) in z.iter().zip(y) {
        let eta = finite(dot(row, b), "linear predictor")?;
        let loss = if method == "logistic" {
            if *target == 1. {
                softplus(-eta)
            } else {
                softplus(eta)
            }
        } else {
            finite(eta.exp() - target * eta, "Poisson likelihood")?
        };
        result = finite(result + loss, "negative log likelihood")?;
    }
    finite(
        result + 0.5 * lambda * dot(&b[1..], &b[1..]),
        "penalized objective",
    )
}

fn poisson_deviance(y: &[f64], mu: &[f64]) -> Result<f64, String> {
    let mut sum = 0.;
    for (target, fit) in y.iter().zip(mu) {
        if *fit <= 0. {
            return Err("Poisson mean underflow".into());
        }
        let contribution = if *target == 0. {
            *fit
        } else {
            target * (target / fit).ln() - (target - fit)
        };
        sum = finite(sum + contribution.max(0.), "Poisson deviance")?;
    }
    finite(2. * sum, "Poisson deviance")
}

fn generalized_linear_model(
    method: &str,
    data: &[Vec<f64>],
    args: &[f64],
) -> Result<String, String> {
    let (lambda, max_iterations, tolerance) = iteration_args(args, 100)?;
    let f = Prepared::new(data, false)?;
    let n = data.len();
    let p = f.scales.len();
    if n.saturating_mul(p + 1)
        .saturating_mul(p + 1)
        .saturating_mul(max_iterations)
        > MAX_WORK
    {
        return Err("GLM exceeds numerical work limit".into());
    }
    for y in &f.y {
        if method == "logistic" && *y != 0. && *y != 1. {
            return Err("logistic response must be 0 or 1".into());
        }
        if method == "poisson" && (*y < 0. || *y > MAX_INTEGER || y.fract() != 0.) {
            return Err("Poisson response must be an exact nonnegative integer".into());
        }
    }
    let ymean = mean(f.y.iter().copied())?;
    if ymean == 0. || (method == "logistic" && ymean == 1.) {
        return Err("response has no finite intercept maximum likelihood estimate".into());
    }
    let initial_intercept = if method == "logistic" {
        (ymean / (1. - ymean)).ln()
    } else {
        ymean.ln()
    };
    let z =
        f.x.iter()
            .map(|row| {
                let mut v = vec![1.];
                v.extend_from_slice(row);
                v
            })
            .collect::<Vec<_>>();
    // An unpenalized singular model has no uniquely identifiable coefficients.
    if lambda == 0. {
        qr_solve(z.clone(), vec![0.; n])?;
    }
    let mut b = vec![0.; p + 1];
    b[0] = initial_intercept;
    let mut objective = glm_objective(method, &z, &f.y, &b, lambda)?;
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..max_iterations {
        let mut weighted = Vec::with_capacity(n + p);
        let mut rhs = Vec::with_capacity(n + p);
        for (row, y) in z.iter().zip(&f.y) {
            let eta = finite(dot(row, &b), "linear predictor")?;
            let (mu, weight) = if method == "logistic" {
                let e = (-eta.abs()).exp();
                (sigmoid(eta), e / (1. + e).powi(2))
            } else {
                let e = finite(eta.exp(), "Poisson mean")?;
                if e == 0. {
                    return Err("Poisson mean underflow".into());
                }
                (e, e)
            };
            if weight <= 0. {
                return Err("GLM information is singular, possibly separated data".into());
            }
            let sqrt_weight = weight.sqrt();
            weighted.push(
                row.iter()
                    .map(|x| finite(x * sqrt_weight, "weighted design"))
                    .collect::<Result<Vec<_>, _>>()?,
            );
            rhs.push(finite((y - mu) / sqrt_weight, "IRLS step")?);
        }
        if lambda > 0. {
            let penalty = lambda.sqrt();
            for j in 1..=p {
                let mut row = vec![0.; p + 1];
                row[j] = penalty;
                weighted.push(row);
                rhs.push(finite(-penalty * b[j], "penalty step")?);
            }
        }
        let step = qr_solve(weighted, rhs)?.coefficients;
        let step_size = norm(step.iter().copied());
        let threshold = tolerance * (1. + norm(b.iter().copied()));
        let mut accepted = None;
        let mut damping = 1.;
        for _ in 0..50 {
            let proposal = b
                .iter()
                .zip(&step)
                .map(|(v, d)| v + damping * d)
                .collect::<Vec<_>>();
            if let Ok(next) = glm_objective(method, &z, &f.y, &proposal, lambda) {
                if next <= objective + f64::EPSILON * 16. * (1. + objective.abs()) {
                    accepted = Some((proposal, next));
                    break;
                }
            }
            damping *= 0.5;
        }
        iterations = iteration + 1;
        let Some((proposal, next)) = accepted else {
            break;
        };
        b = proposal;
        objective = next;
        if step_size <= threshold {
            converged = true;
            break;
        }
    }
    let eta = z
        .iter()
        .map(|row| finite(dot(row, &b), "linear predictor"))
        .collect::<Result<Vec<_>, _>>()?;
    let fitted = eta
        .iter()
        .map(|x| {
            if method == "logistic" {
                Ok(sigmoid(*x))
            } else {
                finite(x.exp(), "Poisson mean")
            }
        })
        .collect::<Result<Vec<_>, _>>()?;
    if method == "logistic"
        && lambda == 0.
        && f.y.iter().zip(&fitted).all(|(y, m)| (y - m).abs() < 1e-8)
    {
        return Err(
            "unpenalized logistic data are completely separated; use positive lambda".into(),
        );
    }
    let deviance = if method == "logistic" {
        finite(2. * glm_objective(method, &z, &f.y, &b, 0.)?, "deviance")?
    } else {
        poisson_deviance(&f.y, &fitted)?
    };
    let null_deviance = if method == "logistic" {
        finite(
            -2. * (n as f64) * (ymean * ymean.ln() + (1. - ymean) * (1. - ymean).ln()),
            "null deviance",
        )?
    } else {
        poisson_deviance(&f.y, &vec![ymean; n])?
    };
    let coefficients = f.original_coefficients(&b[1..])?;
    let intercept = f.intercept(&coefficients, b[0])?;
    let link = if method == "logistic" { 1. } else { 2. };
    let stable = f.stable_parameters(link, b[0], &coefficients)?;
    let mut predict_parameters = vec![link, intercept];
    predict_parameters.extend_from_slice(&coefficients);
    let mut gradient = vec![0.; p + 1];
    for ((row, y), mu) in z.iter().zip(&f.y).zip(&fitted) {
        for j in 0..=p {
            gradient[j] += row[j] * (mu - y);
        }
    }
    for j in 1..=p {
        gradient[j] += lambda * b[j];
    }
    let max_gradient = finite(
        gradient.iter().fold(0_f64, |a, v| a.max(v.abs())),
        "gradient",
    )?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"{method}\",\"n_samples\":{n},\"n_features\":{p},\"lambda\":{lambda},\"intercept\":{intercept},\"coefficients\":{},\"fitted_values\":{},\"linear_predictors\":{},\"converged\":{converged},\"iterations\":{iterations},\"max_iterations\":{max_iterations},\"tolerance\":{tolerance},\"gradient_max_abs\":{max_gradient},\"objective\":{objective},\"deviance\":{deviance},\"null_deviance\":{null_deviance},\"predict_parameters\":{},\"stable_predict_parameters\":{},\"feature_means\":{},\"feature_scales\":{},\"standardization\":\"training_population_standard_deviation\",\"penalty\":\"sum_negative_log_likelihood_plus_lambda_half_squared_standardized_slopes\"}}",numbers(&coefficients),numbers(&fitted),numbers(&eta),numbers(&predict_parameters),numbers(&stable),numbers(&f.means()?),numbers(&f.scales)))
}

fn soft_threshold(x: f64, threshold: f64) -> f64 {
    x.signum() * (x.abs() - threshold).max(0.)
}

fn elastic_net(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 2 || args.len() > 4 {
        return Err(
            "elastic_net requires [lambda, alpha, max_iterations=10000, tolerance=1e-8]".into(),
        );
    }
    let lambda = args[0];
    let alpha = args[1];
    if lambda < 0. || !(0. ..=1.).contains(&alpha) {
        return Err("lambda must be nonnegative and alpha between 0 and 1".into());
    }
    let mut iteration_parameters = vec![lambda];
    iteration_parameters.extend_from_slice(&args[2..]);
    let (_, max_iterations, tolerance) = iteration_args(&iteration_parameters, 10_000)?;
    let f = Prepared::new(data, true)?;
    let n = data.len();
    let p = f.scales.len();
    if n.saturating_mul(p).saturating_mul(max_iterations) > MAX_WORK {
        return Err("elastic net exceeds numerical work limit".into());
    }
    let mut coefficients_scaled = vec![0.; p];
    let mut residuals = f.y.clone();
    let variance = (0..p)
        .map(|j| {
            finite(
                f.x.iter().map(|row| row[j] * row[j]).sum::<f64>() / n as f64,
                "feature variance",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut converged = false;
    let mut iterations = 0;
    let mut kkt = 0.;
    for iteration in 0..max_iterations {
        for j in 0..p {
            let old = coefficients_scaled[j];
            let partial = finite(
                f.x.iter()
                    .zip(&residuals)
                    .map(|(row, r)| row[j] * r)
                    .sum::<f64>()
                    / n as f64
                    + variance[j] * old,
                "coordinate gradient",
            )?;
            let denominator = variance[j] + lambda * (1. - alpha);
            let next = if denominator > 0. {
                finite(
                    soft_threshold(partial, lambda * alpha) / denominator,
                    "coordinate coefficient",
                )?
            } else {
                0.
            };
            let change = next - old;
            for (row, residual) in f.x.iter().zip(&mut residuals) {
                *residual = finite(*residual - row[j] * change, "coordinate residual")?;
            }
            coefficients_scaled[j] = next;
        }
        kkt = 0_f64;
        for j in 0..p {
            let gradient = finite(
                -f.x.iter()
                    .zip(&residuals)
                    .map(|(row, r)| row[j] * r)
                    .sum::<f64>()
                    / n as f64
                    + lambda * (1. - alpha) * coefficients_scaled[j],
                "KKT gradient",
            )?;
            let violation = if coefficients_scaled[j] == 0. {
                (gradient.abs() - lambda * alpha).max(0.)
            } else {
                (gradient + lambda * alpha * coefficients_scaled[j].signum()).abs()
            };
            kkt = kkt.max(violation);
        }
        iterations = iteration + 1;
        if kkt <= tolerance * (1. + norm(f.y.iter().copied()) / (n as f64).sqrt()) {
            converged = true;
            break;
        }
    }
    let sse = finite(norm(residuals.iter().copied()).powi(2), "SSE")?;
    let objective = finite(
        sse / (2. * n as f64)
            + lambda
                * (alpha * coefficients_scaled.iter().map(|b| b.abs()).sum::<f64>()
                    + 0.5 * (1. - alpha) * dot(&coefficients_scaled, &coefficients_scaled)),
        "elastic-net objective",
    )?;
    let coefficients = f.original_coefficients(&coefficients_scaled)?;
    let intercept = f.intercept(&coefficients, 0.)?;
    let fitted =
        f.x.iter()
            .map(|row| {
                finite(
                    f.response_origin + (f.response_offset + dot(row, &coefficients_scaled)),
                    "fitted value",
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
    let stable = f.stable_parameters(0., 0., &coefficients)?;
    let nonzero = coefficients.iter().filter(|x| **x != 0.).count();
    Ok(format!("{{\"version\":1,\"algorithm\":\"elastic_net\",\"n_samples\":{n},\"n_features\":{p},\"lambda\":{lambda},\"alpha\":{alpha},\"intercept\":{intercept},\"coefficients\":{},\"standardized_coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"sse\":{sse},\"objective\":{objective},\"converged\":{converged},\"iterations\":{iterations},\"max_iterations\":{max_iterations},\"tolerance\":{tolerance},\"kkt_max_violation\":{kkt},\"nonzero_coefficients\":{nonzero},\"stable_predict_parameters\":{},\"feature_means\":{},\"feature_scales\":{},\"standardization\":\"training_population_standard_deviation\",\"penalty\":\"sse_over_2n_plus_lambda_alpha_l1_plus_lambda_half_one_minus_alpha_l2\"}}",numbers(&coefficients),numbers(&coefficients_scaled),numbers(&fitted),numbers(&residuals),numbers(&stable),numbers(&f.means()?),numbers(&f.scales)))
}

fn predictions(data: &[Vec<f64>], width: usize, args: &[f64]) -> Result<String, String> {
    let stable = args.len() == 3 * width + 3;
    if args.len() != width + 2 && !stable {
        return Err(
            "prediction requires [link,intercept,slopes...] or stable_predict_parameters".into(),
        );
    }
    let link = args[0];
    if link != 0. && link != 1. && link != 2. {
        return Err("link must be 0 identity, 1 logit, or 2 log".into());
    }
    let mut result = Vec::with_capacity(data.len());
    for row in data {
        let eta = if stable {
            let adjustment = (0..width)
                .map(|j| ((row[j] - args[3 + j]) - args[3 + width + j]) * args[3 + 2 * width + j])
                .sum::<f64>();
            finite(args[1] + (args[2] + adjustment), "prediction")?
        } else {
            finite(args[1] + dot(row, &args[2..]), "prediction")?
        };
        result.push(if link == 1. {
            sigmoid(eta)
        } else if link == 2. {
            finite(eta.exp(), "log-link prediction")?
        } else {
            eta
        });
    }
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"glm_predict\",\"link\":{link},\"predictions\":{}}}",
        numbers(&result)
    ))
}

fn polynomial_features(data: &[Vec<f64>], width: usize, args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 2 {
        return Err("polynomial_features requires [degree, include_bias=0]".into());
    }
    let degree = args[0];
    let bias = args.get(1).copied().unwrap_or(0.);
    if !(1. ..=3.).contains(&degree) || degree.fract() != 0. || (bias != 0. && bias != 1.) {
        return Err("degree must be 1, 2, or 3 and include_bias 0 or 1".into());
    }
    if width > MAX_FEATURES {
        return Err("polynomial expansion accepts at most 256 input features".into());
    }
    let mut count = if bias == 1. { 1_usize } else { 0 };
    let mut degree_count = 1_usize;
    for d in 1..=degree as usize {
        degree_count = degree_count
            .checked_mul(width + d - 1)
            .ok_or("polynomial feature count overflow")?
            / d;
        count = count
            .checked_add(degree_count)
            .ok_or("polynomial feature count overflow")?;
        if count > 4096 || data.len() > MAX_CELLS / count {
            return Err("polynomial expansion exceeds 4096 terms or 1000000 output cells".into());
        }
    }
    let mut powers = Vec::with_capacity(count);
    if bias == 1. {
        powers.push(vec![0.; width]);
    }
    fn combinations(
        width: usize,
        start: usize,
        remaining: usize,
        power: &mut Vec<f64>,
        powers: &mut Vec<Vec<f64>>,
    ) {
        if remaining == 0 {
            powers.push(power.clone());
            return;
        }
        for j in start..width {
            power[j] += 1.;
            combinations(width, j, remaining - 1, power, powers);
            power[j] -= 1.;
        }
    }
    for d in 1..=degree as usize {
        combinations(width, 0, d, &mut vec![0.; width], &mut powers);
    }
    let transformed = data
        .iter()
        .map(|row| {
            powers
                .iter()
                .map(|power| {
                    let mut value = 1.;
                    for (x, p) in row.iter().zip(power) {
                        value = finite(value * x.powi(*p as i32), "polynomial feature")?;
                    }
                    Ok(value)
                })
                .collect::<Result<Vec<_>, String>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"polynomial_features\",\"degree\":{degree},\"include_bias\":{},\"n_input_features\":{width},\"n_output_features\":{count},\"powers\":{},\"transformed_rows\":{}}}",bias==1.,matrix(&powers),matrix(&transformed)))
}

// Student-t and F tails use the regularized incomplete beta function.
fn log_gamma(z: f64) -> f64 {
    let c = [
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
    let mut x = 0.9999999999998099;
    for (i, c) in c.iter().enumerate() {
        x += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * (2. * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

fn beta_fraction(a: f64, b: f64, x: f64) -> Result<f64, String> {
    let minimum = 1e-300;
    let mut c = 1.;
    let mut d = 1. - (a + b) * x / (a + 1.);
    if d.abs() < minimum {
        d = minimum;
    }
    d = 1. / d;
    let mut h = d;
    for m in 1..=10_000 {
        let m = m as f64;
        let m2 = 2. * m;
        let aa = m * (b - m) * x / ((a - 1. + m2) * (a + m2));
        d = 1. + aa * d;
        if d.abs() < minimum {
            d = minimum;
        }
        c = 1. + aa / c;
        if c.abs() < minimum {
            c = minimum;
        }
        d = 1. / d;
        h *= d * c;
        let aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1. + m2));
        d = 1. + aa * d;
        if d.abs() < minimum {
            d = minimum;
        }
        c = 1. + aa / c;
        if c.abs() < minimum {
            c = minimum;
        }
        d = 1. / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.).abs() < 2e-14 {
            return finite(h, "beta fraction");
        }
    }
    Err("incomplete beta failed to converge".into())
}

fn regularized_beta(x: f64, a: f64, b: f64) -> Result<f64, String> {
    if x <= 0. {
        return Ok(0.);
    }
    if x >= 1. {
        return Ok(1.);
    }
    let factor =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    let result = if x < (a + 1.) / (a + b + 2.) {
        factor * beta_fraction(a, b, x)? / a
    } else {
        1. - factor * beta_fraction(b, a, 1. - x)? / b
    };
    Ok(finite(result, "beta probability")?.clamp(0., 1.))
}

fn student_two_sided(t: f64, df: f64) -> Result<f64, String> {
    if t.is_nan() {
        return Ok(f64::NAN);
    }
    if t.is_infinite() {
        return Ok(0.);
    }
    regularized_beta(df / (df + t * t), df / 2., 0.5)
}

fn student_quantile(confidence: f64, df: f64) -> Result<f64, String> {
    // Preserve the tail directly; forming (1+confidence)/2 can round to 1.
    let tail = 1. - confidence;
    let mut lower = 0.;
    let mut upper = 1.;
    while student_two_sided(upper, df)? > tail {
        upper *= 2.;
        if upper > 1e150 {
            return Err("confidence quantile exceeds numerical range".into());
        }
    }
    for _ in 0..100 {
        let midpoint = (lower + upper) / 2.;
        if student_two_sided(midpoint, df)? > tail {
            lower = midpoint;
        } else {
            upper = midpoint;
        }
    }
    Ok((lower + upper) / 2.)
}
