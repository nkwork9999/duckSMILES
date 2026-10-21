//! Gaussian random-intercept mixed models with arbitrary group sizes and fixed covariates.
//!
//! Fit rows: `[fixed_features..., response, group_id]`; an intercept is included.
//! Parameters: `[method=1 (0 ML, 1 REML), max_iterations=200, tolerance=1e-8]`.
//! ML and REML profile over the nonnegative random/residual variance ratio, including
//! the zero boundary. Block eigenscaling plus Householder QR computes GLS without
//! constructing a dense marginal covariance or squaring the design condition number.
//! Standard errors condition on the fitted variance components; no small-sample
//! degrees-of-freedom approximation, random slopes, crossed effects, or GLMM is implied.
//!
//! Prediction rows are `[fixed_features..., group_id]`; parameters are the returned
//! `predict_parameters`. Known groups use fitted BLUPs, new groups use zero random
//! effect. `marginal_predict_parameters` always predicts with zero random effects.
//! Portable layout: `[p,y_origin,y_mean_offset,centered_intercept,feature_origins...,
//! feature_mean_offsets...,original_slopes...,n_groups,group_id,BLUP,...]`.

use std::collections::BTreeMap;

const MAX_CELLS: usize = 1_000_000;
const MAX_FEATURES: usize = 64;
const MAX_WORK: usize = 200_000_000;
const MAX_INTEGER: f64 = 9_007_199_254_740_991.;

/// Runs a bounded mixed-model fit or prediction and returns version-1 JSON.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > MAX_CELLS / width
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
    {
        return Err("data must be nonempty, rectangular, finite, and at most 1000000 cells".into());
    }
    if args.iter().any(|v| !v.is_finite()) {
        return Err("parameters must be finite".into());
    }
    match operation {
        "mixed_model" => fit(data, args),
        "mixed_predict" => predict(data, args),
        _ => Err(format!("unknown mixed-model operation: {operation}")),
    }
}

fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
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
fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0_f64, f64::hypot)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn numbers(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(f64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn matrix(values: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|r| numbers(r))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn group_id(value: f64) -> Result<i64, String> {
    if value.abs() > MAX_INTEGER || value.fract() != 0. {
        Err("group ids must be exact signed integers".into())
    } else {
        Ok(value as i64)
    }
}
fn exact_usize(value: f64, name: &str) -> Result<usize, String> {
    if !(0. ..=MAX_INTEGER).contains(&value) || value.fract() != 0. {
        Err(format!("{name} must be an exact nonnegative integer"))
    } else {
        Ok(value as usize)
    }
}

struct Design {
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    y_origin: f64,
    y_offset: f64,
    groups: BTreeMap<i64, Vec<usize>>,
    original_design_log_scale: f64,
}
impl Design {
    fn new(data: &[Vec<f64>]) -> Result<Self, String> {
        let width = data[0].len();
        if width < 2 {
            return Err("mixed_model requires response and group id columns".into());
        }
        let p = width - 2;
        let n = data.len();
        let q = p + 1;
        if p > MAX_FEATURES || n <= q + 1 {
            return Err(
                "mixed_model requires <=64 features and at least two residual degrees of freedom"
                    .into(),
            );
        }
        let mut groups = BTreeMap::<i64, Vec<usize>>::new();
        for (i, row) in data.iter().enumerate() {
            groups.entry(group_id(row[p + 1])?).or_default().push(i);
        }
        if groups.len() < 2 || groups.len() >= n || groups.len() > 10000 {
            return Err(
                "mixed model requires 2 to 10000 groups and within-group replication".into(),
            );
        }
        let origins = data[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|r| r[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let value = finite(
                    norm(data.iter().map(|r| (r[j] - origins[j]) - offsets[j])) / (n as f64).sqrt(),
                    "feature scale",
                )?;
                if value == 0. {
                    Err(
                        "fixed design contains a constant feature in addition to the intercept"
                            .into(),
                    )
                } else {
                    Ok(value)
                }
            })
            .collect::<Result<Vec<_>, String>>()?;
        let y_origin = data[0][p];
        let y_offset = mean(data.iter().map(|r| r[p] - y_origin))?;
        let y = data
            .iter()
            .map(|r| finite((r[p] - y_origin) - y_offset, "centered response"))
            .collect::<Result<Vec<_>, _>>()?;
        if norm(y.iter().copied()) == 0. {
            return Err("constant response has zero residual variance".into());
        }
        let x = data
            .iter()
            .map(|r| {
                let mut row = vec![1.];
                for j in 0..p {
                    row.push(finite(
                        ((r[j] - origins[j]) - offsets[j]) / scales[j],
                        "scaled feature",
                    )?);
                }
                Ok(row)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let original_design_log_scale = finite(
            2. * scales.iter().map(|s| s.ln()).sum::<f64>(),
            "fixed design log scale",
        )?;
        Ok(Self {
            x,
            y,
            origins,
            offsets,
            scales,
            y_origin,
            y_offset,
            groups,
            original_design_log_scale,
        })
    }
}

struct Qr {
    beta: Vec<f64>,
    r: Vec<Vec<f64>>,
}
#[allow(clippy::needless_range_loop)]
fn qr_solve(mut x: Vec<Vec<f64>>, mut y: Vec<f64>) -> Result<Qr, String> {
    let p = x[0].len();
    if x.len() < p {
        return Err("rank-deficient fixed design".into());
    }
    let threshold =
        f64::EPSILON * 64. * (x.len().max(p) as f64).sqrt() * norm(x.iter().flatten().copied());
    for j in 0..p {
        let length = finite(norm(x[j..].iter().map(|r| r[j])), "QR norm")?;
        if length <= threshold {
            return Err("rank-deficient or numerically singular fixed design".into());
        }
        let diagonal = -length.copysign(x[j][j]);
        let mut reflector = x[j..].iter().map(|r| r[j] / length).collect::<Vec<_>>();
        reflector[0] -= diagonal / length;
        let norm = norm(reflector.iter().copied());
        for v in &mut reflector {
            *v /= norm;
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
    let mut beta = vec![0.; p];
    for j in (0..p).rev() {
        beta[j] = finite(
            (y[j] - dot(&x[j][j + 1..], &beta[j + 1..])) / x[j][j],
            "GLS coefficient",
        )?;
    }
    Ok(Qr {
        beta,
        r: x[..p].to_vec(),
    })
}

struct Profile {
    rho: f64,
    sigma: f64,
    log_likelihood: f64,
    beta: Vec<f64>,
    r: Vec<Vec<f64>>,
}

fn profile(design: &Design, t: f64, reml: bool) -> Result<Profile, String> {
    let n = design.y.len();
    let q = design.x[0].len();
    let rho = t.exp_m1();
    let mut x = design.x.clone();
    let mut y = design.y.clone();
    let mut logdet_w = 0.;
    for indices in design.groups.values() {
        let factor = finite(1. + rho * indices.len() as f64, "group variance ratio")?;
        let shrink = factor.sqrt();
        logdet_w = finite(logdet_w + factor.ln(), "marginal log determinant")?;
        let group_y = mean(indices.iter().map(|i| design.y[*i]))?;
        let group_x = (0..q)
            .map(|j| mean(indices.iter().map(|i| design.x[*i][j])))
            .collect::<Result<Vec<_>, _>>()?;
        for &i in indices {
            y[i] = finite(
                (design.y[i] - group_y) + group_y / shrink,
                "weighted response",
            )?;
            for j in 0..q {
                x[i][j] = finite(
                    (design.x[i][j] - group_x[j]) + group_x[j] / shrink,
                    "weighted design",
                )?;
            }
        }
    }
    let qr = qr_solve(x.clone(), y.clone())?;
    let rss = finite(
        norm(x.iter().zip(&y).map(|(r, v)| v - dot(r, &qr.beta))).powi(2),
        "weighted residual sum of squares",
    )?;
    if rss <= 0. {
        return Err("residual variance collapsed to zero; no finite likelihood optimum".into());
    }
    let df = if reml { n - q } else { n };
    let sigma = finite(rss / df as f64, "residual variance")?;
    if sigma <= 0. {
        return Err("residual variance underflow".into());
    }
    let logdet_x = if reml {
        finite(
            2. * (0..q).map(|j| qr.r[j][j].abs().ln()).sum::<f64>()
                + design.original_design_log_scale,
            "restricted design determinant",
        )?
    } else {
        0.
    };
    let log_likelihood = finite(
        -0.5 * (df as f64 * ((2. * std::f64::consts::PI).ln() + 1. + sigma.ln())
            + logdet_w
            + logdet_x),
        "log likelihood",
    )?;
    Ok(Profile {
        rho,
        sigma,
        log_likelihood,
        beta: qr.beta,
        r: qr.r,
    })
}

fn optimize(
    design: &Design,
    reml: bool,
    max_iterations: usize,
    tolerance: f64,
) -> Result<(Profile, bool, usize, bool), String> {
    // Wide deterministic grid finds a local bracket before a safeguarded golden
    // search. The upper ratio is explicit; hitting it never reports convergence.
    let grid = [0., 0.25, 0.5, 1., 2., 4., 8., 12., 16., 24., 32.];
    let mut candidates = Vec::with_capacity(grid.len());
    for t in grid {
        candidates.push(profile(design, t, reml)?);
    }
    let best_index = (0..candidates.len())
        .max_by(|a, b| {
            candidates[*a]
                .log_likelihood
                .total_cmp(&candidates[*b].log_likelihood)
        })
        .ok_or("no likelihood candidate")?;
    let at_limit = best_index == grid.len() - 1;
    let mut lower = grid[best_index.saturating_sub(1)];
    let mut upper = grid[(best_index + 1).min(grid.len() - 1)];
    let mut best = candidates.swap_remove(best_index);
    let golden = (5_f64.sqrt() - 1.) / 2.;
    let mut a = upper - golden * (upper - lower);
    let mut b = lower + golden * (upper - lower);
    let mut fa = profile(design, a, reml)?;
    let mut fb = profile(design, b, reml)?;
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..max_iterations {
        iterations = iteration + 1;
        if upper - lower <= tolerance * (1. + lower.abs() + upper.abs()) {
            converged = true;
            break;
        }
        if fa.log_likelihood >= fb.log_likelihood {
            upper = b;
            b = a;
            fb = fa;
            a = upper - golden * (upper - lower);
            fa = profile(design, a, reml)?;
        } else {
            lower = a;
            a = b;
            fa = fb;
            b = lower + golden * (upper - lower);
            fb = profile(design, b, reml)?;
        }
    }
    if fa.log_likelihood > best.log_likelihood {
        best = fa;
    }
    if fb.log_likelihood > best.log_likelihood {
        best = fb;
    }
    Ok((best, converged && !at_limit, iterations, at_limit))
}

fn fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 3 {
        return Err("mixed_model requires [method=1,max_iterations=200,tolerance=1e-8]".into());
    }
    let method = args.first().copied().unwrap_or(1.);
    let max = args.get(1).copied().unwrap_or(200.);
    let tolerance = args.get(2).copied().unwrap_or(1e-8);
    if method != 0. && method != 1. {
        return Err("method must be 0 ML or 1 REML".into());
    }
    let max_iterations = exact_usize(max, "max_iterations")?;
    if max_iterations == 0 || max_iterations > 10000 || !(0. < tolerance && tolerance < 1.) {
        return Err("max_iterations must be 1..10000 and tolerance between 0 and 1".into());
    }
    let design = Design::new(data)?;
    let n = data.len();
    let p = design.scales.len();
    let q = p + 1;
    if n.saturating_mul(q)
        .saturating_mul(q)
        .saturating_mul(max_iterations + 13)
        > MAX_WORK
    {
        return Err("mixed model exceeds numerical work limit".into());
    }
    let (model, converged, iterations, search_limit) =
        optimize(&design, method == 1., max_iterations, tolerance)?;
    let coefficients = model.beta[1..]
        .iter()
        .zip(&design.scales)
        .map(|(b, s)| finite(b / s, "original coefficient"))
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        (design.y_origin - dot(&design.origins, &coefficients))
            + (design.y_offset + model.beta[0] - dot(&design.offsets, &coefficients)),
        "intercept",
    )?;
    let random_variance = finite(model.rho * model.sigma, "random intercept variance")?;
    let mut random = BTreeMap::<i64, f64>::new();
    for (group, indices) in &design.groups {
        let residual_mean = mean(
            indices
                .iter()
                .map(|i| design.y[*i] - dot(&design.x[*i], &model.beta)),
        )?;
        let ratio = model.rho * indices.len() as f64;
        random.insert(
            *group,
            finite((ratio / (1. + ratio)) * residual_mean, "group BLUP")?,
        );
    }
    let mut marginal = Vec::with_capacity(n);
    let mut conditional = Vec::with_capacity(n);
    let mut residuals = Vec::with_capacity(n);
    for (i, row) in data.iter().enumerate() {
        let fixed = finite(dot(&design.x[i], &model.beta), "fixed prediction")?;
        let u = random[&group_id(row[p + 1])?];
        marginal.push(finite(
            design.y_origin + (design.y_offset + fixed),
            "marginal fitted value",
        )?);
        conditional.push(finite(
            design.y_origin + (design.y_offset + fixed + u),
            "conditional fitted value",
        )?);
        residuals.push(finite(design.y[i] - fixed - u, "conditional residual")?);
    }
    let (covariance, standard_errors) = fixed_uncertainty(&design, &model)?;
    let mut predict_parameters = vec![p as f64, design.y_origin, design.y_offset, model.beta[0]];
    predict_parameters.extend_from_slice(&design.origins);
    predict_parameters.extend_from_slice(&design.offsets);
    predict_parameters.extend_from_slice(&coefficients);
    let mut marginal_parameters = predict_parameters.clone();
    marginal_parameters.push(0.);
    predict_parameters.push(random.len() as f64);
    for (group, u) in &random {
        predict_parameters.push(*group as f64);
        predict_parameters.push(*u);
    }
    let mut stable = vec![
        0.,
        design.y_origin,
        finite(design.y_offset + model.beta[0], "prediction center")?,
    ];
    stable.extend_from_slice(&design.origins);
    stable.extend_from_slice(&design.offsets);
    stable.extend_from_slice(&coefficients);
    let group_ids = random.keys().map(|g| *g as f64).collect::<Vec<_>>();
    let blups = random.values().copied().collect::<Vec<_>>();
    let aic = if method == 0. {
        finite(-2. * model.log_likelihood + 2. * (q + 2) as f64, "AIC")?.to_string()
    } else {
        "null".into()
    };
    let boundary = model.rho == 0.;
    let log_likelihood = model.log_likelihood;
    let sigma = model.sigma;
    let rho = model.rho;
    Ok(format!("{{\"version\":1,\"algorithm\":\"mixed_model\",\"method\":\"{}\",\"n_samples\":{n},\"n_fixed_features\":{p},\"n_groups\":{},\"intercept\":{intercept},\"coefficients\":{},\"residual_variance\":{sigma},\"random_intercept_variance\":{random_variance},\"variance_ratio\":{rho},\"variance_boundary\":{boundary},\"log_likelihood\":{log_likelihood},\"aic\":{aic},\"aic_defined_for\":\"ML_only\",\"intercept_standard_error\":{},\"coefficient_standard_errors\":{},\"fixed_effect_covariance\":{},\"groups\":{},\"random_intercepts\":{},\"marginal_fitted_values\":{},\"conditional_fitted_values\":{},\"conditional_residuals\":{},\"predict_parameters\":{},\"marginal_predict_parameters\":{},\"stable_predict_parameters\":{},\"converged\":{converged},\"iterations\":{iterations},\"max_iterations\":{max_iterations},\"tolerance\":{tolerance},\"search_limit_reached\":{search_limit},\"maximum_log1p_variance_ratio\":32,\"standard_error_method\":\"GLS_conditional_on_variance_estimates\",\"random_effect_structure\":\"one_grouping_factor_random_intercept\"}}",if method==1.{"REML"}else{"ML"},random.len(),numbers(&coefficients),standard_errors[0],numbers(&standard_errors[1..]),matrix(&covariance),numbers(&group_ids),numbers(&blups),numbers(&marginal),numbers(&conditional),numbers(&residuals),numbers(&predict_parameters),numbers(&marginal_parameters),numbers(&stable)))
}

#[allow(clippy::needless_range_loop)]
fn fixed_uncertainty(
    design: &Design,
    model: &Profile,
) -> Result<(Vec<Vec<f64>>, Vec<f64>), String> {
    let q = model.beta.len();
    let p = q - 1;
    let mut inverse = vec![vec![0.; q]; q];
    for k in 0..q {
        for j in (0..q).rev() {
            let rhs = if j == k { 1. } else { 0. };
            let sum = (j + 1..q)
                .map(|i| model.r[j][i] * inverse[i][k])
                .sum::<f64>();
            inverse[j][k] = finite((rhs - sum) / model.r[j][j], "GLS covariance")?;
        }
    }
    let mut transform = vec![vec![0.; q]; q];
    transform[0][0] = 1.;
    for j in 0..p {
        transform[0][j + 1] = finite(
            -(design.origins[j] + design.offsets[j]) / design.scales[j],
            "intercept covariance",
        )?;
        transform[j + 1][j + 1] = 1. / design.scales[j];
    }
    let mut product = vec![vec![0.; q]; q];
    for i in 0..q {
        for j in 0..q {
            product[i][j] = finite(
                (0..q).map(|k| transform[i][k] * inverse[k][j]).sum(),
                "original-unit covariance",
            )?;
        }
    }
    let mut covariance = vec![vec![0.; q]; q];
    for i in 0..q {
        for j in 0..=i {
            let value = finite(
                model.sigma * dot(&product[i], &product[j]),
                "fixed effect covariance",
            )?;
            covariance[i][j] = value;
            covariance[j][i] = value;
        }
    }
    let standard_errors = (0..q)
        .map(|j| {
            finite(
                model.sigma.sqrt() * norm(product[j].iter().copied()),
                "fixed effect standard error",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok((covariance, standard_errors))
}

fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 5 {
        return Err("mixed_predict requires predict_parameters returned by mixed_model".into());
    }
    let p = exact_usize(args[0], "model feature count")?;
    if p > MAX_FEATURES || data[0].len() != p + 1 || args.len() < 5 + 3 * p {
        return Err("prediction feature count does not match model".into());
    }
    let count = exact_usize(args[4 + 3 * p], "model group count")?;
    if count > 10000 || args.len() != 5 + 3 * p + 2 * count {
        return Err("invalid portable mixed model length".into());
    }
    let mut random = BTreeMap::<i64, f64>::new();
    for pair in args[5 + 3 * p..].chunks_exact(2) {
        if random.insert(group_id(pair[0])?, pair[1]).is_some() {
            return Err("portable model has duplicate group ids".into());
        }
    }
    let mut marginal = Vec::with_capacity(data.len());
    let mut conditional = Vec::with_capacity(data.len());
    let mut known = Vec::with_capacity(data.len());
    for row in data {
        let adjustment = (0..p)
            .map(|j| ((row[j] - args[4 + j]) - args[4 + p + j]) * args[4 + 2 * p + j])
            .sum::<f64>();
        let fixed = finite(args[3] + adjustment, "fixed prediction")?;
        let group = group_id(row[p])?;
        let effect = random.get(&group).copied();
        marginal.push(finite(args[1] + (args[2] + fixed), "marginal prediction")?);
        conditional.push(finite(
            args[1] + (args[2] + fixed + effect.unwrap_or(0.)),
            "conditional prediction",
        )?);
        known.push(effect.is_some().to_string());
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"mixed_predict\",\"marginal_predictions\":{},\"conditional_predictions\":{},\"known_group\":[{}],\"new_group_random_effect\":0}}",numbers(&marginal),numbers(&conditional),known.join(",")))
}
