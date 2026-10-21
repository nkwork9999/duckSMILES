//! Cox proportional hazards with Breslow/Efron partial-likelihood ties.
//!
//! Rows are independent [features..., nonnegative_time, event_0_or_1]. There is
//! no intercept, left truncation, weights, strata, time-varying covariates,
//! frailty, competing-risk model, or automatic PH-assumption test. Confidence
//! intervals are ordinary 95% observed-information Wald intervals. With ridge,
//! curvature SEs are descriptive and calibrated p-values/CIs are suppressed.
//!
//! Risk sets include same-time censors. Baseline jumps always use Breslow's
//! estimator, referenced to the training feature means, even with Efron fitting.

const MAX_ROWS: usize = 100_000;
const MAX_FEATURES: usize = 32;
const MAX_WORK: usize = 200_000_000;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > MAX_ROWS
        || width == 0
        || width > MAX_FEATURES + 2
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|x| !x.is_finite()))
        || args.iter().any(|x| !x.is_finite())
    {
        return Err(
            "Cox requires rectangular finite numeric data, at most 100000 rows and 32 features"
                .into(),
        );
    }
    match operation {
        "cox_ph" => fit(data, args),
        "cox_predict" | "relative_risk_predict" => predict(operation, data, args),
        _ => Err(format!("unknown survival operation: {operation}")),
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
        Err(format!("{name} must be an integer from {min} to {max}"))
    } else {
        Ok(x as usize)
    }
}
fn numbers(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}
fn nulls(n: usize) -> String {
    format!("[{}]", vec!["null"; n].join(","))
}
fn nullable_exp(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter()
            .map(|&v| {
                let e = v.exp();
                if e.is_finite() {
                    e.to_string()
                } else {
                    "null".into()
                }
            })
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn inner(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(
        a.iter().zip(b).map(|(a, b)| a * b).sum(),
        "linear predictor",
    )
}

struct Problem {
    z: Vec<Vec<f64>>,
    times: Vec<f64>,
    events: Vec<bool>,
    order: Vec<usize>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    means: Vec<f64>,
    scales: Vec<f64>,
    deaths: usize,
    efron: bool,
}
struct Evaluation {
    ll: f64,
    score: Vec<f64>,
    info: Vec<Vec<f64>>,
    baseline: Vec<(f64, f64)>,
}
struct Budget {
    used: usize,
}
impl Budget {
    fn charge(&mut self, n: usize, p: usize) -> Result<(), String> {
        let cost = n
            .checked_mul(p)
            .and_then(|v| v.checked_mul(p + 1))
            .and_then(|v| v.checked_mul(2))
            .ok_or("Cox work count overflow")?;
        self.used = self
            .used
            .checked_add(cost)
            .ok_or("Cox work count overflow")?;
        if self.used > MAX_WORK {
            Err("Cox evaluations exceed 200000000 work units".into())
        } else {
            Ok(())
        }
    }
}
impl Problem {
    fn new(data: &[Vec<f64>], efron: bool) -> Result<Self, String> {
        let p = data[0]
            .len()
            .checked_sub(2)
            .filter(|&p| p > 0 && p <= MAX_FEATURES)
            .ok_or("Cox fitting needs 1 to 32 features plus time and event")?;
        if data.len() < 3 || data.len() <= p {
            return Err("Cox needs at least three observations and more rows than features".into());
        }
        if data
            .iter()
            .any(|r| r[p] < 0. || r[p + 1] != 0. && r[p + 1] != 1.)
        {
            return Err("Cox times must be nonnegative and event indicators binary".into());
        }
        let deaths = data.iter().filter(|r| r[p + 1] == 1.).count();
        if deaths == 0 {
            return Err("Cox model requires at least one event".into());
        }
        let mut z = vec![vec![0.; p]; data.len()];
        let mut origins = Vec::new();
        let mut offsets = Vec::new();
        let mut means = Vec::new();
        let mut scales = Vec::new();
        for j in 0..p {
            let origin = data[0][j];
            let mut offset = 0.;
            for (i, r) in data.iter().enumerate() {
                let relative = finite(r[j] - origin, "origin-relative feature")?;
                offset = finite(
                    offset + (relative - offset) / (i + 1) as f64,
                    "feature mean",
                )?;
            }
            let mut centered = Vec::new();
            let mut max = 0_f64;
            for r in data {
                let v = finite((r[j] - origin) - offset, "centered feature")?;
                max = max.max(v.abs());
                centered.push(v);
            }
            if max == 0. {
                return Err("constant Cox features are not identifiable".into());
            }
            let scale = finite(
                max * (centered.iter().map(|v| (v / max).powi(2)).sum::<f64>()
                    / (data.len() - 1) as f64)
                    .sqrt(),
                "feature scale",
            )?;
            if scale == 0. {
                return Err("Cox feature scale underflow".into());
            }
            for (i, v) in centered.into_iter().enumerate() {
                z[i][j] = finite(v / scale, "standardized feature")?;
            }
            origins.push(origin);
            offsets.push(offset);
            means.push(finite(origin + offset, "feature mean")?);
            scales.push(scale);
        }
        let times = data.iter().map(|r| r[p]).collect::<Vec<_>>();
        let events = data.iter().map(|r| r[p + 1] == 1.).collect::<Vec<_>>();
        let mut order = (0..data.len()).collect::<Vec<_>>();
        // Feature ordering within ties makes floating-point accumulation invariant
        // to permutations of otherwise identical rows.
        order.sort_by(|&i, &j| {
            times[j].total_cmp(&times[i]).then_with(|| {
                for (a, b) in z[i].iter().zip(&z[j]) {
                    let c = a.total_cmp(b);
                    if !c.is_eq() {
                        return c;
                    }
                }
                events[i].cmp(&events[j])
            })
        });
        Ok(Self {
            z,
            times,
            events,
            order,
            origins,
            offsets,
            means,
            scales,
            deaths,
            efron,
        })
    }
    fn evaluate(&self, beta: &[f64], budget: &mut Budget) -> Result<Evaluation, String> {
        let p = beta.len();
        let n = self.z.len();
        budget.charge(n, p)?;
        let eta = self
            .z
            .iter()
            .map(|x| inner(x, beta))
            .collect::<Result<Vec<_>, _>>()?;
        let mut ll = 0.;
        let mut score = vec![0.; p];
        let mut info = vec![vec![0.; p]; p];
        let mut baseline = Vec::new();
        let mut risk_max = f64::NEG_INFINITY;
        let mut weight = 0.;
        let mut first = vec![0.; p];
        let mut second = vec![vec![0.; p]; p];
        let mut i = 0;
        while i < n {
            let start = i;
            let time = self.times[self.order[i]];
            while i < n && self.times[self.order[i]] == time {
                i += 1;
            }
            let group = &self.order[start..i];
            let max = group.iter().map(|&idx| eta[idx]).fold(risk_max, f64::max);
            let rescale = (risk_max - max).exp();
            weight *= rescale;
            for j in 0..p {
                first[j] *= rescale;
                for value in &mut second[j] {
                    *value *= rescale;
                }
            }
            risk_max = max;
            let mut event_weight = 0.;
            let mut event_first = vec![0.; p];
            let mut event_second = vec![vec![0.; p]; p];
            let mut d = 0;
            for &idx in group {
                let w = (eta[idx] - risk_max).exp();
                weight += w;
                for j in 0..p {
                    first[j] += w * self.z[idx][j];
                    for (moment, value) in second[j].iter_mut().zip(&self.z[idx]) {
                        *moment += w * self.z[idx][j] * value;
                    }
                }
                if self.events[idx] {
                    d += 1;
                    event_weight += w;
                    ll += eta[idx] - risk_max;
                    for j in 0..p {
                        score[j] += self.z[idx][j];
                        event_first[j] += w * self.z[idx][j];
                        for (moment, value) in event_second[j].iter_mut().zip(&self.z[idx]) {
                            *moment += w * self.z[idx][j] * value;
                        }
                    }
                }
            }
            if d == 0 {
                continue;
            }
            if weight <= 0. || !weight.is_finite() {
                return Err("Cox risk-set weight underflow or overflow".into());
            }
            let increment = finite(
                ((d as f64).ln() - risk_max - weight.ln()).exp(),
                "Breslow baseline hazard jump",
            )?;
            baseline.push((time, increment));
            for l in 0..d {
                let fraction = if self.efron { l as f64 / d as f64 } else { 0. };
                let denominator = weight - fraction * event_weight;
                if denominator <= 0. {
                    return Err("Cox tied risk-set denominator lost precision".into());
                }
                ll -= denominator.ln();
                let average = (0..p)
                    .map(|j| (first[j] - fraction * event_first[j]) / denominator)
                    .collect::<Vec<_>>();
                for j in 0..p {
                    score[j] -= average[j];
                    for k in 0..p {
                        info[j][k] += (second[j][k] - fraction * event_second[j][k]) / denominator
                            - average[j] * average[k];
                    }
                }
            }
        }
        finite(ll, "partial log likelihood")?;
        if score
            .iter()
            .chain(info.iter().flatten())
            .any(|v| !v.is_finite())
        {
            return Err("Cox derivatives exceeded finite range".into());
        }
        Ok(Evaluation {
            ll,
            score,
            info,
            baseline,
        })
    }
}
type PenalizedEvaluation = (f64, Vec<f64>, Vec<Vec<f64>>);

fn penalized(e: &Evaluation, beta: &[f64], lambda: f64) -> Result<PenalizedEvaluation, String> {
    let objective = finite(
        e.ll - 0.5 * lambda * inner(beta, beta)?,
        "penalized partial likelihood",
    )?;
    let gradient = e
        .score
        .iter()
        .zip(beta)
        .map(|(s, b)| s - lambda * b)
        .collect::<Vec<_>>();
    let mut info = e.info.clone();
    for (i, row) in info.iter_mut().enumerate() {
        row[i] += lambda;
    }
    Ok((objective, gradient, info))
}
fn cholesky(matrix: &[Vec<f64>], threshold: f64) -> Result<Vec<Vec<f64>>, String> {
    let p = matrix.len();
    let mut lower = vec![vec![0.; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let value = matrix[i][j] - (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
            if i == j {
                if !value.is_finite() || value <= threshold {
                    return Err(
                        "Cox observed information is singular or numerically unidentified".into(),
                    );
                }
                lower[i][j] = value.sqrt();
            } else {
                lower[i][j] = finite(value / lower[j][j], "information Cholesky factor")?;
            }
        }
    }
    Ok(lower)
}
fn solve(lower: &[Vec<f64>], rhs: &[f64]) -> Result<Vec<f64>, String> {
    let p = rhs.len();
    let mut x = rhs.to_vec();
    for i in 0..p {
        x[i] = (x[i] - (0..i).map(|j| lower[i][j] * x[j]).sum::<f64>()) / lower[i][i];
    }
    for i in (0..p).rev() {
        x[i] = (x[i] - (i + 1..p).map(|j| lower[j][i] * x[j]).sum::<f64>()) / lower[i][i];
        finite(x[i], "information solution")?;
    }
    Ok(x)
}
fn fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 4 {
        return Err("Cox arguments are [ties=1, lambda=0, max_iter=100, tolerance=1e-8]".into());
    }
    let ties = integer(args.first().copied().unwrap_or(1.), 0, 1, "ties method")?;
    let lambda = args.get(1).copied().unwrap_or(0.);
    if !(0. ..=1e12).contains(&lambda) {
        return Err("Cox ridge lambda must be in [0,1e12]".into());
    }
    let maxiter = integer(
        args.get(2).copied().unwrap_or(100.),
        1,
        1000,
        "maximum iterations",
    )?;
    let tolerance = args.get(3).copied().unwrap_or(1e-8);
    if !(1e-12..=1e-3).contains(&tolerance) {
        return Err("Cox tolerance must be in [1e-12,1e-3]".into());
    }
    let problem = Problem::new(data, ties == 1)?;
    let p = problem.scales.len();
    let mut beta = vec![0.; p];
    let mut budget = Budget { used: 0 };
    let mut evaluation = problem.evaluate(&beta, &mut budget)?;
    let null_ll = evaluation.ll;
    let initial_information = (0..p).map(|j| evaluation.info[j][j]).fold(0., f64::max);
    if initial_information <= 0. && lambda == 0. {
        return Err("Cox has no identifiable event information".into());
    }
    let threshold = 1e-12 * (initial_information + lambda).max(1.);
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 0..maxiter {
        iterations = iteration + 1;
        let (objective, gradient, information) = penalized(&evaluation, &beta, lambda)?;
        let lower = cholesky(&information, threshold)?;
        let direction = solve(&lower, &gradient)?;
        let max_gradient = gradient.iter().map(|v| v.abs()).fold(0., f64::max);
        let max_step = direction.iter().map(|v| v.abs()).fold(0., f64::max);
        let max_beta = beta.iter().map(|v| v.abs()).fold(0., f64::max);
        if max_gradient <= tolerance * (1. + problem.deaths as f64)
            && max_step <= tolerance * (1. + max_beta)
        {
            converged = true;
            break;
        }
        if lambda == 0. && max_beta > 50. {
            return Err("Cox unpenalized coefficients diverged; separation is possible".into());
        }
        let slope = inner(&gradient, &direction)?;
        if slope < 0. {
            return Err("Cox Newton direction did not improve likelihood".into());
        }
        let mut step = if max_step > 5. { 5. / max_step } else { 1. };
        let mut accepted = None;
        for _ in 0..40 {
            let next = beta
                .iter()
                .zip(&direction)
                .map(|(b, d)| finite(b + step * d, "Cox coefficient update"))
                .collect::<Result<Vec<_>, _>>()?;
            let next_eval = problem.evaluate(&next, &mut budget)?;
            let next_objective = penalized(&next_eval, &next, lambda)?.0;
            if next_objective >= objective + 1e-4 * step * slope - 1e-12 {
                accepted = Some((next, next_eval));
                break;
            }
            step *= 0.5;
        }
        match accepted {
            Some((next, next_eval)) => {
                beta = next;
                evaluation = next_eval;
            }
            None => return Err("Cox line search failed to improve partial likelihood".into()),
        }
    }
    if !converged {
        return Err("Cox fit did not converge within maximum iterations".into());
    }
    let (_, _, information) = penalized(&evaluation, &beta, lambda)?;
    let lower = cholesky(&information, threshold)?;
    let mut covariance = vec![vec![0.; p]; p];
    for j in 0..p {
        let mut unit = vec![0.; p];
        unit[j] = 1.;
        let column = solve(&lower, &unit)?;
        for i in 0..p {
            covariance[i][j] = finite(
                (column[i] / problem.scales[i]) / problem.scales[j],
                "coefficient covariance",
            )?;
        }
    }
    let coefficients = beta
        .iter()
        .zip(&problem.scales)
        .map(|(b, s)| finite(b / s, "original-unit Cox coefficient"))
        .collect::<Result<Vec<_>, _>>()?;
    let se = (0..p)
        .map(|j| finite(covariance[j][j].sqrt(), "coefficient standard error"))
        .collect::<Result<Vec<_>, _>>()?;
    let z = coefficients
        .iter()
        .zip(&se)
        .map(|(b, s)| finite(b / s, "Wald z statistic"))
        .collect::<Result<Vec<_>, _>>()?;
    let (pvalues, intervals, inference) = if lambda == 0. {
        let pvalues = z
            .iter()
            .map(|v| gamma_sf(0.5, v * v / 2.))
            .collect::<Result<Vec<_>, _>>()?;
        let intervals = coefficients
            .iter()
            .zip(&se)
            .map(|(b, s)| {
                Ok(vec![
                    finite(
                        b - 1.959_963_984_540_054 * s,
                        "lower coefficient confidence limit",
                    )?,
                    finite(
                        b + 1.959_963_984_540_054 * s,
                        "upper coefficient confidence limit",
                    )?,
                ])
            })
            .collect::<Result<Vec<_>, String>>()?;
        (
            numbers(&pvalues),
            matrix(&intervals),
            "observed_information_wald_asymptotic",
        )
    } else {
        (
            nulls(p),
            "null".into(),
            "penalized_curvature_no_calibrated_inference",
        )
    };
    let lr = finite(
        (2. * (evaluation.ll - null_ll)).max(0.),
        "likelihood ratio statistic",
    )?;
    let lr_p = if lambda == 0. {
        gamma_sf(p as f64 / 2., lr / 2.)?.to_string()
    } else {
        "null".into()
    };
    let mut baseline = evaluation.baseline;
    baseline.reverse();
    let mut cumulative = 0.;
    let mut times = Vec::new();
    let mut hazards = Vec::new();
    let mut survival = Vec::new();
    for (time, jump) in baseline {
        cumulative = finite(cumulative + jump, "baseline cumulative hazard")?;
        times.push(time);
        hazards.push(cumulative);
        survival.push((-cumulative).exp());
    }
    let tied = if ties == 0 { "breslow" } else { "efron" };
    Ok(format!("{{\"version\":1,\"algorithm\":\"cox_proportional_hazards_partial_likelihood\",\"ties\":\"{tied}\",\"lambda\":{lambda},\"penalty_scale\":\"sum_partial_log_likelihood_minus_lambda_half_squared_standardized_coefficients\",\"n_samples\":{},\"n_events\":{},\"n_features\":{p},\"coefficients\":{},\"standardized_coefficients\":{},\"feature_origins\":{},\"mean_offsets\":{},\"feature_means\":{},\"feature_scales\":{},\"covariance_matrix\":{},\"standard_errors\":{},\"z_statistics\":{},\"p_values\":{pvalues},\"confidence_level\":0.95,\"confidence_intervals\":{intervals},\"hazard_ratios\":{},\"inference_method\":\"{inference}\",\"partial_log_likelihood\":{},\"null_partial_log_likelihood\":{null_ll},\"likelihood_ratio_statistic\":{lr},\"likelihood_ratio_degrees_of_freedom\":{p},\"likelihood_ratio_p_value\":{lr_p},\"baseline_hazard_estimator\":\"breslow\",\"baseline_reference\":\"training_feature_origins_plus_mean_offsets\",\"baseline_times\":{},\"baseline_cumulative_hazard\":{},\"baseline_survival\":{},\"converged\":true,\"iterations\":{iterations},\"work_units\":{}}}",data.len(),problem.deaths,numbers(&coefficients),numbers(&beta),numbers(&problem.origins),numbers(&problem.offsets),numbers(&problem.means),numbers(&problem.scales),matrix(&covariance),numbers(&se),numbers(&z),nullable_exp(&coefficients),evaluation.ll,numbers(&times),numbers(&hazards),numbers(&survival),budget.used))
}
fn predict(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data[0].len();
    if p > MAX_FEATURES {
        return Err("Cox prediction has at most 32 features".into());
    }
    let offset = usize::from(operation == "cox_predict");
    if args.len() != 3 * p + offset {
        return Err("prediction arguments are [optional baseline_survival, coefficients(p), origins(p), mean_offsets(p)]".into());
    }
    let baseline = if offset == 1 { args[0] } else { 1. };
    if !(0. ..=1.).contains(&baseline) {
        return Err("baseline survival must be in [0,1]".into());
    }
    let coefficients = &args[offset..offset + p];
    let origins = &args[offset + p..offset + 2 * p];
    let means = &args[offset + 2 * p..];
    let mut linear = Vec::new();
    let mut risk = Vec::new();
    let mut survival = Vec::new();
    for row in data {
        let lp = finite(
            (0..p)
                .map(|j| coefficients[j] * ((row[j] - origins[j]) - means[j]))
                .sum(),
            "centered Cox linear predictor",
        )?;
        let relative = finite(lp.exp(), "relative risk")?;
        if relative == 0. {
            return Err(
                "relative risk underflow; linear predictor is outside representable range".into(),
            );
        }
        let probability = if baseline == 0. {
            0.
        } else if baseline == 1. {
            1.
        } else {
            (baseline.ln() * relative).exp()
        };
        linear.push(lp);
        risk.push(relative);
        survival.push(probability);
    }
    let extra = if offset == 1 {
        format!(
            ",\"baseline_survival\":{baseline},\"survival_probabilities\":{}",
            numbers(&survival)
        )
    } else {
        String::new()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"centered_cox_prediction\",\"linear_predictors\":{},\"relative_risks\":{}{extra}}}",numbers(&linear),numbers(&risk)))
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
    0.5 * (2. * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}
fn gamma_sf(a: f64, x: f64) -> Result<f64, String> {
    if x == 0. {
        return Ok(1.);
    }
    if x.is_infinite() {
        return Ok(0.);
    }
    let front = (a * x.ln() - x - log_gamma(a)).exp();
    let q;
    if x < a + 1. {
        let mut term = 1. / a;
        let mut sum = term;
        let mut done = false;
        for i in 1..=20_000 {
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() < sum.abs() * 2e-14 {
                done = true;
                break;
            }
        }
        if !done {
            return Err("chi-square gamma series failed to converge".into());
        }
        q = 1. - sum * front;
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
        let mut done = false;
        for i in 1..=20_000 {
            let i = i as f64;
            let an = -i * (i - a);
            b += 2.;
            d = 1. / nz(an * d + b);
            c = nz(b + an / c);
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-14 {
                done = true;
                break;
            }
        }
        if !done {
            return Err("chi-square gamma fraction failed to converge".into());
        }
        q = front * h;
    }
    if q.is_finite() && (-1e-10..=1. + 1e-10).contains(&q) {
        Ok(q.clamp(0., 1.))
    } else {
        Err("chi-square probability outside numerical range".into())
    }
}
