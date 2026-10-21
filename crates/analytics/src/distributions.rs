//! Bounded unweighted distribution and right-censored lifetime likelihoods.
//! Location is fixed at zero for Gamma, Weibull, exponential and lognormal.
//! No generic fitted-distribution KS p-value is supplied: its null calibration
//! would require accounting for estimated parameters.

use std::f64::consts::PI;

const LIMIT: usize = 100_000;
const LOG_SQRT_TAU: f64 = 0.9189385332046727;

fn check(data: &[Vec<f64>], columns: usize, args: &[f64]) -> Result<(), String> {
    if data.is_empty()
        || data.len() > LIMIT
        || data
            .iter()
            .any(|r| r.len() != columns || r.iter().any(|x| !x.is_finite() || x.abs() > 1e150))
        || args.iter().any(|x| !x.is_finite())
    {
        return Err("expected bounded nonempty finite numeric rows of the specified width".into());
    }
    Ok(())
}

fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("numerical result is not representable".into())
    }
}

fn mean_sd(x: &[f64]) -> Result<(f64, f64), String> {
    let origin = x[0];
    let m = finite(origin + x.iter().map(|v| v - origin).sum::<f64>() / x.len() as f64)?;
    let span = x.iter().map(|v| (v - m).abs()).fold(0.0_f64, f64::max);
    if span == 0.0 {
        return Err("positive identifiable variance is required".into());
    }
    let sd = finite(
        span * (x.iter().map(|v| ((v - m) / span).powi(2)).sum::<f64>() / x.len() as f64).sqrt(),
    )?;
    if sd <= 0.0 {
        return Err("positive identifiable variance is required".into());
    }
    Ok((m, sd))
}

fn vec_json(v: &[f64]) -> Result<String, String> {
    v.iter()
        .map(|x| finite(*x).map(|x| x.to_string()))
        .collect::<Result<Vec<_>, _>>()
        .map(|v| format!("[{}]", v.join(",")))
}

fn fitted(
    operation: &str,
    n: usize,
    params: &[(&str, f64)],
    ll: f64,
    extra: &str,
) -> Result<String, String> {
    let k = params.len() as f64;
    let aic = finite(2.0 * k - 2.0 * ll)?;
    let bic = finite(k * (n as f64).ln() - 2.0 * ll)?;
    finite(ll)?;
    let fields = params
        .iter()
        .map(|(key, value)| finite(*value).map(|v| format!("\"{key}\":{v}")))
        .collect::<Result<Vec<_>, _>>()?
        .join(",");
    Ok(format!("{{\"version\":1,\"operation\":\"{operation}\",\"n\":{n},{fields},\"parameters\":{{{fields}}},\"log_likelihood\":{ll},\"aic\":{aic},\"bic\":{bic},\"converged\":true{extra}}}"))
}

// Lanczos log-gamma for positive arguments, used in count/Gamma likelihoods.
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
        return PI.ln() - (PI * z).sin().ln() - log_gamma(1.0 - z);
    }
    let z = z - 1.0;
    let a = C
        .iter()
        .enumerate()
        .fold(0.9999999999998099, |a, (i, c)| a + c / (z + i as f64 + 1.0));
    let t = z + 7.5;
    LOG_SQRT_TAU + (z + 0.5) * t.ln() - t + a.ln()
}

fn poisson_log_mass(count: f64, mean: f64) -> f64 {
    if count == 0.0 {
        return -mean;
    }
    // log(k!) = (k+1/2)log(k)-k+log(sqrt(2*pi))+error(k).
    // Express the leading cancellation in a nonnegative deviance term.
    let stirling_error = if count > 15.0 {
        let r = 1.0 / count;
        let r2 = r * r;
        r * (1.0 / 12.0 - r2 * (1.0 / 360.0 - r2 * (1.0 / 1260.0 - r2 / 1680.0)))
    } else {
        log_gamma(count + 1.0) - (count + 0.5) * count.ln() + count - LOG_SQRT_TAU
    };
    let deviance = if (count - mean).abs() < 0.1 * (count + mean) {
        let v = (count - mean) / (count + mean);
        let v2 = v * v;
        let mut sum = (count - mean) * v;
        let mut term = 2.0 * count * v;
        for j in 1..100 {
            term *= v2;
            let next = sum + term / (2 * j + 1) as f64;
            if next == sum {
                break;
            }
            sum = next;
        }
        sum
    } else {
        count * (count / mean).ln() + mean - count
    };
    -stirling_error - deviance - LOG_SQRT_TAU - 0.5 * count.ln()
}

fn digamma_trigamma(mut x: f64) -> (f64, f64) {
    let (mut psi, mut tri) = (0.0, 0.0);
    while x < 10.0 {
        psi -= 1.0 / x;
        tri += 1.0 / (x * x);
        x += 1.0;
    }
    let r = 1.0 / x;
    let r2 = r * r;
    psi +=
        x.ln() - 0.5 * r - r2 * (1.0 / 12.0 - r2 * (1.0 / 120.0 - r2 * (1.0 / 252.0 - r2 / 240.0)));
    tri += r + 0.5 * r2 + r * r2 * (1.0 / 6.0 - r2 * (1.0 / 30.0 - r2 * (1.0 / 42.0 - r2 / 30.0)));
    (psi, tri)
}

fn root_decreasing<F>(mut equation: F) -> Result<f64, String>
where
    F: FnMut(f64) -> (f64, f64),
{
    let (mut lo, mut hi) = (1e-10, 1.0);
    while equation(hi).0 > 0.0 && hi < 1e10 {
        hi *= 2.0;
    }
    if equation(lo).0 <= 0.0 || equation(hi).0 >= 0.0 {
        return Err("finite identifiable shape root was not bracketed".into());
    }
    let mut x = (lo + hi) / 2.0;
    for _ in 0..160 {
        let (value, derivative) = equation(x);
        if !value.is_finite() || !derivative.is_finite() {
            return Err("shape likelihood is numerically ill conditioned".into());
        }
        if value > 0.0 {
            lo = x;
        } else {
            hi = x;
        }
        if hi - lo <= 2e-12 * x.max(1e-10) {
            return Ok((lo + hi) / 2.0);
        }
        let next = x - value / derivative;
        x = if next > lo && next < hi {
            next
        } else {
            (lo + hi) / 2.0
        };
    }
    Err("shape likelihood solver did not converge".into())
}

fn gamma_fit(x: &[f64]) -> Result<(f64, f64, f64), String> {
    let max = x.iter().copied().fold(0.0_f64, f64::max);
    let mean_relative = x.iter().map(|v| v / max).sum::<f64>() / x.len() as f64;
    let mean_log_relative = x.iter().map(|v| (v / max).ln()).sum::<f64>() / x.len() as f64;
    let s = mean_relative.ln() - mean_log_relative;
    if !s.is_finite() || s < 1e-12 {
        return Err("Gamma shape is not numerically identifiable".into());
    }
    let shape = root_decreasing(|a| {
        let (psi, tri) = digamma_trigamma(a);
        (a.ln() - psi - s, 1.0 / a - tri)
    })?;
    let scale = finite(max * mean_relative / shape)?;
    if scale <= 0.0 {
        return Err("Gamma scale is not representable".into());
    }
    let mean_log = max.ln() + mean_log_relative;
    let ll =
        x.len() as f64 * ((shape - 1.0) * mean_log - shape * scale.ln() - log_gamma(shape) - shape);
    Ok((shape, scale, finite(ll)?))
}

fn weibull_fit(x: &[f64], events: &[f64]) -> Result<(f64, f64, f64), String> {
    let logs: Vec<_> = x.iter().map(|v| v.ln()).collect();
    let max_log = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let z: Vec<_> = logs.iter().map(|v| v - max_log).collect();
    let d = events.iter().sum::<f64>();
    let event_mean = z.iter().zip(events).map(|(v, e)| v * e).sum::<f64>() / d;
    if !event_mean.is_finite() || event_mean >= 0.0 {
        return Err("finite Weibull shape is not identifiable".into());
    }
    let shape = root_decreasing(|k| {
        let sum = z.iter().map(|v| (k * v).exp()).sum::<f64>();
        let weighted_mean = z.iter().map(|v| v * (k * v).exp()).sum::<f64>() / sum;
        let weighted_var = z
            .iter()
            .map(|v| (v - weighted_mean).powi(2) * (k * v).exp())
            .sum::<f64>()
            / sum;
        (
            1.0 / k + event_mean - weighted_mean,
            -1.0 / (k * k) - weighted_var,
        )
    })?;
    let sum = z.iter().map(|v| (shape * v).exp()).sum::<f64>();
    let delta = (sum / d).ln() / shape;
    let log_scale = max_log + delta;
    let scale = finite(log_scale.exp())?;
    if scale <= 0.0 {
        return Err("Weibull scale is not representable".into());
    }
    let ll = d * (shape.ln() - log_scale + (shape - 1.0) * (event_mean - delta) - 1.0);
    Ok((shape, scale, finite(ll)?))
}

fn fit_distribution(operation: &str, x: &[f64]) -> Result<String, String> {
    let n = x.len();
    let nf = n as f64;
    match operation {
        "fit_normal" | "fit_lognormal" => {
            if n < 2 {
                return Err("at least two observations are required".into());
            }
            let values: Vec<_> = if operation == "fit_lognormal" {
                x.iter().map(|v| v.ln()).collect()
            } else {
                x.to_vec()
            };
            let (mean, sigma) = mean_sd(&values)?;
            let ll = -nf * (LOG_SQRT_TAU + sigma.ln() + 0.5);
            if operation == "fit_lognormal" {
                fitted(
                    operation,
                    n,
                    &[("log_mean", mean), ("log_sigma", sigma)],
                    ll - values.iter().sum::<f64>(),
                    "",
                )
            } else {
                fitted(operation, n, &[("mean", mean), ("sigma", sigma)], ll, "")
            }
        }
        "fit_exponential" => {
            let sum: f64 = x.iter().sum();
            if sum <= 0.0 {
                return Err("positive total exposure is required".into());
            }
            let rate = nf / sum;
            fitted(operation, n, &[("rate", rate)], nf * (rate.ln() - 1.0), "")
        }
        "fit_poisson" => {
            if x.iter()
                .any(|v| v.fract() != 0.0 || *v > 9_007_199_254_740_991.0)
            {
                return Err("Poisson observations must be exact nonnegative integer counts".into());
            }
            let rate = x.iter().sum::<f64>() / nf;
            let ll = if rate == 0.0 {
                0.0
            } else {
                x.iter().map(|v| poisson_log_mass(*v, rate)).sum()
            };
            fitted(operation, n, &[("rate", rate)], ll, "")
        }
        "fit_gamma" => {
            let (shape, scale, ll) = gamma_fit(x)?;
            fitted(operation, n, &[("shape", shape), ("scale", scale)], ll, "")
        }
        "fit_weibull" => {
            let events = vec![1.0; n];
            let (shape, scale, ll) = weibull_fit(x, &events)?;
            fitted(operation, n, &[("shape", shape), ("scale", scale)], ll, "")
        }
        _ => Err("unknown distribution operation".into()),
    }
}

fn normality(x: &[f64]) -> Result<String, String> {
    if x.len() < 3 {
        return Err("at least three observations are required".into());
    }
    let (m, s) = mean_sd(x)?;
    let nf = x.len() as f64;
    let skew = x.iter().map(|v| ((v - m) / s).powi(3)).sum::<f64>() / nf;
    let kurtosis = x.iter().map(|v| ((v - m) / s).powi(4)).sum::<f64>() / nf;
    let statistic = finite(nf / 6.0 * (skew * skew + (kurtosis - 3.0).powi(2) / 4.0))?;
    let p = (-statistic / 2.0).exp();
    Ok(format!("{{\"version\":1,\"operation\":\"normality_test\",\"method\":\"jarque_bera\",\"probability_method\":\"asymptotic_chi_square\",\"large_sample_recommended\":true,\"small_sample_warning\":{},\"n\":{},\"skewness\":{skew},\"kurtosis\":{kurtosis},\"statistic\":{statistic},\"degrees_of_freedom\":2,\"p_value\":{p}}}",x.len()<=2000,x.len()))
}

fn kde(x: &[f64], args: &[f64]) -> Result<String, String> {
    if args.len() < 2
        || args.len() > 4097
        || args[0] <= 0.0
        || x.len() * (args.len() - 1) > 20_000_000
    {
        return Err("KDE requires a positive bandwidth and bounded evaluation points/work".into());
    }
    let h = args[0];
    let factor = finite((1.0 / (2.0 * PI).sqrt()) / h)?;
    let density: Vec<_> = args[1..]
        .iter()
        .map(|at| {
            factor
                * (x.iter()
                    .map(|v| (-0.5 * ((at - v) / h).powi(2)).exp())
                    .sum::<f64>()
                    / x.len() as f64)
        })
        .collect();
    Ok(format!("{{\"version\":1,\"operation\":\"kde_gaussian\",\"bandwidth\":{h},\"evaluation_points\":{},\"densities\":{}}}",vec_json(&args[1..])?,vec_json(&density)?))
}

fn box_cox(x: &[f64], args: &[f64]) -> Result<String, String> {
    if args.len() > 1 || x.len() < 2 {
        return Err("Box-Cox expects at least two values and zero or one lambda".into());
    }
    let logs: Vec<_> = x.iter().map(|v| v.ln()).collect();
    let sum_log: f64 = logs.iter().sum();
    let center_log = sum_log / x.len() as f64;
    let transform = |lambda: f64| -> Vec<f64> {
        logs.iter()
            .map(|l| {
                if lambda == 0.0 {
                    *l
                } else {
                    (lambda * l).exp_m1() / lambda
                }
            })
            .collect()
    };
    let score = |lambda: f64| -> f64 {
        // Center in log units before evaluating variance. The omitted scale
        // exp(lambda*center_log) cancels the lambda-dependent Jacobian exactly.
        // This prevents changed measurement units from destroying the profile.
        let y: Vec<_> = logs
            .iter()
            .map(|l| {
                let z = l - center_log;
                if lambda == 0.0 {
                    z
                } else {
                    (lambda * z).exp_m1() / lambda
                }
            })
            .collect();
        if y.iter().any(|v| !v.is_finite()) {
            return f64::NEG_INFINITY;
        }
        match mean_sd(&y) {
            Ok((_, s)) => -sum_log - x.len() as f64 * s.ln(),
            _ => f64::NEG_INFINITY,
        }
    };
    let lambda = if let Some(lambda) = args.first() {
        *lambda
    } else {
        // Expand a bracket before golden-section maximization; never accept an
        // arbitrary search boundary as a claimed unconstrained optimum.
        let (mut left, mut right) = (-2.0, 2.0);
        for _ in 0..5 {
            let mid = (left + right) / 2.0;
            if score(left) < score(mid) && score(right) < score(mid) {
                break;
            }
            if score(left) > score(right) {
                right = mid;
                left *= 2.0;
            } else {
                left = mid;
                right *= 2.0;
            }
        }
        let (a, b) = (score(left), score(right));
        let mid = score((left + right) / 2.0);
        if !mid.is_finite() || a >= mid || b >= mid {
            return Err("Box-Cox finite optimum was not bracketed".into());
        }
        let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
        let mut c = right - ratio * (right - left);
        let mut d = left + ratio * (right - left);
        let mut fc = score(c);
        let mut fd = score(d);
        for _ in 0..100 {
            if fc > fd {
                right = d;
                d = c;
                fd = fc;
                c = right - ratio * (right - left);
                fc = score(c);
            } else {
                left = c;
                c = d;
                fc = fd;
                d = left + ratio * (right - left);
                fd = score(d);
            }
            if right - left < 1e-9 {
                break;
            }
        }
        (left + right) / 2.0
    };
    let ll = finite(score(lambda))?;
    Ok(format!("{{\"version\":1,\"operation\":\"box_cox\",\"lambda\":{lambda},\"lambda_estimated\":{},\"profile_log_likelihood\":{ll},\"likelihood_constant_omitted\":true,\"transformed\":{},\"converged\":true}}",args.is_empty(),vec_json(&transform(lambda))?))
}

fn survival(operation: &str, data: &[Vec<f64>]) -> Result<String, String> {
    if data
        .iter()
        .any(|r| r[0] < 0.0 || (r[1] != 0.0 && r[1] != 1.0))
    {
        return Err("survival rows require nonnegative times and binary event flags".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    let events: Vec<_> = data.iter().map(|r| r[1]).collect();
    let d = events.iter().sum::<f64>();
    if d == 0.0 {
        return Err("at least one observed event is required".into());
    }
    let extra = format!(",\"events\":{d},\"right_censored\":true");
    if operation == "survival_exponential" {
        let sum = x.iter().sum::<f64>();
        if sum <= 0.0 {
            return Err("positive total exposure is required".into());
        }
        let rate = d / sum;
        let se = rate / d.sqrt();
        finite(se)?;
        fitted(operation,x.len(),&[("rate",rate)],d*(rate.ln()-1.0),&format!("{extra},\"rate_standard_error\":{se},\"standard_error_method\":\"asymptotic_observed_information\""))
    } else {
        if x.iter().any(|v| *v <= 0.0) {
            return Err("Weibull times must be strictly positive".into());
        }
        let (shape, scale, ll) = weibull_fit(&x, &events)?;
        fitted(
            operation,
            x.len(),
            &[("shape", shape), ("scale", scale)],
            ll,
            &extra,
        )
    }
}

// Chi-square(1) upper tail Q(1/2,x/2), with direct continued fraction in
// the tail rather than subtracting a rounded CDF from one.
fn chi1_tail(statistic: f64) -> Result<f64, String> {
    let x = statistic / 2.0;
    if x == 0.0 {
        return Ok(1.0);
    }
    let factor = (-x + 0.5 * x.ln() - 0.5 * PI.ln()).exp();
    if x < 1.5 {
        let mut term = 2.0;
        let mut sum = term;
        for i in 1..1000 {
            term *= x / (0.5 + i as f64);
            sum += term;
            if term.abs() <= sum.abs() * 2e-15 {
                return Ok((1.0 - factor * sum).clamp(0.0, 1.0));
            }
        }
    } else {
        let tiny = 1e-300;
        let mut b = x + 0.5;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..1000 {
            let a = -(i as f64) * (i as f64 - 0.5);
            b += 2.0;
            d = a * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + a / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.0).abs() < 2e-15 {
                return Ok((factor * h).clamp(0.0, 1.0));
            }
        }
    }
    Err("chi-square tail failed to converge".into())
}

fn logrank(data: &[Vec<f64>]) -> Result<String, String> {
    if data
        .iter()
        .any(|r| r[0] < 0.0 || (r[1] != 0.0 && r[1] != 1.0) || (r[2] != 0.0 && r[2] != 1.0))
    {
        return Err("log-rank expects nonnegative time, binary event and binary group".into());
    }
    let mut rows = data.to_vec();
    rows.sort_by(|a, b| a[0].total_cmp(&b[0]));
    let mut risk = [
        rows.iter().filter(|r| r[2] == 0.0).count() as f64,
        rows.iter().filter(|r| r[2] == 1.0).count() as f64,
    ];
    if risk.contains(&0.0) {
        return Err("both comparison groups are required".into());
    }
    let (mut observed, mut expected, mut variance) = (0.0, 0.0, 0.0);
    let mut start = 0;
    while start < rows.len() {
        let mut end = start + 1;
        while end < rows.len() && rows[end][0] == rows[start][0] {
            end += 1;
        }
        let mut departures = [0.0, 0.0];
        let mut events = [0.0, 0.0];
        for row in &rows[start..end] {
            let g = row[2] as usize;
            departures[g] += 1.0;
            events[g] += row[1];
        }
        let n = risk[0] + risk[1];
        let d = events[0] + events[1];
        observed += events[0];
        expected += d * risk[0] / n;
        if n > 1.0 {
            variance += risk[0] * risk[1] * d * (n - d) / (n * n * (n - 1.0));
        }
        risk[0] -= departures[0];
        risk[1] -= departures[1];
        start = end;
    }
    if variance <= 0.0 {
        return Err("log-rank variance must be positive".into());
    }
    let statistic = finite((observed - expected).powi(2) / variance)?;
    let p = chi1_tail(statistic)?;
    Ok(format!("{{\"version\":1,\"operation\":\"logrank_test\",\"observed_group0\":{observed},\"expected_group0\":{expected},\"variance\":{variance},\"statistic\":{statistic},\"degrees_of_freedom\":1,\"p_value\":{p},\"probability_method\":\"asymptotic_chi_square\",\"ties\":\"hypergeometric\"}}"))
}

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let cols = match operation {
        "survival_exponential" | "survival_weibull" => 2,
        "logrank_test" => 3,
        _ => 1,
    };
    check(data, cols, args)?;
    if !matches!(operation, "box_cox" | "kde_gaussian") && !args.is_empty() {
        return Err("unexpected distribution parameters".into());
    }
    if matches!(operation, "survival_exponential" | "survival_weibull") {
        return survival(operation, data);
    }
    if operation == "logrank_test" {
        return logrank(data);
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    if matches!(
        operation,
        "fit_gamma" | "fit_weibull" | "fit_lognormal" | "box_cox"
    ) && x.iter().any(|v| *v <= 0.0)
    {
        return Err("strictly positive observations required".into());
    }
    if matches!(operation, "fit_exponential" | "fit_poisson") && x.iter().any(|v| *v < 0.0) {
        return Err("nonnegative observations required".into());
    }
    match operation {
        "normality_test" => normality(&x),
        "kde_gaussian" => kde(&x, args),
        "box_cox" => box_cox(&x, args),
        _ => fit_distribution(operation, &x),
    }
}
