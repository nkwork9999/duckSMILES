//! Univariate, equally spaced time-series analysis in the caller's row order.
//!
//! This module provides diagnostics, fixed-parameter smoothers, conditional
//! least-squares AR, a direct-DFT periodogram and additive seasonal decomposition.
//! It does not fit MA, ARIMA/SARIMA, automatic smoothing parameters or prediction
//! intervals. Nulls only occur at unidentified decomposition edges.

const MAX_ROWS: usize = 100_000;
const MAX_WORK: usize = 50_000_000;
const MAX_HORIZON: usize = 10_000;

/// All inputs are finite, one-column numeric rows in increasing time order.
/// Outputs use zero-based original input indices; future index n is horizon 1.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data.len() > MAX_ROWS
        || data.iter().any(|r| r.len() != 1 || !r[0].is_finite())
        || args.iter().any(|x| !x.is_finite())
    {
        return Err(
            "time series requires 1 to 100000 finite one-column rows and finite arguments".into(),
        );
    }
    let x = data.iter().map(|r| r[0]).collect::<Vec<_>>();
    let n = x.len();
    let output = match operation {
        "acf" | "pacf" => {
            arguments(args, 1, 1)?;
            let lag = integer(args[0], 0, n.saturating_sub(1).min(10_000), "maximum lag")?;
            if operation == "pacf" && lag > 2048 {
                return Err("PACF is limited to 2048 lags".into());
            }
            let r = acf(&x, lag)?;
            let values = if operation == "pacf" { pacf(&r)? } else { r };
            let (name, key) = if operation == "pacf" {
                (
                    "durbin_levinson_biased_autocovariance",
                    "partial_autocorrelations",
                )
            } else {
                ("biased_sample_autocorrelation", "autocorrelations")
            };
            Ok(format!("{{\"version\":1,\"algorithm\":\"{name}\",\"n_samples\":{n},\"max_lag\":{lag},\"lags\":{},\"{key}\":{}}}",indices(0,lag+1),numbers(&values)))
        }
        "difference" => {
            arguments(args, 0, 2)?;
            let lag = integer(
                args.first().copied().unwrap_or(1.),
                1,
                n.saturating_sub(1),
                "lag",
            )?;
            let order = integer(
                args.get(1).copied().unwrap_or(1.),
                1,
                32,
                "difference order",
            )?;
            if lag.checked_mul(order).is_none_or(|k| k >= n) {
                return Err("differencing must leave at least one observation".into());
            }
            let mut values = x;
            for _ in 0..order {
                values = values[lag..]
                    .iter()
                    .zip(&values)
                    .map(|(a, b)| finite(a - b, "differenced value"))
                    .collect::<Result<Vec<_>, _>>()?;
            }
            Ok(format!("{{\"version\":1,\"algorithm\":\"iterated_lag_difference\",\"lag\":{lag},\"order\":{order},\"original_indices\":{},\"values\":{}}}",indices(lag*order,n),numbers(&values)))
        }
        "moving_average" => {
            arguments(args, 1, 1)?;
            let window = integer(args[0], 1, n, "moving-average window")?;
            bounded_work(n, window, 1)?;
            let values = x.windows(window).map(mean).collect::<Result<Vec<_>, _>>()?;
            Ok(format!("{{\"version\":1,\"algorithm\":\"trailing_simple_moving_average\",\"window\":{window},\"original_indices\":{},\"values\":{}}}",indices(window-1,n),numbers(&values)))
        }
        "exponential_smoothing" => {
            arguments(args, 2, 2)?;
            unit(args[0], "alpha")?;
            let horizon = integer(args[1], 1, MAX_HORIZON, "forecast horizon")?;
            smooth(&x, args[0], None, horizon)
        }
        "holt" => {
            arguments(args, 3, 4)?;
            unit(args[0], "alpha")?;
            unit(args[1], "beta")?;
            let horizon = integer(args[2], 1, MAX_HORIZON, "forecast horizon")?;
            let phi = args.get(3).copied().unwrap_or(1.);
            if !(0. ..=1.).contains(&phi) || phi == 0. {
                return Err("damping phi must be in (0,1]".into());
            }
            if n < 2 {
                return Err("Holt initialization requires at least two observations".into());
            }
            smooth(&x, args[0], Some((args[1], phi)), horizon)
        }
        "ar_model" => {
            arguments(args, 2, 2)?;
            let order = integer(args[0], 1, 32, "AR order")?;
            let horizon = integer(args[1], 1, MAX_HORIZON, "forecast horizon")?;
            ar(&x, order, horizon)
        }
        "periodogram" => {
            arguments(args, 0, 1)?;
            let dt = args.first().copied().unwrap_or(1.);
            if dt <= 0. || !(2..=4096).contains(&n) {
                return Err(
                    "periodogram requires positive sample interval and 2 to 4096 observations"
                        .into(),
                );
            }
            periodogram(&x, dt)
        }
        "seasonal_decompose" => {
            arguments(args, 1, 1)?;
            let period = integer(args[0], 2, n / 2, "seasonal period")?;
            decompose(&x, period)
        }
        "ljung_box" => {
            arguments(args, 1, 2)?;
            let lag = integer(args[0], 1, n.saturating_sub(1).min(10_000), "maximum lag")?;
            let model_df = integer(
                args.get(1).copied().unwrap_or(0.),
                0,
                lag - 1,
                "model degrees of freedom",
            )?;
            let r = acf(&x, lag)?;
            let nf = n as f64;
            let q = finite(
                nf * (nf + 2.) * (1..=lag).map(|k| r[k] * r[k] / (n - k) as f64).sum::<f64>(),
                "Ljung-Box Q",
            )?;
            let df = (lag - model_df) as f64;
            let p = gamma_sf(df / 2., q / 2.)?;
            Ok(format!("{{\"version\":1,\"algorithm\":\"ljung_box_chi_square_approximation\",\"n_samples\":{n},\"max_lag\":{lag},\"model_degrees_of_freedom\":{model_df},\"q_statistic\":{q},\"degrees_of_freedom\":{df},\"p_value\":{p}}}"))
        }
        _ => Err(format!("unknown time series operation: {operation}")),
    }?;
    if output.contains(":NaN")
        || output.contains(":inf")
        || output.contains(":-inf")
        || output.contains("[NaN")
        || output.contains("[inf")
        || output.contains("[-inf")
        || output.contains(",NaN")
        || output.contains(",inf")
        || output.contains(",-inf")
    {
        return Err("time series result exceeds finite numerical range".into());
    }
    Ok(output)
}
fn arguments(args: &[f64], min: usize, max: usize) -> Result<(), String> {
    if (min..=max).contains(&args.len()) {
        Ok(())
    } else {
        Err(format!("operation requires {min} to {max} arguments"))
    }
}
fn integer(x: f64, min: usize, max: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < min as f64 || x > max as f64 {
        Err(format!("{name} must be an integer between {min} and {max}"))
    } else {
        Ok(x as usize)
    }
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn unit(x: f64, name: &str) -> Result<(), String> {
    if (0. ..=1.).contains(&x) {
        Ok(())
    } else {
        Err(format!("{name} must be in [0,1]"))
    }
}
fn bounded_work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.checked_mul(b)
        .and_then(|x| x.checked_mul(c))
        .is_none_or(|w| w > MAX_WORK)
    {
        Err("requested computation exceeds 50000000 work units".into())
    } else {
        Ok(())
    }
}
fn numbers(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn optional_numbers(x: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        x.iter()
            .map(|x| x.map_or_else(|| "null".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn indices(start: usize, end: usize) -> String {
    numbers(&(start..end).map(|i| i as f64).collect::<Vec<_>>())
}
fn mean_offset(x: &[f64]) -> Result<f64, String> {
    let mut m = 0.;
    for (i, &v) in x.iter().enumerate() {
        let d = finite(v - x[0], "centered observation")?;
        m = finite(m + (d - m) / (i + 1) as f64, "sample mean offset")?;
    }
    Ok(m)
}
fn mean(x: &[f64]) -> Result<f64, String> {
    finite(x[0] + mean_offset(x)?, "sample mean")
}
fn centered(x: &[f64]) -> Result<(f64, Vec<f64>), String> {
    let offset = mean_offset(x)?;
    let m = finite(x[0] + offset, "sample mean")?;
    Ok((
        m,
        x.iter()
            .map(|v| finite((v - x[0]) - offset, "centered observation"))
            .collect::<Result<Vec<_>, _>>()?,
    ))
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(a.iter().zip(b).map(|(x, y)| x * y).sum(), "inner product")
}
fn acf(x: &[f64], lag: usize) -> Result<Vec<f64>, String> {
    if x.len() < 2 {
        return Err("correlation requires at least two observations".into());
    }
    bounded_work(x.len(), lag.max(1), 1)?;
    let (_, mut y) = centered(x)?;
    let scale = y.iter().map(|v| v.abs()).fold(0., f64::max);
    if scale == 0. {
        return Err("autocorrelation undefined for a constant series".into());
    }
    // Normalize before forming products to avoid overflow and tiny-square loss.
    for v in &mut y {
        *v /= scale;
    }
    let denominator = dot(&y, &y)?;
    let mut r = vec![1.];
    for k in 1..=lag {
        r.push(dot(&y[k..], &y[..y.len() - k])? / denominator);
    }
    Ok(r)
}
fn pacf(r: &[f64]) -> Result<Vec<f64>, String> {
    bounded_work(r.len(), r.len(), 1)?;
    let mut result = vec![1.];
    let mut phi = Vec::<f64>::new();
    let mut error = 1.;
    for k in 1..r.len() {
        if error <= 1e-14 {
            return Err("PACF recursion has numerically singular covariance".into());
        }
        let reflection = finite(
            (r[k] - (1..k).map(|j| phi[j - 1] * r[k - j]).sum::<f64>()) / error,
            "partial autocorrelation",
        )?;
        if reflection.abs() > 1. + 1e-10 {
            return Err("PACF covariance lost positive definiteness".into());
        }
        let mut next = vec![0.; k];
        for j in 1..k {
            next[j - 1] = phi[j - 1] - reflection * phi[k - j - 1];
        }
        next[k - 1] = reflection;
        phi = next;
        error *= 1. - reflection * reflection;
        result.push(reflection);
    }
    Ok(result)
}
fn smooth(
    x: &[f64],
    alpha: f64,
    trend_parameters: Option<(f64, f64)>,
    horizon: usize,
) -> Result<String, String> {
    let mut level = x[0];
    let mut trend = if trend_parameters.is_some() {
        finite(x[1] - x[0], "initial trend")?
    } else {
        0.
    };
    let mut levels = vec![level];
    let mut trends = vec![trend];
    let mut fitted = Vec::new();
    let mut residuals = Vec::new();
    for &v in &x[1..] {
        let (beta, phi) = trend_parameters.unwrap_or((0., 0.));
        let prediction = finite(level + phi * trend, "one-step prediction")?;
        let residual = finite(v - prediction, "smoothing residual")?;
        let old = level;
        level = finite(prediction + alpha * residual, "smoothed level")?;
        if trend_parameters.is_some() {
            trend = finite(
                beta * (level - old) + (1. - beta) * phi * trend,
                "smoothed trend",
            )?;
        }
        fitted.push(prediction);
        residuals.push(residual);
        levels.push(level);
        trends.push(trend);
    }
    let mut forecasts = Vec::with_capacity(horizon);
    let mut cumulative = 0.;
    let mut damp = 1.;
    for _ in 0..horizon {
        damp *= trend_parameters.map_or(0., |p| p.1);
        cumulative += damp;
        forecasts.push(finite(level + cumulative * trend, "smoothing forecast")?);
    }
    let name = if trend_parameters.is_some() {
        "holt_fixed_parameters_first_difference_initialization"
    } else {
        "simple_exponential_smoothing_first_observation_initialization"
    };
    let extra = trend_parameters.map_or_else(String::new, |(beta, phi)| {
        format!(
            ",\"beta\":{beta},\"phi\":{phi},\"trends\":{},\"final_trend\":{trend}",
            numbers(&trends)
        )
    });
    let sse = dot(&residuals, &residuals)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"{name}\",\"alpha\":{alpha}{extra},\"levels\":{},\"final_level\":{level},\"fitted_indices\":{},\"fitted_values\":{},\"residuals\":{},\"sum_squared_errors\":{sse},\"forecast_indices\":{},\"forecasts\":{}}}",numbers(&levels),indices(1,x.len()),numbers(&fitted),numbers(&residuals),indices(x.len(),x.len()+horizon),numbers(&forecasts)))
}

fn ar(x: &[f64], p: usize, horizon: usize) -> Result<String, String> {
    if x.len() <= 2 * p + 1 {
        return Err(
            "conditional AR fit requires n > 2*order + 1 for residual degrees of freedom".into(),
        );
    }
    let n = x.len() - p;
    bounded_work(n, p, p)?;
    let (ymean, y) = centered(&x[p..])?;
    let origin = x[0];
    let response_offset = finite(
        (x[p] - origin) + mean_offset(&x[p..])?,
        "AR response mean offset",
    )?;
    let mut means = Vec::new();
    let mut mean_offsets = Vec::new();
    let mut scales = Vec::new();
    let mut design = Vec::new();
    for lag in 1..=p {
        let (m, mut column) = centered(&x[p - lag..x.len() - lag])?;
        mean_offsets.push(finite(
            (x[p - lag] - origin) + mean_offset(&x[p - lag..x.len() - lag])?,
            "AR lag mean offset",
        )?);
        let scale = column.iter().map(|x| x.abs()).fold(0., f64::max);
        if scale == 0. {
            return Err("AR lag matrix is rank deficient".into());
        }
        for x in &mut column {
            *x /= scale;
        }
        means.push(m);
        scales.push(scale);
        design.push(column);
    }
    let mut q = Vec::<Vec<f64>>::new();
    let mut r = vec![vec![0.; p]; p];
    for j in 0..p {
        let mut v = design[j].clone();
        let original = dot(&v, &v)?.sqrt();
        // Twice-reorthogonalized modified Gram-Schmidt on centered/scaled lags.
        for _ in 0..2 {
            for i in 0..j {
                let amount = dot(&q[i], &v)?;
                r[i][j] += amount;
                for (v, q) in v.iter_mut().zip(&q[i]) {
                    *v -= amount * q;
                }
            }
        }
        let norm = dot(&v, &v)?.sqrt();
        if norm <= 1e-10 * original.max(1.) {
            return Err("AR lag matrix is numerically rank deficient".into());
        }
        r[j][j] = norm;
        for v in &mut v {
            *v /= norm;
        }
        q.push(v);
    }
    let mut scaled_coeff = q
        .iter()
        .map(|v| dot(v, &y))
        .collect::<Result<Vec<_>, _>>()?;
    for i in (0..p).rev() {
        scaled_coeff[i] = (scaled_coeff[i]
            - (i + 1..p).map(|j| r[i][j] * scaled_coeff[j]).sum::<f64>())
            / r[i][i];
        finite(scaled_coeff[i], "AR coefficient")?;
    }
    let coefficients = scaled_coeff
        .iter()
        .zip(&scales)
        .map(|(b, s)| finite(b / s, "AR coefficient"))
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        origin * (1. - coefficients.iter().sum::<f64>()) + response_offset
            - coefficients
                .iter()
                .zip(&mean_offsets)
                .map(|(b, m)| b * m)
                .sum::<f64>(),
        "AR intercept",
    )?;
    let fitted_centered = (0..n)
        .map(|i| {
            finite(
                (0..p).map(|j| scaled_coeff[j] * design[j][i]).sum(),
                "AR fitted value",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let residuals = y
        .iter()
        .zip(&fitted_centered)
        .map(|(y, f)| finite(y - f, "AR residual"))
        .collect::<Result<Vec<_>, _>>()?;
    let sse = dot(&residuals, &residuals)?;
    // A degenerate Gaussian variance has no finite likelihood/AIC/BIC.
    if sse <= 0. {
        return Err(
            "AR conditional innovation variance must be positive for finite Gaussian likelihood"
                .into(),
        );
    }
    let variance = sse / n as f64;
    if variance == 0. {
        return Err("AR innovation variance underflow".into());
    }
    let unbiased = sse / (n - p - 1) as f64;
    let log_likelihood = finite(
        -0.5 * n as f64 * ((2. * std::f64::consts::PI).ln() + 1. + variance.ln()),
        "conditional Gaussian likelihood",
    )?;
    let parameter_count = p + 2;
    let aic = -2. * log_likelihood + 2. * parameter_count as f64;
    let bic = -2. * log_likelihood + (n as f64).ln() * parameter_count as f64;
    let fitted = fitted_centered
        .iter()
        .map(|f| finite(origin + (response_offset + f), "AR fitted value"))
        .collect::<Result<Vec<_>, _>>()?;
    // Keep recursive state relative to an observed origin. Rounding every
    // future point back into original units would accumulate errors when the
    // origin is large compared with the fluctuations in the series.
    let mut history = x
        .iter()
        .map(|v| finite(v - origin, "AR relative history"))
        .collect::<Result<Vec<_>, _>>()?;
    let mut forecasts = Vec::new();
    for _ in 0..horizon {
        let f = finite(
            response_offset
                + (0..p)
                    .map(|j| coefficients[j] * (history[history.len() - j - 1] - mean_offsets[j]))
                    .sum::<f64>(),
            "AR recursive forecast",
        )?;
        forecasts.push(finite(origin + f, "AR forecast in original units")?);
        history.push(f);
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"conditional_least_squares_ar_qr\",\"order\":{p},\"n_effective\":{n},\"intercept\":{intercept},\"coefficients\":{},\"lag_means\":{},\"response_mean\":{ymean},\"reference_origin\":{origin},\"response_mean_offset\":{response_offset},\"lag_mean_offsets\":{},\"fitted_indices\":{},\"fitted_values\":{},\"residuals\":{},\"innovation_variance\":{variance},\"residual_variance_unbiased\":{unbiased},\"likelihood\":\"conditional_gaussian_mle_variance\",\"log_likelihood\":{log_likelihood},\"parameter_count\":{parameter_count},\"aic\":{aic},\"bic\":{bic},\"forecast_indices\":{},\"forecasts\":{}}}",numbers(&coefficients),numbers(&means),numbers(&mean_offsets),indices(p,x.len()),numbers(&fitted),numbers(&residuals),indices(x.len(),x.len()+horizon),numbers(&forecasts)))
}
fn periodogram(x: &[f64], dt: f64) -> Result<String, String> {
    let n = x.len();
    let (_, y) = centered(x)?;
    let nf = n as f64;
    let duration = finite(nf * dt, "sample duration")?;
    if duration == 0. {
        return Err("sample duration underflow".into());
    }
    let mut frequencies = Vec::new();
    let mut power = Vec::new();
    for k in 0..=n / 2 {
        let mut real = 0.;
        let mut imaginary = 0.;
        for (t, v) in y.iter().enumerate() {
            let angle = 2. * std::f64::consts::PI * k as f64 * t as f64 / nf;
            real += v * angle.cos();
            imaginary -= v * angle.sin();
        }
        let factor = if k == 0 || n & 1 == 0 && k == n / 2 {
            1.
        } else {
            2.
        };
        frequencies.push(finite(k as f64 / duration, "periodogram frequency")?);
        power.push(finite(
            (real * real + imaginary * imaginary) * (dt / nf) * factor,
            "spectral density",
        )?);
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"direct_dft_one_sided_periodogram\",\"detrend\":\"constant\",\"scaling\":\"density\",\"n_samples\":{n},\"sample_interval\":{dt},\"frequencies\":{},\"power_spectral_density\":{}}}",numbers(&frequencies),numbers(&power)))
}
fn decompose(x: &[f64], period: usize) -> Result<String, String> {
    let n = x.len();
    let radius = period / 2;
    let width = 2 * radius + 1;
    bounded_work(n, width, 1)?;
    let mut trend = vec![None; n];
    let mut sums = vec![0.; period];
    let mut counts = vec![0usize; period];
    for i in radius..n - radius {
        let window = &x[i - radius..=i + radius];
        let origin = window[0];
        let mut offset = 0.;
        for (j, &v) in window.iter().enumerate() {
            let w = if period & 1 == 0 && (j == 0 || j == width - 1) {
                0.5
            } else {
                1.
            };
            offset = finite(
                offset + w * finite(v - origin, "trend centered value")? / period as f64,
                "moving trend",
            )?;
        }
        let t = finite(origin + offset, "moving trend")?;
        let d = finite((x[i] - origin) - offset, "detrended observation")?;
        trend[i] = Some(t);
        sums[i % period] = finite(sums[i % period] + d, "seasonal sum")?;
        counts[i % period] += 1;
    }
    if counts.contains(&0) {
        return Err("seasonal phase has no identified detrended observations".into());
    }
    let mut pattern = sums
        .iter()
        .zip(counts)
        .map(|(s, n)| s / n as f64)
        .collect::<Vec<_>>();
    let average = mean(&pattern)?;
    for v in &mut pattern {
        *v = finite(*v - average, "centered seasonal pattern")?;
    }
    let seasonal = (0..n).map(|i| pattern[i % period]).collect::<Vec<_>>();
    let residuals = (0..n)
        .map(|i| {
            trend[i]
                .map(|t| finite((x[i] - t) - seasonal[i], "seasonal residual"))
                .transpose()
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"classical_additive_centered_moving_average\",\"period\":{period},\"original_indices\":{},\"trend\":{},\"seasonal_pattern\":{},\"seasonal\":{},\"residuals\":{},\"trend_extrapolation\":false}}",indices(0,n),optional_numbers(&trend),numbers(&pattern),numbers(&seasonal),optional_numbers(&residuals)))
}
// Regularized upper incomplete gamma for the Ljung-Box chi-square tail.
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
