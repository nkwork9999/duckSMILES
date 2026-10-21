//! Small univariate multiplicative seasonal ARIMA models fitted by conditional
//! sum of squared innovations, not exact state-space maximum likelihood.
//! A(B)(w_t-mu)=C(B)e_t on the differenced series, with zero presample deviations
//! and innovations. CSS discards max(deg(A),deg(C)) initial innovations.

use std::f64::consts::PI;

struct Spec {
    p: usize,
    d: usize,
    q: usize,
    horizon: usize,
    constant: bool,
    iterations: usize,
    tolerance: f64,
    sp: usize,
    sd: usize,
    sq: usize,
    season: usize,
}

fn integer(value: f64, max: usize) -> Result<usize, String> {
    if value.is_finite() && value >= 0.0 && value.fract() == 0.0 && value <= max as f64 {
        Ok(value as usize)
    } else {
        Err("ARIMA orders and limits must be bounded nonnegative integers".into())
    }
}

fn specification(args: &[f64]) -> Result<Spec, String> {
    if !(4..=11).contains(&args.len()) || args.iter().any(|v| !v.is_finite()) {
        return Err(
            "ARIMA expects p,d,q,horizon and optional constant,iterations,tolerance,P,D,Q,season"
                .into(),
        );
    }
    let get = |at: usize, default: f64| args.get(at).copied().unwrap_or(default);
    let spec = Spec {
        p: integer(args[0], 3)?,
        d: integer(args[1], 2)?,
        q: integer(args[2], 3)?,
        horizon: integer(args[3], 1000)?,
        constant: integer(get(4, 1.0), 1)? == 1,
        iterations: integer(get(5, 1000.0), 10_000)?,
        tolerance: get(6, 1e-8),
        sp: integer(get(7, 0.0), 1)?,
        sd: integer(get(8, 0.0), 2)?,
        sq: integer(get(9, 0.0), 1)?,
        season: integer(get(10, 0.0), 52)?,
    };
    if !(1e-12..=1e-3).contains(&spec.tolerance)
        || ((spec.sp + spec.sd + spec.sq) > 0 && spec.season < 2)
        || spec.season == 1
    {
        return Err("invalid seasonal period or optimization tolerance".into());
    }
    Ok(spec)
}

fn finite(value: f64) -> Result<f64, String> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err("ARIMA numerical result is not representable".into())
    }
}

fn json_array(values: &[f64]) -> Result<String, String> {
    values
        .iter()
        .map(|v| finite(*v).map(|v| v.to_string()))
        .collect::<Result<Vec<_>, _>>()
        .map(|v| format!("[{}]", v.join(",")))
}

// The Levinson recursion maps open-unit reflection coefficients to a stationary
// AR polynomial. Negating AR coefficients gives an invertible MA polynomial.
fn stationary(parameters: &[f64]) -> Vec<f64> {
    let mut phi: Vec<f64> = Vec::new();
    for value in parameters {
        let kappa = value.tanh() * (1.0 - 1e-6);
        let mut next = vec![0.0; phi.len() + 1];
        for j in 0..phi.len() {
            next[j] = phi[j] - kappa * phi[phi.len() - j - 1];
        }
        next[phi.len()] = kappa;
        phi = next;
    }
    phi
}

fn multiply(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut out = vec![0.0; a.len() + b.len() - 1];
    for (i, av) in a.iter().enumerate() {
        for (j, bv) in b.iter().enumerate() {
            out[i + j] += av * bv;
        }
    }
    out
}

struct Model {
    ar: Vec<f64>,
    ma: Vec<f64>,
    seasonal_ar: Vec<f64>,
    seasonal_ma: Vec<f64>,
    a: Vec<f64>,
    c: Vec<f64>,
}

fn polynomial(coefficients: &[f64], stride: usize, sign: f64) -> Vec<f64> {
    let mut values = vec![0.0; coefficients.len() * stride + 1];
    values[0] = 1.0;
    for (i, value) in coefficients.iter().enumerate() {
        values[(i + 1) * stride] = sign * value;
    }
    values
}

fn model(parameters: &[f64], s: &Spec) -> Model {
    let ar = stationary(&parameters[..s.p]);
    let ma: Vec<_> = stationary(&parameters[s.p..s.p + s.q])
        .iter()
        .map(|v| -v)
        .collect();
    let seasonal_ar = stationary(&parameters[s.p + s.q..s.p + s.q + s.sp]);
    let seasonal_ma: Vec<_> = stationary(&parameters[s.p + s.q + s.sp..])
        .iter()
        .map(|v| -v)
        .collect();
    let a = multiply(
        &polynomial(&ar, 1, -1.0),
        &polynomial(&seasonal_ar, s.season.max(1), -1.0),
    );
    let c = multiply(
        &polynomial(&ma, 1, 1.0),
        &polynomial(&seasonal_ma, s.season.max(1), 1.0),
    );
    Model {
        ar,
        ma,
        seasonal_ar,
        seasonal_ma,
        a,
        c,
    }
}

struct Evaluation {
    value: f64,
    mean: f64,
    innovations: Vec<f64>,
    sse: f64,
}

fn evaluate(y: &[f64], m: &Model, burn: usize, constant: bool) -> Option<Evaluation> {
    let n = y.len();
    let mut ey = vec![0.0; n];
    let mut unit = vec![0.0; n];
    for i in 0..n {
        ey[i] = y[i];
        unit[i] = 1.0;
        for lag in 1..m.a.len().min(i + 1) {
            ey[i] += m.a[lag] * y[i - lag];
            unit[i] += m.a[lag];
        }
        for lag in 1..m.c.len().min(i + 1) {
            ey[i] -= m.c[lag] * ey[i - lag];
            unit[i] -= m.c[lag] * unit[i - lag];
        }
    }
    let mean = if constant {
        let denominator = unit[burn..].iter().map(|v| v * v).sum::<f64>();
        if denominator < 1e-30 {
            return None;
        }
        ey[burn..]
            .iter()
            .zip(&unit[burn..])
            .map(|(e, u)| e * u)
            .sum::<f64>()
            / denominator
    } else {
        0.0
    };
    let innovations: Vec<_> = ey.iter().zip(unit).map(|(e, u)| e - mean * u).collect();
    let sse = innovations[burn..].iter().map(|v| v * v).sum::<f64>();
    if !sse.is_finite() || sse <= 0.0 || !mean.is_finite() {
        return None;
    }
    Some(Evaluation {
        value: (sse / (n - burn) as f64).ln(),
        mean,
        innovations,
        sse,
    })
}

struct Optimized {
    parameters: Vec<f64>,
    converged: bool,
    iterations: usize,
}

fn optimize<F>(initial: Vec<f64>, max_iterations: usize, tolerance: f64, objective: F) -> Optimized
where
    F: Fn(&[f64]) -> f64,
{
    let size = initial.len();
    if size == 0 {
        return Optimized {
            parameters: initial,
            converged: true,
            iterations: 0,
        };
    }
    let mut simplex = Vec::with_capacity(size + 1);
    simplex.push((objective(&initial), initial.clone()));
    for i in 0..size {
        let mut point = initial.clone();
        point[i] += 0.2;
        simplex.push((objective(&point), point));
    }
    let mut used = 0;
    let mut converged = false;
    while used < max_iterations {
        simplex.sort_by(|a, b| a.0.total_cmp(&b.0));
        let spread = simplex
            .iter()
            .skip(1)
            .flat_map(|(_, x)| x.iter().zip(&simplex[0].1).map(|(a, b)| (a - b).abs()))
            .fold(0.0_f64, f64::max);
        if spread <= tolerance
            && (simplex[size].0 - simplex[0].0).abs() <= tolerance * (1.0 + simplex[0].0.abs())
        {
            converged = true;
            break;
        }
        let mut centroid = vec![0.0; size];
        for (_, point) in &simplex[..size] {
            for (center, v) in centroid.iter_mut().zip(point) {
                *center += v / size as f64;
            }
        }
        let affine = |from: &[f64], factor: f64| -> Vec<f64> {
            centroid
                .iter()
                .zip(from)
                .map(|(c, v)| c + factor * (v - c))
                .collect()
        };
        let reflected = affine(&simplex[size].1, -1.0);
        let fr = objective(&reflected);
        if fr < simplex[0].0 {
            let expanded = affine(&reflected, 2.0);
            let fe = objective(&expanded);
            simplex[size] = if fe < fr {
                (fe, expanded)
            } else {
                (fr, reflected)
            };
        } else if fr < simplex[size - 1].0 {
            simplex[size] = (fr, reflected);
        } else {
            let contracted = if fr < simplex[size].0 {
                affine(&reflected, 0.5)
            } else {
                affine(&simplex[size].1, 0.5)
            };
            let fc = objective(&contracted);
            if fc < fr.min(simplex[size].0) {
                simplex[size] = (fc, contracted);
            } else {
                let best = simplex[0].1.clone();
                for item in simplex.iter_mut().skip(1) {
                    for (x, b) in item.1.iter_mut().zip(&best) {
                        *x = (*x + b) / 2.0;
                    }
                    item.0 = objective(&item.1);
                }
            }
        }
        used += 1;
    }
    simplex.sort_by(|a, b| a.0.total_cmp(&b.0));
    Optimized {
        parameters: simplex.remove(0).1,
        converged,
        iterations: used,
    }
}

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if operation != "arima" {
        return Err("unknown ARIMA operation".into());
    }
    let s = specification(args)?;
    if !(3..=100_000).contains(&data.len())
        || data
            .iter()
            .any(|r| r.len() != 1 || !r[0].is_finite() || r[0].abs() > 1e150)
    {
        return Err("ARIMA expects ordered finite one-column values, 3..100000 rows".into());
    }
    let dim = s.p + s.q + s.sp + s.sq;
    let burn = (s.p + s.sp * s.season).max(s.q + s.sq * s.season);
    let lost = s.d + s.sd * s.season;
    if data.len() <= lost + burn + dim + usize::from(s.constant) + 1 {
        return Err("insufficient observations after differencing and CSS conditioning".into());
    }
    let work = (data.len() as u128)
        * (burn + 1) as u128
        * if dim == 0 {
            1
        } else {
            s.iterations.max(1) as u128 * (dim + 3) as u128
        };
    if work > 500_000_000 {
        return Err("ARIMA requested optimization work exceeds its limit".into());
    }
    let lags: Vec<usize> = std::iter::repeat_n(1, s.d)
        .chain(std::iter::repeat_n(s.season, s.sd))
        .collect();
    let mut layers = vec![data.iter().map(|r| r[0]).collect::<Vec<f64>>()];
    for lag in &lags {
        let prev = layers.last().unwrap();
        let diff: Vec<_> = (*lag..prev.len())
            .map(|i| prev[i] - prev[i - lag])
            .collect();
        if diff.iter().any(|v| !v.is_finite()) {
            return Err("differencing overflow".into());
        }
        layers.push(diff);
    }
    let differenced = layers.last().unwrap();
    let normalization_origin = if s.constant { differenced[0] } else { 0.0 };
    let scale = differenced
        .iter()
        .map(|v| (v - normalization_origin).abs())
        .fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return Err("positive innovation variance is required".into());
    }
    let mut y: Vec<_> = differenced
        .iter()
        .map(|v| (v - normalization_origin) / scale)
        .collect();
    let mut initial = vec![0.0; dim];
    if s.p > 0 {
        let mean = if s.constant {
            y.iter().sum::<f64>() / y.len() as f64
        } else {
            0.0
        };
        let denominator = y[..y.len() - 1]
            .iter()
            .map(|v| (v - mean).powi(2))
            .sum::<f64>();
        if denominator > 0.0 {
            let lag = y[1..]
                .iter()
                .zip(&y)
                .map(|(a, b)| (a - mean) * (b - mean))
                .sum::<f64>()
                / denominator;
            initial[0] = lag.clamp(-0.8, 0.8).atanh();
        }
    }
    let optimized = optimize(initial, s.iterations, s.tolerance, |parameters| {
        if parameters.iter().any(|v| !v.is_finite() || v.abs() > 20.0) {
            return f64::INFINITY;
        }
        let m = model(parameters, &s);
        evaluate(&y, &m, burn, s.constant)
            .map(|e| e.value)
            .unwrap_or(f64::INFINITY)
    });
    let m = model(&optimized.parameters, &s);
    let evaluation = evaluate(&y, &m, burn, s.constant)
        .ok_or("positive finite identifiable conditional innovation variance is required")?;
    let effective_n = y.len() - burn;
    // A deterministic conditional relation has no finite positive-variance
    // Gaussian MLE. Reject a residual below rounding resolution instead of
    // optimizing against an artificial positive SSE floor.
    if evaluation.sse / effective_n as f64 <= 128.0 * f64::EPSILON.powi(2) {
        return Err("innovation variance is below numerical identification precision".into());
    }
    let sigma_squared = finite(evaluation.sse / effective_n as f64 * scale * scale)?;
    if sigma_squared <= 0.0 {
        return Err("innovation variance is not representable".into());
    }
    let log_likelihood =
        finite(-0.5 * effective_n as f64 * ((2.0 * PI).ln() + sigma_squared.ln() + 1.0))?;
    let k = (dim + usize::from(s.constant) + 1) as f64;
    let aic = finite(2.0 * k - 2.0 * log_likelihood)?;
    let bic = finite(k * (effective_n as f64).ln() - 2.0 * log_likelihood)?;
    let mean = evaluation.mean;
    let innovations: Vec<_> = evaluation.innovations.iter().map(|v| v * scale).collect();
    let mut errors = evaluation.innovations;
    let mut forecasts = Vec::with_capacity(s.horizon);
    for _ in 0..s.horizon {
        let t = y.len();
        let mut next = mean;
        for lag in 1..m.a.len() {
            next -= m.a[lag] * (y[t - lag] - mean);
        }
        for lag in 1..m.c.len() {
            next += m.c[lag] * errors[t - lag];
        }
        finite(next)?;
        y.push(next);
        errors.push(0.0);
        let mut level = next * scale + normalization_origin;
        layers.last_mut().unwrap().push(level);
        for stage in (0..lags.len()).rev() {
            let lower = &mut layers[stage];
            level += lower[lower.len() - lags[stage]];
            finite(level)?;
            lower.push(level);
        }
        forecasts.push(finite(level)?);
    }
    let difference_lags: Vec<_> = lags.iter().map(|v| *v as f64).collect();
    let mean_raw = finite(mean * scale + normalization_origin)?;
    let equation_constant = finite(mean_raw * m.a.iter().sum::<f64>())?;
    Ok(format!("{{\"version\":1,\"operation\":\"arima\",\"estimation_method\":\"conditional_sum_of_squares\",\"presample\":\"zero_deviations_and_innovations\",\"order\":[{},{},{}],\"seasonal_order\":[{},{},{},{}],\"effective_n\":{effective_n},\"conditioning_burn\":{burn},\"differencing_lags\":{},\"differenced_mean\":{mean_raw},\"ar_equation_constant\":{equation_constant},\"mean_interpretation\":\"constant_mean_of_differenced_series\",\"ar_coefficients\":{},\"ma_coefficients\":{},\"seasonal_ar_coefficients\":{},\"seasonal_ma_coefficients\":{},\"ar_polynomial\":{},\"ma_polynomial\":{},\"innovations\":{},\"sigma_squared\":{sigma_squared},\"conditional_log_likelihood\":{log_likelihood},\"aic\":{aic},\"bic\":{bic},\"forecasts\":{},\"converged\":{},\"iterations\":{},\"stationarity_enforced\":true,\"invertibility_enforced\":true}}",s.p,s.d,s.q,s.sp,s.sd,s.sq,s.season,json_array(&difference_lags)?,json_array(&m.ar)?,json_array(&m.ma)?,json_array(&m.seasonal_ar)?,json_array(&m.seasonal_ma)?,json_array(&m.a)?,json_array(&m.c)?,json_array(&innovations)?,json_array(&forecasts)?,optimized.converged,optimized.iterations))
}
