//! Time-series transformations and conditional forecast uncertainty.
#![allow(clippy::needless_range_loop)]
use crate::state_space::{confidence, dot, finite, integer, ls, scaled_variance};
fn scalar(data: &[Vec<f64>], min: usize) -> Result<Vec<f64>, String> {
    if data.len() < min
        || data.len() > 10000
        || data.iter().any(|r| r.len() != 1 || !r[0].is_finite())
    {
        return Err("Expected finite scalar rows within resource limits".into());
    }
    Ok(data.iter().map(|r| r[0]).collect())
}
fn standardized(y: &[f64]) -> Result<(f64, f64, Vec<f64>), String> {
    let origin = y[0];
    let dif = y.iter().map(|v| v - origin).collect::<Vec<_>>();
    finite(&dif)?;
    let scale = dif.iter().map(|v| v.abs()).fold(0., f64::max);
    if scale == 0. {
        return Err("Constant series has no variance".into());
    }
    let z = dif.iter().map(|v| v / scale).collect::<Vec<_>>();
    let mean = z.iter().sum::<f64>() / z.len() as f64;
    Ok((origin, scale, z.iter().map(|v| v - mean).collect()))
}
fn intervals(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let means = scalar(data, 1)?;
    if means.len() > 1000 || a.len() < 7 {
        return Err("Forecast interval layout/limit invalid".into());
    }
    finite(a)?;
    let z = confidence(a[0])?;
    let s = a[1];
    if s < 0. {
        return Err("Innovation variance must be nonnegative".into());
    }
    let pa = integer(a[2], 0, 200, "AR degree")?;
    let pc = integer(a[3], 0, 200, "MA degree")?;
    let nd = integer(a[4], 0, 4, "difference factors")?;
    if a.len() != 7 + pa + pc + nd {
        return Err("Forecast polynomial layout length mismatch".into());
    }
    let ar = &a[5..6 + pa];
    let ma = &a[6 + pa..7 + pa + pc];
    if ar[0] != 1. || ma[0] != 1. {
        return Err("Polynomial leading coefficients must be one".into());
    }
    let mut den = ar.to_vec();
    for lag in &a[7 + pa + pc..] {
        let lag = integer(*lag, 1, 100, "difference lag")?;
        if den.len() + lag > 500 {
            return Err("Differencing polynomial degree exceeds limit".into());
        }
        let mut next = vec![0.; den.len() + lag];
        for i in 0..den.len() {
            next[i] += den[i];
            next[i + lag] -= den[i];
        }
        den = next;
    }
    let n = means.len();
    let mut psi = vec![0.; n];
    let mut vars = vec![0.; n];
    let mut cumulative = 0.;
    for i in 0..n {
        let mut v = if i < ma.len() { ma[i] } else { 0. };
        for j in 1..=i.min(den.len() - 1) {
            v -= den[j] * psi[i - j];
        }
        psi[i] = v;
        cumulative += v * v;
        vars[i] = s * cumulative;
    }
    finite(&vars)?;
    let lo = means
        .iter()
        .zip(&vars)
        .map(|(m, v)| m - z * v.sqrt())
        .collect::<Vec<_>>();
    let hi = means
        .iter()
        .zip(&vars)
        .map(|(m, v)| m + z * v.sqrt())
        .collect::<Vec<_>>();
    finite(&lo)?;
    finite(&hi)?;
    Ok(format!("{{\"version\":1,\"operation\":\"arima_forecast_intervals\",\"forecast_means\":{means:?},\"forecast_variances\":{vars:?},\"impulse_response\":{psi:?},\"lower\":{lo:?},\"upper\":{hi:?},\"parameter_uncertainty_included\":false}}"))
}
fn ccf(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if data.len() < 3
        || data.len() > 10000
        || data
            .iter()
            .any(|r| r.len() != 2 || r.iter().any(|x| !x.is_finite()))
        || a.is_empty()
        || a.len() > 2
    {
        return Err("CCF requires two complete series and [maxlag,normalization=0]".into());
    }
    finite(a)?;
    let n = data.len();
    let maxlag = integer(a[0], 0, (n - 2).min(500), "maxlag")?;
    let mode = integer(*a.get(1).unwrap_or(&0.), 0, 1, "normalization")?;
    let x = data.iter().map(|r| r[0]).collect::<Vec<_>>();
    let y = data.iter().map(|r| r[1]).collect::<Vec<_>>();
    let (_, _, x) = standardized(&x)?;
    let (_, _, y) = standardized(&y)?;
    let mut lags = vec![];
    let mut corr = vec![];
    let mut overlaps = vec![];
    let mut se = vec![];
    for lag in -(maxlag as isize)..=maxlag as isize {
        let count = n - lag.unsigned_abs();
        let (xs, ys) = if lag >= 0 {
            (&x[..count], &y[lag as usize..])
        } else {
            (&x[(-lag) as usize..], &y[..count])
        };
        let v = if mode == 0 {
            dot(xs, ys) / (dot(&x, &x).sqrt() * dot(&y, &y).sqrt())
        } else {
            let mx = xs.iter().sum::<f64>() / count as f64;
            let my = ys.iter().sum::<f64>() / count as f64;
            let xx = xs.iter().map(|v| v - mx).collect::<Vec<_>>();
            let yy = ys.iter().map(|v| v - my).collect::<Vec<_>>();
            let den = dot(&xx, &xx).sqrt() * dot(&yy, &yy).sqrt();
            if den == 0. {
                return Err("Constant overlap has undefined Pearson correlation".into());
            }
            dot(&xx, &yy) / den
        };
        lags.push(lag);
        corr.push(v.clamp(-1., 1.));
        overlaps.push(count);
        se.push(1. / (count as f64).sqrt());
    }
    let lower = se.iter().map(|v| -2. * v).collect::<Vec<_>>();
    let upper = se.iter().map(|v| 2. * v).collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"operation\":\"cross_correlation\",\"positive_lag\":\"input_leads_output\",\"normalization\":{mode},\"lags\":{lags:?},\"correlations\":{corr:?},\"overlaps\":{overlaps:?},\"white_noise_standard_errors\":{se:?},\"lower\":{lower:?},\"upper\":{upper:?}}}"))
}
fn trend(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let y = scalar(data, 3)?;
    if a.len() > 1 {
        return Err("Trend takes [n_train=all]".into());
    }
    let n = integer(
        *a.first().unwrap_or(&(y.len() as f64)),
        3,
        y.len(),
        "n_train",
    )?;
    let origin = y[0];
    let scale = y[..n]
        .iter()
        .map(|v| (v - origin).abs())
        .fold(0., f64::max)
        .max(f64::MIN_POSITIVE);
    let yy = y[..n]
        .iter()
        .map(|v| (v - origin) / scale)
        .collect::<Vec<_>>();
    finite(&yy)?;
    let mid = (n - 1) as f64 / 2.;
    let ts = (n - 1) as f64;
    let x = (0..n)
        .map(|i| vec![1., (i as f64 - mid) / ts])
        .collect::<Vec<_>>();
    let (b, _) = ls(&x, &yy)?;
    let slope = scale * b[1] / ts;
    let intercept = origin + scale * b[0] - slope * mid;
    let fit = (0..y.len())
        .map(|i| origin + scale * (b[0] + b[1] * (i as f64 - mid) / ts))
        .collect::<Vec<_>>();
    let res = y.iter().zip(&fit).map(|(y, f)| y - f).collect::<Vec<_>>();
    finite(&fit)?;
    finite(&res)?;
    finite(&[intercept, slope])?;
    Ok(format!("{{\"version\":1,\"operation\":\"remove_linear_trend\",\"n_train\":{n},\"intercept\":{intercept},\"slope\":{slope},\"trend\":{fit:?},\"residuals\":{res:?},\"prediction_model\":[1,{origin},{scale},{},{},{mid},{ts}]}}",b[0],b[1]))
}
fn cycle(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let y = scalar(data, 4)?;
    finite(a)?;
    if a.is_empty() || a.len() > 3 || a[0] <= 2. || a[0] > 1e6 {
        return Err("Cycle takes [period>2,include_constant=1,n_train=all]".into());
    }
    let period = a[0];
    let constant = integer(*a.get(1).unwrap_or(&1.), 0, 1, "include_constant")?;
    let n = integer(
        *a.get(2).unwrap_or(&(y.len() as f64)),
        4,
        y.len(),
        "n_train",
    )?;
    let origin = if constant == 1 { y[0] } else { 0. };
    let scale = y[..n]
        .iter()
        .map(|v| (v - origin).abs())
        .fold(0., f64::max)
        .max(f64::MIN_POSITIVE);
    let yy = y[..n]
        .iter()
        .map(|v| (v - origin) / scale)
        .collect::<Vec<_>>();
    finite(&yy)?;
    let design = (0..y.len())
        .map(|i| {
            let w = 2. * std::f64::consts::PI * (i as f64 / period);
            let mut r = vec![w.cos(), w.sin()];
            if constant == 1 {
                r.push(1.);
            }
            r
        })
        .collect::<Vec<_>>();
    let (b, _) = ls(&design[..n].to_vec(), &yy)?;
    let cos = scale * b[0];
    let sin = scale * b[1];
    let intercept = if constant == 1 {
        origin + scale * b[2]
    } else {
        0.
    };
    let amp = cos.hypot(sin);
    let phase = if amp <= 1e-12 * scale {
        "null".into()
    } else {
        (-sin).atan2(cos).to_string()
    };
    let component = design
        .iter()
        .map(|r| intercept + cos * r[0] + sin * r[1])
        .collect::<Vec<_>>();
    let res = y
        .iter()
        .zip(&component)
        .map(|(y, c)| y - c)
        .collect::<Vec<_>>();
    finite(&res)?;
    finite(&component)?;
    finite(&[cos, sin, intercept, amp])?;
    Ok(format!("{{\"version\":1,\"operation\":\"remove_cycle\",\"period\":{period},\"intercept\":{intercept},\"cosine_coefficient\":{cos},\"sine_coefficient\":{sin},\"amplitude\":{amp},\"phase\":{phase},\"component\":{component:?},\"residuals\":{res:?}}}"))
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    match operation {
        "arima_prewhitening" => prewhite(data, args),
        "transfer_function_model" => transfer(data, args),
        "arima_forecast_intervals" => intervals(data, args),
        "cross_correlation" => ccf(data, args),
        "remove_linear_trend" => trend(data, args),
        "remove_cycle" => cycle(data, args),
        _ => Err("Unknown dynamic-series operation".into()),
    }
}
use crate::state_space::{bounded_opt, chol, eye, mul, mv, trans, zeros};
type Mat = Vec<Vec<f64>>;
#[derive(Clone)]
struct Arma {
    f: Mat,
    h: Vec<f64>,
    p0: Mat,
    ar: Vec<f64>,
    ma: Vec<f64>,
}
fn pacf(u: &[f64]) -> Vec<f64> {
    let k = u.iter().map(|v| 0.98 * (2. * v - 1.)).collect::<Vec<_>>();
    match k.len() {
        0 => vec![],
        1 => k,
        _ => vec![k[0] * (1. - k[1]), k[1]],
    }
}
fn arma(p: usize, q: usize, u: &[f64]) -> Result<Arma, String> {
    let ar = pacf(&u[..p]);
    let ma = pacf(&u[p..p + q]).iter().map(|v| -v).collect::<Vec<_>>();
    let k = p.max(q + 1).max(1);
    let mut f = zeros(k, k);
    f[0][..p].copy_from_slice(&ar);
    for i in 1..k {
        f[i][i - 1] = 1.;
    }
    let mut h = vec![0.; k];
    h[0] = 1.;
    h[1..q + 1].copy_from_slice(&ma);
    let mut cov = zeros(k, k);
    cov[0][0] = 1.;
    let mut power = f.clone();
    let mut done = false;
    for _ in 0..40 {
        let inc = mul(&mul(&power, &cov), &trans(&power));
        let size = cov.iter().flatten().map(|x| x.abs()).fold(0., f64::max);
        let change = inc.iter().flatten().map(|x| x.abs()).fold(0., f64::max);
        for i in 0..k {
            for j in 0..k {
                cov[i][j] += inc[i][j];
            }
        }
        power = mul(&power, &power);
        if change < 1e-14 * size.max(1.) {
            done = true;
            break;
        }
    }
    if !done || cov.iter().flatten().any(|v| !v.is_finite()) {
        return Err("Stationary ARMA covariance failed".into());
    }
    Ok(Arma {
        f,
        h,
        p0: cov,
        ar,
        ma,
    })
}
/// Applies the same exact stationary innovations whitening operator to all columns.
fn whiten(a: &Arma, x: &Mat) -> Result<(Mat, f64, Vec<f64>, Vec<Mat>), String> {
    let n = x.len();
    let cols = x[0].len();
    let k = a.h.len();
    let mut state = zeros(k, cols);
    let mut cov = a.p0.clone();
    let mut out = zeros(n, cols);
    let mut logdet = 0.;
    let mut last_cov = cov.clone();
    let mut last_state = state.clone();
    for t in 0..n {
        let ph = mv(&cov, &a.h);
        let s = dot(&a.h, &ph);
        if !s.is_finite() || s <= 1e-13 {
            return Err("Singular ARMA innovation variance".into());
        }
        logdet += s.ln();
        let gain = ph.iter().map(|v| v / s).collect::<Vec<_>>();
        for j in 0..cols {
            let mean = (0..k).map(|i| a.h[i] * state[i][j]).sum::<f64>();
            let e = x[t][j] - mean;
            out[t][j] = e / s.sqrt();
            for i in 0..k {
                state[i][j] += gain[i] * e;
            }
        }
        let mut ikh = eye(k);
        for i in 0..k {
            for j in 0..k {
                ikh[i][j] -= gain[i] * a.h[j];
            }
        }
        cov = mul(&mul(&ikh, &cov), &trans(&ikh));
        last_cov = cov.clone();
        last_state = state.clone();
        state = mul(&a.f, &state);
        cov = mul(&mul(&a.f, &cov), &trans(&a.f));
        cov[0][0] += 1.;
    }
    finite(&out.iter().flatten().copied().collect::<Vec<_>>())?;
    Ok((
        out,
        logdet,
        last_state.iter().map(|r| r[0]).collect(),
        vec![last_cov],
    ))
}
struct Profile {
    beta: Vec<f64>,
    sigma: f64,
    logdet: f64,
    residual: Vec<f64>,
    a: Arma,
}
fn gls(y: &[f64], x: &Mat, p: usize, q: usize, u: &[f64]) -> Result<Profile, String> {
    let a = arma(p, q, u)?;
    let joined = y
        .iter()
        .zip(x)
        .map(|(y, r)| {
            let mut v = vec![*y];
            v.extend(r);
            v
        })
        .collect::<Mat>();
    let (w, logdet, _, _) = whiten(&a, &joined)?;
    let yy = w.iter().map(|r| r[0]).collect::<Vec<_>>();
    let xx = w.iter().map(|r| r[1..].to_vec()).collect::<Mat>();
    let (beta, residual) = ls(&xx, &yy)?;
    let sigma = dot(&residual, &residual) / y.len() as f64;
    if sigma <= 1e-24 || !sigma.is_finite() {
        return Err("Degenerate Gaussian innovation variance".into());
    }
    Ok(Profile {
        beta,
        sigma,
        logdet,
        residual,
        a,
    })
}
fn gradient_curvature<F: Fn(&[f64]) -> f64>(f: &F, u: &[f64]) -> bool {
    let free = (0..u.len())
        .filter(|i| u[*i] > 1e-4 && u[*i] < 1. - 1e-4)
        .collect::<Vec<_>>();
    if free.is_empty() {
        return true;
    }
    let h = 1e-4;
    let mut info = zeros(free.len(), free.len());
    for (ii, i) in free.iter().enumerate() {
        for (jj, j) in free.iter().enumerate() {
            let mut pp = u.to_vec();
            let mut pm = u.to_vec();
            let mut mp = u.to_vec();
            let mut mm = u.to_vec();
            pp[*i] += h;
            pp[*j] += h;
            pm[*i] += h;
            pm[*j] -= h;
            mp[*i] -= h;
            mp[*j] += h;
            mm[*i] -= h;
            mm[*j] -= h;
            info[ii][jj] = (f(&pp) - f(&pm) - f(&mp) + f(&mm)) / (4. * h * h);
        }
    }
    chol(&info).is_ok()
}
fn bivariate(data: &[Vec<f64>]) -> Result<(), String> {
    if data.is_empty()
        || data.len() > 10000
        || data
            .iter()
            .any(|r| r.len() != 2 || r.iter().any(|v| !v.is_finite()))
    {
        Err("Expected finite [input,output] rows within limits".into())
    } else {
        Ok(())
    }
}
fn settings(a: &[f64], maxiter_at: usize) -> Result<(usize, f64), String> {
    let iterations = integer(*a.get(maxiter_at).unwrap_or(&1000.), 1, 5000, "maxiter")?;
    let tolerance = *a.get(maxiter_at + 1).unwrap_or(&1e-8);
    if !(1e-10..=1e-3).contains(&tolerance) {
        return Err("Tolerance must be [1e-10,1e-3]".into());
    }
    Ok((iterations, tolerance))
}
fn prewhite(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    bivariate(data)?;
    finite(args)?;
    if args.len() < 3 || args.len() > 6 {
        return Err(
            "Prewhitening takes [n_train,p,q,maxlag=10,maxiter=1000,tolerance=1e-8]".into(),
        );
    }
    let n = integer(args[0], 12, data.len().min(2000), "n_train")?;
    let p = integer(args[1], 0, 2, "AR order")?;
    let q = integer(args[2], 0, 2, "MA order")?;
    let maxlag = integer(*args.get(3).unwrap_or(&10.), 0, (n - 2).min(500), "maxlag")?;
    let (iterations, tol) = settings(args, 4)?;
    if n * iterations * (1 + p + q) > 40_000_000 {
        return Err("Prewhitening work limit exceeded".into());
    }
    let raw = data[..n].iter().map(|r| r[0]).collect::<Vec<_>>();
    let (origin, scale, _) = standardized(&raw)?;
    let y = raw.iter().map(|v| (v - origin) / scale).collect::<Vec<_>>();
    let design = vec![vec![1.]; n];
    let objective = |u: &[f64]| {
        gls(&y, &design, p, q, u)
            .map_or(f64::INFINITY, |fit| fit.sigma.ln() + fit.logdet / n as f64)
    };
    let d = p + q;
    let starts = if d == 0 {
        vec![vec![]]
    } else {
        vec![vec![0.5; d], vec![0.25; d], vec![0.75; d]]
    };
    let opt = bounded_opt(objective, &starts, iterations, tol);
    let fit = gls(&y, &design, p, q, &opt.x)?;
    let sigma = scaled_variance(fit.sigma, scale)?;
    let mean = origin + scale * fit.beta[0];
    let ll =
        -0.5 * n as f64 * ((2. * std::f64::consts::PI).ln() + 1. + sigma.ln()) - 0.5 * fit.logdet;
    let input = fit.residual.iter().map(|v| scale * v).collect::<Vec<_>>();
    let rawy = data[..n].iter().map(|r| r[1]).collect::<Vec<_>>();
    let yori = rawy[0];
    let yscale = rawy.iter().map(|v| (v - yori).abs()).fold(0., f64::max);
    if !yscale.is_finite() {
        return Err("Output scale cannot be represented".into());
    }
    let output = if yscale == 0. {
        vec![0.; n]
    } else {
        let z = rawy.iter().map(|v| (v - yori) / yscale).collect::<Vec<_>>();
        let mean = z.iter().sum::<f64>() / n as f64;
        let z = z.iter().map(|v| vec![v - mean]).collect::<Mat>();
        let (w, _, _, _) = whiten(&fit.a, &z)?;
        w.iter().map(|r| r[0] * yscale).collect()
    };
    let pairs = input
        .iter()
        .zip(&output)
        .map(|(x, y)| vec![*x, *y])
        .collect::<Mat>();
    let ccf = ccf(&pairs, &[maxlag as f64]).unwrap_or("null".into());
    let guard = opt.x.iter().any(|v| *v < 1e-6 || *v > 1. - 1e-6);
    let identified = gradient_curvature(&objective, &opt.x);
    let converged = opt.converged && !guard && identified;
    finite(&[mean, sigma, ll, opt.projected_gradient])?;
    finite(&input)?;
    finite(&output)?;
    Ok(format!("{{\"version\":1,\"operation\":\"arima_prewhitening\",\"likelihood\":\"exact_stationary_gaussian\",\"n_train\":{n},\"ar\":{:?},\"ma\":{:?},\"mean\":{mean},\"sigma_squared\":{sigma},\"log_likelihood\":{ll},\"input_innovations\":{input:?},\"output_innovations\":{output:?},\"ccf\":{ccf},\"converged\":{converged},\"optimizer_converged\":{},\"information_identified\":{identified},\"numerical_boundary_guard\":{guard},\"projected_gradient\":{},\"iterations\":{},\"output_used_for_input_fit\":false}}",fit.a.ar,fit.a.ma,opt.converged,opt.projected_gradient,opt.iterations))
}
fn transfer_design(x: &[f64], lag: usize, num: usize, delta: &[f64]) -> Mat {
    let mut design = vec![vec![1.; num + 2]; x.len()];
    for j in 0..=num {
        for t in 0..x.len() {
            let mut v = if t >= lag + j { x[t - lag - j] } else { 0. };
            for k in 1..=delta.len() {
                if t >= k {
                    v += delta[k - 1] * design[t - k][j + 1];
                }
            }
            design[t][j + 1] = v;
        }
    }
    design
}
fn noise_forecast(
    a: &Arma,
    residual: &[f64],
    sigma: f64,
    horizon: usize,
) -> Result<(Vec<f64>, Vec<f64>), String> {
    let rows = residual.iter().map(|v| vec![*v]).collect::<Mat>();
    let (_, _, mut state, covs) = whiten(a, &rows)?;
    let mut cov = covs[0].clone();
    let mut means = vec![];
    let mut vars = vec![];
    for _ in 0..horizon {
        state = mv(&a.f, &state);
        cov = mul(&mul(&a.f, &cov), &trans(&a.f));
        cov[0][0] += 1.;
        let variance = sigma * dot(&a.h, &mv(&cov, &a.h));
        if variance < 0. {
            return Err("Negative ARMA forecast variance".into());
        }
        means.push(dot(&a.h, &state));
        vars.push(variance);
    }
    finite(&means)?;
    finite(&vars)?;
    Ok((means, vars))
}
fn transfer(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    bivariate(data)?;
    finite(a)?;
    if a.len() < 7 || a.len() > 10 {
        return Err("Transfer takes [n_train,lag,numerator_degree,denominator_degree,p,q,horizon,confidence=.95,maxiter=1000,tolerance=1e-8]".into());
    }
    let n = integer(a[0], 20, data.len().min(2000), "n_train")?;
    let lag = integer(a[1], 0, 24, "input delay")?;
    let num = integer(a[2], 0, 2, "numerator degree")?;
    let den = integer(a[3], 0, 2, "denominator degree")?;
    let p = integer(a[4], 0, 2, "AR order")?;
    let q = integer(a[5], 0, 2, "MA order")?;
    let horizon = integer(a[6], 1, 500, "horizon")?;
    if data.len() != n + horizon || n < lag + num + 8 {
        return Err(
            "Transfer requires exactly n_train+horizon rows and enough delayed observations".into(),
        );
    }
    let conf = *a.get(7).unwrap_or(&0.95);
    let z = confidence(conf)?;
    let (iterations, tol) = settings(a, 8)?;
    let d = den + p + q;
    if n * iterations * (d + 1) > 60_000_000 {
        return Err("Transfer optimization work limit exceeded".into());
    }
    let rawx = data[..n].iter().map(|r| r[0]).collect::<Vec<_>>();
    let (xorigin, xscale, xc) = standardized(&rawx)?;
    let xmean = rawx.iter().map(|v| (v - xorigin) / xscale).sum::<f64>() / n as f64;
    let x = data
        .iter()
        .map(|r| (r[0] - xorigin) / xscale - xmean)
        .collect::<Vec<_>>();
    finite(&x)?;
    let _ = xc;
    let rawy = data[..n].iter().map(|r| r[1]).collect::<Vec<_>>();
    let (yorigin, yscale, _) = standardized(&rawy)?;
    let y = rawy
        .iter()
        .map(|v| (v - yorigin) / yscale)
        .collect::<Vec<_>>();
    let objective = |u: &[f64]| {
        let delta = pacf(&u[..den]);
        let design = transfer_design(&x[..n], lag, num, &delta);
        gls(&y, &design, p, q, &u[den..])
            .map_or(f64::INFINITY, |f| f.sigma.ln() + f.logdet / n as f64)
    };
    let starts = if d == 0 {
        vec![vec![]]
    } else {
        vec![
            vec![0.5; d],
            vec![0.25; d],
            vec![0.75; d],
            (0..d).map(|i| if i % 2 == 0 { 0.7 } else { 0.3 }).collect(),
        ]
    };
    let opt = bounded_opt(objective, &starts, iterations, tol);
    if !opt.value.is_finite() {
        return Err("No finite full-rank transfer fit".into());
    }
    let delta = pacf(&opt.x[..den]);
    let design = transfer_design(&x, lag, num, &delta);
    let fit = gls(&y, &design[..n].to_vec(), p, q, &opt.x[den..])?;
    let sigma = scaled_variance(fit.sigma, yscale)?;
    let intercept = yorigin + yscale * fit.beta[0];
    let numerator = fit.beta[1..]
        .iter()
        .map(|v| yscale * v / xscale)
        .collect::<Vec<_>>();
    let ll =
        -0.5 * n as f64 * ((2. * std::f64::consts::PI).ln() + 1. + sigma.ln()) - 0.5 * fit.logdet;
    let residual = y
        .iter()
        .zip(&design)
        .map(|(y, r)| y - dot(r, &fit.beta))
        .collect::<Vec<_>>();
    let (noise, var) = noise_forecast(&fit.a, &residual, fit.sigma, horizon)?;
    let forecasts = (0..horizon)
        .map(|i| yorigin + yscale * (dot(&design[n + i], &fit.beta) + noise[i]))
        .collect::<Vec<_>>();
    let vars = var
        .iter()
        .map(|v| scaled_variance(*v, yscale))
        .collect::<Result<Vec<_>, _>>()?;
    let lower = forecasts
        .iter()
        .zip(&vars)
        .map(|(m, v)| m - z * v.sqrt())
        .collect::<Vec<_>>();
    let upper = forecasts
        .iter()
        .zip(&vars)
        .map(|(m, v)| m + z * v.sqrt())
        .collect::<Vec<_>>();
    let fitted = design[..n]
        .iter()
        .map(|r| yorigin + yscale * dot(r, &fit.beta))
        .collect::<Vec<_>>();
    let innovations = fit.residual.iter().map(|v| yscale * v).collect::<Vec<_>>();
    let guard = opt.x.iter().any(|v| *v < 1e-6 || *v > 1. - 1e-6);
    let identified = gradient_curvature(&objective, &opt.x);
    let converged = opt.converged && !guard && identified;
    let nparams = d + num + 3;
    let aic = if converged {
        (2. * nparams as f64 - 2. * ll).to_string()
    } else {
        "null".into()
    };
    let bic = if converged {
        ((n as f64).ln() * nparams as f64 - 2. * ll).to_string()
    } else {
        "null".into()
    };
    for v in [
        &forecasts,
        &vars,
        &lower,
        &upper,
        &fitted,
        &innovations,
        &numerator,
    ] {
        finite(v)?;
    }
    finite(&[intercept, sigma, ll, opt.projected_gradient])?;
    Ok(format!("{{\"version\":1,\"operation\":\"transfer_function_model\",\"likelihood\":\"exact_stationary_gaussian_noise\",\"input_initialization\":\"zero_centered_presample\",\"n_train\":{n},\"input_delay\":{lag},\"input_center\":{},\"intercept\":{intercept},\"numerator\":{numerator:?},\"denominator\":{delta:?},\"ar\":{:?},\"ma\":{:?},\"sigma_squared\":{sigma},\"log_likelihood\":{ll},\"aic\":{aic},\"bic\":{bic},\"fitted_transfer_means\":{fitted:?},\"innovations\":{innovations:?},\"forecasts\":{forecasts:?},\"forecast_variances\":{vars:?},\"lower\":{lower:?},\"upper\":{upper:?},\"converged\":{converged},\"optimizer_converged\":{},\"information_identified\":{identified},\"numerical_boundary_guard\":{guard},\"projected_gradient\":{},\"iterations\":{},\"future_output_used\":false,\"parameter_uncertainty_included\":false,\"future_input_uncertainty_included\":false}}",xorigin+xscale*xmean,fit.a.ar,fit.a.ma,opt.converged,opt.projected_gradient,opt.iterations))
}
