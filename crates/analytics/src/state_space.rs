//! Known-parameter Gaussian state-space inference and additive innovations ETS.
//! These are distinct disturbance conventions; portable models carry an explicit tag.
#![allow(clippy::needless_range_loop)]
type Mat = Vec<Vec<f64>>;
type Forecast = (Vec<f64>, Vec<f64>, Vec<f64>);
pub(crate) fn zeros(n: usize, m: usize) -> Mat {
    vec![vec![0.; m]; n]
}
pub(crate) fn eye(n: usize) -> Mat {
    let mut a = zeros(n, n);
    for i in 0..n {
        a[i][i] = 1.;
    }
    a
}
pub(crate) fn trans(a: &Mat) -> Mat {
    let mut b = zeros(a[0].len(), a.len());
    for i in 0..a.len() {
        for j in 0..a[0].len() {
            b[j][i] = a[i][j];
        }
    }
    b
}
pub(crate) fn mul(a: &Mat, b: &Mat) -> Mat {
    let mut c = zeros(a.len(), b[0].len());
    for i in 0..a.len() {
        for j in 0..b[0].len() {
            for k in 0..b.len() {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    c
}
pub(crate) fn mv(a: &Mat, x: &[f64]) -> Vec<f64> {
    a.iter().map(|r| dot(r, x)).collect()
}
pub(crate) fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn add(a: &Mat, b: &Mat, sign: f64) -> Mat {
    a.iter()
        .zip(b)
        .map(|(a, b)| a.iter().zip(b).map(|(x, y)| x + sign * y).collect())
        .collect()
}
pub(crate) fn chol(a: &Mat) -> Result<Mat, String> {
    let n = a.len();
    let scale = a.iter().flatten().map(|x| x.abs()).fold(0., f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err("Positive definite covariance required".into());
    }
    let mut l = zeros(n, n);
    for i in 0..n {
        for j in 0..=i {
            let mut v = a[i][j] / scale;
            for k in 0..j {
                v -= l[i][k] * l[j][k];
            }
            if i == j {
                if v <= 1e-14 || !v.is_finite() {
                    return Err("Numerically singular covariance".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let s = scale.sqrt();
    for r in &mut l {
        for v in r {
            *v *= s;
        }
    }
    Ok(l)
}
pub(crate) fn solve_chol(l: &Mat, b: &[f64]) -> Vec<f64> {
    let n = b.len();
    let mut x = b.to_vec();
    for i in 0..n {
        for j in 0..i {
            x[i] -= l[i][j] * x[j];
        }
        x[i] /= l[i][i];
    }
    for i in (0..n).rev() {
        for j in i + 1..n {
            x[i] -= l[j][i] * x[j];
        }
        x[i] /= l[i][i];
    }
    x
}
fn validate_cov(a: &Mat, definite: bool, external: bool) -> Result<(), String> {
    let n = a.len();
    let s = a.iter().flatten().map(|v| v.abs()).fold(0., f64::max);
    if !s.is_finite() {
        return Err("Nonfinite covariance".into());
    }
    if s == 0. {
        return if definite {
            Err("Proper positive definite initial covariance required".into())
        } else {
            Ok(())
        };
    }
    let mut l = zeros(n, n);
    for i in 0..n {
        for j in 0..n {
            if (external && a[i][j] != a[j][i])
                || (!external && (a[i][j] - a[j][i]).abs() > 1e-12 * s)
            {
                return Err("Covariance must be symmetric".into());
            }
        }
        for j in 0..=i {
            let mut v = a[i][j] / s;
            for k in 0..j {
                v -= l[i][k] * l[j][k];
            }
            if !v.is_finite() {
                return Err("Covariance factorization exceeded numeric range".into());
            }
            if i == j {
                if v < if external { 0. } else { -1e-12 } || (definite && v <= 1e-14) {
                    return Err(
                        "Covariance must be positive semidefinite (initial positive definite)"
                            .into(),
                    );
                }
                l[i][j] = v.max(0.).sqrt();
            } else if l[j][j] > if external { 0. } else { 1e-14 } {
                l[i][j] = v / l[j][j];
            } else if (external && v != 0.) || (!external && v.abs() > 1e-10) {
                return Err("Indefinite covariance".into());
            }
        }
    }
    Ok(())
}
fn sym(a: &mut Mat) {
    for i in 0..a.len() {
        for j in 0..i {
            let v = 0.5 * a[i][j] + 0.5 * a[j][i];
            a[i][j] = v;
            a[j][i] = v;
        }
    }
}
pub(crate) fn integer(v: f64, lo: usize, hi: usize, name: &str) -> Result<usize, String> {
    if v.is_finite() && v.fract() == 0. && v >= lo as f64 && v <= hi as f64 {
        Ok(v as usize)
    } else {
        Err(format!("{name} must be an integer in {lo}..={hi}"))
    }
}
pub(crate) fn normal_quantile(p: f64) -> Result<f64, String> {
    if !(0.5..1.).contains(&p) {
        return Err("Normal quantile probability outside (0.5,1)".into());
    }
    let a = [
        -39.69683028665376,
        220.9460984245205,
        -275.9285104469687,
        138.357751867269,
        -30.66479806614716,
        2.506628277459239,
    ];
    let b = [
        -54.47609879822406,
        161.5858368580409,
        -155.6989798598866,
        66.80131188771972,
        -13.28068155288572,
    ];
    let c = [
        -0.007784894002430293,
        -0.3223964580411365,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    let d = [
        0.007784695709041462,
        0.3224671290700398,
        2.445134137142996,
        3.754408661907416,
    ];
    if p > 0.97575 {
        let q = (-2. * (1. - p).ln()).sqrt();
        return Ok(
            -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
                / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.),
        );
    }
    let q = p - 0.5;
    let r = q * q;
    Ok(
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.),
    )
}
pub(crate) fn confidence(c: f64) -> Result<f64, String> {
    if !c.is_finite() || !(0.5..=0.9999).contains(&c) {
        return Err("Confidence must be in [0.5,0.9999]".into());
    }
    normal_quantile(0.5 + 0.5 * c)
}
pub(crate) fn finite(v: &[f64]) -> Result<(), String> {
    if v.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err("Nonfinite input or numeric range exceeded".into())
    }
}
pub(crate) fn scaled_variance(v: f64, scale: f64) -> Result<f64, String> {
    if v < 0. || !v.is_finite() || scale <= 0. || !scale.is_finite() {
        return Err("Invalid variance scaling".into());
    }
    if v == 0. {
        return Ok(0.);
    }
    let sd = v.sqrt() * scale;
    let out = sd * sd;
    if !out.is_finite() || out == 0. {
        return Err("Scaled positive variance is not representable".into());
    }
    Ok(out)
}
pub(crate) fn ls(x: &Mat, y: &[f64]) -> Result<(Vec<f64>, Vec<f64>), String> {
    let n = x.len();
    let p = x[0].len();
    if n < p {
        return Err("Insufficient rows for full-rank model".into());
    }
    let mut q = zeros(p, n);
    let mut r = zeros(p, p);
    for j in 0..p {
        let mut v = x.iter().map(|r| r[j]).collect::<Vec<_>>();
        let norm0 = v.iter().fold(0_f64, |a, b| a.hypot(*b));
        for _ in 0..2 {
            for i in 0..j {
                let vq = dot(&q[i], &v);
                r[i][j] += vq;
                for t in 0..n {
                    v[t] -= vq * q[i][t];
                }
            }
        }
        let norm = v.iter().fold(0_f64, |a, b| a.hypot(*b));
        if !norm.is_finite() || norm <= 1e-11 * norm0.max(1.) {
            return Err("Rank deficient or numerically unidentified design".into());
        }
        r[j][j] = norm;
        q[j] = v.iter().map(|v| v / norm).collect();
    }
    let mut b = q.iter().map(|v| dot(v, y)).collect::<Vec<_>>();
    for i in (0..p).rev() {
        for j in i + 1..p {
            b[i] -= r[i][j] * b[j];
        }
        b[i] /= r[i][i];
    }
    let residual = x
        .iter()
        .zip(y)
        .map(|(r, y)| y - dot(r, &b))
        .collect::<Vec<_>>();
    finite(&b)?;
    finite(&residual)?;
    Ok((b, residual))
}
fn flatten(a: &Mat) -> Vec<f64> {
    a.iter().flatten().copied().collect()
}
fn read_mat(a: &[f64], at: &mut usize, k: usize) -> Mat {
    let m = a[*at..*at + k * k].chunks(k).map(|r| r.to_vec()).collect();
    *at += k * k;
    m
}
#[derive(Clone)]
struct Gaussian {
    m: Vec<f64>,
    p: Mat,
    f: Mat,
    q: Mat,
    h: Vec<f64>,
    r: f64,
}
fn g_forecast(g: &Gaussian, horizon: usize) -> Result<Forecast, String> {
    let mut m = g.m.clone();
    let mut p = g.p.clone();
    let mut means = vec![];
    let mut latent = vec![];
    let mut observation = vec![];
    for _ in 0..horizon {
        m = mv(&g.f, &m);
        p = add(&mul(&mul(&g.f, &p), &trans(&g.f)), &g.q, 1.);
        sym(&mut p);
        let v = dot(&g.h, &mv(&p, &g.h));
        if v < -1e-10 * (g.r + v.abs()) {
            return Err("Negative forecast variance".into());
        }
        means.push(dot(&g.h, &m));
        latent.push(v.max(0.));
        observation.push(v.max(0.) + g.r);
    }
    finite(&means)?;
    finite(&latent)?;
    finite(&observation)?;
    Ok((means, latent, observation))
}
fn bounds(means: &[f64], vars: &[f64], z: f64) -> Result<(Vec<f64>, Vec<f64>), String> {
    let mut lo = vec![];
    let mut hi = vec![];
    for (m, v) in means.iter().zip(vars) {
        if *v < 0. {
            return Err("Negative prediction variance".into());
        }
        let d = z * v.sqrt();
        lo.push(m - d);
        hi.push(m + d);
    }
    finite(&lo)?;
    finite(&hi)?;
    Ok((lo, hi))
}
fn g_carrier(g: &Gaussian) -> Vec<f64> {
    let mut v = vec![1., 0., g.m.len() as f64];
    v.extend(&g.m);
    v.extend(flatten(&g.p));
    v.extend(flatten(&g.f));
    v.extend(flatten(&g.q));
    v.extend(&g.h);
    v.push(g.r);
    v
}
fn gaussian(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() < 3 || data.is_empty() || data.len() > 5000 {
        return Err("Kalman requires 1..5000 rows and explicit model parameters".into());
    }
    finite(a)?;
    let k = integer(a[0], 1, 8, "states")?;
    let horizon = integer(a[1], 0, 500, "horizon")?;
    let z = confidence(a[2])?;
    if a.len() != 4 + 2 * k + 3 * k * k {
        return Err("Kalman parameter layout length mismatch".into());
    }
    let mut at = 3;
    let mut m = a[at..at + k].to_vec();
    at += k;
    let mut p = read_mat(a, &mut at, k);
    let f = read_mat(a, &mut at, k);
    let q = read_mat(a, &mut at, k);
    let h = a[at..at + k].to_vec();
    at += k;
    let r = a[at];
    if r <= 0. {
        return Err("Measurement variance must be positive".into());
    }
    validate_cov(&p, true, true)?;
    validate_cov(&q, false, true)?;
    let mut pm = vec![];
    let mut pc = vec![];
    let mut fm = vec![];
    let mut fc = vec![];
    let mut innov = vec![];
    let mut ivar = vec![];
    let mut stdres = vec![];
    let mut ll = 0.;
    let mut observed = 0;
    for row in data {
        if row.len() != 2 || !row.iter().all(|v| v.is_finite()) || (row[1] != 0. && row[1] != 1.) {
            return Err("Rows must be [finite response, observed flag 0/1]".into());
        }
        pm.push(m.clone());
        pc.push(p.clone());
        let ph = mv(&p, &h);
        let s = r + dot(&h, &ph);
        if !s.is_finite() || s <= 0. {
            return Err("Innovation variance overflow or nonpositive".into());
        }
        if row[1] == 1. {
            let e = row[0] - dot(&h, &m);
            let es = e / s.sqrt();
            ll -= 0.5 * ((2. * std::f64::consts::PI).ln() + s.ln() + es * es);
            observed += 1;
            let gain = ph.iter().map(|v| v / s).collect::<Vec<_>>();
            for i in 0..k {
                m[i] += gain[i] * e;
            }
            let mut ikh = eye(k);
            for i in 0..k {
                for j in 0..k {
                    ikh[i][j] -= gain[i] * h[j];
                }
            }
            p = mul(&mul(&ikh, &p), &trans(&ikh));
            for i in 0..k {
                for j in 0..k {
                    p[i][j] += r * gain[i] * gain[j];
                }
            }
            sym(&mut p);
            innov.push(Some(e));
            ivar.push(Some(s));
            stdres.push(Some(es));
        } else {
            innov.push(None);
            ivar.push(None);
            stdres.push(None);
        }
        finite(&m)?;
        finite(&flatten(&p))?;
        fm.push(m.clone());
        fc.push(p.clone());
        m = mv(&f, &m);
        p = add(&mul(&mul(&f, &p), &trans(&f)), &q, 1.);
        sym(&mut p);
    }
    finite(&[ll])?;
    let n = data.len();
    let mut sm = fm.clone();
    let mut sc = fc.clone();
    for t in (0..n - 1).rev() {
        let l = chol(&pc[t + 1])?;
        let pf = mul(&fc[t], &trans(&f));
        let c = pf.iter().map(|r| solve_chol(&l, r)).collect::<Mat>();
        let diff = sm[t + 1]
            .iter()
            .zip(&pm[t + 1])
            .map(|(x, y)| x - y)
            .collect::<Vec<_>>();
        let update = mv(&c, &diff);
        for i in 0..k {
            sm[t][i] += update[i];
        }
        sc[t] = add(
            &fc[t],
            &mul(&mul(&c, &add(&sc[t + 1], &pc[t + 1], -1.)), &trans(&c)),
            1.,
        );
        sym(&mut sc[t]);
        finite(&sm[t])?;
        finite(&flatten(&sc[t]))?;
        validate_cov(&sc[t], false, false)?;
    }
    let g = Gaussian {
        m: fm[n - 1].clone(),
        p: fc[n - 1].clone(),
        f,
        q,
        h,
        r,
    };
    let (means, latent, obs) = g_forecast(&g, horizon)?;
    let (lo, hi) = bounds(&means, &obs, z)?;
    let opts = |v: &[Option<f64>]| {
        format!(
            "[{}]",
            v.iter()
                .map(|x| x.map_or("null".into(), |v| v.to_string()))
                .collect::<Vec<_>>()
                .join(",")
        )
    };
    Ok(format!("{{\"version\":1,\"operation\":\"kalman_smoother\",\"model_kind\":\"independent_gaussian\",\"observed_count\":{observed},\"log_likelihood\":{ll},\"predicted_means\":{pm:?},\"filtered_means\":{fm:?},\"smoothed_means\":{sm:?},\"predicted_covariances\":{:?},\"filtered_covariances\":{:?},\"smoothed_covariances\":{:?},\"innovations\":{},\"innovation_variances\":{},\"standardized_residuals\":{},\"forecast_means\":{means:?},\"forecast_latent_variances\":{latent:?},\"forecast_observation_variances\":{obs:?},\"lower\":{lo:?},\"upper\":{hi:?},\"prediction_model\":{:?},\"parameter_uncertainty_included\":false}}",pc.iter().map(flatten).collect::<Vec<_>>(),fc.iter().map(flatten).collect::<Vec<_>>(),sc.iter().map(flatten).collect::<Vec<_>>(),opts(&innov),opts(&ivar),opts(&stdres),g_carrier(&g)))
}
fn predict(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if (data.len() > 1
        || data
            .first()
            .is_some_and(|r| r.len() != 1 || !r[0].is_finite()))
        || a.len() < 5
        || a.len() > 16384
    {
        return Err("Prediction takes one finite scalar placeholder (or empty direct data) and [horizon,confidence,model...]".into());
    }
    finite(a)?;
    let horizon = integer(a[0], 1, 500, "horizon")?;
    let z = confidence(a[1])?;
    let v = &a[2..];
    if v[0] != 1. {
        return Err("Unsupported state-space model version".into());
    }
    let k = integer(v[2], 1, 32, "states")?;
    let (means, latent, obs, kind) = if v[1] == 0. {
        if k > 8 || v.len() != 4 + 2 * k + 3 * k * k {
            return Err("Gaussian saved model length mismatch".into());
        }
        let mut at = 3;
        let m = v[at..at + k].to_vec();
        at += k;
        let p = read_mat(v, &mut at, k);
        let f = read_mat(v, &mut at, k);
        let q = read_mat(v, &mut at, k);
        let h = v[at..at + k].to_vec();
        at += k;
        let r = v[at];
        if r <= 0. {
            return Err("Positive measurement variance required".into());
        }
        validate_cov(&p, false, true)?;
        validate_cov(&q, false, true)?;
        let (m, l, o) = g_forecast(&Gaussian { m, p, f, q, h, r }, horizon)?;
        (m, l, o, "independent_gaussian")
    } else if v[1] == 1. {
        if v.len() != 6 + 3 * k + k * k || v[4] <= 0. || v[5] <= 0. {
            return Err("ETS saved model layout/scale invalid".into());
        }
        let origin = v[3];
        let scale = v[4];
        let sigma = v[5];
        let mut at = 6;
        let last = v[at..at + k].to_vec();
        at += k;
        let f = read_mat(v, &mut at, k);
        let w = &v[at..at + k];
        at += k;
        let g = &v[at..at + k];
        let (m, l, o) = ets_forecast(&last, &f, w, g, sigma, horizon)?;
        (
            m.iter().map(|v| origin + scale * v).collect(),
            l.iter()
                .map(|v| scaled_variance(*v, scale))
                .collect::<Result<Vec<_>, _>>()?,
            o.iter()
                .map(|v| scaled_variance(*v, scale))
                .collect::<Result<Vec<_>, _>>()?,
            "additive_innovations_ets",
        )
    } else {
        return Err("Unsupported state-space model kind".into());
    };
    let (lo, hi) = bounds(&means, &obs, z)?;
    Ok(format!("{{\"version\":1,\"operation\":\"state_space_predict\",\"model_kind\":\"{kind}\",\"forecast_means\":{means:?},\"forecast_latent_variances\":{latent:?},\"forecast_observation_variances\":{obs:?},\"lower\":{lo:?},\"upper\":{hi:?},\"parameter_uncertainty_included\":false}}"))
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    match operation {
        "kalman_smoother" => gaussian(data, args),
        "state_space_smoothing" => ets(data, args),
        "state_space_predict" => predict(data, args),
        _ => Err("Unknown state-space operation".into()),
    }
}
#[derive(Clone)]
pub(crate) struct Opt {
    pub x: Vec<f64>,
    pub value: f64,
    pub converged: bool,
    pub iterations: usize,
    pub projected_gradient: f64,
}
/// Deterministic bounded Nelder-Mead with a separate projected-score stopping check.
pub(crate) fn bounded_opt<F: Fn(&[f64]) -> f64>(
    f: F,
    starts: &[Vec<f64>],
    maxiter: usize,
    tol: f64,
) -> Opt {
    let d = starts[0].len();
    if d == 0 {
        return Opt {
            x: vec![],
            value: f(&[]),
            converged: true,
            iterations: 0,
            projected_gradient: 0.,
        };
    }
    let mut best = Opt {
        x: starts[0].clone(),
        value: f64::INFINITY,
        converged: false,
        iterations: 0,
        projected_gradient: f64::INFINITY,
    };
    for start in starts {
        let mut simplex = vec![start.clone()];
        for j in 0..d {
            let mut x = start.clone();
            x[j] = if x[j] < 0.8 { x[j] + 0.15 } else { x[j] - 0.15 };
            simplex.push(x);
        }
        let mut val = simplex.iter().map(|x| f(x)).collect::<Vec<_>>();
        let mut it = 0;
        let mut stopped = false;
        while it < maxiter {
            it += 1;
            let mut order = (0..=d).collect::<Vec<_>>();
            order.sort_by(|i, j| val[*i].total_cmp(&val[*j]));
            simplex = order.iter().map(|i| simplex[*i].clone()).collect();
            val = order.iter().map(|i| val[*i]).collect();
            let spread = simplex[1..]
                .iter()
                .map(|x| {
                    x.iter()
                        .zip(&simplex[0])
                        .map(|(x, y)| (x - y).abs())
                        .fold(0., f64::max)
                })
                .fold(0., f64::max);
            if val[0].is_finite()
                && spread < tol.powf(0.75).max(1e-8)
                && val[d] - val[0] < 0.01 * tol * (1. + val[0].abs())
            {
                stopped = true;
                break;
            }
            let center = (0..d)
                .map(|j| simplex[..d].iter().map(|x| x[j]).sum::<f64>() / d as f64)
                .collect::<Vec<_>>();
            let point = |factor: f64| {
                center
                    .iter()
                    .zip(&simplex[d])
                    .map(|(c, w)| (c + factor * (c - w)).clamp(0., 1.))
                    .collect::<Vec<_>>()
            };
            let xr = point(1.);
            let vr = f(&xr);
            if vr < val[0] {
                let xe = point(2.);
                let ve = f(&xe);
                if ve < vr {
                    simplex[d] = xe;
                    val[d] = ve;
                } else {
                    simplex[d] = xr;
                    val[d] = vr;
                }
            } else if vr < val[d - 1] {
                simplex[d] = xr;
                val[d] = vr;
            } else {
                let factor = if vr < val[d] { 0.5 } else { -0.5 };
                let xc = point(factor);
                let vc = f(&xc);
                if vc < val[d].min(vr) {
                    simplex[d] = xc;
                    val[d] = vc;
                } else {
                    for i in 1..=d {
                        for j in 0..d {
                            simplex[i][j] = 0.5 * (simplex[0][j] + simplex[i][j]);
                        }
                        val[i] = f(&simplex[i]);
                    }
                }
            }
        }
        let idx = (0..=d).min_by(|i, j| val[*i].total_cmp(&val[*j])).unwrap();
        let x = simplex[idx].clone();
        let v = val[idx];
        let mut pg = 0_f64;
        if v.is_finite() {
            for j in 0..d {
                let h = 1e-5;
                let mut xp = x.clone();
                let mut xm = x.clone();
                xp[j] = (xp[j] + h).min(1.);
                xm[j] = (xm[j] - h).max(0.);
                let grad = (f(&xp) - f(&xm)) / (xp[j] - xm[j]);
                let proj = if x[j] < 1e-6 {
                    grad.min(0.)
                } else if x[j] > 1. - 1e-6 {
                    grad.max(0.)
                } else {
                    grad
                };
                pg = pg.max(proj.abs());
            }
        }
        let candidate = Opt {
            x,
            value: v,
            converged: stopped && pg < tol.sqrt().max(1e-7),
            iterations: it,
            projected_gradient: pg,
        };
        if candidate.value < best.value
            || (candidate.value <= best.value + 1e-12 && candidate.converged && !best.converged)
        {
            best = candidate;
        }
    }
    best
}
fn ets_system(code: usize, period: usize, u: &[f64]) -> (Mat, Vec<f64>, Vec<f64>, Mat, Vec<f64>) {
    let trend = [1, 2, 4, 5].contains(&code);
    let seasonal = code >= 3;
    let damped = [2, 5].contains(&code);
    let mut at = 0;
    let alpha = 0.0001 + 0.9998 * u[at];
    at += 1;
    let beta = if trend {
        let b = alpha * u[at];
        at += 1;
        b
    } else {
        0.
    };
    let gamma = if seasonal {
        let g = (1. - alpha) * u[at];
        at += 1;
        g
    } else {
        0.
    };
    let phi = if damped { 0.8 + 0.2 * u[at] } else { 1. };
    let k = 1 + usize::from(trend) + if seasonal { period } else { 0 };
    let mut f = zeros(k, k);
    let mut w = vec![0.; k];
    let mut g = vec![0.; k];
    f[0][0] = 1.;
    w[0] = 1.;
    g[0] = alpha;
    if trend {
        f[0][1] = phi;
        f[1][1] = phi;
        w[1] = phi;
        g[1] = beta;
    }
    if seasonal {
        let s = 1 + usize::from(trend);
        f[s][k - 1] = 1.;
        for i in s + 1..k {
            f[i][i - 1] = 1.;
        }
        w[k - 1] = 1.;
        g[s] = gamma;
    }
    let mut b = eye(k);
    if seasonal {
        b = zeros(k, k - 1);
        for i in 0..k - period {
            b[i][i] = 1.;
        }
        for i in 0..period - 1 {
            b[k - period + i][k - period + i] = 1.;
            b[k - 1][k - period + i] = -1.;
        }
    }
    (f, w, g, b, vec![alpha, beta, gamma, phi])
}
struct EtsFit {
    sigma: f64,
    initial: Vec<f64>,
    last: Vec<f64>,
    residuals: Vec<f64>,
    fitted: Vec<f64>,
}
fn ets_profile(y: &[f64], code: usize, period: usize, u: &[f64]) -> Result<EtsFit, String> {
    let (f, w, g, b, _) = ets_system(code, period, u);
    let k = w.len();
    let mut d = b.clone();
    let mut state = vec![0.; k];
    let mut x = vec![];
    let mut z = vec![];
    let mut a = f.clone();
    for i in 0..k {
        for j in 0..k {
            a[i][j] -= g[i] * w[j];
        }
    }
    for yi in y {
        let fit = dot(&w, &state);
        z.push(yi - fit);
        x.push(
            (0..b[0].len())
                .map(|j| (0..k).map(|i| w[i] * d[i][j]).sum())
                .collect(),
        );
        state = mv(&f, &state);
        for i in 0..k {
            state[i] += g[i] * (yi - fit);
        }
        d = mul(&a, &d);
    }
    let (coef, res) = ls(&x, &z)?;
    let sigma = dot(&res, &res) / y.len() as f64;
    if !sigma.is_finite() || sigma <= 1e-24 {
        return Err("Degenerate ETS residual variance".into());
    }
    let initial = mv(&b, &coef);
    let mut last = initial.clone();
    let mut residuals = vec![];
    let mut fitted = vec![];
    for yi in y {
        let pred = dot(&w, &last);
        let e = yi - pred;
        fitted.push(pred);
        residuals.push(e);
        last = mv(&f, &last);
        for i in 0..k {
            last[i] += g[i] * e;
        }
    }
    finite(&last)?;
    Ok(EtsFit {
        sigma,
        initial,
        last,
        residuals,
        fitted,
    })
}
fn ets_forecast(
    last: &[f64],
    f: &Mat,
    w: &[f64],
    g: &[f64],
    sigma: f64,
    horizon: usize,
) -> Result<Forecast, String> {
    let mut state = last.to_vec();
    let mut impulse = g.to_vec();
    let mut cumulative = 0.;
    let mut means = vec![];
    let mut latent = vec![];
    let mut vars = vec![];
    for i in 0..horizon {
        means.push(dot(w, &state));
        latent.push(sigma * cumulative);
        vars.push(sigma * (1. + cumulative));
        state = mv(f, &state);
        if i + 1 < horizon {
            let weight = dot(w, &impulse);
            cumulative += weight * weight;
            impulse = mv(f, &impulse);
        }
    }
    finite(&means)?;
    finite(&vars)?;
    Ok((means, latent, vars))
}
fn ets(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if data.len() < 8
        || data.len() > 5000
        || data.iter().any(|r| r.len() != 1 || !r[0].is_finite())
        || a.len() > 6
    {
        return Err("ETS requires 8..5000 finite scalar observations".into());
    }
    finite(a)?;
    let code = integer(*a.first().unwrap_or(&0.), 0, 5, "ETS family")?;
    let seasonal = code >= 3;
    let period = if seasonal {
        integer(*a.get(1).unwrap_or(&4.), 2, 24, "period")?
    } else {
        integer(*a.get(1).unwrap_or(&0.), 0, 24, "unused period")?;
        0
    };
    let horizon = integer(*a.get(2).unwrap_or(&10.), 1, 500, "horizon")?;
    let conf = *a.get(3).unwrap_or(&0.95);
    let z = confidence(conf)?;
    let iterations = integer(*a.get(4).unwrap_or(&1000.), 1, 5000, "maxiter")?;
    let tol = *a.get(5).unwrap_or(&1e-8);
    if !(1e-10..=1e-3).contains(&tol) {
        return Err("Tolerance must be in [1e-10,1e-3]".into());
    }
    let trend = [1, 2, 4, 5].contains(&code);
    let damped = [2, 5].contains(&code);
    let k = 1 + usize::from(trend) + if seasonal { period } else { 0 };
    let d = 1 + usize::from(trend) + usize::from(seasonal) + usize::from(damped);
    let n = data.len();
    if n < 2 * k + 4 || n * k * k > 2_000_000 || n * k * k * iterations > 250_000_000 {
        return Err("ETS identification/work resource limit exceeded".into());
    }
    let origin = data[0][0];
    let scale = data
        .iter()
        .map(|r| (r[0] - origin).abs())
        .fold(0., f64::max);
    if !scale.is_finite() || scale == 0. {
        return Err("ETS requires representable nonconstant data".into());
    }
    let y = data
        .iter()
        .map(|r| (r[0] - origin) / scale)
        .collect::<Vec<_>>();
    let objective =
        |u: &[f64]| ets_profile(&y, code, period, u).map_or(f64::INFINITY, |p| p.sigma.ln());
    let starts = [vec![0.2; d], vec![0.5; d], vec![0.8; d]];
    let opt = bounded_opt(objective, &starts, iterations, tol);
    if !opt.value.is_finite() {
        return Err("No finite identifiable ETS fit".into());
    }
    let fit = ets_profile(&y, code, period, &opt.x)?;
    let (f, w, g, _, params) = ets_system(code, period, &opt.x);
    let sigma = scaled_variance(fit.sigma, scale)?;
    finite(&[sigma])?;
    let guard = opt.x[0] < 1e-6 || opt.x[0] > 1. - 1e-6;
    // Curvature is checked on the free smoothing face. No covariance/SE is claimed.
    let free = (0..d)
        .filter(|i| opt.x[*i] > 1e-4 && opt.x[*i] < 1. - 1e-4)
        .collect::<Vec<_>>();
    let mut information = zeros(free.len(), free.len());
    let h = 1e-4;
    for (ii, i) in free.iter().enumerate() {
        for (jj, j) in free.iter().enumerate() {
            let mut pp = opt.x.clone();
            let mut pm = opt.x.clone();
            let mut mp = opt.x.clone();
            let mut mm = opt.x.clone();
            pp[*i] += h;
            pp[*j] += h;
            pm[*i] += h;
            pm[*j] -= h;
            mp[*i] -= h;
            mp[*j] += h;
            mm[*i] -= h;
            mm[*j] -= h;
            information[ii][jj] =
                (objective(&pp) - objective(&pm) - objective(&mp) + objective(&mm)) / (4. * h * h);
        }
    }
    let identified = free.is_empty() || chol(&information).is_ok();
    let converged = opt.converged && !guard && identified;
    let ll = -0.5 * n as f64 * ((2. * std::f64::consts::PI).ln() + sigma.ln() + 1.);
    finite(&[ll])?;
    let nparams = d + k - usize::from(seasonal) + 1;
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
    let (means, latent, obs) = ets_forecast(&fit.last, &f, &w, &g, fit.sigma, horizon)?;
    let means = means.iter().map(|v| origin + scale * v).collect::<Vec<_>>();
    let latent = latent
        .iter()
        .map(|v| scaled_variance(*v, scale))
        .collect::<Result<Vec<_>, _>>()?;
    let obs = obs
        .iter()
        .map(|v| scaled_variance(*v, scale))
        .collect::<Result<Vec<_>, _>>()?;
    let (lo, hi) = bounds(&means, &obs, z)?;
    let residuals = fit.residuals.iter().map(|v| scale * v).collect::<Vec<_>>();
    let fitted = fit
        .fitted
        .iter()
        .map(|v| origin + scale * v)
        .collect::<Vec<_>>();
    let mut carrier = vec![1., 1., k as f64, origin, scale, fit.sigma];
    carrier.extend(&fit.last);
    carrier.extend(flatten(&f));
    carrier.extend(&w);
    carrier.extend(&g);
    finite(&carrier)?;
    finite(&[opt.projected_gradient])?;
    Ok(format!("{{\"version\":1,\"operation\":\"state_space_smoothing\",\"model_kind\":\"additive_innovations_ets\",\"family\":{code},\"period\":{period},\"smoothing_parameters\":{params:?},\"initial_states_standardized\":{:?},\"residuals\":{residuals:?},\"fitted\":{fitted:?},\"sigma_squared\":{sigma},\"log_likelihood\":{ll},\"parameter_count\":{nparams},\"aic\":{aic},\"bic\":{bic},\"optimizer_converged\":{},\"converged\":{converged},\"information_identified\":{identified},\"numerical_boundary_guard\":{guard},\"iterations\":{},\"projected_gradient\":{},\"forecast_means\":{means:?},\"forecast_latent_variances\":{latent:?},\"forecast_observation_variances\":{obs:?},\"lower\":{lo:?},\"upper\":{hi:?},\"prediction_model\":{carrier:?},\"parameter_uncertainty_included\":false}}",fit.initial,opt.converged,opt.iterations,opt.projected_gradient))
}
