//! Explicit CPU curve families, B-spline regression, and local polynomial smoothing.
//! Fits are ordinary unweighted residual least squares. P-splines penalize second
//! coefficient differences, not the integrated curvature used by smoothing splines.
//! No arbitrary formula language, bounds, profile likelihood, or JMP equivalence is implied.

pub fn run(op: &str, data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let w = data.first().map_or(0, Vec::len);
    if w == 0
        || data.len() > 10000
        || w > 2
        || data
            .iter()
            .any(|r| r.len() != w || r.iter().any(|x| !x.is_finite()))
        || a.len() > 16384
        || a.iter().any(|x| !x.is_finite())
    {
        return Err(
            "nonlinear methods require finite rectangular data, <=10000 rows and 2 columns".into(),
        );
    }
    let out = match op {
        "curve_fit" => curve_fit(data, a),
        "curve_predict" => curve_predict(data, a),
        "bspline_basis" => basis_operation(data, a),
        "spline_fit" => spline_fit(data, a),
        "spline_predict" => spline_predict(data, a),
        "loess" => loess(data, a),
        "loess_predict" => loess_predict(data, a),
        _ => Err("unknown nonlinear operation".into()),
    }?;
    if out
        .split([':', ',', '[', ']', '}', ' '])
        .any(|x| matches!(x, "NaN" | "inf" | "-inf"))
    {
        return Err("nonlinear result exceeds finite numeric range".into());
    }
    Ok(out)
}
fn nums(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn mat(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| nums(r)).collect::<Vec<_>>().join(",")
    )
}
fn integer(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("expected integer in [{lo},{hi}]"))
    } else {
        Ok(x as usize)
    }
}
fn width(d: &[Vec<f64>], w: usize) -> Result<(), String> {
    if d[0].len() != w {
        Err(format!("operation requires {w} columns"))
    } else {
        Ok(())
    }
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn norm(x: &[f64]) -> f64 {
    x.iter().fold(0_f64, |a, b| a.hypot(*b))
}
fn inverse(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let scales = (0..n).map(|i| a[i][i].sqrt()).collect::<Vec<_>>();
    if scales.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err("singular information matrix".into());
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v =
                a[i][j] / scales[i] / scales[j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if v <= 1e-12 || !v.is_finite() {
                    return Err("rank deficient information matrix".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut inv = vec![vec![0.; n]; n];
    for col in 0..n {
        let mut b = vec![0.; n];
        for i in 0..n {
            b[i] = ((if i == col { 1. } else { 0. }) - (0..i).map(|j| l[i][j] * b[j]).sum::<f64>())
                / l[i][i];
        }
        let mut x = vec![0.; n];
        for i in (0..n).rev() {
            x[i] = (b[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            inv[i][col] = x[i] / scales[i] / scales[col];
        }
    }
    Ok(inv)
}
fn multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    (0..a.len())
        .map(|i| {
            (0..b[0].len())
                .map(|j| (0..b.len()).map(|k| a[i][k] * b[k][j]).sum())
                .collect()
        })
        .collect()
}
fn information(j: &[Vec<f64>]) -> Vec<Vec<f64>> {
    let p = j[0].len();
    (0..p)
        .map(|i| {
            (0..p)
                .map(|k| j.iter().map(|r| r[i] * r[k]).sum())
                .collect()
        })
        .collect()
}
fn family_p(f: usize) -> usize {
    match f {
        0 | 5 | 7 => 2,
        1 | 4 | 6 | 8 | 9 => 3,
        _ => 4,
    }
}
fn family_name(f: usize) -> &'static str {
    [
        "exponential_2p",
        "exponential_3p",
        "logistic_4p",
        "gompertz_4p",
        "power_3p",
        "michaelis_menten",
        "gaussian_peak_3p",
        "logistic_2p",
        "logistic_3p",
        "gompertz_3p",
    ][f]
}
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let t = x.exp();
        t / (1. + t)
    }
}
fn evaluate(f: usize, x: f64, p: &[f64]) -> Result<(f64, Vec<f64>), String> {
    let (y, g) = match f {
        0 => {
            let e = (p[1] * x).exp();
            (p[0] * e, vec![e, p[0] * x * e])
        }
        1 => {
            let e = (p[2] * x).exp();
            (p[0] + p[1] * e, vec![1., e, p[1] * x * e])
        }
        2 | 7 | 8 => {
            let s = sigmoid(p[0] * (x - p[1]));
            let (lo, hi) = match f {
                2 => (p[2], p[3]),
                7 => (0., 1.),
                _ => (0., p[2]),
            };
            let t = (hi - lo) * s * (1. - s);
            let mut g = vec![t * (x - p[1]), -t * p[0]];
            if f == 2 {
                g.extend([1. - s, s]);
            } else if f == 8 {
                g.push(s);
            }
            (lo + (hi - lo) * s, g)
        }
        3 | 9 => {
            let (lo, hi, r, m) = if f == 3 {
                (p[0], p[1], p[2], p[3])
            } else {
                (0., p[0], p[1], p[2])
            };
            let t = (-r * (x - m)).exp();
            let e = (-t).exp();
            let et = if t.is_finite() { e * t } else { 0. };
            let u = (hi - lo) * et;
            (
                lo + (hi - lo) * e,
                if f == 3 {
                    vec![1. - e, e, u * (x - m), -u * r]
                } else {
                    vec![e, u * (x - m), -u * r]
                },
            )
        }
        4 => {
            if x <= 0. {
                return Err("power family requires positive x".into());
            }
            let e = x.powf(p[2]);
            (p[0] + p[1] * e, vec![1., e, p[1] * e * x.ln()])
        }
        5 => {
            if x <= 0. || p[1] <= 0. {
                return Err("Michaelis-Menten requires positive x and constant".into());
            }
            let t = x / (p[1] + x);
            (p[0] * t, vec![t, -p[0] * t / (p[1] + x)])
        }
        6 => {
            if p[2] <= 0. {
                return Err("Gaussian width must be positive".into());
            }
            let z = (x - p[1]) / p[2];
            let e = (-0.5 * z * z).exp();
            (
                p[0] * e,
                vec![e, p[0] * e * z / p[2], p[0] * e * z * z / p[2]],
            )
        }
        _ => unreachable!(),
    };
    if !y.is_finite() || g.iter().any(|v| !v.is_finite()) {
        Err("curve evaluation exceeds numeric range".into())
    } else {
        Ok((y, g))
    }
}
fn start(f: usize, d: &[Vec<f64>]) -> Result<Vec<f64>, String> {
    let xmin = d.iter().map(|r| r[0]).fold(f64::INFINITY, f64::min);
    let xmax = d.iter().map(|r| r[0]).fold(f64::NEG_INFINITY, f64::max);
    let ymin = d.iter().map(|r| r[1]).fold(f64::INFINITY, f64::min);
    let ymax = d.iter().map(|r| r[1]).fold(f64::NEG_INFINITY, f64::max);
    let range = xmax - xmin;
    if range <= 0. || !range.is_finite() {
        return Err("curve fit requires distinct x values".into());
    }
    let mid = xmin + range / 2.;
    let direction = if d
        .iter()
        .map(|r| (r[0] - mid) * (r[1] - (ymin + (ymax - ymin) / 2.)))
        .sum::<f64>()
        < 0.
    {
        -1.
    } else {
        1.
    };
    Ok(match f {
        0 => {
            if d.iter().all(|r| r[1] > 0.) {
                let xm = d.iter().map(|r| r[0] - xmin).sum::<f64>() / d.len() as f64;
                let ym = d.iter().map(|r| r[1].ln()).sum::<f64>() / d.len() as f64;
                let b = d
                    .iter()
                    .map(|r| (r[0] - xmin - xm) * (r[1].ln() - ym))
                    .sum::<f64>()
                    / d.iter().map(|r| (r[0] - xmin - xm).powi(2)).sum::<f64>();
                vec![(ym - b * (xmin + xm)).exp(), b]
            } else {
                vec![1., direction / range]
            }
        }
        1 => vec![
            ymin - (ymax - ymin) * 0.1,
            (ymax - ymin).max(1e-6),
            direction / range,
        ],
        2 => vec![4. * direction / range, mid, ymin, ymax],
        3 => vec![ymin, ymax, 4. * direction / range, mid],
        4 => vec![ymin, 1., 1.],
        5 => vec![2. * ymax, mid.abs().max(1e-6)],
        6 => vec![
            ymax,
            d.iter().max_by(|a, b| a[1].total_cmp(&b[1])).unwrap()[0],
            range / 4.,
        ],
        7 => vec![4. * direction / range, mid],
        8 => vec![4. * direction / range, mid, ymax],
        9 => vec![ymax, 4. * direction / range, mid],
        _ => unreachable!(),
    })
}
type CurveState = (Vec<f64>, Vec<Vec<f64>>, Vec<f64>);
fn curve_state(f: usize, d: &[Vec<f64>], p: &[f64], scale: f64) -> Result<CurveState, String> {
    let mut fit = Vec::with_capacity(d.len());
    let mut j = Vec::with_capacity(d.len());
    let mut residual = Vec::with_capacity(d.len());
    for row in d {
        let (y, g) = evaluate(f, row[0], p)?;
        fit.push(y);
        j.push(g.iter().map(|x| x / scale).collect());
        residual.push((row[1] - y) / scale);
    }
    Ok((fit, j, residual))
}
fn curve_fit(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 2)?;
    if a.is_empty() {
        return Err("curve_fit requires family".into());
    }
    let f = integer(a[0], 0, 9)?;
    let k = family_p(f);
    if d.len() <= k {
        return Err("curve fit requires n greater than parameter count".into());
    }
    let maxiter = integer(*a.get(1).unwrap_or(&300.), 1, 5000)?;
    let tol = *a.get(2).unwrap_or(&1e-9);
    if !(1e-12..=1e-3).contains(&tol) {
        return Err("curve tolerance must be in [1e-12,1e-3]".into());
    }
    let mut p = if a.len() <= 3 {
        start(f, d)?
    } else if a.len() == 3 + k {
        a[3..].to_vec()
    } else {
        return Err("initial parameters must follow family,max_iter,tolerance".into());
    };
    if d.len().saturating_mul(maxiter).saturating_mul(k * k) > 200_000_000 {
        return Err("curve work limit exceeded".into());
    }
    let scale = d
        .iter()
        .map(|r| r[1].abs())
        .fold(0_f64, f64::max)
        .max(1e-300);
    let mut state = curve_state(f, d, &p, scale)?;
    let mut loss = dot(&state.2, &state.2);
    let mut lambda = 1e-3;
    let mut converged = false;
    let mut iterations = 0;
    for it in 0..maxiter {
        iterations = it + 1;
        let scales = (0..k)
            .map(|c| state.1.iter().fold(0_f64, |v, r| v.hypot(r[c])))
            .collect::<Vec<_>>();
        if scales.iter().any(|s| *s == 0. || !s.is_finite()) {
            break;
        }
        let j = state
            .1
            .iter()
            .map(|r| {
                r.iter()
                    .zip(&scales)
                    .map(|(x, s)| x / s)
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let mut info = information(&j);
        let gradient = (0..k)
            .map(|c| j.iter().zip(&state.2).map(|(r, e)| r[c] * e).sum::<f64>())
            .collect::<Vec<_>>();
        if gradient.iter().fold(0_f64, |v, x| v.max(x.abs())) <= tol * loss.sqrt().max(1e-14) {
            converged = true;
            break;
        }
        for (i, row) in info.iter_mut().enumerate() {
            row[i] += lambda;
        }
        let inv = inverse(&info)?;
        let delta = inv.iter().map(|r| dot(r, &gradient)).collect::<Vec<_>>();
        let trial = p
            .iter()
            .zip(&delta)
            .zip(&scales)
            .map(|((v, step), s)| v + step / s)
            .collect::<Vec<_>>();
        if let Ok(newstate) = curve_state(f, d, &trial, scale) {
            let newloss = dot(&newstate.2, &newstate.2);
            if newloss < loss {
                let improvement = loss - newloss;
                p = trial;
                state = newstate;
                loss = newloss;
                lambda = (lambda / 3.).max(1e-15);
                if improvement <= tol * loss.max(1e-28)
                    || norm(&delta)
                        <= tol
                            * (1.
                                + norm(
                                    &p.iter()
                                        .zip(&scales)
                                        .map(|(x, s)| x * s)
                                        .collect::<Vec<_>>(),
                                ))
                {
                    converged = true;
                    break;
                }
                continue;
            }
        }
        lambda *= 10.;
        if lambda > 1e18 {
            break;
        }
    }
    let df = d.len() - k;
    let sse = loss * scale * scale;
    let sigma = sse / df as f64;
    // Compute information in unit-length Jacobian columns. Forming J'J in
    // response units can overflow even when parameter standard errors are finite.
    let js = (0..k)
        .map(|c| state.1.iter().fold(0_f64, |v, r| v.hypot(r[c])))
        .collect::<Vec<_>>();
    let normalized = state
        .1
        .iter()
        .map(|r| r.iter().zip(&js).map(|(x, s)| x / s).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let covariance = inverse(&information(&normalized)).ok();
    let (covstr, se, identified) = if let Some(inv) = covariance {
        let factor = (loss / df as f64).sqrt();
        let se = (0..k)
            .map(|i| inv[i][i].max(0.).sqrt() * factor / js[i])
            .collect::<Vec<_>>();
        let cov = (0..k)
            .map(|i| {
                (0..k)
                    .map(|j| (inv[i][j] * factor / js[i]) * (factor / js[j]))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        (mat(&cov), nums(&se), true)
    } else {
        ("null".into(), "null".into(), false)
    };
    let residuals = d
        .iter()
        .zip(&state.0)
        .map(|(r, y)| r[1] - y)
        .collect::<Vec<_>>();
    let mut predict = vec![f as f64];
    predict.extend(&p);
    Ok(format!("{{\"version\":1,\"algorithm\":\"nonlinear_least_squares_lm\",\"family\":\"{}\",\"coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"sse\":{},\"residual_df\":{},\"sigma_squared\":{},\"covariance\":{},\"standard_errors\":{},\"locally_identified\":{},\"covariance_method\":\"residual_variance_times_inverse_jacobian_information\",\"converged\":{},\"iterations\":{},\"predict_parameters\":{}}}",family_name(f),nums(&p),nums(&state.0),nums(&residuals),sse,df,sigma,covstr,se,identified,converged,iterations,nums(&predict)))
}
fn curve_predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 1)?;
    if a.is_empty() {
        return Err("curve predictor missing".into());
    }
    let f = integer(a[0], 0, 9)?;
    if a.len() != family_p(f) + 1 {
        return Err("wrong curve parameter count".into());
    }
    let y = d
        .iter()
        .map(|r| evaluate(f, r[0], &a[1..]).map(|v| v.0))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"curve_predict\",\"predictions\":{}}}",
        nums(&y)
    ))
}

fn knots(a: &[f64]) -> Result<(usize, &[f64], usize), String> {
    if a.len() < 5 {
        return Err("B-spline requires degree and full clamped knot vector".into());
    }
    let degree = integer(a[0], 1, 3)?;
    let t = &a[1..];
    if t.len() < 2 * (degree + 1) || t.len() > 69 || t.windows(2).any(|x| x[1] < x[0]) {
        return Err("invalid nondecreasing knot vector".into());
    }
    let p = t.len() - degree - 1;
    if p > 64
        || t[..=degree].iter().any(|x| *x != t[0])
        || t[p..].iter().any(|x| *x != t[p])
        || t[0] >= t[p]
        || !(t[p] - t[0]).is_finite()
        || t[degree + 1..p].iter().any(|x| *x <= t[0] || *x >= t[p])
    {
        return Err(
            "knots must have degree+1 repeated endpoints and interior knots strictly inside".into(),
        );
    }
    for group in t[degree + 1..p].windows(degree + 1) {
        if group[0] == group[degree] {
            return Err("interior knot multiplicity exceeds degree".into());
        }
    }
    Ok((degree, t, p))
}
fn basis(x: f64, degree: usize, t: &[f64]) -> Result<Vec<f64>, String> {
    let p = t.len() - degree - 1;
    if x < t[degree] || x > t[p] {
        return Err("B-spline evaluation lies outside knot domain".into());
    }
    if x == t[p] {
        let mut b = vec![0.; p];
        b[p - 1] = 1.;
        return Ok(b);
    }
    let mut b = (0..t.len() - 1)
        .map(|i| if t[i] <= x && x < t[i + 1] { 1. } else { 0. })
        .collect::<Vec<_>>();
    for k in 1..=degree {
        for i in 0..t.len() - k - 1 {
            let left = if t[i + k] > t[i] {
                (x - t[i]) / (t[i + k] - t[i]) * b[i]
            } else {
                0.
            };
            let right = if t[i + k + 1] > t[i + 1] {
                (t[i + k + 1] - x) / (t[i + k + 1] - t[i + 1]) * b[i + 1]
            } else {
                0.
            };
            b[i] = left + right;
        }
    }
    b.truncate(p);
    Ok(b)
}
fn basis_operation(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 1)?;
    let (degree, t, _) = knots(a)?;
    let b = d
        .iter()
        .map(|r| basis(r[0], degree, t))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"bspline_basis\",\"degree\":{},\"knots\":{},\"basis\":{}}}",
        degree,
        nums(t),
        mat(&b)
    ))
}
fn spline_fit(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 2)?;
    if a.len() < 3 || a[1] < 0. {
        return Err("spline_fit requires degree,nonnegative lambda,knots".into());
    }
    let mut spec = vec![a[0]];
    spec.extend_from_slice(&a[2..]);
    let (degree, t, p) = knots(&spec)?;
    if d.len() <= p || d.len().saturating_mul(p * p) > 100_000_000 {
        return Err("spline fit requires n>basis count and bounded work".into());
    }
    let b = d
        .iter()
        .map(|r| basis(r[0], degree, t))
        .collect::<Result<Vec<_>, _>>()?;
    let raw = information(&b);
    let mut info = raw.clone();
    for i in 0..p - 2 {
        for j in 0..3 {
            for k in 0..3 {
                info[i + j][i + k] += a[1] * [1., -2., 1.][j] * [1., -2., 1.][k];
            }
        }
    }
    let inv = inverse(&info)?;
    let origin = d[0][1];
    let y = d.iter().map(|r| r[1] - origin).collect::<Vec<_>>();
    let rhs = (0..p)
        .map(|c| b.iter().zip(&y).map(|(r, y)| r[c] * y).sum::<f64>())
        .collect::<Vec<_>>();
    let coef = inv.iter().map(|r| dot(r, &rhs)).collect::<Vec<_>>();
    let fit = b.iter().map(|r| origin + dot(r, &coef)).collect::<Vec<_>>();
    let residuals = d
        .iter()
        .zip(&fit)
        .map(|(r, y)| r[1] - y)
        .collect::<Vec<_>>();
    let sse = norm(&residuals).powi(2);
    let h = multiply(&inv, &raw);
    let edf = (0..p).map(|i| h[i][i]).sum::<f64>();
    let trace_square = (0..p)
        .map(|i| (0..p).map(|j| h[i][j] * h[j][i]).sum::<f64>())
        .sum::<f64>();
    let df = d.len() as f64 - 2. * edf + trace_square;
    if df <= 0. {
        return Err("spline variance degrees of freedom is not positive".into());
    }
    let cov = multiply(&multiply(&inv, &raw), &inv)
        .iter()
        .map(|r| r.iter().map(|x| x * sse / df).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let se = (0..p).map(|i| cov[i][i].max(0.).sqrt()).collect::<Vec<_>>();
    let mut predict = vec![degree as f64, t.len() as f64, origin];
    predict.extend(t);
    predict.extend(&coef);
    Ok(format!("{{\"version\":1,\"algorithm\":\"penalized_bspline_regression\",\"penalty\":\"second_coefficient_difference_squared\",\"degree\":{},\"lambda\":{},\"knots\":{},\"response_origin\":{},\"coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"sse\":{},\"effective_df\":{},\"variance_df\":{},\"covariance\":{},\"standard_errors\":{},\"predict_parameters\":{}}}",degree,a[1],nums(t),origin,nums(&coef),nums(&fit),nums(&residuals),sse,edf,df,mat(&cov),nums(&se),nums(&predict)))
}
fn spline_predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 1)?;
    if a.len() < 5 {
        return Err("saved spline predictor missing".into());
    }
    let nt = integer(a[1], 4, 69)?;
    if a.len() < 3 + nt {
        return Err("truncated spline knots".into());
    }
    let mut spec = vec![a[0]];
    spec.extend_from_slice(&a[3..3 + nt]);
    let (degree, t, p) = knots(&spec)?;
    if a.len() != 3 + nt + p {
        return Err("wrong spline coefficient count".into());
    }
    let fit = d
        .iter()
        .map(|r| basis(r[0], degree, t).map(|b| a[2] + dot(&b, &a[3 + nt..])))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"spline_predict\",\"predictions\":{}}}",
        nums(&fit)
    ))
}

fn local(x: f64, d: &[Vec<f64>], robust: &[f64], span: f64, degree: usize) -> Result<f64, String> {
    let k = ((span * d.len() as f64).ceil() as usize)
        .max(degree + 2)
        .min(d.len());
    let mut distances = d.iter().map(|r| (r[0] - x).abs()).collect::<Vec<_>>();
    distances.sort_by(f64::total_cmp);
    let bandwidth = distances[k - 1];
    if bandwidth <= 0. || !bandwidth.is_finite() {
        return Err("local regression has no positive bandwidth".into());
    }
    let mut info = vec![vec![0.; degree + 1]; degree + 1];
    let mut rhs = vec![0.; degree + 1];
    let origin = d[0][1];
    for (row, rw) in d.iter().zip(robust) {
        let u = (row[0] - x) / bandwidth;
        if u.abs() >= 1. {
            continue;
        }
        let weight = (1. - u.abs().powi(3)).powi(3) * rw;
        let b = [1., u, u * u];
        for i in 0..=degree {
            rhs[i] += weight * b[i] * (row[1] - origin);
            for j in 0..=degree {
                info[i][j] += weight * b[i] * b[j];
            }
        }
    }
    let inv = inverse(&info).map_err(|_| {
        "local polynomial is rank deficient; increase span or reduce degree".to_string()
    })?;
    Ok(origin + dot(&inv[0], &rhs))
}
fn loess(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 2)?;
    if a.is_empty() || a.len() > 3 || a[0] <= 0. || a[0] > 1. || d.len() > 2000 {
        return Err("LOESS requires span in (0,1] and <=2000 training rows".into());
    }
    let degree = integer(*a.get(1).unwrap_or(&1.), 0, 2)?;
    let iterations = integer(*a.get(2).unwrap_or(&0.), 0, 10)?;
    if d.len() < degree + 3
        || d.len()
            .saturating_mul(d.len())
            .saturating_mul(iterations + 1)
            > 40_000_000
    {
        return Err("LOESS training count/work limit".into());
    }
    let mut robust = vec![1.; d.len()];
    let mut fit = vec![];
    let mut used = 0;
    for it in 0..=iterations {
        fit = d
            .iter()
            .map(|r| local(r[0], d, &robust, a[0], degree))
            .collect::<Result<Vec<_>, _>>()?;
        used = it;
        if it == iterations {
            break;
        }
        let mut residual = d
            .iter()
            .zip(&fit)
            .map(|(r, y)| (r[1] - y).abs())
            .collect::<Vec<_>>();
        residual.sort_by(f64::total_cmp);
        let median = if residual.len() % 2 == 0 {
            (residual[residual.len() / 2 - 1] + residual[residual.len() / 2]) / 2.
        } else {
            residual[residual.len() / 2]
        };
        if median == 0. {
            break;
        }
        for ((rw, row), y) in robust.iter_mut().zip(d).zip(&fit) {
            let u = (row[1] - y).abs() / (6. * median);
            *rw = if u >= 1. { 0. } else { (1. - u * u).powi(2) };
        }
    }
    let residual = d
        .iter()
        .zip(&fit)
        .map(|(r, y)| r[1] - y)
        .collect::<Vec<_>>();
    let mut predict = vec![a[0], degree as f64, d.len() as f64];
    predict.extend(d.iter().map(|r| r[0]));
    predict.extend(d.iter().map(|r| r[1]));
    predict.extend(&robust);
    Ok(format!("{{\"version\":1,\"algorithm\":\"loess_tricube\",\"span\":{},\"degree\":{},\"robust_iterations_used\":{},\"robust_weights\":{},\"fitted_values\":{},\"residuals\":{},\"sse\":{},\"predict_parameters\":{}}}",a[0],degree,used,nums(&robust),nums(&fit),nums(&residual),norm(&residual).powi(2),nums(&predict)))
}
fn loess_predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    width(d, 1)?;
    if a.len() < 3 || a[0] <= 0. || a[0] > 1. {
        return Err("LOESS saved state missing".into());
    }
    let degree = integer(a[1], 0, 2)?;
    let n = integer(a[2], degree + 3, 2000)?;
    if a.len() != 3 + 3 * n
        || d.len().saturating_mul(n) > 20_000_000
        || a[3 + 2 * n..].iter().any(|v| *v < 0. || *v > 1.)
    {
        return Err("invalid LOESS saved state or work limit".into());
    }
    let train = (0..n)
        .map(|i| vec![a[3 + i], a[3 + n + i]])
        .collect::<Vec<_>>();
    let y = d
        .iter()
        .map(|r| local(r[0], &train, &a[3 + 2 * n..], a[0], degree))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"loess_predict\",\"predictions\":{}}}",
        nums(&y)
    ))
}
