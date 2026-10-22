//! Profile Gaussian ML/REML spatial regression; covariance estimation is distinct
//! from the supplied-covariance operations in `spatial`.
//! Physical variance components are nonnegative. Interior observed-information
//! Wald intervals are explicit approximations, not JMP's bounded Satterthwaite CI.
type Matrix = Vec<Vec<f64>>;
const NMAX: usize = 96;
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("invalid spatial estimation integer".into())
    }
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("spatial estimation result is not representable; rescale inputs".into())
    }
}
fn vector(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn matrix(x: &Matrix) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| vector(r)).collect::<Vec<_>>().join(",")
    )
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn distance(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).fold(0_f64, |a, (x, y)| a.hypot(x - y))
}
fn cholesky(a: &Matrix) -> Result<Matrix, String> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v = a[i][j] - dot(&l[i][..j], &l[j][..j]);
            if i == j {
                if !v.is_finite() || v <= a[i][i].abs() * 1e-12 {
                    return Err("spatial covariance/information is numerically singular".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    Ok(l)
}
fn forward(l: &Matrix, y: &[f64]) -> Vec<f64> {
    let mut x = vec![0.; y.len()];
    for i in 0..y.len() {
        x[i] = (y[i] - dot(&l[i][..i], &x[..i])) / l[i][i];
    }
    x
}
fn solve(l: &Matrix, y: &[f64]) -> Vec<f64> {
    let mut x = forward(l, y);
    for i in (0..y.len()).rev() {
        x[i] = (x[i] - (i + 1..y.len()).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
    }
    x
}
fn inverse(l: &Matrix) -> Matrix {
    let n = l.len();
    let mut a = vec![vec![0.; n]; n];
    for j in 0..n {
        let mut e = vec![0.; n];
        e[j] = 1.;
        for (i, v) in solve(l, &e).into_iter().enumerate() {
            a[i][j] = v;
        }
    }
    for i in 0..n {
        let (before, after) = a.split_at_mut(i);
        let row = &mut after[0];
        for (j, prior) in before.iter_mut().enumerate() {
            let v = (row[j] + prior[i]) * 0.5;
            row[j] = v;
            prior[i] = v;
        }
    }
    a
}
fn correlation(kind: usize, h: f64, range: f64) -> f64 {
    let u = h / range;
    match kind {
        0 if u >= 1. => 0.,
        0 => (1. - u).powi(2) * (1. + 0.5 * u),
        2 => (-u * u).exp(),
        _ => (-u).exp(),
    }
}
struct Design {
    n: usize,
    p: usize,
    x: Matrix,
    y: Vec<f64>,
    points: Matrix,
    subjects: Vec<f64>,
    origins: Vec<f64>,
    centers: Vec<f64>,
    scales: Vec<f64>,
    yo: f64,
    ys: f64,
    distances: Matrix,
    min_distance: f64,
    max_distance: f64,
}
impl Design {
    fn new(data: &Matrix, d: usize, p: usize) -> Result<Self, String> {
        let n = data.len();
        let q = p + 1;
        if d == 0 || n <= q || n > NMAX || data[0].len() != d + p + 2 {
            return Err("spatial fit requires n>p+1,n<=96; rows[coordinates,fixed features,subject,response]".into());
        }
        let points = data.iter().map(|r| r[..d].to_vec()).collect::<Matrix>();
        let subjects = data.iter().map(|r| r[d + p]).collect::<Vec<_>>();
        if subjects
            .iter()
            .any(|v| v.fract() != 0. || v.abs() > 9_007_199_254_740_991.)
        {
            return Err("spatial subjects must be exact integer labels".into());
        }
        let mut distances = vec![vec![0.; n]; n];
        let mut min_distance = f64::INFINITY;
        let mut max_distance = 0_f64;
        for i in 0..n {
            for j in 0..i {
                let h = distance(&points[i], &points[j]);
                distances[i][j] = h;
                distances[j][i] = h;
                if subjects[i] == subjects[j] && h > 0. {
                    min_distance = min_distance.min(h);
                    max_distance = max_distance.max(h);
                }
            }
        }
        let origins = data[0][d..d + p].to_vec();
        let mut centers = Vec::new();
        let mut scales = Vec::new();
        for j in 0..p {
            let scale = data
                .iter()
                .map(|r| (r[d + j] - origins[j]).abs())
                .fold(0_f64, f64::max);
            if scale == 0. {
                return Err("spatial fixed feature is aliased with intercept".into());
            }
            centers.push(
                data.iter()
                    .map(|r| (r[d + j] - origins[j]) / scale / n as f64)
                    .sum(),
            );
            scales.push(scale);
        }
        let x = data
            .iter()
            .map(|r| {
                std::iter::once(1.)
                    .chain((0..p).map(|j| (r[d + j] - origins[j]) / scales[j] - centers[j]))
                    .collect()
            })
            .collect::<Matrix>();
        let yo = data[0][d + p + 1];
        let ys = data
            .iter()
            .map(|r| (r[d + p + 1] - yo).abs())
            .fold(0_f64, f64::max);
        if ys == 0. {
            return Err(
                "estimated spatial variance is zero; no positive finite variance MLE".into(),
            );
        }
        let y = data.iter().map(|r| (r[d + p + 1] - yo) / ys).collect();
        Ok(Self {
            n,
            p,
            x,
            y,
            points,
            subjects,
            origins,
            centers,
            scales,
            yo,
            ys,
            distances,
            min_distance,
            max_distance,
        })
    }
    fn covariance(&self, kind: usize, range: f64, w: f64) -> Matrix {
        (0..self.n)
            .map(|i| {
                (0..self.n)
                    .map(|j| {
                        if i == j {
                            1.
                        } else if self.subjects[i] != self.subjects[j] {
                            0.
                        } else {
                            (1. - w) * correlation(kind, self.distances[i][j], range)
                        }
                    })
                    .collect()
            })
            .collect()
    }
}
struct Profile {
    nll: f64,
    sse: f64,
    beta: Vec<f64>,
    gi: Matrix,
    l: Matrix,
    logdet: f64,
    logdet_x: f64,
}
fn profile(
    design: &Design,
    kind: usize,
    range: f64,
    w: f64,
    reml: bool,
) -> Result<Profile, String> {
    let n = design.n;
    let q = design.p + 1;
    let l = cholesky(&design.covariance(kind, range, w))?;
    let wy = forward(&l, &design.y);
    let cols = (0..q)
        .map(|j| forward(&l, &design.x.iter().map(|r| r[j]).collect::<Vec<_>>()))
        .collect::<Matrix>();
    // Twice-reorthogonalized QR avoids squaring the fixed-design condition.
    let mut axes: Matrix = Vec::new();
    let mut r = vec![vec![0.; q]; q];
    for j in 0..q {
        let mut v = cols[j].clone();
        let original = dot(&v, &v).sqrt();
        for _ in 0..2 {
            for k in 0..j {
                let a = dot(&axes[k], &v);
                r[k][j] += a;
                for i in 0..n {
                    v[i] -= a * axes[k][i];
                }
            }
        }
        let norm = dot(&v, &v).sqrt();
        if !norm.is_finite() || norm <= 1e-11 * original {
            return Err("spatial fixed design is rank deficient".into());
        }
        r[j][j] = norm;
        for a in &mut v {
            *a /= norm;
        }
        axes.push(v);
    }
    let mut beta = axes.iter().map(|v| dot(v, &wy)).collect::<Vec<_>>();
    for i in (0..q).rev() {
        beta[i] = (beta[i] - (i + 1..q).map(|j| r[i][j] * beta[j]).sum::<f64>()) / r[i][i];
    }
    let residual = wy
        .iter()
        .enumerate()
        .map(|(i, y)| y - (0..q).map(|j| cols[j][i] * beta[j]).sum::<f64>())
        .collect::<Vec<_>>();
    let sse = dot(&residual, &residual);
    if !sse.is_finite() || sse <= 128. * f64::EPSILON.powi(2) {
        return Err("spatial residual variance is below numerical identification precision".into());
    }
    let mut inverse_columns = vec![vec![0.; q]; q];
    for (j, column) in inverse_columns.iter_mut().enumerate() {
        for i in (0..q).rev() {
            column[i] = ((if i == j { 1. } else { 0. })
                - (i + 1..q).map(|k| r[i][k] * column[k]).sum::<f64>())
                / r[i][i];
        }
    }
    let ri = (0..q)
        .map(|i| inverse_columns.iter().map(|column| column[i]).collect())
        .collect::<Matrix>();
    let gi = ri
        .iter()
        .map(|a| ri.iter().map(|b| dot(a, b)).collect())
        .collect();
    let logdet = 2. * l.iter().enumerate().map(|(i, r)| r[i].ln()).sum::<f64>();
    let logdet_x = 2. * r.iter().enumerate().map(|(i, r)| r[i].ln()).sum::<f64>();
    let df = n - if reml { q } else { 0 };
    let nll = 0.5
        * (df as f64 * ((2. * std::f64::consts::PI).ln() + 1. + (sse / df as f64).ln())
            + logdet
            + if reml { logdet_x } else { 0. });
    Ok(Profile {
        nll: finite(nll)?,
        sse,
        beta,
        gi,
        l,
        logdet,
        logdet_x,
    })
}
#[derive(Clone)]
struct Candidate {
    x: Vec<f64>,
    value: f64,
    converged: bool,
    iterations: usize,
}
fn simplex<F: FnMut(&[f64]) -> f64>(
    start: Vec<f64>,
    iterations: usize,
    tol: f64,
    eval: &mut F,
) -> Candidate {
    let dim = start.len();
    let mut points = vec![start.clone()];
    for i in 0..dim {
        let mut p = start.clone();
        p[i] = if p[i] > 0.9 { p[i] - 0.05 } else { p[i] + 0.05 };
        points.push(p);
    }
    let mut values = points.iter().map(|p| eval(p)).collect::<Vec<_>>();
    let mut used = 0;
    let mut done = false;
    for it in 0..iterations {
        used = it + 1;
        let mut order = (0..=dim).collect::<Vec<_>>();
        order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
        points = order.iter().map(|&i| points[i].clone()).collect();
        values = order.iter().map(|&i| values[i]).collect();
        let width = points
            .iter()
            .skip(1)
            .flat_map(|p| p.iter().zip(&points[0]).map(|(a, b)| (a - b).abs()))
            .fold(0_f64, f64::max);
        if width < tol.sqrt() * 0.01
            && (values[dim] - values[0]).abs() < tol * (1. + values[0].abs())
        {
            done = true;
            break;
        }
        let center = (0..dim)
            .map(|j| points[..dim].iter().map(|p| p[j] / dim as f64).sum::<f64>())
            .collect::<Vec<_>>();
        let point = |factor: f64| {
            (0..dim)
                .map(|j| (center[j] + factor * (center[j] - points[dim][j])).clamp(0., 1.))
                .collect::<Vec<_>>()
        };
        let reflected = point(1.);
        let rv = eval(&reflected);
        let replacement = if rv < values[0] {
            let expanded = point(2.);
            let ev = eval(&expanded);
            if ev < rv {
                Some((expanded, ev))
            } else {
                Some((reflected, rv))
            }
        } else if rv < values[dim - 1] {
            Some((reflected, rv))
        } else {
            let contracted = if rv < values[dim] {
                point(0.5)
            } else {
                point(-0.5)
            };
            let cv = eval(&contracted);
            if cv < values[dim].min(rv) {
                Some((contracted, cv))
            } else {
                None
            }
        };
        if let Some((p, v)) = replacement {
            points[dim] = p;
            values[dim] = v;
        } else {
            let (best, others) = points.split_first_mut().unwrap();
            for (point, value) in others.iter_mut().zip(values.iter_mut().skip(1)) {
                for (coordinate, origin) in point.iter_mut().zip(best.iter()) {
                    *coordinate = (*coordinate + origin) * 0.5;
                }
                *value = eval(point);
            }
        }
    }
    let best = (0..=dim)
        .min_by(|&a, &b| values[a].total_cmp(&values[b]))
        .unwrap();
    Candidate {
        x: points[best].clone(),
        value: values[best],
        converged: done,
        iterations: used,
    }
}
fn hessian<F: Fn(&[f64]) -> Result<f64, String>>(x: &[f64], f: F) -> Result<Matrix, String> {
    let n = x.len();
    let mut a = vec![vec![0.; n]; n];
    let center = f(x)?;
    let h = 2e-4;
    for i in 0..n {
        let mut p = x.to_vec();
        let mut m = x.to_vec();
        p[i] += h;
        m[i] -= h;
        a[i][i] = (f(&p)? - 2. * center + f(&m)?) / (h * h);
        for j in 0..i {
            let mut pp = x.to_vec();
            let mut pm = x.to_vec();
            let mut mp = x.to_vec();
            let mut mm = x.to_vec();
            pp[i] += h;
            pp[j] += h;
            pm[i] += h;
            pm[j] -= h;
            mp[i] -= h;
            mp[j] += h;
            mm[i] -= h;
            mm[j] -= h;
            a[i][j] = (f(&pp)? - f(&pm)? - f(&mp)? + f(&mm)?) / (4. * h * h);
            a[j][i] = a[i][j];
        }
    }
    if a.iter().flatten().any(|x| !x.is_finite()) {
        return Err("observed information is not representable".into());
    }
    Ok(a)
}
fn raw_coefficients(d: &Design, p: &Profile) -> Result<(f64, Vec<f64>), String> {
    let b = (0..d.p)
        .map(|j| finite(p.beta[j + 1] * d.ys / d.scales[j]))
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        d.yo + p.beta[0] * d.ys
            - (0..d.p)
                .map(|j| b[j] * (d.origins[j] + d.centers[j] * d.scales[j]))
                .sum::<f64>(),
    )?;
    Ok((intercept, b))
}
fn beta_covariance(d: &Design, p: &Profile, total: f64) -> Result<Matrix, String> {
    let q = d.p + 1;
    let mut t = vec![vec![0.; q]; q];
    t[0][0] = 1.;
    for j in 0..d.p {
        t[0][j + 1] = -(d.origins[j] / d.scales[j] + d.centers[j]);
        t[j + 1][j + 1] = 1. / d.scales[j];
    }
    let scale = total.sqrt();
    for r in &mut t {
        for a in r {
            *a *= scale;
        }
    }
    t.iter()
        .map(|a| {
            t.iter()
                .map(|b| finite((0..q).map(|k| a[k] * dot(&p.gi[k], b)).sum()))
                .collect()
        })
        .collect()
}
fn fit(data: &Matrix, args: &[f64]) -> Result<String, String> {
    if !(5..=9).contains(&args.len()) || args.len() == 8 {
        return Err("spatial fit args [dimensions,fixed_features,kind,reml,nugget,maxiter=1000,tol=1e-8,optional_range_lower,range_upper]".into());
    }
    let d = integer(args[0], 3)?;
    let p = integer(args[1], 8)?;
    let kind = integer(args[2], 3)?;
    let reml = integer(args[3], 1)? == 1;
    let nugget = integer(args[4], 1)? == 1;
    let iterations = integer(args.get(5).copied().unwrap_or(1000.), 3000)?;
    let tol = args.get(6).copied().unwrap_or(1e-8);
    if iterations == 0 || !(1e-11..=1e-4).contains(&tol) {
        return Err("invalid spatial iteration limit or tolerance".into());
    }
    let design = Design::new(data, d, p)?;
    let q = p + 1;
    let df = design.n - if reml { q } else { 0 };
    if design.n as u128 * design.n as u128 * design.n as u128 * iterations as u128 > 3_000_000_000 {
        return Err("spatial likelihood requested work exceeds limit".into());
    }
    let lower = args
        .get(7)
        .copied()
        .unwrap_or(if design.min_distance.is_finite() {
            design.min_distance * 0.001
        } else {
            0.001
        });
    let upper = args.get(8).copied().unwrap_or(if design.max_distance > 0. {
        design.max_distance * 1000.
    } else {
        1000.
    });
    if lower <= 0. || upper <= lower || lower < 1e-100 || upper > 1e100 {
        return Err(
            "spatial correlation-length bounds must satisfy 1e-100<=lower<upper<=1e100".into(),
        );
    }
    let log_lower = lower.ln();
    let log_width = upper.ln() - log_lower;
    let mut evaluations = 0usize;
    // Charge every coarse, multistart, boundary-face and score evaluation. A
    // conservative reserve of100 dense likelihood evaluations covers the final
    // fit and the at-most19 observed-information perturbations outside this
    // closure. Exhaustion rejects the fit instead of returning a partial search.
    let dense_work = design.n as u128 * design.n as u128 * design.n as u128;
    let max_evaluations = (3_000_000_000_u128 / dense_work).saturating_sub(100) as usize;
    let mut budget_exhausted = false;
    let mut objective = |x: &[f64]| {
        if evaluations >= max_evaluations {
            budget_exhausted = true;
            return f64::INFINITY;
        }
        evaluations += 1;
        let range = (log_lower + x[0] * log_width).exp();
        let w = if nugget { x[1] } else { 0. };
        profile(&design, kind, range, w, reml)
            .map(|z| z.nll)
            .unwrap_or(f64::INFINITY)
    };
    let mut starts = Vec::new();
    for i in 0..=8 {
        for &w in if nugget {
            &[0., 0.15, 0.5, 0.85, 1.][..]
        } else {
            &[0.][..]
        } {
            let x = if nugget {
                vec![i as f64 / 8., w]
            } else {
                vec![i as f64 / 8.]
            };
            let value = objective(&x);
            starts.push(Candidate {
                x,
                value,
                converged: false,
                iterations: 0,
            });
        }
    }
    starts.sort_by(|a, b| a.value.total_cmp(&b.value));
    if !starts[0].value.is_finite() {
        return Err(
            "no positive definite identifiable spatial likelihood within range bounds".into(),
        );
    }
    let mut best = starts[0].clone();
    for start in starts.into_iter().take(6) {
        let candidate = simplex(start.x, iterations, tol, &mut objective);
        if candidate.value < best.value || (candidate.value == best.value && candidate.converged) {
            best = candidate;
        }
    }
    // Explicit zero-nugget face keeps an attainable variance-component boundary.
    if nugget {
        let mut face = |x: &[f64]| objective(&[x[0], 0.]);
        for i in [0.1, 0.5, 0.9] {
            let mut candidate = simplex(vec![i], iterations, tol, &mut face);
            candidate.x.push(0.);
            if candidate.value < best.value {
                best = candidate;
            }
        }
    }
    // Verify first-order stationarity independently of simplex collapse at a
    // clipped face. Scores use the bounded unit coordinates of the search.
    let mut projected_score_norm = 0_f64;
    for j in 0..best.x.len() {
        let h = 1e-5;
        let mut plus = best.x.clone();
        let mut minus = best.x.clone();
        plus[j] = (plus[j] + h).min(1.);
        minus[j] = (minus[j] - h).max(0.);
        let score = if best.x[j] < h {
            ((objective(&plus) - best.value) / (plus[j] - best.x[j])).min(0.)
        } else if best.x[j] > 1. - h {
            ((best.value - objective(&minus)) / (best.x[j] - minus[j])).max(0.)
        } else {
            (objective(&plus) - objective(&minus)) / (plus[j] - minus[j])
        };
        projected_score_norm = projected_score_norm.max(score.abs() / (1. + best.value.abs()));
    }
    let score_passed = projected_score_norm.is_finite() && projected_score_norm < 1e-4;
    if budget_exhausted {
        return Err(
            "spatial likelihood evaluation budget exhausted; reduce rows or iteration limit".into(),
        );
    }
    if !projected_score_norm.is_finite() {
        projected_score_norm = f64::MAX;
    }
    let range = (log_lower + best.x[0] * log_width).exp();
    let w = if nugget { best.x[1] } else { 0. };
    let result = profile(&design, kind, range, w, reml)?;
    let variance_normalized = result.sse / df as f64;
    let total = finite((variance_normalized.sqrt() * design.ys).powi(2))?;
    if total <= 0. {
        return Err("estimated variance underflows binary64".into());
    }
    let sill = finite(total * (1. - w))?;
    let nugget_variance = finite(total * w)?;
    let ll = finite(
        -result.nll
            - df as f64 * design.ys.ln()
            - if reml {
                design.scales.iter().map(|s| s.ln()).sum::<f64>()
            } else {
                0.
            },
    )?;
    let (intercept, beta) = raw_coefficients(&design, &result)?;
    let beta_cov = beta_covariance(&design, &result, total)?;
    let mut params = vec![
        1.,
        d as f64,
        p as f64,
        design.n as f64,
        kind as f64,
        sill,
        nugget_variance,
        range,
    ];
    params.extend(data.iter().flatten());
    let range_guard = !(1e-6..=1. - 1e-6).contains(&best.x[0]);
    let variance_boundary = nugget && !(1e-6..=1. - 1e-6).contains(&w);
    let mut theta = vec![variance_normalized.ln(), range.ln()];
    if nugget {
        theta.push(w);
    }
    let interior = !range_guard && !variance_boundary && best.converged && score_passed;
    let info = if interior {
        hessian(&theta, |t| {
            let weight = if nugget { t[2] } else { 0. };
            if !(0. ..=1.).contains(&weight) {
                return Err("information perturbation exceeds variance domain".into());
            }
            let z = profile(&design, kind, t[1].exp(), weight, reml)?;
            finite(
                0.5 * (df as f64 * ((2. * std::f64::consts::PI).ln() + t[0])
                    + z.logdet
                    + if reml { z.logdet_x } else { 0. }
                    + z.sse / t[0].exp()),
            )
        })
        .ok()
    } else {
        None
    };
    let info_factor = info.as_ref().and_then(|v| cholesky(v).ok());
    let identified = info_factor.is_some();
    let (cov_json, se_json, ci_json, info_json) =
        if let (Some(information), Some(factor)) = (info.as_ref(), info_factor) {
            let inv = inverse(&factor);
            let dim = theta.len();
            let mut jac = vec![vec![0.; dim]; if nugget { 3 } else { 2 }];
            jac[0][0] = sill;
            if nugget {
                jac[0][2] = -total;
                jac[1][0] = nugget_variance;
                jac[1][2] = total;
                jac[2][1] = range;
            } else {
                jac[1][1] = range;
            }
            let cov = jac
                .iter()
                .map(|a| {
                    jac.iter()
                        .map(|b| finite((0..dim).map(|k| a[k] * dot(&inv[k], b)).sum()))
                        .collect::<Result<Vec<_>, _>>()
                })
                .collect::<Result<Matrix, _>>()?;
            let se = cov
                .iter()
                .enumerate()
                .map(|(i, r)| r[i].max(0.).sqrt())
                .collect::<Vec<_>>();
            let estimates = if nugget {
                vec![sill, nugget_variance, range]
            } else {
                vec![sill, range]
            };
            let z = 1.959963984540054;
            let ci = estimates
                .iter()
                .zip(&se)
                .map(|(v, s)| vec![v - z * s, v + z * s])
                .collect::<Matrix>();
            (matrix(&cov), vector(&se), matrix(&ci), matrix(information))
        } else {
            (
                "null".into(),
                "null".into(),
                "null".into(),
                info.as_ref().map(matrix).unwrap_or_else(|| "null".into()),
            )
        };
    let marginal = design
        .x
        .iter()
        .map(|x| finite(design.yo + design.ys * dot(x, &result.beta)))
        .collect::<Result<Vec<_>, _>>()?;
    let residual = design
        .y
        .iter()
        .zip(&design.x)
        .map(|(y, x)| y - dot(x, &result.beta))
        .collect::<Vec<_>>();
    let alpha = solve(&result.l, &residual);
    let conditional = (0..design.n)
        .map(|i| {
            finite(
                marginal[i]
                    + design.ys
                        * (1. - w)
                        * (0..design.n)
                            .filter(|&j| design.subjects[i] == design.subjects[j])
                            .map(|j| correlation(kind, design.distances[i][j], range) * alpha[j])
                            .sum::<f64>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let kernel_parameter = if kind == 3 {
        (-1. / range).exp()
    } else {
        range
    };
    let names = if nugget {
        "[\"partial_sill\",\"nugget\",\"correlation_length\"]"
    } else {
        "[\"partial_sill\",\"correlation_length\"]"
    };
    Ok(format!("{{\"version\":1,\"method\":\"spatial_fit\",\"estimation\":\"{}\",\"intercept\":{intercept},\"coefficients\":{},\"coefficient_covariance\":{},\"total_variance\":{total},\"partial_sill\":{sill},\"nugget\":{nugget_variance},\"nugget_fraction\":{w},\"range_parameter\":{range},\"correlation_length\":{range},\"kernel_parameter\":{kernel_parameter},\"log_likelihood\":{ll},\"residual_df\":{df},\"marginal_fitted\":{},\"conditional_fitted\":{},\"covariance_parameter_names\":{names},\"covariance_parameter_covariance\":{cov_json},\"covariance_parameter_standard_errors\":{se_json},\"covariance_parameter_intervals\":{ci_json},\"observed_information\":{info_json},\"observed_information_parameters\":\"log_total_variance_log_correlation_length_optional_identity_nugget_fraction\",\"ci_method\":\"interior_delta_method_untruncated_normal_Wald_95pct_not_JMP_bounded_Satterthwaite\",\"covariance_parameters_identified\":{identified},\"variance_boundary\":{variance_boundary},\"range_bound_active\":{range_guard},\"projected_score_norm\":{projected_score_norm},\"projected_score_passed\":{score_passed},\"optimization_converged\":{},\"converged\":{},\"iterations\":{},\"likelihood_evaluations\":{evaluations},\"predict_parameters\":{},\"search\":\"deterministic_multistart_profile_Nelder_Mead_no_global_optimum_certificate\",\"prediction_uncertainty\":\"plug_in_covariance_with_fixed_effect_uncertainty_no_covariance_estimation_uncertainty\"}}",if reml{"REML"}else{"ML"},vector(&beta),matrix(&beta_cov),vector(&marginal),vector(&conditional),best.converged&&score_passed,best.converged&&score_passed&&identified&&!range_guard,best.iterations,vector(&params)))
}
fn predict(data: &Matrix, args: &[f64]) -> Result<String, String> {
    if args.len() < 8 || args[0] != 1. {
        return Err("invalid saved estimated spatial model version".into());
    }
    let d = integer(args[1], 3)?;
    let p = integer(args[2], 8)?;
    let n = integer(args[3], NMAX)?;
    let kind = integer(args[4], 3)?;
    let sill = args[5];
    let nugget = args[6];
    let range = args[7];
    let total = sill + nugget;
    if d == 0
        || n <= p + 1
        || args.len() != 8 + n * (d + p + 2)
        || data[0].len() != d + p + 1
        || data.len() > 1000
        || sill < 0.
        || nugget < 0.
        || total <= 0.
        || !total.is_finite()
        || range <= 0.
        || range > 1e100
    {
        return Err("invalid saved spatial model dimensions or covariance".into());
    }
    let training = args[8..]
        .chunks_exact(d + p + 2)
        .map(<[f64]>::to_vec)
        .collect();
    let design = Design::new(&training, d, p)?;
    let w = nugget / total;
    let result = profile(&design, kind, range, w, false)?;
    let residual = design
        .y
        .iter()
        .zip(&design.x)
        .map(|(y, x)| y - dot(x, &result.beta))
        .collect::<Vec<_>>();
    let alpha = solve(&result.l, &residual);
    let mut means = Vec::new();
    let mut latent = Vec::new();
    let mut observation = Vec::new();
    for row in data {
        let subject = row[d + p];
        if subject.fract() != 0. || subject.abs() > 9_007_199_254_740_991. {
            return Err("invalid query subject integer".into());
        }
        let x =
            std::iter::once(1.)
                .chain((0..p).map(|j| {
                    (row[d + j] - design.origins[j]) / design.scales[j] - design.centers[j]
                }))
                .collect::<Vec<_>>();
        let k = (0..n)
            .map(|i| {
                if subject == design.subjects[i] {
                    (1. - w) * correlation(kind, distance(&row[..d], &design.points[i]), range)
                } else {
                    0.
                }
            })
            .collect::<Vec<_>>();
        means.push(finite(
            design.yo + design.ys * (dot(&x, &result.beta) + dot(&k, &alpha)),
        )?);
        let a = solve(&result.l, &k);
        let v = forward(&result.l, &k);
        let delta = (0..p + 1)
            .map(|j| x[j] - (0..n).map(|i| design.x[i][j] * a[i]).sum::<f64>())
            .collect::<Vec<_>>();
        let correction = (0..p + 1)
            .map(|j| delta[j] * dot(&result.gi[j], &delta))
            .sum::<f64>();
        let removed = dot(&v, &v);
        let variance = finite(total * ((1. - w) - removed + correction))?;
        if variance / total < -1e-9 * (1. + removed.abs() + correction.abs()) {
            return Err("negative plug-in spatial predictive variance".into());
        }
        let exact = nugget == 0.
            && (0..n).any(|i| {
                subject == design.subjects[i] && row[..d] == design.points[i] && x == design.x[i]
            });
        let variance = if exact { 0. } else { variance.max(0.) };
        latent.push(variance);
        observation.push(finite(variance + nugget)?);
    }
    Ok(format!("{{\"version\":1,\"method\":\"spatial_fit_predict\",\"predictions\":{},\"latent_variances\":{},\"observation_variances\":{},\"variance_scope\":\"plug_in_universal_kriging_independent_query_nugget_no_covariance_parameter_estimation_uncertainty\"}}",vector(&means),vector(&latent),vector(&observation)))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 1000
        || data[0].len() > 16
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data.iter().any(|r| {
            r.len() != data[0].len() || r.iter().any(|v| !v.is_finite() || v.abs() > 1e100)
        })
    {
        return Err("spatial estimation requires bounded finite rectangular inputs".into());
    }
    match op {
        "spatial_fit" => fit(&data.to_vec(), args),
        "spatial_fit_predict" => predict(&data.to_vec(), args),
        _ => Err("unknown spatial estimation operation".into()),
    }
}
