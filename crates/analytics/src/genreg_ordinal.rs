//! Bounded unpenalized cumulative-logit maximum likelihood.
//! P(Y <= j) = logistic(alpha_j + X beta), with sorted numeric outcome labels.
//! Observed-information uncertainty is a named native convention; JMP GenReg
//! documents a sandwich construction whose exact details are not established.

// Symmetric triangular solves and threshold blocks require cross-row indexing.
#![allow(clippy::needless_range_loop)]
use std::cell::Cell;
type Matrix = Vec<Vec<f64>>;
const MAX_WORK: usize = 200_000_000;
thread_local! {static WORK:Cell<usize>=const {Cell::new(0)};}
fn charge(n: usize) -> Result<(), String> {
    WORK.with(|w| {
        let n = w
            .get()
            .checked_add(n)
            .ok_or("ordinal numerical work overflow")?;
        w.set(n);
        if n > MAX_WORK {
            Err("ordinal cumulative numerical work exceeds 200000000 units".into())
        } else {
            Ok(())
        }
    })
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("ordinal result is not representable in finite binary64".into())
    }
}
fn integer(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.is_finite() && x.fract() == 0. && x >= lo as f64 && x <= hi as f64 {
        Ok(x as usize)
    } else {
        Err(format!("ordinal integer must be in [{lo},{hi}]"))
    }
}
fn vector(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| vector(r)).collect::<Vec<_>>().join(",")
    )
}
fn optional(x: Option<f64>) -> String {
    x.map_or_else(|| "null".into(), |x| x.to_string())
}
fn optional_matrix(x: &[Vec<Option<f64>>]) -> String {
    format!(
        "[{}]",
        x.iter()
            .map(|r| format!(
                "[{}]",
                r.iter().map(|&x| optional(x)).collect::<Vec<_>>().join(",")
            ))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn logistic(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
fn log_logistic(x: f64) -> f64 {
    if x >= 0. {
        -(-x).exp().ln_1p()
    } else {
        x - x.exp().ln_1p()
    }
}
fn log_interval(lo: f64, hi: f64, gap: f64) -> Result<f64, String> {
    if !gap.is_finite() || gap <= 0. {
        return Err("ordinal cutpoints must be strictly increasing".into());
    }
    // Preserve the physical cutpoint gap before adding a large linear predictor.
    let lg = if gap < std::f64::consts::LN_2 {
        (-(-gap).exp_m1()).ln()
    } else {
        (-(-gap).exp()).ln_1p()
    };
    finite(log_logistic(hi) + log_logistic(-lo) + lg)
}
fn scaled_cholesky(a: &Matrix) -> Result<(Matrix, f64), String> {
    let n = a.len();
    charge(n * n * n)?;
    let scale = (0..n).map(|i| a[i][i]).fold(0_f64, f64::max);
    if scale <= 0. || !scale.is_finite() {
        return Err("ordinal information is not positive definite".into());
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v = a[i][j] / scale - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if !v.is_finite() || v <= 1e-12 {
                    return Err(
                        "ordinal information is singular at relative pivot tolerance1e-12".into(),
                    );
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    Ok((l, scale))
}
fn inverse_spd(a: &Matrix) -> Result<Matrix, String> {
    let n = a.len();
    charge(3 * n * n * n)?;
    let (l, scale) = scaled_cholesky(a)?;
    let mut inv = vec![vec![0.; n]; n];
    for col in 0..n {
        let mut y = vec![0.; n];
        let mut x = vec![0.; n];
        for i in 0..n {
            y[i] =
                ((i == col) as u8 as f64 - (0..i).map(|k| l[i][k] * y[k]).sum::<f64>()) / l[i][i];
        }
        for i in (0..n).rev() {
            x[i] = (y[i] - (i + 1..n).map(|k| l[k][i] * x[k]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            inv[i][col] = finite(x[i] / scale)?;
        }
    }
    for i in 0..n {
        for j in 0..i {
            let v = 0.5 * inv[i][j] + 0.5 * inv[j][i];
            inv[i][j] = v;
            inv[j][i] = v;
        }
    }
    Ok(inv)
}
fn scaled_sd(g: &[f64], cov: &Matrix) -> Option<f64> {
    let gm = g.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    let cm = cov.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
    if gm == 0. || cm <= 0. || !gm.is_finite() || !cm.is_finite() {
        return None;
    }
    let mut q = 0.;
    for i in 0..g.len() {
        for j in 0..g.len() {
            q += (g[i] / gm) * (cov[i][j] / cm) * (g[j] / gm);
        }
    }
    if q <= 0. || !q.is_finite() {
        return None;
    }
    let se = (gm.ln() + 0.5 * cm.ln() + 0.5 * q.ln()).exp();
    (se > 0. && se.is_finite()).then_some(se)
}
fn transformed_covariance(jac: &Matrix, cov: &Matrix) -> Result<Option<Matrix>, String> {
    let n = jac.len();
    charge(4 * n.pow(4))?;
    let cm = cov.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
    let scales = jac
        .iter()
        .map(|r| r.iter().map(|v| v.abs()).fold(0_f64, f64::max))
        .collect::<Vec<_>>();
    if cm <= 0. || scales.iter().any(|s| *s == 0. || !s.is_finite()) {
        return Ok(None);
    }
    let mut c = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut q = 0.;
            for u in 0..n {
                for v in 0..n {
                    q += (jac[i][u] / scales[i]) * (cov[u][v] / cm) * (jac[j][v] / scales[j]);
                }
            }
            let value = if q == 0. {
                0.
            } else {
                q.signum() * (q.abs().ln() + scales[i].ln() + scales[j].ln() + cm.ln()).exp()
            };
            if !value.is_finite() || (i == j && value <= 0.) {
                return Ok(None);
            }
            c[i][j] = value;
            c[j][i] = value;
        }
    }
    Ok(Some(c))
}
fn normal_quantile(p: f64) -> f64 {
    // Acklam rational inverse, adequate to <5e-9 absolute on supported confidence levels.
    let a = [
        -39.69683028665376,
        220.9460984245205,
        -275.9285104469687,
        138.357_751_867_269,
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
        -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.)
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.)
    }
}
struct Problem {
    x: Matrix,
    y: Vec<usize>,
    w: Vec<f64>,
    total: f64,
    labels: Vec<f64>,
    counts: Vec<f64>,
    origin: Vec<f64>,
    offset: Vec<f64>,
    scale: Vec<f64>,
    p: usize,
    d: usize,
}
impl Problem {
    fn new(data: &[Vec<f64>]) -> Result<Self, String> {
        let p = data[0]
            .len()
            .checked_sub(2)
            .ok_or("fit rows need outcome and frequency")?;
        if p > 16 {
            return Err("ordinal supports at most16 predictors".into());
        }
        let mut rows = vec![];
        let mut total = 0.;
        let mut labels = vec![];
        for r in data {
            let w = integer(r[p + 1], 0, 1_000_000)? as f64;
            if w > 0. {
                total += w;
                rows.push((r, w));
                labels.push(r[p]);
            }
        }
        if total > 1e9 {
            return Err("ordinal total frequency exceeds1e9".into());
        }
        labels.sort_by(f64::total_cmp);
        labels.dedup_by(|a, b| *a == *b);
        let k = labels.len();
        if !(2..=16).contains(&k) {
            return Err("ordinal requires2..16 positive-frequency ordered outcome levels".into());
        }
        let d = k - 1 + p;
        if total <= d as f64 || rows.len() < p + 1 {
            return Err("ordinal positive-frequency design is too small".into());
        }
        charge(rows.len() * (p * p + 12 * p + 16))?;
        let mut origin = vec![0.; p];
        let mut offset = vec![0.; p];
        let mut scale = vec![0.; p];
        for j in 0..p {
            origin[j] = rows[0].0[j];
            offset[j] = rows
                .iter()
                .map(|(r, w)| (r[j] - origin[j]) * (w / total))
                .sum();
            let span = rows
                .iter()
                .map(|(r, _)| ((r[j] - origin[j]) - offset[j]).abs())
                .fold(0_f64, f64::max);
            if span <= 0. {
                return Err("ordinal design has a constant/aliased predictor".into());
            }
            scale[j] = span
                * rows
                    .iter()
                    .map(|(r, w)| (((r[j] - origin[j]) - offset[j]) / span).powi(2) * (w / total))
                    .sum::<f64>()
                    .sqrt();
            if !scale[j].is_finite() || scale[j] <= 0. {
                return Err("ordinal feature variation is unrepresentable".into());
            }
        }
        let x = rows
            .iter()
            .map(|(r, _)| {
                (0..p)
                    .map(|j| finite((r[j] - origin[j]) / scale[j] - offset[j] / scale[j]))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Matrix, _>>()?;
        // Weighted centered QR rejects an implicit common intercept and aliased slopes.
        let mut basis: Matrix = vec![];
        for j in 0..p {
            let mut v = x
                .iter()
                .zip(&rows)
                .map(|(r, (_, w))| r[j] * (w / total).sqrt())
                .collect::<Vec<_>>();
            for _ in 0..2 {
                for q in &basis {
                    let t = dot(&v, q);
                    for (a, b) in v.iter_mut().zip(q) {
                        *a -= t * b;
                    }
                }
            }
            let norm = dot(&v, &v).sqrt();
            if norm <= 1e-10 {
                return Err("ordinal design has a constant/aliased predictor".into());
            }
            for a in &mut v {
                *a /= norm;
            }
            basis.push(v);
        }
        let y = rows
            .iter()
            .map(|(r, _)| labels.iter().position(|v| *v == r[p]).unwrap())
            .collect::<Vec<_>>();
        let w = rows.iter().map(|(_, w)| *w).collect::<Vec<_>>();
        let mut counts = vec![0.; k];
        for (&y, &w) in y.iter().zip(&w) {
            counts[y] += w;
        }
        let problem = Self {
            x,
            y,
            w,
            total,
            labels,
            counts,
            origin,
            offset,
            scale,
            p,
            d,
        };
        charge(problem.x.len() * 8 * p * (k - 1))?;
        // A definite univariate complete/quasi-separation direction is cheap to certify.
        for feature in 0..p {
            for sign in [-1., 1.] {
                let mut direction = vec![0.; d];
                direction[k - 1 + feature] = sign;
                let mut possible = true;
                for cut in 0..k - 1 {
                    let lo = problem
                        .x
                        .iter()
                        .zip(&problem.y)
                        .filter(|(_, y)| **y <= cut)
                        .map(|(r, _)| -sign * r[feature])
                        .fold(f64::NEG_INFINITY, f64::max);
                    let hi = problem
                        .x
                        .iter()
                        .zip(&problem.y)
                        .filter(|(_, y)| **y > cut)
                        .map(|(r, _)| -sign * r[feature])
                        .fold(f64::INFINITY, f64::min);
                    if lo > hi {
                        possible = false;
                        break;
                    }
                    direction[cut] = 0.5 * lo + 0.5 * hi;
                }
                if possible && problem.separation(&direction)? {
                    return Err("ordinal complete/quasi separation: certified nonfinite-MLE recession direction".into());
                }
            }
        }
        Ok(problem)
    }
    fn separation(&self, v: &[f64]) -> Result<bool, String> {
        charge(self.x.len() * (self.p + 4))?;
        let m = self.labels.len() - 1;
        if v[..m].windows(2).any(|r| r[1] < r[0]) {
            return Ok(false);
        }
        let mut strict = false;
        for (x, &y) in self.x.iter().zip(&self.y) {
            let eta = dot(x, &v[m..]);
            if y < m {
                let s = v[y] + eta;
                if !s.is_finite() || s < 0. {
                    return Ok(false);
                }
                strict |= s > 0.;
            }
            if y > 0 {
                let s = v[y - 1] + eta;
                if !s.is_finite() || s > 0. {
                    return Ok(false);
                }
                strict |= s < 0.;
            }
        }
        Ok(strict)
    }
    fn evaluate(&self, t: &[f64]) -> Result<(f64, Vec<f64>, Matrix), String> {
        let m = self.labels.len() - 1;
        let d = self.d;
        charge(self.x.len() * (d * d + 12 * d + 64))?;
        if t.iter().any(|x| !x.is_finite() || x.abs() > 1e6)
            || t[..m].windows(2).any(|r| r[1] <= r[0])
        {
            return Err("ordinal parameter/order boundary".into());
        }
        let mut ll = 0.;
        let mut g = vec![0.; d];
        let mut info = vec![vec![0.; d]; d];
        for row in 0..self.x.len() {
            let x = &self.x[row];
            let y = self.y[row];
            let w = self.w[row] / self.total;
            let eta = finite(dot(x, &t[m..]))?;
            let (lp, ga, gb, fa, fb, c, sg) = if y == 0 {
                let b = finite(t[0] + eta)?;
                let f = logistic(b) * logistic(-b);
                (log_logistic(b), 0., logistic(-b), 0., f, 0., logistic(-b))
            } else if y == m {
                let a = finite(t[m - 1] + eta)?;
                let f = logistic(a) * logistic(-a);
                (log_logistic(-a), -logistic(a), 0., f, 0., 0., -logistic(a))
            } else {
                let a = finite(t[y - 1] + eta)?;
                let b = finite(t[y] + eta)?;
                let gap = t[y] - t[y - 1];
                let r = if gap > 700. { 0. } else { 1. / gap.exp_m1() };
                let c = finite(r * (1. + r))?;
                (
                    log_interval(a, b, gap)?,
                    -logistic(a) - r,
                    logistic(-b) + r,
                    logistic(a) * logistic(-a),
                    logistic(b) * logistic(-b),
                    c,
                    logistic(-a) - logistic(b),
                )
            };
            ll += w * lp;
            if y > 0 {
                g[y - 1] += w * ga;
                info[y - 1][y - 1] += w * (fa + c);
            }
            if y < m {
                g[y] += w * gb;
                info[y][y] += w * (fb + c);
            }
            if y > 0 && y < m {
                info[y - 1][y] -= w * c;
                info[y][y - 1] -= w * c;
            }
            for i in 0..self.p {
                g[m + i] += w * sg * x[i];
                if y > 0 {
                    let v = w * fa * x[i];
                    info[y - 1][m + i] += v;
                    info[m + i][y - 1] += v;
                }
                if y < m {
                    let v = w * fb * x[i];
                    info[y][m + i] += v;
                    info[m + i][y] += v;
                }
                for j in 0..self.p {
                    info[m + i][m + j] += w * (fa + fb) * x[i] * x[j];
                }
            }
        }
        finite(ll)?;
        if g.iter()
            .chain(info.iter().flatten())
            .any(|x| !x.is_finite())
        {
            return Err("ordinal score/information is unrepresentable".into());
        }
        Ok((ll, g, info))
    }
}
fn interval_json(values: &[f64], cov: &Option<Matrix>, z: f64) -> (String, String) {
    match cov {
        None => ("null".into(), "null".into()),
        Some(c) => {
            let se = (0..values.len())
                .map(|i| c[i][i].sqrt())
                .collect::<Vec<_>>();
            let ci = values
                .iter()
                .zip(&se)
                .map(|(&x, &s)| {
                    let (a, b) = (x - z * s, x + z * s);
                    if a.is_finite() && b.is_finite() {
                        vector(&[a, b])
                    } else {
                        "null".into()
                    }
                })
                .collect::<Vec<_>>();
            (vector(&se), format!("[{}]", ci.join(",")))
        }
    }
}
fn fitted(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 3 {
        return Err("ordinal fit args are[maxiter,tolerance,confidence]".into());
    }
    let maxiter = integer(args.first().copied().unwrap_or(200.), 0, 2000)?;
    let tol = args.get(1).copied().unwrap_or(1e-7);
    let confidence = args.get(2).copied().unwrap_or(0.95);
    if !(1e-10..=1e-3).contains(&tol) || !(0.5..0.999999).contains(&confidence) {
        return Err("invalid ordinal tolerance/confidence".into());
    }
    let p = Problem::new(data)?;
    let m = p.labels.len() - 1;
    let d = p.d;
    let mut initial = vec![0.; d];
    let mut count = 0.;
    for (j, t) in initial[..m].iter_mut().enumerate() {
        count += p.counts[j];
        let q = count / p.total;
        *t = q.ln() - (-q).ln_1p();
    }
    let mut theta = initial.clone();
    let mut evaluation = p.evaluate(&theta)?;
    let mut converged = false;
    let mut line_failed = false;
    let mut boundary = false;
    let mut singular = false;
    let mut iterations = 0;
    for iteration in 0..maxiter {
        iterations = iteration;
        let score = evaluation.1.iter().map(|x| x.abs()).fold(0_f64, f64::max);
        if score <= tol {
            converged = true;
            break;
        }
        let inverse = match inverse_spd(&evaluation.2) {
            Ok(v) => v,
            Err(_) => {
                singular = true;
                break;
            }
        };
        let step = inverse
            .iter()
            .map(|r| dot(r, &evaluation.1))
            .collect::<Vec<_>>();
        let slope = dot(&step, &evaluation.1);
        if !slope.is_finite() || slope <= 0. {
            line_failed = true;
            break;
        }
        let mut accepted = None;
        let mut fraction = 1.;
        for _ in 0..60 {
            let candidate = theta
                .iter()
                .zip(&step)
                .map(|(a, b)| a + fraction * b)
                .collect::<Vec<_>>();
            if let Ok(e) = p.evaluate(&candidate) {
                if e.0 >= evaluation.0 + 1e-4 * fraction * slope {
                    accepted = Some((candidate, e));
                    break;
                }
            }
            fraction *= 0.5;
        }
        match accepted {
            Some((t, e)) => {
                theta = t;
                evaluation = e;
                iterations = iteration + 1;
            }
            None => {
                line_failed = true;
                break;
            }
        }
    }
    let score = evaluation.1.iter().map(|x| x.abs()).fold(0_f64, f64::max);
    if maxiter > 0 && score <= tol {
        converged = true;
    }
    if theta.iter().any(|t| t.abs() >= 0.999e6) {
        converged = false;
        boundary = true;
    }
    let ray = theta
        .iter()
        .zip(&initial)
        .map(|(a, b)| a - b)
        .collect::<Vec<_>>();
    if p.separation(&ray)? {
        return Err(
            "ordinal complete/quasi separation: certified nonfinite-MLE recession direction".into(),
        );
    }
    let normalized_cov = if converged {
        inverse_spd(&evaluation.2).ok().map(|c| {
            c.iter()
                .map(|r| r.iter().map(|x| x / p.total).collect())
                .collect::<Matrix>()
        })
    } else {
        None
    };
    if converged && normalized_cov.is_none() {
        singular = true;
    }
    let coefficients = (0..p.p)
        .map(|j| finite(theta[m + j] / p.scale[j]))
        .collect::<Result<Vec<_>, _>>()?;
    let shift = finite(
        (0..p.p)
            .map(|j| coefficients[j] * (p.origin[j] + p.offset[j]))
            .sum(),
    )?;
    let cuts = theta[..m]
        .iter()
        .map(|x| finite(x - shift))
        .collect::<Result<Vec<_>, _>>()?;
    if cuts.windows(2).any(|r| r[1] <= r[0]) {
        return Err(
            "original-unit ordinal cutpoint order is unrepresentable; rescale predictor units"
                .into(),
        );
    }
    let mut reported = cuts.clone();
    reported.extend(&coefficients);
    let mut genreg = vec![cuts[0]];
    for j in 1..m {
        genreg.push(theta[j] - theta[j - 1]);
    }
    genreg.extend(&coefficients);
    let mut jac = vec![vec![0.; d]; d];
    for i in 0..m {
        jac[i][i] = 1.;
        for j in 0..p.p {
            jac[i][m + j] = -(p.origin[j] + p.offset[j]) / p.scale[j];
        }
    }
    for j in 0..p.p {
        jac[m + j][m + j] = 1. / p.scale[j];
    }
    let raw_cov = match &normalized_cov {
        Some(c) => transformed_covariance(&jac, c)?,
        None => None,
    };
    for i in 1..m {
        jac[i].fill(0.);
        jac[i][i] = 1.;
        jac[i][i - 1] = -1.;
    }
    let genreg_cov = match &normalized_cov {
        Some(c) => transformed_covariance(&jac, c)?,
        None => None,
    };
    let z = normal_quantile(0.5 * (1. + confidence));
    let (se, ci) = interval_json(&reported, &raw_cov, z);
    let (gse, gci) = interval_json(&genreg, &genreg_cov, z);
    let mut saved = vec![
        1.,
        p.p as f64,
        p.labels.len() as f64,
        confidence,
        normalized_cov.is_some() as u8 as f64,
    ];
    saved.extend(&p.origin);
    saved.extend(&p.offset);
    saved.extend(&p.scale);
    saved.extend(&p.labels);
    saved.extend(&theta);
    saved.extend(
        normalized_cov
            .clone()
            .unwrap_or_else(|| vec![vec![0.; d]; d])
            .into_iter()
            .flatten(),
    );
    let ll = finite(evaluation.0 * p.total)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"genreg_ordinal_fit\",\"n_samples\":{},\"effective_sample_size\":{},\"n_features\":{},\"ordered_labels\":{},\"level_counts\":{},\"cutpoints\":{},\"coefficients\":{},\"genreg_intercepts\":{},\"log_likelihood\":{},\"aic\":{},\"bic\":{},\"covariance\":{},\"standard_errors\":{},\"confidence_intervals\":{},\"genreg_covariance\":{},\"genreg_standard_errors\":{},\"genreg_confidence_intervals\":{},\"covariance_method\":\"observed_information\",\"jmp_sandwich_equivalence_verified\":false,\"converged\":{},\"information_identified\":{},\"inference_valid\":{},\"finite_mle_existence_verified\":false,\"global_optimum_verified\":false,\"boundary_reached\":{},\"information_singular\":{},\"line_search_failed\":{},\"iterations\":{},\"gradient_max\":{},\"confidence\":{},\"work_units\":{},\"predict_parameters\":{}}}",p.x.len(),p.total,p.p,vector(&p.labels),vector(&p.counts),vector(&cuts),vector(&coefficients),vector(&genreg[..m]),ll,finite(2.*d as f64-2.*ll)?,finite((d as f64)*p.total.ln()-2.*ll)?,raw_cov.as_ref().map_or_else(||"null".into(),|c|matrix(c)),se,ci,genreg_cov.as_ref().map_or_else(||"null".into(),|c|matrix(c)),gse,gci,converged,normalized_cov.is_some(),raw_cov.is_some(),boundary,singular,line_failed,iterations,score,confidence,WORK.with(Cell::get),vector(&saved)))
}
fn predicted(data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() < 5 || a[0] != 1. {
        return Err("invalid ordinal model version/header".into());
    }
    let p = integer(a[1], 0, 16)?;
    let k = integer(a[2], 2, 16)?;
    let m = k - 1;
    let d = m + p;
    let confidence = a[3];
    let has_cov = integer(a[4], 0, 1)? == 1;
    if !(0.5..0.999999).contains(&confidence)
        || a.len() != 5 + 3 * p + k + d + d * d
        || data[0].len() != p + 1
    {
        return Err("ordinal saved model dimensions/confidence mismatch".into());
    }
    if data.len() * k > 200000 {
        return Err("ordinal prediction exceeds200000 output probability cells".into());
    }
    let origin = &a[5..5 + p];
    let offset = &a[5 + p..5 + 2 * p];
    let scale = &a[5 + 2 * p..5 + 3 * p];
    if scale.iter().any(|x| *x <= 0.) {
        return Err("ordinal saved feature scales must be positive".into());
    }
    let labels = &a[5 + 3 * p..5 + 3 * p + k];
    let theta = &a[5 + 3 * p + k..5 + 3 * p + k + d];
    let cuts = &theta[..m];
    if labels.windows(2).any(|r| r[1] <= r[0]) || cuts.windows(2).any(|r| r[1] <= r[0]) {
        return Err("ordinal saved labels/cutpoints must be strictly ordered".into());
    }
    let cov = a[5 + 3 * p + k + d..]
        .chunks_exact(d)
        .map(|r| r.to_vec())
        .collect::<Matrix>();
    if !has_cov && cov.iter().flatten().any(|x| *x != 0.) {
        return Err("ordinal disabled covariance must be all zero".into());
    }
    if has_cov {
        for i in 0..d {
            for j in 0..i {
                let s = cov[i][j].abs().max(cov[j][i].abs());
                if s > 0. && (cov[i][j] / s - cov[j][i] / s).abs() > 1e-10 {
                    return Err("ordinal saved covariance must be symmetric".into());
                }
            }
        }
        scaled_cholesky(&cov)?;
    }
    let z = normal_quantile(0.5 * (1. + confidence));
    let mut probs = vec![];
    let mut cdfs = vec![];
    let mut modes = vec![];
    let mut quantiles = vec![];
    let mut eta_values = vec![];
    let mut errors = vec![];
    let mut cis = vec![];
    for row in data {
        charge(32 + p * 8 + k * (d * d + 8 * d + 48))?;
        let q = row[p];
        if !(0. ..1.).contains(&q) || q == 0. {
            return Err("ordinal quantile probability must be strictly(0,1)".into());
        }
        let x = (0..p)
            .map(|j| finite((row[j] - origin[j]) / scale[j] - offset[j] / scale[j]))
            .collect::<Result<Vec<_>, _>>()?;
        let eta = finite(dot(&x, &theta[m..]))?;
        eta_values.push(finite(
            (0..p).map(|j| (theta[m + j] / scale[j]) * row[j]).sum(),
        )?);
        let locations = cuts
            .iter()
            .map(|c| finite(c + eta))
            .collect::<Result<Vec<_>, _>>()?;
        let mut ps = vec![];
        let mut ss = vec![];
        let mut intervals = vec![];
        for j in 0..k {
            let probability = if j == 0 {
                logistic(locations[0])
            } else if j == m {
                logistic(-locations[m - 1])
            } else {
                log_interval(locations[j - 1], locations[j], cuts[j] - cuts[j - 1])?.exp()
            };
            let mut g = vec![0.; d];
            let fa = if j > 0 {
                logistic(locations[j - 1]) * logistic(-locations[j - 1])
            } else {
                0.
            };
            let fb = if j < m {
                logistic(locations[j]) * logistic(-locations[j])
            } else {
                0.
            };
            if j > 0 {
                g[j - 1] = -fa;
            }
            if j < m {
                g[j] = fb;
            }
            let sg = if j == 0 {
                fb
            } else if j == m {
                -fa
            } else {
                probability * (logistic(-locations[j - 1]) - logistic(locations[j]))
            };
            for i in 0..p {
                g[m + i] = sg * x[i];
            }
            let se = if has_cov { scaled_sd(&g, &cov) } else { None };
            let ci = se.and_then(|s| {
                let (a, b) = (probability - z * s, probability + z * s);
                (a.is_finite() && b.is_finite()).then_some(vec![a, b])
            });
            ps.push(probability);
            ss.push(se);
            intervals.push(ci.map_or_else(|| "null".into(), |v| vector(&v)));
        }
        let sum = ps.iter().sum::<f64>();
        if !sum.is_finite() || sum <= 0. {
            return Err("ordinal probability normalization failed".into());
        }
        if (sum - 1.).abs() > 1e-12 {
            return Err("ordinal probabilities failed sum-to-one consistency".into());
        }
        let mut c = locations.iter().map(|&x| logistic(x)).collect::<Vec<_>>();
        c.push(1.);
        let mode = ps
            .iter()
            .enumerate()
            .fold(0, |best, (i, v)| if *v > ps[best] { i } else { best });
        let quantile = c.iter().position(|v| *v >= q).unwrap_or(k - 1);
        modes.push(labels[mode]);
        quantiles.push(labels[quantile]);
        probs.push(ps);
        cdfs.push(c);
        errors.push(ss);
        cis.push(format!("[{}]", intervals.join(",")));
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"genreg_ordinal_predict\",\"ordered_labels\":{},\"probabilities\":{},\"cumulative_probabilities\":{},\"most_likely\":{},\"quantiles\":{},\"linear_predictor\":{},\"probability_standard_errors\":{},\"probability_confidence_intervals\":[{}],\"confidence\":{},\"covariance_method\":\"observed_information_delta_wald\",\"jmp_sandwich_equivalence_verified\":false,\"work_units\":{}}}",vector(labels),matrix(&probs),matrix(&cdfs),vector(&modes),vector(&quantiles),vector(&eta_values),optional_matrix(&errors),cis.join(","),confidence,WORK.with(Cell::get)))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    WORK.with(|w| w.set(0));
    let columns = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > 100000
        || columns == 0
        || columns > 1024
        || data.len().saturating_mul(columns) > 1000000
        || data
            .iter()
            .any(|r| r.len() != columns || r.iter().any(|x| !x.is_finite() || x.abs() > 1e150))
        || args.len() > 16384
        || args.iter().any(|x| !x.is_finite())
    {
        return Err("ordinal input must be bounded finite rectangular numeric data".into());
    }
    charge(8 * data.len() * columns + 8 * args.len())?;
    let result = match op {
        "genreg_ordinal_fit" => fitted(data, args),
        "genreg_ordinal_predict" => predicted(data, args),
        _ => Err("unknown ordinal operation".into()),
    };
    // A numeric singularity may withhold covariance, but must never hide a work
    // budget error caught inside a line search or an optional inversion.
    charge(0)?;
    result
}
