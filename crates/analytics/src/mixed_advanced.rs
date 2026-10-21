//! Gaussian repeated covariance, known-covariance GLS, and correlated random slopes.
//!
//! Repeated rows: [fixed_features...,response,subject,time]. Nine covariance
//! structures use sorted global time levels; omitted visits retain their global
//! index. AR1 uses actual elapsed time (negative rho requires integral gaps);
//! Toeplitz requires a uniform global grid. No simultaneous G-side/R-side fit.
//! Random-slope rows: [other_fixed_features...,slope_covariate,response,subject].
//! The slope is included among both fixed and random effects. Independent groups,
//! a correlated 2x2 G, and independent residual error are fitted by profile ML/REML.
//! Conditional fixed-effect SEs do not use KR/Satterthwaite corrections.
//!
//! Bounds: <=256 fit observations, <=16 numeric fixed features, <=8 repeated
//! levels (builder <=64), 30000 likelihood evaluations and 1e9 nominal work units.
//! Variance-component hypothesis tests and covariance-parameter SEs are omitted.

use std::collections::BTreeMap;
type Mat = Vec<Vec<f64>>;
const MAX_ID: f64 = 9_007_199_254_740_991.;
const NAMES: [&str; 9] = [
    "heterogeneous_diagonal",
    "unstructured",
    "compound_symmetry",
    "heterogeneous_compound_symmetry",
    "ar1",
    "toeplitz",
    "heterogeneous_toeplitz",
    "antedependent",
    "equal_variance_antedependent",
];

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.first().map_or(0, Vec::len);
    if p == 0
        || p > 1024
        || data.len() > 100_000
        || data.len() > 1_000_000 / p
        || data
            .iter()
            .any(|r| r.len() != p || r.iter().any(|x| !x.is_finite()))
        || args.len() > 16384
        || args.iter().any(|x| !x.is_finite())
    {
        return Err(
            "data must be finite, rectangular, nonempty, and bounded; parameters finite".into(),
        );
    }
    match operation {
        "repeated_covariance" => builder(data, args),
        "repeated_model" => fit(data, args, false),
        "random_slope_model" => fit(data, args, true),
        "random_slope_predict" => predict(data, args),
        _ => Err("unknown advanced mixed-model operation".into()),
    }
}
fn nums(a: &[f64]) -> String {
    format!(
        "[{}]",
        a.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn mat(a: &Mat) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| nums(r)).collect::<Vec<_>>().join(",")
    )
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("calculation exceeds finite numeric range".into())
    }
}
fn integer(x: f64) -> Result<usize, String> {
    if (0. ..=MAX_ID).contains(&x) && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("parameter must be an exact nonnegative integer".into())
    }
}
fn id(x: f64) -> Result<i64, String> {
    if x.abs() <= MAX_ID && x.fract() == 0. {
        Ok(x as i64)
    } else {
        Err("subject ids must be exact signed integers".into())
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: impl Iterator<Item = f64>) -> f64 {
    a.fold(0_f64, f64::hypot)
}
fn eye(n: usize) -> Mat {
    let mut a = vec![vec![0.; n]; n];
    for (i, r) in a.iter_mut().enumerate() {
        r[i] = 1.;
    }
    a
}
fn mean(a: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut m = 0.;
    for (i, x) in a.enumerate() {
        m = finite(m + finite(x - m)? / (i + 1) as f64)?;
    }
    Ok(m)
}
// Relative to each original diagonal, so changing variable units does not make
// a well-conditioned correlation matrix appear singular.
fn chol(a: &Mat) -> Result<Mat, String> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        if a[i][i] <= 0. || !a[i][i].is_finite() {
            return Err("covariance must have positive finite variances".into());
        }
        for j in 0..=i {
            let s = finite(a[i][j] - dot(&l[i][..j], &l[j][..j]))?;
            l[i][j] = if i == j {
                if s <= a[i][i] * 1e-12 {
                    return Err(
                        "covariance is singular or not positive definite at relative pivot 1e-12"
                            .into(),
                    );
                }
                s.sqrt()
            } else {
                finite(s / l[j][j])?
            };
        }
    }
    Ok(l)
}
fn forward(l: &Mat, b: &[f64]) -> Result<Vec<f64>, String> {
    let mut x = vec![0.; b.len()];
    for i in 0..b.len() {
        x[i] = finite((b[i] - dot(&l[i][..i], &x[..i])) / l[i][i])?;
    }
    Ok(x)
}
fn solve_chol(l: &Mat, b: &[f64]) -> Result<Vec<f64>, String> {
    let mut x = forward(l, b)?;
    for i in (0..x.len()).rev() {
        x[i] = finite((x[i] - (i + 1..x.len()).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i])?;
    }
    Ok(x)
}
fn cov_count(code: usize, m: usize) -> usize {
    match code {
        0 => m,
        1 => m * (m + 1) / 2,
        2 | 4 => 2,
        3 => m + 1,
        5 | 8 => m,
        6 | 7 => 2 * m - 1,
        _ => 0,
    }
}
fn integral_gaps(t: &[f64]) -> bool {
    t.iter().all(|x| {
        let d = *x - t[0];
        d.abs() <= MAX_ID && (d - d.round()).abs() <= 32. * f64::EPSILON * d.abs().max(1.)
    })
}
fn grid(t: &[f64]) -> Result<(), String> {
    let h = t[1] - t[0];
    if t.iter().enumerate().any(|(i, x)| {
        let d = *x - t[0];
        (d - i as f64 * h).abs() > 64. * f64::EPSILON * d.abs().max(h.abs())
    }) {
        Err("Toeplitz requires a uniformly spaced global time grid".into())
    } else {
        Ok(())
    }
}
fn covariance(code: usize, t: &[f64], p: &[f64]) -> Result<Mat, String> {
    let m = t.len();
    if code > 8 || p.len() != cov_count(code, m) {
        return Err("incorrect covariance structure code or physical parameter count".into());
    }
    if t.windows(2)
        .any(|w| !((w[1] - w[0]).is_finite() && w[1] > w[0]))
    {
        return Err(
            "time levels must be distinct and strictly increasing with finite differences".into(),
        );
    }
    let mut a = vec![vec![0.; m]; m];
    if code == 1 {
        for (k, (i, j)) in (0..m)
            .flat_map(|i| (0..=i).map(move |j| (i, j)))
            .enumerate()
        {
            a[i][j] = p[k];
            a[j][i] = p[k];
        }
    } else {
        let hetero = matches!(code, 0 | 3 | 6 | 7);
        let ns = if hetero { m } else { 1 };
        if p[..ns].iter().any(|x| *x <= 0.) {
            return Err("variances must be strictly positive".into());
        }
        let sd: Vec<_> = (0..m)
            .map(|i| p[if hetero { i } else { 0 }].sqrt())
            .collect();
        let r = &p[ns..];
        if r.iter().any(|x| x.abs() >= 1.) {
            return Err("correlations must be strictly between -1 and 1".into());
        }
        if matches!(code, 2 | 3) && r[0] <= -1. / (m - 1) as f64 {
            return Err("compound symmetry rho must exceed -1/(levels-1)".into());
        }
        if matches!(code, 5 | 6) {
            grid(t)?;
        }
        if code == 4 && r[0] < 0. && !integral_gaps(t) {
            return Err("negative AR1 correlation requires integral actual time gaps".into());
        }
        for i in 0..m {
            for j in 0..=i {
                let c = if i == j {
                    1.
                } else {
                    match code {
                        0 => 0.,
                        2 | 3 => r[0],
                        4 => {
                            let d = t[i] - t[j];
                            if r[0] < 0. {
                                r[0].abs().powf(d.round())
                                    * if d.round().rem_euclid(2.) == 0. {
                                        1.
                                    } else {
                                        -1.
                                    }
                            } else {
                                r[0].powf(d)
                            }
                        }
                        5 | 6 => r[i - j - 1],
                        7 | 8 => r[j..i].iter().product(),
                        _ => unreachable!(),
                    }
                };
                a[i][j] = finite(sd[i] * c * sd[j])?;
                a[j][i] = a[i][j];
            }
        }
    }
    chol(&a)?;
    Ok(a)
}
fn builder(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || !(2..=64).contains(&data.len()) || args.is_empty() {
        return Err(
            "builder requires 2..64 scalar increasing time levels and structure parameters".into(),
        );
    }
    let code = integer(args[0])?;
    let t: Vec<_> = data.iter().map(|r| r[0]).collect();
    let a = covariance(code, &t, &args[1..])?;
    Ok(format!("{{\"version\":1,\"method\":\"{}\",\"time_levels\":{},\"covariance\":{},\"physical_parameters\":{},\"positive_definite\":true}}",NAMES[code],nums(&t),mat(&a),nums(&args[1..])))
}
struct Design {
    x: Mat,
    y: Vec<f64>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    yo: f64,
    ym: f64,
    ys: f64,
    groups: Vec<(i64, Vec<usize>)>,
    times: Vec<f64>,
    level: Vec<usize>,
    log_x_scale: f64,
    work: u64,
}
impl Design {
    fn new(data: &[Vec<f64>], random: bool) -> Result<Self, String> {
        let w = data[0].len();
        let end = if random { 2 } else { 3 };
        if w < end || (random && w < 3) || data.len() > 256 {
            return Err("fit requires <=256 observations and the documented feature/response/subject/time columns".into());
        }
        let p = w - end;
        let n = data.len();
        if p > 16 || n <= p + 2 {
            return Err(
                "fit requires <=16 features and at least two residual degrees of freedom".into(),
            );
        }
        let mut gm = BTreeMap::<i64, Vec<usize>>::new();
        for (i, r) in data.iter().enumerate() {
            gm.entry(id(r[p + 1])?).or_default().push(i);
        }
        if gm.len() < if random { 3 } else { 2 } || gm.len() == n {
            return Err("fit requires independent subjects and within-subject replication".into());
        }
        let groups: Vec<_> = gm.into_iter().collect();
        let mut times = Vec::new();
        let mut level = vec![0; n];
        if !random {
            times = data.iter().map(|r| r[p + 2]).collect();
            times.sort_by(f64::total_cmp);
            times.dedup();
            if !(2..=8).contains(&times.len()) {
                return Err("repeated fit requires 2..8 global time levels".into());
            }
            for (i, r) in data.iter().enumerate() {
                level[i] = times.binary_search_by(|x| x.total_cmp(&r[p + 2])).unwrap();
            }
            for (_, ix) in &groups {
                let mut seen = vec![false; times.len()];
                for &i in ix {
                    if seen[level[i]] {
                        return Err("duplicate subject/time observations are not allowed".into());
                    }
                    seen[level[i]] = true;
                }
            }
        }
        let origins = data[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| {
                mean(
                    data.iter()
                        .map(|r| finite(r[j] - origins[j]))
                        .collect::<Result<Vec<_>, _>>()?
                        .into_iter(),
                )
            })
            .collect::<Result<Vec<_>, String>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = finite(
                    norm(data.iter().map(|r| (r[j] - origins[j]) - offsets[j])) / (n as f64).sqrt(),
                )?;
                if s == 0. {
                    Err("fixed design has a constant feature".into())
                } else {
                    Ok(s)
                }
            })
            .collect::<Result<Vec<_>, String>>()?;
        let x = data
            .iter()
            .map(|r| {
                std::iter::once(1.)
                    .chain((0..p).map(|j| ((r[j] - origins[j]) - offsets[j]) / scales[j]))
                    .collect()
            })
            .collect();
        let yo = data[0][p];
        let ym = mean(
            data.iter()
                .map(|r| finite(r[p] - yo))
                .collect::<Result<Vec<_>, _>>()?
                .into_iter(),
        )?;
        let ys = finite(norm(data.iter().map(|r| (r[p] - yo) - ym)) / (n as f64).sqrt())?;
        if ys == 0. {
            return Err("constant response has unbounded Gaussian variance likelihood".into());
        }
        let y = data.iter().map(|r| ((r[p] - yo) - ym) / ys).collect();
        let k = p + 1;
        let work = groups
            .iter()
            .map(|(_, ix)| {
                let m = ix.len() as u64;
                m * m * m + m * m * (k + 1) as u64 + m * ((k + 1) * (k + 1)) as u64
            })
            .sum::<u64>()
            + (n * k * k) as u64;
        let log_x_scale = scales.iter().map(|v| v.ln()).sum();
        Ok(Self {
            x,
            y,
            origins,
            offsets,
            scales,
            yo,
            ym,
            ys,
            groups,
            times,
            level,
            log_x_scale,
            work,
        })
    }
}
struct Qr {
    beta: Vec<f64>,
    cov: Mat,
    quad: f64,
    logdet: f64,
}
fn qr(mut a: Mat, k: usize) -> Result<Qr, String> {
    let mut largest = 0_f64;
    for j in 0..k {
        let length = norm(a[j..].iter().map(|r| r[j]));
        largest = largest.max(length);
        if length <= largest * 1e-11 || !length.is_finite() {
            return Err("fixed design is rank deficient or numerically singular".into());
        }
        let alpha = -length.copysign(a[j][j]);
        let mut v: Vec<_> = a[j..].iter().map(|r| r[j]).collect();
        v[0] -= alpha;
        let vn = norm(v.iter().copied());
        for x in &mut v {
            *x /= vn;
        }
        for c in j..=k {
            let s = 2. * v.iter().zip(&a[j..]).map(|(x, r)| x * r[c]).sum::<f64>();
            for (r, x) in a[j..].iter_mut().zip(&v) {
                r[c] -= s * x;
            }
        }
        a[j][j] = alpha;
    }
    let mut beta = vec![0.; k];
    for i in (0..k).rev() {
        beta[i] = finite((a[i][k] - dot(&a[i][i + 1..k], &beta[i + 1..k])) / a[i][i])?;
    }
    let mut ri = vec![vec![0.; k]; k];
    for c in (0..k).rev() {
        for i in (0..=c).rev() {
            ri[i][c] = ((if i == c { 1. } else { 0. })
                - (i + 1..=c).map(|j| a[i][j] * ri[j][c]).sum::<f64>())
                / a[i][i];
        }
    }
    let mut cov = vec![vec![0.; k]; k];
    for i in 0..k {
        for j in 0..k {
            cov[i][j] = dot(&ri[i], &ri[j]);
        }
    }
    let quad = norm(a[k..].iter().map(|r| r[k])).powi(2);
    if quad <= 0. || !quad.is_finite() {
        return Err("zero or nonfinite residual variation".into());
    }
    let logdet = 2. * (0..k).map(|i| a[i][i].abs().ln()).sum::<f64>();
    Ok(Qr {
        beta,
        cov,
        quad,
        logdet,
    })
}
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
fn pacf_correlations(theta: &[f64]) -> Vec<f64> {
    let mut r = vec![1.];
    let mut a: Vec<f64> = Vec::new();
    let mut v = 1.;
    for &th in theta {
        let p = th.tanh();
        let lag = r.len();
        r.push((1..lag).map(|j| a[j - 1] * r[lag - j]).sum::<f64>() + p * v);
        let mut b: Vec<_> = (0..a.len())
            .map(|j| a[j] - p * a[a.len() - 1 - j])
            .collect();
        b.push(p);
        a = b;
        v *= 1. - p * p;
    }
    r[1..].to_vec()
}
fn shape(code: usize, t: &[f64], theta: &[f64]) -> Result<Mat, String> {
    if theta.iter().any(|v| v.abs() > 25.) {
        return Err("covariance coordinates exceed supported range".into());
    }
    let m = t.len();
    if code == 1 {
        let mut l = vec![vec![0.; m]; m];
        l[0][0] = 1.;
        let mut k = 0;
        for (i, row) in l.iter_mut().enumerate() {
            for (j, v) in row.iter_mut().enumerate().take(i + 1) {
                if i == 0 && j == 0 {
                    continue;
                }
                *v = if i == j { theta[k].exp() } else { theta[k] };
                k += 1;
            }
        }
        return Ok((0..m)
            .map(|i| (0..m).map(|j| dot(&l[i], &l[j])).collect())
            .collect());
    }
    let hetero = matches!(code, 0 | 3 | 6 | 7);
    let mut p = vec![1.];
    let mut k = 0;
    if hetero {
        for _ in 1..m {
            p.push((2. * theta[k]).exp());
            k += 1;
        }
    }
    match code {
        0 => {}
        2 | 3 => {
            let lo = -1. / (m - 1) as f64;
            p.push(lo + (1. - lo) * sigmoid(theta[k]));
        }
        4 => p.push(if integral_gaps(t) {
            theta[k].tanh()
        } else {
            sigmoid(theta[k])
        }),
        5 | 6 => p.extend(pacf_correlations(&theta[k..])),
        7 | 8 => p.extend(theta[k..].iter().map(|x| x.tanh())),
        _ => return Err("unknown covariance structure".into()),
    }
    covariance(code, t, &p)
}
struct Eval {
    objective: f64,
    qr: Qr,
    scale: f64,
    blocks: Vec<Mat>,
}
struct Engine<'a> {
    d: &'a Design,
    random: bool,
    code: usize,
    method: usize,
    known: Option<Mat>,
    calls: usize,
    work: u64,
}
impl Engine<'_> {
    fn eval(&mut self, theta: &[f64]) -> Result<Eval, String> {
        self.calls += 1;
        self.work = self.work.saturating_add(self.d.work);
        if self.calls > 30000 || self.work > 1_000_000_000 {
            return Err("covariance likelihood work budget exceeded".into());
        }
        let d = self.d;
        let k = d.x[0].len();
        let n = d.y.len();
        let full = if self.random {
            None
        } else {
            Some(if let Some(v) = &self.known {
                v.clone()
            } else {
                shape(self.code, &d.times, theta)?
            })
        };
        let mut a = vec![vec![0.; k + 1]; n];
        let mut logdet = 0.;
        let mut blocks = Vec::new();
        for (_, ix) in &d.groups {
            let m = ix.len();
            let mut v = vec![vec![0.; m]; m];
            if self.random {
                if theta.len() != 3 || theta.iter().any(|v| v.abs() > 1000.) {
                    return Err("random covariance coordinates outside bound".into());
                }
                let g00 = theta[0] * theta[0];
                let g01 = theta[0] * theta[1];
                let g11 = theta[1] * theta[1] + theta[2] * theta[2];
                for i in 0..m {
                    for j in 0..m {
                        let zi = d.x[ix[i]][k - 1];
                        let zj = d.x[ix[j]][k - 1];
                        v[i][j] =
                            g00 + g01 * (zi + zj) + g11 * zi * zj + if i == j { 1. } else { 0. };
                    }
                }
            } else {
                for i in 0..m {
                    for j in 0..m {
                        v[i][j] = full.as_ref().unwrap()[d.level[ix[i]]][d.level[ix[j]]];
                    }
                }
            }
            let l = chol(&v)?;
            logdet += 2. * (0..m).map(|i| l[i][i].ln()).sum::<f64>();
            for c in (0..=k).rev() {
                let b: Vec<_> = ix
                    .iter()
                    .map(|&i| if c == k { d.y[i] } else { d.x[i][c] })
                    .collect();
                let w = forward(&l, &b)?;
                for (j, &i) in ix.iter().enumerate() {
                    a[i][c] = w[j];
                }
            }
            blocks.push(l);
        }
        let qr = qr(a, k)?;
        let df = if self.method == 1 { n - k } else { n };
        let scale = if self.method == 2 {
            1.
        } else {
            qr.quad / df as f64
        };
        let objective = if self.method == 2 {
            0.5 * (n as f64 * (2. * std::f64::consts::PI).ln() + logdet + qr.quad)
        } else {
            0.5 * (df as f64 * ((2. * std::f64::consts::PI).ln() + 1. + scale.ln())
                + logdet
                + if self.method == 1 { qr.logdet } else { 0. })
        };
        finite(objective)?;
        Ok(Eval {
            objective,
            qr,
            scale,
            blocks,
        })
    }
    fn value(&mut self, x: &[f64]) -> f64 {
        self.eval(x).map_or(f64::INFINITY, |r| r.objective)
    }
    fn gradient(&mut self, x: &[f64]) -> Result<Vec<f64>, String> {
        let mut g = vec![0.; x.len()];
        for i in 0..x.len() {
            let h = 1e-5 * x[i].abs().max(1.);
            let mut p = x.to_vec();
            p[i] += h;
            let fp = self.eval(&p)?.objective;
            p[i] = x[i] - h;
            let fm = self.eval(&p)?.objective;
            g[i] = finite((fp - fm) / (2. * h))?;
        }
        Ok(g)
    }
}
struct Opt {
    theta: Vec<f64>,
    iterations: usize,
    converged: bool,
    gradient: f64,
}
fn optimize(e: &mut Engine<'_>, start: Vec<f64>, maxiter: usize, tol: f64) -> Result<Opt, String> {
    if start.is_empty() {
        return Ok(Opt {
            theta: start,
            iterations: 0,
            converged: true,
            gradient: 0.,
        });
    }
    let mut x = start;
    let d = x.len();
    let mut h = eye(d);
    let mut f = e.eval(&x)?.objective;
    let mut g = e.gradient(&x)?;
    let threshold = tol.sqrt().clamp(1e-7, 1e-5);
    for iter in 0..maxiter {
        let gn = norm(g.iter().copied());
        if gn < threshold {
            return Ok(Opt {
                theta: x,
                iterations: iter,
                converged: true,
                gradient: gn,
            });
        }
        let mut direction: Vec<_> = h.iter().map(|r| -dot(r, &g)).collect();
        if dot(&direction, &g) >= -1e-12 * gn * norm(direction.iter().copied()) {
            h = eye(d);
            direction = g.iter().map(|v| -v).collect();
        }
        let slope = dot(&g, &direction);
        let dn = norm(direction.iter().copied());
        let mut step = if dn > 4. { 4. / dn } else { 1. };
        let mut accepted = None;
        for _ in 0..40 {
            let z: Vec<_> = x
                .iter()
                .zip(&direction)
                .map(|(a, b)| a + step * b)
                .collect();
            let fz = e.value(&z);
            if fz.is_finite() && fz <= f + 1e-4 * step * slope {
                accepted = Some((z, fz));
                break;
            }
            step *= 0.5;
        }
        let Some((z, fz)) = accepted else {
            return Ok(Opt {
                theta: x,
                iterations: iter,
                converged: false,
                gradient: gn,
            });
        };
        let gz = e.gradient(&z)?;
        let s: Vec<_> = z.iter().zip(&x).map(|(a, b)| a - b).collect();
        let y: Vec<_> = gz.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        if sy > 1e-12 * norm(s.iter().copied()) * norm(y.iter().copied()) {
            let hy: Vec<_> = h.iter().map(|r| dot(r, &y)).collect();
            let factor = (sy + dot(&y, &hy)) / (sy * sy);
            for i in 0..d {
                for j in 0..d {
                    h[i][j] += factor * s[i] * s[j] - (hy[i] * s[j] + s[i] * hy[j]) / sy;
                }
            }
        } else {
            h = eye(d);
        }
        x = z;
        f = fz;
        g = gz;
    }
    Ok(Opt {
        theta: x,
        iterations: maxiter,
        converged: false,
        gradient: norm(g.iter().copied()),
    })
}
fn initial(code: usize, m: usize, t: &[f64]) -> Vec<f64> {
    let mut v = vec![0.; cov_count(code, m) - 1];
    if matches!(code, 2 | 3) {
        let i = v.len() - 1;
        v[i] = -((m - 1) as f64).ln();
    }
    if code == 4 {
        v[0] = if integral_gaps(t) {
            0.2_f64.atanh()
        } else {
            (0.2_f64 / 0.8).ln()
        };
    }
    v
}
fn physical(code: usize, v: &Mat, t: &[f64]) -> Vec<f64> {
    let m = v.len();
    if code == 0 {
        return (0..m).map(|i| v[i][i]).collect();
    }
    if code == 1 {
        return (0..m).flat_map(|i| (0..=i).map(move |j| v[i][j])).collect();
    }
    let hetero = matches!(code, 3 | 6 | 7);
    let mut p = if hetero {
        (0..m).map(|i| v[i][i]).collect()
    } else {
        vec![v[0][0]]
    };
    let cor = |i: usize, j: usize| v[i][j] / v[i][i].sqrt() / v[j][j].sqrt();
    match code {
        2 | 3 => p.push(cor(0, 1)),
        4 => {
            let c = cor(0, 1);
            let h = t[1] - t[0];
            p.push(if c < 0. {
                -(-c).powf(1. / h)
            } else {
                c.powf(1. / h)
            });
        }
        5 | 6 => p.extend((1..m).map(|j| cor(0, j))),
        7 | 8 => p.extend((1..m).map(|j| cor(j - 1, j))),
        _ => {}
    }
    p
}
fn transform_beta(d: &Design, b: &[f64], c: &Mat) -> Result<(Vec<f64>, Mat), String> {
    let k = b.len();
    let mut tr = vec![vec![0.; k]; k];
    tr[0][0] = d.ys;
    let mut out = vec![0.; k];
    out[0] = finite(d.yo + d.ym + d.ys * b[0])?;
    for j in 1..k {
        out[j] = finite(d.ys * b[j] / d.scales[j - 1])?;
        out[0] = finite(out[0] - out[j] * finite(d.origins[j - 1] + d.offsets[j - 1])?)?;
        tr[0][j] = -d.ys * (d.origins[j - 1] + d.offsets[j - 1]) / d.scales[j - 1];
        tr[j][j] = d.ys / d.scales[j - 1];
    }
    let mut cv = vec![vec![0.; k]; k];
    for i in 0..k {
        for j in 0..k {
            cv[i][j] = finite(
                (0..k)
                    .map(|a| (0..k).map(|b| tr[i][a] * c[a][b] * tr[j][b]).sum::<f64>())
                    .sum(),
            )?;
        }
    }
    Ok((out, cv))
}
fn fit(data: &[Vec<f64>], args: &[f64], random: bool) -> Result<String, String> {
    let shift = usize::from(!random);
    let code = if random {
        0
    } else {
        integer(
            *args
                .first()
                .ok_or("repeated_model requires structure code")?,
        )?
    };
    if code > 8 {
        return Err("unknown covariance structure code".into());
    }
    let method = integer(args.get(shift).copied().unwrap_or(1.))?;
    if method > if random { 1 } else { 2 } {
        return Err("method must be 0 ML, 1 REML, or 2 known repeated covariance".into());
    }
    let maxiter = integer(args.get(shift + 1).copied().unwrap_or(300.))?;
    let tol = args.get(shift + 2).copied().unwrap_or(1e-8);
    if !(1..=1000).contains(&maxiter) || !(1e-12..=1e-3).contains(&tol) {
        return Err("max iterations must be 1..1000 and tolerance 1e-12..1e-3".into());
    }
    let d = Design::new(data, random)?;
    let m = d.times.len();
    let k = d.x[0].len();
    let n = data.len();
    if !random && matches!(code, 5 | 6) {
        grid(&d.times)?;
    }
    let known = if method == 2 {
        if args.len() != shift + 3 + cov_count(code, m) {
            return Err("known covariance GLS requires complete physical covariance parameters after controls".into());
        }
        let a = covariance(code, &d.times, &args[shift + 3..])?;
        Some(
            a.iter()
                .map(|r| r.iter().map(|x| finite(x / d.ys / d.ys)).collect())
                .collect::<Result<Mat, _>>()?,
        )
    } else {
        if args.len() > shift + 3 {
            return Err("unexpected covariance fit parameters".into());
        }
        None
    };
    let dshape = if random { 3 } else { cov_count(code, m) - 1 };
    if method != 2 && n <= k + dshape + 1 {
        return Err("insufficient observations for fixed and covariance parameters".into());
    }
    let mut engine = Engine {
        d: &d,
        random,
        code,
        method,
        known,
        calls: 0,
        work: 0,
    };
    let opt = if method == 2 {
        Opt {
            theta: vec![],
            iterations: 0,
            converged: true,
            gradient: 0.,
        }
    } else {
        optimize(
            &mut engine,
            if random {
                vec![1., 0., 1.]
            } else {
                initial(code, m, &d.times)
            },
            maxiter,
            tol,
        )?
    };
    let ev = engine.eval(&opt.theta)?;
    let scale = ev.scale;
    let actual_scale = finite(scale * d.ys * d.ys)?;
    let bc: Mat = ev
        .qr
        .cov
        .iter()
        .map(|r| r.iter().map(|v| v * scale).collect())
        .collect();
    let (beta, bc) = transform_beta(&d, &ev.qr.beta, &bc)?;
    let se: Vec<_> = (0..k).map(|i| bc[i][i].max(0.).sqrt()).collect();
    let normalized_fitted: Vec<_> = d.x.iter().map(|r| dot(r, &ev.qr.beta)).collect();
    let fitted: Vec<_> = normalized_fitted
        .iter()
        .map(|v| finite(d.yo + d.ym + d.ys * v))
        .collect::<Result<_, _>>()?;
    let residuals: Vec<_> =
        d.y.iter()
            .zip(&normalized_fitted)
            .map(|(a, b)| (a - b) * d.ys)
            .collect();
    let ll = finite(
        -ev.objective
            - (if method == 1 { n - k } else { n }) as f64 * d.ys.ln()
            - if method == 1 { d.log_x_scale } else { 0. },
    )?;
    let np = if method == 2 { k } else { k + dshape + 1 };
    let mut boundary = false;
    let mut identified = method == 2;
    if method != 2 {
        // Central second differences of the independently profiled objective.
        // This is only an identification diagnostic, not calibrated inference.
        let dim = opt.theta.len();
        let mut h = vec![vec![0.; dim]; dim];
        let base = ev.objective;
        for i in 0..dim {
            let hi = 2e-4 * opt.theta[i].abs().max(1.);
            let mut p = opt.theta.clone();
            p[i] += hi;
            let fp = engine.value(&p);
            p[i] -= 2. * hi;
            let fm = engine.value(&p);
            h[i][i] = (fp + fm - 2. * base) / (hi * hi);
            for j in 0..i {
                let hj = 2e-4 * opt.theta[j].abs().max(1.);
                let mut z = opt.theta.clone();
                z[i] += hi;
                z[j] += hj;
                let pp = engine.value(&z);
                z[j] -= 2. * hj;
                let pm = engine.value(&z);
                z[i] -= 2. * hi;
                let mm = engine.value(&z);
                z[j] += 2. * hj;
                let mp = engine.value(&z);
                h[i][j] = (pp - pm + mm - mp) / (4. * hi * hj);
                h[j][i] = h[i][j];
            }
        }
        identified = dim == 0 || chol(&h).is_ok();
        boundary = opt.theta.iter().any(|x| x.abs() > 24.9);
    }
    let extra;
    let mut conditional = fitted.clone();
    if random {
        let th = &opt.theta;
        let gg = [
            [th[0] * th[0], th[0] * th[1]],
            [th[0] * th[1], th[1] * th[1] + th[2] * th[2]],
        ];
        let slope_scale = d.scales[k - 2];
        let slope_center = d.origins[k - 2] + d.offsets[k - 2];
        let g11 = actual_scale * gg[1][1] / slope_scale / slope_scale;
        let g01 = actual_scale * gg[0][1] / slope_scale - g11 * slope_center;
        let g00 = actual_scale * gg[0][0]
            - 2. * slope_center * actual_scale * gg[0][1] / slope_scale
            + g11 * slope_center * slope_center;
        let g = vec![
            vec![finite(g00)?, finite(g01)?],
            vec![finite(g01)?, finite(g11)?],
        ];
        let det = gg[0][0] * gg[1][1] - gg[0][1] * gg[0][1];
        boundary |= det <= 1e-8 * (gg[0][0] * gg[1][1]).max(1.) || gg[0][0].min(gg[1][1]) <= 1e-8;
        let mut blups = Vec::new();
        let mut ids = Vec::new();
        // Portable normalized training reference, then group/BLUP triples.
        let p = k - 1;
        let mut portable = vec![p as f64, d.yo, d.ym, d.ys, ev.qr.beta[0]];
        portable.extend(&d.origins);
        portable.extend(&d.offsets);
        portable.extend(&d.scales);
        portable.extend(&ev.qr.beta[1..]);
        portable.push(d.groups.len() as f64);
        for ((gid, ix), l) in d.groups.iter().zip(&ev.blocks) {
            let r: Vec<_> = ix.iter().map(|&i| d.y[i] - normalized_fitted[i]).collect();
            let vi = solve_chol(l, &r)?;
            let b0 = dot(&vec![1.; ix.len()], &vi);
            let b1 = ix
                .iter()
                .zip(&vi)
                .map(|(&i, v)| d.x[i][k - 1] * v)
                .sum::<f64>();
            let u0 = gg[0][0] * b0 + gg[0][1] * b1;
            let u1 = gg[1][0] * b0 + gg[1][1] * b1;
            ids.push(*gid as f64);
            blups.push(vec![
                d.ys * u0 - d.ys * u1 / slope_scale * slope_center,
                d.ys * u1 / slope_scale,
            ]);
            portable.extend([*gid as f64, u0, u1]);
            for &i in ix {
                conditional[i] =
                    finite(d.yo + d.ym + d.ys * (normalized_fitted[i] + u0 + u1 * d.x[i][k - 1]))?;
            }
        }
        let cr: Vec<_> = data
            .iter()
            .zip(&conditional)
            .map(|(r, p)| r[k - 1] - p)
            .collect();
        extra=format!(",\"residual_variance\":{actual_scale},\"random_effect_covariance\":{},\"subject_ids\":{},\"blups\":{},\"conditional_fitted_values\":{},\"conditional_residuals\":{},\"predict_parameters\":{}",mat(&g),nums(&ids),mat(&blups),nums(&conditional),nums(&cr),nums(&portable));
    } else {
        let mut v = if method == 2 {
            covariance(code, &d.times, &args[shift + 3..])?
        } else {
            shape(code, &d.times, &opt.theta)?
        };
        if method != 2 {
            for r in &mut v {
                for x in r {
                    *x = finite(*x * actual_scale)?;
                }
            }
        }
        if method != 2 {
            let l = chol(&v)?;
            let largest = (0..m).map(|i| v[i][i]).fold(0_f64, f64::max);
            boundary |=
                (0..m).any(|i| v[i][i] <= largest * 1e-12 || l[i][i] * l[i][i] / v[i][i] <= 1e-8);
        }
        let mut params = if method == 2 {
            args[shift + 3..].to_vec()
        } else {
            physical(code, &v, &d.times)
        };
        if code == 4 && method != 2 {
            params[1] = if integral_gaps(&d.times) {
                opt.theta[0].tanh()
            } else {
                sigmoid(opt.theta[0])
            };
        }
        extra = format!(
            ",\"time_levels\":{},\"covariance_by_level\":{},\"covariance_parameters\":{}",
            nums(&d.times),
            mat(&v),
            nums(&params)
        );
    }
    let converged = opt.converged && identified && !boundary;
    let method_label = if method == 0 {
        "ML"
    } else if method == 1 {
        "REML"
    } else {
        "known_covariance_GLS"
    };
    let estimate = if method == 2 { 1. } else { actual_scale };
    let aic = if method == 2 {
        String::from("null")
    } else {
        (2. * np as f64 - 2. * ll).to_string()
    };
    let bic = if method == 2 {
        String::from("null")
    } else {
        ((np as f64) * (n as f64).ln() - 2. * ll).to_string()
    };
    Ok(format!("{{\"version\":1,\"method\":\"{method_label}\",\"covariance_structure\":\"{}\",\"coefficients\":{},\"coefficient_covariance\":{},\"standard_errors\":{},\"fitted_values\":{},\"residuals\":{},\"log_likelihood\":{ll},\"estimated_scale\":{estimate},\"aic\":{aic},\"bic\":{bic},\"converged\":{converged},\"optimizer_converged\":{},\"covariance_parameters_identified\":{identified},\"variance_boundary\":{boundary},\"iterations\":{},\"gradient_norm\":{},\"likelihood_evaluations\":{},\"inference\":\"conditional GLS fixed-effect covariance; no small-sample df correction or covariance-parameter tests\"{extra}}}",if random{"random_intercept_and_slope"}else{NAMES[code]},nums(&beta),mat(&bc),nums(&se),nums(&fitted),nums(&residuals),opt.converged,opt.iterations,opt.gradient,engine.calls))
}
fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 6 {
        return Err("prediction requires exact returned predict_parameters".into());
    }
    let p = integer(args[0])?;
    if !(1..=16).contains(&p) || data[0].len() != p + 1 || args.len() < 6 + 4 * p || args[3] <= 0. {
        return Err("invalid prediction parameter layout".into());
    }
    let ng = integer(args[5 + 4 * p])?;
    if ng > 256 || args.len() != 6 + 4 * p + 3 * ng {
        return Err("invalid prediction group table".into());
    }
    let origins = &args[5..5 + p];
    let offsets = &args[5 + p..5 + 2 * p];
    let scales = &args[5 + 2 * p..5 + 3 * p];
    let beta = &args[5 + 3 * p..5 + 4 * p];
    if scales.iter().any(|v| *v <= 0.) {
        return Err("prediction scales must be positive".into());
    }
    let mut groups = BTreeMap::new();
    for r in args[6 + 4 * p..].chunks_exact(3) {
        if groups.insert(id(r[0])?, (r[1], r[2])).is_some() {
            return Err("duplicate prediction subject".into());
        }
    }
    let mut out = Vec::new();
    let mut marginal = Vec::new();
    let mut known = Vec::new();
    for r in data {
        let z: Vec<_> = (0..p)
            .map(|j| ((r[j] - origins[j]) - offsets[j]) / scales[j])
            .collect();
        let f = args[4] + dot(&z, beta);
        let g = groups.get(&id(r[p])?);
        known.push(g.is_some().to_string());
        marginal.push(finite(args[1] + args[2] + args[3] * f)?);
        out.push(finite(
            args[1] + args[2] + args[3] * (f + g.map_or(0., |(a, b)| a + b * z[p - 1])),
        )?);
    }
    Ok(format!(
        "{{\"version\":1,\"predictions\":{},\"marginal_predictions\":{},\"known_subjects\":[{}]}}",
        nums(&out),
        nums(&marginal),
        known.join(",")
    ))
}
