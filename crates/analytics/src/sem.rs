//! Bounded Gaussian simple-structure CFA,summary fit indices,factor scores,RAM.
//! Covariance inputs use divisor N (ML),not N-1;conversion is caller responsibility.
//! Latent variances are1,first loading per factor positive;free means count toward k.
//! No general SEM fitting,FIML,multigroup,robust correction,or causal guarantee.
type Mat = Vec<Vec<f64>>;
fn nums(a: &[f64]) -> String {
    format!(
        "[{}]",
        a.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn mat(a: &Mat) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| nums(r)).collect::<Vec<_>>().join(",")
    )
}
fn option(a: Option<f64>) -> String {
    a.filter(|x| x.is_finite())
        .map_or("null".into(), |v| v.to_string())
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("SEM numeric range exceeded".into())
    }
}
fn uint(x: f64) -> Result<usize, String> {
    if (0. ..=1e9).contains(&x) && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("parameter must be an exact nonnegative integer<=1e9".into())
    }
}
fn eye(n: usize) -> Mat {
    let mut a = vec![vec![0.; n]; n];
    for (i, r) in a.iter_mut().enumerate() {
        r[i] = 1.;
    }
    a
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(a: impl Iterator<Item = f64>) -> f64 {
    a.fold(0_f64, f64::hypot)
}
fn transpose(a: &Mat) -> Mat {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j]).collect())
        .collect()
}
fn mul(a: &Mat, b: &Mat) -> Result<Mat, String> {
    let bt = transpose(b);
    a.iter()
        .map(|r| bt.iter().map(|c| finite(dot(r, c))).collect())
        .collect()
}
fn chol(a: &Mat) -> Result<Mat, String> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        if a[i][i] <= 0. || !a[i][i].is_finite() {
            return Err("covariance diagonal must be positive".into());
        }
        for j in 0..=i {
            if (a[i][j] - a[j][i]).abs() > 1e-10 * a[i][i].sqrt() * a[j][j].sqrt() {
                return Err("covariance must be symmetric".into());
            }
            let x = finite(a[i][j] - dot(&l[i][..j], &l[j][..j]))?;
            l[i][j] = if i == j {
                if x <= 1e-12 * a[i][i] {
                    return Err("matrix not positive definite at relative pivot1e-12".into());
                }
                x.sqrt()
            } else {
                finite(x / l[j][j])?
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
fn solve(l: &Mat, b: &[f64]) -> Result<Vec<f64>, String> {
    let mut x = forward(l, b)?;
    for i in (0..x.len()).rev() {
        x[i] = finite((x[i] - (i + 1..x.len()).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i])?;
    }
    Ok(x)
}
fn spd_inverse(a: &Mat) -> Result<(Mat, f64), String> {
    let l = chol(a)?;
    let n = a.len();
    let mut out = vec![vec![0.; n]; n];
    for j in 0..n {
        let mut e = vec![0.; n];
        e[j] = 1.;
        let x = solve(&l, &e)?;
        for i in 0..n {
            out[i][j] = x[i];
        }
    }
    Ok((out, 2. * (0..n).map(|i| l[i][i].ln()).sum::<f64>()))
}
fn inverse(a: &Mat) -> Result<Mat, String> {
    let n = a.len();
    let mut w: Mat = a
        .iter()
        .enumerate()
        .map(|(i, r)| {
            r.iter()
                .copied()
                .chain((0..n).map(|j| f64::from(i == j)))
                .collect()
        })
        .collect();
    for c in 0..n {
        let pivot = (c..n)
            .max_by(|&i, &j| {
                let si = norm(w[i][..n].iter().copied());
                let sj = norm(w[j][..n].iter().copied());
                (w[i][c].abs() / si).total_cmp(&(w[j][c].abs() / sj))
            })
            .unwrap();
        let scale = norm(w[pivot][..n].iter().copied());
        if w[pivot][c].abs() <= scale * 1e-12 || !scale.is_finite() {
            return Err("linear system singular or ill conditioned".into());
        }
        w.swap(c, pivot);
        let d = w[c][c];
        for x in &mut w[c] {
            *x = finite(*x / d)?;
        }
        let row = w[c].clone();
        for (i, r) in w.iter_mut().enumerate() {
            if i != c {
                let f = r[c];
                for (j, x) in r.iter_mut().enumerate() {
                    *x = finite(*x - f * row[j])?;
                }
            }
        }
    }
    Ok(w.into_iter().map(|r| r[n..].to_vec()).collect())
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let w = data.first().map_or(0, Vec::len);
    if w == 0
        || w > 256
        || data.len() > 100000
        || data.len() > 1_000_000 / w
        || data
            .iter()
            .any(|r| r.len() != w || r.iter().any(|x| !x.is_finite()))
        || args.len() > 4096
        || args.iter().any(|x| !x.is_finite())
    {
        return Err(
            "SEM data and parameters must be finite rectangular nonempty and bounded".into(),
        );
    }
    match op {
        "cfa_model" => cfa(data, args),
        "cfa_scores" => scores(data, args),
        "sem_fit_indices" => indices(data, args),
        "sem_ram" => ram(data, args),
        _ => Err("unknown SEM operation".into()),
    }
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
    let zz = z - 1.;
    let x = C
        .iter()
        .enumerate()
        .fold(0.9999999999998099, |s, (i, c)| s + c / (zz + i as f64 + 1.));
    let t = zz + 7.5;
    0.5 * (2. * std::f64::consts::PI).ln() + (zz + 0.5) * t.ln() - t + x.ln()
}
fn gamma_pq(a: f64, x: f64) -> Result<(f64, f64), String> {
    if x == 0. {
        return Ok((0., 1.));
    }
    let pref = (a * x.ln() - x - log_gamma(a)).exp();
    if x < a + 1. {
        let mut sum = 1. / a;
        let mut term = sum;
        for i in 1..200000 {
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() < sum.abs() * 2e-15 {
                let p = (sum * pref).clamp(0., 1.);
                return Ok((p, 1. - p));
            }
        }
    } else {
        let tiny = 1e-300;
        let mut b = x + 1. - a;
        let mut c = 1. / tiny;
        let mut d = 1. / b;
        let mut h = d;
        for i in 1..200000 {
            let fi = i as f64;
            let an = -fi * (fi - a);
            b += 2.;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1. / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-15 {
                let q = (pref * h).clamp(0., 1.);
                return Ok((1. - q, q));
            }
        }
    }
    Err("gamma evaluation convergence limit".into())
}
fn nc_cdf(x: f64, df: f64, lambda: f64) -> Result<f64, String> {
    if x == 0. {
        return Ok(0.);
    }
    if lambda == 0. {
        return Ok(gamma_pq(df / 2., x / 2.)?.0);
    }
    if !(0. ..=1e7).contains(&lambda) {
        return Err("noncentrality exceeds supported1e7".into());
    }
    let mean = lambda / 2.;
    let k = mean.floor() as usize;
    let xx = x / 2.;
    let a = df / 2. + k as f64;
    let weight = (-mean + k as f64 * mean.ln() - log_gamma(k as f64 + 1.)).exp();
    let p = gamma_pq(a, xx)?.0;
    let term = (a * xx.ln() - xx - log_gamma(a + 1.)).exp();
    let mut total = weight * p;
    let mut mass = weight;
    let (mut wu, mut wd) = (weight, weight);
    let (mut pu, mut pd) = (p, p);
    let (mut tu, mut td) = (term, term);
    let (mut au, mut ad) = (a, a);
    for step in 1..=50000 {
        wu *= mean / (k + step) as f64;
        pu = (pu - tu).clamp(0., 1.);
        tu *= xx / (au + 1.);
        au += 1.;
        mass += wu;
        total += wu * pu;
        if step <= k {
            wd *= (k - step + 1) as f64 / mean;
            td *= ad / xx;
            ad -= 1.;
            pd = (pd + td).clamp(0., 1.);
            mass += wd;
            total += wd * pd;
        } else {
            wd = 0.;
        }
        if step > 20 && wu + wd < 1e-16 * mass {
            return Ok((total / mass).clamp(0., 1.));
        }
    }
    Err("noncentral chi square mixture limit50000".into())
}
fn nc_interval(x: f64, df: f64, n: f64) -> Result<[f64; 2], String> {
    let central = gamma_pq(df / 2., x / 2.)?.0;
    let mut result = [0.; 2];
    for (i, q) in [0.95, 0.05].into_iter().enumerate() {
        if central <= q {
            continue;
        }
        let mut hi = x.clamp(1., 1e7);
        while nc_cdf(x, df, hi)? > q {
            if hi >= 1e7 {
                return Err("noncentrality CI beyond supported1e7".into());
            }
            hi = (2. * hi).min(1e7);
        }
        let mut lo = 0.;
        for _ in 0..75 {
            let mid = (lo + hi) / 2.;
            if nc_cdf(x, df, mid)? > q {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        result[i] = ((lo + hi) / (2. * n * df)).sqrt();
    }
    Ok(result)
}
struct FitStats {
    body: String,
}
fn fit_stats(
    s: &Mat,
    v: &Mat,
    observed: &[f64],
    means: &[f64],
    n: usize,
    k: usize,
    inference: bool,
) -> Result<FitStats, String> {
    let p = s.len();
    let moments = p * (p + 3) / 2;
    if n <= p || k > moments {
        return Err("sample size must exceed p;free parameter count cannot exceed unique mean/covariance moments".into());
    }
    let (si, lds) = spd_inverse(s)?;
    let _ = si;
    let (vi, ldv) = spd_inverse(v)?;
    let delta: Vec<_> = observed.iter().zip(means).map(|(a, b)| a - b).collect();
    let quad = dot(
        &delta,
        &vi.iter().map(|r| dot(r, &delta)).collect::<Vec<_>>(),
    );
    let tr = (0..p)
        .map(|i| (0..p).map(|j| vi[i][j] * s[j][i]).sum::<f64>())
        .sum::<f64>();
    let f = finite(ldv - lds + tr + quad - p as f64)?;
    if f < -1e-8 {
        return Err("negative Gaussian discrepancy beyond rounding tolerance".into());
    }
    let chi = finite(n as f64 * f.max(0.))?;
    let df = moments - k;
    let baseline_discrepancy = (0..p).map(|i| s[i][i].ln()).sum::<f64>() - lds;
    let baseline = if baseline_discrepancy.abs() <= 64. * f64::EPSILON * p as f64 {
        0.
    } else {
        finite(n as f64 * baseline_discrepancy.max(0.))?
    };
    let df0 = p * (p - 1) / 2;
    let ll =
        finite(-0.5 * n as f64 * (p as f64 * (2. * std::f64::consts::PI).ln() + ldv + tr + quad))?;
    let cp = if inference && df > 0 {
        Some(gamma_pq(df as f64 / 2., chi / 2.)?.1)
    } else {
        None
    };
    let excess = (baseline - df0 as f64).max(0.);
    let cfi = if excess > 0. {
        Some((1. - (chi - df as f64).max(0.) / excess).clamp(0., 1.))
    } else {
        None
    };
    let rni = if baseline != df0 as f64 {
        Some(1. - (chi - df as f64) / (baseline - df0 as f64))
    } else {
        None
    };
    let tli = if df > 0 && df0 > 0 && (baseline / df0 as f64 - 1.).abs() > 1e-14 {
        Some(
            ((baseline / df0 as f64 - chi / df as f64) / (baseline / df0 as f64 - 1.))
                .clamp(0., 1.),
        )
    } else {
        None
    };
    let nfi = if baseline > 0. {
        Some((1. - chi / baseline).clamp(0., 1.))
    } else {
        None
    };
    let rmsea = if df > 0 {
        Some(((chi - df as f64).max(0.) / (n as f64 * df as f64)).sqrt())
    } else {
        None
    };
    let (ci, status) = if !inference {
        ("null".into(), "unconverged_boundary_or_unidentified")
    } else if df > 0 {
        match nc_interval(chi, df as f64, n as f64) {
            Ok(a) => (nums(&a), "computed"),
            Err(_) => ("null".into(), "noncentrality_or_iteration_limit"),
        }
    } else {
        ("null".into(), "undefined_zero_df")
    };
    let mut residual = Vec::new();
    let mut standard = Vec::new();
    for i in 0..p {
        for j in 0..=i {
            residual.push(s[i][j] - v[i][j]);
            standard.push((s[i][j] - v[i][j]) / s[i][i].sqrt() / s[j][j].sqrt());
        }
        residual.push(delta[i]);
        standard.push(delta[i] / s[i][i].sqrt());
    }
    let rmr = finite(norm(residual.into_iter()) / (moments as f64).sqrt())?;
    let srmr = finite(norm(standard.into_iter()) / (moments as f64).sqrt())?;
    let aic = finite(-2. * ll + 2. * k as f64)?;
    let aicc = if inference && n > k + 1 {
        Some(aic + 2. * k as f64 * (k + 1) as f64 / (n - k - 1) as f64)
    } else {
        None
    };
    Ok(FitStats{body:format!("\"n_free_parameters\":{k},\"degrees_of_freedom\":{df},\"chi_square\":{chi},\"chi_square_p_value\":{},\"baseline_chi_square\":{baseline},\"baseline_degrees_of_freedom\":{df0},\"cfi\":{},\"rni\":{},\"tli\":{},\"nfi\":{},\"rmsea\":{},\"rmsea_confidence_interval\":{ci},\"rmsea_ci_status\":\"{status}\",\"rmr\":{rmr},\"srmr\":{srmr},\"log_likelihood\":{ll},\"aic\":{},\"aicc\":{},\"bic\":{},\"bicu\":{},\"fit_index_convention\":\"ML divisorN;means included;SRMR lower covariance plus means over p(p+3)/2;RMSEA formula has no arbitrary upper clipping\"",option(cp),option(cfi),option(rni),option(tli),option(nfi),option(rmsea),option(if inference{Some(aic)}else{None}),option(aicc),option(if inference{Some(-2.*ll+k as f64*(n as f64).ln())}else{None}),option(if inference{Some(chi-df as f64*(n as f64).ln())}else{None}))})
}
fn indices(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.len();
    if p > 20 || data[0].len() != 2 * p + 2 || args.len() != 2 {
        return Err("fit indices require p<=20rows[Srow(p),Sigma row(p),samplemean,modelmean];args[N,kfree]".into());
    }
    let n = uint(args[0])?;
    let k = uint(args[1])?;
    let s = data.iter().map(|r| r[..p].to_vec()).collect();
    let v = data.iter().map(|r| r[p..2 * p].to_vec()).collect();
    let observed: Vec<_> = data.iter().map(|r| r[2 * p]).collect();
    let means: Vec<_> = data.iter().map(|r| r[2 * p + 1]).collect();
    let stats = fit_stats(&s, &v, &observed, &means, n, k, true)?;
    Ok(format!(
        "{{\"version\":1,\"method\":\"supplied_Gaussian_mean_covariance_fit\",{}}}",
        stats.body
    ))
}
struct CfaSpec {
    p: usize,
    q: usize,
    n: usize,
    s: Mat,
    assign: Vec<usize>,
    anchor: Vec<usize>,
    correlated: bool,
    calls: usize,
}
struct CfaModel {
    loading: Vec<f64>,
    phi: Mat,
    unique: Vec<f64>,
    cov: Mat,
}
impl CfaSpec {
    fn decode(&self, t: &[f64]) -> Result<CfaModel, String> {
        let p = self.p;
        let q = self.q;
        let mut loading = t[..p].to_vec();
        for &i in &self.anchor {
            if t[i].abs() > 25. {
                return Err("anchor log loading coordinate outside[-25,25]".into());
            }
            loading[i] = t[i].exp();
        }
        if loading.iter().any(|x| x.abs() > 100.) || t[p..2 * p].iter().any(|x| x.abs() > 25.) {
            return Err("loading/uniqueness coordinate outside supported range".into());
        }
        let unique: Vec<_> = t[p..2 * p].iter().map(|x| x.exp()).collect();
        let mut l = eye(q);
        if self.correlated {
            let mut k = 2 * p;
            for (i, j) in (0..q).flat_map(|i| (0..i).map(move |j| (i, j))) {
                if t[k].abs() > 1000. {
                    return Err("factor correlation coordinate exceeds1000".into());
                }
                l[i][j] = t[k];
                k += 1;
            }
        }
        let mut phi = mul(&l, &transpose(&l))?;
        let sd: Vec<_> = (0..q).map(|i| phi[i][i].sqrt()).collect();
        for i in 0..q {
            for j in 0..q {
                phi[i][j] /= sd[i] * sd[j];
            }
        }
        let cov = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        loading[i] * loading[j] * phi[self.assign[i]][self.assign[j]]
                            + if i == j { unique[i] } else { 0. }
                    })
                    .collect()
            })
            .collect();
        Ok(CfaModel {
            loading,
            phi,
            unique,
            cov,
        })
    }
    fn value(&mut self, t: &[f64]) -> Result<f64, String> {
        self.calls += 1;
        if self.calls > 100000 || self.calls * self.p * self.p * self.p > 1_000_000_000 {
            return Err("CFA likelihood evaluation budget exceeded".into());
        }
        let model = self.decode(t)?;
        let (vi, ld) = spd_inverse(&model.cov)?;
        let (_, lds) = spd_inverse(&self.s)?;
        let tr = (0..self.p).map(|i| dot(&vi[i], &self.s[i])).sum::<f64>();
        finite(0.5 * self.n as f64 * (ld - lds + tr - self.p as f64))
    }
    fn gradient(&mut self, t: &[f64]) -> Result<Vec<f64>, String> {
        self.calls += 1;
        if self.calls > 100000 || self.calls * self.p * self.p * self.p > 1_000_000_000 {
            return Err("CFA likelihood/gradient work budget exceeded".into());
        }
        let model = self.decode(t)?;
        let (vi, _) = spd_inverse(&model.cov)?;
        let vsvi = mul(&mul(&vi, &self.s)?, &vi)?;
        let p = self.p;
        let q = self.q;
        let e: Mat = (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| {
                        self.n as f64 * 0.25 * ((vi[i][j] - vsvi[i][j]) + (vi[j][i] - vsvi[j][i]))
                    })
                    .collect()
            })
            .collect();
        let mut gradient = vec![0.; t.len()];
        for i in 0..p {
            let loading_derivative = 2.
                * (0..p)
                    .map(|j| e[i][j] * model.loading[j] * model.phi[self.assign[i]][self.assign[j]])
                    .sum::<f64>();
            gradient[i] = loading_derivative
                * if self.anchor.contains(&i) {
                    model.loading[i]
                } else {
                    1.
                };
            gradient[p + i] = e[i][i] * model.unique[i];
        }
        if self.correlated {
            let mut l = eye(q);
            let mut k = 2 * p;
            for (a, b) in (0..q).flat_map(|a| (0..a).map(move |b| (a, b))) {
                l[a][b] = t[k];
                k += 1;
            }
            let var: Vec<_> = l.iter().map(|r| dot(r, r)).collect();
            let mut k = 2 * p;
            for (a, b) in (0..q).flat_map(|a| (0..a).map(move |b| (a, b))) {
                let mut g = 0.;
                for (i, row) in e.iter().enumerate() {
                    for (j, &entry) in row.iter().enumerate() {
                        let f = self.assign[i];
                        let h = self.assign[j];
                        let derivative = ((if f == a { l[h][b] } else { 0. })
                            + (if h == a { l[f][b] } else { 0. }))
                            / var[f].sqrt()
                            / var[h].sqrt()
                            - model.phi[f][h]
                                * ((if f == a { l[a][b] / var[a] } else { 0. })
                                    + (if h == a { l[a][b] / var[a] } else { 0. }));
                        g += entry * model.loading[i] * model.loading[j] * derivative;
                    }
                }
                gradient[k] = g;
                k += 1;
            }
        }
        gradient.into_iter().map(finite).collect()
    }
}
struct Opt {
    t: Vec<f64>,
    converged: bool,
    iterations: usize,
    gradient: f64,
}
fn optimize(spec: &mut CfaSpec, start: Vec<f64>, maxiter: usize, tol: f64) -> Result<Opt, String> {
    let mut t = start;
    let d = t.len();
    let mut h = eye(d);
    let mut f = spec.value(&t)?;
    let mut g = spec.gradient(&t)?;
    for iter in 0..maxiter {
        let gn = norm(g.iter().copied());
        if gn < tol.sqrt().clamp(1e-7, 1e-5) {
            return Ok(Opt {
                t,
                converged: true,
                iterations: iter,
                gradient: gn,
            });
        }
        let mut direction: Vec<_> = h.iter().map(|r| -dot(r, &g)).collect();
        if dot(&direction, &g) >= -1e-12 * gn * norm(direction.iter().copied()) {
            h = eye(d);
            direction = g.iter().map(|x| -x).collect();
        }
        let slope = dot(&direction, &g);
        let mut step = 1_f64.min(4. / norm(direction.iter().copied()));
        let mut accepted = None;
        for _ in 0..50 {
            let z: Vec<_> = t
                .iter()
                .zip(&direction)
                .map(|(a, b)| a + step * b)
                .collect();
            if let Ok(fz) = spec.value(&z) {
                if fz <= f + 1e-4 * step * slope {
                    accepted = Some((z, fz));
                    break;
                }
            }
            step *= 0.5;
        }
        let Some((z, fz)) = accepted else {
            return Ok(Opt {
                t,
                converged: false,
                iterations: iter,
                gradient: gn,
            });
        };
        let gz = spec.gradient(&z)?;
        let s: Vec<_> = z.iter().zip(&t).map(|(a, b)| a - b).collect();
        let y: Vec<_> = gz.iter().zip(&g).map(|(a, b)| a - b).collect();
        let sy = dot(&s, &y);
        if sy > 1e-12 * norm(s.iter().copied()) * norm(y.iter().copied()) {
            let hy: Vec<_> = h.iter().map(|r| dot(r, &y)).collect();
            let fac = (sy + dot(&y, &hy)) / (sy * sy);
            for i in 0..d {
                for j in 0..d {
                    h[i][j] += fac * s[i] * s[j] - (hy[i] * s[j] + s[i] * hy[j]) / sy;
                }
            }
        } else {
            h = eye(d);
        }
        t = z;
        g = gz;
        f = fz;
    }
    Ok(Opt {
        t,
        converged: false,
        iterations: maxiter,
        gradient: norm(g.iter().copied()),
    })
}
fn physical(spec: &CfaSpec, t: &[f64], scale: &[f64]) -> Result<Vec<f64>, String> {
    let m = spec.decode(t)?;
    let mut out: Vec<_> = m.loading.iter().zip(scale).map(|(a, b)| a * b).collect();
    out.extend(m.unique.iter().zip(scale).map(|(a, b)| a * b * b));
    if spec.correlated {
        for (i, j) in (0..spec.q).flat_map(|i| (0..i).map(move |j| (i, j))) {
            out.push(m.phi[i][j]);
        }
    }
    Ok(out)
}
fn cfa(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.len();
    if !(3..=20).contains(&p) || data[0].len() != p + 1 || args.len() != 5 + p {
        return Err("CFA requires p3..20summary rows[S_MLrow(p),mean],args[N,q,correlated,maxiter,tol,assignments(p)]".into());
    }
    let n = uint(args[0])?;
    let q = uint(args[1])?;
    let correlated = uint(args[2])?;
    let maxiter = uint(args[3])?;
    let tol = args[4];
    if !(1..=4).contains(&q)
        || correlated > 1
        || !(1..=1000).contains(&maxiter)
        || !(1e-12..=1e-3).contains(&tol)
        || n <= p
    {
        return Err(
            "invalid CFA sample count,factors,correlation flag,iteration/tolerance controls".into(),
        );
    }
    let assign = args[5..]
        .iter()
        .map(|x| uint(*x))
        .collect::<Result<Vec<_>, _>>()?;
    let mut anchors = Vec::new();
    for f in 0..q {
        let indices: Vec<_> = assign
            .iter()
            .enumerate()
            .filter(|(_, v)| **v == f)
            .map(|(i, _)| i)
            .collect();
        if indices.len() < 3 {
            return Err(
                "simple-structure CFA requires at least3indicators for every factor".into(),
            );
        }
        anchors.push(indices[0]);
    }
    if assign.iter().any(|i| *i >= q) {
        return Err("CFA assignments must be0..q-1".into());
    }
    let observed: Mat = data.iter().map(|r| r[..p].to_vec()).collect();
    chol(&observed)?;
    let scales: Vec<_> = (0..p).map(|i| observed[i][i].sqrt()).collect();
    let s: Mat = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| observed[i][j] / scales[i] / scales[j])
                .collect()
        })
        .collect();
    let means: Vec<_> = data.iter().map(|r| r[p]).collect();
    let d = 2 * p + if correlated == 1 { q * (q - 1) / 2 } else { 0 };
    let mut spec = CfaSpec {
        p,
        q,
        n,
        s,
        assign,
        anchor: anchors,
        correlated: correlated == 1,
        calls: 0,
    };
    let mut start = vec![0.; d];
    for i in 0..p {
        start[i] = if spec.anchor.contains(&i) {
            0.7_f64.ln()
        } else {
            0.7_f64.copysign(spec.s[i][spec.anchor[spec.assign[i]]])
        };
        start[p + i] = 0.51_f64.ln();
    }
    let opt = optimize(&mut spec, start, maxiter, tol)?;
    let m = spec.decode(&opt.t)?;
    let heywood = m.unique.iter().any(|x| *x <= 1e-8);
    let phi_chol = chol(&m.phi)?;
    let latent_boundary = spec.anchor.iter().any(|i| m.loading[*i] <= 1e-8)
        || (0..q).any(|i| phi_chol[i][i] * phi_chol[i][i] <= 1e-8)
        || opt.t[2 * p..].iter().any(|x| x.abs() > 999.);
    let boundary = heywood || latent_boundary;
    let mut h = vec![vec![0.; d]; d];
    // Observed information differentiates the analytic Gaussian score, with
    // Richardson extrapolation. No expected-information substitution is made.
    if !boundary {
        for j in 0..d {
            let step = 1e-4 * opt.t[j].abs().max(1.);
            let mut z = opt.t.clone();
            z[j] += step;
            let plus = spec.gradient(&z)?;
            z[j] -= 2. * step;
            let minus = spec.gradient(&z)?;
            z[j] += 1.5 * step;
            let half_plus = spec.gradient(&z)?;
            z[j] -= step;
            let half_minus = spec.gradient(&z)?;
            for i in 0..d {
                h[i][j] = (4. * (half_plus[i] - half_minus[i]) / step
                    - (plus[i] - minus[i]) / (2. * step))
                    / 3.;
            }
        }
        for (i, j) in (0..d).flat_map(|i| (0..i).map(move |j| (i, j))) {
            let average = (h[i][j] + h[j][i]) / 2.;
            h[i][j] = average;
            h[j][i] = average;
        }
    }
    let identified = !boundary && chol(&h).is_ok();
    let converged = opt.converged && identified && !boundary;
    let (se, pcov) = if converged {
        let (hi, _) = spd_inverse(&h)?;
        let mut jac = vec![vec![0.; d]; d];
        for j in 0..d {
            let step = 1e-5 * opt.t[j].abs().max(1.);
            let mut t = opt.t.clone();
            t[j] += step;
            let a = physical(&spec, &t, &scales)?;
            t[j] -= 2. * step;
            let b = physical(&spec, &t, &scales)?;
            for i in 0..d {
                jac[i][j] = (a[i] - b[i]) / (2. * step);
            }
        }
        let cov = mul(&mul(&jac, &hi)?, &transpose(&jac))?;
        let se: Vec<_> = (0..d).map(|i| cov[i][i].max(0.).sqrt()).collect();
        (nums(&se), mat(&cov))
    } else {
        ("null".into(), "null".into())
    };
    let loading: Mat = (0..p)
        .map(|i| {
            (0..q)
                .map(|f| {
                    if spec.assign[i] == f {
                        m.loading[i] * scales[i]
                    } else {
                        0.
                    }
                })
                .collect()
        })
        .collect();
    let unique: Vec<_> = m
        .unique
        .iter()
        .zip(&scales)
        .map(|(a, b)| a * b * b)
        .collect();
    let covariance: Mat = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| m.cov[i][j] * scales[i] * scales[j])
                .collect()
        })
        .collect();
    let stats = fit_stats(&observed, &covariance, &means, &means, n, d + p, converged)?;
    let mean_cov: Mat = covariance
        .iter()
        .map(|r| r.iter().map(|v| v / n as f64).collect())
        .collect();
    let mut saved = vec![p as f64, q as f64];
    saved.extend(&means);
    saved.extend(loading.iter().flatten());
    saved.extend(m.phi.iter().flatten());
    saved.extend(&unique);
    Ok(format!("{{\"version\":1,\"method\":\"simple_structure_Gaussian_CFA_ML\",\"loadings\":{},\"factor_covariance\":{},\"uniquenesses\":{},\"means\":{},\"implied_covariance\":{},\"standard_errors\":{se},\"parameter_covariance\":{pcov},\"mean_covariance\":{},\"score_parameters\":{},\"converged\":{converged},\"optimizer_converged\":{},\"identified\":{identified},\"heywood_boundary\":{heywood},\"latent_boundary\":{latent_boundary},\"iterations\":{},\"gradient_norm\":{},\"likelihood_evaluations\":{},\"inference\":\"observed information;SE null unless identified converged interior;latentvariance1 and positive first loading;means free\",{}}}",mat(&loading),mat(&m.phi),nums(&unique),nums(&means),mat(&covariance),mat(&mean_cov),nums(&saved),opt.converged,opt.iterations,opt.gradient,spec.calls,stats.body))
}
fn scores(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 3 {
        return Err("factor scores require[method,p,q,means,Lambda,Phi,theta]".into());
    }
    let method = uint(args[0])?;
    let p = uint(args[1])?;
    let q = uint(args[2])?;
    if method > 1
        || !(3..=20).contains(&p)
        || !(1..=4).contains(&q)
        || data[0].len() != p
        || args.len() != 3 + 2 * p + p * q + q * q
    {
        return Err("invalid saved factor-score parameter layout".into());
    }
    let means = &args[3..3 + p];
    let load: Mat = args[3 + p..3 + p + p * q]
        .chunks(q)
        .map(|r| r.to_vec())
        .collect();
    let phi: Mat = args[3 + p + p * q..3 + p + p * q + q * q]
        .chunks(q)
        .map(|r| r.to_vec())
        .collect();
    chol(&phi)?;
    let unique = &args[3 + p + p * q + q * q..];
    if unique.iter().any(|v| *v <= 0.) {
        return Err("score uniquenesses must be positive".into());
    }
    let lt = transpose(&load);
    let (weights, uncertainty) = if method == 0 {
        let mut v = mul(&mul(&load, &phi)?, &lt)?;
        for i in 0..p {
            v[i][i] += unique[i];
        }
        let (vi, _) = spd_inverse(&v)?;
        let weights = mul(&mul(&phi, &lt)?, &vi)?;
        let adjustment = mul(&mul(&weights, &load)?, &phi)?;
        let uc = (0..q)
            .map(|i| (0..q).map(|j| phi[i][j] - adjustment[i][j]).collect())
            .collect();
        (weights, uc)
    } else {
        let weighted: Mat = lt
            .iter()
            .map(|r| r.iter().zip(unique).map(|(v, u)| v / u).collect())
            .collect();
        let info = mul(&weighted, &load)?;
        let (inv, _) = spd_inverse(&info)?;
        (mul(&inv, &weighted)?, inv)
    };
    let rows: Mat = data
        .iter()
        .map(|r| {
            let x: Vec<_> = r.iter().zip(means).map(|(x, m)| x - m).collect();
            weights.iter().map(|w| finite(dot(w, &x))).collect()
        })
        .collect::<Result<_, _>>()?;
    Ok(format!("{{\"version\":1,\"method\":\"{}\",\"scores\":{},\"{}\":{},\"scope\":\"saved model factor scores only;does not propagate fitted parameter uncertainty\"}}",if method==0{"regression_factor_scores"}else{"Bartlett_factor_scores"},mat(&rows),if method==0{"conditional_factor_covariance"}else{"measurement_error_covariance"},mat(&uncertainty)))
}
fn ram(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.len();
    if p > 16 || data[0].len() != 2 * p + 1 || !args.is_empty() {
        return Err("RAM requires<=16rows[Brow(p),Omega row(p),intercept],no args".into());
    }
    let b: Mat = data.iter().map(|r| r[..p].to_vec()).collect();
    if (0..p).any(|i| b[i][i] != 0.) {
        return Err("RAM self paths must be zero".into());
    }
    let omega: Mat = data.iter().map(|r| r[p..2 * p].to_vec()).collect();
    chol(&omega)?;
    let a: Mat = (0..p)
        .map(|i| (0..p).map(|j| f64::from(i == j) - b[i][j]).collect())
        .collect();
    let ai = inverse(&a)?;
    let intercept: Vec<_> = data.iter().map(|r| r[2 * p]).collect();
    let means: Vec<_> = ai
        .iter()
        .map(|r| finite(dot(r, &intercept)))
        .collect::<Result<_, _>>()?;
    let cov = mul(&mul(&ai, &omega)?, &transpose(&ai))?;
    chol(&cov)?;
    let total: Mat = (0..p)
        .map(|i| (0..p).map(|j| ai[i][j] - f64::from(i == j)).collect())
        .collect();
    let indirect: Mat = (0..p)
        .map(|i| (0..p).map(|j| total[i][j] - b[i][j]).collect())
        .collect();
    // Strict discrete dynamic stability is separate from static invertibility.
    // B P B' - P + I=0 has SPD P exactly when B is strictly stable.
    let k: Mat = (0..p * p)
        .map(|r| {
            let (i, j) = (r / p, r % p);
            (0..p * p)
                .map(|c| {
                    let (a, d) = (c / p, c % p);
                    f64::from(r == c) - b[i][a] * b[j][d]
                })
                .collect()
        })
        .collect();
    let (stable, status) = match inverse(&k) {
        Ok(ki) => {
            let mut solution = vec![vec![0.; p]; p];
            for i in 0..p {
                for j in 0..p {
                    solution[i][j] = (0..p).map(|a| ki[i * p + j][a * p + a]).sum();
                }
            }
            let yes = chol(&solution).is_ok();
            (
                yes.to_string(),
                if yes {
                    "Lyapunov_SPD"
                } else {
                    "Lyapunov_not_SPD"
                },
            )
        }
        Err(_) => ("null".into(), "Lyapunov_singular_or_numerically_unresolved"),
    };
    Ok(format!("{{\"version\":1,\"method\":\"supplied_RAM_linear_system\",\"means\":{},\"covariance\":{},\"direct_effects\":{},\"total_effects\":{},\"indirect_effects\":{},\"invertible\":true,\"dynamically_stable\":{stable},\"stability_status\":\"{status}\",\"scope\":\"supplied parameters only;invertibility or positive Omega does not establish model identification,fit,or causality\"}}",nums(&means),mat(&cov),mat(&b),mat(&total),mat(&indirect)))
}
