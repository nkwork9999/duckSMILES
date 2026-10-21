//! Bounded frequentist degradation paths, stability, independent destructive
//! degradation and cumulative-exposure step-stress reliability. These contracts
//! do not implement JMP's Bayesian repeated-measures degradation platform.
const MAX_WORK: usize = 200_000_000;
const LN_2PI: f64 = 1.8378770664093453;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let w = data.first().map_or(0, Vec::len);
    if w == 0
        || data.len() > 10000
        || w > 128
        || data.len() > 1_000_000 / w
        || args.len() > 16384
        || args.iter().any(|x| !x.is_finite())
        || data
            .iter()
            .any(|r| r.len() != w || r.iter().any(|x| !x.is_finite()))
    {
        return Err(
            "degradation requires finite nonempty rectangular bounded data/arguments".into(),
        );
    }
    match operation {
        "degradation_path" => path_fit(data, args),
        "degradation_predict" => path_predict(data, args),
        "stability_analysis" => stability(data, args),
        "destructive_degradation" => destructive_fit(data, args),
        "destructive_predict" => destructive_predict(data, args),
        "cumulative_damage_fit" => damage_fit(data, args),
        "cumulative_damage_predict" => damage_predict(data, args),
        _ => Err("unknown degradation operation".into()),
    }
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numeric range"))
    }
}
fn int(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.fract() == 0. && x >= lo as f64 && x <= hi as f64 {
        Ok(x as usize)
    } else {
        Err(format!("integer argument must be in [{lo},{hi}]"))
    }
}
fn id(x: f64) -> Result<(), String> {
    if x.fract() == 0. && x.abs() <= 9_007_199_254_740_991. {
        Ok(())
    } else {
        Err("unit/pattern identifier must be an exact integer within 2^53-1".into())
    }
}
fn num(x: f64) -> String {
    if x.is_finite() {
        x.to_string()
    } else {
        "null".into()
    }
}
fn arr(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(|x| num(*x)).collect::<Vec<_>>().join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn transform(x: f64, kind: usize, time: bool) -> Result<f64, String> {
    if time && x < 0. {
        return Err("degradation time must be nonnegative".into());
    }
    finite(
        match kind {
            0 => x,
            1 => {
                if x <= 0. {
                    return Err("log transform requires positive input".into());
                }
                x.ln()
            }
            2 => x.sqrt(),
            _ => return Err("invalid transform".into()),
        },
        "transformed value",
    )
}
fn inverse(x: f64, kind: usize) -> f64 {
    match kind {
        0 => x,
        1 => x.exp(),
        2 => x * x,
        _ => f64::NAN,
    }
}
struct Ols {
    beta: Vec<f64>,
    inverse: Vec<Vec<f64>>,
    sse: f64,
    df: usize,
}
#[allow(clippy::needless_range_loop)]
fn ols(x: &[Vec<f64>], y: &[f64]) -> Result<Ols, String> {
    let n = x.len();
    let d = x[0].len();
    if d == 0 || n <= d || n.saturating_mul(d).saturating_mul(d) > MAX_WORK {
        return Err("unidentified or excessive linear path design".into());
    }
    let scales = (0..d)
        .map(|j| x.iter().map(|r| r[j].abs()).fold(0_f64, f64::max))
        .collect::<Vec<_>>();
    if scales.contains(&0.) {
        return Err("linear path design is rank deficient".into());
    }
    let mut q: Vec<Vec<f64>> = Vec::new();
    let mut r = vec![vec![0.; d]; d];
    for j in 0..d {
        let mut v = x.iter().map(|row| row[j] / scales[j]).collect::<Vec<_>>();
        for _ in 0..2 {
            for k in 0..j {
                let a = dot(&q[k], &v);
                r[k][j] += a;
                for i in 0..n {
                    v[i] -= a * q[k][i];
                }
            }
        }
        let length = v.iter().fold(0_f64, |a, b| a.hypot(*b));
        if length <= 1e-11 * (n as f64).sqrt() {
            return Err("linear path design is rank deficient".into());
        }
        r[j][j] = length;
        q.push(v.iter().map(|v| v / length).collect());
    }
    let mut beta = q.iter().map(|col| dot(col, y)).collect::<Vec<_>>();
    for j in (0..d).rev() {
        for k in j + 1..d {
            beta[j] -= r[j][k] * beta[k];
        }
        beta[j] /= r[j][j];
    }
    for j in 0..d {
        beta[j] = finite(beta[j] / scales[j], "linear coefficient")?;
    }
    let sse = finite(
        x.iter()
            .zip(y)
            .map(|(x, y)| (y - dot(x, &beta)).powi(2))
            .sum(),
        "residual sum of squares",
    )?;
    let mut ri = vec![vec![0.; d]; d];
    for col in 0..d {
        for j in (0..d).rev() {
            let rhs = if j == col { 1. } else { 0. };
            ri[j][col] = (rhs - (j + 1..d).map(|k| r[j][k] * ri[k][col]).sum::<f64>()) / r[j][j];
        }
    }
    let inverse = (0..d)
        .map(|j| {
            (0..d)
                .map(|k| {
                    finite(
                        dot(&ri[j], &ri[k]) / scales[j] / scales[k],
                        "linear information inverse",
                    )
                })
                .collect()
        })
        .collect::<Result<Vec<Vec<f64>>, String>>()?;
    Ok(Ols {
        beta,
        inverse,
        sse,
        df: n - d,
    })
}
#[derive(Clone)]
struct Path {
    id: f64,
    a: f64,
    b: f64,
    caa: f64,
    cab: f64,
    cbb: f64,
    mse: f64,
    df: f64,
}
struct Paths {
    paths: Vec<Path>,
    sse: f64,
    df: usize,
    response: usize,
    time: usize,
    rank: usize,
}
fn fit_paths(
    data: &[Vec<f64>],
    im: usize,
    sm: usize,
    response: usize,
    time: usize,
) -> Result<Paths, String> {
    if data[0].len() != 3 {
        return Err("path rows must be [time,response,unit_id]".into());
    }
    let mut ids = Vec::new();
    let mut indices = Vec::new();
    let mut u = Vec::new();
    let mut y = Vec::new();
    for row in data {
        id(row[2])?;
        let g = if let Some(i) = ids.iter().position(|x| *x == row[2]) {
            i
        } else {
            ids.push(row[2]);
            ids.len() - 1
        };
        indices.push(g);
        u.push(transform(row[0], time, true)?);
        y.push(transform(row[1], response, false)?);
    }
    let groups = ids.len();
    if groups > 128 {
        return Err("at most 128 path units are supported".into());
    }
    let na = if im == 0 {
        groups
    } else if im == 1 {
        1
    } else {
        0
    };
    let nb = if sm == 0 { groups } else { 1 };
    let d = na + nb;
    let no_center = im == 2 || im == 1 && sm == 0;
    let origin = if no_center { 0. } else { u[0] };
    let offset = if no_center {
        0.
    } else {
        finite(
            u.iter().map(|v| (v - origin) / u.len() as f64).sum(),
            "time center",
        )?
    };
    let center = finite(origin + offset, "time center")?;
    let scale = u.iter().map(|v| (v - center).abs()).fold(0_f64, f64::max);
    if scale == 0. {
        return Err("time variation is required".into());
    }
    let mut x = vec![vec![0.; d]; data.len()];
    for i in 0..data.len() {
        if na > 0 {
            x[i][if im == 0 { indices[i] } else { 0 }] = 1.;
        }
        x[i][na + if sm == 0 { indices[i] } else { 0 }] = (u[i] - center) / scale;
    }
    let fit = ols(&x, &y)?;
    let mse = fit.sse / fit.df as f64;
    let mut paths = Vec::new();
    for (g, &identifier) in ids.iter().enumerate() {
        let si = na + if sm == 0 { g } else { 0 };
        let mut la = vec![0.; d];
        let mut lb = vec![0.; d];
        if na > 0 {
            la[if im == 0 { g } else { 0 }] = 1.;
        }
        la[si] = -center / scale;
        lb[si] = 1. / scale;
        let covariance = |a: &[f64], b: &[f64]| -> Result<f64, String> {
            finite(
                a.iter()
                    .enumerate()
                    .map(|(j, a)| a * dot(&fit.inverse[j], b))
                    .sum::<f64>()
                    * mse,
                "path covariance",
            )
        };
        paths.push(Path {
            id: identifier,
            a: finite(dot(&la, &fit.beta), "path intercept")?,
            b: finite(dot(&lb, &fit.beta), "path slope")?,
            caa: covariance(&la, &la)?,
            cab: covariance(&la, &lb)?,
            cbb: covariance(&lb, &lb)?,
            mse,
            df: fit.df as f64,
        });
    }
    Ok(Paths {
        paths,
        sse: fit.sse,
        df: fit.df,
        response,
        time,
        rank: d,
    })
}
fn path_model(paths: &Paths) -> Vec<f64> {
    let mut v = vec![
        1.,
        paths.response as f64,
        paths.time as f64,
        paths.paths.len() as f64,
    ];
    for p in &paths.paths {
        v.extend([p.id, p.a, p.b, p.caa, p.cab, p.cbb, p.mse, p.df]);
    }
    v
}
fn path_fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 4 {
        return Err(
            "path arguments are [intercept_mode,slope_mode,response_transform,time_transform]"
                .into(),
        );
    }
    let im = int(args.first().copied().unwrap_or(0.), 0, 2)?;
    let sm = int(args.get(1).copied().unwrap_or(0.), 0, 1)?;
    let response = int(args.get(2).copied().unwrap_or(0.), 0, 1)?;
    let time = int(args.get(3).copied().unwrap_or(0.), 0, 2)?;
    let p = fit_paths(data, im, sm, response, time)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"transformed_fixed_unit_linear_paths\",\"n_samples\":{},\"rank\":{},\"residual_df\":{},\"sse\":{},\"unit_ids\":{},\"intercepts\":{},\"slopes\":{},\"predict_parameters\":{},\"uncertainty_contract\":\"conditional_independent_homoskedastic_errors_not_serial_correlation_adjusted\",\"unsupported_bayesian_repeated_measures\":true}}",data.len(),p.rank,p.df,p.sse,arr(&p.paths.iter().map(|x|x.id).collect::<Vec<_>>()),arr(&p.paths.iter().map(|x|x.a).collect::<Vec<_>>()),arr(&p.paths.iter().map(|x|x.b).collect::<Vec<_>>()),arr(&path_model(&p))))
}
fn parse_paths(args: &[f64]) -> Result<Paths, String> {
    if args.len() < 4 || args[0] != 1. {
        return Err("path prediction requires version1 saved parameters".into());
    }
    let response = int(args[1], 0, 1)?;
    let time = int(args[2], 0, 2)?;
    let g = int(args[3], 1, 128)?;
    if args.len() != 4 + 8 * g {
        return Err("saved path model length is invalid".into());
    }
    let mut paths: Vec<Path> = Vec::new();
    for row in args[4..].chunks_exact(8) {
        id(row[0])?;
        if paths.iter().any(|p| p.id == row[0])
            || row[3] < 0.
            || row[5] < 0.
            || row[6] < 0.
            || row[7] < 1.
            || row[7] > 10000.
            || row[7].fract() != 0.
        {
            return Err("saved path covariance, residual variance, df or IDs are invalid".into());
        }
        let corr = if row[3] == 0. || row[5] == 0. {
            if row[4] == 0. {
                0.
            } else {
                f64::INFINITY
            }
        } else {
            row[4] / row[3].sqrt() / row[5].sqrt()
        };
        if corr.abs() > 1. + 1e-10 {
            return Err("saved path covariance must be positive semidefinite".into());
        }
        paths.push(Path {
            id: row[0],
            a: row[1],
            b: row[2],
            caa: row[3],
            cab: row[4],
            cbb: row[5],
            mse: row[6],
            df: row[7],
        });
    }
    Ok(Paths {
        paths,
        sse: 0.,
        df: 0,
        response,
        time,
        rank: 0,
    })
}
fn log_gamma(z: f64) -> f64 {
    let c = [
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
    let mut x = 0.9999999999998099;
    for (i, c) in c.iter().enumerate() {
        x += c / (zz + i as f64 + 1.);
    }
    let t = zz + 7.5;
    0.5 * LN_2PI + (zz + 0.5) * t.ln() - t + x.ln()
}
fn beta_cf(a: f64, b: f64, x: f64) -> f64 {
    let qab = a + b;
    let qap = a + 1.;
    let qam = a - 1.;
    let mut c = 1.;
    let mut d = 1. - qab * x / qap;
    if d.abs() < 1e-300 {
        d = 1e-300;
    }
    d = 1. / d;
    let mut h = d;
    for m in 1..=400 {
        let mf = m as f64;
        let m2 = 2. * mf;
        for aa in [
            mf * (b - mf) * x / ((qam + m2) * (a + m2)),
            -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2)),
        ] {
            d = 1. + aa * d;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = 1. + aa / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1. / d;
            let change = d * c;
            h *= change;
            if aa < 0. && (change - 1.).abs() < 2e-15 {
                return h;
            }
        }
    }
    h
}
fn beta_i(a: f64, b: f64, x: f64) -> f64 {
    if x <= 0. {
        return 0.;
    }
    if x >= 1. {
        return 1.;
    }
    let f = (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    if x < (a + 1.) / (a + b + 2.) {
        f * beta_cf(a, b, x) / a
    } else {
        1. - f * beta_cf(b, a, 1. - x) / b
    }
}
fn student_critical(confidence: f64, two_sided: bool, df: f64) -> f64 {
    let tail = (1. - confidence) / if two_sided { 2. } else { 1. };
    let sf = |t: f64| 0.5 * beta_i(df / 2., 0.5, df / (df + t * t));
    let mut hi = 1.;
    while sf(hi) > tail && hi < 1e150 {
        hi *= 2.;
    }
    let mut lo = 0.;
    for _ in 0..120 {
        let mid = lo / 2. + hi / 2.;
        if sf(mid) > tail {
            lo = mid
        } else {
            hi = mid
        }
    }
    lo / 2. + hi / 2.
}
fn f_p(f: f64, a: f64, b: f64) -> f64 {
    if f.is_infinite() {
        0.
    } else {
        beta_i(b / 2., a / 2., b / (b + a * f.max(0.)))
    }
}
fn roots(a: f64, b: f64, c: f64) -> Result<Vec<f64>, String> {
    let scale = a.abs().max(b.abs()).max(c.abs());
    if scale == 0. {
        return Ok(vec![]);
    }
    let (a, b, c) = (a / scale, b / scale, c / scale);
    if a == 0. {
        return if b == 0. {
            Ok(vec![])
        } else {
            Ok(vec![finite(-c / b, "linear inverse root")?])
        };
    }
    let disc = b * b - 4. * a * c;
    if disc < 0. {
        return Ok(vec![]);
    }
    if disc == 0. {
        return Ok(vec![finite(-b / (2. * a), "inverse root")?]);
    }
    let q = -0.5 * (b + if b >= 0. { disc.sqrt() } else { -disc.sqrt() });
    let mut out = vec![
        finite(q / a, "inverse root")?,
        finite(c / q, "inverse root")?,
    ];
    out.sort_by(f64::total_cmp);
    out.dedup();
    Ok(out)
}
fn variance(p: &Path, u: f64) -> Result<f64, String> {
    let v = finite(p.caa + 2. * p.cab * u + p.cbb * u * u, "mean variance")?;
    let scale = p.caa.abs() + (2. * p.cab * u).abs() + (p.cbb * u * u).abs();
    if v < -1e-12 * scale {
        return Err("negative mean prediction variance".into());
    }
    Ok(v.max(0.))
}
fn inverse_sets(p: &Path, target: f64, q: f64, side: usize, time: usize) -> Result<String, String> {
    let z = p.a - target;
    let (a, b, c) = (
        finite(p.b * p.b - q * q * p.cbb, "inverse quadratic")?,
        finite(2. * (z * p.b - q * q * p.cab), "inverse quadratic")?,
        finite(z * z - q * q * p.caa, "inverse quadratic")?,
    );
    let mut cuts = roots(a, b, c)?;
    if p.b != 0. {
        cuts.push(finite(-z / p.b, "inverse sign root")?);
    }
    let lower = if time == 1 { f64::NEG_INFINITY } else { 0. };
    cuts.retain(|x| *x > lower);
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    cuts.insert(0, lower);
    cuts.push(f64::INFINITY);
    let accepts = |u: f64| -> Result<bool, String> {
        let d = finite(z + p.b * u, "inverse mean difference")?;
        let s = q * variance(p, u)?.sqrt();
        let e = 64. * f64::EPSILON * (1. + d.abs() + s);
        Ok(match side {
            0 => d.abs() <= s + e,
            1 => d <= s + e,
            _ => d >= -s - e,
        })
    };
    let mut intervals: Vec<(f64, f64)> = Vec::new();
    for i in 0..cuts.len() - 1 {
        let (l, h) = (cuts[i], cuts[i + 1]);
        let mid = if l.is_infinite() {
            h - (1. + h.abs())
        } else if h.is_infinite() {
            l + (1. + l.abs())
        } else {
            l / 2. + h / 2.
        };
        if accepts(mid)? {
            if let Some(last) = intervals.last_mut() {
                if last.1 == l {
                    last.1 = h;
                    continue;
                }
            }
            intervals.push((l, h));
        }
    }
    for &x in &cuts {
        if x.is_finite() && accepts(x)? && !intervals.iter().any(|(l, h)| x >= *l && x <= *h) {
            intervals.push((x, x));
        }
    }
    intervals.sort_by(|a, b| a.0.total_cmp(&b.0));
    Ok(format!("[{}]",intervals.iter().map(|(l,h)|{let low=if l.is_infinite(){0.}else{inverse(*l,time)};let high=if h.is_infinite(){f64::INFINITY}else{inverse(*h,time)};format!("{{\"lower\":{},\"upper\":{},\"lower_closed\":{},\"upper_closed\":{},\"upper_unbounded\":{}}}",num(low),num(high),l.is_finite(),h.is_finite(),h.is_infinite())}).collect::<Vec<_>>().join(",")))
}
fn path_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 6 || data[0].len() != 3 {
        return Err(
            "path prediction rows [time,id,threshold], arguments [confidence,side,model...]".into(),
        );
    }
    let conf = args[0];
    if !(0.5..1.).contains(&conf) {
        return Err("confidence must be in (0.5,1)".into());
    }
    let side = int(args[1], 0, 2)?;
    let model = parse_paths(&args[2..])?;
    let mut fits = Vec::new();
    let mut lo = Vec::new();
    let mut hi = Vec::new();
    let mut pl = Vec::new();
    let mut ph = Vec::new();
    let mut times = Vec::new();
    let mut sets = Vec::new();
    for row in data {
        let p = model
            .paths
            .iter()
            .find(|p| p.id == row[1])
            .ok_or("unknown saved path ID")?;
        let u = transform(row[0], model.time, true)?;
        let c = transform(row[2], model.response, false)?;
        let q = student_critical(conf, side == 0, p.df);
        let mu = finite(p.a + p.b * u, "path prediction")?;
        let v = variance(p, u)?;
        let se = v.sqrt();
        let spe = finite(v + p.mse, "prediction variance")?.sqrt();
        fits.push(inverse(mu, model.response));
        lo.push(if side == 2 {
            f64::NEG_INFINITY
        } else {
            inverse(mu - q * se, model.response)
        });
        hi.push(if side == 1 {
            f64::INFINITY
        } else {
            inverse(mu + q * se, model.response)
        });
        pl.push(if side == 2 {
            f64::NEG_INFINITY
        } else {
            inverse(mu - q * spe, model.response)
        });
        ph.push(if side == 1 {
            f64::INFINITY
        } else {
            inverse(mu + q * spe, model.response)
        });
        let crossing = if p.b == 0. {
            f64::NAN
        } else {
            let x = (c - p.a) / p.b;
            if model.time != 1 && x < 0. {
                f64::NAN
            } else {
                inverse(x, model.time)
            }
        };
        times.push(crossing);
        sets.push(inverse_sets(p, c, q, side, model.time)?);
    }
    Ok(format!("{{\"version\":1,\"fitted\":{},\"mean_lower\":{},\"mean_upper\":{},\"prediction_lower\":{},\"prediction_upper\":{},\"threshold_times\":{},\"inverse_confidence_sets\":[{}],\"inverse_contract\":\"pointwise_mean_band_inversion_not_simultaneous_or_first_passage\",\"uncertainty_contract\":\"conditional_independent_errors_not_serial_correlation_adjusted\"}}",arr(&fits),arr(&lo),arr(&hi),arr(&pl),arr(&ph),arr(&times),sets.join(",")))
}
fn shelf_life(
    p: &Path,
    lower: Option<f64>,
    upper: Option<f64>,
    q: f64,
    horizon: f64,
) -> Result<Option<f64>, String> {
    let mut cuts = vec![0., horizon];
    for target in [lower, upper].into_iter().flatten() {
        let z = p.a - target;
        cuts.extend(
            roots(
                p.b * p.b - q * q * p.cbb,
                2. * (z * p.b - q * q * p.cab),
                z * z - q * q * p.caa,
            )?
            .into_iter()
            .filter(|x| *x > 0. && *x < horizon),
        );
    }
    cuts.sort_by(f64::total_cmp);
    cuts.dedup();
    let valid = |t: f64| -> Result<bool, String> {
        let m = finite(p.a + p.b * t, "stability mean")?;
        let e = q * variance(p, t)?.sqrt();
        Ok(lower.is_none_or(|l| m - e >= l) && upper.is_none_or(|u| m + e <= u))
    };
    if !valid(0.)? {
        return Ok(Some(0.));
    }
    for pair in cuts.windows(2) {
        if !valid(pair[0] / 2. + pair[1] / 2.)? {
            return Ok(Some(pair[0]));
        }
    }
    if !valid(horizon)? {
        return Ok(Some(horizon));
    }
    Ok(None)
}
fn stability(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 8 {
        return Err("stability arguments [lower_enabled,lower,upper_enabled,upper,pooling_alpha,confidence,pool_separate_MSE,horizon]".into());
    }
    let le = int(args[0], 0, 1)?;
    let ue = int(args[2], 0, 1)?;
    let pool = int(args[6], 0, 1)?;
    if le + ue == 0
        || le == 1 && ue == 1 && args[1] >= args[3]
        || !(0.01..0.5).contains(&args[4])
        || !(0.5..1.).contains(&args[5])
        || args[7] <= 0.
    {
        return Err("invalid stability specifications, alpha, confidence or horizon".into());
    }
    let full = fit_paths(data, 0, 0, 0, 0)?;
    let common_slope = fit_paths(data, 0, 1, 0, 0)?;
    let common = fit_paths(data, 1, 1, 0, 0)?;
    let g = full.paths.len();
    if g < 2 {
        return Err("stability pooling requires at least two batches".into());
    }
    let compare = |restricted: &Paths, unrestricted: &Paths| {
        let df = (unrestricted.rank - restricted.rank) as f64;
        let delta = (restricted.sse - unrestricted.sse).max(0.);
        let denominator = unrestricted.sse / unrestricted.df as f64;
        let f = if denominator == 0. {
            if delta == 0. {
                0.
            } else {
                f64::INFINITY
            }
        } else {
            delta / df / denominator
        };
        (f, f_p(f, df, unrestricted.df as f64))
    };
    let (fs, ps) = compare(&common_slope, &full);
    let (fi, pi) = compare(&common, &common_slope);
    let selected = if ps < args[4] {
        0
    } else if pi < args[4] {
        1
    } else {
        2
    };
    let mut separate = full.paths.clone();
    if pool == 0 {
        for p in &mut separate {
            let rows = data
                .iter()
                .filter(|r| r[2] == p.id)
                .cloned()
                .collect::<Vec<_>>();
            let own = fit_paths(&rows, 0, 0, 0, 0)?;
            *p = own.paths[0].clone();
        }
    }
    let reports = [separate, common_slope.paths.clone(), common.paths.clone()];
    let mut shelf = Vec::new();
    let lower = if le == 1 { Some(args[1]) } else { None };
    let upper = if ue == 1 { Some(args[3]) } else { None };
    for paths in &reports {
        let mut first = None;
        for p in paths {
            let q = student_critical(args[5], le + ue == 2, p.df);
            if let Some(t) = shelf_life(p, lower, upper, q, args[7])? {
                first = Some(first.map_or(t, |f: f64| f.min(t)));
            }
        }
        shelf.push(first.unwrap_or(f64::INFINITY));
    }
    let names = [
        "separate_intercepts_separate_slopes",
        "separate_intercepts_common_slope",
        "common_intercept_common_slope",
    ];
    Ok(format!("{{\"version\":1,\"selected_model\":\"{}\",\"test_order\":[\"slopes\",\"intercepts\"],\"pooling_alpha\":{},\"slope_pooling_f\":{},\"slope_pooling_p\":{},\"intercept_pooling_f\":{},\"intercept_pooling_p\":{},\"model_sse\":{},\"model_residual_df\":{},\"model_shelf_lives\":{},\"shelf_life\":{},\"beyond_horizon\":{},\"confidence_sides\":{},\"separate_model_pooled_mse\":{},\"uncertainty_contract\":\"normal_independent_batch_regressions_not_correlated_repeated_measure_errors\",\"regulatory_extrapolation_rules_applied\":false}}",names[selected],args[4],num(fs),ps,num(fi),pi,arr(&[full.sse,common_slope.sse,common.sse]),arr(&[full.df as f64,common_slope.df as f64,common.df as f64]),arr(&shelf),num(shelf[selected]),!shelf[selected].is_finite(),if le+ue==2{2}else{1},pool==1))
}

fn gamma_q_half(x: f64) -> f64 {
    if x == 0. {
        return 1.;
    }
    let a = 0.5;
    if x < a + 1. {
        let mut sum = 1. / a;
        let mut term = sum;
        for k in 1..400 {
            term *= x / (a + k as f64);
            sum += term;
            if term.abs() < sum.abs() * 2e-16 {
                break;
            }
        }
        (1. - sum * (-x + a * x.ln() - log_gamma(a)).exp()).clamp(0., 1.)
    } else {
        let mut b = x + 1. - a;
        let mut c = 1e300;
        let mut d = 1. / b;
        let mut h = d;
        for i in 1..400 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.;
            d = an * d + b;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = b + an / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1. / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-15 {
                break;
            }
        }
        (h * (-x + a * x.ln() - log_gamma(a)).exp()).clamp(0., 1.)
    }
}
fn log_phi_cdf(z: f64) -> f64 {
    if z <= -8. {
        let x = -z;
        let mut tail = 0.;
        for k in (1..=200).rev() {
            tail = k as f64 / (x + tail);
        }
        return -0.5 * z * z - 0.5 * LN_2PI - (x + tail).ln();
    }
    if z >= 8. {
        return (-log_phi_cdf(-z).exp()).ln_1p();
    }
    let q = 0.5 * gamma_q_half(0.5 * z * z);
    if z < 0. {
        q.ln()
    } else {
        (-q).ln_1p()
    }
}
fn normal_quantile(p: f64) -> f64 {
    let mut lo = -40.;
    let mut hi = 40.;
    for _ in 0..120 {
        let m = (lo + hi) / 2.;
        if log_phi_cdf(m) < p.ln() {
            lo = m
        } else {
            hi = m
        }
    }
    (lo + hi) / 2.
}
fn log_diff(a: f64, b: f64) -> Result<f64, String> {
    if b == f64::NEG_INFINITY {
        return Ok(a);
    }
    if a <= b {
        return Err("censor interval probability is numerically unresolved".into());
    }
    finite(
        a + (-(b - a).exp_m1()).ln(),
        "log censor interval probability",
    )
}
// Normal response interval probability and scores for location and log scale.
fn normal_term(
    lo: f64,
    hi: f64,
    kind: usize,
    mu: f64,
    sigma: f64,
) -> Result<(f64, f64, f64), String> {
    let pdf = |z: f64| -0.5 * z * z - 0.5 * LN_2PI;
    let (l, m, s) = match kind {
        0 => {
            let z = (lo - mu) / sigma;
            (-sigma.ln() + pdf(z), z / sigma, z * z - 1.)
        }
        1 => {
            let z = (hi - mu) / sigma;
            let l = log_phi_cdf(z);
            let r = (pdf(z) - l).exp();
            (l, -r / sigma, -z * r)
        }
        2 => {
            let z = (lo - mu) / sigma;
            let l = log_phi_cdf(-z);
            let r = (pdf(z) - l).exp();
            (l, r / sigma, z * r)
        }
        _ => {
            let a = (lo - mu) / sigma;
            let b = (hi - mu) / sigma;
            let l = if a > 0. {
                log_diff(log_phi_cdf(-a), log_phi_cdf(-b))?
            } else {
                log_diff(log_phi_cdf(b), log_phi_cdf(a))?
            };
            let ra = (pdf(a) - l).exp();
            let rb = (pdf(b) - l).exp();
            (l, (ra - rb) / sigma, a * ra - b * rb)
        }
    };
    Ok((
        finite(l, "normal log likelihood")?,
        finite(m, "location score")?,
        finite(s, "scale score")?,
    ))
}
struct Opt {
    theta: Vec<f64>,
    converged: bool,
    guard_bound: bool,
    iterations: usize,
    gradient: f64,
    line_failed: bool,
    value: f64,
}
#[allow(clippy::needless_range_loop)]
fn optimize<F>(
    initial: Vec<f64>,
    bounds: &[(f64, f64)],
    maxiter: usize,
    tol: f64,
    eval: &F,
) -> Result<Opt, String>
where
    F: Fn(&[f64]) -> Result<(f64, Vec<f64>), String>,
{
    let d = initial.len();
    let mut x = initial;
    let (mut f, mut g) = eval(&x)?;
    let mut h = vec![vec![0.; d]; d];
    for i in 0..d {
        h[i][i] = 1.;
    }
    let mut iterations = 0;
    let mut failed = false;
    for k in 0..maxiter {
        iterations = k;
        let projected = g
            .iter()
            .enumerate()
            .map(|(i, g)| {
                if x[i] <= bounds[i].0 && *g > 0. || x[i] >= bounds[i].1 && *g < 0. {
                    0.
                } else {
                    g.abs()
                }
            })
            .fold(0_f64, f64::max);
        if projected <= tol {
            break;
        }
        let mut direction = h.iter().map(|r| -dot(r, &g)).collect::<Vec<_>>();
        if dot(&direction, &g) >= -1e-16 * dot(&g, &g) {
            direction = g.iter().map(|v| -v).collect();
            h = vec![vec![0.; d]; d];
            for i in 0..d {
                h[i][i] = 1.;
            }
        }
        let mut accepted = None;
        let mut step = 1.;
        for _ in 0..50 {
            let candidate = x
                .iter()
                .enumerate()
                .map(|(j, v)| (v + step * direction[j]).clamp(bounds[j].0, bounds[j].1))
                .collect::<Vec<_>>();
            let displacement = candidate
                .iter()
                .zip(&x)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            if let Ok((nf, ng)) = eval(&candidate) {
                if nf <= f + 1e-4 * dot(&g, &displacement) && dot(&g, &displacement) < 0. {
                    accepted = Some((candidate, nf, ng, displacement));
                    break;
                }
            }
            step *= 0.5;
        }
        let Some((nx, nf, ng, s)) = accepted else {
            failed = true;
            break;
        };
        let y = ng.iter().zip(&g).map(|(a, b)| a - b).collect::<Vec<_>>();
        let ys = dot(&y, &s);
        if ys
            > 1e-12
                * s.iter().fold(0_f64, |a, b| a.hypot(*b))
                * y.iter().fold(0_f64, |a, b| a.hypot(*b))
        {
            let hy = h.iter().map(|r| dot(r, &y)).collect::<Vec<_>>();
            let yhy = dot(&y, &hy);
            for i in 0..d {
                for j in 0..d {
                    h[i][j] +=
                        (1. + yhy / ys) * s[i] * s[j] / ys - (s[i] * hy[j] + hy[i] * s[j]) / ys;
                }
            }
        }
        x = nx;
        f = nf;
        g = ng;
        iterations = k + 1;
    }
    let gradient = g
        .iter()
        .enumerate()
        .map(|(i, g)| {
            if x[i] <= bounds[i].0 && *g > 0. || x[i] >= bounds[i].1 && *g < 0. {
                0.
            } else {
                g.abs()
            }
        })
        .fold(0_f64, f64::max);
    let guard = x
        .iter()
        .zip(bounds)
        .any(|(x, (l, u))| *x <= *l + 1e-8 || *x >= *u - 1e-8);
    Ok(Opt {
        theta: x,
        converged: gradient <= tol && !guard,
        guard_bound: guard,
        iterations,
        gradient,
        line_failed: failed,
        value: f,
    })
}
#[allow(clippy::needless_range_loop)]
fn inverse_spd(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let d = a.len();
    let mut l = vec![vec![0.; d]; d];
    let scale = (0..d).map(|i| a[i][i].abs()).fold(0_f64, f64::max);
    for i in 0..d {
        for j in 0..=i {
            let v = (a[i][j] + a[j][i]) / 2. - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if v <= 1e-12 * scale {
                    return Err("observed information is not identified positive definite".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut out = vec![vec![0.; d]; d];
    for col in 0..d {
        let mut y = vec![0.; d];
        for i in 0..d {
            y[i] = ((if i == col { 1. } else { 0. }) - (0..i).map(|j| l[i][j] * y[j]).sum::<f64>())
                / l[i][i];
        }
        for i in (0..d).rev() {
            out[i][col] = (y[i] - (i + 1..d).map(|j| l[j][i] * out[j][col]).sum::<f64>()) / l[i][i];
            finite(out[i][col], "information inverse")?;
        }
    }
    Ok(out)
}
fn observed_covariance<F>(x: &[f64], n: usize, eval: &F) -> Result<Vec<Vec<f64>>, String>
where
    F: Fn(&[f64]) -> Result<(f64, Vec<f64>), String>,
{
    let d = x.len();
    let mut h = vec![vec![0.; d]; d];
    for j in 0..d {
        let step = 1e-4 * (1. + x[j].abs());
        let mut a = x.to_vec();
        let mut b = x.to_vec();
        a[j] += step;
        b[j] -= step;
        let ga = eval(&a)?.1;
        let gb = eval(&b)?.1;
        for i in 0..d {
            h[i][j] = finite(
                (ga[i] - gb[i]) / (2. * step) * n as f64,
                "observed information",
            )?;
        }
    }
    inverse_spd(&h)
}
fn mat(a: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| arr(r)).collect::<Vec<_>>().join(",")
    )
}
fn stress(x: f64, kind: usize) -> Result<f64, String> {
    finite(
        match kind {
            0 | 1 => x,
            2 => {
                if x <= -273.15 {
                    return Err("Arrhenius Celsius stress must exceed absolute zero".into());
                }
                11604.5181215503 / (x + 273.15)
            }
            3 => {
                if x <= 0. {
                    return Err("log stress must be positive".into());
                }
                x.ln()
            }
            _ => return Err("unknown stress relation".into()),
        },
        "stress transform",
    )
}
fn destructive_features(t: f64, s: f64, time: usize, mode: usize) -> Result<Vec<f64>, String> {
    let u = transform(t, time, true)?;
    if mode == 0 {
        Ok(vec![u])
    } else {
        let z = stress(s, mode)?;
        Ok(vec![u, z, finite(u * z, "time stress interaction")?])
    }
}
fn censor_bounds(row: &[f64], family: usize) -> Result<(f64, f64, usize), String> {
    let kind = int(row[4], 0, 3)?;
    let (lo, hi) = (row[2], row[3]);
    if kind == 0 && lo != hi
        || kind == 1 && lo != 0.
        || kind == 2 && hi != 0.
        || kind == 3 && lo >= hi
    {
        return Err("invalid exact/left/right/interval response bounds".into());
    }
    let tr = |x: f64| -> Result<f64, String> {
        if family == 1 {
            if x <= 0. {
                return Err("lognormal finite response bounds must be positive".into());
            }
            Ok(x.ln())
        } else {
            Ok(x)
        }
    };
    Ok((
        if kind == 1 { 0. } else { tr(lo)? },
        if kind == 2 { 0. } else { tr(hi)? },
        kind,
    ))
}
fn destructive_fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 5 || args.len() > 5 || data.len() > 4000 {
        return Err(
            "destructive rows [time,stress,lower,upper,censor_type] with <=4000 observations"
                .into(),
        );
    }
    let family = int(args.first().copied().unwrap_or(0.), 0, 1)?;
    let time = int(args.get(1).copied().unwrap_or(0.), 0, 2)?;
    let mode = int(args.get(2).copied().unwrap_or(0.), 0, 3)?;
    let maxiter = int(args.get(3).copied().unwrap_or(300.), 1, 1000)?;
    let tol = args.get(4).copied().unwrap_or(1e-8);
    if !(1e-12..=1e-3).contains(&tol) {
        return Err(
            "destructive args [family,time_transform,stress_mode,maxiter,tolerance]".into(),
        );
    }
    let raw = data
        .iter()
        .map(|r| destructive_features(r[0], r[1], time, mode))
        .collect::<Result<Vec<_>, _>>()?;
    let p = raw[0].len();
    let d = p + 1;
    let origin = raw[0].clone();
    let offsets = (0..p)
        .map(|j| {
            finite(
                raw.iter()
                    .map(|r| (r[j] - origin[j]) / raw.len() as f64)
                    .sum(),
                "feature mean",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let scales = (0..p)
        .map(|j| {
            let s = raw
                .iter()
                .fold(0_f64, |a, r| a.hypot((r[j] - origin[j]) - offsets[j]))
                / (data.len() as f64).sqrt();
            if s == 0. {
                Err("destructive fixed design contains constant feature".into())
            } else {
                finite(s, "feature scale")
            }
        })
        .collect::<Result<Vec<_>, String>>()?;
    let x = raw
        .iter()
        .map(|r| {
            let mut row = vec![1.];
            for j in 0..p {
                row.push(((r[j] - origin[j]) - offsets[j]) / scales[j]);
            }
            row
        })
        .collect::<Vec<_>>();
    let bounds = data
        .iter()
        .map(|r| censor_bounds(r, family))
        .collect::<Result<Vec<_>, _>>()?;
    let exact = bounds
        .iter()
        .enumerate()
        .filter(|(_, b)| b.2 == 0)
        .map(|(i, _)| i)
        .collect::<Vec<_>>();
    if exact.len() <= d {
        return Err(
            "destructive fit requires more exact responses than location coefficients".into(),
        );
    }
    let ex = exact.iter().map(|i| x[*i].clone()).collect::<Vec<_>>();
    let ey = exact.iter().map(|i| bounds[*i].0).collect::<Vec<_>>();
    let init = ols(&ex, &ey)?;
    let sigma = (init.sse / exact.len() as f64).sqrt().max(1e-3);
    let mut theta = init.beta;
    theta.push(sigma.ln());
    let budget = std::cell::Cell::new(0_usize);
    let eval = |b: &[f64]| -> Result<(f64, Vec<f64>), String> {
        budget.set(budget.get().saturating_add(data.len().saturating_mul(d)));
        if budget.get() > MAX_WORK {
            return Err("destructive likelihood evaluation budget exhausted".into());
        }
        let sigma = b[d].exp();
        let mut ll = 0.;
        let mut g = vec![0.; d + 1];
        for (i, (l, h, kind)) in bounds.iter().enumerate() {
            let (v, mu, ls) = normal_term(*l, *h, *kind, dot(&x[i], &b[..d]), sigma)?;
            ll += v - if family == 1 && *kind == 0 {
                data[i][2].ln()
            } else {
                0.
            };
            for j in 0..d {
                g[j] -= mu * x[i][j] / data.len() as f64;
            }
            g[d] -= ls / data.len() as f64;
        }
        if g.iter().any(|v| !v.is_finite()) {
            return Err("destructive likelihood score exceeds finite range".into());
        }
        Ok((
            finite(-ll / data.len() as f64, "negative log likelihood")?,
            g,
        ))
    };
    let mut search = vec![(-1e6, 1e6); d + 1];
    search[d] = (-30., 30.);
    let result = optimize(theta, &search, maxiter, tol, &eval)?;
    let sigma = result.theta[d].exp();
    let slopes = (0..p)
        .map(|j| result.theta[j + 1] / scales[j])
        .collect::<Vec<_>>();
    let intercept = finite(
        result.theta[0]
            - (0..p)
                .map(|j| slopes[j] * (origin[j] + offsets[j]))
                .sum::<f64>(),
        "destructive raw intercept",
    )?;
    let mut coeff = vec![intercept];
    coeff.extend(&slopes);
    let covariance = observed_covariance(&result.theta, data.len(), &eval).ok();
    if budget.get() > MAX_WORK {
        return Err("destructive likelihood evaluation budget exhausted".into());
    }
    let information_identified = covariance.is_some();
    let optimizer_converged = result.gradient <= tol;
    let parameter_guard_bound = result.guard_bound;
    let inference_valid = result.converged && information_identified;
    let converged = inference_valid;
    let covariance_json = if inference_valid {
        covariance.as_ref().map_or("null".into(), |c| mat(c))
    } else {
        "null".into()
    };
    let aic = if inference_valid {
        num(2. * result.value * data.len() as f64 + 2. * (d + 1) as f64)
    } else {
        "null".into()
    };
    let mut portable = vec![
        1.,
        family as f64,
        time as f64,
        mode as f64,
        sigma,
        result.theta[0],
        p as f64,
    ];
    portable.extend(&origin);
    portable.extend(&offsets);
    portable.extend(&slopes);
    if portable.iter().any(|v| !v.is_finite()) {
        return Err("destructive model is not finite".into());
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"independent_destructive_normal_location_scale_ml\",\"family\":\"{}\",\"coefficients\":{},\"sigma\":{sigma},\"log_likelihood\":{},\"aic\":{},\"n_exact\":{},\"n_samples\":{},\"converged\":{},\"optimizer_converged\":{optimizer_converged},\"parameter_guard_bound\":{parameter_guard_bound},\"information_identified\":{information_identified},\"inference_valid\":{inference_valid},\"iterations\":{},\"gradient_max\":{},\"line_search_failed\":{},\"standardized_parameter_covariance\":{},\"covariance_axes\":\"centered_intercept,standardized_slopes,log_sigma\",\"predict_parameters\":{},\"uncertainty_contract\":\"independent_units_not_repeated_unit_errors\"}}",if family==0{"normal"}else{"lognormal"},arr(&coeff),-result.value*data.len() as f64,aic,exact.len(),data.len(),converged,result.iterations,result.gradient,result.line_failed,covariance_json,arr(&portable)))
}
fn destructive_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 7 || args[0] != 1. || data[0].len() != 3 {
        return Err(
            "destructive prediction rows [time,stress,threshold] with unchanged saved parameters"
                .into(),
        );
    }
    let family = int(args[1], 0, 1)?;
    let time = int(args[2], 0, 2)?;
    let mode = int(args[3], 0, 3)?;
    let sigma = args[4];
    let p = int(args[6], 1, 3)?;
    if p != if mode == 0 { 1 } else { 3 } || args.len() != 7 + 3 * p || sigma <= 0. {
        return Err("invalid destructive saved model shape/scale".into());
    }
    let origins = &args[7..7 + p];
    let offsets = &args[7 + p..7 + 2 * p];
    let slopes = &args[7 + 2 * p..];
    let mut medians = Vec::new();
    let mut means = Vec::new();
    let mut lower = Vec::new();
    let mut upper = Vec::new();
    let mut probs = Vec::new();
    for row in data {
        let features = destructive_features(row[0], row[1], time, mode)?;
        let mu = finite(
            args[5]
                + (0..p)
                    .map(|j| slopes[j] * ((features[j] - origins[j]) - offsets[j]))
                    .sum::<f64>(),
            "destructive mean predictor",
        )?;
        let z = if family == 1 {
            if row[2] <= 0. {
                f64::NEG_INFINITY
            } else {
                (row[2].ln() - mu) / sigma
            }
        } else {
            (row[2] - mu) / sigma
        };
        probs.push(log_phi_cdf(z).exp());
        medians.push(if family == 0 { mu } else { mu.exp() });
        means.push(if family == 0 {
            mu
        } else {
            (mu + sigma * sigma / 2.).exp()
        });
        lower.push(if family == 0 {
            mu - 1.959963984540054 * sigma
        } else {
            (mu - 1.959963984540054 * sigma).exp()
        });
        upper.push(if family == 0 {
            mu + 1.959963984540054 * sigma
        } else {
            (mu + 1.959963984540054 * sigma).exp()
        });
    }
    Ok(format!("{{\"version\":1,\"medians\":{},\"means\":{},\"distribution_lower_025\":{},\"distribution_upper_975\":{},\"lower_failure_probabilities\":{},\"failure_contract\":\"marginal_response_below_threshold_not_first_passage_lifetime\"}}",arr(&medians),arr(&means),arr(&lower),arr(&upper),arr(&probs)))
}

fn damage_stress(x: f64, kind: usize) -> Result<f64, String> {
    match kind {
        0 => stress(x, 1),
        1 => stress(x, 3),
        2 => stress(x, 2),
        3 => {
            if x <= 0. {
                Err("reciprocal stress must be positive".into())
            } else {
                finite(1. / x, "reciprocal stress")
            }
        }
        _ => Err("unknown damage stress relation".into()),
    }
}
#[derive(Clone)]
struct Pattern {
    id: f64,
    segments: Vec<(f64, f64)>,
    duration: f64,
}
fn parse_patterns(args: &[f64], start: usize, relation: usize) -> Result<Vec<Pattern>, String> {
    if start >= args.len() {
        return Err("stress patterns are missing".into());
    }
    let count = int(args[start], 1, 128)?;
    let mut cursor = start + 1;
    let mut patterns: Vec<Pattern> = Vec::new();
    let mut total = 0;
    for _ in 0..count {
        if cursor + 2 > args.len() {
            return Err("truncated stress pattern header".into());
        }
        let identifier = args[cursor];
        id(identifier)?;
        if patterns.iter().any(|p| p.id == identifier) {
            return Err("stress pattern IDs must be unique".into());
        }
        let n = int(args[cursor + 1], 1, 128)?;
        cursor += 2;
        total += n;
        if total > 1024 || cursor + 2 * n > args.len() {
            return Err("truncated or excessive stress segments".into());
        }
        let mut segments = Vec::new();
        let mut duration = 0.;
        for _ in 0..n {
            let dt = args[cursor];
            let st = args[cursor + 1];
            if dt <= 0. {
                return Err("stress segment durations must be positive".into());
            }
            damage_stress(st, relation)?;
            duration = finite(duration + dt, "stress pattern duration")?;
            segments.push((dt, st));
            cursor += 2;
        }
        patterns.push(Pattern {
            id: identifier,
            segments,
            duration,
        });
    }
    if cursor != args.len() {
        return Err("unexpected trailing stress pattern arguments".into());
    }
    Ok(patterns)
}
fn append_patterns(v: &mut Vec<f64>, patterns: &[Pattern]) {
    v.push(patterns.len() as f64);
    for p in patterns {
        v.extend([p.id, p.segments.len() as f64]);
        for (d, s) in &p.segments {
            v.extend([*d, *s]);
        }
    }
}
// Cumulative exposure A(t)=integral exp[-mu(stress(s))] ds. A(t) is
// continuous; the rate at an interior jump uses its right-hand stress level.
#[allow(clippy::too_many_arguments)]
fn exposure(
    pattern: &Pattern,
    t: f64,
    b0: f64,
    b1: f64,
    relation: usize,
    continuation: usize,
    center: f64,
    scale: f64,
) -> Result<(f64, f64, f64), String> {
    if t < 0. {
        return Err("event/prediction time must be nonnegative".into());
    }
    if continuation == 0 && t > pattern.duration {
        return Err("event time exceeds a stopped stress pattern".into());
    }
    let mut logs = Vec::new();
    let mut zs = Vec::new();
    let mut remaining = t;
    let mut cycles = 0.;
    if continuation == 2 && remaining >= pattern.duration {
        cycles = (remaining / pattern.duration).floor();
        remaining = remaining.rem_euclid(pattern.duration);
        finite(cycles, "stress pattern repetitions")?;
    }
    let zvalue = |s: f64| -> Result<f64, String> {
        finite(
            (damage_stress(s, relation)? - center) / scale,
            "standardized stress",
        )
    };
    if cycles > 0. {
        for (dt, s) in &pattern.segments {
            let z = zvalue(*s)?;
            logs.push(finite(
                dt.ln() + cycles.ln() - b0 - b1 * z,
                "cumulative log exposure term",
            )?);
            zs.push(z);
        }
    }
    let mut current = pattern.segments.last().unwrap().1;
    for (dt, s) in &pattern.segments {
        if remaining < *dt {
            current = *s;
        }
        let used = remaining.min(*dt);
        if used > 0. {
            let z = zvalue(*s)?;
            logs.push(finite(
                used.ln() - b0 - b1 * z,
                "cumulative log exposure term",
            )?);
            zs.push(z);
        }
        remaining = (remaining - *dt).max(0.);
        if used < *dt {
            break;
        }
    }
    if remaining > 0. {
        let z = zvalue(current)?;
        logs.push(finite(
            remaining.ln() - b0 - b1 * z,
            "extended log exposure term",
        )?);
        zs.push(z);
    }
    // A zero remaining time at an interior boundary belongs to the next segment.
    if continuation == 2 && t >= pattern.duration && t.rem_euclid(pattern.duration) == 0. {
        current = pattern.segments[0].1;
    }
    let current_z = zvalue(current)?;
    if logs.is_empty() {
        return Ok((f64::NEG_INFINITY, 0., current_z));
    }
    let maximum = logs.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let weights = logs.iter().map(|x| (x - maximum).exp()).collect::<Vec<_>>();
    let sum = weights.iter().sum::<f64>();
    Ok((
        finite(maximum + sum.ln(), "log cumulative exposure")?,
        dot(&weights, &zs) / sum,
        current_z,
    ))
}
struct Dist {
    lf: f64,
    ls: f64,
    ld: f64,
    fa: f64,
    fs: f64,
    sa: f64,
    ss: f64,
    da: f64,
    ds: f64,
}
fn baseline(log_age: f64, distribution: usize, log_sigma: f64) -> Result<Dist, String> {
    if log_age == f64::NEG_INFINITY {
        let density = if distribution == 0 || distribution == 1 && log_sigma == 0. {
            0.
        } else if distribution == 1 && log_sigma > 0. {
            f64::INFINITY
        } else {
            f64::NEG_INFINITY
        };
        return Ok(Dist {
            lf: f64::NEG_INFINITY,
            ls: 0.,
            ld: density,
            fa: 0.,
            fs: 0.,
            sa: 0.,
            ss: 0.,
            da: 0.,
            ds: 0.,
        });
    }
    let sigma = if distribution == 0 {
        1.
    } else {
        log_sigma.exp()
    };
    let z = log_age / sigma;
    let result = if distribution < 2 {
        let h = finite(z.exp(), "cumulative hazard")?;
        let lf = if h == 0. { z } else { (-(-h).exp_m1()).ln() };
        let ratio = if h < 1e-5 {
            1. - h / 2. + h * h / 12.
        } else if h > 700. {
            0.
        } else {
            h / h.exp_m1()
        };
        Dist {
            lf,
            ls: -h,
            ld: -sigma.ln() + z - log_age - h,
            fa: ratio / sigma,
            fs: -z * ratio,
            sa: -h / sigma,
            ss: h * z,
            da: (1. - h) / sigma - 1.,
            ds: -1. - z * (1. - h),
        }
    } else {
        let pdf = -0.5 * z * z - 0.5 * LN_2PI;
        let lf = log_phi_cdf(z);
        let ls = log_phi_cdf(-z);
        let rf = (pdf - lf).exp();
        let rs = (pdf - ls).exp();
        Dist {
            lf,
            ls,
            ld: -sigma.ln() - log_age + pdf,
            fa: rf / sigma,
            fs: -z * rf,
            sa: -rs / sigma,
            ss: z * rs,
            da: -z / sigma - 1.,
            ds: z * z - 1.,
        }
    };
    if [
        result.ls, result.ld, result.fa, result.fs, result.sa, result.ss, result.da, result.ds,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err("cumulative-exposure density/score is not representable".into());
    }
    Ok(result)
}
fn damage_bounds(row: &[f64]) -> Result<(f64, f64, usize), String> {
    if row.len() != 4 {
        return Err("damage rows [lower_time,upper_time,censor_type,pattern_id]".into());
    }
    let k = int(row[2], 0, 3)?;
    let (l, h) = (row[0], row[1]);
    id(row[3])?;
    if l < 0.
        || h < 0.
        || k == 0 && (l <= 0. || l != h)
        || k == 1 && (l != 0. || h <= 0.)
        || k == 2 && (l <= 0. || h != 0.)
        || k == 3 && (l >= h || h <= 0.)
    {
        return Err("invalid cumulative-damage exact/left/right/interval bounds".into());
    }
    Ok((l, h, k))
}
#[allow(clippy::too_many_arguments)]
fn damage_term(
    pattern: &Pattern,
    l: f64,
    h: f64,
    kind: usize,
    theta: &[f64],
    distribution: usize,
    relation: usize,
    continuation: usize,
    center: f64,
    scale: f64,
) -> Result<(f64, Vec<f64>), String> {
    let logs = if distribution == 0 { 0. } else { theta[2] };
    let d = theta.len();
    let evaluate = |t: f64| -> Result<(Dist, f64, f64), String> {
        let (age, average, current) = exposure(
            pattern,
            t,
            theta[0],
            theta[1],
            relation,
            continuation,
            center,
            scale,
        )?;
        Ok((baseline(age, distribution, logs)?, average, current))
    };
    let mut score = vec![0.; d];
    let ll = match kind {
        0 => {
            let (v, average, current) = evaluate(l)?;
            score[0] = -v.da - 1.;
            score[1] = -v.da * average - current;
            if d == 3 {
                score[2] = v.ds;
            }
            v.ld - theta[0] - theta[1] * current
        }
        1 => {
            let (v, average, _) = evaluate(h)?;
            score[0] = -v.fa;
            score[1] = -v.fa * average;
            if d == 3 {
                score[2] = v.fs;
            }
            v.lf
        }
        2 => {
            let (v, average, _) = evaluate(l)?;
            score[0] = -v.sa;
            score[1] = -v.sa * average;
            if d == 3 {
                score[2] = v.ss;
            }
            v.ls
        }
        _ => {
            let (a, za, _) = evaluate(l)?;
            let (b, zb, _) = evaluate(h)?;
            if a.lf > -(2.0_f64).ln() {
                let lp = log_diff(a.ls, b.ls)?;
                let wa = (a.ls - lp).exp();
                let wb = (b.ls - lp).exp();
                score[0] = -wa * a.sa + wb * b.sa;
                score[1] = -wa * a.sa * za + wb * b.sa * zb;
                if d == 3 {
                    score[2] = wa * a.ss - wb * b.ss;
                }
                lp
            } else {
                let lp = log_diff(b.lf, a.lf)?;
                let wa = (a.lf - lp).exp();
                let wb = (b.lf - lp).exp();
                score[0] = -wb * b.fa + wa * a.fa;
                score[1] = -wb * b.fa * zb + wa * a.fa * za;
                if d == 3 {
                    score[2] = wb * b.fs - wa * a.fs;
                }
                lp
            }
        }
    };
    if score.iter().any(|v| !v.is_finite()) {
        return Err("damage likelihood score exceeds numeric range".into());
    }
    Ok((finite(ll, "damage log likelihood")?, score))
}
fn damage_fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 6 || data[0].len() != 4 || data.len() > 4000 {
        return Err("damage fit args [distribution,relation,continuation,maxiter,tolerance,pattern_count,patterns...]".into());
    }
    let distribution = int(args[0], 0, 2)?;
    let relation = int(args[1], 0, 3)?;
    let continuation = int(args[2], 0, 2)?;
    let maxiter = int(args[3], 1, 1000)?;
    let tol = args[4];
    if !(1e-12..=1e-3).contains(&tol) {
        return Err("damage tolerance must be in [1e-12,1e-3]".into());
    }
    let patterns = parse_patterns(args, 5, relation)?;
    let bounds = data
        .iter()
        .map(|r| damage_bounds(r))
        .collect::<Result<Vec<_>, _>>()?;
    let groups = data
        .iter()
        .map(|r| {
            patterns
                .iter()
                .position(|p| p.id == r[3])
                .ok_or("unknown damage pattern ID")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let exact = bounds.iter().filter(|b| b.2 == 0).count();
    let d = if distribution == 0 { 2 } else { 3 };
    if data.len() <= d || exact < d {
        return Err(
            "damage fit requires n>parameters and at least parameters exact failures".into(),
        );
    }
    let mut transformed = Vec::new();
    for p in &patterns {
        for (_, s) in &p.segments {
            transformed.push(damage_stress(*s, relation)?);
        }
    }
    let center = transformed.iter().sum::<f64>() / transformed.len() as f64;
    let scale = transformed
        .iter()
        .map(|s| (s - center).abs())
        .fold(0_f64, f64::max);
    if !center.is_finite() || scale == 0. || !scale.is_finite() {
        return Err("at least two finite distinct stress levels are required".into());
    }
    let units = data
        .len()
        .saturating_mul(patterns.iter().map(|p| p.segments.len()).max().unwrap_or(1))
        .saturating_mul(maxiter)
        .saturating_mul(60)
        .saturating_mul(4);
    if units > MAX_WORK {
        return Err("damage fit exceeds nominal likelihood work bound".into());
    }
    let mut calls = 0_usize;
    let budget = std::cell::Cell::new(0_usize);
    let eval = |theta: &[f64]| -> Result<(f64, Vec<f64>), String> {
        let cost = data.len() * patterns.iter().map(|p| p.segments.len()).max().unwrap_or(1);
        budget.set(budget.get().saturating_add(cost));
        if budget.get() > MAX_WORK {
            return Err("damage likelihood evaluation budget exhausted".into());
        }
        let mut ll = 0.;
        let mut g = vec![0.; d];
        for (i, (l, h, k)) in bounds.iter().enumerate() {
            let (v, s) = damage_term(
                &patterns[groups[i]],
                *l,
                *h,
                *k,
                theta,
                distribution,
                relation,
                continuation,
                center,
                scale,
            )?;
            ll += v;
            for j in 0..d {
                g[j] -= s[j] / data.len() as f64;
            }
        }
        Ok((
            finite(
                -ll / data.len() as f64,
                "damage average negative likelihood",
            )?,
            g,
        ))
    };
    let mean_time = bounds
        .iter()
        .map(|(l, h, k)| if *k == 1 { *h } else { *l })
        .sum::<f64>()
        / data.len() as f64;
    let mut search = vec![(-1e4, 1e4); d];
    if d == 3 {
        search[2] = (-10., 10.);
    }
    let mut best = None;
    for slope in [0., -0.5, 0.5] {
        let mut initial = vec![mean_time.ln(), slope];
        if d == 3 {
            initial.push(0.);
        }
        if let Ok(result) = optimize(initial, &search, maxiter, tol, &eval) {
            calls += 1;
            if best.as_ref().is_none_or(|b: &Opt| result.value < b.value) {
                best = Some(result);
            }
        }
    }
    let best = best.ok_or("no finite cumulative-damage starting fit")?;
    if budget.get() > MAX_WORK {
        return Err("damage fit exhausted likelihood evaluation budget".into());
    }
    let raw_slope = finite(best.theta[1] / scale, "damage stress slope")?;
    let raw_intercept = finite(
        best.theta[0] - raw_slope * center,
        "damage location intercept",
    )?;
    let sigma = if distribution == 0 {
        1.
    } else {
        best.theta[2].exp()
    };
    let covariance = observed_covariance(&best.theta, data.len(), &eval).ok();
    if budget.get() > MAX_WORK {
        return Err("damage fit exhausted likelihood evaluation budget".into());
    }
    let information_identified = covariance.is_some();
    let optimizer_converged = best.gradient <= tol;
    let parameter_guard_bound = best.guard_bound;
    let inference_valid = best.converged && information_identified;
    let converged = inference_valid;
    let covariance_json = if inference_valid {
        covariance.as_ref().map_or("null".into(), |c| mat(c))
    } else {
        "null".into()
    };
    let aic = if inference_valid {
        num(2. * best.value * data.len() as f64 + 2. * d as f64)
    } else {
        "null".into()
    };
    let mut portable = vec![
        1.,
        distribution as f64,
        relation as f64,
        continuation as f64,
        raw_intercept,
        raw_slope,
        sigma,
    ];
    append_patterns(&mut portable, &patterns);
    Ok(format!("{{\"version\":1,\"algorithm\":\"step_stress_cumulative_exposure_ml\",\"coefficients\":{},\"sigma\":{sigma},\"log_likelihood\":{},\"aic\":{},\"n_samples\":{},\"n_exact\":{exact},\"converged\":{},\"optimizer_converged\":{optimizer_converged},\"parameter_guard_bound\":{parameter_guard_bound},\"information_identified\":{information_identified},\"inference_valid\":{inference_valid},\"iterations\":{},\"gradient_max\":{},\"line_search_failed\":{},\"successful_starts\":{calls},\"likelihood_work_units\":{},\"stress_center\":{center},\"stress_scale\":{scale},\"centered_parameter_covariance\":{},\"covariance_axes\":\"location_at_stress_center,scaled_stress_slope,log_sigma_if_estimated\",\"predict_parameters\":{},\"exposure_contract\":\"A(t)=integral_exp(-mu(stress(s)))ds;density=f_A(A(t))*A_prime(t)\",\"jump_density_convention\":\"right_continuous_stress_at_internal_segment_boundary\",\"unsupported\":[\"Bayesian_repeated_measures\",\"non_step_stress\",\"stress_history_effects_beyond_cumulative_exposure\"]}}",arr(&[raw_intercept,raw_slope]),-best.value*data.len() as f64,aic,data.len(),converged,best.iterations,best.gradient,best.line_failed,budget.get(),covariance_json,arr(&portable)))
}
#[allow(clippy::too_many_arguments)]
fn exposure_quantile(
    p: &Pattern,
    prob: f64,
    b0: f64,
    b1: f64,
    distribution: usize,
    sigma: f64,
    relation: usize,
    continuation: usize,
) -> Result<f64, String> {
    let target = match distribution {
        0 => (-(-prob).ln_1p()).ln(),
        1 => sigma * (-(-prob).ln_1p()).ln(),
        _ => sigma * normal_quantile(prob),
    };
    // Search log time so quantiles many orders below a segment duration retain
    // relative accuracy. The representable positive lifetime domain is explicit.
    let minimum = f64::from_bits(1).ln();
    let maximum = f64::MAX.ln();
    let at = |log_t: f64| -> Result<f64, String> {
        let t = if continuation == 0 {
            log_t.exp().min(p.duration)
        } else {
            log_t.exp()
        };
        Ok(exposure(p, t, b0, b1, relation, continuation, 0., 1.)?.0)
    };
    let mut hi = p.duration.ln();
    let mut brackets = 0;
    while at(hi)? < target {
        if continuation == 0 {
            return Ok(f64::NAN);
        }
        if hi >= maximum {
            return Ok(f64::INFINITY);
        }
        hi = (hi + 0.5 * hi.abs().max(1.)).min(maximum);
        brackets += 1;
        if brackets > 1024 {
            return Err("exposure quantile bracket work limit exceeded".into());
        }
    }
    let mut lo = minimum;
    if hi < lo || at(lo)? > target {
        return Ok(f64::NAN);
    }
    for _ in 0..120 {
        let mid = lo / 2. + hi / 2.;
        if at(mid)? < target {
            lo = mid
        } else {
            hi = mid
        }
    }
    Ok((lo / 2. + hi / 2.).exp())
}
fn damage_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 8 || args[0] != 1. || data[0].len() != 3 {
        return Err("damage prediction rows [time,pattern_id,quantile_probability] and saved version1 parameters".into());
    }
    let distribution = int(args[1], 0, 2)?;
    let relation = int(args[2], 0, 3)?;
    let continuation = int(args[3], 0, 2)?;
    let (b0, b1, sigma) = (args[4], args[5], args[6]);
    if sigma <= 0. || distribution == 0 && sigma != 1. {
        return Err("invalid cumulative-damage saved scale".into());
    }
    let patterns = parse_patterns(args, 7, relation)?;
    if data
        .len()
        .saturating_mul(patterns.iter().map(|p| p.segments.len()).max().unwrap_or(1))
        .saturating_mul(1200)
        > MAX_WORK
    {
        return Err("damage prediction exceeds bounded quantile work".into());
    }
    let mut ages = Vec::new();
    let mut logages = Vec::new();
    let mut survivals = Vec::new();
    let mut failure = Vec::new();
    let mut hazards = Vec::new();
    let mut densities = Vec::new();
    let mut quantiles = Vec::new();
    let mut loghazards = Vec::new();
    for row in data {
        if !(0. ..1.).contains(&row[2]) || row[2] == 0. {
            return Err("quantile probability must be in (0,1)".into());
        }
        let p = patterns
            .iter()
            .find(|p| p.id == row[1])
            .ok_or("unknown saved stress pattern")?;
        let (age, _, current) = exposure(p, row[0], b0, b1, relation, continuation, 0., 1.)?;
        let v = baseline(age, distribution, sigma.ln())?;
        let rate = -b0 - b1 * current;
        let loghazard = if age == f64::NEG_INFINITY {
            if distribution == 0 {
                rate
            } else if distribution == 1 {
                if sigma == 1. {
                    rate
                } else if sigma < 1. {
                    f64::NEG_INFINITY
                } else {
                    f64::INFINITY
                }
            } else {
                f64::NEG_INFINITY
            }
        } else {
            finite(v.ld + rate - v.ls, "damage log hazard")?
        };
        ages.push(age.exp());
        logages.push(age);
        survivals.push(v.ls.exp());
        failure.push(v.lf.exp());
        loghazards.push(loghazard);
        hazards.push(loghazard.exp());
        densities.push((v.ld + rate).exp());
        quantiles.push(exposure_quantile(
            p,
            row[2],
            b0,
            b1,
            distribution,
            sigma,
            relation,
            continuation,
        )?);
    }
    Ok(format!("{{\"version\":1,\"effective_ages\":{},\"log_effective_ages\":{},\"survivals\":{},\"failure_probabilities\":{},\"densities\":{},\"hazards\":{},\"log_hazards\":{},\"quantiles\":{},\"null_quantile_contract\":\"beyond_stopped_pattern_or_nonrepresentable_lifetime\",\"exposure_contract\":\"cumulative_exposure_not_Miner_fatigue_damage\"}}",arr(&ages),arr(&logages),arr(&survivals),arr(&failure),arr(&densities),arr(&hazards),arr(&loghazards),arr(&quantiles)))
}
