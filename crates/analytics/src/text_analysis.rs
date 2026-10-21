//! Numeric Text Explorer kernels. Tokenization and term-ID dictionaries belong
//! to the caller; no language-specific parsing or sentiment lexicon is implied.
#![allow(clippy::needless_range_loop)]
type Matrix = Vec<Vec<f64>>;
fn array(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| array(r)).collect::<Vec<_>>().join(",")
    )
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded nonnegative integer".into())
    }
}
fn counts(data: &[Vec<f64>]) -> Result<(), String> {
    if data
        .iter()
        .flatten()
        .any(|&v| !(0. ..=1e9).contains(&v) || v.fract() != 0.)
    {
        Err("DTM counts must be integer0..1e9".into())
    } else {
        Ok(())
    }
}
fn summary(data: &[Vec<f64>]) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut total = vec![0.; data[0].len()];
    let mut df = total.clone();
    let mut lengths = vec![0.; data.len()];
    for (i, row) in data.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            total[j] += v;
            df[j] += f64::from(v > 0.);
            lengths[i] += v;
        }
    }
    (total, df, lengths)
}
fn dtm(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 2 || data[0].len() != 3 {
        return Err("text_dtm rows[document_id,term_id,count]; args[n_documents,n_terms]".into());
    }
    let n = integer(args[0], 4096)?;
    let p = integer(args[1], 4096)?;
    if n == 0 || p == 0 || n.saturating_mul(p) > 131072 {
        return Err("DTM output limited to131072cells and4096dimensions".into());
    }
    let mut out = vec![vec![0.; p]; n];
    for r in data {
        let i = integer(r[0], n - 1)?;
        let j = integer(r[1], p - 1)?;
        if !(0. ..=1e9).contains(&r[2]) || r[2].fract() != 0. {
            return Err("token counts must be integer0..1e9".into());
        }
        out[i][j] += r[2];
        if out[i][j] > 1e9 {
            return Err("aggregated cell count exceeds1e9".into());
        }
    }
    let (total, df, lengths) = summary(&out);
    Ok(format!("{{\"version\":1,\"method\":\"text_dtm\",\"counts\":{},\"term_totals\":{},\"document_frequency\":{},\"document_lengths\":{},\"tokenization\":\"caller_supplied_numeric_term_ids\"}}",matrix(&out),array(&total),array(&df),array(&lengths)))
}
fn weighted(data: &[Vec<f64>], mode: usize, selected: &[usize], idf: &[f64]) -> Matrix {
    data.iter()
        .map(|r| {
            selected
                .iter()
                .zip(idf)
                .map(|(&j, &weight)| match mode {
                    0 => f64::from(r[j] > 0.),
                    1 => r[j].min(2.),
                    2 => r[j],
                    3 => r[j].ln_1p() / std::f64::consts::LN_10,
                    _ => r[j] * weight,
                })
                .collect()
        })
        .collect()
}
fn weight(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    counts(data)?;
    if args.len() > 3 {
        return Err("text_weight args[mode0binary/1ternary/2frequency/3log10/4TFIDF=4,min_total=1,max_terms=all]".into());
    }
    let p = data[0].len();
    let mode = integer(*args.first().unwrap_or(&4.), 4)?;
    let minimum = *args.get(1).unwrap_or(&1.);
    let limit = integer(*args.get(2).unwrap_or(&(p as f64)), p)?;
    if minimum < 1. || minimum.fract() != 0. || minimum > 1e12 || limit == 0 {
        return Err("invalid minimum corpus count or maximum terms".into());
    }
    let (total, df, lengths) = summary(data);
    let mut selected = (0..p).filter(|&j| total[j] >= minimum).collect::<Vec<_>>();
    selected.sort_by(|&a, &b| total[b].total_cmp(&total[a]).then(a.cmp(&b)));
    selected.truncate(limit);
    if selected.is_empty() {
        return Err("no terms survive training corpus frequency filter".into());
    }
    let idf = selected
        .iter()
        .map(|&j| (data.len() as f64 / df[j]).log10())
        .collect::<Vec<_>>();
    let out = weighted(data, mode, &selected, &idf);
    let mut parameters = vec![
        1.,
        p as f64,
        selected.len() as f64,
        mode as f64,
        data.len() as f64,
    ];
    parameters.extend(selected.iter().map(|&j| j as f64));
    parameters.extend(&idf);
    Ok(format!("{{\"version\":1,\"method\":\"text_weight\",\"matrix\":{},\"selected_terms\":{},\"term_totals\":{},\"document_frequency\":{},\"document_lengths\":{},\"inverse_document_frequency\":{},\"weight_parameters\":{},\"idf_convention\":\"unsmoothed_log10_n_documents_over_training_document_frequency\"}}",matrix(&out),array(&selected.iter().map(|&j|j as f64).collect::<Vec<_>>()),array(&total),array(&df),array(&lengths),array(&idf),array(&parameters)))
}
fn weight_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    counts(data)?;
    if args.len() < 5 || args[0] != 1. {
        return Err("invalid saved weighting version".into());
    }
    let p = integer(args[1], 4096)?;
    let k = integer(args[2], p)?;
    let mode = integer(args[3], 4)?;
    let n = integer(args[4], 131072)?;
    if p != data[0].len() || k == 0 || n == 0 || args.len() != 5 + 2 * k {
        return Err("invalid saved weighting dimensions".into());
    }
    let selected = args[5..5 + k]
        .iter()
        .map(|&v| integer(v, p - 1))
        .collect::<Result<Vec<_>, _>>()?;
    let mut unique = selected.clone();
    unique.sort_unstable();
    unique.dedup();
    let idf = &args[5 + k..];
    if unique.len() != k
        || idf
            .iter()
            .any(|&v| v < 0. || v > (n as f64).log10() + 1e-12)
    {
        return Err("invalid saved dictionary indices or IDF".into());
    }
    Ok(format!("{{\"version\":1,\"method\":\"text_weight_predict\",\"matrix\":{},\"training_documents\":{},\"dictionary_refitted\":false}}",matrix(&weighted(data,mode,&selected,idf)),n))
}
struct Svd {
    a: Matrix,
    u: Matrix,
    v: Matrix,
    s: Vec<f64>,
    means: Vec<f64>,
    scales: Vec<f64>,
    divisor: f64,
    error: f64,
    rank: usize,
}
fn svd(data: &[Vec<f64>], k: usize, preprocess: usize) -> Result<Svd, String> {
    let n = data.len();
    let p = data[0].len();
    if n > 1024
        || p > 128
        || n.saturating_mul(p).saturating_mul(p) > 1_000_000
        || k == 0
        || k > 32
        || k > n.min(p)
        || data.iter().flatten().any(|&v| v < 0.)
    {
        return Err("LSA requires nonnegative supplied weighted DTM,n<=1024,p<=128,k<=32,min(n,p),n*p*p<=1e6".into());
    }
    if preprocess > 0 && n < 2 {
        return Err("centered LSA needs at least two documents".into());
    }
    let divisor = if preprocess == 0 { n } else { n - 1 } as f64;
    let mut means = vec![0.; p];
    let mut scales = vec![1.; p];
    let mut a = data.to_vec();
    for j in 0..p {
        if preprocess > 0 {
            let origin = data[0][j];
            means[j] = origin + data.iter().map(|r| (r[j] - origin) / n as f64).sum::<f64>();
        }
        if preprocess == 2 {
            let maximum = data
                .iter()
                .map(|r| (r[j] - means[j]).abs())
                .fold(0., f64::max);
            if maximum > 0. {
                scales[j] = maximum
                    * (data
                        .iter()
                        .map(|r| ((r[j] - means[j]) / maximum).powi(2))
                        .sum::<f64>()
                        / (n - 1) as f64)
                        .sqrt();
            }
        }
        for i in 0..n {
            a[i][j] = (data[i][j] - means[j]) / scales[j] / divisor;
        }
    }
    let global = a.iter().flatten().map(|v| v.abs()).fold(0., f64::max);
    if global == 0. || !global.is_finite() {
        return Err("transformed DTM has zero/nonfinite numerical rank".into());
    }
    let mut b = a
        .iter()
        .map(|r| r.iter().map(|v| v / global).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut v = vec![vec![0.; p]; p];
    for j in 0..p {
        v[j][j] = 1.;
    }
    let mut converged = false;
    for _ in 0..100 {
        let mut changed = false;
        for j in 0..p {
            for l in j + 1..p {
                let alpha = b.iter().map(|r| r[j] * r[j]).sum::<f64>();
                let beta = b.iter().map(|r| r[l] * r[l]).sum::<f64>();
                let gamma = b.iter().map(|r| r[j] * r[l]).sum::<f64>();
                if alpha == 0. || beta == 0. || gamma.abs() <= 1e-13 * alpha.sqrt() * beta.sqrt() {
                    continue;
                }
                changed = true;
                let tau = (beta - alpha) / (2. * gamma);
                let t = if tau == 0. {
                    1.
                } else {
                    tau.signum() / (tau.abs() + tau.hypot(1.))
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                for row in b.iter_mut().chain(v.iter_mut()) {
                    let x = row[j];
                    let y = row[l];
                    row[j] = c * x - s * y;
                    row[l] = s * x + c * y;
                }
            }
        }
        if !changed {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("one-sided Jacobi SVD did not converge in100sweeps".into());
    }
    let norms = (0..p)
        .map(|j| b.iter().map(|r| r[j] * r[j]).sum::<f64>().sqrt())
        .collect::<Vec<_>>();
    let mut order = (0..p).collect::<Vec<_>>();
    order.sort_by(|&j, &l| norms[l].total_cmp(&norms[j]).then(j.cmp(&l)));
    let largest = norms[order[0]];
    let rank = norms.iter().filter(|&&v| v > largest * 1e-10).count();
    if k > rank {
        return Err("requested singular vectors exceed relative1e-10 numerical rank".into());
    }
    let mut term = vec![vec![0.; k]; p];
    let mut u = vec![vec![0.; k]; n];
    let mut singular = vec![0.; k];
    for (c, &j) in order[..k].iter().enumerate() {
        let anchor = (0..p)
            .max_by(|&i, &l| v[i][j].abs().total_cmp(&v[l][j].abs()).then(l.cmp(&i)))
            .unwrap();
        let sign = if v[anchor][j] < 0. { -1. } else { 1. };
        singular[c] = norms[j] * global;
        if singular[c] == 0. || !singular[c].is_finite() {
            return Err("singular value outside representable units; rescale supplied DTM".into());
        }
        for i in 0..p {
            term[i][c] = sign * v[i][j];
        }
        for i in 0..n {
            u[i][c] = sign * b[i][j] / norms[j];
        }
    }
    let normalized_loadings = term
        .iter()
        .map(|r| {
            r.iter()
                .zip(&singular)
                .map(|(&v, &s)| v * (s / global))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let error = reconstruction_error(&a, &u, &normalized_loadings, global);
    Ok(Svd {
        a,
        u,
        v: term,
        s: singular,
        means,
        scales,
        divisor,
        error,
        rank,
    })
}
fn projection_parameters(f: &Svd, t: &[Vec<f64>]) -> Vec<f64> {
    let mut p = vec![1., t.len() as f64, t[0].len() as f64, f.divisor];
    p.extend(&f.means);
    p.extend(&f.scales);
    p.extend(t.iter().flatten());
    p
}
fn project(
    data: &[Vec<f64>],
    means: &[f64],
    scales: &[f64],
    divisor: f64,
    t: &[Vec<f64>],
) -> Matrix {
    data.iter()
        .map(|r| {
            (0..t[0].len())
                .map(|j| {
                    (0..t.len())
                        .map(|l| ((r[l] - means[l]) / scales[l] / divisor) * t[l][j])
                        .sum()
                })
                .collect()
        })
        .collect()
}
fn projection(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 4 || args[0] != 1. {
        return Err("invalid saved projection version".into());
    }
    let p = integer(args[1], 128)?;
    let k = integer(args[2], 32)?;
    if p != data[0].len()
        || k == 0
        || p == 0
        || data.len().saturating_mul(k) > 131072
        || args.len() != 4 + 2 * p + p * k
        || args[3] <= 0.
        || args[4 + p..4 + 2 * p].iter().any(|&v| v <= 0.)
        || data.iter().flatten().any(|&v| v < 0.)
    {
        return Err("invalid saved projection dimensions/scales or nonnegative DTM".into());
    }
    let t = args[4 + 2 * p..]
        .chunks(k)
        .map(|r| r.to_vec())
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"method\":\"text_projection\",\"scores\":{},\"preprocessing_refitted\":false}}",matrix(&project(data,&args[4..4+p],&args[4+p..4+2*p],args[3],&t))))
}
fn lsa(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 2 {
        return Err("text_lsa args[k,preprocess0uncentered/1centered/2centered_scaled=2]".into());
    }
    let k = integer(args[0], 32)?;
    let preprocess = integer(*args.get(1).unwrap_or(&2.), 2)?;
    let f = svd(data, k, preprocess)?;
    let scores =
        f.u.iter()
            .map(|r| r.iter().zip(&f.s).map(|(a, b)| a * b).collect())
            .collect::<Vec<Vec<f64>>>();
    let params = projection_parameters(&f, &f.v);
    Ok(format!("{{\"version\":1,\"method\":\"text_lsa\",\"singular_values\":{},\"document_vectors\":{},\"document_scores\":{},\"term_vectors\":{},\"means\":{},\"scales\":{},\"divisor\":{},\"numerical_rank\":{},\"relative_reconstruction_error\":{},\"projection_parameters\":{},\"normalization\":\"transformed_DTM_divided_by_n_or_n_minus_one_literal_JMP_help_convention\"}}",array(&f.s),matrix(&f.u),matrix(&scores),matrix(&f.v),array(&f.means),array(&f.scales),f.divisor,f.rank,f.error,array(&params)))
}
fn multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Matrix {
    a.iter()
        .map(|r| {
            (0..b[0].len())
                .map(|j| r.iter().enumerate().map(|(l, &v)| v * b[l][j]).sum())
                .collect()
        })
        .collect()
}
fn reconstruction_error(
    a: &[Vec<f64>],
    left: &[Vec<f64>],
    normalized_right: &[Vec<f64>],
    scale: f64,
) -> f64 {
    let mut residual = 0.;
    let mut total = 0.;
    for (i, row) in a.iter().enumerate() {
        for (j, &a) in row.iter().enumerate() {
            let target = a / scale;
            let fitted = left[i]
                .iter()
                .zip(&normalized_right[j])
                .map(|(a, b)| a * b)
                .sum::<f64>();
            residual += (target - fitted).powi(2);
            total += target * target;
        }
    }
    (residual / total).sqrt()
}
fn varimax_objective(l: &[Vec<f64>]) -> f64 {
    let p = l.len() as f64;
    (0..l[0].len())
        .map(|j| {
            l.iter().map(|r| r[j].powi(4)).sum::<f64>()
                - l.iter().map(|r| r[j] * r[j]).sum::<f64>().powi(2) / p
        })
        .sum()
}
fn topic(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 5 {
        return Err(
            "text_topic_svd args[k,preprocess=2,max_sweeps=200,tol=1e-8,kaiser0/1=1]".into(),
        );
    }
    let k = integer(args[0], 16)?;
    let preprocess = integer(*args.get(1).unwrap_or(&2.), 2)?;
    let maximum = integer(*args.get(2).unwrap_or(&200.), 2000)?;
    let tolerance = *args.get(3).unwrap_or(&1e-8);
    let kaiser = integer(*args.get(4).unwrap_or(&1.), 1)? == 1;
    let f = svd(data, k, preprocess)?;
    let p = f.v.len();
    if maximum == 0
        || !(1e-12..=1e-3).contains(&tolerance)
        || p.saturating_mul(k)
            .saturating_mul(k)
            .saturating_mul(maximum)
            > 50_000_000
    {
        return Err("invalid varimax controls or50Mwork cap".into());
    }
    let loadings =
        f.v.iter()
            .map(|r| r.iter().zip(&f.s).map(|(&v, &s)| v * s).collect::<Vec<_>>())
            .collect::<Vec<_>>();
    // A common positive scaling leaves the maximizer unchanged, while making
    // fourth-power objectives and stopping criteria independent of DTM units.
    let loading_scale = loadings
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(0., f64::max);
    let mut l = loadings
        .iter()
        .map(|r| r.iter().map(|v| v / loading_scale).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let norms = l
        .iter()
        .map(|r| r.iter().map(|x| x * x).sum::<f64>().sqrt())
        .collect::<Vec<_>>();
    if kaiser {
        for (i, r) in l.iter_mut().enumerate() {
            if norms[i] > 0. {
                for x in r {
                    *x /= norms[i];
                }
            }
        }
    }
    let mut rotation = vec![vec![0.; k]; k];
    for i in 0..k {
        rotation[i][i] = 1.;
    }
    let mut converged = k == 1;
    let mut iterations = 0;
    let mut objective = varimax_objective(&l);
    for sweep in 0..maximum {
        if converged {
            break;
        }
        let before = objective;
        for j in 0..k {
            for m in j + 1..k {
                let (mut su, mut sv, mut aa, mut bb) = (0., 0., 0., 0.);
                for r in &l {
                    let u = r[j] * r[j] - r[m] * r[m];
                    let v = 2. * r[j] * r[m];
                    su += u;
                    sv += v;
                    aa += u * u - v * v;
                    bb += 2. * u * v;
                }
                let angle = 0.25
                    * (bb - 2. * su * sv / p as f64).atan2(aa - (su * su - sv * sv) / p as f64);
                let (c, s) = (angle.cos(), angle.sin());
                for r in l.iter_mut().chain(rotation.iter_mut()) {
                    let x = r[j];
                    let y = r[m];
                    r[j] = c * x + s * y;
                    r[m] = -s * x + c * y;
                }
            }
        }
        iterations = sweep + 1;
        objective = varimax_objective(&l);
        if objective < before - 1e-10 * (1. + before.abs()) {
            return Err("varimax objective decreased beyond floating tolerance".into());
        }
        converged = (objective - before).abs() <= tolerance * (1. + before.abs());
    }
    let mut rotated = multiply(&loadings, &rotation);
    let mut scores = multiply(&f.u, &rotation);
    for j in 0..k {
        let anchor = (0..p)
            .max_by(|&i, &l| {
                rotated[i][j]
                    .abs()
                    .total_cmp(&rotated[l][j].abs())
                    .then(l.cmp(&i))
            })
            .unwrap();
        if rotated[anchor][j] < 0. {
            for r in &mut rotated {
                r[j] = -r[j];
            }
            for r in &mut scores {
                r[j] = -r[j];
            }
            for r in &mut rotation {
                r[j] = -r[j];
            }
        }
    }
    let projection =
        f.v.iter()
            .map(|r| {
                (0..k)
                    .map(|j| (0..k).map(|m| r[m] / f.s[m] * rotation[m][j]).sum())
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
    let energy = (0..k)
        .map(|j| rotated.iter().map(|r| r[j] * r[j]).sum())
        .collect::<Vec<f64>>();
    if energy.contains(&0.) {
        return Err(
            "positive topic energy underflows physical binary64 units; rescale supplied DTM".into(),
        );
    }
    let normalized_loadings = rotated
        .iter()
        .map(|r| r.iter().map(|v| v / loading_scale).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let error = reconstruction_error(&f.a, &scores, &normalized_loadings, loading_scale);
    Ok(format!("{{\"version\":1,\"method\":\"text_topic_svd\",\"topic_loadings\":{},\"document_scores\":{},\"rotation\":{},\"topic_energy\":{},\"varimax_objective\":{},\"varimax_loading_scale\":{},\"kaiser_normalized\":{},\"rotation_sweeps\":{},\"rotation_converged\":{},\"relative_reconstruction_error\":{},\"projection_parameters\":{},\"definition\":\"loadings_VS_rotated_by_R_document_vectors_UR\"}}",matrix(&rotated),matrix(&scores),matrix(&rotation),array(&energy),objective,loading_scale,kaiser,iterations,converged,error,array(&projection_parameters(&f,&projection))))
}
struct Random(u64);
impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn index(&mut self, n: usize) -> usize {
        let cap = u64::MAX - u64::MAX % n as u64;
        loop {
            let x = self.next();
            if x < cap {
                return (x % n as u64) as usize;
            }
        }
    }
}
fn expectation(data: &[Vec<f64>], mix: &[f64], prob: &[Vec<f64>]) -> (f64, Matrix, Vec<f64>) {
    let mut ll = 0.;
    let mut out = Vec::with_capacity(data.len());
    let mut document_ll = Vec::with_capacity(data.len());
    for row in data {
        let scores = prob
            .iter()
            .zip(mix)
            .map(|(theta, &m)| {
                m.ln()
                    + theta
                        .iter()
                        .zip(row)
                        .map(|(&p, &x)| if x > 0. { p.ln() } else { (-p).ln_1p() })
                        .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let high = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let normal = high + scores.iter().map(|v| (v - high).exp()).sum::<f64>().ln();
        out.push(scores.iter().map(|v| (v - normal).exp()).collect());
        ll += normal;
        document_ll.push(normal);
    }
    (ll, out, document_ll)
}
fn class_parameters(mix: &[f64], prob: &[Vec<f64>], floor: f64) -> Vec<f64> {
    let mut args = vec![1., mix.len() as f64, prob[0].len() as f64, floor];
    args.extend(mix);
    args.extend(prob.iter().flatten());
    args
}
fn class_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    counts(data)?;
    if args.len() < 4 || args[0] != 1. {
        return Err("invalid saved latent-class version".into());
    }
    let k = integer(args[1], 16)?;
    let p = integer(args[2], 128)?;
    let floor = args[3];
    if k == 0
        || p != data[0].len()
        || !(1e-12..=0.1).contains(&floor)
        || data.len().saturating_mul(k) > 131072
        || args.len() != 4 + k + k * p
    {
        return Err("invalid saved latent-class dimensions/floor".into());
    }
    let mix = &args[4..4 + k];
    let prob = args[4 + k..]
        .chunks(p)
        .map(|r| r.to_vec())
        .collect::<Vec<_>>();
    if mix.iter().any(|&v| v <= 0. || v > 1.)
        || (mix.iter().sum::<f64>() - 1.).abs() > 1e-10
        || prob.iter().flatten().any(|&v| v < floor || v > 1. - floor)
    {
        return Err("invalid saved mixture/term probabilities".into());
    }
    let (ll, responsibilities, document_ll) = expectation(data, mix, &prob);
    let predicted = responsibilities
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .max_by(|(i, a), (j, b)| a.total_cmp(b).then(j.cmp(i)))
                .unwrap()
                .0 as f64
        })
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"method\":\"text_latent_class_predict\",\"probabilities\":{},\"predicted_class\":{},\"document_log_probability\":{},\"log_likelihood\":{}}}",matrix(&responsibilities),array(&predicted),array(&document_ll),ll))
}
fn class(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    counts(data)?;
    if args.is_empty() || args.len() > 5 {
        return Err(
            "text_latent_class args[k,seed=0,max_iter=200,tol=1e-8,probability_floor=1e-9]".into(),
        );
    }
    let n = data.len();
    let p = data[0].len();
    let k = integer(args[0], 16)?;
    let seed = *args.get(1).unwrap_or(&0.);
    let maximum = integer(*args.get(2).unwrap_or(&200.), 2000)?;
    let tolerance = *args.get(3).unwrap_or(&1e-8);
    let floor = *args.get(4).unwrap_or(&1e-9);
    if !(2..=2048).contains(&n)
        || p > 128
        || k == 0
        || k > n
        || maximum == 0
        || !(0. ..=9007199254740991.).contains(&seed)
        || seed.fract() != 0.
        || !(1e-12..=1e-3).contains(&tolerance)
        || !(1e-12..=0.1).contains(&floor)
        || n.saturating_mul(p)
            .saturating_mul(k)
            .saturating_mul(maximum)
            > 100_000_000
    {
        return Err("invalid bounded Bernoulli mixture controls or100Mwork cap".into());
    }
    let mut rng = Random(seed as u64);
    let mut indices = (0..n).collect::<Vec<_>>();
    for i in (1..n).rev() {
        let j = rng.index(i + 1);
        indices.swap(i, j);
    }
    let initial = indices[..k].to_vec();
    let mut probability = initial
        .iter()
        .map(|&i| {
            data[i]
                .iter()
                .map(|&v| if v > 0. { 0.75 } else { 0.25 })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut mix = vec![1. / k as f64; k];
    let (mut ll, mut responsibility, _) = expectation(data, &mix, &probability);
    let mut trace = vec![ll];
    let mut optimizer_converged = false;
    let mut iterations = 0;
    for iteration in 0..maximum {
        let mass = (0..k)
            .map(|j| responsibility.iter().map(|r| r[j]).sum::<f64>())
            .collect::<Vec<_>>();
        if mass.iter().any(|&v| v <= 1e-12) {
            return Err("latent class component has numerically zero posterior mass; change seed or reduce k".into());
        }
        let next_mix = mass.iter().map(|v| v / n as f64).collect::<Vec<_>>();
        let mut next = vec![vec![0.; p]; k];
        for j in 0..k {
            for t in 0..p {
                next[j][t] = (data
                    .iter()
                    .zip(&responsibility)
                    .filter(|(r, _)| r[t] > 0.)
                    .map(|(_, r)| r[j])
                    .sum::<f64>()
                    / mass[j])
                    .clamp(floor, 1. - floor);
            }
        }
        let (next_ll, next_responsibility, _) = expectation(data, &next_mix, &next);
        if next_ll < ll - 1e-10 * (1. + ll.abs()) {
            return Err("Bernoulli EM likelihood decreased beyond roundoff".into());
        }
        let change = next
            .iter()
            .flatten()
            .zip(probability.iter().flatten())
            .map(|(a, b)| (a - b).abs())
            .chain(next_mix.iter().zip(&mix).map(|(a, b)| (a - b).abs()))
            .fold(0., f64::max);
        optimizer_converged =
            (next_ll - ll).abs() <= tolerance * (1. + ll.abs()) && change <= tolerance.sqrt();
        ll = next_ll;
        mix = next_mix;
        probability = next;
        responsibility = next_responsibility;
        trace.push(ll);
        iterations = iteration + 1;
        if optimizer_converged {
            break;
        }
    }
    let mut order = (0..k).collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        mix[b].total_cmp(&mix[a]).then_with(|| {
            for t in 0..p {
                let c = probability[b][t].total_cmp(&probability[a][t]);
                if !c.is_eq() {
                    return c;
                }
            }
            a.cmp(&b)
        })
    });
    let mix = order.iter().map(|&j| mix[j]).collect::<Vec<_>>();
    let probability = order
        .iter()
        .map(|&j| probability[j].clone())
        .collect::<Vec<_>>();
    let (ll, responsibility, _) = expectation(data, &mix, &probability);
    let distinct = (0..k).all(|j| {
        (j + 1..k).all(|l| {
            probability[j]
                .iter()
                .zip(&probability[l])
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max)
                > 1e-6
        })
    });
    let parameters = k * p + k - 1;
    let dimension_bound = parameters as f64 <= 2_f64.powi(p as i32) - 1.;
    let converged = optimizer_converged && distinct && dimension_bound;
    let aic = if converged {
        (2. * parameters as f64 - 2. * ll).to_string()
    } else {
        "null".into()
    };
    let bic = if converged {
        (parameters as f64 * (n as f64).ln() - 2. * ll).to_string()
    } else {
        "null".into()
    };
    let predicted = responsibility
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .max_by(|(i, a), (j, b)| a.total_cmp(b).then(j.cmp(i)))
                .unwrap()
                .0 as f64
        })
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"method\":\"text_latent_class\",\"mixture_probabilities\":{},\"term_probabilities\":{},\"probabilities\":{},\"predicted_class\":{},\"log_likelihood\":{},\"log_likelihood_trace\":{},\"parameters\":{},\"aic\":{},\"bic\":{},\"optimizer_converged\":{},\"components_distinct\":{},\"necessary_dimension_bound\":{},\"identification_verified\":false,\"converged\":{},\"iterations\":{},\"initial_document_indices\":{},\"predict_parameters\":{},\"estimator\":\"bounded_conditional_independence_Bernoulli_mixture_local_EM\"}}",array(&mix),matrix(&probability),matrix(&responsibility),array(&predicted),ll,array(&trace),parameters,aic,bic,optimizer_converged,distinct,dimension_bound,converged,iterations,array(&initial.iter().map(|&i|i as f64).collect::<Vec<_>>()),array(&class_parameters(&mix,&probability,floor))))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 131072
        || data[0].len() > 4096
        || data.len().saturating_mul(data[0].len()) > 131072
        || data
            .iter()
            .any(|r| r.len() != data[0].len() || r.iter().any(|v| !v.is_finite() || v.abs() > 1e12))
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
    {
        return Err(
            "text kernels require bounded finite rectangular numeric input,<=131072cells,abs<=1e12"
                .into(),
        );
    }
    let result = match op {
        "text_dtm" => dtm(data, args),
        "text_weight" => weight(data, args),
        "text_weight_predict" => weight_predict(data, args),
        "text_lsa" => lsa(data, args),
        "text_projection" => projection(data, args),
        "text_topic_svd" => topic(data, args),
        "text_latent_class" => class(data, args),
        "text_latent_class_predict" => class_predict(data, args),
        _ => Err("unknown text-analysis operation".into()),
    }?;
    if result
        .split(['[', ']', ',', ':', '{', '}'])
        .any(|s| matches!(s, "NaN" | "inf" | "-inf"))
    {
        Err(
            "derived text model result outside finite binary64 units; rescale supplied matrix"
                .into(),
        )
    } else {
        Ok(result)
    }
}
