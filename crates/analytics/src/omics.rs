//! Count normalization, quantile maps and independent-locus diploid admixture.
use crate::resampling::Rng;
const WORK: usize = 200_000_000;
fn a(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn m(x: &[Vec<f64>]) -> String {
    format!("[{}]", x.iter().map(|r| a(r)).collect::<Vec<_>>().join(","))
}
fn j(op: &str, body: String) -> String {
    format!("{{\"version\":1,\"method\":\"{op}\",{body}}}")
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn normalized(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() {
        return Err("assay_normalize requires a method code".into());
    }
    let method = integer(args[0], 6)?;
    let p = data[0].len();
    if (method == 0 || method == 4) && args.len() > 2
        || [1, 5, 6].contains(&method) && args.len() != 1
        || [2, 3].contains(&method) && args.len() != p + 1
    {
        return Err("normalization parameters do not match selected method".into());
    }
    let ddof = if method == 0 {
        integer(*args.get(1).unwrap_or(&1.), 1)?
    } else {
        0
    };
    let pseudo = if method == 4 {
        *args.get(1).unwrap_or(&0.)
    } else {
        0.
    };
    if pseudo < 0. || [2, 3].contains(&method) && args[1..].iter().any(|&v| v <= 0.) {
        return Err("pseudocount must be nonnegative and lengths positive".into());
    }
    let mut out = Vec::new();
    for row in data {
        let values: Vec<f64> = if method == 0 {
            if p <= ddof || row.iter().any(|v| v.abs() > 1e150) {
                return Err("row standardization needs variation and bounded values".into());
            }
            let origin = row[0];
            let scale = row.iter().map(|&v| (v - origin).abs()).fold(0., f64::max);
            if scale == 0. {
                return Err("row variance is zero".into());
            }
            // Standardize dimensionless differences so centering and the SD
            // remain representable even when physical units are subnormal.
            let offset = row
                .iter()
                .map(|&v| ((v - origin) / scale) / p as f64)
                .sum::<f64>();
            let centered: Vec<_> = row.iter().map(|&v| (v - origin) / scale - offset).collect();
            let norm = centered.iter().fold(0_f64, |sum, &v| sum.hypot(v));
            let multiplier = ((p - ddof) as f64).sqrt();
            centered.iter().map(|v| (v / norm) * multiplier).collect()
        } else {
            if row.iter().any(|&v| v < 0.) {
                return Err("count normalization requires nonnegative observations".into());
            }
            if method == 4 {
                let mut logs = Vec::new();
                for &x in row {
                    let high = x.max(pseudo);
                    let low = x.min(pseudo);
                    if high == 0. {
                        return Err(
                            "CLR requires positive counts after explicit pseudocount".into()
                        );
                    }
                    logs.push(high.ln() + (low / high).ln_1p());
                }
                let origin = logs[0];
                let mean_offset = logs.iter().map(|&v| (v - origin) / p as f64).sum::<f64>();
                logs.iter().map(|&v| (v - origin) - mean_offset).collect()
            } else {
                let maximum = row.iter().copied().fold(0., f64::max);
                if maximum == 0. {
                    return Err("count row total is zero".into());
                }
                let scaled_total = row.iter().map(|v| v / maximum).sum::<f64>();
                if method == 2 {
                    let log_total = maximum.ln() + scaled_total.ln();
                    row.iter()
                        .zip(&args[1..])
                        .map(|(&x, &length)| {
                            if x == 0. {
                                0.
                            } else {
                                (x.ln() - log_total - length.ln() + 1e9_f64.ln()).exp()
                            }
                        })
                        .collect()
                } else if method == 3 {
                    let log_rates: Vec<_> = row
                        .iter()
                        .zip(&args[1..])
                        .map(|(&x, &length)| {
                            if x == 0. {
                                f64::NEG_INFINITY
                            } else {
                                x.ln() - length.ln()
                            }
                        })
                        .collect();
                    let high = log_rates.iter().copied().fold(f64::NEG_INFINITY, f64::max);
                    let rates: Vec<_> = log_rates.iter().map(|&v| (v - high).exp()).collect();
                    let total = rates.iter().sum::<f64>();
                    rates.iter().map(|v| 1e6 * v / total).collect()
                } else {
                    let scale = if method == 1 {
                        1e6
                    } else if method == 5 {
                        100.
                    } else {
                        1.
                    };
                    row.iter()
                        .map(|v| scale * (v / maximum) / scaled_total)
                        .collect()
                }
            }
        };
        if values.iter().any(|v: &f64| !v.is_finite()) {
            return Err("normalized output is outside finite binary64 range".into());
        }
        out.push(values);
    }
    Ok(j("assay_normalize",format!("\"rows\":{},\"method_code\":{},\"orientation\":\"samples_rows_features_columns\",\"length_units\":\"base_pairs_for_RPKM_and_TPM\",\"pseudocount\":{}",m(&out),method,pseudo)))
}
fn quantile_map(row: &[f64], reference: &[f64]) -> Vec<f64> {
    let mut order: Vec<_> = (0..row.len()).collect();
    order.sort_by(|&i, &j| row[i].total_cmp(&row[j]).then(i.cmp(&j)));
    let mut out = vec![0.; row.len()];
    let mut start = 0;
    while start < order.len() {
        let mut end = start + 1;
        while end < order.len() && row[order[end]] == row[order[start]] {
            end += 1;
        }
        let origin = reference[start];
        let mean = origin
            + reference[start..end]
                .iter()
                .map(|&v| (v - origin) / (end - start) as f64)
                .sum::<f64>();
        for &i in &order[start..end] {
            out[i] = mean;
        }
        start = end;
    }
    out
}
fn quantiles(data: &[Vec<f64>], args: &[f64], predict: bool) -> Result<String, String> {
    let p = data[0].len();
    if data.iter().flatten().any(|v| v.abs() > 1e150) {
        return Err("quantile values exceed numerical bounds".into());
    }
    let reference = if predict {
        if args.len() != p
            || args.windows(2).any(|v| v[1] < v[0])
            || args.iter().any(|v| v.abs() > 1e150)
        {
            return Err(
                "quantile reference must be finite sorted values matching row width".into(),
            );
        }
        args.to_vec()
    } else {
        if !args.is_empty() {
            return Err("quantile_normalize takes no parameters".into());
        }
        let sorted: Vec<_> = data
            .iter()
            .map(|r| {
                let mut r = r.clone();
                r.sort_by(f64::total_cmp);
                r
            })
            .collect();
        (0..p)
            .map(|i| {
                let origin = sorted[0][i];
                origin
                    + sorted
                        .iter()
                        .map(|r| (r[i] - origin) / data.len() as f64)
                        .sum::<f64>()
            })
            .collect()
    };
    let out: Vec<_> = data.iter().map(|r| quantile_map(r, &reference)).collect();
    Ok(j(if predict{"quantile_normalize_predict"}else{"quantile_normalize"},format!("\"rows\":{},\"reference_quantiles\":{},\"tie_policy\":\"mean_reference_quantiles_over_equal_value_ranks\",\"fit_scope\":\"{}\"",m(&out),a(&reference),if predict{"saved_reference_only"}else{"all_supplied_samples"})))
}
fn validate_genotypes(data: &[Vec<f64>]) -> Result<Vec<f64>, String> {
    if data
        .iter()
        .flatten()
        .any(|v| ![-1., 0., 1., 2.].contains(v))
        || data.iter().any(|r| r.iter().all(|&v| v == -1.))
    {
        return Err(
            "admixture requires diploid0/1/2,-1missing and >=1observed call per row".into(),
        );
    }
    let mut frequency = vec![0.; data[0].len()];
    for k in 0..frequency.len() {
        let mut called = 0;
        for r in data {
            if r[k] >= 0. {
                frequency[k] += r[k];
                called += 1;
            }
        }
        if called == 0 {
            return Err("all-missing marker is unidentified".into());
        }
        frequency[k] /= 2. * called as f64;
    }
    Ok(frequency)
}
fn likelihood(data: &[Vec<f64>], q: &[Vec<f64>], freq: &[Vec<f64>]) -> f64 {
    let mut value = 0.;
    for (i, row) in data.iter().enumerate() {
        for (k, &g) in row.iter().enumerate() {
            if g < 0. {
                continue;
            }
            let p = q[i]
                .iter()
                .enumerate()
                .map(|(c, &weight)| weight * freq[c][k])
                .sum::<f64>();
            value += g * p.ln() + (2. - g) * (-p).ln_1p() + if g == 1. { 2_f64.ln() } else { 0. };
        }
    }
    value
}
struct Admixture {
    q: Vec<Vec<f64>>,
    frequency: Vec<Vec<f64>>,
    history: Vec<f64>,
    iterations: usize,
    converged: bool,
}
fn start(
    data: &[Vec<f64>],
    global: &[f64],
    k: usize,
    rng: &mut Rng,
    floor: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let n = data.len();
    let p = global.len();
    let mut chosen = vec![rng.index(n)];
    while chosen.len() < k {
        let mut best = (0, f64::NEG_INFINITY);
        for i in 0..n {
            if chosen.contains(&i) {
                continue;
            }
            let distance = chosen
                .iter()
                .map(|&c| {
                    (0..p)
                        .map(|j| {
                            let a = if data[i][j] < 0. {
                                global[j]
                            } else {
                                data[i][j] / 2.
                            };
                            let b = if data[c][j] < 0. {
                                global[j]
                            } else {
                                data[c][j] / 2.
                            };
                            (a - b).abs()
                        })
                        .sum::<f64>()
                })
                .fold(f64::INFINITY, f64::min);
            if distance > best.1 {
                best = (i, distance);
            }
        }
        chosen.push(best.0);
    }
    let freq: Vec<_> = chosen
        .iter()
        .map(|&i| {
            (0..p)
                .map(|j| {
                    let observed = if data[i][j] < 0. {
                        global[j]
                    } else {
                        data[i][j] / 2.
                    };
                    (0.95 * observed + 0.05 * global[j]).clamp(floor, 1. - floor)
                })
                .collect()
        })
        .collect();
    let q = (0..n)
        .map(|_| {
            let mut row: Vec<_> = (0..k).map(|_| 0.1 + rng.uniform()).collect();
            let sum = row.iter().sum::<f64>();
            for x in &mut row {
                *x /= sum;
            }
            row
        })
        .collect();
    (q, freq)
}
fn fit_once(
    data: &[Vec<f64>],
    global: &[f64],
    k: usize,
    limit: usize,
    tolerance: f64,
    floor: f64,
    rng: &mut Rng,
) -> Result<Admixture, String> {
    let n = data.len();
    let p = global.len();
    if k == 1 {
        let q = vec![vec![1.]; n];
        let freq = vec![global.iter().map(|&v| v.clamp(floor, 1. - floor)).collect()];
        let ll = likelihood(data, &q, &freq);
        return Ok(Admixture {
            q,
            frequency: freq,
            history: vec![ll],
            iterations: 1,
            converged: true,
        });
    }
    let (mut q, mut frequency) = start(data, global, k, rng, floor);
    let mut history = vec![likelihood(data, &q, &frequency)];
    let mut converged = false;
    let mut iterations = 0;
    for iteration in 1..=limit {
        let mut q_count = vec![vec![0.; k]; n];
        let mut alt = vec![vec![0.; p]; k];
        let mut total = vec![vec![0.; p]; k];
        for (i, row) in data.iter().enumerate() {
            for (j, &g) in row.iter().enumerate() {
                if g < 0. {
                    continue;
                }
                let prob = (0..k).map(|c| q[i][c] * frequency[c][j]).sum::<f64>();
                for c in 0..k {
                    let a = g * q[i][c] * frequency[c][j] / prob;
                    let b = (2. - g) * q[i][c] * (1. - frequency[c][j]) / (1. - prob);
                    q_count[i][c] += a + b;
                    alt[c][j] += a;
                    total[c][j] += a + b;
                }
            }
        }
        for row in &mut q_count {
            let sum = row.iter().sum::<f64>();
            for v in row {
                *v /= sum;
            }
        }
        let mut next_frequency = frequency.clone();
        for c in 0..k {
            for j in 0..p {
                if total[c][j] > 0. {
                    next_frequency[c][j] = (alt[c][j] / total[c][j]).clamp(floor, 1. - floor);
                }
            }
        }
        let next = likelihood(data, &q_count, &next_frequency);
        let previous = *history.last().unwrap();
        if !next.is_finite() || next < previous - 1e-9 * (1. + previous.abs()) {
            return Err("admixture EM lost finite monotone likelihood".into());
        }
        q = q_count;
        frequency = next_frequency;
        history.push(next);
        iterations = iteration;
        if (next - previous).abs() <= tolerance * (1. + previous.abs()) {
            converged = true;
            break;
        }
    }
    Ok(Admixture {
        q,
        frequency,
        history,
        iterations,
        converged,
    })
}
fn admixture(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 6 {
        return Err("marker_admixture expects [populations,maxiter=1000,tolerance=1e-8,seed=0,n_init=3,frequency_floor=1e-8]".into());
    }
    let k = integer(args[0], 16)?;
    let limit = integer(*args.get(1).unwrap_or(&1000.), 10000)?;
    let tol = *args.get(2).unwrap_or(&1e-8);
    let seed = *args.get(3).unwrap_or(&0.);
    let n_init = integer(*args.get(4).unwrap_or(&3.), 8)?;
    let floor = *args.get(5).unwrap_or(&1e-8);
    let mut rng = Rng::new(seed)?;
    let n = data.len();
    let p = data[0].len();
    let work = n
        .saturating_mul(p)
        .saturating_mul(k)
        .saturating_mul(if k == 1 { 1 } else { limit + k + 1 })
        .saturating_mul(n_init);
    if k == 0
        || k > n
        || limit == 0
        || n_init == 0
        || !(1e-12..=1e-3).contains(&tol)
        || !(1e-12..=0.01).contains(&floor)
        || 5 + k + k * p > 16384
        || work > WORK
    {
        return Err("invalid or excessive admixture parameters".into());
    }
    let global = validate_genotypes(data)?;
    let mut best: Option<Admixture> = None;
    for _ in 0..n_init {
        let candidate = fit_once(data, &global, k, limit, tol, floor, &mut rng)?;
        if best
            .as_ref()
            .is_none_or(|b| candidate.history.last().unwrap() > b.history.last().unwrap())
        {
            best = Some(candidate);
        }
    }
    let mut fit = best.unwrap();
    let mut order: Vec<_> = (0..k).collect();
    order.sort_by(|&a, &b| {
        for j in 0..p {
            let compare = fit.frequency[a][j].total_cmp(&fit.frequency[b][j]);
            if compare != std::cmp::Ordering::Equal {
                return compare;
            }
        }
        a.cmp(&b)
    });
    fit.frequency = order.iter().map(|&i| fit.frequency[i].clone()).collect();
    fit.q = fit
        .q
        .iter()
        .map(|r| order.iter().map(|&i| r[i]).collect())
        .collect();
    let mut initial = vec![0.; k];
    for row in &fit.q {
        for c in 0..k {
            initial[c] += row[c] / n as f64;
        }
    }
    let mut params = vec![k as f64, p as f64, limit as f64, tol, floor];
    params.extend_from_slice(&initial);
    params.extend(fit.frequency.iter().flatten());
    Ok(j("marker_admixture",format!("\"ancestry\":{},\"allele_frequencies\":{},\"log_likelihood\":{},\"log_likelihood_history\":{},\"predict_parameters\":{},\"iterations\":{},\"converged\":{},\"n_init\":{},\"seed\":{},\"frequency_floor\":{},\"convergence_criterion\":\"relative_observed_log_likelihood\",\"identifiability\":\"population_labels_canonicalized_but_parameters_can_be_nonunique\",\"parameter_uncertainty\":null,\"assumptions\":\"unrelated_diploid_independent_loci_HWE_binomial_admixture\"",m(&fit.q),m(&fit.frequency),fit.history.last().unwrap(),a(&fit.history),a(&params),fit.iterations,fit.converged,n_init,seed,floor)))
}
fn admixture_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 7 {
        return Err("invalid saved admixture model".into());
    }
    let k = integer(args[0], 16)?;
    let p = integer(args[1], 1024)?;
    let limit = integer(args[2], 10000)?;
    let tol = args[3];
    let floor = args[4];
    if k == 0
        || p != data[0].len()
        || limit == 0
        || args.len() != 5 + k + k * p
        || !(1e-12..=1e-3).contains(&tol)
        || !(1e-12..=0.01).contains(&floor)
        || data
            .len()
            .saturating_mul(p)
            .saturating_mul(k)
            .saturating_mul(limit)
            > WORK
    {
        return Err("saved admixture dimensions/limits invalid".into());
    }
    if data
        .iter()
        .flatten()
        .any(|v| ![-1., 0., 1., 2.].contains(v))
        || data.iter().any(|r| r.iter().all(|&v| v < 0.))
    {
        return Err("query genotypes must be observed at >=1marker".into());
    }
    let initial = &args[5..5 + k];
    let sum = initial.iter().sum::<f64>();
    if initial.iter().any(|&v| v <= 0.) || !sum.is_finite() {
        return Err("invalid ancestry starting weights".into());
    }
    let frequency: Vec<Vec<f64>> = args[5 + k..].chunks_exact(p).map(<[f64]>::to_vec).collect();
    if frequency
        .iter()
        .flatten()
        .any(|&v| v < floor || v > 1. - floor)
    {
        return Err("saved allele frequency outside fitted floor".into());
    }
    let mut ancestry = Vec::new();
    let mut row_converged = Vec::new();
    let mut row_iterations = Vec::new();
    let mut row_likelihood = Vec::new();
    for row in data {
        let mut q: Vec<_> = initial.iter().map(|v| v / sum).collect();
        let one = vec![row.clone()];
        let mut ll = likelihood(&one, &[q.clone()], &frequency);
        let mut converged = k == 1;
        let mut iterations = 0;
        if k > 1 {
            for iteration in 1..=limit {
                let mut count = vec![0.; k];
                for (j, &g) in row.iter().enumerate() {
                    if g < 0. {
                        continue;
                    }
                    let prob = (0..k).map(|c| q[c] * frequency[c][j]).sum::<f64>();
                    for c in 0..k {
                        count[c] += g * q[c] * frequency[c][j] / prob
                            + (2. - g) * q[c] * (1. - frequency[c][j]) / (1. - prob);
                    }
                }
                let total = count.iter().sum::<f64>();
                for v in &mut count {
                    *v /= total;
                }
                let next = likelihood(&one, &[count.clone()], &frequency);
                if !next.is_finite() {
                    return Err("admixture query likelihood invalid".into());
                }
                let difference = (next - ll).abs();
                q = count;
                ll = next;
                iterations = iteration;
                if difference <= tol * (1. + ll.abs()) {
                    converged = true;
                    break;
                }
            }
        }
        ancestry.push(q);
        row_converged.push(converged);
        row_iterations.push(iterations as f64);
        row_likelihood.push(ll);
    }
    Ok(j("admixture_predict",format!("\"ancestry\":{},\"sample_log_likelihood\":{},\"iterations\":{},\"row_converged\":[{}],\"converged\":{},\"fitting_scope\":\"ancestry_only_saved_allele_frequencies_fixed\"",m(&ancestry),a(&row_likelihood),a(&row_iterations),row_converged.iter().map(bool::to_string).collect::<Vec<_>>().join(","),row_converged.iter().all(|&v|v))))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 10000
        || data[0].len() > 1024
        || data.len().saturating_mul(data[0].len()) > 1_000_000
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != data[0].len() || r.iter().any(|v| !v.is_finite()))
    {
        return Err("omics requires bounded finite rectangular inputs".into());
    }
    match op {
        "assay_normalize" => normalized(data, args),
        "quantile_normalize" => quantiles(data, args, false),
        "quantile_normalize_predict" => quantiles(data, args, true),
        "marker_admixture" => admixture(data, args),
        "admixture_predict" => admixture_predict(data, args),
        _ => Err("unknown omics operation".into()),
    }
}
