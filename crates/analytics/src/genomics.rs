//! Diploid marker QC, imputation, relatedness, simulation and fixed-ratio GBLUP.
use crate::{optimize::normal_cdf, resampling::Rng};
fn f(x: f64) -> String {
    if x.is_finite() {
        x.to_string()
    } else {
        "null".into()
    }
}
fn a(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(|&v| f(v)).collect::<Vec<_>>().join(",")
    )
}
fn m(x: &[Vec<f64>]) -> String {
    format!("[{}]", x.iter().map(|r| a(r)).collect::<Vec<_>>().join(","))
}
fn json(op: &str, body: String) -> String {
    format!("{{\"version\":1,\"method\":\"{op}\",{body}}}")
}
fn int(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn dosage(x: f64) -> bool {
    x == -1. || (0. ..=2.).contains(&x)
}
fn calls(data: &[Vec<f64>]) -> Result<(), String> {
    if data
        .iter()
        .flatten()
        .any(|v| ![-1., 0., 1., 2.].contains(v))
    {
        Err("diploid calls must be 0,1,2 or missing code -1".into())
    } else {
        Ok(())
    }
}
fn means(data: &[Vec<f64>], p: usize) -> Result<Vec<f64>, String> {
    let mut sums = vec![0.; p];
    let mut counts = vec![0; p];
    for r in data {
        for j in 0..p {
            if !dosage(r[j]) {
                return Err("invalid marker dosage".into());
            }
            if r[j] != -1. {
                sums[j] += r[j];
                counts[j] += 1;
            }
        }
    }
    for j in 0..p {
        if counts[j] == 0 {
            return Err("all-missing marker cannot be inferred".into());
        }
        sums[j] /= counts[j] as f64;
    }
    Ok(sums)
}
fn genotype_counts(data: &[Vec<f64>]) -> Vec<[usize; 3]> {
    let mut c = vec![[0; 3]; data[0].len()];
    for r in data {
        for (j, &v) in r.iter().enumerate() {
            if v >= 0. {
                c[j][v as usize] += 1;
            }
        }
    }
    c
}
fn statistics(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !args.is_empty() {
        return Err("marker_statistics takes no parameters".into());
    }
    calls(data)?;
    let counts = genotype_counts(data);
    let p = counts.len();
    let mut n = vec![0.; p];
    let mut freq = vec![f64::NAN; p];
    let mut observed = vec![f64::NAN; p];
    let mut diversity = vec![f64::NAN; p];
    let mut pic = vec![f64::NAN; p];
    let mut chi = vec![f64::NAN; p];
    let mut prob = vec![f64::NAN; p];
    for j in 0..p {
        n[j] = counts[j].iter().sum::<usize>() as f64;
        if n[j] == 0. {
            continue;
        }
        freq[j] = (counts[j][1] as f64 + 2. * counts[j][2] as f64) / (2. * n[j]);
        let q = 1. - freq[j];
        observed[j] = counts[j][1] as f64 / n[j];
        diversity[j] = 2. * freq[j] * q;
        pic[j] = diversity[j] - 2. * (freq[j] * q).powi(2);
        if freq[j] > 0. && q > 0. {
            let expected = [
                n[j] * q * q,
                n[j] * 2. * freq[j] * q,
                n[j] * freq[j] * freq[j],
            ];
            chi[j] = (0..3)
                .map(|k| (counts[j][k] as f64 - expected[k]).powi(2) / expected[k])
                .sum();
            prob[j] = 2. * normal_cdf(-chi[j].sqrt());
        }
    }
    let mut order: Vec<_> = (0..p).filter(|&j| prob[j].is_finite()).collect();
    order.sort_by(|&i, &j| prob[i].total_cmp(&prob[j]));
    let mut adjusted = vec![f64::NAN; p];
    let mut last = 1_f64;
    for rank in (0..order.len()).rev() {
        let j = order[rank];
        last = last.min(prob[j] * order.len() as f64 / (rank + 1) as f64);
        adjusted[j] = last;
    }
    Ok(json("marker_statistics",format!("\"called_count\":{},\"call_rate\":{},\"genotype_counts\":{},\"allele_frequency\":{},\"minor_allele_frequency\":{},\"observed_heterozygosity\":{},\"expected_heterozygosity\":{},\"polymorphism_information_content\":{},\"hwe_chi_square\":{},\"hwe_p_value\":{},\"hwe_bh_p_value\":{},\"hwe_test\":\"Pearson_chi_square_df1_asymptotic\",\"missing_code\":-1",a(&n),a(&n.iter().map(|&v|v/data.len() as f64).collect::<Vec<_>>()),m(&counts.iter().map(|c|c.iter().map(|&v|v as f64).collect()).collect::<Vec<_>>()),a(&freq),a(&freq.iter().map(|&v|if v.is_nan(){v}else{v.min(1.-v)}).collect::<Vec<_>>()),a(&observed),a(&diversity),a(&pic),a(&chi),a(&prob),a(&adjusted))))
}
fn mode(c: &[usize; 3]) -> f64 {
    let best = *c.iter().max().unwrap();
    c.iter().position(|&v| v == best).unwrap() as f64
}
fn pair_correlation(data: &[Vec<f64>], j: usize, k: usize) -> Option<f64> {
    let pairs: Vec<_> = data
        .iter()
        .filter(|r| r[j] >= 0. && r[k] >= 0.)
        .map(|r| (r[j], r[k]))
        .collect();
    if pairs.len() < 3 {
        return None;
    }
    let n = pairs.len() as f64;
    let x = pairs.iter().map(|v| v.0 / n).sum::<f64>();
    let y = pairs.iter().map(|v| v.1 / n).sum::<f64>();
    let xx = pairs.iter().map(|v| (v.0 - x).powi(2)).sum::<f64>();
    let yy = pairs.iter().map(|v| (v.1 - y).powi(2)).sum::<f64>();
    if xx == 0. || yy == 0. {
        None
    } else {
        Some(
            (pairs.iter().map(|v| (v.0 - x) * (v.1 - y)).sum::<f64>() / xx.sqrt() / yy.sqrt())
                .powi(2)
                .min(1.),
        )
    }
}
fn impute(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 5 {
        return Err("marker_impute expects [method0mean/1mode/2LDknn/3HWE/4observed/5uniform/6specified,nearest_markers=10,nearest_samples=5,seed=0,specified=0]".into());
    }
    calls(data)?;
    let method = int(*args.first().unwrap_or(&0.), 6)?;
    let km = int(*args.get(1).unwrap_or(&10.), 1024)?;
    let ks = int(*args.get(2).unwrap_or(&5.), 100000)?;
    let mut rng = Rng::new(*args.get(3).unwrap_or(&0.))?;
    let specified = int(*args.get(4).unwrap_or(&0.), 2)?;
    if km == 0 || ks == 0 {
        return Err("nearest counts must be positive".into());
    }
    let p = data[0].len();
    let center = means(data, p)?;
    let counts = genotype_counts(data);
    let modes: Vec<_> = counts.iter().map(mode).collect();
    let mut out = data.to_vec();
    let mut fallback = 0;
    let missing = data.iter().flatten().filter(|&&v| v < 0.).count();
    if method == 2
        && (p > 128
            || data.len() > 2000
            || data.len().saturating_mul(data.len()).saturating_mul(p) > 100_000_000
            || data.len().saturating_mul(p).saturating_mul(p) > 100_000_000)
    {
        return Err("LD-kNN work exceeds bounded complete-pair budget".into());
    }
    let mut neighbors = vec![Vec::new(); p];
    if method == 2 {
        for (j, near) in neighbors.iter_mut().enumerate() {
            let mut candidates: Vec<_> = (0..p)
                .filter(|&k| k != j)
                .filter_map(|k| pair_correlation(data, j, k).map(|r| (k, r)))
                .collect();
            candidates.sort_by(|a, b| b.1.total_cmp(&a.1).then(a.0.cmp(&b.0)));
            *near = candidates.into_iter().take(km).map(|v| v.0).collect();
        }
    }
    for (i, r) in out.iter_mut().enumerate() {
        for j in 0..p {
            if r[j] >= 0. {
                continue;
            }
            r[j] = match method {
                0 => center[j],
                1 => modes[j],
                2 => {
                    let mut candidates = Vec::new();
                    for (k, other) in data.iter().enumerate() {
                        if other[j] < 0. {
                            continue;
                        }
                        let mut distance = 0.;
                        let mut common = 0;
                        for &l in &neighbors[j] {
                            if data[i][l] >= 0. && other[l] >= 0. {
                                distance += (data[i][l] - other[l]).abs();
                                common += 1;
                            }
                        }
                        if common > 0 {
                            candidates.push((k, distance / common as f64));
                        }
                    }
                    candidates.sort_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
                    let mut votes = [0; 3];
                    for &(k, _) in candidates.iter().take(ks) {
                        votes[data[k][j] as usize] += 1;
                    }
                    if candidates.is_empty() {
                        fallback += 1;
                        modes[j]
                    } else {
                        mode(&votes)
                    }
                }
                6 => specified as f64,
                _ => {
                    let p = center[j] / 2.;
                    let total = counts[j].iter().sum::<usize>() as f64;
                    let probabilities = match method {
                        3 => [(1. - p) * (1. - p), 2. * p * (1. - p), p * p],
                        4 => [
                            counts[j][0] as f64 / total,
                            counts[j][1] as f64 / total,
                            counts[j][2] as f64 / total,
                        ],
                        _ => [1. / 3.; 3],
                    };
                    let uniform = rng.uniform();
                    let mut cumulative = 0.;
                    let mut selected = 2;
                    for (call, &probability) in probabilities.iter().enumerate() {
                        cumulative += probability;
                        if uniform < cumulative {
                            selected = call;
                            break;
                        }
                    }
                    selected as f64
                }
            };
        }
    }
    let saved = if method == 0 {
        a(&center)
    } else if method == 1 {
        a(&modes)
    } else {
        "null".into()
    };
    Ok(json("marker_impute",format!("\"rows\":{},\"imputed_cells\":{},\"method_code\":{},\"predict_parameters\":{},\"fallback_to_mode_count\":{},\"ld_metric\":\"squared_dosage_correlation_and_mean_taxicab_sample_distance\",\"tie_policy\":\"lowest_marker_sample_genotype_index\",\"training_scope\":\"all_supplied_rows_original_observed_calls_only\"",m(&out),missing,method,saved,fallback)))
}
fn impute_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != data[0].len()
        || args.iter().any(|v| !(0. ..=2.).contains(v))
        || data.iter().flatten().any(|&v| !dosage(v))
    {
        return Err("saved imputation values must match marker columns in [0,2]".into());
    }
    let out: Vec<_> = data
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .map(|(j, &v)| if v == -1. { args[j] } else { v })
                .collect()
        })
        .collect();
    Ok(json(
        "marker_impute_predict",
        format!("\"rows\":{}", m(&out)),
    ))
}
type CodedMarkers = (Vec<Vec<f64>>, Vec<f64>, f64);
fn coded(data: &[Vec<f64>], method: usize) -> Result<CodedMarkers, String> {
    let p = data[0].len();
    let center = means(data, p)?;
    let variance: Vec<_> = center.iter().map(|&v| v * (1. - v / 2.)).collect();
    if method == 2 {
        calls(data)?;
    }
    let denom = if method == 0 {
        variance.iter().sum()
    } else if method == 1 {
        variance.iter().filter(|&&v| v > 0.).count() as f64
    } else {
        variance.iter().map(|v| v * v).sum()
    };
    if denom <= 0. {
        return Err("relatedness requires at least one polymorphic marker".into());
    }
    let out = data
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .map(|(j, &v)| {
                    if v == -1. || variance[j] == 0. {
                        0.
                    } else if method == 0 {
                        v - center[j]
                    } else if method == 1 {
                        (v - center[j]) / variance[j].sqrt()
                    } else {
                        let p = center[j] / 2.;
                        match v as usize {
                            0 => -2. * p * p,
                            1 => 2. * p * (1. - p),
                            _ => -2. * (1. - p).powi(2),
                        }
                    }
                })
                .collect()
        })
        .collect();
    Ok((out, center, denom))
}
fn relatedness(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 1 {
        return Err("marker_relatedness expects [method0VanRaden1/1VanRaden2/2dominance]".into());
    }
    let method = int(*args.first().unwrap_or(&0.), 2)?;
    let n = data.len();
    let p = data[0].len();
    if n > 1000 || n.saturating_mul(n).saturating_mul(p) > 100_000_000 {
        return Err("pairwise relatedness exceeds work budget".into());
    }
    let (z, center, denom) = coded(data, method)?;
    let mut kernel = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v = z[i].iter().zip(&z[j]).map(|(&a, &b)| a * b / denom).sum();
            kernel[i][j] = v;
            kernel[j][i] = v;
        }
    }
    Ok(json("marker_relatedness",format!("\"matrix\":{},\"allele_frequency\":{},\"denominator\":{},\"method_code\":{},\"missing_policy\":\"zero_centered_additive_or_HWE_expected_dominance\"",m(&kernel),a(&center.iter().map(|v|v/2.).collect::<Vec<_>>()),f(denom),method)))
}
fn distance(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 2 {
        return Err("distance_matrix expects [method0Euclidean/1Jaccard/2BrayCurtis/3numericGower/4Hamming,optional_missing_code]".into());
    }
    let method = int(*args.first().unwrap_or(&0.), 4)?;
    let missing = args.get(1).copied();
    let n = data.len();
    let p = data[0].len();
    if n > 1000
        || n.saturating_mul(n).saturating_mul(p) > 100_000_000
        || data.iter().flatten().any(|v| v.abs() > 1e150)
    {
        return Err("distance matrix exceeds work or numerical bounds".into());
    }
    let valid = |v: f64| missing != Some(v);
    let mut ranges = vec![0.; p];
    for j in 0..p {
        let mut low = f64::INFINITY;
        let mut high = f64::NEG_INFINITY;
        for r in data {
            if valid(r[j]) {
                if method == 2 && r[j] < 0. {
                    return Err("Bray-Curtis requires nonnegative observations".into());
                }
                low = low.min(r[j]);
                high = high.max(r[j]);
            }
        }
        if high >= low {
            ranges[j] = high - low;
        }
    }
    let mut out = vec![vec![0.; n]; n];
    let mut common_count = out.clone();
    for i in 0..n {
        if !data[i].iter().any(|&v| valid(v)) {
            return Err("all-missing distance row".into());
        }
        common_count[i][i] = data[i].iter().filter(|&&v| valid(v)).count() as f64;
        for j in 0..i {
            let mut norm = 0_f64;
            let mut sum = 0.;
            let mut denominator = 0.;
            let mut common = 0;
            for (k, &range) in ranges.iter().enumerate() {
                let (x, y) = (data[i][k], data[j][k]);
                if !valid(x) || !valid(y) {
                    continue;
                }
                common += 1;
                match method {
                    0 => norm = norm.hypot(x - y),
                    1 => {
                        sum += f64::from((x != 0.) != (y != 0.));
                        denominator += f64::from(x != 0. || y != 0.);
                    }
                    2 => {
                        sum += (x - y).abs();
                        denominator += x + y;
                    }
                    3 => {
                        if range > 0. {
                            sum += (x - y).abs() / range;
                        }
                    }
                    _ => sum += f64::from(x != y),
                }
            }
            if common == 0 {
                return Err("distance pair has no common observed coordinates".into());
            }
            let value = match method {
                0 => norm,
                1 | 2 => {
                    if denominator > 0. {
                        sum / denominator
                    } else {
                        0.
                    }
                }
                _ => sum / common as f64,
            };
            out[i][j] = value;
            out[j][i] = value;
            common_count[i][j] = common as f64;
            common_count[j][i] = common as f64;
        }
    }
    Ok(json("distance_matrix",format!("\"matrix\":{},\"common_count\":{},\"method_code\":{},\"gower_policy\":\"numeric_ranges_constant_columns_contribute_zero\"",m(&out),m(&common_count),method)))
}
fn simulate(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() < 2 || args.len() > 3 || data.len() > 1024 {
        return Err(
            "marker_simulate requires allele-frequency scalar rows and [samples,seed,inbreeding=0]"
                .into(),
        );
    }
    let n = int(args[0], 100000)?;
    let p = data.len();
    let mut rng = Rng::new(args[1])?;
    let inbreeding = *args.get(2).unwrap_or(&0.);
    if n == 0
        || n.saturating_mul(p) > 1_000_000
        || !(0. ..=1.).contains(&inbreeding)
        || data.iter().any(|r| !(0. ..=1.).contains(&r[0]))
    {
        return Err("invalid genotype simulation parameters".into());
    }
    let mut out = vec![vec![0.; p]; n];
    for r in &mut out {
        for j in 0..p {
            let p = data[j][0];
            let q = 1. - p;
            let pr0 = q * q + inbreeding * p * q;
            let pr1 = 2. * p * q * (1. - inbreeding);
            let u = rng.uniform();
            r[j] = if u < pr0 {
                0.
            } else if u < pr0 + pr1 {
                1.
            } else {
                2.
            };
        }
    }
    Ok(json("marker_simulate",format!("\"rows\":{},\"seed\":{},\"inbreeding\":{},\"assumption\":\"independent_biallelic_diploid_loci\"",m(&out),f(args[1]),f(inbreeding))))
}
fn cross(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    calls(data)?;
    if args.len() != 4 || data.iter().flatten().any(|&v| v < 0.) {
        return Err(
            "marker_cross requires complete parental calls and [parent0,parent1,offspring,seed]"
                .into(),
        );
    }
    let p0 = int(args[0], data.len() - 1)?;
    let p1 = int(args[1], data.len() - 1)?;
    let n = int(args[2], 100000)?;
    let mut rng = Rng::new(args[3])?;
    let p = data[0].len();
    if n == 0 || n.saturating_mul(p) > 1_000_000 {
        return Err("cross output budget exceeded".into());
    }
    let mut out = vec![vec![0.; p]; n];
    for r in &mut out {
        for (j, v) in r.iter_mut().enumerate() {
            *v = f64::from(rng.uniform() < data[p0][j] / 2.)
                + f64::from(rng.uniform() < data[p1][j] / 2.);
        }
    }
    Ok(json("marker_cross",format!("\"rows\":{},\"seed\":{},\"parent_indices\":{},\"assumption\":\"one_generation_unlinked_loci_no_selection_or_recombination_map\"",m(&out),f(args[3]),a(&[p0 as f64,p1 as f64]))))
}
fn cholesky(matrix: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let p = matrix.len();
    let mut l = vec![vec![0.; p]; p];
    for i in 0..p {
        for j in 0..=i {
            let mut v = matrix[i][j];
            for (&left, &right) in l[i][..j].iter().zip(&l[j][..j]) {
                v -= left * right;
            }
            if i == j {
                if v <= 0. || !v.is_finite() {
                    return Err("genomic ridge information not positive definite".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    Ok(l)
}
fn solve(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let p = b.len();
    let mut x = b.to_vec();
    for i in 0..p {
        for j in 0..i {
            x[i] -= l[i][j] * x[j];
        }
        x[i] /= l[i][i];
    }
    for i in (0..p).rev() {
        for j in i + 1..p {
            x[i] -= l[j][i] * x[j];
        }
        x[i] /= l[i][i];
    }
    x
}
fn blup(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data[0].len() - 1;
    if p == 0 || p > 64 || data.len() > 10000 || args.is_empty() || args.len() > 3 {
        return Err("genomic_blup requires <=64markers plus outcome and [n_train,lambda>=1e-8=1,method0VanRaden1/1VanRaden2=0]".into());
    }
    let n = int(args[0], data.len())?;
    let lambda = *args.get(1).unwrap_or(&1.);
    let method = int(*args.get(2).unwrap_or(&0.), 1)?;
    if n < 2
        || !(1e-8..=1e12).contains(&lambda)
        || data[..n].iter().any(|r| r[p].abs() > 1e150)
        || data.iter().any(|r| r[..p].iter().any(|&v| !dosage(v)))
    {
        return Err("invalid genomic regression parameters".into());
    }
    let markers: Vec<_> = data[..n].iter().map(|r| r[..p].to_vec()).collect();
    let (z, center, c) = coded(&markers, method)?;
    let origin = data[0][p];
    let offset = data[..n]
        .iter()
        .map(|r| (r[p] - origin) / n as f64)
        .sum::<f64>();
    let y: Vec<_> = data[..n].iter().map(|r| (r[p] - origin) - offset).collect();
    let mut gram = vec![vec![0.; p]; p];
    let mut rhs = vec![0.; p];
    for (i, row) in z.iter().enumerate() {
        for j in 0..p {
            rhs[j] += row[j] * y[i];
            for (k, &value) in row.iter().enumerate() {
                gram[j][k] += row[j] * value;
            }
        }
    }
    for (j, row) in gram.iter_mut().enumerate() {
        row[j] += lambda * c;
    }
    let b = solve(&cholesky(&gram)?, &rhs);
    let effects: Vec<_> = (0..p)
        .map(|j| {
            if method == 0 {
                b[j]
            } else {
                let variance = center[j] * (1. - center[j] / 2.);
                if variance > 0. {
                    b[j] / variance.sqrt()
                } else {
                    0.
                }
            }
        })
        .collect();
    let mut params = vec![p as f64, origin, offset];
    params.extend_from_slice(&center);
    params.extend_from_slice(&effects);
    let mut predictions = Vec::new();
    let mut fitted = Vec::new();
    for (i, r) in data.iter().enumerate() {
        let value = origin
            + (offset
                + (0..p)
                    .map(|j| {
                        if r[j] == -1. {
                            0.
                        } else {
                            (r[j] - center[j]) * effects[j]
                        }
                    })
                    .sum::<f64>());
        if !value.is_finite() {
            return Err("genomic prediction overflow".into());
        }
        if i < n {
            fitted.push(value);
        } else {
            predictions.push(value);
        }
    }
    Ok(json("genomic_blup",format!("\"fitted\":{},\"predictions\":{},\"marker_effects\":{},\"feature_means\":{},\"predict_parameters\":{},\"variance_ratio\":{},\"relationship_method\":{},\"relationship_denominator\":{},\"variance_estimation\":\"fixed_user_supplied_ratio_not_REML\",\"training_rows\":{}",a(&fitted),a(&predictions),a(&effects),a(&center),a(&params),f(lambda),method,f(c),n)))
}
fn blup_predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 5 {
        return Err("invalid genomic prediction model".into());
    }
    let p = int(args[0], 64)?;
    if p == 0
        || data[0].len() != p
        || args.len() != 3 + 2 * p
        || args[3..3 + p].iter().any(|v| !(0. ..=2.).contains(v))
        || data.iter().flatten().any(|&v| !dosage(v))
    {
        return Err("genomic model does not match marker columns".into());
    }
    let out: Vec<_> = data
        .iter()
        .map(|r| {
            args[1]
                + (args[2]
                    + (0..p)
                        .map(|j| {
                            if r[j] == -1. {
                                0.
                            } else {
                                (r[j] - args[3 + j]) * args[3 + p + j]
                            }
                        })
                        .sum::<f64>())
        })
        .collect();
    if out.iter().any(|v| !v.is_finite()) {
        return Err("genomic prediction overflow".into());
    }
    Ok(json(
        "genomic_blup_predict",
        format!("\"predictions\":{}", a(&out)),
    ))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 100000
        || data[0].len() > 1024
        || data.len().saturating_mul(data[0].len()) > 1_000_000
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != data[0].len() || r.iter().any(|v| !v.is_finite()))
    {
        return Err("invalid or excessive genomics input".into());
    }
    match op {
        "marker_statistics" => statistics(data, args),
        "marker_impute" => impute(data, args),
        "marker_impute_predict" => impute_predict(data, args),
        "marker_relatedness" => relatedness(data, args),
        "distance_matrix" => distance(data, args),
        "marker_simulate" => simulate(data, args),
        "marker_cross" => cross(data, args),
        "genomic_blup" => blup(data, args),
        "genomic_blup_predict" => blup_predict(data, args),
        _ => Err("unknown genomics operation".into()),
    }
}
