//! Binary market baskets and complete-case item reliability.
use std::collections::{BTreeMap, BTreeSet};
const WORK: usize = 100_000_000;
fn a(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn nullable(x: Option<f64>) -> String {
    x.filter(|v| v.is_finite())
        .map_or("null".into(), |v| v.to_string())
}
fn optional_array(x: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|&v| nullable(v)).collect::<Vec<_>>().join(",")
    )
}
fn result(op: &str, body: String) -> String {
    format!("{{\"version\":1,\"method\":\"{op}\",{body}}}")
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn indices(mask: u32) -> String {
    format!(
        "[{}]",
        (0..32)
            .filter(|&i| mask & (1 << i) != 0)
            .map(|i| i.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn baskets(data: &[Vec<f64>]) -> Result<Vec<u32>, String> {
    if data[0].len() > 32 || data.iter().flatten().any(|v| *v != 0. && *v != 1.) {
        return Err("market baskets require <=32 binary item columns".into());
    }
    Ok(data
        .iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .fold(0, |mask, (i, &v)| mask | if v == 1. { 1 << i } else { 0 })
        })
        .collect())
}
fn association(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 6 {
        return Err("association_rules expects [support=.1,confidence=.4,lift=1.2,max_antecedents=3,max_rule_size=4,max_rules=10000]".into());
    }
    let support = *args.first().unwrap_or(&0.1);
    let confidence = *args.get(1).unwrap_or(&0.4);
    let lift = *args.get(2).unwrap_or(&1.2);
    let antecedents = integer(*args.get(3).unwrap_or(&3.), 5)?;
    let size = integer(*args.get(4).unwrap_or(&4.), 6)?;
    let limit = integer(*args.get(5).unwrap_or(&10000.), 100000)?;
    if !(1e-12..=1.).contains(&support)
        || !(0. ..=1.).contains(&confidence)
        || lift < 0.
        || antecedents == 0
        || size < 2
        || limit == 0
    {
        return Err("invalid rule thresholds or size limits".into());
    }
    let rows = baskets(data)?;
    let n = data.len();
    let p = data[0].len();
    let mut all = BTreeMap::new();
    let mut level = Vec::new();
    let mut work = 0_usize;
    let mut generated = 0_usize;
    for i in 0..p {
        let mask = 1 << i;
        let count = rows.iter().filter(|&&r| r & mask == mask).count();
        if count as f64 / n as f64 >= support {
            all.insert(mask, count);
            level.push(mask);
        }
    }
    for width in 2..=size.min(p) {
        let previous: BTreeSet<_> = level.iter().copied().collect();
        let mut candidates = BTreeSet::new();
        for (i, &left) in level.iter().enumerate() {
            for &right in &level[i + 1..] {
                generated += 1;
                if generated > 1_000_000 {
                    return Err("association candidate join work limit reached".into());
                }
                let mask = left | right;
                if mask.count_ones() as usize != width {
                    continue;
                }
                if (0..p)
                    .filter(|&j| mask & (1 << j) != 0)
                    .all(|j| previous.contains(&(mask & !(1 << j))))
                {
                    candidates.insert(mask);
                }
            }
        }
        level.clear();
        for mask in candidates {
            work = work.saturating_add(n);
            if work > WORK {
                return Err("association transaction counting work limit reached".into());
            }
            let count = rows.iter().filter(|&&r| r & mask == mask).count();
            if count as f64 / n as f64 >= support {
                all.insert(mask, count);
                level.push(mask);
                if all.len() > 100000 {
                    return Err("frequent itemset output limit reached".into());
                }
            }
        }
        if level.is_empty() {
            break;
        }
    }
    let mut rules = Vec::new();
    let mut considered = 0;
    for (&joint, &count) in &all {
        if joint.count_ones() < 2 {
            continue;
        }
        let mut lhs = (joint - 1) & joint;
        while lhs > 0 {
            let rhs = joint ^ lhs;
            if lhs.count_ones() as usize <= antecedents && rhs > 0 {
                considered += 1;
                if considered > 1_000_000 {
                    return Err("association rule generation work limit reached".into());
                }
                let nx = all[&lhs];
                let ny = all[&rhs];
                let px = nx as f64 / n as f64;
                let py = ny as f64 / n as f64;
                let pxy = count as f64 / n as f64;
                let conf = count as f64 / nx as f64;
                let ratio = conf / py;
                if conf >= confidence && ratio >= lift {
                    let conviction = if count == nx {
                        None
                    } else {
                        Some((1. - py) / (1. - conf))
                    };
                    rules.push((lhs,rhs,format!("{{\"antecedent\":{},\"consequent\":{},\"count\":{},\"support\":{},\"confidence\":{},\"lift\":{},\"leverage\":{},\"conviction\":{},\"conviction_infinite\":{}}}",indices(lhs),indices(rhs),count,pxy,conf,ratio,pxy-px*py,nullable(conviction),count==nx&&ny<n)));
                    if rules.len() > limit {
                        return Err("association rule output limit reached; raise support/confidence/lift or max_rules".into());
                    }
                }
            }
            lhs = (lhs - 1) & joint;
        }
    }
    rules.sort_by_key(|r| (r.0, r.1));
    let sets: Vec<_> = all
        .iter()
        .map(|(&mask, &count)| {
            format!(
                "{{\"items\":{},\"count\":{},\"support\":{}}}",
                indices(mask),
                count,
                count as f64 / n as f64
            )
        })
        .collect();
    Ok(result("association_rules",format!("\"transaction_count\":{},\"item_count\":{},\"frequent_itemsets\":[{}],\"rules\":[{}],\"rule_count\":{},\"algorithm\":\"Apriori_complete_within_explicit_size_and_work_limits\",\"empty_transactions\":\"included_in_support_denominator\",\"truncated\":false",n,p,sets.join(","),rules.iter().map(|r|r.2.as_str()).collect::<Vec<_>>().join(","),rules.len())))
}
fn affinity(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !args.is_empty() {
        return Err("basket_affinity takes no parameters".into());
    }
    let rows = baskets(data)?;
    let n = rows.len();
    let p = data[0].len();
    if n.saturating_mul(p).saturating_mul(p) > WORK {
        return Err("affinity pair work limit reached".into());
    }
    let counts: Vec<_> = (0..p)
        .map(|i| rows.iter().filter(|&&r| r & (1 << i) != 0).count())
        .collect();
    let mut pairs = Vec::new();
    for i in 0..p {
        for j in i + 1..p {
            let mask = (1 << i) | (1 << j);
            let joint = rows.iter().filter(|&&r| r & mask == mask).count();
            let a = counts[i];
            let b = counts[j];
            let union = a + b - joint;
            let support = joint as f64 / n as f64;
            let px = a as f64 / n as f64;
            let py = b as f64 / n as f64;
            pairs.push(format!("{{\"item0\":{},\"item1\":{},\"joint_count\":{},\"support\":{},\"confidence0_to1\":{},\"confidence1_to0\":{},\"lift\":{},\"jaccard\":{},\"leverage\":{}}}",i,j,joint,support,nullable((a>0).then(||joint as f64/a as f64)),nullable((b>0).then(||joint as f64/b as f64)),nullable((a>0&&b>0).then(||support/(px*py))),nullable((union>0).then(||joint as f64/union as f64)),support-px*py));
        }
    }
    Ok(result("basket_affinity",format!("\"transaction_count\":{},\"item_counts\":{},\"pairs\":[{}],\"interpretation\":\"descriptive_pair_affinity_without_multiple_testing_or_causal_inference\"",n,a(&counts.iter().map(|&v|v as f64).collect::<Vec<_>>()),pairs.join(","))))
}
fn reliability(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 1 {
        return Err("item_reliability expects [standardized0/1=0]".into());
    }
    let standardized = integer(*args.first().unwrap_or(&0.), 1)? == 1;
    let n = data.len();
    let p = data[0].len();
    if n < 2
        || !(2..=64).contains(&p)
        || n.saturating_mul(p).saturating_mul(p) > WORK
        || data.iter().flatten().any(|v| v.abs() > 1e150)
    {
        return Err("item reliability needs >=2rows,2..64bounded columns and <=100Mwork".into());
    }
    let origins = &data[0];
    let scale = data
        .iter()
        .flat_map(|r| r.iter().zip(origins).map(|(&v, &o)| (v - o).abs()))
        .fold(0., f64::max)
        .max(f64::MIN_POSITIVE);
    let scales: Vec<_> = (0..p)
        .map(|j| {
            if standardized {
                data.iter()
                    .map(|r| (r[j] - origins[j]).abs())
                    .fold(0., f64::max)
                    .max(f64::MIN_POSITIVE)
            } else {
                scale
            }
        })
        .collect();
    let means: Vec<_> = (0..p)
        .map(|j| {
            data.iter()
                .map(|r| ((r[j] - origins[j]) / scales[j]) / n as f64)
                .sum::<f64>()
        })
        .collect();
    let mut centered: Vec<Vec<f64>> = data
        .iter()
        .map(|r| {
            (0..p)
                .map(|j| (r[j] - origins[j]) / scales[j] - means[j])
                .collect()
        })
        .collect();
    let variance: Vec<_> = (0..p)
        .map(|j| {
            centered
                .iter()
                .map(|r| r[j] * r[j] / (n - 1) as f64)
                .sum::<f64>()
        })
        .collect();
    if standardized {
        if variance.iter().any(|&v| v <= 0.) {
            return Err("standardized alpha needs positive variance for every item".into());
        }
        for row in &mut centered {
            for j in 0..p {
                row[j] /= variance[j].sqrt();
            }
        }
    }
    let diagonal: Vec<_> = (0..p)
        .map(|j| {
            centered
                .iter()
                .map(|r| r[j] * r[j] / (n - 1) as f64)
                .sum::<f64>()
        })
        .collect();
    let totals: Vec<_> = centered.iter().map(|r| r.iter().sum::<f64>()).collect();
    let total_variance = totals.iter().map(|v| v * v / (n - 1) as f64).sum::<f64>();
    let trace = diagonal.iter().sum::<f64>();
    let alpha =
        (total_variance > 0.).then(|| p as f64 / (p - 1) as f64 * (1. - trace / total_variance));
    let mut deleted = Vec::new();
    let mut correlations = Vec::new();
    for j in 0..p {
        let rest: Vec<_> = centered
            .iter()
            .map(|r| {
                r.iter()
                    .enumerate()
                    .filter(|&(i, _)| i != j)
                    .map(|(_, &v)| v)
                    .sum::<f64>()
            })
            .collect();
        let rest_var = rest.iter().map(|v| v * v / (n - 1) as f64).sum::<f64>();
        let rest_trace = diagonal
            .iter()
            .enumerate()
            .filter(|&(i, _)| i != j)
            .map(|(_, &v)| v)
            .sum::<f64>();
        deleted.push(
            (p > 2 && rest_var > 0.)
                .then(|| (p - 1) as f64 / (p - 2) as f64 * (1. - rest_trace / rest_var)),
        );
        let cov = centered
            .iter()
            .zip(&rest)
            .map(|(r, &v)| r[j] * v / (n - 1) as f64)
            .sum::<f64>();
        correlations.push(
            (diagonal[j] > 0. && rest_var > 0.)
                .then(|| (cov / diagonal[j].sqrt() / rest_var.sqrt()).clamp(-1., 1.)),
        );
    }
    Ok(result("item_reliability",format!("\"alpha\":{},\"alpha_if_deleted\":{},\"corrected_item_total_correlations\":{},\"standardized\":{},\"sample_count\":{},\"item_count\":{},\"missing_policy\":\"complete_finite_rows_only\",\"interpretation\":\"internal_consistency_not_unidimensionality_or_validity\"",nullable(alpha),optional_array(&deleted),optional_array(&correlations),standardized,n,p)))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 100000
        || data[0].len() > 1024
        || data.len().saturating_mul(data[0].len()) > 1_000_000
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != data[0].len() || r.iter().any(|v| !v.is_finite()))
    {
        return Err("market analytics requires bounded finite rectangular inputs".into());
    }
    match op {
        "association_rules" => association(data, args),
        "basket_affinity" => affinity(data, args),
        "item_reliability" => reliability(data, args),
        _ => Err("unknown market operation".into()),
    }
}
