//! Numeric covariance, clustering, distances, and prefix-trained classifiers.
//!
//! Classification/regression rows are `[features..., label]`; `n_train` selects
//! a training prefix, and predictions cover the suffix only. Query labels are
//! ignored. Classification labels must be exactly represented signed integers.
//! LDA uses pooled unbiased within-class covariance; QDA uses unbiased per-class
//! covariance. Regularization adds an absolute diagonal ridge in original units.
//! Gaussian naive Bayes uses population class variances with an absolute floor.
//! KNN standardization is fitted on the training prefix alone; clustering and
//! covariance use original feature units without automatic standardization.

use std::collections::{BTreeMap, VecDeque};

const MAX_CELLS: usize = 1_000_000;
const MAX_FEATURES: usize = 256;
const MAX_WORK: usize = 100_000_000;
const MAX_INTEGER: f64 = 9_007_199_254_740_991.;

/// Runs a bounded numeric operation and returns version-1 JSON.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = validate(data)?;
    if args.iter().any(|x| !x.is_finite()) {
        return Err("parameters must be finite".into());
    }
    match operation {
        "covariance_matrix" | "correlation_matrix" => covariance(operation, data, args),
        "hierarchical_cluster" => hierarchical(data, args),
        "mahalanobis" => mahalanobis(data, args),
        "discriminant_analysis" => discriminant(data, args),
        "gaussian_naive_bayes" => naive_bayes(data, args),
        "knn_predict" => knn(data, args),
        "dbscan" => dbscan(data, p, args),
        _ => Err(format!("unknown multivariate operation: {operation}")),
    }
}

fn validate(data: &[Vec<f64>]) -> Result<usize, String> {
    let p = data.first().map_or(0, Vec::len);
    if p == 0 || p > MAX_FEATURES + 1 || data.len() > MAX_CELLS / p {
        return Err("data must have 1 to 257 columns and at most 1000000 cells".into());
    }
    if data
        .iter()
        .any(|row| row.len() != p || row.iter().any(|x| !x.is_finite()))
    {
        return Err("data must be finite and rectangular".into());
    }
    Ok(p)
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn integer(x: f64, name: &str) -> Result<usize, String> {
    if !(0. ..=MAX_INTEGER).contains(&x) || x.fract() != 0. {
        return Err(format!("{name} must be an exact nonnegative integer"));
    }
    Ok(x as usize)
}
fn mean(values: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut m = 0.;
    for (i, x) in values.enumerate() {
        m = finite(
            m + finite(x - m, "mean difference")? / (i + 1) as f64,
            "mean",
        )?;
    }
    Ok(m)
}
fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0_f64, f64::hypot)
}
fn numbers(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(v: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}

struct Center {
    origins: Vec<f64>,
    offsets: Vec<f64>,
    rows: Vec<Vec<f64>>,
}
impl Center {
    fn new(data: &[Vec<f64>], p: usize) -> Result<Self, String> {
        let origins = data[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|r| r[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let rows = data
            .iter()
            .map(|r| {
                (0..p)
                    .map(|j| finite((r[j] - origins[j]) - offsets[j], "centered feature"))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            origins,
            offsets,
            rows,
        })
    }
    fn mean(&self) -> Result<Vec<f64>, String> {
        self.origins
            .iter()
            .zip(&self.offsets)
            .map(|(o, m)| finite(o + m, "feature mean"))
            .collect()
    }
    fn difference(&self, row: &[f64]) -> Result<Vec<f64>, String> {
        (0..self.origins.len())
            .map(|j| {
                finite(
                    (row[j] - self.origins[j]) - self.offsets[j],
                    "feature difference",
                )
            })
            .collect()
    }
    fn covariance(&self, denominator: usize) -> Result<Vec<Vec<f64>>, String> {
        let p = self.origins.len();
        if denominator == 0 {
            return Err("covariance requires positive degrees of freedom".into());
        }
        if self.rows.len().saturating_mul(p).saturating_mul(p) > MAX_WORK {
            return Err("covariance exceeds numerical work limit".into());
        }
        let mut cov = vec![vec![0.; p]; p];
        for i in 0..p {
            for j in 0..=i {
                let c = finite(
                    self.rows.iter().map(|r| r[i] * r[j]).sum::<f64>() / denominator as f64,
                    "covariance",
                )?;
                cov[i][j] = c;
                cov[j][i] = c;
            }
        }
        Ok(cov)
    }
}

fn covariance(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() > MAX_FEATURES {
        return Err("covariance supports at most 256 features".into());
    }
    let sample = if operation == "correlation_matrix" {
        if !args.is_empty() {
            return Err("correlation_matrix takes no parameters".into());
        }
        true
    } else {
        if args.len() > 1 || args.first().is_some_and(|x| *x != 0. && *x != 1.) {
            return Err("covariance requires [sample=1], 0 population or 1 sample".into());
        }
        args.first().copied().unwrap_or(1.) == 1.
    };
    if sample && data.len() < 2 {
        return Err("sample covariance requires at least two rows".into());
    }
    let center = Center::new(data, data[0].len())?;
    let cov = center.covariance(data.len() - usize::from(sample))?;
    if operation == "correlation_matrix" {
        let p = cov.len();
        let mut correlation = vec![vec![0.; p]; p];
        for i in 0..p {
            if cov[i][i] <= 0. {
                return Err("correlation is undefined for constant features".into());
            }
            for j in 0..p {
                if cov[j][j] <= 0. {
                    return Err("correlation is undefined for constant features".into());
                }
                correlation[i][j] = finite(
                    cov[i][j] / cov[i][i].sqrt() / cov[j][j].sqrt(),
                    "correlation",
                )?
                .clamp(-1., 1.);
            }
        }
        Ok(format!("{{\"version\":1,\"algorithm\":\"correlation_matrix\",\"n_samples\":{},\"mean\":{},\"correlation\":{}}}",data.len(),numbers(&center.mean()?),matrix(&correlation)))
    } else {
        Ok(format!("{{\"version\":1,\"algorithm\":\"covariance_matrix\",\"n_samples\":{},\"sample\":{sample},\"mean\":{},\"covariance\":{}}}",data.len(),numbers(&center.mean()?),matrix(&cov)))
    }
}

fn distance(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(
        norm(a.iter().zip(b).map(|(x, y)| x - y)),
        "Euclidean distance",
    )
}

fn hierarchical(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 2 {
        return Err("hierarchical_cluster requires [clusters, linkage=2]".into());
    }
    let n = data.len();
    let clusters = integer(args[0], "clusters")?;
    let method = args.get(1).copied().unwrap_or(2.);
    if clusters == 0 || clusters > n || method.fract() != 0. || !(0. ..=3.).contains(&method) {
        return Err(
            "clusters must be 1..n and linkage 0 single,1 complete,2 average,3 Ward".into(),
        );
    }
    if n > 500 || n.saturating_mul(n).saturating_mul(n) > MAX_WORK {
        return Err("hierarchical clustering exceeds 500 rows or numerical work limit".into());
    }
    let size = 2 * n - 1;
    let mut distances = vec![vec![0.; size]; size];
    for i in 0..n {
        for j in 0..i {
            let d = distance(&data[i], &data[j])?;
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }
    let mut active = vec![false; size];
    active[..n].fill(true);
    let mut counts = vec![1_usize; size];
    let mut members = (0..size)
        .map(|i| if i < n { vec![i] } else { vec![] })
        .collect::<Vec<_>>();
    let mut merges = Vec::with_capacity(n.saturating_sub(1));
    let mut labels = vec![0.; n];
    fn assign_labels(active: &[bool], members: &[Vec<usize>], labels: &mut [f64]) {
        let mut groups = active
            .iter()
            .enumerate()
            .filter(|(_, a)| **a)
            .map(|(i, _)| i)
            .collect::<Vec<_>>();
        groups.sort_by_key(|i| members[*i].iter().min().copied().unwrap_or(usize::MAX));
        for (label, i) in groups.iter().enumerate() {
            for j in &members[*i] {
                labels[*j] = label as f64;
            }
        }
    }
    if clusters == n {
        assign_labels(&active, &members, &mut labels);
    }
    for step in 0..n.saturating_sub(1) {
        let mut best = None;
        for i in 0..n + step {
            if !active[i] {
                continue;
            }
            for j in i + 1..n + step {
                if !active[j] {
                    continue;
                }
                let candidate = (distances[i][j], i, j);
                if best.is_none_or(|old: (f64, usize, usize)| candidate < old) {
                    best = Some(candidate);
                }
            }
        }
        let (d, a, b) = best.ok_or("no active cluster pair")?;
        let u = n + step;
        counts[u] = counts[a] + counts[b];
        members[u] = members[a].iter().chain(&members[b]).copied().collect();
        for v in 0..u {
            if !active[v] || v == a || v == b {
                continue;
            }
            let da = distances[a][v];
            let db = distances[b][v];
            let updated = match method as usize {
                0 => da.min(db),
                1 => da.max(db),
                2 => (counts[a] as f64 * da + counts[b] as f64 * db) / counts[u] as f64,
                _ => {
                    // Scale before squaring in the Ward Lance-Williams update.
                    let scale = da.max(db).max(d);
                    if scale == 0. {
                        0.
                    } else {
                        let total = (counts[v] + counts[u]) as f64;
                        let squared = ((counts[v] + counts[a]) as f64 * (da / scale).powi(2)
                            + (counts[v] + counts[b]) as f64 * (db / scale).powi(2)
                            - counts[v] as f64 * (d / scale).powi(2))
                            / total;
                        scale * squared.max(0.).sqrt()
                    }
                }
            };
            distances[u][v] = finite(updated, "linkage distance")?;
            distances[v][u] = updated;
        }
        active[a] = false;
        active[b] = false;
        active[u] = true;
        merges.push(vec![a as f64, b as f64, d, counts[u] as f64]);
        if n - step - 1 == clusters {
            assign_labels(&active, &members, &mut labels);
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"hierarchical_cluster\",\"linkage\":{method},\"n_clusters\":{clusters},\"labels\":{},\"merges\":{},\"tie_break\":\"smallest_cluster_id_pair\"}}",numbers(&labels),matrix(&merges)))
}

struct Cholesky {
    lower: Vec<Vec<f64>>,
    scales: Vec<f64>,
    log_determinant: f64,
}
impl Cholesky {
    #[allow(clippy::needless_range_loop)]
    fn new(covariance: &[Vec<f64>]) -> Result<Self, String> {
        let p = covariance.len();
        let mut scales = Vec::with_capacity(p);
        for (j, row) in covariance.iter().enumerate() {
            if row[j] <= 0. {
                return Err("covariance is singular; use positive regularization".into());
            }
            scales.push(row[j].sqrt());
        }
        let mut lower = vec![vec![0.; p]; p];
        for i in 0..p {
            for j in 0..=i {
                let c = finite(
                    covariance[i][j] / scales[i] / scales[j],
                    "scaled covariance",
                )?;
                let sum = (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
                if i == j {
                    let pivot = c - sum;
                    if pivot <= f64::EPSILON * 64. * p as f64 {
                        return Err(
                            "covariance is numerically singular; use positive regularization"
                                .into(),
                        );
                    }
                    lower[i][j] = pivot.sqrt();
                } else {
                    lower[i][j] = finite((c - sum) / lower[j][j], "Cholesky entry")?;
                }
            }
        }
        let log_determinant = finite(
            2. * (0..p)
                .map(|i| scales[i].ln() + lower[i][i].ln())
                .sum::<f64>(),
            "log determinant",
        )?;
        Ok(Self {
            lower,
            scales,
            log_determinant,
        })
    }
    fn squared_distance(&self, difference: &[f64]) -> Result<f64, String> {
        let p = difference.len();
        let mut whitened = vec![0.; p];
        for i in 0..p {
            let rhs = finite(difference[i] / self.scales[i], "scaled difference")?;
            let sum = (0..i).map(|j| self.lower[i][j] * whitened[j]).sum::<f64>();
            whitened[i] = finite((rhs - sum) / self.lower[i][i], "Mahalanobis coordinate")?;
        }
        finite(
            norm(whitened.into_iter()).powi(2),
            "squared Mahalanobis distance",
        )
    }
}

fn regularize(cov: &mut [Vec<f64>], lambda: f64) -> Result<(), String> {
    for (i, row) in cov.iter_mut().enumerate() {
        row[i] = finite(row[i] + lambda, "regularized covariance")?;
    }
    Ok(())
}

fn training_count(args: &[f64], data: &[Vec<f64>]) -> Result<usize, String> {
    let n = integer(args[0], "training row count")?;
    if n == 0 || n >= data.len() {
        return Err("training row count must leave a nonempty query suffix".into());
    }
    Ok(n)
}

fn mahalanobis(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 2 {
        return Err("mahalanobis requires [n_reference,lambda=0]".into());
    }
    let n = training_count(args, data)?;
    let p = data[0].len();
    let lambda = args.get(1).copied().unwrap_or(0.);
    if n < 2 || p > MAX_FEATURES || lambda < 0. {
        return Err("reference requires at least two rows, <=256 features, lambda>=0".into());
    }
    let center = Center::new(&data[..n], p)?;
    let mut cov = center.covariance(n - 1)?;
    regularize(&mut cov, lambda)?;
    let solve = Cholesky::new(&cov)?;
    let distances = data[n..]
        .iter()
        .map(|row| solve.squared_distance(&center.difference(row)?))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"mahalanobis\",\"n_reference\":{n},\"lambda\":{lambda},\"reference_mean\":{},\"covariance\":{},\"squared_distances\":{},\"preprocessing\":\"reference_prefix_only\"}}",numbers(&center.mean()?),matrix(&cov),numbers(&distances)))
}

fn class_groups(
    data: &[Vec<f64>],
    n: usize,
    p: usize,
) -> Result<BTreeMap<i64, Vec<usize>>, String> {
    let mut classes = BTreeMap::<i64, Vec<usize>>::new();
    for (i, row) in data[..n].iter().enumerate() {
        let label = row[p];
        if label.abs() > MAX_INTEGER || label.fract() != 0. {
            return Err("training classification labels must be exact signed integers".into());
        }
        classes.entry(label as i64).or_default().push(i);
    }
    if classes.len() < 2 || classes.len() > 256 {
        return Err("classification requires 2 to 256 training classes".into());
    }
    Ok(classes)
}

fn softmax(scores: &[f64]) -> Result<Vec<f64>, String> {
    let max = scores.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let exponentials = scores.iter().map(|s| (s - max).exp()).collect::<Vec<_>>();
    let total = finite(exponentials.iter().sum(), "posterior normalization")?;
    if total <= 0. {
        return Err("posterior normalization underflow".into());
    }
    Ok(exponentials.iter().map(|x| x / total).collect())
}

fn class_result(
    operation: &str,
    classes: &[f64],
    priors: &[f64],
    posterior: &[Vec<f64>],
    extra: &str,
) -> String {
    let predictions = posterior
        .iter()
        .map(|row| {
            let mut best = 0;
            for i in 1..row.len() {
                if row[i] > row[best] {
                    best = i;
                }
            }
            classes[best]
        })
        .collect::<Vec<_>>();
    format!("{{\"version\":1,\"algorithm\":\"{operation}\",\"classes\":{},\"class_priors\":{},\"predictions\":{},\"posterior_probabilities\":{},\"preprocessing\":\"training_prefix_only\",{extra}}}",numbers(classes),numbers(priors),numbers(&predictions),matrix(posterior))
}

fn discriminant(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 3 || data[0].len() < 2 {
        return Err("discriminant_analysis requires [n_train,method=0,regularization=0] and features plus label".into());
    }
    let n = training_count(args, data)?;
    let p = data[0].len() - 1;
    let method = args.get(1).copied().unwrap_or(0.);
    let lambda = args.get(2).copied().unwrap_or(0.);
    if (method != 0. && method != 1.) || lambda < 0. {
        return Err("method must be 0 LDA or 1 QDA and regularization nonnegative".into());
    }
    let groups = class_groups(data, n, p)?;
    let k = groups.len();
    if data
        .len()
        .saturating_mul(p)
        .saturating_mul(p)
        .saturating_mul(k)
        > MAX_WORK
    {
        return Err("discriminant analysis exceeds numerical work limit".into());
    }
    let classes = groups.keys().map(|x| *x as f64).collect::<Vec<_>>();
    let priors = groups
        .values()
        .map(|rows| rows.len() as f64 / n as f64)
        .collect::<Vec<_>>();
    let mut centers = Vec::with_capacity(k);
    let mut class_covariances = Vec::with_capacity(k);
    for indices in groups.values() {
        let rows = indices
            .iter()
            .map(|i| data[*i][..p].to_vec())
            .collect::<Vec<_>>();
        let center = Center::new(&rows, p)?;
        if method == 1. && indices.len() < 2 {
            return Err("QDA needs at least two training observations in every class".into());
        }
        let covariance = center.covariance(if method == 1. { indices.len() - 1 } else { 1 })?;
        centers.push(center);
        class_covariances.push(covariance);
    }
    let covariances = if method == 0. {
        if n <= k {
            return Err("LDA needs positive pooled covariance degrees of freedom".into());
        }
        let mut pooled = vec![vec![0.; p]; p];
        for covariance in &class_covariances {
            for i in 0..p {
                for j in 0..p {
                    pooled[i][j] = finite(
                        pooled[i][j] + covariance[i][j] / (n - k) as f64,
                        "pooled covariance",
                    )?;
                }
            }
        }
        regularize(&mut pooled, lambda)?;
        vec![pooled; k]
    } else {
        for covariance in &mut class_covariances {
            regularize(covariance, lambda)?;
        }
        class_covariances
    };
    let solvers = covariances
        .iter()
        .map(|cov| Cholesky::new(cov))
        .collect::<Result<Vec<_>, _>>()?;
    let mut posterior = Vec::with_capacity(data.len() - n);
    for row in &data[n..] {
        let scores = (0..k)
            .map(|c| {
                let d = solvers[c].squared_distance(&centers[c].difference(&row[..p])?)?;
                finite(
                    priors[c].ln() - 0.5 * (solvers[c].log_determinant + d),
                    "discriminant score",
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        posterior.push(softmax(&scores)?);
    }
    let means = centers
        .iter()
        .map(Center::mean)
        .collect::<Result<Vec<_>, _>>()?;
    let covariance_json = format!(
        "[{}]",
        covariances
            .iter()
            .map(|v| matrix(v))
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok(class_result("discriminant_analysis",&classes,&priors,&posterior,&format!("\"n_train\":{n},\"method\":{method},\"regularization\":{lambda},\"class_means\":{},\"class_covariances\":{covariance_json},\"covariance_estimator\":\"unbiased_within_class\"",matrix(&means))))
}

fn naive_bayes(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 2 || data[0].len() < 2 {
        return Err(
            "gaussian_naive_bayes requires [n_train,variance_floor=1e-9] and features plus label"
                .into(),
        );
    }
    let n = training_count(args, data)?;
    let p = data[0].len() - 1;
    let floor = args.get(1).copied().unwrap_or(1e-9);
    if floor <= 0. {
        return Err("variance_floor must be positive".into());
    }
    let groups = class_groups(data, n, p)?;
    let k = groups.len();
    if data.len().saturating_mul(p).saturating_mul(k) > MAX_WORK {
        return Err("naive Bayes exceeds numerical work limit".into());
    }
    let classes = groups.keys().map(|x| *x as f64).collect::<Vec<_>>();
    let priors = groups
        .values()
        .map(|g| g.len() as f64 / n as f64)
        .collect::<Vec<_>>();
    let mut centers = Vec::with_capacity(k);
    let mut variances = Vec::with_capacity(k);
    for indices in groups.values() {
        let rows = indices
            .iter()
            .map(|i| data[*i][..p].to_vec())
            .collect::<Vec<_>>();
        let center = Center::new(&rows, p)?;
        let variance = (0..p)
            .map(|j| {
                finite(
                    norm(center.rows.iter().map(|r| r[j])).powi(2) / rows.len() as f64,
                    "class variance",
                )
                .map(|v| v.max(floor))
            })
            .collect::<Result<Vec<_>, _>>()?;
        centers.push(center);
        variances.push(variance);
    }
    let mut posterior = Vec::with_capacity(data.len() - n);
    for row in &data[n..] {
        let scores = (0..k)
            .map(|c| {
                let difference = centers[c].difference(&row[..p])?;
                let terms = (0..p)
                    .map(|j| {
                        variances[c][j].ln() + (difference[j] / variances[c][j].sqrt()).powi(2)
                    })
                    .sum::<f64>();
                finite(priors[c].ln() - 0.5 * terms, "naive Bayes score")
            })
            .collect::<Result<Vec<_>, _>>()?;
        posterior.push(softmax(&scores)?);
    }
    let means = centers
        .iter()
        .map(Center::mean)
        .collect::<Result<Vec<_>, _>>()?;
    Ok(class_result("gaussian_naive_bayes",&classes,&priors,&posterior,&format!("\"n_train\":{n},\"variance_floor\":{floor},\"class_means\":{},\"class_variances\":{},\"variance_estimator\":\"population\"",matrix(&means),matrix(&variances))))
}

fn knn(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 3 || data[0].len() < 2 {
        return Err(
            "knn_predict requires [n_train,k,task] with task 0 regression or 1 classification"
                .into(),
        );
    }
    let n = training_count(args, data)?;
    let k = integer(args[1], "k")?;
    let task = args[2];
    let p = data[0].len() - 1;
    if k == 0 || k > n || (task != 0. && task != 1.) {
        return Err("k must be 1..n_train and task 0 regression or 1 classification".into());
    }
    if n.saturating_mul(data.len() - n).saturating_mul(p) > MAX_WORK {
        return Err("KNN exceeds numerical work limit".into());
    }
    if task == 1. {
        for row in &data[..n] {
            if row[p].abs() > MAX_INTEGER || row[p].fract() != 0. {
                return Err("training classification labels must be exact signed integers".into());
            }
        }
    }
    let center = Center::new(&data[..n], p)?;
    let scales = (0..p)
        .map(|j| {
            let s = finite(
                norm(center.rows.iter().map(|r| r[j])) / (n as f64).sqrt(),
                "feature scale",
            )?;
            Ok(if s == 0. { 1. } else { s })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let training = center
        .rows
        .iter()
        .map(|row| {
            row.iter()
                .zip(&scales)
                .map(|(v, s)| finite(v / s, "scaled training feature"))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut predictions = Vec::with_capacity(data.len() - n);
    let mut neighbor_indices = Vec::with_capacity(data.len() - n);
    for row in &data[n..] {
        let scaled = center
            .difference(&row[..p])?
            .iter()
            .zip(&scales)
            .map(|(v, s)| finite(v / s, "scaled query feature"))
            .collect::<Result<Vec<_>, _>>()?;
        let mut distances = training
            .iter()
            .enumerate()
            .map(|(i, r)| distance(r, &scaled).map(|d| (d, i)))
            .collect::<Result<Vec<_>, _>>()?;
        distances.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        let neighbors = &distances[..k];
        let prediction = if task == 0. {
            mean(neighbors.iter().map(|(_, i)| data[*i][p]))?
        } else {
            let mut votes = BTreeMap::<i64, usize>::new();
            for (_, i) in neighbors {
                *votes.entry(data[*i][p] as i64).or_default() += 1;
            }
            let mut best = (i64::MAX, 0);
            for (class, count) in votes {
                if count > best.1 {
                    best = (class, count);
                }
            }
            best.0 as f64
        };
        predictions.push(prediction);
        neighbor_indices.push(neighbors.iter().map(|(_, i)| *i as f64).collect::<Vec<_>>());
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"knn_predict\",\"n_train\":{n},\"k\":{k},\"task\":{task},\"predictions\":{},\"neighbor_indices\":{},\"feature_means\":{},\"feature_scales\":{},\"preprocessing\":\"training_prefix_only_population_standard_deviation\",\"distance_tie_break\":\"training_row_index\",\"vote_tie_break\":\"smallest_class_label\"}}",numbers(&predictions),matrix(&neighbor_indices),numbers(&center.mean()?),numbers(&scales)))
}

fn dbscan(data: &[Vec<f64>], p: usize, args: &[f64]) -> Result<String, String> {
    if args.len() != 2 || args[0] < 0. {
        return Err("dbscan requires [eps>=0,min_samples]".into());
    }
    let eps = args[0];
    let min_samples = integer(args[1], "min_samples")?;
    let n = data.len();
    if min_samples == 0 {
        return Err("min_samples must be positive and includes the point itself".into());
    }
    if n > 5000 || n.saturating_mul(n).saturating_mul(p) > MAX_WORK {
        return Err("DBSCAN exceeds 5000 rows or numerical work limit".into());
    }
    let mut neighbors = vec![Vec::new(); n];
    for i in 0..n {
        for j in 0..=i {
            if distance(&data[i], &data[j])? <= eps {
                neighbors[i].push(j);
                if i != j {
                    neighbors[j].push(i);
                }
            }
        }
    }
    let core = neighbors
        .iter()
        .map(|row| row.len() >= min_samples)
        .collect::<Vec<_>>();
    let mut labels = vec![-1.; n];
    let mut visited = vec![false; n];
    let mut cluster = 0;
    for i in 0..n {
        if visited[i] {
            continue;
        }
        visited[i] = true;
        if !core[i] {
            continue;
        }
        labels[i] = cluster as f64;
        let mut queue = VecDeque::new();
        let mut queued = vec![false; n];
        for &j in &neighbors[i] {
            queue.push_back(j);
            queued[j] = true;
        }
        while let Some(j) = queue.pop_front() {
            if !visited[j] {
                visited[j] = true;
                if core[j] {
                    for &v in &neighbors[j] {
                        if !queued[v] {
                            queue.push_back(v);
                            queued[v] = true;
                        }
                    }
                }
            }
            if labels[j] < 0. {
                labels[j] = cluster as f64;
            }
        }
        cluster += 1;
    }
    let core_indices = core
        .iter()
        .enumerate()
        .filter_map(|(i, c)| if *c { Some(i as f64) } else { None })
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"dbscan\",\"eps\":{eps},\"min_samples\":{min_samples},\"n_clusters\":{cluster},\"labels\":{},\"core_indices\":{},\"border_assignment\":\"first_discovered_cluster\"}}",numbers(&labels),numbers(&core_indices)))
}
