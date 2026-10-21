//! Bounded classical inference for numeric SQL data, without external runtimes.
//!
//! Normal-theory tests assume independent observations and the documented normal
//! sampling model. Rank tests use tie-corrected asymptotic distributions, with
//! continuity correction for Mann–Whitney and Wilcoxon. Binomial and Fisher
//! tests enumerate probability-ordered outcomes, rather than doubling one tail.
//! Confidence levels are bounded to [1e-10, 1-1e-10]; exact enumerations have at
//! most 100,000 trials/counts. Missing/nonfinite observations are rejected.

use std::collections::BTreeMap;

const MAX_ROWS: usize = 100_000;
const MAX_CELLS: usize = 1_000_000;
const MIN_TAIL: f64 = 1e-10;
const MAX_EXACT_INTEGER: f64 = 9_007_199_254_740_991.;

/// Dispatch a classical test; arguments and matrix layouts are operation specific.
/// Every successful result is a version-1 JSON object with explicit method labels.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0 || width > 1024 || data.len() > MAX_ROWS || data.len() > MAX_CELLS / width {
        return Err(
            "data must be a nonempty matrix within 100000 rows, 1024 columns and 1000000 cells"
                .into(),
        );
    }
    if data
        .iter()
        .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("data must be rectangular and all data and arguments finite".into());
    }
    let out = match operation {
        "one_sample_t" => {
            columns(width, 1)?;
            arguments(args, 1, 2)?;
            mean_test(
                operation,
                summary(data.iter().map(|r| r[0]))?,
                args[0],
                confidence(args.get(1))?,
            )
        }
        "paired_t" => {
            columns(width, 2)?;
            arguments(args, 0, 2)?;
            let differences = data
                .iter()
                .map(|r| checked(r[0] - r[1], "paired difference"))
                .collect::<Result<Vec<_>, _>>()?;
            mean_test(
                operation,
                summary(differences.into_iter())?,
                args.first().copied().unwrap_or(0.),
                confidence(args.get(1))?,
            )
        }
        "pooled_t" => {
            columns(width, 2)?;
            arguments(args, 0, 2)?;
            pooled_t(
                data,
                args.first().copied().unwrap_or(0.),
                confidence(args.get(1))?,
            )
        }
        "variance_test" => {
            columns(width, 1)?;
            arguments(args, 1, 2)?;
            variance_test(data, args[0], confidence(args.get(1))?)
        }
        "variance_ratio_test" => {
            columns(width, 2)?;
            arguments(args, 0, 2)?;
            variance_ratio(
                data,
                args.first().copied().unwrap_or(1.),
                confidence(args.get(1))?,
            )
        }
        "proportion_test" => {
            columns(width, 1)?;
            arguments(args, 1, 2)?;
            proportion_test(data, args[0], confidence(args.get(1))?)
        }
        "contingency_test" => {
            arguments(args, 0, 0)?;
            contingency(data)
        }
        "fisher_exact" => {
            arguments(args, 0, 0)?;
            fisher(data)
        }
        "mann_whitney" => {
            columns(width, 2)?;
            arguments(args, 0, 0)?;
            mann_whitney(data)
        }
        "wilcoxon_signed_rank" => {
            columns(width, 2)?;
            arguments(args, 0, 0)?;
            wilcoxon(data)
        }
        "kruskal_wallis" => {
            columns(width, 2)?;
            arguments(args, 0, 0)?;
            kruskal(data)
        }
        "levene_test" => {
            columns(width, 2)?;
            arguments(args, 0, 1)?;
            levene(data, args.first().copied().unwrap_or(1.))
        }
        "equivalence_t" => {
            columns(width, 1)?;
            arguments(args, 2, 3)?;
            equivalence(data, args[0], args[1], args.get(2).copied().unwrap_or(0.05))
        }
        _ => Err(format!("unknown inference operation: {operation}")),
    }?;
    // A numerical overflow must be an error, never an invalid JSON token.
    if out.contains(":NaN")
        || out.contains(":inf")
        || out.contains(":-inf")
        || out.contains("[NaN")
        || out.contains("[inf")
        || out.contains("[-inf")
        || out.contains(",NaN")
        || out.contains(",inf")
        || out.contains(",-inf")
    {
        return Err("inference result exceeds finite numerical range".into());
    }
    Ok(out)
}

fn columns(actual: usize, expected: usize) -> Result<(), String> {
    if actual == expected {
        Ok(())
    } else {
        Err(format!("operation requires {expected} columns"))
    }
}
fn arguments(args: &[f64], min: usize, max: usize) -> Result<(), String> {
    if (min..=max).contains(&args.len()) {
        Ok(())
    } else {
        Err(format!("operation requires {min} to {max} arguments"))
    }
}
fn confidence(value: Option<&f64>) -> Result<f64, String> {
    let value = value.copied().unwrap_or(0.95);
    if (MIN_TAIL..=1. - MIN_TAIL).contains(&value) {
        Ok(value)
    } else {
        Err("confidence must be between 1e-10 and 1-1e-10".into())
    }
}
fn checked(value: f64, what: &str) -> Result<f64, String> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{what} exceeds finite numerical range"))
    }
}
fn numbers(values: impl IntoIterator<Item = f64>) -> String {
    format!(
        "[{}]",
        values
            .into_iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn matrix(values: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|r| numbers(r.iter().copied()))
            .collect::<Vec<_>>()
            .join(",")
    )
}

#[derive(Clone, Copy, Default)]
struct Summary {
    n: usize,
    origin: f64,
    offset: f64,
    m2: f64,
}
impl Summary {
    fn push(&mut self, x: f64) -> Result<(), String> {
        if self.n == 0 {
            self.n = 1;
            self.origin = x;
            return Ok(());
        }
        self.n += 1;
        let relative = checked(x - self.origin, "centered value")?;
        let delta = checked(relative - self.offset, "mean change")?;
        self.offset = checked(self.offset + delta / self.n as f64, "mean")?;
        self.m2 = checked(self.m2 + delta * (relative - self.offset), "sum of squares")?;
        if self.m2 < 0. {
            return Err("negative numerical sum of squares".into());
        }
        Ok(())
    }
    fn mean(self) -> Result<f64, String> {
        checked(self.origin + self.offset, "sample mean")
    }
    fn variance(self) -> Result<f64, String> {
        if self.n < 2 {
            return Err("at least two observations required".into());
        }
        Ok(self.m2 / (self.n - 1) as f64)
    }
    fn positive_variance(self) -> Result<f64, String> {
        let v = self.variance()?;
        if v > 0. {
            Ok(v)
        } else {
            Err("test undefined with zero sample variance".into())
        }
    }
    fn centered_at(self, mu: f64) -> Result<f64, String> {
        checked((self.origin - mu) + self.offset, "mean difference")
    }
}
fn summary(values: impl Iterator<Item = f64>) -> Result<Summary, String> {
    let mut s = Summary::default();
    for x in values {
        s.push(x)?;
    }
    Ok(s)
}
fn binary_groups(data: &[Vec<f64>]) -> Result<[Summary; 2], String> {
    let mut g = [Summary::default(); 2];
    for r in data {
        if r[1] != 0. && r[1] != 1. {
            return Err("group labels must be 0 or 1".into());
        }
        g[r[1] as usize].push(r[0])?;
    }
    if g.iter().any(|s| s.n < 2) {
        return Err("each group requires at least two observations".into());
    }
    Ok(g)
}
fn mean_test(op: &str, s: Summary, null: f64, conf: f64) -> Result<String, String> {
    let se = (s.positive_variance()? / s.n as f64).sqrt();
    if se == 0. {
        return Err("standard error underflow".into());
    }
    let difference = s.centered_at(null)?;
    let t = checked(difference / se, "t statistic")?;
    let df = (s.n - 1) as f64;
    let p = 2. * student_sf(t.abs(), df)?;
    let q = student_quantile((1. - conf) / 2., df)?;
    let mean = s.mean()?;
    let ci = numbers([
        checked(
            s.origin + (s.offset - q * se),
            "lower mean confidence limit",
        )?,
        checked(
            s.origin + (s.offset + q * se),
            "upper mean confidence limit",
        )?,
    ]);
    let key = if op == "paired_t" {
        "mean_difference"
    } else {
        "mean"
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"{op}_two_sided\",\"n_samples\":{},\"{key}\":{mean},\"null_value\":{null},\"standard_error\":{se},\"t_statistic\":{t},\"degrees_of_freedom\":{df},\"p_value\":{p},\"confidence_level\":{conf},\"confidence_interval\":{ci}}}",s.n))
}
fn pooled_t(data: &[Vec<f64>], null: f64, conf: f64) -> Result<String, String> {
    let g = binary_groups(data)?;
    let df = (g[0].n + g[1].n - 2) as f64;
    let pooled = checked(g[0].m2 + g[1].m2, "pooled sum of squares")? / df;
    if pooled <= 0. {
        return Err("pooled test undefined with zero pooled variance".into());
    }
    let se = (pooled * (1. / g[0].n as f64 + 1. / g[1].n as f64)).sqrt();
    if se == 0. {
        return Err("standard error underflow".into());
    }
    let difference = checked(
        (g[0].origin - g[1].origin) + (g[0].offset - g[1].offset),
        "mean difference",
    )?;
    let t = checked((difference - null) / se, "t statistic")?;
    let p = 2. * student_sf(t.abs(), df)?;
    let q = student_quantile((1. - conf) / 2., df)?;
    let ci = numbers([
        checked(difference - q * se, "lower confidence limit")?,
        checked(difference + q * se, "upper confidence limit")?,
    ]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"pooled_two_sample_t_two_sided\",\"n_group0\":{},\"n_group1\":{},\"mean_difference\":{difference},\"null_difference\":{null},\"pooled_variance\":{pooled},\"standard_error\":{se},\"t_statistic\":{t},\"degrees_of_freedom\":{df},\"p_value\":{p},\"confidence_level\":{conf},\"confidence_interval\":{ci}}}",g[0].n,g[1].n))
}
fn variance_test(data: &[Vec<f64>], sigma: f64, conf: f64) -> Result<String, String> {
    if sigma <= 0. {
        return Err("null standard deviation must be positive".into());
    }
    let s = summary(data.iter().map(|r| r[0]))?;
    let variance = s.positive_variance()?;
    let df = (s.n - 1) as f64;
    let statistic = checked(
        (variance.sqrt() / sigma).powi(2) * df,
        "chi-square statistic",
    )?;
    let (cdf, sf) = gamma_tails(df / 2., statistic / 2.)?;
    let p = (2. * cdf.min(sf)).min(1.);
    let alpha = (1. - conf) / 2.;
    let ci = numbers([
        checked(
            s.m2 / chi_quantile(alpha, df, true)?,
            "lower variance limit",
        )?,
        checked(
            s.m2 / chi_quantile(alpha, df, false)?,
            "upper variance limit",
        )?,
    ]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"normal_variance_chi_square_equal_tail\",\"n_samples\":{},\"sample_variance\":{variance},\"null_standard_deviation\":{sigma},\"chi_square_statistic\":{statistic},\"degrees_of_freedom\":{df},\"p_value\":{p},\"confidence_level\":{conf},\"variance_confidence_interval\":{ci}}}",s.n))
}
fn variance_ratio(data: &[Vec<f64>], null: f64, conf: f64) -> Result<String, String> {
    if null <= 0. {
        return Err("null variance ratio must be positive".into());
    }
    let g = binary_groups(data)?;
    let ratio = checked(
        g[0].positive_variance()? / g[1].positive_variance()?,
        "sample variance ratio",
    )?;
    let statistic = checked(ratio / null, "F statistic")?;
    let df0 = (g[0].n - 1) as f64;
    let df1 = (g[1].n - 1) as f64;
    let (cdf, sf) = f_tails(statistic, df0, df1)?;
    let p = (2. * cdf.min(sf)).min(1.);
    let alpha = (1. - conf) / 2.;
    let ci = numbers([
        checked(
            ratio / f_quantile(alpha, df0, df1, true)?,
            "lower variance ratio limit",
        )?,
        checked(
            ratio / f_quantile(alpha, df0, df1, false)?,
            "upper variance ratio limit",
        )?,
    ]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"normal_variance_ratio_f_equal_tail\",\"variance_ratio\":{ratio},\"null_ratio\":{null},\"f_statistic\":{statistic},\"df_numerator\":{df0},\"df_denominator\":{df1},\"p_value\":{p},\"confidence_level\":{conf},\"variance_ratio_confidence_interval\":{ci}}}"))
}

fn log_factorials(n: usize) -> Vec<f64> {
    let mut values = vec![0.; n + 1];
    for i in 2..=n {
        values[i] = values[i - 1] + (i as f64).ln();
    }
    values
}
fn log_choose(n: usize, k: usize, factorials: &[f64]) -> f64 {
    factorials[n] - factorials[k] - factorials[n - k]
}
fn probability_ordered(log_probabilities: impl Iterator<Item = f64>, observed: f64) -> f64 {
    // Include equal-probability outcomes despite accumulated log-factorial error.
    let tolerance = 1e-10;
    let mut sum = 0.;
    let mut correction = 0.;
    let mut total = 0.;
    let mut total_correction = 0.;
    for lp in log_probabilities {
        let mass = lp.exp();
        let y = mass - total_correction;
        let t = total + y;
        total_correction = (t - total) - y;
        total = t;
        if lp <= observed + tolerance {
            let y = mass - correction;
            let t = sum + y;
            correction = (t - sum) - y;
            sum = t;
        }
    }
    // Enumerating the complete support lets us normalize accumulated
    // log-factorial error, which otherwise shifts large-table p-values.
    (sum / total).clamp(0., 1.)
}
fn proportion_test(data: &[Vec<f64>], null: f64, conf: f64) -> Result<String, String> {
    if !(0. ..=1.).contains(&null) || data.iter().any(|r| r[0] != 0. && r[0] != 1.) {
        return Err("proportion null must be in [0,1] and observations binary 0/1".into());
    }
    let n = data.len();
    let k = data.iter().filter(|r| r[0] == 1.).count();
    let p = if null == 0. {
        f64::from(k == 0)
    } else if null == 1. {
        f64::from(k == n)
    } else {
        let factorials = log_factorials(n);
        let lp = |x: usize| {
            log_choose(n, x, &factorials) + x as f64 * null.ln() + (n - x) as f64 * (-null).ln_1p()
        };
        probability_ordered((0..=n).map(lp), lp(k))
    };
    let phat = k as f64 / n as f64;
    let z = normal_quantile((1. - conf) / 2.)?;
    let z2 = z * z;
    let den = 1. + z2 / n as f64;
    let center = (phat + z2 / (2. * n as f64)) / den;
    let half = z * ((phat * (1. - phat) + z2 / (4. * n as f64)) / n as f64).sqrt() / den;
    let ci = numbers([(center - half).max(0.), (center + half).min(1.)]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"exact_binomial_probability_ordered_two_sided\",\"n_trials\":{n},\"successes\":{k},\"sample_proportion\":{phat},\"null_proportion\":{null},\"p_value\":{p},\"confidence_level\":{conf},\"confidence_interval_method\":\"wilson_score\",\"confidence_interval\":{ci}}}"))
}
fn count_table(data: &[Vec<f64>], cap: f64) -> Result<(Vec<f64>, Vec<f64>, f64), String> {
    if data.len() < 2 || data[0].len() < 2 {
        return Err("contingency tables require at least two rows and columns".into());
    }
    let mut rows = vec![0.; data.len()];
    let mut cols = vec![0.; data[0].len()];
    let mut total = 0.;
    for (i, r) in data.iter().enumerate() {
        for (j, &x) in r.iter().enumerate() {
            if x < 0. || x.fract() != 0. || x > cap {
                return Err("counts must be nonnegative integers within count limits".into());
            }
            rows[i] += x;
            cols[j] += x;
            total += x;
            if total > cap {
                return Err(format!("total table count exceeds {cap}"));
            }
        }
    }
    if rows.contains(&0.) || cols.contains(&0.) {
        return Err("all row and column marginal totals must be positive".into());
    }
    Ok((rows, cols, total))
}
fn contingency(data: &[Vec<f64>]) -> Result<String, String> {
    let (rows, cols, total) = count_table(data, 1e12)?;
    let expected = rows
        .iter()
        .map(|r| cols.iter().map(|c| r / total * c).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut statistic = 0.;
    let residuals = data
        .iter()
        .zip(&expected)
        .map(|(r, e)| {
            r.iter()
                .zip(e)
                .map(|(x, m)| {
                    let z = (x - m) / m.sqrt();
                    statistic += z * z;
                    z
                })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    checked(statistic, "Pearson chi-square statistic")?;
    let df = ((rows.len() - 1) * (cols.len() - 1)) as f64;
    let p = gamma_tails(df / 2., statistic / 2.)?.1;
    let v = (statistic / (total * (rows.len() - 1).min(cols.len() - 1) as f64)).sqrt();
    let small = expected.iter().flatten().filter(|&&x| x < 5.).count();
    Ok(format!("{{\"version\":1,\"algorithm\":\"pearson_chi_square_uncorrected_asymptotic\",\"total_count\":{total},\"chi_square_statistic\":{statistic},\"degrees_of_freedom\":{df},\"p_value\":{p},\"cramers_v\":{v},\"expected_counts\":{},\"pearson_residuals\":{},\"expected_cells_below_5\":{small},\"yates_correction\":false}}",matrix(&expected),matrix(&residuals)))
}
fn fisher(data: &[Vec<f64>]) -> Result<String, String> {
    if data.len() != 2 || data[0].len() != 2 {
        return Err("Fisher exact requires a 2x2 count table".into());
    }
    let (rows, cols, total) = count_table(data, MAX_ROWS as f64)?;
    let n = total as usize;
    let row = rows[0] as usize;
    let col = cols[0] as usize;
    let lo = row.saturating_sub(n - col);
    let hi = row.min(col);
    let factorials = log_factorials(n);
    let lp = |x: usize| {
        log_choose(col, x, &factorials) + log_choose(n - col, row - x, &factorials)
            - log_choose(n, row, &factorials)
    };
    let p = probability_ordered((lo..=hi).map(lp), lp(data[0][0] as usize));
    let num = data[0][0] * data[1][1];
    let den = data[0][1] * data[1][0];
    let odds = if den == 0. {
        "null".into()
    } else {
        (num / den).to_string()
    };
    let infinite = den == 0. && num > 0.;
    Ok(format!("{{\"version\":1,\"algorithm\":\"fisher_fixed_margins_probability_ordered_two_sided\",\"p_value\":{p},\"odds_ratio\":{odds},\"odds_ratio_infinite\":{infinite},\"odds_ratio_estimator\":\"sample_cross_product\",\"enumerated_tables\":{}}}",hi-lo+1))
}

fn group_id(x: f64) -> Result<u64, String> {
    if !(0. ..=MAX_EXACT_INTEGER).contains(&x) || x.fract() != 0. {
        Err("group labels must be exact nonnegative integers".into())
    } else {
        Ok(x as u64)
    }
}
fn groups(data: &[Vec<f64>]) -> Result<BTreeMap<u64, Vec<f64>>, String> {
    let mut g = BTreeMap::<u64, Vec<f64>>::new();
    for r in data {
        g.entry(group_id(r[1])?).or_default().push(r[0]);
    }
    if g.len() < 2 {
        return Err("at least two groups required".into());
    }
    Ok(g)
}
// Ranks in original order, together with sum(t^3-t) across tie blocks.
fn ranks(values: &[f64]) -> (Vec<f64>, f64) {
    let mut order = (0..values.len()).collect::<Vec<_>>();
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let mut rank = vec![0.; values.len()];
    let mut ties = 0.;
    let mut i = 0;
    while i < order.len() {
        let start = i;
        while i < order.len() && values[order[i]] == values[order[start]] {
            i += 1;
        }
        let t = (i - start) as f64;
        ties += t * t * t - t;
        let r = (start + i + 1) as f64 / 2.;
        for &j in &order[start..i] {
            rank[j] = r;
        }
    }
    (rank, ties)
}
fn mann_whitney(data: &[Vec<f64>]) -> Result<String, String> {
    let mut n = [0usize; 2];
    for r in data {
        if r[1] != 0. && r[1] != 1. {
            return Err("Mann-Whitney groups must be 0 or 1".into());
        }
        n[r[1] as usize] += 1;
    }
    if n.contains(&0) {
        return Err("both groups must be nonempty".into());
    }
    let values = data.iter().map(|r| r[0]).collect::<Vec<_>>();
    let (rank, ties) = ranks(&values);
    let rank0 = data
        .iter()
        .zip(rank)
        .filter(|(r, _)| r[1] == 0.)
        .map(|(_, v)| v)
        .sum::<f64>();
    let n0 = n[0] as f64;
    let n1 = n[1] as f64;
    let total = (n[0] + n[1]) as f64;
    let u = rank0 - n0 * (n0 + 1.) / 2.;
    let mean = n0 * n1 / 2.;
    let variance = n0 * n1 / 12. * ((total + 1.) - ties / (total * (total - 1.)));
    if variance <= 0. {
        return Err("Mann-Whitney undefined when all observations are tied".into());
    }
    let z = ((u - mean).abs() - 0.5).max(0.) / variance.sqrt();
    let p = 2. * normal_sf(z)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"mann_whitney_tie_corrected_normal_approximation\",\"n_group0\":{},\"n_group1\":{},\"u_statistic\":{u},\"z_absolute\":{z},\"p_value\":{p},\"continuity_correction\":true}}",n[0],n[1]))
}
fn wilcoxon(data: &[Vec<f64>]) -> Result<String, String> {
    let diff = data
        .iter()
        .map(|r| checked(r[0] - r[1], "paired difference"))
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .filter(|&v| v != 0.)
        .collect::<Vec<_>>();
    if diff.is_empty() {
        return Err("Wilcoxon requires at least one nonzero difference".into());
    }
    let n = diff.len() as f64;
    let absolute = diff.iter().map(|x| x.abs()).collect::<Vec<_>>();
    let (rank, ties) = ranks(&absolute);
    let plus = diff
        .iter()
        .zip(rank.iter())
        .filter(|(d, _)| **d > 0.)
        .map(|(_, r)| r)
        .sum::<f64>();
    let minus = n * (n + 1.) / 2. - plus;
    let statistic = plus.min(minus);
    let variance = n * (n + 1.) * (2. * n + 1.) / 24. - ties / 48.;
    if variance <= 0. {
        return Err("Wilcoxon variance is zero".into());
    }
    let z = ((plus - n * (n + 1.) / 4.).abs() - 0.5).max(0.) / variance.sqrt();
    let p = 2. * normal_sf(z)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"wilcoxon_tie_corrected_normal_approximation\",\"zero_method\":\"wilcox_drop\",\"n_nonzero\":{},\"n_zero_dropped\":{},\"positive_rank_sum\":{plus},\"negative_rank_sum\":{minus},\"rank_statistic\":{statistic},\"z_absolute\":{z},\"p_value\":{p},\"continuity_correction\":true}}",diff.len(),data.len()-diff.len()))
}
fn kruskal(data: &[Vec<f64>]) -> Result<String, String> {
    let g = groups(data)?;
    let values = data.iter().map(|r| r[0]).collect::<Vec<_>>();
    let (rank, ties) = ranks(&values);
    let mut sums = BTreeMap::<u64, f64>::new();
    for (r, v) in data.iter().zip(rank) {
        *sums.entry(group_id(r[1])?).or_default() += v;
    }
    let n = data.len() as f64;
    let correction = 1. - ties / (n * n * n - n);
    if correction <= 0. {
        return Err("Kruskal-Wallis undefined when all observations are tied".into());
    }
    let h = ((12. / (n * (n + 1.))
        * sums
            .iter()
            .map(|(id, s)| s * s / g[id].len() as f64)
            .sum::<f64>()
        - 3. * (n + 1.))
        / correction)
        .max(0.);
    let df = (g.len() - 1) as f64;
    let p = gamma_tails(df / 2., h / 2.)?.1;
    Ok(format!("{{\"version\":1,\"algorithm\":\"kruskal_wallis_tie_corrected_chi_square_approximation\",\"n_samples\":{},\"n_groups\":{},\"h_statistic\":{h},\"degrees_of_freedom\":{df},\"tie_correction\":{correction},\"p_value\":{p}}}",data.len(),g.len()))
}
fn levene(data: &[Vec<f64>], center: f64) -> Result<String, String> {
    if center != 0. && center != 1. {
        return Err("Levene center must be 0 (mean) or 1 (median)".into());
    }
    let g = groups(data)?;
    if data.len() <= g.len() {
        return Err("Levene test requires residual degrees of freedom".into());
    }
    let mut deviations = BTreeMap::new();
    let mut overall = Summary::default();
    for (id, values) in &g {
        if values.len() < 2 {
            return Err("Levene requires at least two observations in every group".into());
        }
        let c = if center == 0. {
            summary(values.iter().copied())?.mean()?
        } else {
            let mut v = values.clone();
            v.sort_by(f64::total_cmp);
            let m = v.len() / 2;
            if v.len() % 2 == 0 {
                v[m - 1] / 2. + v[m] / 2.
            } else {
                v[m]
            }
        };
        let mut s = Summary::default();
        for &x in values {
            let d = checked(x - c, "absolute deviation")?.abs();
            s.push(d)?;
            overall.push(d)?;
        }
        deviations.insert(*id, s);
    }
    let within = checked(
        deviations.values().map(|s| s.m2).sum::<f64>(),
        "within-deviation sum of squares",
    )?;
    if within <= 0. {
        return Err("Levene F undefined with zero within-group deviation variance".into());
    }
    let between = checked(
        deviations
            .values()
            .map(|s| {
                s.n as f64 * ((s.origin - overall.origin) + (s.offset - overall.offset)).powi(2)
            })
            .sum::<f64>(),
        "between-deviation sum of squares",
    )?;
    let df0 = (g.len() - 1) as f64;
    let df1 = (data.len() - g.len()) as f64;
    let statistic = checked((between / df0) / (within / df1), "Levene F statistic")?;
    let p = f_tails(statistic, df0, df1)?.1;
    let name = if center == 0. { "mean" } else { "median" };
    Ok(format!("{{\"version\":1,\"algorithm\":\"levene_absolute_deviation_anova\",\"center\":\"{name}\",\"n_samples\":{},\"n_groups\":{},\"f_statistic\":{statistic},\"df_numerator\":{df0},\"df_denominator\":{df1},\"p_value\":{p}}}",data.len(),g.len()))
}
fn equivalence(data: &[Vec<f64>], lower: f64, upper: f64, alpha: f64) -> Result<String, String> {
    if lower >= upper || !(MIN_TAIL..0.5).contains(&alpha) {
        return Err("equivalence requires lower < upper and 1e-10 <= alpha < 0.5".into());
    }
    let s = summary(data.iter().map(|r| r[0]))?;
    let se = (s.positive_variance()? / s.n as f64).sqrt();
    if se == 0. {
        return Err("standard error underflow".into());
    }
    let df = (s.n - 1) as f64;
    let tl = checked(s.centered_at(lower)? / se, "lower t statistic")?;
    let tu = checked(s.centered_at(upper)? / se, "upper t statistic")?;
    let pl = student_signed_sf(tl, df)?;
    let pu = student_signed_sf(-tu, df)?;
    let p = pl.max(pu);
    let equivalent = p < alpha;
    let q = student_quantile(alpha, df)?;
    let ci = numbers([
        checked(s.origin + (s.offset - q * se), "lower equivalence CI")?,
        checked(s.origin + (s.offset + q * se), "upper equivalence CI")?,
    ]);
    let conf = 1. - 2. * alpha;
    let mean = s.mean()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"one_sample_tost\",\"mean\":{mean},\"lower_equivalence_bound\":{lower},\"upper_equivalence_bound\":{upper},\"alpha\":{alpha},\"t_lower\":{tl},\"t_upper\":{tu},\"degrees_of_freedom\":{df},\"p_lower\":{pl},\"p_upper\":{pu},\"p_value\":{p},\"equivalent\":{equivalent},\"confidence_level\":{conf},\"confidence_interval\":{ci}}}"))
}

// Special functions are private to this module. All continued fractions and
// series have iteration caps and return errors rather than unconverged values.
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
    let z = z - 1.;
    let mut a = 0.999_999_999_999_809_9;
    for (i, c) in C.iter().enumerate() {
        a += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * (2. * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + a.ln()
}
fn beta_fraction(a: f64, b: f64, x: f64) -> Result<f64, String> {
    let nz = |v: f64| {
        if v.abs() < 1e-300 {
            1e-300_f64.copysign(v)
        } else {
            v
        }
    };
    let mut c = 1.;
    let mut d = 1. / nz(1. - (a + b) * x / (a + 1.));
    let mut h = d;
    for i in 1..=10_000 {
        let m = i as f64;
        let even = m * (b - m) * x / ((a + 2. * m - 1.) * (a + 2. * m));
        d = 1. / nz(1. + even * d);
        c = nz(1. + even / c);
        h *= d * c;
        let odd = -(a + m) * (a + b + m) * x / ((a + 2. * m) * (a + 2. * m + 1.));
        d = 1. / nz(1. + odd * d);
        c = nz(1. + odd / c);
        let delta = d * c;
        h *= delta;
        if (delta - 1.).abs() < 2e-14 {
            return checked(h, "incomplete beta fraction");
        }
    }
    Err("incomplete beta fraction failed to converge".into())
}
fn beta(x: f64, a: f64, b: f64) -> Result<f64, String> {
    if x <= 0. {
        return Ok(0.);
    }
    if x >= 1. {
        return Ok(1.);
    }
    let front =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    let p = if x < (a + 1.) / (a + b + 2.) {
        front * beta_fraction(a, b, x)? / a
    } else {
        1. - front * beta_fraction(b, a, 1. - x)? / b
    };
    probability(p)
}
fn probability(p: f64) -> Result<f64, String> {
    if p.is_finite() && (-1e-10..=1. + 1e-10).contains(&p) {
        Ok(p.clamp(0., 1.))
    } else {
        Err("probability outside numerical range".into())
    }
}
fn gamma_tails(a: f64, x: f64) -> Result<(f64, f64), String> {
    if x == 0. {
        return Ok((0., 1.));
    }
    if x < 0. || !x.is_finite() || a <= 0. || !a.is_finite() {
        return Err("invalid gamma parameters".into());
    }
    let front = (a * x.ln() - x - log_gamma(a)).exp();
    if x < a + 1. {
        let mut term = 1. / a;
        let mut sum = term;
        for i in 1..=20_000 {
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() < sum.abs() * 2e-14 {
                let p = probability(sum * front)?;
                return Ok((p, 1. - p));
            }
        }
    } else {
        let nz = |v: f64| {
            if v.abs() < 1e-300 {
                1e-300_f64.copysign(v)
            } else {
                v
            }
        };
        let mut b = x + 1. - a;
        let mut c = 1e300;
        let mut d = 1. / nz(b);
        let mut h = d;
        for i in 1..=20_000 {
            let i = i as f64;
            let an = -i * (i - a);
            b += 2.;
            d = 1. / nz(an * d + b);
            c = nz(b + an / c);
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-14 {
                let q = probability(front * h)?;
                return Ok((1. - q, q));
            }
        }
    }
    Err("incomplete gamma failed to converge".into())
}
fn student_sf(t: f64, df: f64) -> Result<f64, String> {
    let r = t / df.sqrt();
    let r2 = r * r;
    if r <= 1. {
        Ok(0.5 - 0.5 * beta(r2 / (1. + r2), 0.5, df / 2.)?)
    } else {
        Ok(0.5 * beta(1. / (1. + r2), df / 2., 0.5)?)
    }
}
fn student_signed_sf(t: f64, df: f64) -> Result<f64, String> {
    let tail = student_sf(t.abs(), df)?;
    Ok(if t >= 0. { tail } else { 1. - tail })
}
fn f_tails(x: f64, df0: f64, df1: f64) -> Result<(f64, f64), String> {
    if x == 0. {
        return Ok((0., 1.));
    }
    let r = x * (df0 / df1);
    if r <= 1. {
        let p = beta(r / (1. + r), df0 / 2., df1 / 2.)?;
        Ok((p, 1. - p))
    } else {
        let q = beta(1. / (1. + r), df1 / 2., df0 / 2.)?;
        Ok((1. - q, q))
    }
}
fn normal_sf(z: f64) -> Result<f64, String> {
    Ok(0.5 * gamma_tails(0.5, z * z / 2.)?.1)
}
// Inversion brackets the positive axis and then bisects; stopping uses the
// bracket width rather than a probability threshold (which loses tail digits).
fn quantile(
    target: f64,
    start: f64,
    decreasing: bool,
    prob: impl Fn(f64) -> Result<f64, String>,
) -> Result<f64, String> {
    let mut lo = 0.;
    let mut hi = start.max(1.);
    let bracketed = |p: f64| if decreasing { p <= target } else { p >= target };
    let mut found = false;
    for _ in 0..512 {
        if bracketed(prob(hi)?) {
            found = true;
            break;
        }
        hi *= 2.;
        if !hi.is_finite() {
            break;
        }
    }
    if !found {
        return Err("quantile could not be bracketed in finite range".into());
    }
    for _ in 0..512 {
        let mid = lo + (hi - lo) / 2.;
        if mid == lo || mid == hi || hi - lo <= 2e-14 * hi.max(f64::MIN_POSITIVE) {
            return Ok(lo + (hi - lo) / 2.);
        }
        if bracketed(prob(mid)?) {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    Err("quantile bisection failed to converge".into())
}
fn student_quantile(tail: f64, df: f64) -> Result<f64, String> {
    quantile(tail, 1., true, |x| student_sf(x, df))
}
fn chi_quantile(tail: f64, df: f64, upper: bool) -> Result<f64, String> {
    quantile(tail, df, upper, |x| {
        let (p, q) = gamma_tails(df / 2., x / 2.)?;
        Ok(if upper { q } else { p })
    })
}
fn f_quantile(tail: f64, df0: f64, df1: f64, upper: bool) -> Result<f64, String> {
    quantile(tail, 1., upper, |x| {
        let (p, q) = f_tails(x, df0, df1)?;
        Ok(if upper { q } else { p })
    })
}
fn normal_quantile(tail: f64) -> Result<f64, String> {
    quantile(tail, 1., true, normal_sf)
}
