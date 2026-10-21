//! Statistical summaries and bounded inference for numerical SQL workflows.
//!
//! Welch and ANOVA p-values use a converged continued fraction for regularized
//! incomplete beta (NIST DLMF 8.17); p-values are numerical approximations.
//! Random intercept inference is restricted to a balanced one-way ANOVA variance
//! estimator. Capability uses overall sample SD, not a fitted subgroup model.

use std::collections::BTreeMap;

const MAX_CELLS: usize = 1_000_000;
const MAX_ROWS: usize = 100_000;
const MAX_INTEGER: f64 = 9_007_199_254_740_991.;

#[derive(Clone, Copy, Default)]
struct Moments {
    n: usize,
    mean: f64,
    origin: f64,
    mean_offset: f64,
    m2: f64,
    min: f64,
    max: f64,
}

impl Moments {
    fn push(&mut self, value: f64) -> Result<(), String> {
        if self.n == 0 {
            *self = Self {
                n: 1,
                mean: value,
                origin: value,
                mean_offset: 0.,
                m2: 0.,
                min: value,
                max: value,
            };
            return Ok(());
        }
        self.n += 1;
        let offset = finite(value - self.origin, "origin-relative value")?;
        let delta = finite(offset - self.mean_offset, "centered value")?;
        self.mean_offset = finite(self.mean_offset + delta / self.n as f64, "mean offset")?;
        self.mean = finite(self.origin + self.mean_offset, "mean")?;
        self.m2 = finite(
            self.m2 + delta * (offset - self.mean_offset),
            "sum of squares",
        )?;
        self.min = self.min.min(value);
        self.max = self.max.max(value);
        Ok(())
    }

    fn variance(&self) -> f64 {
        self.m2 / (self.n - 1) as f64
    }
}

/// Numeric rows and parameters follow the operation-specific SQL documentation.
pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0 || width > 1024 || data.len() > MAX_ROWS || data.len() > MAX_CELLS / width {
        return Err("data must be a nonempty matrix within row, column, and cell limits".into());
    }
    if data
        .iter()
        .any(|row| row.len() != width || row.iter().any(|value| !value.is_finite()))
        || parameters.iter().any(|value| !value.is_finite())
    {
        return Err("data must be rectangular; all data and parameters must be finite".into());
    }
    match operation {
        "summarize" => {
            require_parameters(parameters, 0)?;
            summarize(data)
        }
        "welch_t_test" => {
            require_columns(width, 2)?;
            require_parameters(parameters, 0)?;
            welch(data)
        }
        "one_way_anova" => {
            require_columns(width, 2)?;
            require_parameters(parameters, 0)?;
            anova(data, false)
        }
        "adjust_pvalues" => {
            require_columns(width, 1)?;
            require_parameters(parameters, 1)?;
            adjust(data, parameters[0])
        }
        "process_capability" => {
            require_columns(width, 1)?;
            if parameters.len() != 2 && parameters.len() != 3 {
                return Err(
                    "capability requires [lower limit, upper limit, optional within-process sigma]"
                        .into(),
                );
            }
            capability(data, parameters)
        }
        "kaplan_meier" => {
            require_columns(width, 2)?;
            require_parameters(parameters, 0)?;
            kaplan_meier(data)
        }
        "random_intercept" => {
            require_columns(width, 2)?;
            require_parameters(parameters, 0)?;
            anova(data, true)
        }
        _ => Err(format!("unknown statistics operation: {operation}")),
    }
}

fn require_columns(actual: usize, expected: usize) -> Result<(), String> {
    if actual == expected {
        Ok(())
    } else {
        Err(format!("operation requires {expected} columns"))
    }
}

fn require_parameters(parameters: &[f64], expected: usize) -> Result<(), String> {
    if parameters.len() == expected {
        Ok(())
    } else {
        Err(format!("operation requires {expected} parameters"))
    }
}

fn finite(value: f64, name: &str) -> Result<f64, String> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}

fn integer(value: f64) -> Result<u64, String> {
    if !(0. ..=MAX_INTEGER).contains(&value) || value.fract() != 0. {
        Err("group ids must be exact nonnegative integers".into())
    } else {
        Ok(value as u64)
    }
}

fn moments(values: impl Iterator<Item = f64>) -> Result<Moments, String> {
    let mut output = Moments::default();
    for value in values {
        output.push(value)?;
    }
    Ok(output)
}

fn numbers(values: impl IntoIterator<Item = f64>) -> String {
    format!(
        "[{}]",
        values
            .into_iter()
            .map(|value| value.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn optional(value: Option<f64>) -> String {
    value.map_or_else(|| "null".to_string(), |value| value.to_string())
}

fn summarize(data: &[Vec<f64>]) -> Result<String, String> {
    let statistics = (0..data[0].len())
        .map(|j| moments(data.iter().map(|row| row[j])))
        .collect::<Result<Vec<_>, _>>()?;
    let variances = statistics
        .iter()
        .map(|s| optional((s.n > 1).then(|| s.variance())))
        .collect::<Vec<_>>()
        .join(",");
    let deviations = statistics
        .iter()
        .map(|s| optional((s.n > 1).then(|| s.variance().sqrt())))
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!("{{\"version\":1,\"algorithm\":\"sample_summary\",\"n_samples\":{},\"n_features\":{},\"means\":{},\"sample_variances\":[{variances}],\"standard_deviations\":[{deviations}],\"minima\":{},\"maxima\":{}}}", data.len(), data[0].len(), numbers(statistics.iter().map(|s| s.mean)), numbers(statistics.iter().map(|s| s.min)), numbers(statistics.iter().map(|s| s.max))))
}

fn welch(data: &[Vec<f64>]) -> Result<String, String> {
    let mut groups = [Moments::default(); 2];
    for row in data {
        if row[1] != 0. && row[1] != 1. {
            return Err("Welch groups must be 0 or 1".into());
        }
        groups[row[1] as usize].push(row[0])?;
    }
    if groups.iter().any(|group| group.n < 2) {
        return Err("Welch test requires at least two rows in each group".into());
    }
    let a = groups[0].variance() / groups[0].n as f64;
    let b = groups[1].variance() / groups[1].n as f64;
    let variance = finite(a + b, "variance of mean difference")?;
    if variance == 0. {
        return Err("Welch test is undefined with zero variance in both groups".into());
    }
    let difference = finite(
        (groups[0].origin - groups[1].origin) + (groups[0].mean_offset - groups[1].mean_offset),
        "mean difference",
    )?;
    let t = finite(difference / variance.sqrt(), "t statistic")?;
    // Normalize before squaring to avoid underflow/overflow in Satterthwaite df.
    let df = finite(
        1. / ((a / variance).powi(2) / (groups[0].n - 1) as f64
            + (b / variance).powi(2) / (groups[1].n - 1) as f64),
        "Welch degrees of freedom",
    )?;
    let ratio = t.abs() / df.sqrt();
    // Evaluate the complementary representation near t=0 to retain small tails.
    let ratio_squared = ratio * ratio;
    let p_value = if ratio < 1. {
        1. - beta_regularized(ratio_squared / (1. + ratio_squared), 0.5, df / 2.)?
    } else {
        beta_regularized(1. / (1. + ratio_squared), df / 2., 0.5)?
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"welch_two_sided\",\"n_group0\":{},\"n_group1\":{},\"mean_group0\":{},\"mean_group1\":{},\"mean_difference\":{difference},\"t_statistic\":{t},\"degrees_of_freedom\":{df},\"p_value\":{p_value}}}", groups[0].n, groups[1].n, groups[0].mean, groups[1].mean))
}

fn anova(data: &[Vec<f64>], random_intercept: bool) -> Result<String, String> {
    let mut groups = BTreeMap::<u64, Moments>::new();
    let overall = moments(data.iter().map(|row| row[0]))?;
    for row in data {
        groups.entry(integer(row[1])?).or_default().push(row[0])?;
    }
    let k = groups.len();
    if k < 2 || data.len() <= k {
        return Err("ANOVA requires at least two groups and residual degrees of freedom".into());
    }
    let df_between = (k - 1) as f64;
    let df_within = (data.len() - k) as f64;
    let within = finite(
        groups.values().map(|group| group.m2).sum::<f64>(),
        "within-group sum of squares",
    )?;
    let between = finite(
        groups
            .values()
            .map(|group| {
                group.n as f64
                    * ((group.origin - overall.origin) + (group.mean_offset - overall.mean_offset))
                        .powi(2)
            })
            .sum::<f64>(),
        "between-group sum of squares",
    )?;
    let ms_between = between / df_between;
    let ms_within = within / df_within;
    if random_intercept {
        let size = groups.values().next().map_or(0, |group| group.n);
        if size < 2 || groups.values().any(|group| group.n != size) {
            return Err(
                "random intercept ANOVA requires balanced groups with at least two replicates"
                    .into(),
            );
        }
        let between_variance = finite(
            ((ms_between - ms_within) / size as f64).max(0.),
            "between-group variance",
        )?;
        let total = finite(between_variance + ms_within, "total variance")?;
        let icc = optional((total > 0.).then(|| between_variance / total));
        return Ok(format!("{{\"version\":1,\"algorithm\":\"balanced_random_intercept_anova\",\"n_samples\":{},\"n_groups\":{k},\"replicates_per_group\":{size},\"grand_mean\":{},\"within_group_variance\":{ms_within},\"between_group_variance\":{between_variance},\"icc\":{icc},\"negative_between_estimates_truncated_to_zero\":true}}", data.len(), overall.mean));
    }
    if ms_within == 0. {
        return Err("ANOVA F statistic is undefined with zero within-group variance".into());
    }
    let f = finite(ms_between / ms_within, "F statistic")?;
    let x = 1. / (1. + (df_between / df_within) * f);
    let p_value = beta_regularized(x, df_within / 2., df_between / 2.)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"one_way_anova\",\"n_samples\":{},\"n_groups\":{k},\"grand_mean\":{},\"sum_squares_between\":{between},\"sum_squares_within\":{within},\"df_between\":{df_between},\"df_within\":{df_within},\"f_statistic\":{f},\"p_value\":{p_value}}}", data.len(), overall.mean))
}

fn adjust(data: &[Vec<f64>], method: f64) -> Result<String, String> {
    if method != 0. && method != 1. {
        return Err("adjustment method must be 0 (Holm) or 1 (Benjamini-Hochberg)".into());
    }
    if data.iter().any(|row| !(0. ..=1.).contains(&row[0])) {
        return Err("p-values must be between zero and one".into());
    }
    let n = data.len();
    let mut ordered = data
        .iter()
        .enumerate()
        .map(|(i, row)| (i, row[0]))
        .collect::<Vec<_>>();
    ordered.sort_by(|(i, a), (j, b)| a.total_cmp(b).then(i.cmp(j)));
    let mut adjusted = vec![0.; n];
    if method == 0. {
        let mut previous = 0_f64;
        for (rank, &(i, p)) in ordered.iter().enumerate() {
            previous = previous.max(((n - rank) as f64 * p).min(1.));
            adjusted[i] = previous;
        }
    } else {
        let mut previous = 1_f64;
        for (rank, &(i, p)) in ordered.iter().enumerate().rev() {
            previous = previous.min(n as f64 * p / (rank + 1) as f64);
            adjusted[i] = previous;
        }
    }
    let name = if method == 0. {
        "holm"
    } else {
        "benjamini_hochberg"
    };
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"{name}\",\"adjusted_p_values\":{}}}",
        numbers(adjusted)
    ))
}

fn capability(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    if data.len() < 2 || parameters[0] >= parameters[1] {
        return Err(
            "capability requires at least two rows and lower < upper specification limit".into(),
        );
    }
    let summary = moments(data.iter().map(|row| row[0]))?;
    let sd = summary.variance().sqrt();
    if sd == 0. {
        return Err("capability is undefined with zero sample standard deviation".into());
    }
    if parameters.len() == 3 && parameters[2] <= 0. {
        return Err("supplied within-process sigma must be positive".into());
    }
    let width = finite(parameters[1] - parameters[0], "specification width")?;
    let upper_distance = finite(
        (parameters[1] - summary.origin) - summary.mean_offset,
        "upper specification distance",
    )?;
    let lower_distance = finite(
        (summary.origin - parameters[0]) + summary.mean_offset,
        "lower specification distance",
    )?;
    let pp = finite((width / sd) / 6., "Pp")?;
    let ppk = finite((upper_distance / sd).min(lower_distance / sd) / 3., "Ppk")?;
    let (cp, cpk) = if parameters.len() == 3 {
        let within_sigma = parameters[2];
        (
            optional(Some(finite((width / within_sigma) / 6., "Cp")?)),
            optional(Some(finite(
                (upper_distance / within_sigma).min(lower_distance / within_sigma) / 3.,
                "Cpk",
            )?)),
        )
    } else {
        ("null".to_string(), "null".to_string())
    };
    let within_sigma = optional(parameters.get(2).copied());
    let lower = finite(
        summary.origin + (summary.mean_offset - 3. * sd),
        "lower reference limit",
    )?;
    let upper = finite(
        summary.origin + (summary.mean_offset + 3. * sd),
        "upper reference limit",
    )?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"overall_sample_capability\",\"n_samples\":{},\"mean\":{},\"sample_standard_deviation\":{sd},\"lower_specification_limit\":{},\"upper_specification_limit\":{},\"cp\":{cp},\"cpk\":{cpk},\"pp\":{pp},\"ppk\":{ppk},\"within_process_sigma\":{within_sigma},\"sigma_estimator\":\"overall_sample_standard_deviation\",\"individual_lower_reference_limit\":{lower},\"individual_upper_reference_limit\":{upper},\"reference_limit_model\":\"descriptive_mean_plus_minus_3_overall_sample_sd\"}}", data.len(), summary.mean, parameters[0], parameters[1]))
}

fn kaplan_meier(data: &[Vec<f64>]) -> Result<String, String> {
    if data
        .iter()
        .any(|row| row[0] < 0. || (row[1] != 0. && row[1] != 1.))
    {
        return Err("survival times must be nonnegative; event indicators must be 0 or 1".into());
    }
    let mut ordered = data
        .iter()
        .map(|row| (row[0], row[1] == 1.))
        .collect::<Vec<_>>();
    ordered.sort_by(|a, b| a.0.total_cmp(&b.0));
    let mut times = Vec::new();
    let mut risks = Vec::new();
    let mut events = Vec::new();
    let mut censors = Vec::new();
    let mut survivals = Vec::new();
    let mut variances = Vec::new();
    let mut risk = ordered.len();
    let mut survival = 1.;
    let mut greenwood_sum = 0.;
    let mut median = None;
    let mut i = 0;
    while i < ordered.len() {
        let time = ordered[i].0;
        let start = i;
        let mut deaths = 0;
        while i < ordered.len() && ordered[i].0 == time {
            deaths += usize::from(ordered[i].1);
            i += 1;
        }
        let censored = i - start - deaths;
        if deaths > 0 {
            survival *= 1. - deaths as f64 / risk as f64;
            if deaths < risk {
                greenwood_sum += deaths as f64 / (risk as f64 * (risk - deaths) as f64);
            }
            if survival <= 0.5 && median.is_none() {
                median = Some(time);
            }
        }
        times.push(time);
        risks.push(risk as f64);
        events.push(deaths as f64);
        censors.push(censored as f64);
        survivals.push(survival);
        variances.push(if survival == 0. {
            0.
        } else {
            finite(survival * survival * greenwood_sum, "Greenwood variance")?
        });
        risk -= deaths + censored;
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"kaplan_meier\",\"n_samples\":{},\"times\":{},\"at_risk\":{},\"events\":{},\"censored\":{},\"survival\":{},\"greenwood_variance\":{},\"median_survival\":{}}}", data.len(), numbers(times), numbers(risks), numbers(events), numbers(censors), numbers(survivals), numbers(variances), optional(median)))
}

// Lanczos approximation with g=7, used on positive beta shape parameters.
fn log_gamma(z: f64) -> f64 {
    const COEFFICIENTS: [f64; 8] = [
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
    let mut series = 0.999_999_999_999_809_9;
    for (i, value) in COEFFICIENTS.iter().enumerate() {
        series += value / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * (2. * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + series.ln()
}

fn beta_regularized(x: f64, a: f64, b: f64) -> Result<f64, String> {
    if x <= 0. {
        return Ok(0.);
    }
    if x >= 1. {
        return Ok(1.);
    }
    let logarithm = log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p();
    let front = logarithm.exp();
    let value = if x < (a + 1.) / (a + b + 2.) {
        front * beta_fraction(a, b, x)? / a
    } else {
        1. - front * beta_fraction(b, a, 1. - x)? / b
    };
    let value = finite(value, "incomplete beta probability")?;
    if !(-1e-12..=1. + 1e-12).contains(&value) {
        return Err("incomplete beta probability is outside numerical range".into());
    }
    Ok(value.clamp(0., 1.))
}

fn beta_fraction(a: f64, b: f64, x: f64) -> Result<f64, String> {
    const TINY: f64 = 1e-300;
    let nonzero = |value: f64| {
        if value.abs() < TINY {
            TINY.copysign(value)
        } else {
            value
        }
    };
    let mut c = 1.;
    let mut d = 1. / nonzero(1. - (a + b) * x / (a + 1.));
    let mut fraction = d;
    for m in 1..=2000 {
        let m = m as f64;
        let even = m * (b - m) * x / ((a + 2. * m - 1.) * (a + 2. * m));
        d = 1. / nonzero(1. + even * d);
        c = nonzero(1. + even / c);
        fraction *= d * c;
        let odd = -(a + m) * (a + b + m) * x / ((a + 2. * m) * (a + 2. * m + 1.));
        d = 1. / nonzero(1. + odd * d);
        c = nonzero(1. + odd / c);
        let delta = d * c;
        fraction *= delta;
        if (delta - 1.).abs() < 2e-14 {
            return finite(fraction, "incomplete beta continued fraction");
        }
    }
    Err("incomplete beta continued fraction did not converge".into())
}
