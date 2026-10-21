//! Reproducible resampling, simulation and explicitly scoped planning formulas.
use crate::{optimize::normal_cdf, quality::log_gamma};
use std::f64::consts::PI;
fn array(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn confidence(x: f64) -> Result<f64, String> {
    if (1e-10..=1. - 1e-10).contains(&x) {
        Ok(x)
    } else {
        Err("confidence must lie in [1e-10,1-1e-10]".into())
    }
}
fn quantile_normal(p: f64) -> f64 {
    let (mut lo, mut hi) = (-40., 40.);
    for _ in 0..100 {
        let mid = (lo + hi) / 2.;
        if normal_cdf(mid) < p {
            lo = mid
        } else {
            hi = mid
        }
    }
    (lo + hi) / 2.
}
fn mean(x: &[f64]) -> f64 {
    let origin = x[0];
    origin + x.iter().map(|v| (v - origin) / x.len() as f64).sum::<f64>()
}
fn sd(x: &[f64]) -> f64 {
    let origin = x[0];
    let offset = x.iter().map(|v| (v - origin) / x.len() as f64).sum::<f64>();
    (x.iter()
        .map(|v| ((v - origin) - offset).powi(2) / (x.len() - 1) as f64)
        .sum::<f64>())
    .sqrt()
}
pub(crate) struct Rng(u64);
impl Rng {
    pub(crate) fn new(seed: f64) -> Result<Self, String> {
        if !(0. ..=9_007_199_254_740_991.).contains(&seed) || seed.fract() != 0. {
            return Err("seed must be an exact nonnegative integer".into());
        }
        Ok(Self(seed as u64))
    }
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    pub(crate) fn index(&mut self, n: usize) -> usize {
        let n = n as u64;
        let threshold = n.wrapping_neg() % n;
        loop {
            let x = self.next();
            if x >= threshold {
                return (x % n) as usize;
            }
        }
    }
    pub(crate) fn uniform(&mut self) -> f64 {
        // The upper midpoint can round to 1.0 in binary64. Preserve the open
        // interval required by logarithms, jitter and certain Bernoulli events.
        let midpoint = ((self.next() >> 11) as f64 + 0.5) / 9_007_199_254_740_992.;
        midpoint.min(f64::from_bits(1_f64.to_bits() - 1))
    }
    pub(crate) fn shuffle(&mut self, x: &mut [usize]) {
        for i in (1..x.len()).rev() {
            let j = self.index(i + 1);
            x.swap(i, j);
        }
    }
}
fn percentile(sorted: &[f64], p: f64) -> f64 {
    let pos = p * (sorted.len() - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    sorted[lo] + (sorted[hi] - sorted[lo]) * (pos - lo as f64)
}
fn bootstrap(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || data.len() < 2 || args.is_empty() || args.len() > 3 {
        return Err(
            "bootstrap_mean requires scalar rows and [replicates,confidence=.95,seed=0]".into(),
        );
    }
    let count = integer(args[0], 100000)?;
    let conf = confidence(*args.get(1).unwrap_or(&0.95))?;
    let mut rng = Rng::new(*args.get(2).unwrap_or(&0.))?;
    if count < 2 || count.saturating_mul(data.len()) > 20_000_000 {
        return Err("invalid or excessive bootstrap work".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    let mut samples = Vec::with_capacity(count);
    let origin = x[0];
    for _ in 0..count {
        samples.push(
            (0..x.len())
                .map(|_| (x[rng.index(x.len())] - origin) / x.len() as f64)
                .sum::<f64>(),
        );
    }
    samples.sort_by(f64::total_cmp);
    let alpha = (1. - conf) / 2.;
    let ci = [
        origin + percentile(&samples, alpha),
        origin + percentile(&samples, 1. - alpha),
    ];
    Ok(format!("{{\"version\":1,\"algorithm\":\"iid_percentile_bootstrap_mean\",\"estimate\":{},\"bootstrap_standard_error\":{},\"confidence_interval\":{},\"replicates\":{count},\"seed\":{}}}",mean(&x),sd(&samples),array(&ci),args.get(2).unwrap_or(&0.)))
}
fn permutation(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 2 || args.is_empty() || args.len() > 3 {
        return Err("permutation_test needs [value,group0/1] and [replicates0exact,seed=0,alternative0two/1greater/2less]".into());
    }
    let count = integer(args[0], 1_000_000)?;
    let alt = integer(*args.get(2).unwrap_or(&0.), 2)?;
    let mut rng = Rng::new(*args.get(1).unwrap_or(&0.))?;
    let n = data.len();
    let mut a = Vec::new();
    let mut b = Vec::new();
    let origin = data[0][0];
    for r in data {
        match r[1] {
            0. => a.push(r[0] - origin),
            1. => b.push(r[0] - origin),
            _ => return Err("permutation groups must be 0/1".into()),
        }
    }
    if a.is_empty() || b.is_empty() {
        return Err("both groups are required".into());
    }
    let observed = mean(&a) - mean(&b);
    let x: Vec<_> = data.iter().map(|r| r[0] - origin).collect();
    let total = x.iter().sum::<f64>();
    let na = a.len();
    let tie_tolerance =
        64. * f64::EPSILON * x.iter().map(|v| v.abs()).fold(0_f64, f64::max) * n as f64;
    let exceeds = |v: f64| match alt {
        0 => v.abs() + tie_tolerance >= observed.abs(),
        1 => v + tie_tolerance >= observed,
        _ => v - tie_tolerance <= observed,
    };
    let statistic = |ix: &[usize]| {
        let sum = ix.iter().map(|&i| x[i]).sum::<f64>();
        sum / na as f64 - (total - sum) / (n - na) as f64
    };
    let (mut assignments, mut extreme) = (0usize, 0usize);
    if count == 0 {
        let mut ways = 1.;
        for k in 1..=na.min(n - na) {
            ways *= ((n + 1 - k) as f64) / (k as f64);
            if ways > 1_000_000. || ways * n as f64 > 20_000_000. {
                return Err("exact permutation work exceeds limit; request Monte Carlo".into());
            }
        }
        let mut ix: Vec<_> = (0..na).collect();
        loop {
            assignments += 1;
            if exceeds(statistic(&ix)) {
                extreme += 1;
            }
            let pivot = (0..na).rev().find(|&i| ix[i] < n - na + i);
            match pivot {
                None => break,
                Some(i) => {
                    ix[i] += 1;
                    for j in i + 1..na {
                        ix[j] = ix[j - 1] + 1;
                    }
                }
            }
        }
    } else {
        if count.saturating_mul(n) > 20_000_000 {
            return Err("permutation work exceeds limit".into());
        }
        let mut ix: Vec<_> = (0..n).collect();
        for _ in 0..count {
            rng.shuffle(&mut ix);
            assignments += 1;
            if exceeds(statistic(&ix[..na])) {
                extreme += 1;
            }
        }
    }
    let p = if count == 0 {
        extreme as f64 / assignments as f64
    } else {
        (extreme + 1) as f64 / (assignments + 1) as f64
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"label_permutation_mean_difference\",\"statistic\":{observed},\"p_value\":{p},\"assignments\":{assignments},\"exact\":{}}}",count==0))
}
fn sample(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 3 {
        return Err("sample_rows requires [count,replacement0/1,seed]".into());
    }
    let count = integer(args[0], 100000)?;
    let replacement = integer(args[1], 1)?;
    let mut rng = Rng::new(args[2])?;
    if count == 0
        || count.saturating_mul(data[0].len()) > 1000000
        || (replacement == 0 && count > data.len())
    {
        return Err("invalid sample shape".into());
    }
    let indices = if replacement == 1 {
        (0..count)
            .map(|_| rng.index(data.len()))
            .collect::<Vec<_>>()
    } else {
        let mut indices: Vec<_> = (0..data.len()).collect();
        rng.shuffle(&mut indices);
        indices.truncate(count);
        indices
    };
    let rows = indices
        .iter()
        .map(|&i| array(&data[i]))
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!(
        "{{\"version\":1,\"indices\":{},\"rows\":[{rows}]}}",
        array(&indices.iter().map(|&v| v as f64).collect::<Vec<_>>())
    ))
}
fn simulate(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 2 || args.len() != 2 {
        return Err("simulate_normal requires [mean,sd] rows and [count_per_row,seed]".into());
    }
    let count = integer(args[0], 1000000)?;
    let mut rng = Rng::new(args[1])?;
    if count == 0 || count.saturating_mul(data.len()) > 1000000 || data.iter().any(|r| r[1] < 0.) {
        return Err("invalid simulation size or standard deviation".into());
    }
    let mut x = Vec::new();
    for r in data {
        for _ in 0..count {
            let z = (-2. * rng.uniform().ln()).sqrt() * (2. * PI * rng.uniform()).cos();
            x.push(r[0] + r[1] * z);
        }
    }
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"box_muller_normal\",\"shape\":[{},{}],\"samples\":{}}}",
        data.len(),
        count,
        array(&x)
    ))
}
fn planning(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if operation == "power_z" {
        if data[0].len() != 2 || args.is_empty() || args.len() > 3 {
            return Err("power_z requires [effect,known_sigma] and [n,alpha=.05,alternative0two/1greater/2less]".into());
        }
        let n = integer(args[0], 1_000_000_000)?;
        let alpha = *args.get(1).unwrap_or(&0.05);
        confidence(alpha)?;
        if n == 0 || alpha >= 0.5 || data.iter().any(|r| r[1] <= 0.) {
            return Err("invalid power design".into());
        }
        let alt = integer(*args.get(2).unwrap_or(&0.), 2)?;
        let z = quantile_normal(1. - if alt == 0 { alpha / 2. } else { alpha });
        let power: Vec<_> = data
            .iter()
            .map(|r| {
                let shift = r[0] / r[1] * (n as f64).sqrt();
                match alt {
                    0 => normal_cdf(-z - shift) + normal_cdf(shift - z),
                    1 => normal_cdf(shift - z),
                    _ => normal_cdf(-shift - z),
                }
            })
            .collect();
        return Ok(format!("{{\"version\":1,\"algorithm\":\"one_sample_known_sigma_normal_power\",\"power\":{},\"n\":{n}}}",array(&power)));
    }
    if operation == "reliability_demonstration" {
        if data[0].len() != 1 || args.len() > 1 {
            return Err(
                "reliability demonstration requires reliability rows and [confidence=.95]".into(),
            );
        }
        let conf = confidence(*args.first().unwrap_or(&0.95))?;
        let mut sizes = Vec::new();
        for r in data {
            if !(0. < r[0] && r[0] < 1.) {
                return Err("reliability target must lie in (0,1)".into());
            }
            let n = ((1. - conf).ln() / r[0].ln()).ceil();
            if n > 1e12 {
                return Err("sample size exceeds limit".into());
            }
            sizes.push(n.max(1.));
        }
        return Ok(format!("{{\"version\":1,\"algorithm\":\"independent_zero_failure_binomial_demonstration\",\"sample_sizes\":{}}}",array(&sizes)));
    }
    if data[0].len() != 1 || args.is_empty() || args.len() > 2 || args[0] <= 0. {
        return Err(
            "margin planning requires scalar sigma/proportion rows and [margin,confidence=.95]"
                .into(),
        );
    }
    let conf = confidence(*args.get(1).unwrap_or(&0.95))?;
    let z = quantile_normal(0.5 + conf / 2.);
    let mut sizes = Vec::new();
    for r in data {
        let ratio = if operation == "sample_size_mean" {
            if r[0] <= 0. {
                return Err("sigma must be positive".into());
            }
            r[0] / args[0]
        } else {
            if !(0. ..=1.).contains(&r[0]) {
                return Err("p must lie in [0,1]".into());
            }
            (r[0] * (1. - r[0])).sqrt() / args[0]
        };
        let size = (z * ratio).powi(2);
        if !size.is_finite() || size > 1e12 {
            return Err("sample size exceeds limit".into());
        }
        sizes.push(size.ceil().max(1.));
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"normal_approximation_margin_planning\",\"sample_sizes\":{}}}",array(&sizes)))
}
fn effect(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if op == "cohen_d" {
        if data[0].len() != 2 || !args.is_empty() {
            return Err("cohen_d requires [value,group0/1], []".into());
        }
        let origin = data[0][0];
        let mut a = Vec::new();
        let mut b = Vec::new();
        for r in data {
            match r[1] {
                0. => a.push(r[0] - origin),
                1. => b.push(r[0] - origin),
                _ => return Err("groups must be 0/1".into()),
            }
        }
        if a.len() < 2 || b.len() < 2 {
            return Err("both groups require at least two observations".into());
        }
        let df = (a.len() + b.len() - 2) as f64;
        let pooled =
            (((a.len() - 1) as f64 * sd(&a).powi(2) + (b.len() - 1) as f64 * sd(&b).powi(2)) / df)
                .sqrt();
        if pooled == 0. {
            return Err("pooled variance is zero".into());
        }
        let d = (mean(&a) - mean(&b)) / pooled;
        let correction =
            (log_gamma(df / 2.) - 0.5 * (df / 2.).ln() - log_gamma((df - 1.) / 2.)).exp();
        return Ok(format!("{{\"version\":1,\"cohen_d\":{d},\"hedges_g\":{},\"pooled_sd\":{pooled},\"residual_df\":{df}}}",d*correction));
    }
    if data.len() != 2 || data[0].len() != 2 || args.len() > 2 {
        return Err(
            "risk effect requires a 2x2 table and [confidence=.95,continuity_correction=0]".into(),
        );
    }
    let conf = confidence(*args.first().unwrap_or(&0.95))?;
    let correction = *args.get(1).unwrap_or(&0.);
    if correction < 0. {
        return Err("correction must be nonnegative".into());
    }
    for v in data.iter().flatten() {
        integer(*v, 1_000_000_000)?;
    }
    let a = data[0][0] + correction;
    let b = data[0][1] + correction;
    let c = data[1][0] + correction;
    let d = data[1][1] + correction;
    if [a, b, c, d].iter().any(|v| *v <= 0.) {
        return Err("positive cells required; specify correction explicitly for zeros".into());
    }
    let (estimate, se) = if op == "odds_ratio" {
        (a * d / (b * c), (1. / a + 1. / b + 1. / c + 1. / d).sqrt())
    } else {
        (
            (a / (a + b)) / (c / (c + d)),
            (1. / a - 1. / (a + b) + 1. / c - 1. / (c + d)).sqrt(),
        )
    };
    let z = quantile_normal(0.5 + conf / 2.);
    let ci = [
        (estimate.ln() - z * se).exp(),
        (estimate.ln() + z * se).exp(),
    ];
    Ok(format!("{{\"version\":1,\"algorithm\":\"log_scale_wald_{op}\",\"estimate\":{estimate},\"log_standard_error\":{se},\"confidence_interval\":{}}}",array(&ci)))
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 100000
        || data.len().saturating_mul(width) > 1000000
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite() || v.abs() > 1e100))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("resampling requires bounded finite rectangular data".into());
    }
    let output = match operation {
        "bootstrap_mean" => bootstrap(data, args),
        "permutation_test" => permutation(data, args),
        "sample_rows" => sample(data, args),
        "simulate_normal" => simulate(data, args),
        "power_z" | "sample_size_mean" | "sample_size_proportion" | "reliability_demonstration" => {
            planning(operation, data, args)
        }
        "odds_ratio" | "relative_risk" | "cohen_d" => effect(operation, data, args),
        _ => Err("unknown sampling operation".into()),
    }?;
    if output
        .split([':', '[', ']', ',', '}', ' '])
        .any(|token| matches!(token, "inf" | "-inf" | "NaN"))
    {
        Err("resampling calculation exceeds numeric range".into())
    } else {
        Ok(output)
    }
}
