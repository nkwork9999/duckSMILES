//! Control-chart calculations and explicitly scoped measurement-system models.
//! Rows must be in time order. Estimated limits use the supplied Phase-I sample;
//! monitoring future data should use saved, externally fixed center/sigma.
use std::collections::{BTreeMap, BTreeSet};
use std::f64::consts::PI;

fn arr(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn int(v: f64, max: usize) -> Result<usize, String> {
    if v.is_finite() && v >= 0. && v <= max as f64 && v.fract() == 0. {
        Ok(v as usize)
    } else {
        Err("expected bounded integer".into())
    }
}
fn moments(x: &[f64]) -> (f64, f64) {
    let origin = x[0];
    let offset = x.iter().map(|v| (v - origin) / x.len() as f64).sum::<f64>();
    let variance = if x.len() > 1 {
        x.iter()
            .map(|v| ((v - origin) - offset).powi(2) / (x.len() - 1) as f64)
            .sum::<f64>()
    } else {
        0.
    };
    (origin + offset, variance.sqrt())
}
fn signals(x: &[f64], lo: &[f64], hi: &[f64]) -> Vec<f64> {
    x.iter()
        .zip(lo)
        .zip(hi)
        .map(|((x, l), h)| if x < l || x > h { 1. } else { 0. })
        .collect()
}
fn chart(
    name: &str,
    x: &[f64],
    center: f64,
    sigma: f64,
    lo: &[f64],
    hi: &[f64],
    extra: &str,
) -> String {
    format!("{{\"version\":1,\"algorithm\":\"{name}\",\"center\":{center},\"sigma\":{sigma},\"statistics\":{},\"lower\":{},\"upper\":{},\"signals\":{}{extra}}}",arr(x),arr(lo),arr(hi),arr(&signals(x,lo,hi)))
}
fn ir(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || data.len() < 2 || !matches!(args.len(), 0 | 2) {
        return Err("IR requires ordered scalar rows and optional [center,sigma]".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    let mr: Vec<_> = x.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let mrbar = moments(&mr).0;
    let center = if args.is_empty() {
        moments(&x).0
    } else {
        args[0]
    };
    let sigma = if args.is_empty() {
        mrbar * PI.sqrt() / 2.
    } else {
        args[1]
    };
    if sigma <= 0. {
        return Err("positive within-process sigma required".into());
    }
    let lo = vec![center - 3. * sigma; x.len()];
    let hi = vec![center + 3. * sigma; x.len()];
    let d4 = 1. + 3. * (PI / 2. - 1.).sqrt();
    Ok(chart(
        "individuals_moving_range",
        &x,
        center,
        sigma,
        &lo,
        &hi,
        &format!(
            ",\"moving_ranges\":{},\"moving_range_center\":{mrbar},\"moving_range_upper\":{}",
            arr(&mr),
            d4 * mrbar
        ),
    ))
}
pub(crate) fn log_gamma(z: f64) -> f64 {
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
    let z = z - 1.;
    let mut x = 0.9999999999998099;
    for (i, c) in c.iter().enumerate() {
        x += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * (2. * PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}
fn subgroup(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let n = data[0].len();
    if data.len() < 2 || !(2..=1000).contains(&n) || !args.is_empty() {
        return Err(
            "subgroup chart requires >=2 equal subgroups of size 2..1000, no parameters".into(),
        );
    }
    let means: Vec<_> = data.iter().map(|r| moments(r).0).collect();
    let center = moments(&means).0;
    let spread: Vec<_> = data
        .iter()
        .map(|r| {
            if operation == "xbar_r_chart" {
                r.iter().copied().reduce(f64::max).unwrap()
                    - r.iter().copied().reduce(f64::min).unwrap()
            } else {
                moments(r).1
            }
        })
        .collect();
    let spreadbar = moments(&spread).0;
    let (sigma, spread_sigma) = if operation == "xbar_r_chart" {
        if n > 10 {
            return Err("Xbar-R currently supports subgroup sizes 2..10".into());
        }
        let d2 = [
            1.128, 1.693, 2.059, 2.326, 2.534, 2.704, 2.847, 2.970, 3.078,
        ][n - 2];
        let d3 = [
            0.853, 0.888, 0.880, 0.864, 0.848, 0.833, 0.820, 0.808, 0.797,
        ][n - 2];
        let sigma = spreadbar / d2;
        (sigma, d3 * sigma)
    } else {
        let c4 = (2. / (n - 1) as f64).sqrt()
            * (log_gamma(n as f64 / 2.) - log_gamma((n - 1) as f64 / 2.)).exp();
        let sigma = spreadbar / c4;
        (sigma, sigma * (1. - c4 * c4).max(0.).sqrt())
    };
    if sigma <= 0. {
        return Err("subgroups have no estimable within-process variation".into());
    }
    let delta = 3. * sigma / (n as f64).sqrt();
    let lo = vec![center - delta; data.len()];
    let hi = vec![center + delta; data.len()];
    Ok(chart(operation,&means,center,sigma,&lo,&hi,&format!(",\"subgroup_size\":{n},\"spread\":{},\"spread_center\":{spreadbar},\"spread_lower\":{},\"spread_upper\":{}",arr(&spread),(spreadbar-3.*spread_sigma).max(0.),spreadbar+3.*spread_sigma)))
}
fn attribute(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let count_only = op == "c_chart";
    let width = if count_only { 1 } else { 2 };
    if data[0].len() != width || !args.is_empty() {
        return Err("attribute rows are [count] for C, [count,exposure] for P/NP/U".into());
    }
    for r in data {
        int(r[0], 1_000_000_000)?;
        if !count_only && (r[1] <= 0. || (op != "u_chart" && (r[0] > r[1] || r[1].fract() != 0.))) {
            return Err("invalid count or sample size".into());
        }
    }
    if op == "np_chart" && data.iter().any(|r| r[1] != data[0][1]) {
        return Err("NP chart requires constant sample size".into());
    }
    let total_count = data.iter().map(|r| r[0]).sum::<f64>();
    let exposure = if count_only {
        data.len() as f64
    } else {
        data.iter().map(|r| r[1]).sum::<f64>()
    };
    let rate = total_count / exposure;
    let center = if op == "np_chart" {
        rate * data[0][1]
    } else {
        rate
    };
    let mut x = Vec::new();
    let mut lo = Vec::new();
    let mut hi = Vec::new();
    for r in data {
        let (stat, sd) = match op {
            "p_chart" => (r[0] / r[1], (rate * (1. - rate) / r[1]).sqrt()),
            "np_chart" => (r[0], (r[1] * rate * (1. - rate)).sqrt()),
            "c_chart" => (r[0], rate.sqrt()),
            _ => (r[0] / r[1], (rate / r[1]).sqrt()),
        };
        x.push(stat);
        lo.push((center - 3. * sd).max(0.));
        hi.push(if op == "p_chart" {
            (center + 3. * sd).min(1.)
        } else if op == "np_chart" {
            (center + 3. * sd).min(r[1])
        } else {
            center + 3. * sd
        });
    }
    Ok(chart(
        op,
        &x,
        center,
        0.,
        &lo,
        &hi,
        ",\"sigma_policy\":\"point_specific_count_model\"",
    ))
}

fn laney(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 2 || data.len() < 2 || !args.is_empty() {
        return Err("Laney charts require >=2 [count,exposure] rows".into());
    }
    let is_p = op == "laney_p_chart";
    for row in data {
        int(row[0], 1_000_000_000)?;
        if row[1] <= 0. || (is_p && (row[1].fract() != 0. || row[0] > row[1])) {
            return Err("invalid Laney count/exposure".into());
        }
    }
    let center = data.iter().map(|r| r[0]).sum::<f64>() / data.iter().map(|r| r[1]).sum::<f64>();
    if center <= 0. || (is_p && center >= 1.) {
        return Err("Laney variation is not identified".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0] / r[1]).collect();
    let sd: Vec<_> = data
        .iter()
        .map(|r| (center * if is_p { 1. - center } else { 1. } / r[1]).sqrt())
        .collect();
    let z: Vec<_> = x.iter().zip(&sd).map(|(x, sd)| (x - center) / sd).collect();
    let mr: Vec<_> = z.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    let sigma_z = moments(&mr).0 * PI.sqrt() / 2.;
    if sigma_z <= 0. {
        return Err("Laney standardized moving-range variance is zero".into());
    }
    let lo: Vec<_> = sd
        .iter()
        .map(|s| (center - 3. * sigma_z * s).max(0.))
        .collect();
    let hi: Vec<_> = sd
        .iter()
        .map(|s| {
            if is_p {
                (center + 3. * sigma_z * s).min(1.)
            } else {
                center + 3. * sigma_z * s
            }
        })
        .collect();
    Ok(chart(
        op,
        &x,
        center,
        0.,
        &lo,
        &hi,
        &format!(",\"sigma_z\":{sigma_z}"),
    ))
}

fn levey(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() != 2 || args[1] <= 0. {
        return Err("Levey-Jennings requires [reference_mean,reference_sigma>0]".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    let lo = vec![args[0] - 3. * args[1]; data.len()];
    let hi = vec![args[0] + 3. * args[1]; data.len()];
    Ok(chart("levey_jennings", &x, args[0], args[1], &lo, &hi, ""))
}
fn sequential(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() != 4 || args[1] <= 0. || args[3] <= 0. {
        return Err("sequential chart requires [tuning,limit,known_center,known_sigma>0]".into());
    }
    let center = args[2];
    let sigma = args[3];
    let mut x = Vec::new();
    let mut lo = Vec::new();
    let mut hi = Vec::new();
    match op {
        "ewma_chart" => {
            let lambda = args[0];
            if !(0. < lambda && lambda <= 1.) {
                return Err("EWMA lambda must lie in (0,1]".into());
            }
            let mut z = center;
            for (i, r) in data.iter().enumerate() {
                z = lambda * r[0] + (1. - lambda) * z;
                let delta = args[1]
                    * sigma
                    * (lambda / (2. - lambda)
                        * (-((2. * (i + 1) as f64) * (-lambda).ln_1p()).exp_m1()))
                    .sqrt();
                x.push(z);
                lo.push(center - delta);
                hi.push(center + delta);
            }
        }
        "moving_average_chart" => {
            let window = int(args[0], 100000)?;
            if window == 0 {
                return Err("window must be positive".into());
            }
            let origin = data[0][0];
            let mut sum = 0.;
            for (i, r) in data.iter().enumerate() {
                sum += r[0] - origin;
                if i >= window {
                    sum -= data[i - window][0] - origin;
                }
                let n = (i + 1).min(window);
                x.push(origin + sum / n as f64);
                let delta = args[1] * sigma / (n as f64).sqrt();
                lo.push(center - delta);
                hi.push(center + delta);
            }
        }
        "cusum_chart" => {
            if args[0] < 0. {
                return Err("CUSUM reference value must be nonnegative".into());
            }
            let mut plus = 0_f64;
            let mut minus = 0_f64;
            let mut negative = Vec::new();
            let mut flags = Vec::new();
            for r in data {
                let z = (r[0] - center) / sigma;
                plus = (plus + z - args[0]).max(0.);
                minus = (minus - z - args[0]).max(0.);
                x.push(plus);
                negative.push(minus);
                flags.push(if plus > args[1] || minus > args[1] {
                    1.
                } else {
                    0.
                });
            }
            return Ok(format!("{{\"version\":1,\"algorithm\":\"standardized_two_sided_cusum\",\"positive\":{},\"negative\":{},\"signals\":{},\"reference\":{},\"decision_limit\":{}}}",arr(&x),arr(&negative),arr(&flags),args[0],args[1]));
        }
        _ => unreachable!(),
    }
    Ok(chart(op, &x, center, sigma, &lo, &hi, ""))
}
fn western(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() != 2 || args[1] <= 0. {
        return Err("Western Electric requires [center,sigma>0]".into());
    }
    let z: Vec<_> = data.iter().map(|r| (r[0] - args[0]) / args[1]).collect();
    let mut masks = Vec::new();
    for i in 0..z.len() {
        let mut mask = if z[i].abs() > 3. { 1 } else { 0 };
        for (window, threshold, needed, bit) in [(3, 2., 2, 2), (5, 1., 4, 4), (8, 0., 8, 8)] {
            if i + 1 >= window {
                let w = &z[i + 1 - window..=i];
                if w.iter().filter(|&&v| v > threshold).count() >= needed
                    || w.iter().filter(|&&v| v < -threshold).count() >= needed
                {
                    mask |= bit;
                }
            }
        }
        masks.push(mask as f64);
    }
    Ok(format!("{{\"version\":1,\"rule_masks\":{},\"rule_bits\":{{\"beyond_3sigma\":1,\"two_of_three_beyond_2sigma\":2,\"four_of_five_beyond_1sigma\":4,\"eight_same_side\":8}}}}",arr(&masks)))
}
fn gage(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 3 || !args.is_empty() {
        return Err("gage_rr requires [value,part_id,operator_id] rows, no parameters".into());
    }
    let mut cells = BTreeMap::<(usize, usize), Vec<f64>>::new();
    let mut parts = BTreeSet::new();
    let mut operators = BTreeSet::new();
    for r in data {
        let p = int(r[1], 1_000_000_000)?;
        let o = int(r[2], 1_000_000_000)?;
        parts.insert(p);
        operators.insert(o);
        cells.entry((p, o)).or_default().push(r[0]);
    }
    let a = parts.len();
    let b = operators.len();
    let repeats = cells.values().next().unwrap().len();
    if a < 2
        || b < 2
        || repeats < 2
        || a.saturating_mul(b) > 10000
        || cells.len() != a * b
        || cells.values().any(|v| v.len() != repeats)
    {
        return Err(
            "gage_rr needs a complete balanced crossed design with >=2 parts/operators/repeats"
                .into(),
        );
    }
    let origin = data[0][0];
    let all: Vec<_> = data.iter().map(|r| r[0] - origin).collect();
    let grand = moments(&all).0;
    let mut partmeans = BTreeMap::new();
    let mut opmeans = BTreeMap::new();
    let mut cellmeans = BTreeMap::new();
    let mut ss_error = 0.;
    for (key, values) in &cells {
        let shifted: Vec<_> = values.iter().map(|v| v - origin).collect();
        let mean = moments(&shifted).0;
        ss_error += shifted.iter().map(|v| (v - mean).powi(2)).sum::<f64>();
        cellmeans.insert(*key, mean);
    }
    for &p in &parts {
        partmeans.insert(
            p,
            operators
                .iter()
                .map(|&o| cellmeans[&(p, o)] / b as f64)
                .sum::<f64>(),
        );
    }
    for &o in &operators {
        opmeans.insert(
            o,
            parts
                .iter()
                .map(|&p| cellmeans[&(p, o)] / a as f64)
                .sum::<f64>(),
        );
    }
    let ms_part = (b * repeats) as f64
        * partmeans.values().map(|v| (v - grand).powi(2)).sum::<f64>()
        / (a - 1) as f64;
    let ms_op = (a * repeats) as f64 * opmeans.values().map(|v| (v - grand).powi(2)).sum::<f64>()
        / (b - 1) as f64;
    let ms_interaction = repeats as f64
        * cellmeans
            .iter()
            .map(|(&(p, o), v)| (v - partmeans[&p] - opmeans[&o] + grand).powi(2))
            .sum::<f64>()
        / ((a - 1) * (b - 1)) as f64;
    let repeatability = ss_error / (a * b * (repeats - 1)) as f64;
    let interaction = ((ms_interaction - repeatability) / repeats as f64).max(0.);
    let operator = ((ms_op - ms_interaction) / (a * repeats) as f64).max(0.);
    let part = ((ms_part - ms_interaction) / (b * repeats) as f64).max(0.);
    let rr = repeatability + interaction + operator;
    let total = rr + part;
    let percent = if total > 0. {
        (100. * rr / total).to_string()
    } else {
        "null".into()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"balanced_crossed_random_effects_anova_no_pooling\",\"repeatability_variance\":{repeatability},\"interaction_variance\":{interaction},\"operator_variance\":{operator},\"part_variance\":{part},\"gage_variance\":{rr},\"total_variance\":{total},\"gage_percent_contribution\":{percent},\"n_parts\":{a},\"n_operators\":{b},\"n_repeats\":{repeats}}}"))
}
fn type1(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || data.len() < 2 || args.len() != 2 || args[1] <= 0. {
        return Err("type1_gage requires [reference,tolerance_width>0]".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[0]).collect();
    let (mean, sd) = moments(&x);
    if sd <= 0. {
        return Err("gage variance undefined".into());
    }
    let bias = mean - args[0];
    let cg = 0.2 * args[1] / (6. * sd);
    let cgk = (0.1 * args[1] - bias.abs()) / (3. * sd);
    Ok(format!("{{\"version\":1,\"mean\":{mean},\"sd\":{sd},\"bias\":{bias},\"cg\":{cg},\"cgk\":{cgk},\"tolerance_fraction\":0.2}}"))
}
fn kappa(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let k = data.len();
    if k < 2 || k != data[0].len() || !args.is_empty() {
        return Err("kappa requires a square >=2 contingency count matrix".into());
    }
    for v in data.iter().flatten() {
        int(*v, 1_000_000_000)?;
    }
    let n = data.iter().flatten().sum::<f64>();
    if n == 0. {
        return Err("empty contingency table".into());
    }
    let po = (0..k).map(|i| data[i][i] / n).sum::<f64>();
    let pe = (0..k)
        .map(|i| data[i].iter().sum::<f64>() / n * data.iter().map(|r| r[i]).sum::<f64>() / n)
        .sum::<f64>();
    if pe >= 1. {
        return Err("chance agreement is one; kappa undefined".into());
    }
    let kappa = (po - pe) / (1. - pe);
    Ok(format!(
        "{{\"version\":1,\"observed_agreement\":{po},\"chance_agreement\":{pe},\"kappa\":{kappa}}}"
    ))
}
fn oc(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() != 2 {
        return Err("acceptance_oc requires probabilities and [n,acceptance_count]".into());
    }
    let n = int(args[0], 100000)?;
    let c = int(args[1], n)?;
    if n == 0 || data.len().saturating_mul(c + 1) > 20_000_000 {
        return Err("invalid or oversized sampling plan".into());
    }
    let mut probs = Vec::new();
    for r in data {
        let p = r[0];
        if !(0. ..=1.).contains(&p) {
            return Err("probability outside [0,1]".into());
        }
        let sum = if p == 0. || c == n {
            1.
        } else if p == 1. {
            0.
        } else {
            let mut mass = n as f64 * (-p).ln_1p();
            let mut sum = mass.exp();
            for j in 0..c {
                mass += ((n - j) as f64).ln() - ((j + 1) as f64).ln() + p.ln() - (-p).ln_1p();
                sum += mass.exp();
            }
            sum.clamp(0., 1.)
        };
        probs.push(sum);
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"single_binomial_acceptance_plan\",\"n\":{n},\"acceptance_count\":{c},\"acceptance_probability\":{}}}",arr(&probs)))
}
pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 100000
        || data.len().saturating_mul(width) > 1000000
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite() || v.abs() > 1e100))
        || parameters.iter().any(|v| !v.is_finite() || v.abs() > 1e100)
    {
        return Err("quality requires bounded finite rectangular data".into());
    }
    let result = match operation {
        "ir_chart" => ir(data, parameters),
        "xbar_r_chart" | "xbar_s_chart" => subgroup(operation, data, parameters),
        "p_chart" | "np_chart" | "c_chart" | "u_chart" => attribute(operation, data, parameters),
        "laney_p_chart" | "laney_u_chart" => laney(operation, data, parameters),
        "levey_jennings" => levey(data, parameters),
        "ewma_chart" | "cusum_chart" | "moving_average_chart" => {
            sequential(operation, data, parameters)
        }
        "western_electric" => western(data, parameters),
        "gage_rr" => gage(data, parameters),
        "type1_gage" => type1(data, parameters),
        "agreement_kappa" => kappa(data, parameters),
        "acceptance_oc" => oc(data, parameters),
        _ => Err("unknown quality operation".into()),
    }?;
    // All strings in this schema are fixed labels. A bare nonfinite numeric
    // token therefore means a derived calculation overflowed, never valid JSON.
    if result
        .split([':', '[', ']', ',', '}', ' '])
        .any(|token| matches!(token, "inf" | "-inf" | "NaN"))
    {
        Err("quality calculation exceeds finite numeric range".into())
    } else {
        Ok(result)
    }
}
