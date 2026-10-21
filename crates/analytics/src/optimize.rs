//! Candidate-pool Bayesian optimization with an explicit RBF Gaussian process.
//! Hyperparameters are supplied by the caller; no hidden optimizer or fitting
//! outside the training rows is performed. Posterior variances are latent.
use std::f64::consts::PI;

fn array(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|v| v.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn integer(value: f64, max: usize) -> Result<usize, String> {
    if value.is_finite() && value >= 0. && value.fract() == 0. && value <= max as f64 {
        Ok(value as usize)
    } else {
        Err("expected a bounded nonnegative integer".into())
    }
}

fn validate(data: &[Vec<f64>], parameters: &[f64]) -> Result<usize, String> {
    let p = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > 100_000
        || p == 0
        || p > 130
        || data.len().saturating_mul(p) > 1_000_000
        || data
            .iter()
            .any(|r| r.len() != p || r.iter().any(|v| !v.is_finite()))
        || parameters.iter().any(|v| !v.is_finite())
    {
        return Err("expected a bounded nonempty rectangular finite matrix".into());
    }
    Ok(p)
}

// Phi(x) = Q(1/2, x²/2)/2 in the lower tail. A convergent gamma
// series near zero and a continued fraction in the tail avoid subtracting an
// erf rounded to 1, which otherwise destroys small expected improvements.
pub(crate) fn normal_cdf(x: f64) -> f64 {
    if x == 0. {
        return 0.5;
    }
    if x.abs() >= 38. {
        return if x > 0. { 1. } else { 0. };
    }
    let z = x * x / 2.;
    let factor = (-z + 0.5 * z.ln() - 0.5723649429247001).exp();
    let lower_tail = if z < 1.5 {
        let mut term = 2.;
        let mut sum = term;
        for i in 1..=500 {
            term *= z / (0.5 + i as f64);
            sum += term;
            if term.abs() <= sum.abs() * f64::EPSILON {
                break;
            }
        }
        0.5 * (1. - sum * factor)
    } else {
        let mut b = z + 0.5;
        let mut c = 1e300;
        let mut d = 1. / b;
        let mut h = d;
        for i in 1..=500 {
            let a = -(i as f64) * (i as f64 - 0.5);
            b += 2.;
            d = a * d + b;
            if d.abs() < 1e-300 {
                d = 1e-300;
            }
            c = b + a / c;
            if c.abs() < 1e-300 {
                c = 1e-300;
            }
            d = 1. / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() <= 4. * f64::EPSILON {
                break;
            }
        }
        0.5 * factor * h
    };
    if x < 0. {
        lower_tail
    } else {
        1. - lower_tail
    }
}

fn improvement(mean: f64, sigma: f64, best: f64, minimize: bool, xi: f64) -> f64 {
    let delta = if minimize { best - mean } else { mean - best } - xi;
    if sigma == 0. {
        return delta.max(0.);
    }
    let z = delta / sigma;
    (delta * normal_cdf(z) + sigma * (-0.5 * z * z).exp() / (2. * PI).sqrt()).max(0.)
}

fn acquisition(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 3
        || !matches!(args[1], 0. | 1.)
        || args[2] < 0.
        || data[0].len() != 2
        || data.iter().any(|r| r[1] < 0.)
    {
        return Err("acquisition requires [mean,sigma] rows and [best,minimize,xi>=0]".into());
    }
    let ei: Vec<_> = data
        .iter()
        .map(|r| improvement(r[0], r[1], args[0], args[1] == 1., args[2]))
        .collect();
    if ei.iter().any(|v| !v.is_finite()) {
        return Err("acquisition numerical range exceeded".into());
    }
    Ok(format!(
        "{{\"version\":1,\"expected_improvement\":{}}}",
        array(&ei)
    ))
}

fn forward(l: &[Vec<f64>], b: &[f64]) -> Vec<f64> {
    let mut x = vec![0.; b.len()];
    for i in 0..b.len() {
        x[i] = (b[i] - (0..i).map(|j| l[i][j] * x[j]).sum::<f64>()) / l[i][i];
    }
    x
}

fn gp(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 6
        || args[1] <= 0.
        || args[2] <= 0.
        || args[3] < 0.
        || !matches!(args[4], 0. | 1.)
        || args[5] < 0.
        || data[0].len() < 2
    {
        return Err("GP requires [n_train,length_scale>0,signal_variance>0,noise_variance>=0,minimize,xi>=0]".into());
    }
    let n = integer(args[0], 512)?;
    let p = data[0].len() - 1;
    if n == 0 || n >= data.len() || data.len() - n > 10_000 {
        return Err("GP needs 1..512 training rows followed by 1..10000 candidates".into());
    }
    if n.saturating_mul(data.len()).saturating_mul(p) > 100_000_000
        || (data.len() - n).saturating_mul(n).saturating_mul(n) > 100_000_000
    {
        return Err("GP numerical work limit exceeded".into());
    }
    // All preprocessing is trained exclusively on the observed prefix.
    let mut means = vec![0.; p];
    let origins = data[0][..p].to_vec();
    let mut offsets = vec![0.; p];
    let mut scales = vec![1.; p];
    for j in 0..p {
        offsets[j] = data[..n]
            .iter()
            .map(|r| (r[j] - origins[j]) / n as f64)
            .sum();
        means[j] = origins[j] + offsets[j];
        let v = data[..n]
            .iter()
            .map(|r| ((r[j] - origins[j]) - offsets[j]).powi(2) / n as f64)
            .sum::<f64>();
        if v > 0. {
            scales[j] = v.sqrt();
        }
    }
    let x: Vec<Vec<_>> = data
        .iter()
        .map(|r| {
            (0..p)
                .map(|j| ((r[j] - origins[j]) - offsets[j]) / scales[j])
                .collect()
        })
        .collect();
    let kernel = |a: &[f64], b: &[f64]| -> f64 {
        let distance = a
            .iter()
            .zip(b)
            .map(|(a, b)| ((a - b) / args[1]).powi(2))
            .sum::<f64>();
        args[2] * (-0.5 * distance).exp()
    };
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let mut value = kernel(&x[i], &x[j]);
            if i == j {
                value += args[3];
            }
            value -= (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if !value.is_finite() || value <= f64::EPSILON * args[2] {
                    return Err("GP covariance is singular; use positive observation noise".into());
                }
                l[i][j] = value.sqrt();
            } else {
                l[i][j] = value / l[j][j];
            }
        }
    }
    let y_mean = data[..n].iter().map(|r| r[p] / n as f64).sum::<f64>();
    let y: Vec<_> = data[..n].iter().map(|r| r[p] - y_mean).collect();
    let z = forward(&l, &y);
    let mut alpha = vec![0.; n];
    for i in (0..n).rev() {
        alpha[i] = (z[i] - ((i + 1)..n).map(|j| l[j][i] * alpha[j]).sum::<f64>()) / l[i][i];
    }
    let minimize = args[4] == 1.;
    let best = data[..n]
        .iter()
        .map(|r| r[p])
        .reduce(|a, b| if minimize { a.min(b) } else { a.max(b) })
        .unwrap();
    let mut predictions = Vec::new();
    let mut variances = Vec::new();
    let mut ei = Vec::new();
    for row in &x[n..] {
        let k: Vec<_> = x[..n].iter().map(|r| kernel(row, r)).collect();
        let mean = y_mean + k.iter().zip(&alpha).map(|(a, b)| a * b).sum::<f64>();
        let v = forward(&l, &k);
        let variance = (args[2] - v.iter().map(|v| v * v).sum::<f64>()).max(0.);
        predictions.push(mean);
        variances.push(variance);
        ei.push(improvement(mean, variance.sqrt(), best, minimize, args[5]));
    }
    let mut order: Vec<_> = (0..ei.len()).collect();
    order.sort_by(|&a, &b| ei[b].total_cmp(&ei[a]).then(a.cmp(&b)));
    let indices: Vec<_> = order.iter().map(|&v| v as f64).collect();
    let result=format!("{{\"version\":1,\"algorithm\":\"rbf_gaussian_process\",\"n_train\":{n},\"predictions\":{},\"variances\":{},\"expected_improvement\":{},\"ranked_candidate_indices\":{},\"feature_means\":{},\"feature_scales\":{},\"target_mean\":{y_mean}}}",array(&predictions),array(&variances),array(&ei),array(&indices),array(&means),array(&scales));
    if predictions
        .iter()
        .chain(&variances)
        .chain(&ei)
        .chain(&means)
        .chain(&scales)
        .any(|v| !v.is_finite())
        || !y_mean.is_finite()
    {
        Err("GP numerical range exceeded".into())
    } else {
        Ok(result)
    }
}

fn diversity(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 2 {
        return Err("diversity_select requires [n_select,n_existing]".into());
    }
    if data.iter().flatten().any(|v| v.abs() > 1e150) {
        return Err("diversity distance exceeds numerical range".into());
    }
    let count = integer(args[0], 10_000)?;
    let existing = integer(args[1], data.len())?;
    if count == 0
        || existing >= data.len()
        || count > data.len() - existing
        || (count + existing)
            .saturating_mul(data.len())
            .saturating_mul(data[0].len())
            > 20_000_000
    {
        return Err("invalid or oversized selection request".into());
    }
    let distance = |a: usize, b: usize| {
        data[a]
            .iter()
            .zip(&data[b])
            .map(|(a, b)| (a - b).powi(2))
            .sum::<f64>()
    };
    let mut nearest = vec![f64::INFINITY; data.len()];
    for (i, minimum) in nearest.iter_mut().enumerate().skip(existing) {
        for j in 0..existing {
            *minimum = minimum.min(distance(i, j));
        }
    }
    let mut selected: Vec<usize> = Vec::new();
    let mut used = vec![false; data.len()];
    for _ in 0..count {
        let next = (existing..data.len())
            .filter(|&i| !used[i])
            .max_by(|&a, &b| nearest[a].total_cmp(&nearest[b]).then(b.cmp(&a)))
            .unwrap();
        if nearest[next] <= 0. {
            return Err("not enough distinct unobserved candidates".into());
        }
        selected.push(next);
        used[next] = true;
        for (i, minimum) in nearest.iter_mut().enumerate().skip(existing) {
            *minimum = minimum.min(distance(i, next));
        }
    }
    Ok(format!(
        "{{\"version\":1,\"selected_indices\":{}}}",
        array(&selected.iter().map(|v| *v as f64).collect::<Vec<_>>())
    ))
}

pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    validate(data, parameters)?;
    match operation {
        "acquisition" => acquisition(data, parameters),
        "gp_predict" => gp(data, parameters),
        "diversity_select" => diversity(data, parameters),
        _ => Err("unknown optimization operation".into()),
    }
}
