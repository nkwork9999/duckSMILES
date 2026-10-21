//! Bounded binary-treatment propensity, balance, matching and effect estimators.
//! Consistency/no interference, conditional exchangeability and positivity are
//! assumptions, not conclusions of this code. Each row has one factual outcome;
//! supplied outcome regressions are nuisance predictions, not observed counterfactuals.
//! IPW SE conditions on fixed nuisance scores. AIPW inference additionally needs
//! valid nuisance fitting/cross-fitting and rate/regularity conditions.

pub fn run(op: &str, d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let w = d.first().map_or(0, Vec::len);
    if w == 0
        || w > 66
        || d.len() > 100000
        || d.len().saturating_mul(w) > 1_000_000
        || d.iter()
            .any(|r| r.len() != w || r.iter().any(|v| !v.is_finite()))
        || a.len() > 16384
        || a.iter().any(|v| !v.is_finite())
    {
        return Err("causal methods require bounded finite rectangular data".into());
    }
    let out = match op {
        "propensity_fit" => fit(d, a),
        "propensity_predict" => predict(d, a),
        "causal_balance" => balance(d, a),
        "propensity_match" => matching(d, a),
        "ipw_effect" => effect(d, a, false),
        "aipw_effect" => effect(d, a, true),
        "uplift_curve" => uplift(d, a),
        _ => Err("unknown causal operation".into()),
    }?;
    if out
        .split([':', ',', '[', ']', '}', ' '])
        .any(|s| matches!(s, "NaN" | "inf" | "-inf"))
    {
        Err("causal result exceeds finite numeric range".into())
    } else {
        Ok(out)
    }
}
fn nums(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn opts(v: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|v| v.map_or_else(|| "null".into(), |v| v.to_string()))
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn norm(v: &[f64]) -> f64 {
    v.iter().fold(0_f64, |v, x| v.hypot(*x))
}
fn int(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("integer in [{lo},{hi}] required"))
    } else {
        Ok(x as usize)
    }
}
fn tolerance(x: f64) -> Result<f64, String> {
    if (1e-12..=1e-2).contains(&x) {
        Ok(x)
    } else {
        Err("tolerance must be in [1e-12,1e-2]".into())
    }
}
fn mode(a: &[f64], maximum: usize) -> Result<usize, String> {
    int(*a.first().unwrap_or(&0.), 0, maximum)
}
fn label(x: f64) -> Result<bool, String> {
    match x {
        0. => Ok(false),
        1. => Ok(true),
        _ => Err("binary treatment must be 0/1".into()),
    }
}
fn propensity(p: f64, clip: f64) -> Result<(f64, bool), String> {
    if p <= 0. || p >= 1. {
        return Err(
            "propensity must be strictly in (0,1); clipping does not restore failed positivity"
                .into(),
        );
    }
    let q = p.clamp(clip, 1. - clip);
    Ok((q, q != p))
}
fn clip(a: &[f64]) -> Result<f64, String> {
    let x = *a.get(1).unwrap_or(&0.);
    if !(0. ..0.5).contains(&x) {
        Err("explicit propensity clipping must be in [0,.5)".into())
    } else {
        Ok(x)
    }
}
fn groups(d: &[Vec<f64>], tcol: usize) -> Result<Vec<bool>, String> {
    let t = d
        .iter()
        .map(|r| label(r[tcol]))
        .collect::<Result<Vec<_>, _>>()?;
    if !t.contains(&true) || !t.contains(&false) {
        Err("both treatment groups are required".into())
    } else {
        Ok(t)
    }
}
fn normalized(w: &[f64]) -> Result<Vec<f64>, String> {
    let scale = w.iter().copied().fold(0_f64, f64::max);
    if scale <= 0. || !scale.is_finite() || w.iter().any(|w| *w < 0. || !w.is_finite()) {
        return Err("invalid or overflowing causal weights".into());
    }
    let total = w.iter().map(|w| w / scale).sum::<f64>();
    Ok(w.iter().map(|w| (w / scale) / total).collect())
}
fn inverse(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let scales = (0..n).map(|i| a[i][i].sqrt()).collect::<Vec<_>>();
    if scales.iter().any(|x| *x <= 0. || !x.is_finite()) {
        return Err("singular propensity design information".into());
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v =
                a[i][j] / scales[i] / scales[j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if v <= 1e-12 || !v.is_finite() {
                    return Err("rank deficient propensity information".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut out = vec![vec![0.; n]; n];
    for c in 0..n {
        let mut b = vec![0.; n];
        for i in 0..n {
            b[i] = ((if i == c { 1. } else { 0. }) - (0..i).map(|j| l[i][j] * b[j]).sum::<f64>())
                / l[i][i];
        }
        let mut x = vec![0.; n];
        for i in (0..n).rev() {
            x[i] = (b[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            out[i][c] = x[i] / scales[i] / scales[c];
        }
    }
    Ok(out)
}
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let a = x.exp();
        a / (1. + a)
    }
}
fn softplus(x: f64) -> f64 {
    if x >= 0. {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}
struct Frame {
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
}
impl Frame {
    fn new(d: &[Vec<f64>], p: usize) -> Result<Self, String> {
        let origins = d[0][..p].to_vec();
        let mut offsets = vec![];
        let mut scales = vec![];
        for c in 0..p {
            let off = d
                .iter()
                .map(|r| (r[c] - origins[c]) / d.len() as f64)
                .sum::<f64>();
            let sd = d.iter().fold(0_f64, |v, r| {
                v.hypot(((r[c] - origins[c]) - off) / (d.len() as f64).sqrt())
            });
            if !off.is_finite() || !sd.is_finite() {
                return Err("propensity scaling exceeds range".into());
            }
            offsets.push(off);
            scales.push(if sd > 0. { sd } else { 1. });
        }
        Ok(Self {
            origins,
            offsets,
            scales,
        })
    }
    fn row(&self, r: &[f64]) -> Result<Vec<f64>, String> {
        let mut v = vec![1.];
        v.extend(
            (0..self.origins.len())
                .map(|i| ((r[i] - self.origins[i]) - self.offsets[i]) / self.scales[i]),
        );
        if v.iter().any(|x| !x.is_finite()) {
            Err("propensity query exceeds range".into())
        } else {
            Ok(v)
        }
    }
}
fn loss(x: &[Vec<f64>], y: &[bool], b: &[f64], lambda: f64) -> f64 {
    x.iter()
        .zip(y)
        .map(|(x, y)| softplus(if *y { -dot(x, b) } else { dot(x, b) }))
        .sum::<f64>()
        + lambda * dot(&b[1..], &b[1..]) / 2.
}
fn fit(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.is_empty() || a.len() > 4 || d[0].len() > 33 {
        return Err(
            "propensity_fit requires [n_train,lambda=0,maxiter=100,tol=1e-8], <=32features".into(),
        );
    }
    let n = int(a[0], 3, d.len().min(10000))?;
    let p = d[0].len() - 1;
    if n <= p + 1 {
        return Err("propensity fit needs n_train>parameter count".into());
    }
    let lambda = *a.get(1).unwrap_or(&0.);
    if !(0. ..=1e12).contains(&lambda) {
        return Err("nonnegative bounded propensity ridge required".into());
    }
    let maxiter = int(*a.get(2).unwrap_or(&100.), 1, 2000)?;
    let tol = tolerance(*a.get(3).unwrap_or(&1e-8))?;
    if n.saturating_mul((p + 1) * (p + 1)).saturating_mul(maxiter) > 200_000_000 {
        return Err("propensity work limit exceeded".into());
    }
    let y = groups(&d[..n], p)?;
    let frame = Frame::new(&d[..n], p)?;
    let x = d
        .iter()
        .map(|r| frame.row(r))
        .collect::<Result<Vec<_>, _>>()?;
    let fraction = y.iter().filter(|y| **y).count() as f64 / n as f64;
    let mut b = vec![0.; p + 1];
    b[0] = (fraction / (1. - fraction)).ln();
    let mut objective = loss(&x[..n], &y, &b, lambda);
    let mut iterations = 0;
    for it in 0..maxiter {
        iterations = it + 1;
        let prob = x[..n]
            .iter()
            .map(|x| sigmoid(dot(x, &b)))
            .collect::<Vec<_>>();
        let mut g = (0..=p)
            .map(|c| {
                (0..n)
                    .map(|i| x[i][c] * (prob[i] - if y[i] { 1. } else { 0. }))
                    .sum::<f64>()
                    + if c > 0 { lambda * b[c] } else { 0. }
            })
            .collect::<Vec<_>>();
        let gradient_max = g.iter().fold(0_f64, |v, g| v.max(g.abs()));
        if gradient_max <= tol * n as f64 {
            break;
        }
        let h = (0..=p)
            .map(|c| {
                (0..=p)
                    .map(|k| {
                        (0..n)
                            .map(|i| x[i][c] * x[i][k] * prob[i] * (1. - prob[i]))
                            .sum::<f64>()
                            + if c == k && c > 0 { lambda } else { 0. }
                    })
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        let inverse = inverse(&h)?;
        g = inverse.iter().map(|r| dot(r, &g)).collect();
        let mut step = 1.;
        let mut accepted = false;
        for _ in 0..50 {
            let next = b
                .iter()
                .zip(&g)
                .map(|(b, g)| b - step * g)
                .collect::<Vec<_>>();
            let value = loss(&x[..n], &y, &next, lambda);
            if value.is_finite() && value < objective {
                b = next;
                objective = value;
                accepted = true;
                break;
            }
            step /= 2.;
        }
        if !accepted {
            break;
        }
    }
    let prob = x.iter().map(|x| sigmoid(dot(x, &b))).collect::<Vec<_>>();
    let gradient_max = (0..=p)
        .map(|c| {
            (0..n)
                .map(|i| x[i][c] * (prob[i] - if y[i] { 1. } else { 0. }))
                .sum::<f64>()
                + if c > 0 { lambda * b[c] } else { 0. }
        })
        .fold(0_f64, |v, g| v.max(g.abs()));
    let converged = gradient_max <= tol * n as f64;

    if lambda == 0. {
        let margins = x[..n]
            .iter()
            .zip(&y)
            .map(|(x, y)| dot(x, &b) * if *y { 1. } else { -1. })
            .collect::<Vec<_>>();
        if margins.iter().all(|m| *m >= -1e-12) && margins.iter().any(|m| *m > 1e-6) {
            return Err("unpenalized propensity design is separated; no finite MLE".into());
        }
    }
    let coef = b[1..]
        .iter()
        .zip(&frame.scales)
        .map(|(b, s)| b / s)
        .collect::<Vec<_>>();
    let intercept = b[0]
        - coef
            .iter()
            .enumerate()
            .map(|(i, b)| b * (frame.origins[i] + frame.offsets[i]))
            .sum::<f64>();
    let mut saved = vec![p as f64, b[0]];
    saved.extend(&frame.origins);
    saved.extend(&frame.offsets);
    saved.extend(&frame.scales);
    saved.extend(&b[1..]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"prefix_logistic_propensity\",\"n_train\":{n},\"lambda\":{lambda},\"intercept\":{intercept},\"coefficients\":{},\"propensities\":{},\"objective\":{objective},\"gradient_max\":{gradient_max},\"converged\":{converged},\"iterations\":{iterations},\"predict_parameters\":{},\"scope\":\"nuisance_score_model_not_causal_identification\"}}",nums(&coef),nums(&prob),nums(&saved)))
}
fn predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() < 2 {
        return Err("propensity model missing".into());
    }
    let p = int(a[0], 0, 32)?;
    if a.len() != 2 + 4 * p
        || (d[0].len() != p && !(p == 0 && d[0].len() == 1))
        || a[2 + 2 * p..2 + 3 * p].iter().any(|s| *s <= 0.)
    {
        return Err("propensity saved model dimension/scale mismatch".into());
    }
    let frame = Frame {
        origins: a[2..2 + p].to_vec(),
        offsets: a[2 + p..2 + 2 * p].to_vec(),
        scales: a[2 + 2 * p..2 + 3 * p].to_vec(),
    };
    let prob = d
        .iter()
        .map(|r| {
            frame
                .row(r)
                .map(|x| sigmoid(a[1] + dot(&x[1..], &a[2 + 3 * p..])))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"saved_logistic_propensity\",\"propensities\":{}}}",
        nums(&prob)
    ))
}
fn weight_vectors(t: &[bool], e: &[f64], att: usize) -> Result<(Vec<f64>, Vec<f64>), String> {
    let a = t
        .iter()
        .zip(e)
        .map(|(t, e)| {
            if *t {
                if att == 1 {
                    1.
                } else {
                    1. / e
                }
            } else {
                0.
            }
        })
        .collect::<Vec<_>>();
    let b = t
        .iter()
        .zip(e)
        .map(|(t, e)| {
            if *t {
                0.
            } else if att == 1 {
                e / (1. - e)
            } else {
                1. / (1. - e)
            }
        })
        .collect::<Vec<_>>();
    Ok((normalized(&a)?, normalized(&b)?))
}
fn balance(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() > 2 || d[0].len() < 3 {
        return Err(
            "balance requires [features...,treatment,propensity] and [estimand0ATE/1ATT,clip=0]"
                .into(),
        );
    }
    let p = d[0].len() - 2;
    let t = groups(d, p)?;
    let nt = t.iter().filter(|t| **t).count();
    let nc = d.len() - nt;
    if nt < 2 || nc < 2 {
        return Err("balance SMD requires >=2 rows in each treatment group".into());
    }
    let att = mode(a, 1)?;
    let threshold = clip(a)?;
    let raw = d
        .iter()
        .map(|r| propensity(r[p + 1], threshold))
        .collect::<Result<Vec<_>, _>>()?;
    let e = raw.iter().map(|x| x.0).collect::<Vec<_>>();
    let (w1, w0) = weight_vectors(&t, &e, att)?;
    let mut unweighted = vec![];
    let mut weighted = vec![];
    let mut means1 = vec![];
    let mut means0 = vec![];
    let mut baseline = vec![];
    for c in 0..p {
        let origin = d[0][c];
        let values = d.iter().map(|r| r[c] - origin).collect::<Vec<_>>();
        let mt = values
            .iter()
            .zip(&t)
            .filter(|(_, t)| **t)
            .map(|(x, _)| x / nt as f64)
            .sum::<f64>();
        let mc = values
            .iter()
            .zip(&t)
            .filter(|(_, t)| !**t)
            .map(|(x, _)| x / nc as f64)
            .sum::<f64>();
        let sd1 = values
            .iter()
            .zip(&t)
            .filter(|(_, t)| **t)
            .fold(0_f64, |v, (x, _)| {
                v.hypot((x - mt) / ((nt - 1) as f64).sqrt())
            });
        let sd0 = values
            .iter()
            .zip(&t)
            .filter(|(_, t)| !**t)
            .fold(0_f64, |v, (x, _)| {
                v.hypot((x - mc) / ((nc - 1) as f64).sqrt())
            });
        let sd = sd1.hypot(sd0) / 2_f64.sqrt();
        let wt = dot(&values, &w1);
        let wc = dot(&values, &w0);
        let smd = |diff: f64| {
            if sd > 0. {
                let z = diff / sd;
                z.is_finite().then_some(z)
            } else if diff == 0. {
                Some(0.)
            } else {
                None
            }
        };
        unweighted.push(smd(mt - mc));
        weighted.push(smd(wt - wc));
        means1.push(origin + wt);
        means0.push(origin + wc);
        baseline.push(sd);
    }
    let range = |group: bool| {
        let v = d
            .iter()
            .zip(&t)
            .filter(|(_, t)| **t == group)
            .map(|(row, _)| row[p + 1])
            .collect::<Vec<_>>();
        (
            v.iter().copied().fold(f64::INFINITY, f64::min),
            v.iter().copied().fold(f64::NEG_INFINITY, f64::max),
        )
    };
    let r1 = range(true);
    let r0 = range(false);
    let lo = r1.0.max(r0.0);
    let hi = r1.1.min(r0.1);
    let support = if lo <= hi {
        nums(&[lo, hi])
    } else {
        "null".into()
    };
    let outside = d
        .iter()
        .filter(|row| row[p + 1] < lo || row[p + 1] > hi)
        .count();
    let ess = [1. / dot(&w1, &w1), 1. / dot(&w0, &w0)];
    Ok(format!("{{\"version\":1,\"algorithm\":\"weighted_causal_balance\",\"estimand\":\"{}\",\"unweighted_smd\":{},\"weighted_smd\":{},\"weighted_treated_means\":{},\"weighted_control_means\":{},\"baseline_pooled_sample_sd\":{},\"effective_sample_sizes\":{},\"common_support\":{},\"outside_common_support\":{},\"n_clipped\":{},\"denominator\":\"fixed_unweighted_pooled_sample_SD\",\"assumptions_verified\":false}}",if att==0{"ATE"}else{"ATT"},opts(&unweighted),opts(&weighted),nums(&means1),nums(&means0),nums(&baseline),nums(&ess),support,outside,raw.iter().filter(|r|r.1).count()))
}
fn inference(estimate: f64, values: &[f64]) -> (f64, Vec<f64>, Vec<f64>) {
    let center = values.iter().map(|v| v / values.len() as f64).sum::<f64>();
    let centered = values.iter().map(|v| v - center).collect::<Vec<_>>();
    let se = norm(&centered) / ((values.len() * (values.len() - 1)) as f64).sqrt();
    (
        se,
        vec![
            estimate - 1.959963984540054 * se,
            estimate + 1.959963984540054 * se,
        ],
        centered,
    )
}
fn effect(d: &[Vec<f64>], a: &[f64], augmented: bool) -> Result<String, String> {
    if a.len() > 2 || d[0].len() != if augmented { 5 } else { 3 } || d.len() < 2 {
        return Err(
            "IPW rows[y,treatment,propensity]; AIPW adds[m0,m1]; args[estimand0ATE/1ATT,clip=0]"
                .into(),
        );
    }
    let att = mode(a, 1)?;
    let threshold = clip(a)?;
    let t = groups(d, 1)?;
    let pe = d
        .iter()
        .map(|r| propensity(r[2], threshold))
        .collect::<Result<Vec<_>, _>>()?;
    let e = pe.iter().map(|r| r.0).collect::<Vec<_>>();
    let n = d.len();
    let nt = t.iter().filter(|t| **t).count();
    let prevalence = nt as f64 / n as f64;
    let (estimate, influence, mu1, mu0) = if augmented {
        let score = d
            .iter()
            .zip(&t)
            .zip(&e)
            .map(|((r, t), e)| {
                if att == 0 {
                    r[4] - r[3]
                        + if *t {
                            (r[0] - r[4]) / e
                        } else {
                            -(r[0] - r[3]) / (1. - e)
                        }
                } else if *t {
                    r[0] - r[3]
                } else {
                    -e / (1. - e) * (r[0] - r[3])
                }
            })
            .collect::<Vec<_>>();
        let estimate = score
            .iter()
            .map(|v| v / if att == 0 { n as f64 } else { nt as f64 })
            .sum::<f64>();
        let influence = score
            .iter()
            .zip(&t)
            .map(|(g, t)| {
                if att == 0 {
                    g - estimate
                } else {
                    (g - if *t { estimate } else { 0. }) / prevalence
                }
            })
            .collect::<Vec<_>>();
        (estimate, influence, None, None)
    } else {
        let origin = d[0][0];
        let y = d.iter().map(|r| r[0] - origin).collect::<Vec<_>>();
        let (w1, w0) = weight_vectors(&t, &e, att)?;
        let m1 = dot(&w1, &y);
        let m0 = dot(&w0, &y);
        let influence = (0..n)
            .map(|i| n as f64 * (w1[i] * (y[i] - m1) - w0[i] * (y[i] - m0)))
            .collect();
        (m1 - m0, influence, Some(origin + m1), Some(origin + m0))
    };
    if !estimate.is_finite() || influence.iter().any(|v| !v.is_finite()) {
        return Err("causal score exceeds numerical range".into());
    }
    let (se, ci, influence) = inference(estimate, &influence);
    let means = opts(&[mu0, mu1]);
    Ok(format!("{{\"version\":1,\"algorithm\":\"{}\",\"estimand\":\"{}\",\"estimate\":{},\"standard_error\":{},\"confidence_interval_95\":{},\"influence_values\":{},\"mean_outcomes_control_treated\":{},\"n_clipped\":{},\"nuisance_inference\":\"{}\",\"assumptions\":\"consistency_no_interference_conditional_exchangeability_positivity_independent_rows\",\"assumptions_verified\":false}}",if augmented{"canonical_AIPW_score"}else{"Hajek_IPW_ratios"},if att==0{"ATE"}else{"ATT"},estimate,se,nums(&ci),nums(&influence),means,pe.iter().filter(|r|r.1).count(),if augmented{"external_or_crossfit_nuisance_under_orthogonality_rate_conditions_no_automatic_provenance_check"}else{"fixed_nuisance_delta_method_does_not_include_estimated_propensity_uncertainty"}))
}
fn matching(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if d[0].len() != 3 || a.len() > 3 || d.len() > 2000 {
        return Err("matching rows[y,treatment,propensity], <=2000rows; args[caliper=1,replacement=1,metric0probability/1logit=0]".into());
    }
    let caliper = *a.first().unwrap_or(&1.);
    if caliper < 0. {
        return Err("nonnegative matching caliper required".into());
    }
    let replacement = int(*a.get(1).unwrap_or(&1.), 0, 1)? == 1;
    let logit = int(*a.get(2).unwrap_or(&0.), 0, 1)? == 1;
    let t = groups(d, 1)?;
    let e = d
        .iter()
        .map(|r| propensity(r[2], 0.).map(|x| x.0))
        .collect::<Result<Vec<_>, _>>()?;
    let score = e
        .iter()
        .map(|p| if logit { p.ln() - (-p).ln_1p() } else { *p })
        .collect::<Vec<_>>();
    let mut used = vec![false; d.len()];
    let mut treated = vec![];
    let mut controls = vec![];
    let mut unmatched = vec![];
    let mut distance = vec![];
    let mut effects = vec![];
    for i in 0..d.len() {
        if !t[i] {
            continue;
        }
        let best = (0..d.len())
            .filter(|j| !t[*j] && (replacement || !used[*j]))
            .map(|j| (j, (score[i] - score[j]).abs()))
            .filter(|(_, v)| *v <= caliper)
            .min_by(|a, b| a.1.total_cmp(&b.1).then(a.0.cmp(&b.0)));
        if let Some((j, dist)) = best {
            treated.push(i as f64);
            controls.push(j as f64);
            distance.push(dist);
            effects.push(d[i][0] - d[j][0]);
            used[j] = true;
        } else {
            unmatched.push(i as f64);
        }
    }
    if effects.is_empty() {
        return Err("no propensity matches inside explicit caliper".into());
    }
    let estimate = effects
        .iter()
        .map(|e| e / effects.len() as f64)
        .sum::<f64>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"greedy_nearest_propensity_match\",\"estimand\":\"ATT_on_matched_treated_subset\",\"estimate\":{},\"treated_indices\":{},\"control_indices\":{},\"unmatched_treated_indices\":{},\"distances\":{},\"pair_differences\":{},\"caliper\":{},\"replacement\":{},\"metric\":\"{}\",\"order\":\"treated_input_order_exact_distance_ties_lowest_control_index\",\"standard_error\":null,\"inference\":\"not_supplied_control_reuse_and_estimated_matching_require_valid_matching_variance\",\"assumptions_verified\":false}}",estimate,nums(&treated),nums(&controls),nums(&unmatched),nums(&distance),nums(&effects),caliper,replacement,if logit{"logit"}else{"probability"}))
}
fn uplift(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if d[0].len() != 4 || !a.is_empty() {
        return Err("uplift evaluation rows[factual_outcome,treatment,predicted_effect_score,assignment_propensity], args[]".into());
    }
    let t = groups(d, 1)?;
    let p = d
        .iter()
        .map(|r| propensity(r[3], 0.).map(|v| v.0))
        .collect::<Result<Vec<_>, _>>()?;
    let z = d
        .iter()
        .enumerate()
        .map(|(i, r)| {
            if t[i] {
                r[0] / p[i]
            } else {
                -r[0] / (1. - p[i])
            }
        })
        .collect::<Vec<_>>();
    if z.iter().any(|v| !v.is_finite()) {
        return Err("uplift assignment-weighted outcome exceeds range".into());
    }
    let mut order = (0..d.len()).collect::<Vec<_>>();
    order.sort_by(|i, j| d[*j][2].total_cmp(&d[*i][2]).then(i.cmp(j)));
    let mut fraction = vec![0.];
    let mut gain = vec![0.];
    let mut thresholds = vec![None];
    let mut end = 0;
    let mut total = 0.;
    while end < order.len() {
        let score = d[order[end]][2];
        let mut next = end;
        while next < order.len() && d[order[next]][2] == score {
            total += z[order[next]] / d.len() as f64;
            next += 1;
        }
        fraction.push(next as f64 / d.len() as f64);
        gain.push(total);
        thresholds.push(Some(score));
        end = next;
    }
    let qini = fraction
        .iter()
        .zip(&gain)
        .map(|(f, g)| g - f * total)
        .collect::<Vec<_>>();
    let area = |v: &[f64]| {
        fraction
            .windows(2)
            .zip(v.windows(2))
            .map(|(f, g)| (f[1] - f[0]) * (g[0] / 2. + g[1] / 2.))
            .sum::<f64>()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"assignment_weighted_uplift_evaluation\",\"selection_fraction\":{},\"thresholds\":{},\"gain\":{},\"qini\":{},\"auuc\":{},\"qini_area\":{},\"population_HT_effect\":{},\"gain_units\":\"cumulative_HT_effect_divided_by_total_sample_size\",\"qini_definition\":\"gain_minus_selection_fraction_times_final_gain\",\"assumptions\":\"known_random_assignment_or_valid_propensity_with_exchangeability_positivity_consistency\",\"scope\":\"supplied_score_evaluation_no_uplift_tree_or_individual_counterfactual_truth\",\"assumptions_verified\":false}}",nums(&fraction),opts(&thresholds),nums(&gain),nums(&qini),area(&gain),area(&qini),total))
}
