//! Unidimensional 1PL/2PL/3PL marginal item-response models with fixed N(0,1).
//! Finite quadrature, constrained local optimization and identification diagnostics
//! are explicit. Posterior EAP ability scores are not JMP's person MLE formula.
use std::collections::BTreeMap;
const MAX_WORK: usize = 500_000_000;
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || width > 128
        || data.len() > 10000
        || data.len() > 1_000_000 / width
        || args.len() > 16384
        || data.iter().flatten().chain(args).any(|x| !x.is_finite())
        || data.iter().any(|r| r.len() != width)
    {
        return Err("IRT needs finite nonempty rectangular bounded data/parameters".into());
    }
    match operation {
        "irt_fit" => fit(data, args),
        "irt_predict" => predict(data, args),
        "irt_score" => score(data, args),
        _ => Err("unknown IRT operation".into()),
    }
}
struct Model {
    kind: usize,
    items: usize,
    missing: usize,
    order: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    c: Vec<f64>,
}
impl Model {
    fn unpack(args: &[f64]) -> Result<Self, String> {
        if args.len() < 5 || args[0] != 1. {
            return Err("IRT needs unchanged version1 model_parameters".into());
        }
        let kind = integer(args[1], 1, 3, "IRT model")?;
        let items = integer(args[2], 1, 12, "item count")?;
        let missing = integer(args[3], 0, 2, "missing code")?;
        let order = integer(args[4], 11, 89, "saved quadrature order")?;
        if order % 2 != 1 || args.len() != 5 + 3 * items {
            return Err("IRT saved model length/order is invalid".into());
        }
        let mut a = Vec::new();
        let mut b = Vec::new();
        let mut c = Vec::new();
        for item in args[5..].chunks_exact(3) {
            if !(0.05..=5.).contains(&item[0])
                || item[1].abs() > 8.
                || !(0. ..=0.8).contains(&item[2])
                || (kind == 1 && item[0] != 1.)
                || (kind < 3 && item[2] != 0.)
            {
                return Err("IRT model needs discrimination[0.05,5],difficulty[-8,8],guessing[0,0.8],1PL a1,and1/2PL c0".into());
            }
            a.push(item[0]);
            b.push(item[1]);
            c.push(item[2]);
        }
        Ok(Self {
            kind,
            items,
            missing,
            order,
            a,
            b,
            c,
        })
    }
    fn parameters(&self) -> Vec<f64> {
        let mut v = vec![
            1.,
            self.kind as f64,
            self.items as f64,
            self.missing as f64,
            self.order as f64,
        ];
        for j in 0..self.items {
            v.extend([self.a[j], self.b[j], self.c[j]]);
        }
        v
    }
}
fn response(row: &[f64], missing: usize) -> Result<Vec<f64>, String> {
    row.iter()
        .map(|v| match *v {
            0. | 1. => Ok(*v),
            -1. if missing == 1 => Ok(0.),
            -1. if missing == 2 => Ok(-1.),
            _ => Err("IRT responses must be0/1; -1requires explicit missing policy".into()),
        })
        .collect()
}
fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let model = Model::unpack(args)?;
    if data[0].len() != 1 {
        return Err("irt_predict accepts scalar ability rows".into());
    }
    work(data.len(), model.items, 1)?;
    let mut probabilities = Vec::new();
    let mut information = Vec::new();
    let mut test_info = Vec::new();
    let mut expected = Vec::new();
    let mut variances = Vec::new();
    for row in data {
        let mut p = Vec::new();
        let mut info = Vec::new();
        for j in 0..model.items {
            let s = sigmoid(model.a[j] * (row[0] - model.b[j]));
            let probability = model.c[j] + (1. - model.c[j]) * s;
            let value = if model.c[j] == 0. {
                model.a[j].powi(2) * s * (1. - s)
            } else {
                model.a[j].powi(2) * (1. - model.c[j]) * s * s * (1. - s) / probability
            };
            p.push(probability);
            info.push(value);
        }
        test_info.push(info.iter().sum());
        expected.push(p.iter().sum());
        variances.push(p.iter().map(|v| v * (1. - v)).sum());
        probabilities.push(p);
        information.push(info);
    }
    Ok(format!("{{\"version\":1,\"probabilities\":{},\"item_information\":{},\"test_information\":{},\"expected_total_score\":{},\"conditional_score_variance\":{},\"parameterization\":\"c+(1-c)*sigmoid(a*(theta-b));D=1\"}}",matrix(&probabilities),matrix(&information),numbers(&test_info),numbers(&expected),numbers(&variances)))
}
struct Rule {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}
fn normal_rule(order: usize) -> Result<Rule, String> {
    let mut jacobi = vec![vec![0.; order]; order];
    for j in 1..order {
        jacobi[j - 1][j] = (j as f64).sqrt();
        jacobi[j][j - 1] = jacobi[j - 1][j];
    }
    let e = eigen(&jacobi)?;
    let mut nodes = e.values;
    nodes.reverse();
    // Derivative/orthonormal Hermite weights avoid using tiny eigenvector entries
    // whose relative accuracy is lost for extreme high-order tail nodes.
    let mut weights = Vec::new();
    for &z in &nodes {
        let mut previous = 1.;
        let mut current = z;
        for k in 2..order {
            let next =
                z / (k as f64).sqrt() * current - ((k - 1) as f64 / k as f64).sqrt() * previous;
            previous = current;
            current = next;
        }
        weights.push(finite(
            1. / (order as f64 * current * current),
            "normal quadrature weight",
        )?);
    }
    let total: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= total;
        if *w <= 0. {
            return Err("normal quadrature weight underflow".into());
        }
    }
    Ok(Rule { nodes, weights })
}
struct Moment {
    log_probability: f64,
    mean: f64,
    variance: f64,
    weights: Vec<f64>,
}
struct NodeItem {
    success: f64,
    failure: f64,
    eta: f64,
    deta_success: f64,
    deta_failure: f64,
    dc_success: f64,
    dc_failure: f64,
}
fn node_items(model: &Model, rule: &Rule) -> Vec<Vec<NodeItem>> {
    rule.nodes
        .iter()
        .map(|theta| {
            (0..model.items)
                .map(|j| {
                    let eta = model.a[j] * (theta - model.b[j]);
                    let s = sigmoid(eta);
                    let log_s = -softplus(-eta);
                    let log_failure = -softplus(eta);
                    let lc = (-model.c[j]).ln_1p();
                    let success = if model.c[j] == 0. {
                        log_s
                    } else {
                        logadd(model.c[j].ln(), lc + log_s)
                    };
                    NodeItem {
                        success,
                        failure: lc + log_failure,
                        eta,
                        deta_success: if model.c[j] == 0. {
                            1. - s
                        } else {
                            (lc + log_s + log_failure - success).exp()
                        },
                        deta_failure: -s,
                        dc_success: (log_failure - success).exp(),
                        dc_failure: -1. / (1. - model.c[j]),
                    }
                })
                .collect()
        })
        .collect()
}
fn moment(row: &[f64], nodes: &[Vec<NodeItem>], rule: &Rule) -> Result<Moment, String> {
    if row.iter().all(|x| *x < 0.) {
        return Ok(Moment {
            log_probability: 0.,
            mean: 0.,
            variance: 1.,
            weights: rule.weights.clone(),
        });
    }
    let logarithms = nodes
        .iter()
        .zip(&rule.weights)
        .map(|(items, w)| {
            w.ln()
                + row
                    .iter()
                    .zip(items)
                    .filter(|(x, _)| **x >= 0.)
                    .map(|(x, item)| if *x == 1. { item.success } else { item.failure })
                    .sum::<f64>()
        })
        .collect::<Vec<_>>();
    let maximum = logarithms.iter().copied().fold(-f64::INFINITY, f64::max);
    let mut weights = logarithms
        .iter()
        .map(|x| (x - maximum).exp())
        .collect::<Vec<_>>();
    let sum: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= sum;
    }
    let mean = finite(dot(&rule.nodes, &weights), "posterior ability mean")?;
    let variance = finite(
        rule.nodes
            .iter()
            .zip(&weights)
            .map(|(x, w)| w * (x - mean).powi(2))
            .sum(),
        "posterior ability variance",
    )?;
    Ok(Moment {
        log_probability: finite(maximum + sum.ln(), "marginal pattern log probability")?,
        mean,
        variance,
        weights,
    })
}
fn score(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let model = Model::unpack(args)?;
    if data[0].len() != model.items {
        return Err("irt_score needs one response per saved item".into());
    }
    work(data.len(), model.items, 2 * (model.order + 8))?;
    let low_rule = normal_rule(model.order)?;
    let rule = normal_rule(model.order + 8)?;
    let low_nodes = node_items(&model, &low_rule);
    let nodes = node_items(&model, &rule);
    let mut means = Vec::new();
    let mut variances = Vec::new();
    let mut logs = Vec::new();
    let mut counts = Vec::new();
    let mut incorrect = Vec::new();
    let mut correct = Vec::new();
    let mut missing = Vec::new();
    let mut differences = Vec::new();
    let mut converged = true;
    for raw in data {
        let row = response(raw, model.missing)?;
        let low = moment(&row, &low_nodes, &low_rule)?;
        let high = moment(&row, &nodes, &rule)?;
        let delta = (low.log_probability - high.log_probability)
            .abs()
            .max((low.mean - high.mean).abs())
            .max((low.variance - high.variance).abs());
        converged &= delta <= 1e-7 * (1. + high.log_probability.abs());
        differences.push(delta);
        let observed = row.iter().filter(|x| **x >= 0.).count();
        counts.push(observed as f64);
        missing.push(observed == 0);
        incorrect.push(observed > 0 && row.iter().filter(|x| **x >= 0.).all(|x| *x == 0.));
        correct.push(observed > 0 && row.iter().filter(|x| **x >= 0.).all(|x| *x == 1.));
        means.push(high.mean);
        variances.push(high.variance);
        logs.push(high.log_probability);
    }
    let sd = variances.iter().map(|v| v.sqrt()).collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"method\":\"posterior_EAP_fixed_standard_normal\",\"posterior_means\":{},\"posterior_variances\":{},\"posterior_standard_deviations\":{},\"marginal_log_probabilities\":{},\"observed_items\":{},\"all_incorrect\":{},\"all_correct\":{},\"all_missing\":{},\"quadrature_orders\":[{},{}],\"quadrature_order_deltas\":{},\"quadrature_converged\":{converged},\"uncertainty\":\"conditional_on_supplied_item_parameters\",\"JMP_difference\":\"EAP_returns_finite_all0_all1_scores_instead_of_person_MLE\"}}",numbers(&means),numbers(&variances),numbers(&sd),numbers(&logs),numbers(&counts),booleans(&incorrect),booleans(&correct),booleans(&missing),model.order,model.order+8,numbers(&differences)))
}
struct Dataset {
    patterns: Vec<Vec<f64>>,
    frequencies: Vec<f64>,
    total: f64,
    items: usize,
    kind: usize,
    missing: usize,
    penalty: f64,
    guess_upper: f64,
}
impl Dataset {
    fn new(
        data: &[Vec<f64>],
        kind: usize,
        frequency: bool,
        missing: usize,
        penalty: f64,
        guess_upper: f64,
    ) -> Result<Self, String> {
        let items = data[0].len() - usize::from(frequency);
        let minimum = if kind == 1 {
            2
        } else if kind == 2 {
            3
        } else {
            5
        };
        if data.len() > 2000 || !(minimum..=12).contains(&items) {
            return Err(
                "IRT fit requires<=2000rows and2..12items for1PL,3..12for2PL,5..12for3PL".into(),
            );
        }
        let mut collapsed = BTreeMap::<u32, (Vec<f64>, f64)>::new();
        for raw in data {
            let row = response(&raw[..items], missing)?;
            let freq = if frequency {
                let value = raw[items];
                integer(value, 0, 1_000_000, "integer frequency")?;
                value
            } else {
                1.
            };
            if freq == 0. || row.iter().all(|x| *x < 0.) {
                continue;
            }
            let key = row.iter().enumerate().fold(0_u32, |key, (j, x)| {
                key | ((if *x < 0. { 2 } else { *x as u32 }) << (2 * j))
            });
            let entry = collapsed.entry(key).or_insert((row, 0.));
            entry.1 += freq;
        }
        let (patterns, frequencies): (Vec<_>, Vec<_>) = collapsed.into_values().unzip();
        let total: f64 = frequencies.iter().sum();
        if total > 10_000_000. || total <= ((kind * items) + 1) as f64 {
            return Err(
                "IRT effective frequency must exceed parameter_count+1 and be<=10000000".into(),
            );
        }
        for j in 0..items {
            let positives: f64 = patterns
                .iter()
                .zip(&frequencies)
                .filter(|(r, _)| r[j] == 1.)
                .map(|(_, w)| w)
                .sum();
            let negatives: f64 = patterns
                .iter()
                .zip(&frequencies)
                .filter(|(r, _)| r[j] == 0.)
                .map(|(_, w)| w)
                .sum();
            if positives == 0. || negatives == 0. {
                return Err("each fitted item needs both correct and incorrect observed responses; infinite/unidentified item estimates rejected".into());
            }
        }
        Ok(Self {
            patterns,
            frequencies,
            total,
            items,
            kind,
            missing,
            penalty,
            guess_upper,
        })
    }
    fn model(&self, theta: &[f64], order: usize) -> Model {
        let k = self.items;
        Model {
            kind: self.kind,
            items: k,
            missing: self.missing,
            order,
            b: theta[..k].to_vec(),
            a: if self.kind >= 2 {
                theta[k..2 * k].iter().map(|x| x.exp()).collect()
            } else {
                vec![1.; k]
            },
            c: if self.kind == 3 {
                theta[2 * k..].to_vec()
            } else {
                vec![0.; k]
            },
        }
    }
    fn bounds(&self) -> (Vec<f64>, Vec<f64>) {
        let mut lower = vec![-8.; self.items];
        let mut upper = vec![8.; self.items];
        if self.kind >= 2 {
            lower.extend(vec![0.05_f64.ln(); self.items]);
            upper.extend(vec![5_f64.ln(); self.items]);
        }
        if self.kind == 3 {
            lower.extend(vec![0.; self.items]);
            upper.extend(vec![self.guess_upper; self.items]);
        }
        (lower, upper)
    }
}
struct Evaluation {
    objective: f64,
    gradient: Vec<f64>,
    log_likelihood: f64,
    moments: Vec<Moment>,
}
fn evaluate(dataset: &Dataset, theta: &[f64], rule: &Rule) -> Result<Evaluation, String> {
    let model = dataset.model(theta, rule.nodes.len());
    let nodes = node_items(&model, rule);
    let mut gradient = vec![0.; theta.len()];
    let mut log_likelihood = 0.;
    let mut moments = Vec::new();
    let k = dataset.items;
    for (row, &frequency) in dataset.patterns.iter().zip(&dataset.frequencies) {
        let posterior = moment(row, &nodes, rule)?;
        log_likelihood += frequency * posterior.log_probability;
        for (items, &weight) in nodes.iter().zip(&posterior.weights) {
            if weight == 0. {
                continue;
            }
            let factor = -frequency / dataset.total * weight;
            for j in 0..k {
                if row[j] < 0. {
                    continue;
                }
                let deta = if row[j] == 1. {
                    items[j].deta_success
                } else {
                    items[j].deta_failure
                };
                gradient[j] += factor * (-model.a[j] * deta);
                if dataset.kind >= 2 {
                    gradient[k + j] += factor * items[j].eta * deta;
                }
                if dataset.kind == 3 {
                    gradient[2 * k + j] += factor
                        * if row[j] == 1. {
                            items[j].dc_success
                        } else {
                            items[j].dc_failure
                        };
                }
            }
        }
        moments.push(posterior);
    }
    let penalty = dataset.penalty * model.c.iter().map(|v| v * v).sum::<f64>() / 2.;
    if dataset.kind == 3 {
        for j in 0..k {
            gradient[2 * k + j] += dataset.penalty / dataset.total * model.c[j];
        }
    }
    if gradient.iter().any(|v| !v.is_finite()) {
        return Err("IRT integrated score exceeds finite range".into());
    }
    Ok(Evaluation {
        objective: finite(
            (-log_likelihood + penalty) / dataset.total,
            "normalized penalized marginal likelihood",
        )?,
        gradient,
        log_likelihood: finite(log_likelihood, "marginal log likelihood")?,
        moments,
    })
}
struct Optimum {
    theta: Vec<f64>,
    evaluation: Evaluation,
    iterations: usize,
    converged: bool,
    gradient_max: f64,
    line_failed: bool,
}
fn active(x: f64, g: f64, lo: f64, hi: f64) -> bool {
    let slack = 1e-12 * (hi - lo).min(1.);
    (x <= lo + slack && g > 0.) || (x >= hi - slack && g < 0.)
}
fn optimize(
    dataset: &Dataset,
    mut theta: Vec<f64>,
    rule: &Rule,
    maxiter: usize,
    tol: f64,
) -> Result<Optimum, String> {
    let d = theta.len();
    let (lower, upper) = dataset.bounds();
    let mut current = evaluate(dataset, &theta, rule)?;
    let mut inverse = identity(d);
    let mut iterations = 0;
    let mut line_failed = false;
    for iteration in 0..maxiter {
        let active_set = (0..d)
            .map(|j| active(theta[j], current.gradient[j], lower[j], upper[j]))
            .collect::<Vec<_>>();
        let projected = (0..d)
            .map(|j| {
                if active_set[j] {
                    0.
                } else {
                    current.gradient[j]
                }
            })
            .collect::<Vec<_>>();
        if projected.iter().map(|v| v.abs()).fold(0_f64, f64::max) <= tol {
            break;
        }
        let mut direction = inverse
            .iter()
            .map(|r| -dot(r, &projected))
            .collect::<Vec<_>>();
        for j in 0..d {
            if active_set[j] {
                direction[j] = 0.;
            }
        }
        if dot(&direction, &projected) >= 0. {
            direction = projected.iter().map(|v| -v).collect();
            inverse = identity(d);
        }
        let mut alpha = 1.;
        let mut accepted = None;
        for _ in 0..50 {
            let candidate = (0..d)
                .map(|j| (theta[j] + alpha * direction[j]).clamp(lower[j], upper[j]))
                .collect::<Vec<_>>();
            let step = (0..d).map(|j| candidate[j] - theta[j]).collect::<Vec<_>>();
            let slope = dot(&current.gradient, &step);
            if slope < 0. {
                if let Ok(next) = evaluate(dataset, &candidate, rule) {
                    if next.objective
                        <= current.objective
                            + 1e-4 * slope
                            + 16. * f64::EPSILON * (1. + current.objective.abs())
                    {
                        accepted = Some((candidate, step, next));
                        break;
                    }
                }
            }
            alpha *= 0.5;
        }
        let Some((candidate, step, next)) = accepted else {
            line_failed = true;
            break;
        };
        let changed = (0..d)
            .any(|j| active_set[j] != active(candidate[j], next.gradient[j], lower[j], upper[j]));
        let change = (0..d)
            .map(|j| {
                if active_set[j] {
                    0.
                } else {
                    next.gradient[j] - current.gradient[j]
                }
            })
            .collect::<Vec<_>>();
        let curvature = dot(&step, &change);
        if !changed
            && curvature > 0.
            && curvature > 1e-12 * norm(step.iter().copied()) * norm(change.iter().copied())
        {
            let hy = inverse.iter().map(|r| dot(r, &change)).collect::<Vec<_>>();
            let coefficient = (1. + dot(&change, &hy) / curvature) / curvature;
            for i in 0..d {
                for j in 0..d {
                    inverse[i][j] = finite(
                        inverse[i][j] + coefficient * step[i] * step[j]
                            - (step[i] * hy[j] + hy[i] * step[j]) / curvature,
                        "IRT BFGS inverse Hessian",
                    )?;
                }
            }
        } else {
            inverse = identity(d);
        }
        theta = candidate;
        current = next;
        iterations = iteration + 1;
    }
    let gradient_max = (0..d)
        .map(|j| {
            if active(theta[j], current.gradient[j], lower[j], upper[j]) {
                0.
            } else {
                current.gradient[j].abs()
            }
        })
        .fold(0_f64, f64::max);
    Ok(Optimum {
        theta,
        evaluation: current,
        iterations,
        converged: gradient_max <= tol,
        gradient_max,
        line_failed,
    })
}
struct Information {
    covariance: Vec<Vec<Option<f64>>>,
    standard_errors: Vec<Option<f64>>,
    free: Vec<usize>,
    minimum_eigenvalue: f64,
}
#[allow(clippy::needless_range_loop)]
fn information(dataset: &Dataset, theta: &[f64], rule: &Rule) -> Result<Information, String> {
    let d = theta.len();
    let (lower, upper) = dataset.bounds();
    let base = evaluate(dataset, theta, rule)?;
    let free = (0..d)
        .filter(|&j| {
            !(dataset.kind == 3
                && j >= 2 * dataset.items
                && active(theta[j], base.gradient[j], lower[j], upper[j]))
        })
        .collect::<Vec<_>>();
    let mut hessian = vec![vec![0.; free.len()]; free.len()];
    for (column, &j) in free.iter().enumerate() {
        // A small positive guessing parameter is free, not an active zero.
        // Keep its central-score stencil inside the actual feasible interval.
        let margin = (theta[j] - lower[j]).min(upper[j] - theta[j]);
        let h = (1e-4 * (1. + theta[j].abs())).min(0.25 * margin);
        if h <= 1e-10 * (1. + theta[j].abs()) {
            return Err("information free-parameter stencil is numerically unresolved".into());
        }
        let mut plus = theta.to_vec();
        let mut minus = theta.to_vec();
        plus[j] += h;
        minus[j] -= h;
        let gp = evaluate(dataset, &plus, rule)?.gradient;
        let gm = evaluate(dataset, &minus, rule)?.gradient;
        for (row, &i) in free.iter().enumerate() {
            hessian[row][column] = finite(
                (gp[i] - gm[i]) / (2. * h) * dataset.total,
                "IRT observed information",
            )?;
        }
    }
    for i in 0..free.len() {
        for j in 0..i {
            let value = (hessian[i][j] + hessian[j][i]) / 2.;
            hessian[i][j] = value;
            hessian[j][i] = value;
        }
    }
    let e = eigen(&hessian)?;
    let minimum_eigenvalue = *e.values.last().unwrap();
    if minimum_eigenvalue <= e.values[0].abs() * 1e-10 {
        return Err("IRT free observed information is not positive identified".into());
    }
    let model = dataset.model(theta, rule.nodes.len());
    let mut covariance = vec![vec![None; d]; d];
    let mut standard_errors = vec![None; d];
    for (i, &index) in free.iter().enumerate() {
        let left = if dataset.kind >= 2 && (dataset.items..2 * dataset.items).contains(&index) {
            model.a[index - dataset.items]
        } else {
            1.
        };
        for (j, &other) in free.iter().enumerate() {
            let right = if dataset.kind >= 2 && (dataset.items..2 * dataset.items).contains(&other)
            {
                model.a[other - dataset.items]
            } else {
                1.
            };
            covariance[index][other] = Some(finite(
                (0..free.len())
                    .map(|k| e.vectors[i][k] * e.vectors[j][k] / e.values[k])
                    .sum::<f64>()
                    * left
                    * right,
                "IRT parameter covariance",
            )?);
        }
        standard_errors[index] = Some(covariance[index][index].unwrap().sqrt());
    }
    Ok(Information {
        covariance,
        standard_errors,
        free,
        minimum_eigenvalue,
    })
}
fn fit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 8 {
        return Err("irt_fit args [model1/2/3=2,frequency0/1=0,missing0reject/1wrong/2omit=0,maxiter=300,tol=1e-6,oddq=41,guess_upper=.5,guess_penalty=0]".into());
    }
    let kind = integer(args.first().copied().unwrap_or(2.), 1, 3, "IRT model")?;
    let frequency = integer(args.get(1).copied().unwrap_or(0.), 0, 1, "frequency column")? == 1;
    let missing = integer(args.get(2).copied().unwrap_or(0.), 0, 2, "missing policy")?;
    let maxiter = integer(
        args.get(3).copied().unwrap_or(300.),
        1,
        2000,
        "iteration limit",
    )?;
    let tol = args.get(4).copied().unwrap_or(1e-6);
    let order = integer(
        args.get(5).copied().unwrap_or(41.),
        11,
        81,
        "quadrature order",
    )?;
    let guess_upper = args.get(6).copied().unwrap_or(0.5);
    let penalty = args.get(7).copied().unwrap_or(0.);
    if !(1e-10..=1e-3).contains(&tol)
        || order % 2 != 1
        || !(0.01..=0.8).contains(&guess_upper)
        || !(0. ..=1e6).contains(&penalty)
        || (kind < 3 && penalty != 0.)
    {
        return Err(
            "IRT needs tol[1e-10,1e-3],oddq11..81,guess_upper[.01,.8],penalty[0,1e6]onlyfor3PL"
                .into(),
        );
    }
    let dataset = Dataset::new(data, kind, frequency, missing, penalty, guess_upper)?;
    let k = dataset.items;
    let dimensions = kind * k;
    let starts = if kind == 1 { 1 } else { 3 };
    work(dataset.patterns.len(), k * order, starts * maxiter)?;
    let rule = normal_rule(order)?;
    let high_rule = normal_rule(order + 8)?;
    let proportions = (0..k)
        .map(|j| {
            let observed: f64 = dataset
                .patterns
                .iter()
                .zip(&dataset.frequencies)
                .filter(|(r, _)| r[j] >= 0.)
                .map(|(_, w)| w)
                .sum();
            dataset
                .patterns
                .iter()
                .zip(&dataset.frequencies)
                .filter(|(r, _)| r[j] == 1.)
                .map(|(_, w)| w)
                .sum::<f64>()
                / observed
        })
        .collect::<Vec<_>>();
    let mut best: Option<Optimum> = None;
    let mut converged_starts = 0;
    for index in 0..starts {
        let mut theta = proportions
            .iter()
            .map(|p| ((1. - p) / p).ln().clamp(-7., 7.))
            .collect::<Vec<_>>();
        if kind >= 2 {
            theta.extend(vec![[0., -0.3, 0.3][index]; k]);
        }
        if kind == 3 {
            theta.extend(vec![[0.1_f64, 0., 0.2][index].min(guess_upper); k]);
        }
        if let Ok(candidate) = optimize(&dataset, theta, &rule, maxiter, tol) {
            converged_starts += usize::from(candidate.converged);
            if best.as_ref().is_none_or(|b| {
                candidate.evaluation.objective < b.evaluation.objective - 1e-12
                    || ((candidate.evaluation.objective - b.evaluation.objective).abs() <= 1e-12
                        && candidate.converged
                        && !b.converged)
            }) {
                best = Some(candidate);
            }
        }
    }
    let best = best.ok_or("all IRT marginal likelihood starts failed")?;
    let high = evaluate(&dataset, &best.theta, &high_rule)?;
    let model = dataset.model(&best.theta, order + 8);
    let ll_delta = (high.log_likelihood - best.evaluation.log_likelihood).abs();
    let pattern_deltas = high
        .moments
        .iter()
        .zip(&best.evaluation.moments)
        .map(|(a, b)| (a.log_probability - b.log_probability).abs())
        .collect::<Vec<_>>();
    let quadrature_converged = ll_delta <= tol.max(1e-9) * (1. + high.log_likelihood.abs())
        && high
            .moments
            .iter()
            .zip(&pattern_deltas)
            .all(|(p, d)| *d <= tol.max(1e-9) * (1. + p.log_probability.abs()));
    let (lower, upper) = dataset.bounds();
    let numerical_guard_bound = (0..if kind == 1 { k } else { 2 * k })
        .any(|j| best.theta[j] <= lower[j] + 1e-7 || best.theta[j] >= upper[j] - 1e-7);
    let guess_zero = (0..k)
        .map(|j| kind == 3 && model.c[j] <= 1e-10)
        .collect::<Vec<_>>();
    let guess_top = (0..k)
        .map(|j| kind == 3 && model.c[j] >= guess_upper - 1e-10)
        .collect::<Vec<_>>();
    let info = information(&dataset, &best.theta, &high_rule).ok();
    let information_identified = info.is_some();
    let converged =
        best.converged && quadrature_converged && !numerical_guard_bound && information_identified;
    let log_likelihood = high.log_likelihood;
    let penalized = log_likelihood - penalty * model.c.iter().map(|v| v * v).sum::<f64>() / 2.;
    let aic = if converged && penalty == 0. {
        Some(-2. * log_likelihood + 2. * dimensions as f64)
    } else {
        None
    };
    let bic = if converged && penalty == 0. {
        Some(-2. * log_likelihood + dimensions as f64 * dataset.total.ln())
    } else {
        None
    };
    let reported_info = info.as_ref().filter(|_| converged);
    let inference_valid = converged;
    let uncertainty_kind = if penalty > 0. {
        "penalized_objective_curvature_not_MLE_sampling_covariance"
    } else {
        "observed_ML_curvature_conditional_active_guessing_face"
    };
    let reported_order = order + 8;
    let covariance = reported_info.map_or("null".into(), |x| {
        format!(
            "[{}]",
            x.covariance
                .iter()
                .map(|r| optional_numbers(r))
                .collect::<Vec<_>>()
                .join(",")
        )
    });
    let se = reported_info.map_or("null".into(), |x| optional_numbers(&x.standard_errors));
    let free = info
        .as_ref()
        .map_or(vec![], |x| x.free.iter().map(|j| *j as f64).collect());
    let means = high.moments.iter().map(|m| m.mean).collect::<Vec<_>>();
    let variances = high.moments.iter().map(|m| m.variance).collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"fixed_standard_normal_logistic_IRT_marginal_likelihood\",\"model\":{kind},\"items\":{k},\"effective_frequency\":{},\"response_patterns\":{},\"pattern_frequencies\":{},\"difficulties\":{},\"discriminations\":{},\"guessing\":{},\"log_likelihood\":{log_likelihood},\"guessing_penalty\":{penalty},\"penalized_log_likelihood\":{penalized},\"penalty_definition\":\"lambda/2*sum(c_j^2)_on_total_loglikelihood\",\"aic\":{},\"bic\":{},\"parameter_count\":{dimensions},\"optimizer_converged\":{},\"converged\":{converged},\"iterations\":{},\"projected_gradient_max\":{},\"line_search_failed\":{},\"converged_starts\":{converged_starts},\"numerical_guard_bound\":{numerical_guard_bound},\"guessing_zero_bound\":{},\"guessing_upper_bound\":{},\"quadrature_orders\":[{order},{}],\"log_likelihood_order_delta\":{ll_delta},\"pattern_log_probability_order_deltas\":{},\"quadrature_converged\":{quadrature_converged},\"information_identified\":{information_identified},\"inference_valid\":{inference_valid},\"uncertainty_kind\":\"{uncertainty_kind}\",\"optimizer_gradient_quadrature_order\":{order},\"reported_likelihood_quadrature_order\":{reported_order},\"information_quadrature_order\":{reported_order},\"minimum_information_eigenvalue\":{},\"parameter_covariance\":{covariance},\"parameter_standard_errors\":{se},\"parameter_order\":\"difficulties_then_discriminations_if2PL3PL_then_guessing_if3PL\",\"free_parameter_indices\":{},\"uncertainty\":\"inverse_observed_penalized_information_conditional_on_active_guessing_face\",\"pattern_posterior_means\":{},\"pattern_posterior_variances\":{},\"model_parameters\":{},\"quadrature_error_convention\":\"q_versus_q_plus_8_consistency_not_rigorous_bound\"}}",dataset.total,matrix(&dataset.patterns),numbers(&dataset.frequencies),numbers(&model.b),numbers(&model.a),numbers(&model.c),optional(aic),optional(bic),best.converged,best.iterations,best.gradient_max,best.line_failed,booleans(&guess_zero),booleans(&guess_top),order+8,numbers(&pattern_deltas),optional(info.as_ref().map(|x|x.minimum_eigenvalue)),numbers(&free),numbers(&means),numbers(&variances),numbers(&model.parameters())))
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn integer(x: f64, lo: usize, hi: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("{name} must be integer[{lo},{hi}]"))
    } else {
        Ok(x as usize)
    }
}
fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.saturating_mul(b).saturating_mul(c) > MAX_WORK {
        Err("IRT nominal computation exceeds work limit".into())
    } else {
        Ok(())
    }
}
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
fn softplus(x: f64) -> f64 {
    if x > 0. {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
}
fn logadd(x: f64, y: f64) -> f64 {
    let m = x.max(y);
    m + ((x - m).exp() + (y - m).exp()).ln()
}
fn optional(x: Option<f64>) -> String {
    x.map_or("null".into(), |v| v.to_string())
}
fn optional_numbers(x: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|v| optional(*v)).collect::<Vec<_>>().join(",")
    )
}
fn booleans(x: &[bool]) -> String {
    format!(
        "[{}]",
        x.iter().map(bool::to_string).collect::<Vec<_>>().join(",")
    )
}
fn numbers(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}
fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn norm(x: impl Iterator<Item = f64>) -> f64 {
    x.fold(0_f64, f64::hypot)
}
fn identity(p: usize) -> Vec<Vec<f64>> {
    (0..p)
        .map(|i| (0..p).map(|j| f64::from(i == j)).collect())
        .collect()
}

struct Eigen {
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}
#[allow(clippy::needless_range_loop)]
fn eigen(matrix: &[Vec<f64>]) -> Result<Eigen, String> {
    let p = matrix.len();
    work(p, p, 64 * p)?;
    let scale = matrix
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(0_f64, f64::max);
    if scale == 0. {
        return Ok(Eigen {
            values: vec![0.; p],
            vectors: identity(p),
        });
    }
    let mut a = matrix
        .iter()
        .map(|r| r.iter().map(|v| v / scale).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut v = identity(p);
    let mut converged = p == 1;
    for _ in 0..64 {
        let mut maximum = 0_f64;
        for i in 0..p {
            for j in i + 1..p {
                let off = a[i][j];
                maximum = maximum.max(off.abs());
                if off.abs() <= 2e-15 {
                    continue;
                }
                let tau = (a[j][j] - a[i][i]) / (2. * off);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / ((-tau) + tau.hypot(1.))
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                a[i][i] -= t * off;
                a[j][j] += t * off;
                a[i][j] = 0.;
                a[j][i] = 0.;
                for k in 0..p {
                    if k != i && k != j {
                        let x = a[k][i];
                        let y = a[k][j];
                        a[k][i] = c * x - s * y;
                        a[i][k] = a[k][i];
                        a[k][j] = s * x + c * y;
                        a[j][k] = a[k][j];
                    }
                    let x = v[k][i];
                    let y = v[k][j];
                    v[k][i] = c * x - s * y;
                    v[k][j] = s * x + c * y;
                }
            }
        }
        if maximum <= 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("symmetric eigensolver did not converge".into());
    }
    let mut order = (0..p).collect::<Vec<_>>();
    order.sort_by(|i, j| a[*j][*j].total_cmp(&a[*i][*i]).then(i.cmp(j)));
    let values = order
        .iter()
        .map(|j| finite(a[*j][*j] * scale, "eigenvalue"))
        .collect::<Result<Vec<_>, _>>()?;
    let mut vectors = (0..p)
        .map(|i| order.iter().map(|j| v[i][*j]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for j in 0..p {
        let anchor = (0..p)
            .max_by(|i, k| {
                vectors[*i][j]
                    .abs()
                    .total_cmp(&vectors[*k][j].abs())
                    .then(k.cmp(i))
            })
            .unwrap();
        if vectors[anchor][j] < 0. {
            for row in &mut vectors {
                row[j] = -row[j];
            }
        }
    }
    Ok(Eigen { values, vectors })
}
