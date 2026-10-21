//! Bounded CPU learners. Training rows are a prefix; query labels are ignored.
use crate::resampling::Rng;
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
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
struct Node {
    feature: Option<usize>,
    threshold: f64,
    value: f64,
    probabilities: Vec<f64>,
    left: Option<Box<Node>>,
    right: Option<Box<Node>>,
}
impl Node {
    fn leaf<'a>(&'a self, x: &[f64]) -> &'a Self {
        match self.feature {
            None => self,
            Some(j) => {
                if x[j] <= self.threshold {
                    self.left.as_ref().unwrap().leaf(x)
                } else {
                    self.right.as_ref().unwrap().leaf(x)
                }
            }
        }
    }
    fn json(&self) -> String {
        if let Some(feature) = self.feature {
            format!(
                "{{\"feature\":{feature},\"threshold\":{},\"left\":{},\"right\":{}}}",
                self.threshold,
                self.left.as_ref().unwrap().json(),
                self.right.as_ref().unwrap().json()
            )
        } else {
            format!(
                "{{\"value\":{},\"class_probabilities\":{}}}",
                self.value,
                array(&self.probabilities)
            )
        }
    }
}
fn mean(x: impl Iterator<Item = f64>, n: usize) -> f64 {
    x.map(|v| v / n as f64).sum()
}
fn population_scale(y: &[f64]) -> f64 {
    let center = mean(y.iter().copied(), y.len());
    let scale = y.iter().map(|v| v - center).fold(0_f64, f64::hypot) / (y.len() as f64).sqrt();
    if scale > 0. {
        scale
    } else {
        1.
    }
}
fn losses(y: &[f64]) -> Vec<f64> {
    let (mut mu, mut m2) = (0., 0.);
    y.iter()
        .enumerate()
        .map(|(i, &v)| {
            let d = v - mu;
            mu += d / (i + 1) as f64;
            m2 += d * (v - mu);
            m2.max(0.)
        })
        .collect()
}
struct Builder<'a> {
    x: &'a [Vec<f64>],
    y: &'a [f64],
    classes: usize,
    min_leaf: usize,
    max_depth: usize,
    mtry: usize,
    importance: Vec<f64>,
    rng: Rng,
}
impl Builder<'_> {
    fn build(&mut self, ids: &[usize], depth: usize) -> Node {
        let ys: Vec<_> = ids.iter().map(|&i| self.y[i]).collect();
        let value = mean(ys.iter().copied(), ys.len());
        let mut probs = vec![0.; self.classes];
        for &y in &ys {
            if self.classes > 0 {
                probs[y as usize] += 1. / ys.len() as f64;
            }
        }
        let parent = if self.classes == 0 {
            *losses(&ys).last().unwrap()
        } else {
            ys.len() as f64 * (1. - probs.iter().map(|v| v * v).sum::<f64>())
        };
        let mut node = Node {
            feature: None,
            threshold: 0.,
            value,
            probabilities: probs,
            left: None,
            right: None,
        };
        if depth >= self.max_depth || ids.len() < 2 * self.min_leaf || parent <= 0. {
            return node;
        }
        let mut features: Vec<_> = (0..self.x[0].len()).collect();
        self.rng.shuffle(&mut features);
        features.truncate(self.mtry);
        features.sort_unstable();
        let mut best: Option<(f64, usize, f64)> = None;
        for j in features {
            let mut sorted = ids.to_vec();
            sorted.sort_by(|&a, &b| self.x[a][j].total_cmp(&self.x[b][j]).then(a.cmp(&b)));
            let values: Vec<_> = sorted.iter().map(|&i| self.y[i]).collect();
            let prefix = if self.classes == 0 {
                losses(&values)
            } else {
                Vec::new()
            };
            let suffix = if self.classes == 0 {
                let mut v = values.clone();
                v.reverse();
                losses(&v)
            } else {
                Vec::new()
            };
            let mut left_counts = vec![0.; self.classes];
            let mut right_counts = vec![0.; self.classes];
            for &v in &values {
                if self.classes > 0 {
                    right_counts[v as usize] += 1.;
                }
            }
            for split in 1..sorted.len() {
                if self.classes > 0 {
                    left_counts[values[split - 1] as usize] += 1.;
                    right_counts[values[split - 1] as usize] -= 1.;
                }
                if split < self.min_leaf
                    || sorted.len() - split < self.min_leaf
                    || self.x[sorted[split - 1]][j] == self.x[sorted[split]][j]
                {
                    continue;
                }
                let loss = if self.classes == 0 {
                    prefix[split - 1] + suffix[sorted.len() - split - 1]
                } else {
                    let a = split as f64;
                    let b = (sorted.len() - split) as f64;
                    a - left_counts.iter().map(|v| v * v / a).sum::<f64>() + b
                        - right_counts.iter().map(|v| v * v / b).sum::<f64>()
                };
                let gain = parent - loss;
                let low = self.x[sorted[split - 1]][j];
                let high = self.x[sorted[split]][j];
                let mut threshold = low + (high - low) / 2.;
                if threshold >= high {
                    threshold = low;
                }
                if gain > parent.abs() * 1e-14
                    && best
                        .as_ref()
                        .is_none_or(|b| gain > b.0 + parent.abs() * 1e-14)
                {
                    best = Some((gain, j, threshold));
                }
            }
        }
        if let Some((gain, feature, threshold)) = best {
            let (left, right): (Vec<_>, Vec<_>) = ids
                .iter()
                .copied()
                .partition(|&i| self.x[i][feature] <= threshold);
            self.importance[feature] += gain;
            node.feature = Some(feature);
            node.threshold = threshold;
            node.left = Some(Box::new(self.build(&left, depth + 1)));
            node.right = Some(Box::new(self.build(&right, depth + 1)));
        }
        node
    }
}
fn dimensions(data: &[Vec<f64>], n: f64) -> Result<(usize, usize), String> {
    let train = integer(n, data.len())?;
    let p = data[0].len().saturating_sub(1);
    if train < 2 || p == 0 || p > 64 {
        return Err("learners need >=2 training rows and1..64 features".into());
    }
    Ok((train, p))
}
fn labels(
    data: &[Vec<f64>],
    n: usize,
    p: usize,
    task: usize,
) -> Result<(Vec<f64>, Vec<f64>, f64, f64), String> {
    if task == 0 {
        let origin = data[0][p];
        let y: Vec<_> = data[..n].iter().map(|r| r[p] - origin).collect();
        let scale = population_scale(&y);
        return Ok((
            y.iter().map(|v| v / scale).collect(),
            Vec::new(),
            origin,
            scale,
        ));
    }
    let mut classes: Vec<_> = data[..n].iter().map(|r| r[p]).collect();
    for &v in &classes {
        integer(v, 1_000_000_000)?;
    }
    classes.sort_by(f64::total_cmp);
    classes.dedup();
    if classes.len() < 2 || classes.len() > 64 {
        return Err("classification needs2..64 observed classes".into());
    }
    let y = data[..n]
        .iter()
        .map(|r| classes.binary_search_by(|v| v.total_cmp(&r[p])).unwrap() as f64)
        .collect();
    Ok((y, classes, 0., 1.))
}
fn trees(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let is_forest = op == "random_forest";
    let expected = if is_forest { 6 } else { 4 };
    if args.len() != expected {
        return Err("tree params [n_train,depth,min_leaf,task]; forest [n_train,trees,depth,min_leaf,task,seed]".into());
    }
    let (n, p) = dimensions(data, args[0])?;
    let count = if is_forest { integer(args[1], 200)? } else { 1 };
    let d = usize::from(is_forest);
    let depth = integer(args[1 + d], 16)?;
    let min_leaf = integer(args[2 + d], n)?;
    let task = integer(args[3 + d], 1)?;
    let seed = if is_forest { args[5] } else { 0. };
    if depth == 0
        || min_leaf == 0
        || count == 0
        || count
            .saturating_mul(n)
            .saturating_mul(p)
            .saturating_mul(depth)
            > 100_000_000
        || count.saturating_mul(n) > 100000
    {
        return Err("invalid or excessive tree work".into());
    }
    let (y, classes, origin, target_scale) = labels(data, n, p, task)?;
    if count.saturating_mul(n).saturating_mul(classes.len().max(1)) > 1_000_000 {
        return Err("forest model size exceeds limit".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[..p].to_vec()).collect();
    let mtry = if is_forest {
        if task == 1 {
            (p as f64).sqrt().ceil() as usize
        } else {
            (p / 3).max(1)
        }
    } else {
        p
    };
    let mut builder = Builder {
        x: &x,
        y: &y,
        classes: classes.len(),
        min_leaf,
        max_depth: depth,
        mtry,
        importance: vec![0.; p],
        rng: Rng::new(seed)?,
    };
    let mut models = Vec::new();
    for _ in 0..count {
        let ids: Vec<usize> = if is_forest {
            (0..n).map(|_| builder.rng.index(n)).collect()
        } else {
            (0..n).collect()
        };
        models.push(builder.build(&ids, 0));
    }
    let mut predictions = Vec::new();
    let mut posterior = Vec::new();
    for row in &x[n..] {
        if task == 0 {
            predictions.push(
                origin
                    + target_scale
                        * models
                            .iter()
                            .map(|m| m.leaf(row).value / count as f64)
                            .sum::<f64>(),
            );
        } else {
            let mut prob = vec![0.; classes.len()];
            for model in &models {
                for (j, &v) in model.leaf(row).probabilities.iter().enumerate() {
                    prob[j] += v / count as f64;
                }
            }
            let class = (0..prob.len())
                .max_by(|&a, &b| prob[a].total_cmp(&prob[b]).then(b.cmp(&a)))
                .unwrap();
            predictions.push(classes[class]);
            posterior.push(array(&prob));
        }
    }
    let total = builder.importance.iter().sum::<f64>();
    let importance: Vec<_> = builder
        .importance
        .iter()
        .map(|v| if total > 0. { v / total } else { 0. })
        .collect();
    Ok(format!("{{\"version\":1,\"algorithm\":\"{op}\",\"n_train\":{n},\"task\":{task},\"classes\":{},\"target_origin\":{origin},\"target_scale\":{target_scale},\"predictions\":{},\"posterior_probabilities\":[{}],\"feature_importance\":{},\"trees\":[{}],\"seed\":{seed}}}",array(&classes),array(&predictions),posterior.join(","),array(&importance),models.iter().map(Node::json).collect::<Vec<_>>().join(",")))
}
fn boosting(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 7 {
        return Err("boosting requires [n_train,iterations,learning_rate,depth,min_leaf,task0squared/1binary,seed]".into());
    }
    let (n, p) = dimensions(data, args[0])?;
    let count = integer(args[1], 500)?;
    let rate = args[2];
    let depth = integer(args[3], 8)?;
    let min_leaf = integer(args[4], n)?;
    let task = integer(args[5], 1)?;
    Rng::new(args[6])?;
    if count == 0
        || rate <= 0.
        || rate > 1.
        || depth == 0
        || min_leaf == 0
        || count
            .saturating_mul(n)
            .saturating_mul(p)
            .saturating_mul(depth)
            > 100_000_000
        || count.saturating_mul(n) > 100000
    {
        return Err("invalid or excessive boosting work".into());
    }
    let origin = if task == 0 { data[0][p] } else { 0. };
    let y: Vec<_> = data[..n].iter().map(|r| r[p] - origin).collect();
    let target_scale = if task == 0 { population_scale(&y) } else { 1. };
    let y: Vec<_> = y.iter().map(|v| v / target_scale).collect();
    if task == 1 && (!y.iter().all(|&v| v == 0. || v == 1.) || y.iter().all(|&v| v == y[0])) {
        return Err("binary boosting needs both0/1 outcomes".into());
    }
    let x: Vec<_> = data.iter().map(|r| r[..p].to_vec()).collect();
    let base_mean = mean(y.iter().copied(), n);
    let base = if task == 0 {
        base_mean
    } else {
        (base_mean / (1. - base_mean)).ln()
    };
    let mut scores = vec![base; data.len()];
    let mut models = Vec::new();
    let mut history = Vec::new();
    let mut importance = vec![0.; p];
    for _ in 0..count {
        let residual: Vec<_> = y
            .iter()
            .enumerate()
            .map(|(i, &y)| {
                y - if task == 0 {
                    scores[i]
                } else {
                    sigmoid(scores[i])
                }
            })
            .collect();
        let mut builder = Builder {
            x: &x,
            y: &residual,
            classes: 0,
            min_leaf,
            max_depth: depth,
            mtry: p,
            importance: vec![0.; p],
            rng: Rng::new(args[6])?,
        };
        let model = builder.build(&(0..n).collect::<Vec<_>>(), 0);
        for (i, row) in x.iter().enumerate() {
            scores[i] += rate * model.leaf(row).value;
        }
        for (j, v) in builder.importance.iter().enumerate() {
            importance[j] += v;
        }
        let loss = if task == 0 {
            y.iter()
                .enumerate()
                .map(|(i, &y)| (y - scores[i]).powi(2) / n as f64)
                .sum()
        } else {
            y.iter()
                .enumerate()
                .map(|(i, &y)| {
                    let z = scores[i];
                    (z.max(0.) - y * z + (-z.abs()).exp().ln_1p()) / n as f64
                })
                .sum()
        };
        history.push(loss);
        models.push(model);
    }
    let predictions: Vec<_> = scores[n..]
        .iter()
        .map(|&v| {
            if task == 0 {
                origin + target_scale * v
            } else {
                if v >= 0. {
                    1.
                } else {
                    0.
                }
            }
        })
        .collect();
    let probabilities: Vec<_> = if task == 1 {
        scores[n..].iter().map(|&v| sigmoid(v)).collect()
    } else {
        Vec::new()
    };
    let total = importance.iter().sum::<f64>();
    if total > 0. {
        for v in &mut importance {
            *v /= total;
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"gradient_boosted_regression_trees\",\"loss\":\"{}\",\"loss_units\":\"standardized_training_target\",\"base_score\":{base},\"target_origin\":{origin},\"target_scale\":{target_scale},\"learning_rate\":{rate},\"predictions\":{},\"probabilities\":{},\"training_loss\":{},\"feature_importance\":{},\"trees\":[{}]}}",if task==0{"squared_error"}else{"binomial_gradient"},array(&predictions),array(&probabilities),array(&history),array(&importance),models.iter().map(Node::json).collect::<Vec<_>>().join(",")))
}
fn mlp(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !(6..=7).contains(&args.len()) {
        return Err("MLP needs [n_train,hidden,epochs,learning_rate,task0regression/1binary,seed,tolerance=1e-8]".into());
    }
    let (n, p) = dimensions(data, args[0])?;
    let h = integer(args[1], 128)?;
    let epochs = integer(args[2], 10000)?;
    let rate = args[3];
    let task = integer(args[4], 1)?;
    let mut rng = Rng::new(args[5])?;
    let tolerance = *args.get(6).unwrap_or(&1e-8);
    if h == 0
        || epochs == 0
        || rate <= 0.
        || rate > 1.
        || tolerance <= 0.
        || tolerance > 0.1
        || n.saturating_mul(p).saturating_mul(h).saturating_mul(epochs) > 100_000_000
        || data.len().saturating_mul(p).saturating_mul(h) > 100_000_000
    {
        return Err("invalid MLP settings or excessive work".into());
    }
    let origins = data[0][..p].to_vec();
    let mut offsets = vec![0.; p];
    let mut scales = vec![1.; p];
    for j in 0..p {
        offsets[j] = mean(data[..n].iter().map(|r| r[j] - origins[j]), n);
        let deviation = data[..n]
            .iter()
            .map(|r| (r[j] - origins[j]) - offsets[j])
            .fold(0_f64, f64::hypot)
            / (n as f64).sqrt();
        if deviation > 0. {
            scales[j] = deviation;
        }
    }
    let x: Vec<Vec<f64>> = data
        .iter()
        .map(|r| {
            (0..p)
                .map(|j| ((r[j] - origins[j]) - offsets[j]) / scales[j])
                .collect()
        })
        .collect();
    let target_origin = if task == 0 { data[0][p] } else { 0. };
    let target_offset = if task == 0 {
        mean(data[..n].iter().map(|r| r[p] - target_origin), n)
    } else {
        0.
    };
    let deviation = if task == 0 {
        data[..n]
            .iter()
            .map(|r| (r[p] - target_origin) - target_offset)
            .fold(0_f64, f64::hypot)
            / (n as f64).sqrt()
    } else {
        1.
    };
    let target_scale = if deviation > 0. { deviation } else { 1. };
    let y: Vec<_> = data[..n]
        .iter()
        .map(|r| ((r[p] - target_origin) - target_offset) / target_scale)
        .collect();
    if task == 1 && (!y.iter().all(|&v| v == 0. || v == 1.) || y.iter().all(|&v| v == y[0])) {
        return Err("binary MLP needs both0/1 outcomes".into());
    }
    let b1 = h * p;
    let w2 = b1 + h;
    let b2 = w2 + h;
    let mut weights = vec![0.; b2 + 1];
    let a = (6. / (p + h) as f64).sqrt();
    for w in &mut weights[..b1] {
        *w = (2. * rng.uniform() - 1.) * a;
    }
    let a = (6. / (h + 1) as f64).sqrt();
    for w in &mut weights[w2..b2] {
        *w = (2. * rng.uniform() - 1.) * a;
    }
    let mut first = vec![0.; weights.len()];
    let mut second = vec![0.; weights.len()];
    let mut iterations = 0;
    let mut converged = false;
    let mut training_loss: f64;
    for epoch in 1..=epochs {
        let mut gradient = vec![0.; weights.len()];
        training_loss = 0.;
        for (i, row) in x[..n].iter().enumerate() {
            let hidden: Vec<_> = (0..h)
                .map(|k| {
                    (weights[b1 + k] + (0..p).map(|j| weights[k * p + j] * row[j]).sum::<f64>())
                        .tanh()
                })
                .collect();
            let score = weights[b2] + (0..h).map(|k| weights[w2 + k] * hidden[k]).sum::<f64>();
            let prediction = if task == 0 { score } else { sigmoid(score) };
            let g = (prediction - y[i]) / n as f64;
            training_loss += if task == 0 {
                0.5 * (prediction - y[i]).powi(2) / n as f64
            } else {
                (score.max(0.) - y[i] * score + (-score.abs()).exp().ln_1p()) / n as f64
            };
            gradient[b2] += g;
            for k in 0..h {
                gradient[w2 + k] += g * hidden[k];
                let local = g * weights[w2 + k] * (1. - hidden[k] * hidden[k]);
                gradient[b1 + k] += local;
                for j in 0..p {
                    gradient[k * p + j] += local * row[j];
                }
            }
        }
        iterations = epoch;
        if gradient.iter().any(|v| !v.is_finite()) || !training_loss.is_finite() {
            return Err("MLP numerical range exceeded".into());
        }
        let grad_norm = gradient.iter().fold(0_f64, |n, &v| n.hypot(v));
        if grad_norm < tolerance {
            converged = true;
            break;
        }
        for j in 0..weights.len() {
            first[j] = 0.9 * first[j] + 0.1 * gradient[j];
            second[j] = 0.999 * second[j] + 0.001 * gradient[j] * gradient[j];
            let m = first[j] / (1. - 0.9_f64.powi(epoch as i32));
            let v = second[j] / (1. - 0.999_f64.powi(epoch as i32));
            weights[j] -= rate * m / (v.sqrt() + 1e-8);
        }
    }
    // Diagnostics describe the returned weights, including the final update.
    training_loss = x[..n]
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let hidden: Vec<_> = (0..h)
                .map(|k| {
                    (weights[b1 + k] + (0..p).map(|j| weights[k * p + j] * row[j]).sum::<f64>())
                        .tanh()
                })
                .collect();
            let score = weights[b2] + (0..h).map(|k| weights[w2 + k] * hidden[k]).sum::<f64>();
            if task == 0 {
                0.5 * (score - y[i]).powi(2) / n as f64
            } else {
                (score.max(0.) - y[i] * score + (-score.abs()).exp().ln_1p()) / n as f64
            }
        })
        .sum();
    let mut predictions = Vec::new();
    let mut probabilities = Vec::new();
    for row in &x[n..] {
        let hidden: Vec<_> = (0..h)
            .map(|k| {
                (weights[b1 + k] + (0..p).map(|j| weights[k * p + j] * row[j]).sum::<f64>()).tanh()
            })
            .collect();
        let score = weights[b2] + (0..h).map(|k| weights[w2 + k] * hidden[k]).sum::<f64>();
        if task == 0 {
            predictions.push(target_origin + (target_offset + target_scale * score));
        } else {
            let prob = sigmoid(score);
            probabilities.push(prob);
            predictions.push(if prob >= 0.5 { 1. } else { 0. });
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"cpu_single_hidden_tanh_adam\",\"task\":{task},\"n_train\":{n},\"n_features\":{p},\"hidden_units\":{h},\"iterations\":{iterations},\"converged\":{converged},\"training_loss\":{training_loss},\"predictions\":{},\"probabilities\":{},\"feature_origins\":{},\"feature_mean_offsets\":{},\"feature_scales\":{},\"target_origin\":{target_origin},\"target_mean_offset\":{target_offset},\"target_scale\":{target_scale},\"weights_input_hidden\":{},\"bias_hidden\":{},\"weights_hidden_output\":{},\"bias_output\":{},\"seed\":{}}}",array(&predictions),array(&probabilities),array(&origins),array(&offsets),array(&scales),array(&weights[..b1]),array(&weights[b1..w2]),array(&weights[w2..b2]),weights[b2],args[5]))
}

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width < 2
        || data.len() > 100000
        || data.len().saturating_mul(width) > 1000000
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite() || v.abs() > 1e100))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("learning requires bounded finite rectangular data".into());
    }
    let out = match operation {
        "decision_tree" | "random_forest" => trees(operation, data, args),
        "gradient_boosting" => boosting(data, args),
        "mlp" => mlp(data, args),
        _ => Err("learning method not implemented".into()),
    }?;
    if out
        .split([':', '[', ']', ',', '}', ' '])
        .any(|token| matches!(token, "inf" | "-inf" | "NaN"))
    {
        Err("learner exceeded finite numeric range".into())
    } else {
        Ok(out)
    }
}
