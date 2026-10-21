//! Prediction evaluation and leakage-conscious partition indices.
use crate::resampling::Rng;
use std::collections::BTreeMap;

fn num(x: f64) -> String {
    if x.is_finite() {
        x.to_string()
    } else {
        "null".into()
    }
}
fn array(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(|&v| num(v)).collect::<Vec<_>>().join(",")
    )
}
fn indices(x: &[usize]) -> String {
    format!(
        "[{}]",
        x.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|v| array(v)).collect::<Vec<_>>().join(",")
    )
}
fn index_matrix(x: &[Vec<usize>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|v| indices(v)).collect::<Vec<_>>().join(",")
    )
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded nonnegative integer".into())
    }
}
fn label(x: f64) -> Result<i64, String> {
    if x.abs() <= 9_007_199_254_740_991. && x.fract() == 0. {
        Ok(x as i64)
    } else {
        Err("labels must be exact integers".into())
    }
}
fn ratio(a: f64, b: f64) -> f64 {
    if b > 0. {
        a / b
    } else {
        f64::NAN
    }
}
fn json(op: &str, body: String) -> String {
    format!("{{\"version\":1,\"method\":\"{op}\",{body}}}")
}
fn weights(data: &[Vec<f64>], base: usize) -> Result<(Vec<f64>, f64), String> {
    if ![base, base + 1].contains(&data[0].len()) {
        return Err("wrong observation width".into());
    }
    let w: Vec<_> = data
        .iter()
        .map(|r| if r.len() == base { 1. } else { r[base] })
        .collect();
    let total = w.iter().sum::<f64>();
    if w.iter().any(|&v| v < 0.) || total <= 0. || !total.is_finite() {
        return Err("weights require a finite positive total".into());
    }
    Ok((w, total))
}
fn regression_values(data: &[Vec<f64>]) -> Result<(Vec<f64>, [f64; 6]), String> {
    let (w, total) = weights(data, 2)?;
    if data
        .iter()
        .any(|r| r[0].abs() > 1e150 || r[1].abs() > 1e150)
    {
        return Err("regression metric observations exceed numerical bounds".into());
    }
    let origin = data[0][0];
    let offset = data
        .iter()
        .zip(&w)
        .map(|(r, &v)| (v / total) * (r[0] - origin))
        .sum::<f64>();
    let residuals: Vec<_> = data.iter().map(|r| r[0] - r[1]).collect();
    let mut rmse = 0_f64;
    let mut target_sd = 0_f64;
    let mut mae = 0.;
    let mut bias = 0.;
    for ((r, &e), &v) in data.iter().zip(&residuals).zip(&w) {
        let q = v / total;
        rmse = rmse.hypot(q.sqrt() * e);
        target_sd = target_sd.hypot(q.sqrt() * ((r[0] - origin) - offset));
        mae += q * e.abs();
        bias += q * e;
    }
    let r2 = if target_sd > 0. {
        1. - (rmse / target_sd).powi(2)
    } else {
        f64::NAN
    };
    Ok((residuals, [rmse, mae, bias, r2, origin + offset, total]))
}
fn regression(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !args.is_empty() {
        return Err("regression_metrics takes no parameters".into());
    }
    let (residuals, v) = regression_values(data)?;
    Ok(json("regression_metrics",format!("\"residuals\":{},\"rmse\":{},\"mse\":{},\"mae\":{},\"bias\":{},\"r_squared\":{},\"target_mean\":{},\"weight_sum\":{},\"bias_convention\":\"observed_minus_prediction\"",array(&residuals),num(v[0]),num(v[0]*v[0]),num(v[1]),num(v[2]),num(v[3]),num(v[4]),num(v[5]))))
}
fn confusion(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !args.is_empty() {
        return Err("confusion_matrix takes no parameters".into());
    }
    let (w, total) = weights(data, 2)?;
    let mut labels = Vec::new();
    for (r, &v) in data.iter().zip(&w) {
        let a = label(r[0])?;
        let b = label(r[1])?;
        if v > 0. {
            labels.push(a);
            labels.push(b);
        }
    }
    labels.sort_unstable();
    labels.dedup();
    if labels.len() > 256 {
        return Err("too many classification levels".into());
    }
    let k = labels.len();
    let mut m = vec![vec![0.; k]; k];
    for (r, &v) in data.iter().zip(&w) {
        if v > 0. {
            let a = labels.binary_search(&label(r[0])?).unwrap();
            let b = labels.binary_search(&label(r[1])?).unwrap();
            m[a][b] += v;
        }
    }
    let actual: Vec<_> = m.iter().map(|r| r.iter().sum::<f64>()).collect();
    let predicted: Vec<_> = (0..k)
        .map(|j| m.iter().map(|r| r[j]).sum::<f64>())
        .collect();
    let precision: Vec<_> = (0..k).map(|j| ratio(m[j][j], predicted[j])).collect();
    let recall: Vec<_> = (0..k).map(|j| ratio(m[j][j], actual[j])).collect();
    let f1: Vec<_> = (0..k)
        .map(|j| {
            let scale = actual[j].max(predicted[j]);
            2. * (m[j][j] / scale) / (actual[j] / scale + predicted[j] / scale)
        })
        .collect();
    let accuracy = (0..k).map(|j| m[j][j] / total).sum::<f64>();
    let balanced = recall.iter().filter(|v| v.is_finite()).sum::<f64>()
        / actual.iter().filter(|&&v| v > 0.).count() as f64;
    Ok(json("confusion_matrix",format!("\"labels\":{},\"matrix\":{},\"actual_weight\":{},\"predicted_weight\":{},\"precision\":{},\"recall\":{},\"f1\":{},\"accuracy\":{},\"balanced_accuracy\":{},\"macro_f1\":{},\"weight_sum\":{},\"matrix_orientation\":\"actual_rows_predicted_columns\"",array(&labels.iter().map(|&v|v as f64).collect::<Vec<_>>()),matrix(&m),array(&actual),array(&predicted),array(&precision),array(&recall),array(&f1),num(accuracy),num(balanced),num(f1.iter().sum::<f64>()/k as f64),num(total))))
}
fn binary_rows(data: &[Vec<f64>], probability: bool) -> Result<(Vec<f64>, f64), String> {
    let (w, total) = weights(data, 2)?;
    if data
        .iter()
        .any(|r| ![0., 1.].contains(&r[0]) || (probability && !(0. ..=1.).contains(&r[1])))
    {
        return Err("requires binary outcomes and valid probabilities/scores".into());
    }
    Ok((w, total))
}
fn binary_f1(tp: f64, fp: f64, fn_: f64) -> f64 {
    let scale = tp.max(fp).max(fn_);
    if scale == 0. {
        f64::NAN
    } else {
        2. * (tp / scale) / (2. * (tp / scale) + fp / scale + fn_ / scale)
    }
}
fn binary(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 2 {
        return Err("binary_metrics expects [threshold=.5,clip=1e-15]".into());
    }
    let threshold = *args.first().unwrap_or(&0.5);
    let clip = *args.get(1).unwrap_or(&1e-15);
    if !(0. ..=1.).contains(&threshold) || !(1e-16..=1e-2).contains(&clip) {
        return Err("invalid threshold or probability clip".into());
    }
    let (w, total) = binary_rows(data, true)?;
    let mut counts = [0.; 4];
    let mut brier = 0.;
    let mut loss = 0.;
    for (r, &v) in data.iter().zip(&w) {
        let y = r[0];
        let p = r[1];
        let predicted = usize::from(p >= threshold);
        counts[(y as usize) * 2 + predicted] += v;
        let q = p.clamp(clip, 1. - clip);
        let a = v / total;
        brier += a * (p - y).powi(2);
        loss -= a * if y == 1. { q.ln() } else { (-q).ln_1p() };
    }
    let [tn, fp, fn_, tp] = counts;
    Ok(json("binary_metrics",format!("\"counts\":{},\"count_order\":\"TN_FP_FN_TP\",\"threshold\":{},\"accuracy\":{},\"precision\":{},\"recall\":{},\"specificity\":{},\"f1\":{},\"brier_score\":{},\"log_loss\":{},\"probability_clip\":{},\"weight_sum\":{}",array(&[tn,fp,fn_,tp]),num(threshold),num((tp+tn)/total),num(ratio(tp,tp+fp)),num(ratio(tp,tp+fn_)),num(ratio(tn,tn+fp)),num(binary_f1(tp,fp,fn_)),num(brier),num(loss),num(clip),num(total))))
}
fn curve(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if !args.is_empty() {
        return Err("prediction curves take no parameters".into());
    }
    let (w, total) = binary_rows(data, false)?;
    let positives = data.iter().zip(&w).map(|(r, &v)| r[0] * v).sum::<f64>();
    let negatives = data
        .iter()
        .zip(&w)
        .map(|(r, &v)| (1. - r[0]) * v)
        .sum::<f64>();
    if positives <= 0. || negatives <= 0. {
        return Err("curves require positive weight in both classes".into());
    }
    let mut sorted: Vec<_> = (0..data.len()).filter(|&i| w[i] > 0.).collect();
    sorted.sort_by(|&i, &j| data[j][1].total_cmp(&data[i][1]).then(i.cmp(&j)));
    let mut thresholds = vec![f64::NAN];
    let mut tpr = vec![0.];
    let mut fpr = vec![0.];
    let mut selection = vec![0.];
    let mut precision = vec![1.];
    let mut lift = vec![f64::NAN];
    let (mut tp, mut fp) = (0., 0.);
    let mut i = 0;
    while i < sorted.len() {
        let threshold = data[sorted[i]][1];
        let mut j = i;
        while j < sorted.len() && data[sorted[j]][1] == threshold {
            let row = sorted[j];
            if data[row][0] == 1. {
                tp += w[row];
            } else {
                fp += w[row];
            }
            j += 1;
        }
        thresholds.push(threshold);
        tpr.push((tp / positives).min(1.));
        fpr.push((fp / negatives).min(1.));
        selection.push(((tp + fp) / total).min(1.));
        precision.push(tp / (tp + fp));
        lift.push((tp / (tp + fp)) / (positives / total));
        i = j;
    }
    let mut roc = 0.;
    let mut ap = 0.;
    let mut pr = 0.;
    for k in 1..tpr.len() {
        roc += (fpr[k] - fpr[k - 1]) * (tpr[k] + tpr[k - 1]) / 2.;
        ap += (tpr[k] - tpr[k - 1]) * precision[k];
        pr += (tpr[k] - tpr[k - 1]) * (precision[k] + precision[k - 1]) / 2.;
    }
    let body=match op {
        "roc_curve"=>format!("\"false_positive_rate\":{},\"true_positive_rate\":{},\"auc\":{}",array(&fpr),array(&tpr),num(roc)),
        "precision_recall_curve"=>format!("\"precision\":{},\"recall\":{},\"average_precision\":{},\"trapezoid_area\":{}",array(&precision),array(&tpr),num(ap),num(pr)),
        _=>format!("\"selection_fraction\":{},\"cumulative_gain\":{},\"response_rate\":{},\"lift\":{},\"population_response_rate\":{}",array(&selection),array(&tpr),array(&precision.iter().enumerate().map(|(i,&v)|if i==0 {f64::NAN}else{v}).collect::<Vec<_>>()),array(&lift),num(positives/total)),
    };
    Ok(json(op,format!("\"thresholds\":{},\"threshold_policy\":\"score_greater_or_equal_grouped_ties\",\"weight_sum\":{},{}",array(&thresholds),num(total),body)))
}
fn calibration(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() > 1 {
        return Err("calibration_curve expects [bins=10]".into());
    }
    let k = integer(*args.first().unwrap_or(&10.), 1000)?;
    if k == 0 {
        return Err("bins must be positive".into());
    }
    let (w, total) = binary_rows(data, true)?;
    let mut mass = vec![0.; k];
    let mut sp = vec![0.; k];
    let mut sy = vec![0.; k];
    // Normalize within bins before averaging so arbitrary finite weights are safe.
    for (r, &v) in data.iter().zip(&w) {
        let b = ((r[1] * k as f64) as usize).min(k - 1);
        mass[b] += v;
    }
    for (r, &v) in data.iter().zip(&w) {
        let b = ((r[1] * k as f64) as usize).min(k - 1);
        if mass[b] > 0. {
            sp[b] += (v / mass[b]) * r[1];
            sy[b] += (v / mass[b]) * r[0];
        }
    }
    let ece = (0..k)
        .map(|i| (mass[i] / total) * (sp[i] - sy[i]).abs())
        .sum::<f64>();
    let max_error = (0..k)
        .filter(|&i| mass[i] > 0.)
        .map(|i| (sp[i] - sy[i]).abs())
        .fold(0., f64::max);
    for i in 0..k {
        if mass[i] == 0. {
            sp[i] = f64::NAN;
            sy[i] = f64::NAN;
        }
    }
    Ok(json("calibration_curve",format!("\"bin_edges\":{},\"bin_weight\":{},\"mean_probability\":{},\"observed_fraction\":{},\"expected_calibration_error\":{},\"maximum_calibration_error\":{},\"weight_sum\":{}",array(&(0..=k).map(|i|i as f64/k as f64).collect::<Vec<_>>()),array(&mass),array(&sp),array(&sy),num(ece),num(max_error),num(total))))
}
fn profit(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 4 {
        return Err(
            "profit_decision requires row-major [profit00,profit01,profit10,profit11]".into(),
        );
    }
    let (w, total) = binary_rows(data, true)?;
    let mut decisions = Vec::new();
    let mut expected = Vec::new();
    let mut actual = Vec::new();
    let mut mean_actual = 0.;
    let mut mean_expected = 0.;
    for (r, &v) in data.iter().zip(&w) {
        let e0 = (1. - r[1]) * args[0] + r[1] * args[2];
        let e1 = (1. - r[1]) * args[1] + r[1] * args[3];
        if !e0.is_finite() || !e1.is_finite() {
            return Err("profit arithmetic overflow".into());
        }
        let decision = usize::from(e1 > e0);
        let e = if decision == 0 { e0 } else { e1 };
        let a = args[(r[0] as usize) * 2 + decision];
        decisions.push(decision);
        expected.push(e);
        actual.push(a);
        mean_actual += (v / total) * a;
        mean_expected += (v / total) * e;
    }
    Ok(json("profit_decision",format!("\"decisions\":{},\"actual_profit\":{},\"expected_profit\":{},\"mean_actual_profit\":{},\"total_actual_profit\":{},\"mean_expected_profit\":{},\"total_expected_profit\":{},\"weight_sum\":{},\"matrix_orientation\":\"actual_rows_decision_columns\"",indices(&decisions),array(&actual),array(&expected),num(mean_actual),num(mean_actual*total),num(mean_expected),num(mean_expected*total),num(total))))
}
fn ensemble(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data[0].len();
    if p > 1024 || (!args.is_empty() && args.len() != p) {
        return Err("ensemble weights must match member columns".into());
    }
    let mut w = if args.is_empty() {
        vec![1.; p]
    } else {
        args.to_vec()
    };
    let max = w.iter().copied().fold(0., f64::max);
    if max <= 0. || w.iter().any(|&v| v < 0.) {
        return Err("ensemble needs nonnegative weights and positive total".into());
    }
    for v in &mut w {
        *v /= max;
    }
    let total = w.iter().sum::<f64>();
    for v in &mut w {
        *v /= total;
    }
    let predictions: Vec<_> = data
        .iter()
        .map(|r| r[0] + r.iter().zip(&w).map(|(&x, &q)| q * (x - r[0])).sum::<f64>())
        .collect();
    if predictions.iter().any(|v| !v.is_finite()) {
        return Err("ensemble arithmetic overflow".into());
    }
    Ok(json(
        "ensemble_predict",
        format!(
            "\"predictions\":{},\"normalized_weights\":{}",
            array(&predictions),
            array(&w)
        ),
    ))
}
fn compare(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() < 2 || data[0].len() > 65 || args.len() > 2 {
        return Err("model_compare expects outcome plus <=64 prediction columns and [task0regression/1binary,metric0RMSEorLogLoss/1MAEorBrier]".into());
    }
    let task = integer(*args.first().unwrap_or(&0.), 1)?;
    let metric = integer(*args.get(1).unwrap_or(&0.), 1)?;
    let mut scores = Vec::new();
    for j in 1..data[0].len() {
        let rows: Vec<_> = data.iter().map(|r| vec![r[0], r[j]]).collect();
        let score = if task == 0 {
            let (_, v) = regression_values(&rows)?;
            v[metric]
        } else {
            binary_rows(&rows, true)?;
            rows.iter()
                .map(|r| {
                    if metric == 0 {
                        let p = r[1].clamp(1e-15, 1. - 1e-15);
                        -if r[0] == 1. { p.ln() } else { (-p).ln_1p() }
                    } else {
                        (r[1] - r[0]).powi(2)
                    }
                })
                .sum::<f64>()
                / rows.len() as f64
        };
        scores.push(score);
    }
    let mut ranked: Vec<_> = (0..scores.len()).collect();
    ranked.sort_by(|&i, &j| scores[i].total_cmp(&scores[j]).then(i.cmp(&j)));
    let name = match (task, metric) {
        (0, 0) => "rmse",
        (0, _) => "mae",
        (_, 0) => "clipped_log_loss",
        _ => "brier_score",
    };
    Ok(json("model_compare",format!("\"scores\":{},\"ranked_model_indices\":{},\"selection_metric\":\"{}\",\"evaluation_scope\":\"caller_supplied_predictions_no_refit\"",array(&scores),indices(&ranked),name)))
}
fn holdout(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.len() < 2 || args.len() < 2 || args.len() > 3 {
        return Err("holdout_split requires [train_fraction,seed,stratified0/1=0]".into());
    }
    let frac = args[0];
    if !(0. ..1.).contains(&frac) || frac == 0. {
        return Err("train fraction must lie strictly between zero and one".into());
    }
    let stratified = integer(*args.get(2).unwrap_or(&0.), 1)? == 1;
    let mut rng = Rng::new(args[1])?;
    let mut groups = BTreeMap::<i64, Vec<usize>>::new();
    if stratified && data[0].len() != 1 {
        return Err("stratified split requires scalar integer labels".into());
    }
    for (i, r) in data.iter().enumerate() {
        groups
            .entry(if stratified { label(r[0])? } else { 0 })
            .or_default()
            .push(i);
    }
    let mut train = Vec::new();
    let mut test = Vec::new();
    for group in groups.values_mut() {
        if group.len() < 2 {
            return Err("every stratum needs >=2 observations".into());
        }
        rng.shuffle(group);
        let n = ((frac * group.len() as f64).floor() as usize).clamp(1, group.len() - 1);
        train.extend_from_slice(&group[..n]);
        test.extend_from_slice(&group[n..]);
    }
    train.sort_unstable();
    test.sort_unstable();
    let mut fold = vec![1; data.len()];
    for &i in &train {
        fold[i] = 0;
    }
    Ok(json("holdout_split",format!("\"train_indices\":{},\"test_indices\":{},\"partition\":{},\"train_fraction_actual\":{},\"seed\":{},\"stratified\":{},\"rounding_policy\":\"per_group_floor_clamped_to_keep_both_nonempty\"",indices(&train),indices(&test),indices(&fold),num(train.len() as f64/data.len() as f64),num(args[1]),stratified)))
}
fn rolling(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 2 || args.len() > 4 {
        return Err(
            "rolling_origin requires [initial_train,horizon,step=1,window0expanding=0]".into(),
        );
    }
    let n = data.len();
    let initial = integer(args[0], n)?;
    let horizon = integer(args[1], n)?;
    let step = integer(*args.get(2).unwrap_or(&1.), n)?;
    let window = integer(*args.get(3).unwrap_or(&0.), n)?;
    if initial == 0 || horizon == 0 || step == 0 || initial + horizon > n {
        return Err("invalid rolling origin bounds".into());
    }
    let mut origins = Vec::new();
    let mut train = Vec::new();
    let mut test = Vec::new();
    let mut count = 0;
    for origin in (initial..=n - horizon).step_by(step) {
        let start = if window == 0 {
            0
        } else {
            origin.saturating_sub(window)
        };
        count += origin - start + horizon;
        if count > 1_000_000 {
            return Err("rolling split output exceeds one million indices".into());
        }
        origins.push(origin);
        train.push((start..origin).collect::<Vec<_>>());
        test.push((origin..origin + horizon).collect::<Vec<_>>());
    }
    Ok(json(
        "rolling_origin",
        format!(
            "\"origins\":{},\"train_indices\":{},\"test_indices\":{},\"horizon\":{},\"window\":{}",
            indices(&origins),
            index_matrix(&train),
            index_matrix(&test),
            horizon,
            window
        ),
    ))
}
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data.is_empty()
        || data[0].is_empty()
        || data.len() > 100_000
        || data[0].len() > 1024
        || data.len().saturating_mul(data[0].len()) > 1_000_000
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != data[0].len() || r.iter().any(|v| !v.is_finite()))
    {
        return Err("invalid or excessive evaluation input".into());
    }
    match op {
        "regression_metrics" => regression(data, args),
        "confusion_matrix" => confusion(data, args),
        "binary_metrics" => binary(data, args),
        "roc_curve" | "precision_recall_curve" | "lift_curve" => curve(op, data, args),
        "calibration_curve" => calibration(data, args),
        "profit_decision" => profit(data, args),
        "ensemble_predict" => ensemble(data, args),
        "model_compare" => compare(data, args),
        "holdout_split" => holdout(data, args),
        "rolling_origin" => rolling(data, args),
        _ => Err("unknown evaluation operation".into()),
    }
}
