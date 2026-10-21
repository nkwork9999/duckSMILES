use ducksmiles_analytics::execute as run;

fn numbers(json: &str, key: &str) -> Vec<Option<f64>> {
    let tail = json.split_once(&format!("\"{key}\":")).expect(key).1;
    let end = if tail.starts_with('[') {
        let mut depth = 0;
        tail.char_indices()
            .find_map(|(i, c)| {
                if c == '[' {
                    depth += 1;
                }
                if c == ']' {
                    depth -= 1;
                    if depth == 0 {
                        return Some(i + 1);
                    }
                }
                None
            })
            .unwrap()
    } else {
        tail.find([',', '}']).unwrap_or(tail.len())
    };
    tail[..end]
        .split(['[', ']', ',', ' '])
        .filter(|v| !v.is_empty())
        .map(|v| {
            if v == "null" {
                None
            } else {
                Some(v.parse().unwrap())
            }
        })
        .collect()
}
fn scalar(j: &str, k: &str) -> f64 {
    numbers(j, k)[0].unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() <= 2e-10 * b.abs().max(1.), "{a} != {b}");
}
fn vector(j: &str, k: &str, want: &[f64]) {
    let got = numbers(j, k);
    assert_eq!(got.len(), want.len());
    for (a, &b) in got.iter().zip(want) {
        near(a.unwrap(), b);
    }
}
fn scored() -> Vec<Vec<f64>> {
    vec![vec![0., 0.1], vec![1., 0.4], vec![0., 0.4], vec![1., 0.8]]
}

#[test]
fn weighted_regression_metrics_preserve_residual_and_target_units() {
    let j = run(
        "regression_metrics",
        &[vec![1., 2., 1.], vec![2., 2., 2.], vec![3., 1., 1.]],
        &[],
    )
    .unwrap();
    near(scalar(&j, "rmse"), (1.25_f64).sqrt());
    near(scalar(&j, "mae"), 0.75);
    near(scalar(&j, "r_squared"), -1.5);
    near(scalar(&j, "bias"), 0.25);
    vector(&j, "residuals", &[-1., 0., 2.]);
    let tiny = run(
        "regression_metrics",
        &[vec![1e-200, 2e-200], vec![3e-200, 2e-200]],
        &[],
    )
    .unwrap();
    near(scalar(&tiny, "rmse") / 1e-200, 1.);
    near(scalar(&tiny, "r_squared"), 0.);
    let shifted = run(
        "regression_metrics",
        &[
            vec![1e16, 1e16],
            vec![1e16 + 2., 1e16],
            vec![1e16 + 4., 1e16 + 2.],
        ],
        &[],
    )
    .unwrap();
    near(scalar(&shifted, "r_squared"), 0.);
}

#[test]
fn confusion_matrix_counts_multiclass_actual_rows_and_predicted_columns() {
    let j = run(
        "confusion_matrix",
        &[
            vec![10., 10.],
            vec![10., 20.],
            vec![20., 20.],
            vec![30., 20.],
        ],
        &[],
    )
    .unwrap();
    vector(&j, "labels", &[10., 20., 30.]);
    vector(&j, "matrix", &[1., 1., 0., 0., 1., 0., 0., 1., 0.]);
    near(scalar(&j, "accuracy"), 0.5);
    near(scalar(&j, "balanced_accuracy"), 0.5);
    near(scalar(&j, "macro_f1"), 7. / 18.);
    assert_eq!(numbers(&j, "precision")[2], None);
}

#[test]
fn binary_metrics_use_probability_loss_and_inclusive_threshold() {
    let j = run("binary_metrics", &scored(), &[0.4]).unwrap();
    vector(&j, "counts", &[1., 1., 0., 2.]); // [TN,FP,FN,TP]
    near(scalar(&j, "accuracy"), 0.75);
    near(scalar(&j, "brier_score"), 0.1425);
    near(scalar(&j, "log_loss"), 0.4389051056530454);
}

#[test]
fn roc_groups_equal_scores_before_integrating_and_handles_weights() {
    let j = run("roc_curve", &scored(), &[]).unwrap();
    vector(&j, "false_positive_rate", &[0., 0., 0.5, 1.]);
    vector(&j, "true_positive_rate", &[0., 0.5, 1., 1.]);
    near(scalar(&j, "auc"), 0.875);
    assert_eq!(numbers(&j, "thresholds")[0], None);
    let weighted = run(
        "roc_curve",
        &[vec![0., 0.4, 2.], vec![1., 0.4, 1.], vec![1., 0.8, 3.]],
        &[],
    )
    .unwrap();
    near(scalar(&weighted, "auc"), 0.875);
}

#[test]
fn precision_recall_reports_step_average_precision_separately_from_trapezoids() {
    let j = run("precision_recall_curve", &scored(), &[]).unwrap();
    vector(&j, "precision", &[1., 1., 2. / 3., 0.5]);
    vector(&j, "recall", &[0., 0.5, 1., 1.]);
    near(scalar(&j, "average_precision"), 5. / 6.);
    near(scalar(&j, "trapezoid_area"), 11. / 12.);
}

#[test]
fn lift_uses_cumulative_positive_rate_over_population_positive_rate() {
    let j = run("lift_curve", &scored(), &[]).unwrap();
    vector(&j, "selection_fraction", &[0., 0.25, 0.75, 1.]);
    assert_eq!(numbers(&j, "lift")[0], None);
    let got = numbers(&j, "lift");
    near(got[1].unwrap(), 2.);
    near(got[2].unwrap(), 4. / 3.);
    near(got[3].unwrap(), 1.);
}

#[test]
fn calibration_includes_probability_one_and_keeps_empty_bins_explicit() {
    let j = run(
        "calibration_curve",
        &[vec![0., 0.], vec![1., 0.5], vec![1., 1.]],
        &[4.],
    )
    .unwrap();
    vector(&j, "bin_weight", &[1., 0., 1., 1.]);
    assert_eq!(
        numbers(&j, "mean_probability"),
        vec![Some(0.), None, Some(0.5), Some(1.)]
    );
    near(scalar(&j, "expected_calibration_error"), 1. / 6.);
}

#[test]
fn profit_decision_maximizes_expected_profit_with_actual_by_decision_matrix() {
    let j = run(
        "profit_decision",
        &[vec![0., 0.1], vec![1., 0.6], vec![1., 0.2]],
        &[2., -1., -5., 10.],
    )
    .unwrap();
    vector(&j, "decisions", &[0., 1., 1.]);
    vector(&j, "actual_profit", &[2., 10., 10.]);
    near(scalar(&j, "total_actual_profit"), 22.);
    near(scalar(&j, "mean_expected_profit"), 2.7);
}

#[test]
fn model_comparison_ranks_supplied_predictions_without_refitting() {
    let j = run(
        "model_compare",
        &[vec![0., 0., 1.], vec![1., 1., 0.4], vec![2., 2., 1.]],
        &[],
    )
    .unwrap();
    vector(&j, "ranked_model_indices", &[0., 1.]);
    near(numbers(&j, "scores")[0].unwrap(), 0.);
    near(numbers(&j, "scores")[1].unwrap(), (2.36_f64 / 3.).sqrt());
}

#[test]
fn ensemble_normalizes_member_weights_without_losing_large_origin_spreads() {
    let j = run("ensemble_predict", &[vec![1., 3.], vec![3., 7.]], &[1., 3.]).unwrap();
    vector(&j, "predictions", &[2.5, 6.]);
    let j = run(
        "ensemble_predict",
        &[vec![1e16, 1e16 + 4.]],
        &[1e300, 1e300],
    )
    .unwrap();
    assert_eq!(scalar(&j, "predictions"), 1e16 + 2.);
}

#[test]
fn stratified_holdout_is_seeded_lossless_and_retains_both_levels() {
    let data = vec![
        vec![0.],
        vec![0.],
        vec![0.],
        vec![0.],
        vec![1.],
        vec![1.],
        vec![1.],
        vec![1.],
    ];
    let j = run("holdout_split", &data, &[0.75, 42., 1.]).unwrap();
    assert_eq!(j, run("holdout_split", &data, &[0.75, 42., 1.]).unwrap());
    let train = numbers(&j, "train_indices");
    let test = numbers(&j, "test_indices");
    assert_eq!(train.len(), 6);
    assert_eq!(test.len(), 2);
    let mut all: Vec<_> = train
        .iter()
        .chain(&test)
        .map(|v| v.unwrap() as usize)
        .collect();
    all.sort_unstable();
    assert_eq!(all, (0..8).collect::<Vec<_>>());
    assert_ne!(
        data[test[0].unwrap() as usize],
        data[test[1].unwrap() as usize]
    );
}

#[test]
fn rolling_origin_never_trains_on_current_or_future_validation_rows() {
    let j = run(
        "rolling_origin",
        &(0..8).map(|i| vec![i as f64]).collect::<Vec<_>>(),
        &[3., 2., 2., 3.],
    )
    .unwrap();
    vector(&j, "train_indices", &[0., 1., 2., 2., 3., 4.]);
    vector(&j, "test_indices", &[3., 4., 5., 6.]);
    vector(&j, "origins", &[3., 5.]);
}

#[test]
fn evaluation_rejects_undefined_auc_invalid_probabilities_and_excessive_fold_output() {
    assert!(run("roc_curve", &[vec![1., 0.5], vec![1., 0.8]], &[]).is_err());
    assert!(run("binary_metrics", &[vec![0., 2.]], &[]).is_err());
    assert!(run("confusion_matrix", &[vec![0.1, 1.]], &[]).is_err());
    assert!(run("regression_metrics", &[vec![1., 1., 0.]], &[]).is_err());
    assert!(run("holdout_split", &[vec![0.], vec![1.]], &[0.5, 0., 1.]).is_err());
    assert!(run("ensemble_predict", &[vec![1., 2.]], &[-1., 1.]).is_err());
    assert!(run("rolling_origin", &vec![vec![1.]; 2000], &[1., 1.]).is_err());
}

#[test]
fn weighted_curves_do_not_drop_the_minority_class_through_total_subtraction() {
    let j = run("roc_curve", &[vec![1., 0.8, 1e300], vec![0., 0.2, 1.]], &[]).unwrap();
    near(scalar(&j, "auc"), 1.);
}

#[test]
fn f1_remains_finite_when_doubling_valid_count_weights_would_overflow() {
    let j = run("confusion_matrix", &[vec![1., 1., 1e308]], &[]).unwrap();
    near(scalar(&j, "macro_f1"), 1.);
    let j = run("binary_metrics", &[vec![1., 1., 1e308]], &[]).unwrap();
    near(scalar(&j, "f1"), 1.);
}
