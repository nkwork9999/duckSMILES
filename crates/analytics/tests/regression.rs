use ducksmiles_analytics::regression::run;

// These readers inspect the public numeric JSON payload without a test-only dependency.
fn numeric(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = json[start..].trim_start();
    let end = tail.find([',', '}', ']']).unwrap_or(tail.len());
    tail[..end].parse().expect("finite JSON number")
}

fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let end = json[start..].find(']').unwrap() + start;
    if json[start..end].trim().is_empty() {
        return Vec::new();
    }
    json[start..end]
        .split(',')
        .map(|n| n.trim().parse().unwrap())
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
}

#[test]
fn ridge_penalty_uses_standardized_features_and_returns_original_units() {
    // x mean=1, population variance=2/3. With lambda=3 the slope 2 shrinks to 1.
    let json = run("ridge", &[vec![0., 1.], vec![1., 3.], vec![2., 5.]], &[3.]).unwrap();
    close(numeric(&json, "version"), 1.);
    close(numeric(&json, "intercept"), 2.);
    close(array(&json, "coefficients")[0], 1.);
    close(numeric(&json, "rmse"), (2_f64 / 3.).sqrt());
    close(numeric(&json, "r2"), 0.75);
}

#[test]
fn ridge_centering_handles_feature_and_response_shifts() {
    let json = run(
        "ridge",
        &[vec![1000., 11.], vec![1001., 13.], vec![1002., 15.]],
        &[0.],
    )
    .unwrap();
    close(numeric(&json, "intercept"), -1989.);
    close(array(&json, "coefficients")[0], 2.);
}

#[test]
fn ridge_regularization_handles_constant_and_duplicate_columns() {
    let json = run(
        "ridge",
        &[
            vec![0., 0., 7., 1.],
            vec![1., 1., 7., 3.],
            vec![2., 2., 7., 5.],
        ],
        &[3.],
    )
    .unwrap();
    let coefficients = array(&json, "coefficients");
    close(coefficients[0], 2. / 3.);
    close(coefficients[1], 2. / 3.);
    close(coefficients[2], 0.);
    assert!(run(
        "ridge",
        &[vec![0., 0., 1.], vec![1., 1., 3.], vec![2., 2., 5.]],
        &[0.]
    )
    .is_err());
}

#[test]
fn pls_deflates_multiple_components_to_recover_a_linear_relation() {
    let data = vec![
        vec![0., 0., 3.],
        vec![1., 0., 5.],
        vec![0., 2., -5.],
        vec![2., 1., 3.],
        vec![3., 5., -11.],
    ];
    let json = run("pls", &data, &[2.]).unwrap();
    close(numeric(&json, "intercept"), 3.);
    let coefficients = array(&json, "coefficients");
    close(coefficients[0], 2.);
    close(coefficients[1], -4.);
    close(numeric(&json, "rmse"), 0.);
    close(numeric(&json, "components_used"), 2.);
}

#[test]
fn pls_stops_at_predictive_rank_and_supports_constant_targets() {
    let json = run(
        "pls",
        &[
            vec![0., 0., 1.],
            vec![1., 1., 4.],
            vec![2., 2., 7.],
            vec![3., 3., 10.],
        ],
        &[2.],
    )
    .unwrap();
    close(numeric(&json, "components_used"), 1.);
    close(numeric(&json, "rmse"), 0.);
    let constant = run("pls", &[vec![0., 9.], vec![1., 9.], vec![2., 9.]], &[1.]).unwrap();
    close(numeric(&constant, "intercept"), 9.);
    close(array(&constant, "coefficients")[0], 0.);
    assert!(constant.contains("\"r2\":null"));
}

#[test]
fn predict_consumes_intercept_then_original_unit_coefficients() {
    let json = run("predict", &[vec![1., 2.], vec![3., -1.]], &[3., 2., -4.]).unwrap();
    assert_eq!(array(&json, "predictions"), vec![-3., 13.]);
    assert!(run("predict", &[vec![1.]], &[3., 2., -4.]).is_err());
}

#[test]
fn cross_validation_fits_scaling_on_each_training_fold_only() {
    // lambda=2 halves a slope trained on two samples. Each fold has a different variance.
    let data = vec![
        vec![0., 0., 0.],
        vec![2., 2., 0.],
        vec![10., 10., 1.],
        vec![14., 14., 1.],
    ];
    let json = run("cross_validate", &data, &[0., 2.]).unwrap();
    let predictions = array(&json, "out_of_fold_predictions");
    for (actual, expected) in predictions.iter().zip([6., 7., 5.5, 7.5]) {
        close(*actual, expected);
    }
    close(numeric(&json, "rmse"), 30.875_f64.sqrt());
    close(numeric(&json, "r2"), 1. - 123.5 / 131.);
}

#[test]
fn cross_validation_runs_pls_and_rejects_bad_fold_assignments() {
    let data = vec![
        vec![0., 1., 0.],
        vec![1., 3., 1.],
        vec![2., 5., 0.],
        vec![3., 7., 1.],
    ];
    let json = run("cross_validate", &data, &[1., 1.]).unwrap();
    let predictions = array(&json, "out_of_fold_predictions");
    for (actual, expected) in predictions.iter().zip([1., 3., 5., 7.]) {
        close(*actual, expected);
    }
    for invalid in [
        vec![vec![1., 2., 0.], vec![2., 3., 0.]],
        vec![vec![1., 2., -1.], vec![2., 3., 0.]],
        vec![vec![1., 2., 0.5], vec![2., 3., 1.]],
    ] {
        assert!(run("cross_validate", &invalid, &[0., 1.]).is_err());
    }
}

#[test]
fn group_folds_keep_groups_indivisible_and_membership_permutation_invariant() {
    let data = vec![
        vec![10.],
        vec![10.],
        vec![20.],
        vec![30.],
        vec![40.],
        vec![40.],
    ];
    let json = run("group_folds", &data, &[3., 42.]).unwrap();
    let folds = array(&json, "folds");
    assert_eq!(folds[0], folds[1]);
    assert_eq!(folds[4], folds[5]);
    let mut counts = [0; 3];
    for f in &folds {
        counts[*f as usize] += 1;
    }
    assert_eq!(counts, [2, 2, 2]);
    let reordered = vec![
        vec![40.],
        vec![30.],
        vec![10.],
        vec![20.],
        vec![40.],
        vec![10.],
    ];
    let reordered_folds = array(
        &run("group_folds", &reordered, &[3., 42.]).unwrap(),
        "folds",
    );
    for (index, original) in [4, 3, 0, 2, 5, 1].iter().enumerate() {
        assert_eq!(reordered_folds[index], folds[*original]);
    }
    assert!(run("group_folds", &data, &[5., 42.]).is_err());
}

#[test]
fn rejects_nonfinite_ragged_missing_and_invalid_parameters() {
    for operation in ["ridge", "pls", "predict", "cross_validate", "group_folds"] {
        assert!(run(operation, &[], &[]).is_err());
        assert!(run(operation, &[vec![f64::NAN, 1.]], &[1.]).is_err());
        assert!(run(operation, &[vec![1., 2.], vec![3.]], &[1.]).is_err());
    }
    let data = [vec![0., 1.], vec![1., 3.]];
    assert!(run("ridge", &data, &[-1.]).is_err());
    assert!(run("ridge", &data, &[f64::INFINITY]).is_err());
    assert!(run("pls", &data, &[0.]).is_err());
    assert!(run("pls", &data, &[1.5]).is_err());
    assert!(run("unknown", &data, &[1.]).is_err());
    assert!(run("group_folds", &[vec![0.5], vec![1.]], &[2., 1.]).is_err());
}

#[test]
fn decimal_constant_features_do_not_create_spurious_predictive_rank() {
    // Repeated decimal values must stay exactly centered, despite binary representation.
    let data = vec![vec![0.1, 0.1]; 21];
    let json = run("pls", &data, &[1.]).unwrap();
    close(numeric(&json, "components_used"), 0.);
    assert!(json.contains("\"r2\":null"));
}

#[test]
fn numerical_overflow_returns_errors_instead_of_invalid_json() {
    assert!(run("predict", &[vec![f64::MAX]], &[0., 2.]).is_err());
    assert!(run(
        "ridge",
        &[vec![-f64::MAX, -f64::MAX], vec![f64::MAX, f64::MAX]],
        &[1.]
    )
    .is_err());
}

#[test]
fn regression_preserves_small_spreads_at_large_feature_and_target_origins() {
    let origin = 1e16;
    for operation in ["ridge", "pls"] {
        let lambda_or_components = if operation == "ridge" { 0. } else { 1. };
        let data = vec![
            vec![origin, 1.],
            vec![origin + 2., 3.],
            vec![origin + 4., 5.],
            vec![origin + 6., 7.],
        ];
        let json = run(operation, &data, &[lambda_or_components]).unwrap();
        close(array(&json, "coefficients")[0], 1.);
        close(numeric(&json, "rmse"), 0.);
        let portable = array(&json, "stable_predict_parameters");
        let predictions = array(
            &run(
                "predict",
                &[
                    vec![origin],
                    vec![origin + 2.],
                    vec![origin + 4.],
                    vec![origin + 6.],
                ],
                &portable,
            )
            .unwrap(),
            "predictions",
        );
        for (actual, expected) in predictions.iter().zip([1., 3., 5., 7.]) {
            close(*actual, expected);
        }
        let shifted_target = vec![
            vec![0., origin],
            vec![2., origin + 2.],
            vec![4., origin + 4.],
            vec![6., origin + 6.],
        ];
        let json = run(operation, &shifted_target, &[lambda_or_components]).unwrap();
        close(array(&json, "coefficients")[0], 1.);
        close(numeric(&json, "rmse"), 0.);
    }
}
