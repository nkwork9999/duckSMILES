use ducksmiles_analytics::mixed::run;

fn numeric(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = &json[start..];
    let end = tail.find([',', '}', ']']).unwrap_or(tail.len());
    tail[..end].parse().unwrap()
}
fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let end = json[start..].find(']').unwrap() + start;
    if json[start..end].is_empty() {
        return vec![];
    }
    json[start..end]
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect()
}
fn close(a: f64, e: f64, t: f64) {
    assert!((a - e).abs() <= t * (1. + e.abs()), "{a} != {e}");
}
fn balanced() -> Vec<Vec<f64>> {
    vec![
        vec![1., 0.],
        vec![2., 0.],
        vec![1., 0.],
        vec![2., 0.],
        vec![4., 1.],
        vec![5., 1.],
        vec![4., 1.],
        vec![5., 1.],
        vec![7., 2.],
        vec![8., 2.],
        vec![7., 2.],
        vec![8., 2.],
    ]
}

#[test]
fn balanced_reml_recovers_anova_variance_components_and_blups() {
    // Mutation caught: using ML residual df or omitting fixed-design REML determinant.
    let json = run("mixed_model", &balanced(), &[]).unwrap();
    assert!(json.contains("\"converged\":true"));
    close(numeric(&json, "intercept"), 4.5, 1e-9);
    close(numeric(&json, "residual_variance"), 1. / 3., 2e-6);
    close(
        numeric(&json, "random_intercept_variance"),
        107. / 12.,
        2e-6,
    );
    close(
        numeric(&json, "intercept_standard_error"),
        3_f64.sqrt(),
        2e-6,
    );
    close(numeric(&json, "log_likelihood"), -15.490540829595014, 1e-8);
    let random = array(&json, "random_intercepts");
    close(random[0], -107. / 36., 2e-6);
    close(random[1], 0., 1e-8);
    close(random[2], 107. / 36., 2e-6);
    assert!(json.contains("\"aic\":null"));
}

#[test]
fn balanced_ml_uses_full_likelihood_and_aic_parameter_count() {
    let json = run("mixed_model", &balanced(), &[0.]).unwrap();
    close(numeric(&json, "residual_variance"), 1. / 3., 2e-6);
    close(numeric(&json, "random_intercept_variance"), 71. / 12., 2e-6);
    close(
        numeric(&json, "intercept_standard_error"),
        2_f64.sqrt(),
        2e-6,
    );
    close(numeric(&json, "log_likelihood"), -16.850587844971496, 1e-8);
    close(numeric(&json, "aic"), 39.70117568994299, 1e-8);
}

#[test]
fn unbalanced_fixed_covariate_fit_matches_independent_dense_scipy_likelihood() {
    let rows = [
        vec![-2., -1.4, 0.],
        vec![0., 1.1, 0.],
        vec![1., 2.3, 0.],
        vec![0., 3.3, 1.],
        vec![2., 5.8, 1.],
        vec![-1., 0.2, 2.],
        vec![0., 1.9, 2.],
        vec![1., 3.1, 2.],
        vec![3., 5.4, 2.],
        vec![-2., -2.1, 3.],
        vec![2., 3.2, 3.],
    ];
    for (method, se2, sb2, intercept, slope, ll) in [
        (
            1.,
            0.024989341562286042,
            1.3639029341277962,
            1.6492975136512928,
            1.283270756743652,
            -6.3848188101185706,
        ),
        (
            0.,
            0.021430315590386825,
            1.020866382458063,
            1.6491705050937222,
            1.2834755425573339,
            -4.151452235322462,
        ),
    ] {
        let fit = run("mixed_model", &rows, &[method, 200., 1e-10]).unwrap();
        close(numeric(&fit, "residual_variance"), se2, 3e-6);
        close(numeric(&fit, "random_intercept_variance"), sb2, 3e-6);
        close(numeric(&fit, "intercept"), intercept, 3e-6);
        close(array(&fit, "coefficients")[0], slope, 3e-6);
        close(numeric(&fit, "log_likelihood"), ll, 1e-8);
    }
}

#[test]
fn zero_random_variance_boundary_is_considered_explicitly() {
    let rows = [
        vec![-1., 0.],
        vec![0., 0.],
        vec![1., 0.],
        vec![-1., 1.],
        vec![0., 1.],
        vec![1., 1.],
        vec![-1., 2.],
        vec![0., 2.],
        vec![1., 2.],
    ];
    let json = run("mixed_model", &rows, &[1.]).unwrap();
    close(numeric(&json, "random_intercept_variance"), 0., 0.);
    close(numeric(&json, "residual_variance"), 0.75, 1e-10);
    assert!(json.contains("\"variance_boundary\":true"));
}

#[test]
fn mixed_prediction_distinguishes_known_and_new_groups() {
    let fit = run("mixed_model", &balanced(), &[]).unwrap();
    let model = array(&fit, "predict_parameters");
    let json = run(
        "mixed_predict",
        &[vec![0.], vec![1.], vec![2.], vec![99.]],
        &model,
    )
    .unwrap();
    let marginal = array(&json, "marginal_predictions");
    for y in marginal {
        close(y, 4.5, 1e-9);
    }
    let conditional = array(&json, "conditional_predictions");
    close(conditional[0], 4.5 - 107. / 36., 2e-6);
    close(conditional[2], 4.5 + 107. / 36., 2e-6);
    close(conditional[3], 4.5, 1e-9);
}

#[test]
fn finite_iteration_budget_does_not_report_false_convergence() {
    let fit = run("mixed_model", &balanced(), &[1., 1., 1e-12]).unwrap();
    assert!(fit.contains("\"converged\":false"));
}

#[test]
fn rejects_invalid_groups_rank_degeneracy_nonfinite_and_bad_arguments() {
    for data in [
        vec![],
        vec![vec![f64::NAN, 0.]],
        vec![vec![1., 0.], vec![2.]],
        vec![vec![1., 0.], vec![2., 0.], vec![3., 0.]],
        vec![vec![1., 0.], vec![2., 1.], vec![3., 2.]],
        vec![vec![1., 0.5], vec![2., 1.], vec![3., 1.]],
        vec![
            vec![2., 1., 0.],
            vec![2., 2., 0.],
            vec![2., 3., 1.],
            vec![2., 4., 1.],
        ],
    ] {
        assert!(run("mixed_model", &data, &[]).is_err());
    }
    for args in [
        vec![2.],
        vec![1., 0.],
        vec![1., 1.5],
        vec![1., 2., 0.],
        vec![1., 2., f64::NAN],
    ] {
        assert!(run("mixed_model", &balanced(), &args).is_err());
    }
    assert!(run("mixed_predict", &[vec![0.]], &[]).is_err());
    assert!(run(
        "mixed_model",
        &[vec![1., 0.], vec![1., 0.], vec![1., 1.], vec![1., 1.]],
        &[]
    )
    .is_err());
}

#[test]
fn exactly_determined_fixed_response_does_not_claim_a_finite_variance_optimum() {
    let rows = [
        vec![0., 1., 0.],
        vec![1., 3., 0.],
        vec![2., 5., 1.],
        vec![3., 7., 1.],
    ];
    assert!(run("mixed_model", &rows, &[]).is_err());
}

#[test]
fn reml_likelihood_accounts_for_original_fixed_feature_units() {
    let rows = [
        vec![-2., -1.4, 0.],
        vec![0., 1.1, 0.],
        vec![1., 2.3, 0.],
        vec![0., 3.3, 1.],
        vec![2., 5.8, 1.],
        vec![-1., 0.2, 2.],
        vec![0., 1.9, 2.],
        vec![1., 3.1, 2.],
        vec![3., 5.4, 2.],
        vec![-2., -2.1, 3.],
        vec![2., 3.2, 3.],
    ];
    let shifted = rows
        .iter()
        .map(|r| vec![r[0] * 10. + 100., r[1] + 50., r[2]])
        .collect::<Vec<_>>();
    for method in [0., 1.] {
        let base = run("mixed_model", &rows, &[method]).unwrap();
        let transformed = run("mixed_model", &shifted, &[method]).unwrap();
        close(
            numeric(&transformed, "residual_variance"),
            numeric(&base, "residual_variance"),
            1e-6,
        );
        close(
            numeric(&transformed, "random_intercept_variance"),
            numeric(&base, "random_intercept_variance"),
            1e-6,
        );
        close(
            array(&transformed, "coefficients")[0] * 10.,
            array(&base, "coefficients")[0],
            1e-7,
        );
        let ll_shift = if method == 1. { -10_f64.ln() } else { 0. };
        close(
            numeric(&transformed, "log_likelihood"),
            numeric(&base, "log_likelihood") + ll_shift,
            1e-9,
        );
    }
}

#[test]
fn mixed_prediction_keeps_large_absolute_feature_origins_separate() {
    let rows = [
        vec![0., -1.4, 0.],
        vec![2., 1.1, 0.],
        vec![4., 2.3, 0.],
        vec![2., 3.3, 1.],
        vec![6., 5.8, 1.],
        vec![0., 0.2, 2.],
        vec![2., 1.9, 2.],
        vec![4., 3.1, 2.],
        vec![8., 5.4, 2.],
        vec![0., -2.1, 3.],
        vec![6., 3.2, 3.],
    ];
    let shifted = rows
        .iter()
        .map(|r| vec![r[0] + 1e16, r[1], r[2]])
        .collect::<Vec<_>>();
    let base = run("mixed_model", &rows, &[]).unwrap();
    let fit = run("mixed_model", &shifted, &[]).unwrap();
    let queries = shifted.iter().map(|r| vec![r[0], r[2]]).collect::<Vec<_>>();
    let predictions = run(
        "mixed_predict",
        &queries,
        &array(&fit, "predict_parameters"),
    )
    .unwrap();
    for (a, e) in array(&predictions, "conditional_predictions")
        .iter()
        .zip(array(&base, "conditional_fitted_values"))
    {
        close(*a, e, 1e-7);
    }
}
