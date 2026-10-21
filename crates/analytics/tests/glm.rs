use ducksmiles_analytics::glm::run;

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

fn matrix(json: &str, key: &str) -> Vec<Vec<f64>> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &json[start..];
    let end = tail.find("]]").unwrap() + 1;
    tail[..end]
        .split("],")
        .map(|s| {
            s.trim_matches(['[', ']'])
                .split(',')
                .map(|x| x.parse().unwrap())
                .collect()
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1. + expected.abs()),
        "{actual} != {expected}"
    );
}

#[test]
fn ordinary_least_squares_inference_matches_independent_scipy_fixture() {
    // Mutation caught: ridge penalty, normal approximation, or wrong residual df in OLS.
    let data = [
        vec![0., 1.],
        vec![1., 2.],
        vec![2., 1.],
        vec![3., 4.],
        vec![4., 4.],
    ];
    let json = run("linear_model", &data, &[]).unwrap();
    close(numeric(&json, "intercept"), 0.8, 1e-11);
    close(array(&json, "coefficients")[0], 0.8, 1e-11);
    close(numeric(&json, "sse"), 2.8, 1e-11);
    close(numeric(&json, "residual_df"), 3., 0.);
    close(numeric(&json, "sigma_squared"), 0.9333333333333333, 1e-11);
    close(numeric(&json, "r2"), 0.6956521739130435, 1e-11);
    close(numeric(&json, "adjusted_r2"), 0.5942028985507246, 1e-11);
    close(
        numeric(&json, "intercept_standard_error"),
        0.7483314773547882,
        1e-10,
    );
    close(
        array(&json, "coefficient_standard_errors")[0],
        0.30550504633038933,
        1e-10,
    );
    close(
        array(&json, "coefficient_t_statistics")[0],
        2.618614682831908,
        1e-10,
    );
    close(
        array(&json, "coefficient_p_values")[0],
        0.07909569938205474,
        1e-9,
    );
    close(numeric(&json, "overall_f"), 6.857142857142858, 1e-10);
    close(numeric(&json, "overall_p_value"), 0.07909569938205474, 1e-9);
    let ci = matrix(&json, "coefficient_confidence_intervals");
    close(ci[0][0], -0.17225341398776328, 1e-8);
    close(ci[0][1], 1.7722534139877633, 1e-8);
    for (a, e) in array(&json, "fitted_values")
        .iter()
        .zip([0.8, 1.6, 2.4, 3.2, 4.])
    {
        close(*a, e, 1e-11);
    }
}

#[test]
fn linear_inference_rejects_rank_deficiency_and_insufficient_residual_df() {
    assert!(run(
        "linear_model",
        &[
            vec![0., 0., 1.],
            vec![1., 1., 3.],
            vec![2., 2., 4.],
            vec![3., 3., 9.]
        ],
        &[]
    )
    .is_err());
    assert!(run("linear_model", &[vec![0., 1.], vec![1., 3.]], &[]).is_err());
    assert!(run(
        "linear_model",
        &[vec![1., 1.], vec![1., 2.], vec![1., 3.]],
        &[]
    )
    .is_err());
}

#[test]
fn stable_model_predictions_preserve_large_offsets_and_feature_units() {
    // Mutation caught: rounding an absolute mean before subtracting feature origins.
    let origin = 1e16;
    let data = [
        vec![origin, 1.],
        vec![origin + 2., 3.],
        vec![origin + 4., 5.],
        vec![origin + 6., 7.],
    ];
    for operation in ["linear_model", "elastic_net"] {
        let args = if operation == "elastic_net" {
            vec![0., 0.5]
        } else {
            vec![]
        };
        let fit = run(operation, &data, &args).unwrap();
        close(array(&fit, "coefficients")[0], 1., 1e-10);
        let model = array(&fit, "stable_predict_parameters");
        let prediction = run(
            "glm_predict",
            &[
                vec![origin],
                vec![origin + 2.],
                vec![origin + 4.],
                vec![origin + 6.],
            ],
            &model,
        )
        .unwrap();
        for (a, e) in array(&prediction, "predictions")
            .iter()
            .zip([1., 3., 5., 7.])
        {
            close(*a, e, 1e-10);
        }
    }
}

#[test]
fn logistic_likelihood_returns_exact_two_group_probabilities_and_deviance() {
    let data = [
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 1.],
        vec![1., 0.],
        vec![1., 1.],
        vec![1., 1.],
        vec![1., 1.],
    ];
    let json = run("logistic", &data, &[]).unwrap();
    assert!(json.contains("\"converged\":true"));
    close(numeric(&json, "intercept"), -3_f64.ln(), 1e-8);
    close(array(&json, "coefficients")[0], 2. * 3_f64.ln(), 1e-8);
    close(numeric(&json, "deviance"), 8.997362313900933, 1e-9);
    for (p, want) in array(&json, "fitted_values")
        .iter()
        .zip([0.25, 0.25, 0.25, 0.25, 0.75, 0.75, 0.75, 0.75])
    {
        close(*p, want, 1e-8);
    }
    let portable = array(&json, "stable_predict_parameters");
    let prediction = run("glm_predict", &[vec![0.], vec![1.]], &portable).unwrap();
    close(array(&prediction, "predictions")[0], 0.25, 1e-8);
    close(array(&prediction, "predictions")[1], 0.75, 1e-8);
}

#[test]
fn logistic_multivariable_solution_matches_independent_scipy_optimizer() {
    let data = [
        vec![-2., 1., 0.],
        vec![-1., 0., 0.],
        vec![0., 1., 1.],
        vec![1., -1., 0.],
        vec![2., 2., 1.],
        vec![3., 0., 1.],
        vec![4., 1., 1.],
        vec![1., 3., 0.],
    ];
    let json = run("logistic", &data, &[0., 200., 1e-10]).unwrap();
    close(numeric(&json, "intercept"), -1.2623082993160446, 2e-6);
    close(array(&json, "coefficients")[0], 1.1586255710775668, 2e-6);
    close(array(&json, "coefficients")[1], 0.1023925502064955, 2e-6);
}

#[test]
fn logistic_separation_and_exhausted_iterations_are_not_false_successes() {
    let data = [vec![-2., 0.], vec![-1., 0.], vec![1., 1.], vec![2., 1.]];
    let result = run("logistic", &data, &[]);
    assert!(result.is_err() || result.unwrap().contains("\"converged\":false"));
    let regularized = run("logistic", &data, &[1.]).unwrap();
    assert!(regularized.contains("\"converged\":true"));
    let early = run("logistic", &data, &[1., 1., 1e-14]).unwrap();
    assert!(early.contains("\"converged\":false"));
    assert!(run("logistic", &[vec![0., 1.], vec![1., 1.]], &[]).is_err());
}

#[test]
fn poisson_log_link_matches_group_means_and_exact_deviance() {
    let data = [
        vec![0., 0.],
        vec![0., 1.],
        vec![0., 2.],
        vec![0., 1.],
        vec![1., 2.],
        vec![1., 4.],
        vec![1., 6.],
        vec![1., 4.],
    ];
    let json = run("poisson", &data, &[]).unwrap();
    assert!(json.contains("\"converged\":true"));
    close(numeric(&json, "intercept"), 0., 1e-8);
    close(array(&json, "coefficients")[0], 4_f64.ln(), 1e-8);
    close(numeric(&json, "deviance"), 4.865581297297973, 1e-8);
    for (a, e) in array(&json, "fitted_values")
        .iter()
        .zip([1., 1., 1., 1., 4., 4., 4., 4.])
    {
        close(*a, e, 1e-8);
    }
    assert!(run("poisson", &[vec![0., 0.5], vec![1., 1.]], &[]).is_err());
    assert!(run("poisson", &[vec![0., -1.], vec![1., 1.]], &[]).is_err());
    assert!(run("poisson", &[vec![0., 0.], vec![1., 0.]], &[]).is_err());
}

#[test]
fn elastic_net_soft_thresholds_standardized_features_and_unpenalized_intercept() {
    let data = [
        vec![-1., -1., 0.8],
        vec![-1., 1., 1.2],
        vec![1., -1., 4.8],
        vec![1., 1., 5.2],
    ];
    let fit = run("elastic_net", &data, &[0.5, 1.]).unwrap();
    assert!(fit.contains("\"converged\":true"));
    close(numeric(&fit, "intercept"), 3., 1e-10);
    close(array(&fit, "coefficients")[0], 1.5, 1e-10);
    close(array(&fit, "coefficients")[1], 0., 1e-10);
    close(numeric(&fit, "objective"), 0.895, 1e-10);
    let ridge = run("elastic_net", &data, &[0.5, 0.]).unwrap();
    close(array(&ridge, "coefficients")[0], 4. / 3., 1e-10);
    close(array(&ridge, "coefficients")[1], 2. / 15., 1e-10);
}

#[test]
fn elastic_net_regularizes_duplicate_and_constant_features_without_false_convergence() {
    let data = [
        vec![0., 0., 7., 1.],
        vec![1., 1., 7., 3.],
        vec![2., 2., 7., 5.],
    ];
    let fit = run("elastic_net", &data, &[0.2, 0.5]).unwrap();
    close(array(&fit, "coefficients")[2], 0., 1e-10);
    assert!(fit.contains("\"converged\":true"));
    let early = run("elastic_net", &data, &[0.2, 0.5, 1., 1e-14]).unwrap();
    assert!(early.contains("\"converged\":false"));
}

#[test]
fn polynomial_features_include_all_total_degree_terms_in_stable_order() {
    let json = run(
        "polynomial_features",
        &[vec![2., 3.], vec![-1., 2.]],
        &[3., 1.],
    )
    .unwrap();
    assert_eq!(
        matrix(&json, "powers"),
        vec![
            vec![0., 0.],
            vec![1., 0.],
            vec![0., 1.],
            vec![2., 0.],
            vec![1., 1.],
            vec![0., 2.],
            vec![3., 0.],
            vec![2., 1.],
            vec![1., 2.],
            vec![0., 3.]
        ]
    );
    assert_eq!(
        matrix(&json, "transformed_rows"),
        vec![
            vec![1., 2., 3., 4., 6., 9., 8., 12., 18., 27.],
            vec![1., -1., 2., 1., -2., 4., -1., 2., -4., 8.]
        ]
    );
    let no_bias = run("polynomial_features", &[vec![2., 3.]], &[1., 0.]).unwrap();
    assert_eq!(matrix(&no_bias, "transformed_rows"), vec![vec![2., 3.]]);
}

#[test]
fn prediction_link_functions_are_stable_and_reject_poisson_overflow() {
    assert_eq!(
        array(
            &run("glm_predict", &[vec![1.], vec![2.]], &[0., 3., 2.]).unwrap(),
            "predictions"
        ),
        vec![5., 7.]
    );
    assert_eq!(
        array(
            &run("glm_predict", &[vec![-1000.], vec![1000.]], &[1., 0., 1.]).unwrap(),
            "predictions"
        ),
        vec![0., 1.]
    );
    let log = run("glm_predict", &[vec![0.], vec![1.]], &[2., 0., 2_f64.ln()]).unwrap();
    close(array(&log, "predictions")[1], 2., 1e-12);
    assert!(run("glm_predict", &[vec![1000.]], &[2., 0., 1.]).is_err());
    assert!(run("glm_predict", &[vec![1.]], &[3., 0., 1.]).is_err());
}

#[test]
fn invalid_matrices_parameters_and_bounded_work_return_errors() {
    for operation in [
        "linear_model",
        "logistic",
        "poisson",
        "elastic_net",
        "polynomial_features",
        "glm_predict",
    ] {
        assert!(run(operation, &[], &[]).is_err());
        assert!(run(operation, &[vec![1., 2.], vec![3.]], &[]).is_err());
        assert!(run(operation, &[vec![f64::NAN, 2.]], &[]).is_err());
        assert!(run(operation, &[vec![1., 2.]], &[f64::INFINITY]).is_err());
    }
    let data = [vec![0., 0.], vec![1., 1.], vec![2., 1.]];
    for args in [
        vec![-1.],
        vec![0., 0.],
        vec![0., 1.5],
        vec![0., 2., 0.],
        vec![0., 2., -1.],
    ] {
        assert!(run("logistic", &data, &args).is_err());
    }
    for args in [vec![-1., 0.5], vec![1., -0.1], vec![1., 1.1]] {
        assert!(run("elastic_net", &data, &args).is_err());
    }
    for args in [vec![0.], vec![4.], vec![2., 0.5]] {
        assert!(run("polynomial_features", &data, &args).is_err());
    }
    assert!(run("polynomial_features", &[vec![f64::MAX]], &[2.]).is_err());
    assert!(run("polynomial_features", &[vec![1.; 256]], &[3.]).is_err());
    assert!(run("linear_model", &data, &[1.]).is_err());
    assert!(run("linear_model", &data, &[0.]).is_err());
    assert!(run("unknown", &data, &[]).is_err());
}

#[test]
fn extremely_high_representable_confidence_retains_nonzero_tail_probability() {
    // Mutation caught: (1+confidence)/2 rounds to 1 before obtaining the t tail.
    let data = [
        vec![0., 1.],
        vec![1., 2.],
        vec![2., 1.],
        vec![3., 4.],
        vec![4., 4.],
    ];
    let json = run("linear_model", &data, &[1. - f64::EPSILON / 2.]).unwrap();
    let interval = matrix(&json, "coefficient_confidence_intervals");
    close(interval[0][0], -82737.2397048031, 1e-8);
    close(interval[0][1], 82738.8397048031, 1e-8);
}

#[test]
fn generalized_models_retain_predictions_after_large_feature_translation() {
    let origin = 1e16;
    let logistic = [
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 1.],
        vec![2., 0.],
        vec![2., 1.],
        vec![2., 1.],
        vec![2., 1.],
    ];
    let poisson = [
        vec![0., 0.],
        vec![0., 1.],
        vec![0., 2.],
        vec![0., 1.],
        vec![2., 2.],
        vec![2., 4.],
        vec![2., 6.],
        vec![2., 4.],
    ];
    for (operation, data) in [("logistic", logistic), ("poisson", poisson)] {
        let base = run(operation, &data, &[0.3]).unwrap();
        let shifted = data
            .iter()
            .map(|r| vec![r[0] + origin, r[1]])
            .collect::<Vec<_>>();
        let fit = run(operation, &shifted, &[0.3]).unwrap();
        for (a, e) in array(&fit, "fitted_values")
            .iter()
            .zip(array(&base, "fitted_values"))
        {
            close(*a, e, 1e-9);
        }
        let predictions = run(
            "glm_predict",
            &shifted.iter().map(|r| vec![r[0]]).collect::<Vec<_>>(),
            &array(&fit, "stable_predict_parameters"),
        )
        .unwrap();
        for (a, e) in array(&predictions, "predictions")
            .iter()
            .zip(array(&base, "fitted_values"))
        {
            close(*a, e, 1e-9);
        }
    }
}

#[test]
fn standardized_penalties_are_invariant_to_feature_units() {
    let data = [
        vec![-1., -1., 0.8],
        vec![-1., 1., 1.2],
        vec![1., -1., 4.8],
        vec![1., 1., 5.2],
    ];
    let changed = data
        .iter()
        .map(|r| vec![r[0] * 1e-6 + 100., r[1] * 1e3 - 40., r[2]])
        .collect::<Vec<_>>();
    let base = run("elastic_net", &data, &[0.5, 0.3]).unwrap();
    let fit = run("elastic_net", &changed, &[0.5, 0.3]).unwrap();
    for (a, e) in array(&fit, "fitted_values")
        .iter()
        .zip(array(&base, "fitted_values"))
    {
        close(*a, e, 2e-8);
    }
    close(
        array(&fit, "coefficients")[0] * 1e-6,
        array(&base, "coefficients")[0],
        2e-8,
    );
    close(
        array(&fit, "coefficients")[1] * 1e3,
        array(&base, "coefficients")[1],
        2e-8,
    );
}

#[test]
fn multivariable_linear_inference_matches_independent_matrix_fixture() {
    let data = [
        vec![-2., 1., 1.],
        vec![-1., 0., 0.],
        vec![0., 1., 2.],
        vec![1., -1., 3.],
        vec![2., 2., 1.],
        vec![3., 0., 4.],
        vec![4., 1., 2.],
        vec![1., 3., 5.],
    ];
    let fit = run("linear_model", &data, &[]).unwrap();
    close(numeric(&fit, "intercept"), 1.6688632619439865, 1e-9);
    close(array(&fit, "coefficients")[0], 0.3476112026359146, 1e-9);
    close(array(&fit, "coefficients")[1], 0.26688632619439845, 1e-9);
    close(numeric(&fit, "sse"), 15.156507413509061, 1e-9);
    close(
        numeric(&fit, "intercept_standard_error"),
        0.8271428196124827,
        1e-9,
    );
    close(
        array(&fit, "coefficient_standard_errors")[0],
        0.32957165362695556,
        1e-9,
    );
    close(
        array(&fit, "coefficient_standard_errors")[1],
        0.5288277468926229,
        1e-9,
    );
}
