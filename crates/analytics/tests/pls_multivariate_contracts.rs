use ducksmiles_analytics::execute;
fn values(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let text = &j[start..];
    let end = if text.starts_with('[') {
        let mut depth = 0;
        text.char_indices()
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
        text.find([',', '}']).unwrap_or(text.len())
    };
    text[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * (1. + b.abs()), "{a} != {b}");
}
fn data() -> Vec<Vec<f64>> {
    ROWS.into_iter().map(Vec::from).collect()
}
const ROWS: [[f64; 5]; 8] = [
    [0.0, 1.0, 2.0, 1.0, 4.0],
    [1.0, 0.0, 2.0, 2.0, 3.0],
    [2.0, 1.0, 0.0, 4.0, 1.0],
    [3.0, 2.0, 1.0, 5.0, 2.0],
    [4.0, 0.0, 3.0, 6.0, 6.0],
    [5.0, 3.0, 2.0, 8.0, 3.0],
    [6.0, 2.0, 4.0, 7.0, 9.0],
    [7.0, 4.0, 1.0, 11.0, 4.0],
];
#[test]
fn nipals_multiresponse_matches_independent_coefficients_vip_and_saved_scores() {
    let fit = execute("pls2_fit", &data(), &[2., 2., 0.]).unwrap();
    assert_eq!(values(&fit, "coefficients").len(), 6);
    assert_eq!(values(&fit, "vip").len(), 3);
    for (v, w) in values(&fit, "coefficients").into_iter().zip([
        0.7439616005501649,
        0.33621822710050764,
        1.0786347018276106,
        -0.15922830703454988,
        0.08284842205126429,
        1.6271586786547896,
    ]) {
        near(v, w);
    }
    for (v, w) in values(&fit, "vip").into_iter().zip([
        1.0808985261599124,
        0.8563459868205708,
        1.0480123696796866,
    ]) {
        near(v, w);
    }
    let query = [[1.5, 1.0, 2.0], [8.0, 3.0, 1.0]]
        .into_iter()
        .map(Vec::from)
        .collect::<Vec<_>>();
    let pred = execute("pls2_predict", &query, &values(&fit, "predict_parameters")).unwrap();
    assert_eq!(values(&pred, "predictions").len(), 4);
    for (v, w) in values(&pred, "predictions").into_iter().zip([
        3.3482861630138214,
        3.630476072527427,
        10.25845754819385,
        3.8702792559568375,
    ]) {
        near(v, w);
    }
    near(values(&fit, "components_used")[0], 2.);
}
#[test]
fn simpls_multiresponse_matches_independent_coefficients_vip_and_saved_scores() {
    let fit = execute("pls2_fit", &data(), &[2., 2., 1.]).unwrap();
    assert_eq!(values(&fit, "coefficients").len(), 6);
    assert_eq!(values(&fit, "vip").len(), 3);
    for (v, w) in values(&fit, "coefficients").into_iter().zip([
        0.7458064184373452,
        0.33240160839608535,
        1.0761300274244947,
        -0.15267961501072572,
        0.07935357763570298,
        1.6312288338972363,
    ]) {
        near(v, w);
    }
    for (v, w) in values(&fit, "vip").into_iter().zip([
        1.0749598796042106,
        0.8147426728276999,
        1.0865797873671665,
    ]) {
        near(v, w);
    }
    let query = [[1.5, 1.0, 2.0], [8.0, 3.0, 1.0]]
        .into_iter()
        .map(Vec::from)
        .collect::<Vec<_>>();
    let pred = execute("pls2_predict", &query, &values(&fit, "predict_parameters")).unwrap();
    assert_eq!(values(&pred, "predictions").len(), 4);
    for (v, w) in values(&pred, "predictions").into_iter().zip([
        3.345725093189463,
        3.6345251468266873,
        10.266373290245493,
        3.8585475374825546,
    ]) {
        near(v, w);
    }
    near(values(&fit, "components_used")[0], 2.);
}
#[test]
fn orthogonal_scores_reconstruct_centered_responses_and_vip_norm() {
    for method in [0., 1.] {
        let fit = execute("pls2_fit", &data(), &[2., 3., method]).unwrap();
        let scores = values(&fit, "x_scores");
        let k = values(&fit, "components_used")[0] as usize;
        for a in 0..k {
            for b in 0..a {
                near(
                    (0..8).map(|i| scores[i * k + a] * scores[i * k + b]).sum(),
                    0.,
                );
            }
        }
        near(values(&fit, "vip").iter().map(|v| v * v).sum(), 3.);
        let pred = execute(
            "pls2_predict",
            &data().iter().map(|r| r[..3].to_vec()).collect::<Vec<_>>(),
            &values(&fit, "predict_parameters"),
        )
        .unwrap();
        for (a, b) in values(&fit, "fitted")
            .into_iter()
            .zip(values(&pred, "predictions"))
        {
            near(a, b);
        }
    }
}
#[test]
fn one_predictor_multivariate_pls_has_closed_form_ols_mean_and_individual_limits() {
    let rows = vec![
        vec![-2., -3., 4.],
        vec![-1., -1., 2.],
        vec![0., 2., 2.],
        vec![1., 3., 1.],
        vec![2., 6., 0.],
    ];
    for method in [0., 1.] {
        let fit = execute("pls2_fit", &rows, &[2., 1., method]).unwrap();
        let pred = execute(
            "pls2_predict",
            &[vec![1.]],
            &values(&fit, "predict_parameters"),
        )
        .unwrap();
        near(values(&pred, "predictions")[0], 3.6);
        near(values(&pred, "predictions")[1], 0.9);
        near(values(&pred, "mean_standard_error")[0], 0.08_f64.sqrt());
        near(values(&pred, "mean_standard_error")[1], 0.07_f64.sqrt());
        near(
            values(&pred, "individual_standard_error")[0],
            (0.8_f64 / 3. * 1.3).sqrt(),
        );
        near(
            values(&pred, "mean_ci_upper")[0],
            3.6 + 3.182446305284263 * 0.08_f64.sqrt(),
        );
    }
}
#[test]
fn pls_discriminant_uses_sorted_numeric_labels_and_raw_indicator_scores() {
    let rows = vec![
        vec![-3., 10.],
        vec![-2., 10.],
        vec![-1., 10.],
        vec![1., 30.],
        vec![2., 30.],
        vec![3., 30.],
    ];
    for method in [0., 1.] {
        let fit = execute("plsda_fit", &rows, &[1., method]).unwrap();
        assert_eq!(values(&fit, "class_labels"), vec![10., 30.]);
        let pred = execute(
            "plsda_predict",
            &[vec![-7.], vec![0.], vec![7.]],
            &values(&fit, "predict_parameters"),
        )
        .unwrap();
        assert_eq!(values(&pred, "predicted_labels"), vec![10., 10., 30.]);
        let scores = values(&pred, "predictions");
        near(scores[0], 2.);
        near(scores[1], -1.);
        assert!(pred.contains("not_calibrated_probabilities"));
    }
}
#[test]
fn predictive_rank_exhaustion_and_constant_responses_remain_valid() {
    let rows = vec![
        vec![0., 0., 2., 5.],
        vec![1., 2., 4., 5.],
        vec![2., 4., 6., 5.],
    ];
    for method in [0., 1.] {
        let fit = execute("pls2_fit", &rows, &[2., 2., method]).unwrap();
        near(values(&fit, "components_used")[0], 1.);
        let pred = execute(
            "pls2_predict",
            &[vec![3., 6.]],
            &values(&fit, "predict_parameters"),
        )
        .unwrap();
        near(values(&pred, "predictions")[0], 8.);
        near(values(&pred, "predictions")[1], 5.);
    }
}
#[test]
fn large_offsets_reuse_training_transforms_without_cancellation() {
    let rows = (0..8)
        .map(|i| {
            let x = i as f64;
            vec![1e12 + x, -1e12 + 2. * x, 1e12 + 3. * x]
        })
        .collect::<Vec<_>>();
    let fit = execute("pls2_fit", &rows, &[2., 1.]).unwrap();
    let pred = execute(
        "pls2_predict",
        &[vec![1e12 + 10.]],
        &values(&fit, "predict_parameters"),
    )
    .unwrap();
    assert!((values(&pred, "predictions")[0] - (-1e12 + 20.)).abs() < 0.001);
    assert!((values(&pred, "predictions")[1] - (1e12 + 30.)).abs() < 0.001);
}
#[test]
fn malformed_inputs_carriers_and_wrong_model_kind_are_rejected() {
    let fit = execute("pls2_fit", &data(), &[2., 2.]).unwrap();
    let saved = values(&fit, "predict_parameters");
    for bad in [
        vec![],
        vec![0.],
        vec![2., 0.],
        vec![2., 2., 2.],
        vec![2., 2., 0., 2.],
        vec![2., 2., 0., 1., 1., 1.],
    ] {
        assert!(execute("pls2_fit", &data(), &bad).is_err());
    }
    assert!(execute(
        "pls2_predict",
        &[vec![1., 2., 3.]],
        &saved[..saved.len() - 1]
    )
    .is_err());
    assert!(execute("plsda_predict", &[vec![1., 2., 3.]], &saved).is_err());
    assert!(execute("plsda_fit", &[vec![0., 1.], vec![1., 1.]], &[1.]).is_err());
    let mut bad = saved.clone();
    bad[0] = 2.;
    assert!(execute("pls2_predict", &[vec![1., 2., 3.]], &bad).is_err());
}

#[test]
fn no_predictive_covariance_returns_intercept_with_undefined_vip() {
    for method in [0., 1.] {
        let fit = execute(
            "pls2_fit",
            &[vec![0., 2.], vec![0., 4.], vec![0., 6.]],
            &[1., 1., method],
        )
        .unwrap();
        near(values(&fit, "components_used")[0], 0.);
        assert!(fit.contains("\"vip\":null"));
        let pred = execute(
            "pls2_predict",
            &[vec![100.]],
            &values(&fit, "predict_parameters"),
        )
        .unwrap();
        near(values(&pred, "predictions")[0], 4.);
        near(values(&pred, "mean_standard_error")[0], (4_f64 / 3.).sqrt());
    }
}

#[test]
fn saturated_two_row_model_withholds_undefined_residual_variance() {
    let fit = execute("pls2_fit", &[vec![0., 1.], vec![1., 3.]], &[1., 1.]).unwrap();
    let pred = execute(
        "pls2_predict",
        &[vec![2.]],
        &values(&fit, "predict_parameters"),
    )
    .unwrap();
    near(values(&pred, "predictions")[0], 5.);
    assert!(pred.contains("\"mean_standard_error\":null"));
}
