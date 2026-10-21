use ducksmiles_analytics::execute;

fn number(j: &str, key: &str) -> f64 {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    j[start..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn values(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &j[start..];
    let end = tail.find(']').unwrap();
    if end == 0 {
        return vec![];
    }
    tail[..end].split(',').map(|s| s.parse().unwrap()).collect()
}
fn matrix(j: &str, key: &str) -> Vec<Vec<f64>> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &j[start..];
    if tail.starts_with(']') {
        return vec![];
    }
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
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
fn check_array(a: &[f64], b: &[f64], t: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y, t);
    }
}

#[test]
fn bernoulli_random_intercept_zero_variance_matches_glm_likelihood_and_information() {
    let rows = (0..6)
        .flat_map(|g| (0..6).map(move |i| vec![f64::from(i >= 3), g as f64]))
        .collect::<Vec<_>>();
    let j = execute("glmm_binomial", &rows, &[]).unwrap();
    close(number(&j, "intercept"), 0., 1e-7);
    close(number(&j, "random_intercept_variance"), 0., 1e-10);
    close(number(&j, "log_likelihood"), 36. * 0.5_f64.ln(), 1e-8);
    close(number(&j, "intercept_standard_error"), 1. / 3., 2e-5);
    check_array(&values(&j, "posterior_means"), &[0.; 6], 1e-10);
    check_array(&values(&j, "posterior_variances"), &[0.; 6], 1e-10);
    assert!(j.contains("\"boundary\":true"));
    assert!(j.contains("\"converged\":true"), "{j}");
}

#[test]
fn poisson_random_intercept_zero_variance_includes_factorial_constants_and_glm_se() {
    let rows = (0..6)
        .flat_map(|g| [1., 2., 3., 2.].into_iter().map(move |y| vec![y, g as f64]))
        .collect::<Vec<_>>();
    let j = execute("glmm_poisson", &rows, &[]).unwrap();
    close(number(&j, "intercept"), 2.0_f64.ln(), 1e-7);
    close(number(&j, "random_intercept_variance"), 0., 1e-10);
    let expected = 6. * (8. * 2.0_f64.ln() - 8. - (2. * 2.0_f64.ln() + 6.0_f64.ln()));
    close(number(&j, "log_likelihood"), expected, 1e-8);
    close(
        number(&j, "intercept_standard_error"),
        1. / 48.0_f64.sqrt(),
        2e-5,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
}

#[test]
fn saved_glmm_prediction_separates_conditional_posterior_and_population_means() {
    let d = 0.3_f64.sqrt();
    // Header/version,family,p,G,Q,variance,eta_center, then groups/posterior moments/nodes/weights.
    let model = [
        1.,
        1.,
        0.,
        1.,
        3.,
        0.5,
        2.0_f64.ln(),
        7.,
        0.2,
        0.2,
        0.1,
        0.2 - d,
        0.2,
        0.2 + d,
        1. / 6.,
        2. / 3.,
        1. / 6.,
    ];
    let j = execute("glmm_predict", &[vec![7.], vec![99.]], &model).unwrap();
    check_array(
        &values(&j, "conditional_predictions"),
        &[2. * 0.2_f64.exp(), 2.],
        1e-10,
    );
    check_array(
        &values(&j, "posterior_predictions"),
        &[
            2. * 0.2_f64.exp() * (2. / 3. + d.cosh() / 3.),
            2. * 0.25_f64.exp(),
        ],
        1e-10,
    );
    check_array(
        &values(&j, "population_predictions"),
        &[2. * 0.25_f64.exp(); 2],
        1e-10,
    );
    assert!(j.contains("\"known_groups\":[true,false]"));
}

#[test]
fn glmm_rejects_bad_inputs_and_never_claims_a_finite_separated_estimate() {
    for op in ["glmm_binomial", "glmm_poisson", "glmm_predict"] {
        assert!(execute(op, &[], &[]).is_err());
        assert!(execute(op, &[vec![f64::NAN]], &[]).is_err());
        assert!(execute(op, &[vec![1., 2.], vec![3.]], &[]).is_err());
    }
    let separated = (0..4)
        .flat_map(|g| {
            [-1., 1.]
                .into_iter()
                .map(move |x| vec![x, f64::from(x > 0.), g as f64])
        })
        .collect::<Vec<_>>();
    if let Ok(j) = execute("glmm_binomial", &separated, &[]) {
        assert!(j.contains("\"converged\":false"), "{j}");
    }
    let all_zero = (0..3)
        .flat_map(|g| vec![vec![0., g as f64]; 3])
        .collect::<Vec<_>>();
    if let Ok(j) = execute("glmm_poisson", &all_zero, &[]) {
        assert!(j.contains("\"converged\":false"), "{j}");
    }
}

fn positive_poisson_rows() -> Vec<Vec<f64>> {
    vec![
        vec![-1.0, 0.0, 0.0],
        vec![-0.5, 0.0, 0.0],
        vec![0.0, 0.0, 0.0],
        vec![0.5, 0.0, 0.0],
        vec![1.0, 0.0, 0.0],
        vec![-1.0, 0.0, 1.0],
        vec![-0.5, 1.0, 1.0],
        vec![0.0, 0.0, 1.0],
        vec![0.5, 1.0, 1.0],
        vec![1.0, 0.0, 1.0],
        vec![-1.0, 0.0, 2.0],
        vec![-0.5, 3.0, 2.0],
        vec![0.0, 1.0, 2.0],
        vec![0.5, 0.0, 2.0],
        vec![1.0, 0.0, 2.0],
        vec![-1.0, 3.0, 3.0],
        vec![-0.5, 1.0, 3.0],
        vec![0.0, 3.0, 3.0],
        vec![0.5, 1.0, 3.0],
        vec![1.0, 2.0, 3.0],
        vec![-1.0, 0.0, 4.0],
        vec![-0.5, 0.0, 4.0],
        vec![0.0, 4.0, 4.0],
        vec![0.5, 3.0, 4.0],
        vec![1.0, 3.0, 4.0],
        vec![-1.0, 3.0, 5.0],
        vec![-0.5, 5.0, 5.0],
        vec![0.0, 4.0, 5.0],
        vec![0.5, 7.0, 5.0],
        vec![1.0, 6.0, 5.0],
    ]
}
#[test]
fn bernoulli_positive_variance_matches_independent_real_line_integration() {
    let rows = [0, 1, 2, 4, 5, 6]
        .iter()
        .enumerate()
        .flat_map(|(g, &successes)| (0..6).map(move |i| vec![f64::from(i < successes), g as f64]))
        .collect::<Vec<_>>();
    let j = execute("glmm_binomial", &rows, &[300., 1e-7, 23., 100.]).unwrap();
    close(number(&j, "intercept"), 0., 2e-5);
    close(
        number(&j, "random_intercept_variance"),
        4.325165492488126,
        2e-4,
    );
    close(number(&j, "log_likelihood"), -20.5708709482773, 2e-7);
    check_array(
        &values(&j, "posterior_means"),
        &[
            -2.6721283663129514,
            -1.4731592942402087,
            -0.6709416279569924,
            0.6709416662048762,
            1.4731593301110053,
            2.6721283954485546,
        ],
        2e-4,
    );
    check_array(
        &values(&j, "posterior_modes"),
        &[
            -2.320673256994725,
            -1.28748484224813,
            -0.5921139663252365,
            0.5921140054969823,
            1.2874848795275393,
            2.320673288051971,
        ],
        2e-4,
    );
    check_array(
        &values(&j, "posterior_variances"),
        &[
            1.5709758604225454,
            0.9342977824027511,
            0.7095916558931925,
            0.7095916546019672,
            0.9342977785582625,
            1.5709758496995272,
        ],
        2e-4,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(j.contains("\"quadrature_converged\":true"));
}
#[test]
fn poisson_positive_variance_fixed_slope_and_eb_moments_match_quad_oracle() {
    let j = execute(
        "glmm_poisson",
        &positive_poisson_rows(),
        &[300., 1e-7, 23., 100.],
    )
    .unwrap();
    close(number(&j, "intercept"), -0.004986696855306316, 3e-4);
    check_array(&values(&j, "coefficients"), &[0.23671924408457024], 3e-4);
    close(
        number(&j, "random_intercept_variance"),
        1.2844466299244024,
        3e-4,
    );
    close(number(&j, "log_likelihood"), -46.55605211708837, 2e-7);
    check_array(
        &values(&j, "posterior_means"),
        &[
            -1.625409436275414,
            -0.8019060922938205,
            -0.28040836326528584,
            0.5882898634863929,
            0.5882898634863929,
            1.5311442806602342,
        ],
        3e-4,
    );
    check_array(
        &values(&j, "posterior_modes"),
        &[
            -1.4780334576493455,
            -0.6879787029113644,
            -0.19489222614506102,
            0.6335818083462794,
            0.6335818083462794,
            1.5509530086153533,
        ],
        3e-4,
    );
    check_array(
        &values(&j, "posterior_variances"),
        &[
            0.5300936978649464,
            0.3186387975758752,
            0.21435114322137372,
            0.10096456878385252,
            0.10096456878385253,
            0.04145635766779107,
        ],
        3e-4,
    );
    let covariance = matrix(&j, "fixed_effect_covariance");
    assert_eq!(covariance.len(), 2);
    assert!(covariance[0][0] > 0. && covariance[1][1] > 0.);
    assert!(covariance[0][0] * covariance[1][1] > covariance[0][1] * covariance[1][0]);
    assert!(j.contains("\"converged\":true"), "{j}");
}
#[test]
fn glmm_reports_iteration_and_quadrature_limits_without_false_success() {
    let j = execute(
        "glmm_poisson",
        &positive_poisson_rows(),
        &[1., 1e-8, 15., 100.],
    )
    .unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    let rows = [0, 1, 2, 4, 5, 6]
        .iter()
        .enumerate()
        .flat_map(|(g, &successes)| (0..6).map(move |i| vec![f64::from(i < successes), g as f64]))
        .collect::<Vec<_>>();
    let j = execute("glmm_binomial", &rows, &[300., 1e-8, 3., 100.]).unwrap();
    assert!(j.contains("\"quadrature_converged\":false"));
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(number(&j, "log_likelihood_order_delta") > 1e-6);
}
#[test]
fn glmm_saved_model_preserves_predictions_after_large_feature_origin_and_units() {
    let rows = positive_poisson_rows();
    let shifted = rows
        .iter()
        .map(|r| vec![1e12 + 2. * r[0], r[1], r[2]])
        .collect::<Vec<_>>();
    let a = execute("glmm_poisson", &rows, &[300., 1e-7, 23., 100.]).unwrap();
    let b = execute("glmm_poisson", &shifted, &[300., 1e-7, 23., 100.]).unwrap();
    close(
        number(&a, "random_intercept_variance"),
        number(&b, "random_intercept_variance"),
        2e-6,
    );
    let pa = execute(
        "glmm_predict",
        &[vec![0., 0.], vec![0.5, 99.]],
        &values(&a, "predict_parameters"),
    )
    .unwrap();
    let pb = execute(
        "glmm_predict",
        &[vec![1e12, 0.], vec![1e12 + 1., 99.]],
        &values(&b, "predict_parameters"),
    )
    .unwrap();
    for key in [
        "conditional_predictions",
        "posterior_predictions",
        "population_predictions",
    ] {
        check_array(&values(&pa, key), &values(&pb, key), 2e-6);
    }
}

#[test]
fn glmm_rejects_unrepresentable_saved_posterior_second_moment() {
    let model = [
        1., 0., 0., 1., 3., 1., 0., 7., 0., 0., 0., -1e308, 0., 1e308, 0.25, 0.5, 0.25,
    ];
    assert!(execute("glmm_predict", &[vec![7.]], &model).is_err());
}

#[test]
fn tiny_variance_search_bound_does_not_turn_a_zero_optimum_into_an_upper_bound_hit() {
    let rows = (0..6)
        .flat_map(|g| (0..6).map(move |i| vec![f64::from(i >= 3), g as f64]))
        .collect::<Vec<_>>();
    let j = execute("glmm_binomial", &rows, &[200., 1e-6, 15., 1e-14]).unwrap();
    close(number(&j, "random_intercept_variance"), 0., 1e-15);
    assert!(j.contains("\"variance_upper_bound\":false"), "{j}");
    assert!(j.contains("\"converged\":true"), "{j}");
}

#[test]
fn poisson_saved_posterior_mean_uses_weighted_log_sum_before_exponentiation() {
    let model = [
        1., 1., 0., 1., 3., 1., 700., 7., 0., 0., 0.0001, -10., 0., 10., 0.0000005, 0.999999,
        0.0000005,
    ];
    let j = execute("glmm_predict", &[vec![7.]], &model).unwrap();
    close(
        values(&j, "conditional_predictions")[0] / 700_f64.exp(),
        1.,
        1e-12,
    );
    close(
        values(&j, "posterior_predictions")[0] / 700_f64.exp(),
        0.999999 + 0.000001 * 10_f64.cosh(),
        1e-12,
    );
    close(
        values(&j, "population_predictions")[0] / 700.5_f64.exp(),
        1.,
        1e-12,
    );
}

#[test]
fn accepted_rounded_posterior_weights_preserve_bernoulli_probability_bounds() {
    let model = [1., 0., 0., 1., 1., 1., 40., 7., 0., 0., 0., 0., 1.000000005];
    let j = execute("glmm_predict", &[vec![7.]], &model).unwrap();
    for x in values(&j, "posterior_predictions") {
        assert!((0. ..=1.).contains(&x));
    }
}

#[test]
fn rare_control_separation_is_not_hidden_by_loose_optimizer_tolerance() {
    let mut rows = Vec::new();
    for g in 0..20 {
        for _ in 0..12 {
            for x in [-1., 1.] {
                for z in [-1., 1.] {
                    rows.push(vec![x, z, f64::from(x > 0.), g as f64]);
                }
            }
        }
    }
    rows.push(vec![0.1, -1., 0., 0.]);
    if let Ok(j) = execute("glmm_binomial", &rows, &[200., 1e-3, 3., 100.]) {
        assert!(j.contains("\"converged\":false"), "{j}");
    }
}

#[test]
fn poisson_zero_only_direction_has_no_finite_extended_mle() {
    let rows = (0..6)
        .flat_map(|g| {
            (0..4).flat_map(move |_| vec![vec![-1., 1., g as f64], vec![1., 0., g as f64]])
        })
        .collect::<Vec<_>>();
    if let Ok(j) = execute("glmm_poisson", &rows, &[]) {
        assert!(j.contains("\"converged\":false"), "{j}");
    }
}
