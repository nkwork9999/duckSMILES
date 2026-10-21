use ducksmiles_analytics::execute;

fn number(j: &str, key: &str) -> f64 {
    let s = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    j[s..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn values(j: &str, key: &str) -> Vec<f64> {
    let s = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let e = j[s..].find(']').unwrap() + s;
    if e == s {
        return vec![];
    }
    j[s..e].split(',').map(|x| x.parse().unwrap()).collect()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1. + b.abs()), "{a} != {b}");
}
fn paths() -> Vec<Vec<f64>> {
    vec![
        vec![0., 10.1, 1.],
        vec![1., 7.9, 1.],
        vec![2., 5.9, 1.],
        vec![3., 4.1, 1.],
        vec![0., 12.1, 2.],
        vec![1., 9.9, 2.],
        vec![2., 7.9, 2.],
        vec![3., 6.1, 2.],
    ]
}
#[test]
fn linear_unit_paths_recover_separate_intercepts_and_common_slope() {
    let j = execute("degradation_path", &paths(), &[0., 1., 0., 0.]).unwrap();
    assert_eq!(values(&j, "intercepts").len(), 2);
    close(values(&j, "intercepts")[0], 10., 1e-12);
    close(values(&j, "intercepts")[1], 12., 1e-12);
    close(values(&j, "slopes")[0], -2., 1e-12);
    close(number(&j, "sse"), 0.08, 1e-10);
    close(number(&j, "residual_df"), 5., 1e-12);
}
#[test]
fn common_intercept_is_shared_at_time_zero_when_unit_slopes_differ() {
    let data = vec![
        vec![0., 10.1, 1.],
        vec![1., 8.9, 1.],
        vec![2., 7.9, 1.],
        vec![3., 7.1, 1.],
        vec![0., 10.1, 2.],
        vec![1., 7.9, 2.],
        vec![2., 5.9, 2.],
        vec![3., 4.1, 2.],
    ];
    let j = execute("degradation_path", &data, &[1., 0., 0., 0.]).unwrap();
    close(number(&j, "sse"), 0.08, 1e-10);
    close(values(&j, "intercepts")[0], 10., 1e-10);
    close(values(&j, "intercepts")[1], 10., 1e-10);
    close(values(&j, "slopes")[0], -1., 1e-10);
    close(values(&j, "slopes")[1], -2., 1e-10);
}
#[test]
fn saved_path_prediction_has_exact_threshold_and_mean_uncertainty() {
    // v1,response-transform,time-transform,G; ID,a,b,covaa,covab,covbb,MSE,df.
    let model = [1., 0., 0., 1., 7., 10., -2., 0.04, -0.01, 0.01, 0.16, 8.];
    let mut args = vec![0.95, 0.];
    args.extend(model);
    let j = execute("degradation_predict", &[vec![2., 7., 6.]], &args).unwrap();
    close(values(&j, "fitted")[0], 6., 1e-12);
    close(values(&j, "threshold_times")[0], 2., 1e-12);
    assert!(j.contains("\"inverse_confidence_sets\":[["));
}
#[test]
fn inverse_path_confidence_sets_preserve_empty_two_piece_and_unbounded_cases() {
    let mut args = vec![0.95, 0., 1., 0., 0., 1., 7., 0., 0., 1., 0., 0., 1., 8.];
    let empty = execute("degradation_predict", &[vec![1., 7., 10.]], &args).unwrap();
    assert!(
        empty.contains("\"inverse_confidence_sets\":[[]]"),
        "{empty}"
    );
    let all = execute("degradation_predict", &[vec![1., 7., 0.]], &args).unwrap();
    assert!(all.contains("\"upper\":null"), "{all}");
    args = vec![0.95, 0., 1., 0., 1., 1., 7., 2., 0., 0., 0., 1., 1., 8.];
    let split = execute("degradation_predict", &[vec![1., 7., 0.]], &args).unwrap();
    assert_eq!(split.matches("\"lower\":").count(), 2, "{split}");
    assert!(split.contains("\"lower_closed\":false"), "{split}");
    assert!(split.contains("\"upper\":null"), "{split}");
}
#[test]
fn stability_pooling_tests_slopes_before_intercepts_and_defaults_to_separate_mse() {
    let j = execute(
        "stability_analysis",
        &paths(),
        &[1., 3., 0., 0., 0.25, 0.95, 0., 20.],
    )
    .unwrap();
    assert!(
        j.contains("\"selected_model\":\"separate_intercepts_common_slope\""),
        "{j}"
    );
    assert!(j.contains("\"test_order\":[\"slopes\",\"intercepts\"]"));
    assert!(number(&j, "slope_pooling_p") > 0.99);
    assert!(number(&j, "intercept_pooling_p") < 0.001);
    assert!(number(&j, "shelf_life") > 3. && number(&j, "shelf_life") < 3.5);
}
#[test]
fn independent_destructive_gaussian_fit_keeps_density_constant_and_mle_scale() {
    let rows = vec![
        vec![0., 0., 10.1, 10.1, 0.],
        vec![1., 0., 7.9, 7.9, 0.],
        vec![2., 0., 5.9, 5.9, 0.],
        vec![3., 0., 4.1, 4.1, 0.],
    ];
    let j = execute("destructive_degradation", &rows, &[]).unwrap();
    close(values(&j, "coefficients")[0], 10., 1e-7);
    close(values(&j, "coefficients")[1], -2., 1e-7);
    close(number(&j, "sigma"), 0.1, 1e-7);
    close(number(&j, "log_likelihood"), 3.534586239157491, 1e-7);
}
#[test]
fn saved_destructive_lognormal_prediction_differentiates_mean_median_and_soft_failure() {
    // v1,family,time-transform,stress-mode,sigma,centered-intercept,p;
    // feature origins,mean offsets,raw slopes (one time feature).
    let model = [1., 1., 0., 0., 0.5, 2., 1., 0., 0., -0.2];
    let j = execute(
        "destructive_predict",
        &[vec![5., 0., 1.0_f64.exp()]],
        &model,
    )
    .unwrap();
    close(values(&j, "medians")[0], 1.0_f64.exp(), 1e-12);
    close(values(&j, "means")[0], 1.125_f64.exp(), 1e-12);
    close(values(&j, "lower_failure_probabilities")[0], 0.5, 1e-12);
}
#[test]
fn constant_stress_cumulative_exposure_reduces_to_exponential_total_exposure() {
    // distribution=exponential,relation=linear,continuation=extend; two patterns
    // with distinct stress to identify slope. Each lifetime data [lo,hi,type,pattern].
    let data = vec![
        vec![1., 1., 0., 1.],
        vec![3., 3., 0., 1.],
        vec![2., 2., 0., 2.],
        vec![6., 6., 0., 2.],
    ];
    let args = [0., 0., 1., 300., 1e-9, 2., 1., 1., 10., 0., 2., 1., 10., 1.];
    let j = execute("cumulative_damage_fit", &data, &args).unwrap();
    close(values(&j, "coefficients")[0], 2.0_f64.ln(), 1e-6);
    close(values(&j, "coefficients")[1], 2.0_f64.ln(), 1e-6);
    close(number(&j, "log_likelihood"), -8.158883083359672, 1e-6);
}
#[test]
fn saved_step_stress_cumulative_prediction_integrates_scale_and_jacobian() {
    // v1,dist(exp),relation(linear),continuation(extend),b0,b1,sigma,patterncount;
    // id,nsegments,duration,stress,...
    let model = [
        1.,
        0.,
        0.,
        1.,
        0.,
        2.0_f64.ln(),
        1.,
        1.,
        7.,
        2.,
        2.,
        0.,
        3.,
        1.,
    ];
    let j = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &model).unwrap();
    close(values(&j, "effective_ages")[0], 3., 1e-12);
    close(values(&j, "survivals")[0], (-3.0_f64).exp(), 1e-12);
    close(values(&j, "hazards")[0], 0.5, 1e-12);
    close(values(&j, "quantiles")[0], 2.0_f64.ln(), 1e-12);
}
#[test]
fn cumulative_saved_weibull_and_lognormal_use_the_exposure_distribution() {
    let mut model = vec![
        1.,
        1.,
        0.,
        1.,
        0.,
        2.0_f64.ln(),
        0.5,
        1.,
        7.,
        2.,
        2.,
        0.,
        3.,
        1.,
    ];
    let j = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &model).unwrap();
    close(values(&j, "survivals")[0], (-9.0_f64).exp(), 1e-12);
    close(values(&j, "hazards")[0], 3., 1e-12);
    close(values(&j, "quantiles")[0], 2.0_f64.ln().sqrt(), 1e-12);
    model[1] = 2.;
    model[6] = 1.;
    let j = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &model).unwrap();
    close(values(&j, "survivals")[0], 0.13596860764142443, 1e-11);
    close(values(&j, "quantiles")[0], 1., 1e-11);
}
#[test]
fn cumulative_zero_time_weibull_density_has_correct_shape_boundary() {
    let mut model = vec![1., 1., 0., 1., 0., 0., 1., 1., 7., 1., 2., 0.];
    let j = execute("cumulative_damage_predict", &[vec![0., 7., 0.5]], &model).unwrap();
    close(values(&j, "densities")[0], 1., 1e-12);
    model[6] = 2.;
    let j = execute("cumulative_damage_predict", &[vec![0., 7., 0.5]], &model).unwrap();
    assert!(j.contains("\"densities\":[null]"), "{j}");
}
#[test]
fn unused_stress_patterns_cannot_create_identified_damage_fit() {
    let data = vec![
        vec![1., 1., 0., 1.],
        vec![2., 2., 0., 1.],
        vec![3., 3., 0., 1.],
        vec![4., 4., 0., 1.],
    ];
    let args = [0., 0., 1., 300., 1e-8, 2., 1., 1., 10., 0., 2., 1., 10., 1.];
    if let Ok(j) = execute("cumulative_damage_fit", &data, &args) {
        assert!(j.contains("\"converged\":false"), "{j}");
    }
}
#[test]
fn destructive_likelihood_uses_left_right_and_interval_response_probabilities() {
    let rows = vec![
        vec![0., 0., 10.1, 10.1, 0.],
        vec![1., 0., 7.9, 7.9, 0.],
        vec![2., 0., 5.9, 5.9, 0.],
        vec![3., 0., 4.1, 4.1, 0.],
        vec![4., 0., 0., 2.2, 1.],
        vec![0., 0., 9.8, 0., 2.],
        vec![2., 0., 5.5, 6.5, 3.],
    ];
    let j = execute("destructive_degradation", &rows, &[]).unwrap();
    close(number(&j, "log_likelihood"), 3.4947580021678655, 1e-8);
    close(number(&j, "sigma"), 0.09784559481146529, 1e-6);
    assert!(j.contains("\"converged\":true"), "{j}");
}
#[test]
fn cumulative_log_and_reciprocal_stress_relations_integrate_transformed_location() {
    let mut model = vec![
        1.,
        0.,
        1.,
        1.,
        0.,
        1.,
        1.,
        1.,
        7.,
        2.,
        2.,
        1.,
        3.,
        1.0_f64.exp(),
    ];
    let j = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &model).unwrap();
    close(
        values(&j, "effective_ages")[0],
        2. + 2. * (-1.0_f64).exp(),
        1e-12,
    );
    model[2] = 3.;
    model[13] = 2.;
    let j = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &model).unwrap();
    close(
        values(&j, "effective_ages")[0],
        2. * (-1.0_f64).exp() + 2. * (-0.5_f64).exp(),
        1e-12,
    );
}
#[test]
fn cumulative_quantiles_remain_relative_accurate_far_below_first_segment_duration() {
    let model = [1., 0., 0., 1., -200., 0., 1., 1., 7., 1., 1., 0.];
    let j = execute("cumulative_damage_predict", &[vec![1e-90, 7., 0.5]], &model).unwrap();
    let actual = values(&j, "quantiles")[0];
    let expected = (-200.0_f64).exp() * 2.0_f64.ln();
    assert!(
        (actual / expected - 1.).abs() < 1e-11,
        "{actual} / {expected}"
    );
}
#[test]
fn stopped_stress_quantiles_do_not_round_pattern_duration_outside_its_domain() {
    let model = [1., 0., 0., 0., 0., 0., 1., 1., 7., 1., 10., 0.];
    let j = execute("cumulative_damage_predict", &[vec![1., 7., 0.5]], &model).unwrap();
    close(values(&j, "quantiles")[0], 2.0_f64.ln(), 1e-12);
}
#[test]
fn tiny_repeating_patterns_bracket_quantiles_without_overflowing_cycle_count() {
    let model = [
        1.,
        0.,
        0.,
        2.,
        1e-300_f64.ln(),
        0.,
        1.,
        1.,
        7.,
        1.,
        1e-300,
        0.,
    ];
    let j = execute(
        "cumulative_damage_predict",
        &[vec![1e-301, 7., 0.9]],
        &model,
    )
    .unwrap();
    let actual = values(&j, "quantiles")[0];
    let expected = 1e-300 * 10.0_f64.ln();
    assert!(
        (actual / expected - 1.).abs() < 1e-11,
        "{actual} / {expected}"
    );
}
#[test]
fn saved_path_rejects_impossible_residual_df_instead_of_collapsing_uncertainty() {
    let model = [0.95, 0., 1., 0., 0., 1., 7., 0., 1., 1., 0., 0., 1., 1e308];
    assert!(execute("degradation_predict", &[vec![1., 7., 0.]], &model).is_err());
}
#[test]
fn near_unit_confidence_uses_the_small_student_tail_without_rounding_to_one() {
    let confidence = f64::from_bits(1.0_f64.to_bits() - 1);
    let model = [
        confidence, 0., 1., 0., 0., 1., 7., 0., 1., 1., 0., 0., 1., 1.,
    ];
    let j = execute("degradation_predict", &[vec![0., 7., 0.]], &model).unwrap();
    let expected = (std::f64::consts::PI * ((1. - confidence) / 2.))
        .tan()
        .recip();
    let actual = values(&j, "mean_upper")[0];
    assert!(
        (actual / expected - 1.).abs() < 1e-11,
        "{actual} / {expected}"
    );
}
#[test]
fn malformed_carriers_censor_intervals_and_excessive_resources_reject_without_panics() {
    for op in [
        "degradation_path",
        "degradation_predict",
        "stability_analysis",
        "destructive_degradation",
        "destructive_predict",
        "cumulative_damage_fit",
        "cumulative_damage_predict",
    ] {
        assert!(execute(op, &[], &[]).is_err());
        assert!(execute(op, &[vec![f64::NAN]], &[]).is_err());
        assert!(execute(op, &[vec![1.], vec![1., 2.]], &[]).is_err());
        assert!(execute(op, &[vec![1.]], &[f64::INFINITY]).is_err());
    }
    let excessive = vec![vec![0., 0., 1.]; 10001];
    assert!(execute("degradation_path", &excessive, &[]).is_err());
    assert!(execute("degradation_path", &paths(), &[0., 0., 0., 0., 0.]).is_err());
    let bad_cov = [0.95, 0., 1., 0., 0., 1., 7., 0., 1., 1., 2., 1., 1., 8.];
    assert!(execute("degradation_predict", &[vec![1., 7., 0.]], &bad_cov).is_err());
    assert!(execute("destructive_degradation", &[vec![1., 0., 2., 1., 3.]], &[]).is_err());
    assert!(execute(
        "destructive_predict",
        &[vec![1., 0., 2.]],
        &[1., 0., 0., 0., -1., 0., 1., 0., 0., 1.]
    )
    .is_err());
    let model = [1., 0., 0., 0., 0., 0., 1., 1., 7., 1., 2., 0.];
    assert!(execute("cumulative_damage_predict", &[vec![3., 7., 0.5]], &model).is_err());
    assert!(execute("cumulative_damage_predict", &[vec![1., 7., 1.]], &model).is_err());
    assert!(execute("cumulative_damage_predict", &[vec![1., 8., 0.5]], &model).is_err());
    let missing = [1., 0., 0., 0., 0., 0., 1., 1., 7., 2., 2., 0.];
    assert!(execute("cumulative_damage_predict", &[vec![1., 7., 0.5]], &missing).is_err());
}
#[test]
fn changing_time_units_preserves_cumulative_survival_and_rescales_density() {
    let original = [
        1.,
        0.,
        0.,
        1.,
        0.,
        2.0_f64.ln(),
        1.,
        1.,
        7.,
        2.,
        2.,
        0.,
        3.,
        1.,
    ];
    let scaled = [
        1.,
        0.,
        0.,
        1.,
        10.0_f64.ln(),
        2.0_f64.ln(),
        1.,
        1.,
        7.,
        2.,
        20.,
        0.,
        30.,
        1.,
    ];
    let a = execute("cumulative_damage_predict", &[vec![4., 7., 0.5]], &original).unwrap();
    let b = execute("cumulative_damage_predict", &[vec![40., 7., 0.5]], &scaled).unwrap();
    close(
        values(&a, "survivals")[0],
        values(&b, "survivals")[0],
        1e-12,
    );
    close(
        values(&a, "densities")[0],
        10. * values(&b, "densities")[0],
        1e-12,
    );
    close(
        values(&a, "hazards")[0],
        10. * values(&b, "hazards")[0],
        1e-12,
    );
    close(
        10. * values(&a, "quantiles")[0],
        values(&b, "quantiles")[0],
        1e-12,
    );
}
#[test]
fn nonconverged_destructive_iterate_with_positive_information_has_no_mle_inference() {
    let rows = vec![
        vec![0., 0., 10.1, 10.1, 0.],
        vec![1., 0., 7.9, 7.9, 0.],
        vec![2., 0., 5.9, 5.9, 0.],
        vec![3., 0., 4.1, 4.1, 0.],
        vec![4., 0., 0., 2.2, 1.],
        vec![0., 0., 9.8, 0., 2.],
        vec![2., 0., 5.5, 6.5, 3.],
    ];
    let j = execute("destructive_degradation", &rows, &[0., 0., 0., 1., 1e-12]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(
        j.contains("\"standardized_parameter_covariance\":null"),
        "{j}"
    );
    assert!(j.contains("\"aic\":null"), "{j}");
    assert!(j.contains("\"inference_valid\":false"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
    assert!(number(&j, "log_likelihood").is_finite());
    assert_eq!(values(&j, "coefficients").len(), 2);
    let prediction = execute(
        "destructive_predict",
        &[vec![2., 0., 6.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert!(values(&prediction, "means")[0].is_finite());
}
#[test]
fn nonconverged_damage_iterate_with_positive_information_has_no_mle_inference() {
    let data = vec![
        vec![1., 1., 0., 1.],
        vec![3., 3., 0., 1.],
        vec![2., 2., 0., 2.],
        vec![6., 6., 0., 2.],
    ];
    let args = [0., 0., 1., 1., 1e-12, 2., 1., 1., 10., 0., 2., 1., 10., 1.];
    let j = execute("cumulative_damage_fit", &data, &args).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"centered_parameter_covariance\":null"), "{j}");
    assert!(j.contains("\"aic\":null"), "{j}");
    assert!(j.contains("\"inference_valid\":false"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
    assert!(number(&j, "log_likelihood").is_finite());
    assert_eq!(values(&j, "coefficients").len(), 2);
    let prediction = execute(
        "cumulative_damage_predict",
        &[vec![2., 1., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert!(values(&prediction, "survivals")[0].is_finite());
}
#[test]
fn destructive_scale_guard_never_allows_mle_inference_from_a_projected_stationary_point() {
    let data = vec![
        vec![0., 0., 10., 10., 0.],
        vec![1., 0., 8., 8., 0.],
        vec![2., 0., 6., 6., 0.],
        vec![3., 0., 4., 4., 0.],
    ];
    let j = execute("destructive_degradation", &data, &[0., 0., 0., 300., 1e-8]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"aic\":null"), "{j}");
    assert!(j.contains("\"parameter_guard_bound\":true"), "{j}");
    assert!(j.contains("\"inference_valid\":false"), "{j}");
    assert!(
        j.contains("\"standardized_parameter_covariance\":null"),
        "{j}"
    );
}
