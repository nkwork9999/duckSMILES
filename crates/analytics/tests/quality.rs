use ducksmiles_analytics::quality::run;
fn value(json: &str, key: &str) -> f64 {
    json.split_once(&format!("\"{key}\":"))
        .unwrap()
        .1
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn array(json: &str, key: &str) -> Vec<f64> {
    json.split_once(&format!("\"{key}\":["))
        .unwrap()
        .1
        .split(']')
        .next()
        .unwrap()
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-9, "{a} != {b}");
}
#[test]
fn individuals_chart_estimates_within_sigma_from_moving_ranges() {
    let j = run("ir_chart", &[vec![1.], vec![2.], vec![3.], vec![4.]], &[]).unwrap();
    close(value(&j, "center"), 2.5);
    close(value(&j, "sigma"), std::f64::consts::PI.sqrt() / 2.);
    close(
        array(&j, "upper")[0],
        2.5 + 1.5 * std::f64::consts::PI.sqrt(),
    );
}
#[test]
fn subgroup_charts_use_ranges_or_unbiased_standard_deviations() {
    let x = vec![vec![1., 3.], vec![2., 4.], vec![3., 5.]];
    let j = run("xbar_r_chart", &x, &[]).unwrap();
    close(value(&j, "center"), 3.);
    close(value(&j, "sigma"), 2. / 1.128);
    let j = run("xbar_s_chart", &x, &[]).unwrap();
    close(value(&j, "sigma"), std::f64::consts::PI.sqrt());
}
#[test]
fn attribute_charts_use_binomial_and_poisson_limits_with_exposures() {
    let j = run("p_chart", &[vec![2., 100.], vec![4., 100.]], &[]).unwrap();
    close(value(&j, "center"), 0.03);
    close(
        array(&j, "upper")[0],
        0.03 + 3. * (0.03_f64 * 0.97 / 100.).sqrt(),
    );
    let j = run("np_chart", &[vec![2., 100.], vec![4., 100.]], &[]).unwrap();
    close(value(&j, "center"), 3.);
    let j = run("c_chart", &[vec![4.], vec![6.]], &[]).unwrap();
    close(value(&j, "center"), 5.);
    let j = run("u_chart", &[vec![2., 2.], vec![6., 2.]], &[]).unwrap();
    close(value(&j, "center"), 2.);
    close(array(&j, "upper")[0], 5.);
}
#[test]
fn sequential_charts_track_time_varying_uncertainty_and_accumulated_shift() {
    let x = vec![vec![1.], vec![2.], vec![3.]];
    let j = run("ewma_chart", &x, &[0.5, 3., 0., 1.]).unwrap();
    assert_eq!(array(&j, "statistics"), vec![0.5, 1.25, 2.125]);
    close(array(&j, "upper")[0], 1.5);
    let j = run("cusum_chart", &x, &[0.5, 2., 0., 1.]).unwrap();
    assert_eq!(array(&j, "positive"), vec![0.5, 2., 4.5]);
    assert_eq!(array(&j, "signals"), vec![0., 0., 1.]);
    let j = run("moving_average_chart", &x, &[2., 3., 0., 1.]).unwrap();
    assert_eq!(array(&j, "statistics"), vec![1., 1.5, 2.5]);
    close(array(&j, "upper")[1], 3. / 2_f64.sqrt());
}
#[test]
fn western_electric_rules_detect_each_rule_without_future_data() {
    let j = run(
        "western_electric",
        &[vec![2.1], vec![2.2], vec![0.]],
        &[0., 1.],
    )
    .unwrap();
    assert_eq!(array(&j, "rule_masks"), vec![0., 0., 2.]);
    let j = run("western_electric", &vec![vec![0.1]; 8], &[0., 1.]).unwrap();
    assert_eq!(array(&j, "rule_masks")[7], 8.);
    let j = run("western_electric", &[vec![3.1]], &[0., 1.]).unwrap();
    assert_eq!(array(&j, "rule_masks"), vec![1.]);
}
#[test]
fn crossed_gage_rr_separates_repeatability_operator_and_part_variance() {
    let x = vec![
        vec![1., 0., 0.],
        vec![3., 0., 0.],
        vec![2., 0., 1.],
        vec![4., 0., 1.],
        vec![5., 1., 0.],
        vec![7., 1., 0.],
        vec![6., 1., 1.],
        vec![8., 1., 1.],
    ];
    let j = run("gage_rr", &x, &[]).unwrap();
    close(value(&j, "repeatability_variance"), 2.);
    close(value(&j, "operator_variance"), 0.5);
    close(value(&j, "part_variance"), 8.);
    close(value(&j, "gage_variance"), 2.5);
}
#[test]
fn type1_gage_and_kappa_account_for_tolerance_and_chance_agreement() {
    let j = run("type1_gage", &[vec![9.], vec![10.], vec![11.]], &[10., 6.]).unwrap();
    close(value(&j, "cg"), 0.2);
    close(value(&j, "cgk"), 0.2);
    let j = run("agreement_kappa", &[vec![20., 5.], vec![5., 20.]], &[]).unwrap();
    close(value(&j, "kappa"), 0.6);
    let j = run(
        "acceptance_oc",
        &[vec![0.], vec![0.1], vec![1.]],
        &[10., 0.],
    )
    .unwrap();
    let p = array(&j, "acceptance_probability");
    close(p[0], 1.);
    close(p[1], 0.9_f64.powi(10));
    close(p[2], 0.);
}
#[test]
fn malformed_or_unidentified_quality_requests_fail() {
    assert!(run("p_chart", &[vec![4., 3.]], &[]).is_err());
    assert!(run("np_chart", &[vec![1., 10.], vec![1., 20.]], &[]).is_err());
    assert!(run("ir_chart", &[vec![1.], vec![1.]], &[]).is_err());
    assert!(run("ewma_chart", &[vec![1.]], &[0., 3., 0., 1.]).is_err());
    assert!(run("gage_rr", &[vec![1., 0., 0.], vec![2., 1., 1.]], &[]).is_err());
    assert!(run("agreement_kappa", &[vec![1., 0.], vec![0., 0.]], &[]).is_err());
}

#[test]
fn quality_rejects_derived_overflow_and_estimated_limits_are_translation_equivariant() {
    assert!(run("cusum_chart", &[vec![1e100]], &[0.5, 5., 0., 1e-300]).is_err());
    let a = run("ir_chart", &[vec![1.], vec![2.], vec![4.]], &[]).unwrap();
    let b = run("ir_chart", &[vec![101.], vec![102.], vec![104.]], &[]).unwrap();
    close(value(&a, "sigma"), value(&b, "sigma"));
    close(value(&b, "center") - value(&a, "center"), 100.);
}

#[test]
fn levey_jennings_preserves_fixed_reference_limits_and_laney_adjusts_overdispersion() {
    let j = run("levey_jennings", &[vec![0.], vec![3.1]], &[0., 1.]).unwrap();
    assert_eq!(array(&j, "signals"), vec![0., 1.]);
    let j = run("laney_u_chart", &[vec![1., 1.], vec![9., 1.]], &[]).unwrap();
    // Standardized difference 8/sqrt(5), corrected process sigma = 8/d2.
    close(value(&j, "center"), 5.);
    close(
        array(&j, "upper")[0],
        5. + 3. * 8. * std::f64::consts::PI.sqrt() / 2.,
    );
}

#[test]
fn ewma_transient_limits_preserve_tiny_valid_smoothing_parameters() {
    let j = run("ewma_chart", &[vec![2.], vec![4.]], &[1e-20, 3., 0., 1.]).unwrap();
    let upper = array(&j, "upper");
    assert!((upper[0] / 3e-20 - 1.).abs() < 1e-12);
    assert_eq!(array(&j, "signals"), vec![0., 1.]);
}
