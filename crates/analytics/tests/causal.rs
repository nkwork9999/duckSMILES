use ducksmiles_analytics::execute;
fn num(j: &str, k: &str) -> f64 {
    let p = j.find(&format!("\"{k}\":")).expect(k) + k.len() + 3;
    j[p..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn vals(j: &str, k: &str) -> Vec<f64> {
    let p = j.find(&format!("\"{k}\":[")).expect(k) + k.len() + 4;
    let s = &j[p..p + j[p..].find(']').unwrap()];
    if s.is_empty() {
        vec![]
    } else {
        s.split(',').map(|x| x.parse().unwrap()).collect()
    }
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
fn check(a: &[f64], b: &[f64], t: f64) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        near(*a, *b, t);
    }
}
fn ipw_data() -> Vec<Vec<f64>> {
    vec![
        vec![4., 1., 0.5],
        vec![6., 1., 0.75],
        vec![1., 0., 0.25],
        vec![3., 0., 0.5],
    ]
}
#[test]
fn propensity_logistic_fits_training_only_and_reuses_portable_preprocessing() {
    // Mutation caught: pooled treatment prevalence or suffix-label/preprocessing leakage.
    let d = [
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 1.],
        vec![1., 0.],
        vec![1., 1.],
        vec![1., 1.],
        vec![1., 1.],
        vec![2., -999.],
    ];
    let j = execute("propensity_fit", &d, &[8., 0., 200., 1e-10]).unwrap();
    near(num(&j, "intercept"), -3_f64.ln(), 1e-8);
    check(&vals(&j, "coefficients"), &[2. * 3_f64.ln()], 1e-8);
    assert!(j.contains("\"converged\":true"));
    let q = execute(
        "propensity_predict",
        &[vec![0.], vec![1.], vec![2.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    check(&vals(&q, "propensities"), &[0.25, 0.75, 27. / 28.], 1e-8);
    let mut poisoned = d.clone();
    poisoned[8][1] = 1e150;
    assert_eq!(
        j,
        execute("propensity_fit", &poisoned, &[8., 0., 200., 1e-10]).unwrap()
    );
}
#[test]
fn weighted_balance_uses_fixed_pretreatment_pooled_sd_and_ess() {
    // Mutation caught: recomputing SMD denominator after weighting or ignoring estimand weights.
    let d = [
        vec![0., 1., 0.5],
        vec![2., 1., 0.75],
        vec![1., 0., 0.25],
        vec![3., 0., 0.5],
    ];
    let j = execute("causal_balance", &d, &[0.]).unwrap();
    check(&vals(&j, "unweighted_smd"), &[-1. / 2_f64.sqrt()], 1e-12);
    check(&vals(&j, "weighted_smd"), &[-1.4 / 2_f64.sqrt()], 1e-12);
    check(
        &vals(&j, "effective_sample_sizes"),
        &[25. / 13., 25. / 13.],
        1e-12,
    );
    check(&vals(&j, "common_support"), &[0.5, 0.5], 0.);
    near(num(&j, "outside_common_support"), 2., 0.);
}
#[test]
fn hajek_ipw_ate_and_att_use_distinct_ratios_and_fixed_nuisance_influence_se() {
    // Mutation caught: raw Horvitz-Thompson mean or ATE weights used for ATT.
    let j = execute("ipw_effect", &ipw_data(), &[0.]).unwrap();
    near(num(&j, "estimate"), 2.6, 1e-12);
    near(num(&j, "standard_error"), 1.2288_f64.sqrt(), 1e-12);
    check(
        &vals(&j, "influence_values"),
        &[-1.92, 1.92, 1.92, -1.92],
        1e-12,
    );
    let j = execute("ipw_effect", &ipw_data(), &[1.]).unwrap();
    near(num(&j, "estimate"), 2.5, 1e-12);
    near(num(&j, "standard_error"), (12.5_f64 / 12.).sqrt(), 1e-12);
}
#[test]
fn aipw_orthogonal_ate_att_scores_do_not_treat_observed_outcomes_as_both_potential_outcomes() {
    // Supplied nuisance predictions differ by treatment; one factual outcome is observed per row.
    let d = [
        vec![4., 1., 0.5, 2., 4.],
        vec![6., 1., 0.75, 3., 6.],
        vec![1., 0., 0.25, 1., 3.],
        vec![3., 0., 0.5, 3., 5.],
    ];
    let j = execute("aipw_effect", &d, &[0.]).unwrap();
    near(num(&j, "estimate"), 2.25, 1e-12);
    near(num(&j, "standard_error"), 0.25, 1e-12);
    let j = execute("aipw_effect", &d, &[1.]).unwrap();
    near(num(&j, "estimate"), 2.5, 1e-12);
    near(num(&j, "standard_error"), (1_f64 / 6.).sqrt(), 1e-12);
    check(&vals(&j, "influence_values"), &[-1., 1., 0., 0.], 1e-12);
}
#[test]
fn propensity_matching_has_exact_tie_replacement_caliper_and_subset_estimand_contract() {
    // Mutation caught: hidden control reuse, wrong caliper metric, or naive independent-pair SE.
    let d = [
        vec![5., 1., 0.25],
        vec![8., 1., 0.75],
        vec![9., 1., 0.875],
        vec![1., 0., 0.125],
        vec![2., 0., 0.375],
        vec![4., 0., 0.875],
    ];
    let j = execute("propensity_match", &d, &[0.2, 1., 0.]).unwrap();
    check(&vals(&j, "treated_indices"), &[0., 1., 2.], 0.);
    check(&vals(&j, "control_indices"), &[3., 5., 5.], 0.);
    near(num(&j, "estimate"), 13. / 3., 1e-12);
    assert!(j.contains("\"standard_error\":null"));
    let j = execute("propensity_match", &d, &[0.2, 0., 0.]).unwrap();
    check(&vals(&j, "unmatched_treated_indices"), &[2.], 0.);
    near(num(&j, "estimate"), 4., 1e-12);
}
#[test]
fn uplift_curve_uses_assignment_weighted_actual_outcomes_and_groups_score_ties() {
    // Mutation caught: ranking treated outcomes as true individual treatment effects.
    let d = [
        vec![3., 1., 0.9, 0.5],
        vec![1., 0., 0.8, 0.5],
        vec![2., 1., 0.8, 0.5],
        vec![2., 0., 0.1, 0.5],
    ];
    let j = execute("uplift_curve", &d, &[]).unwrap();
    check(
        &vals(&j, "selection_fraction"),
        &[0., 0.25, 0.75, 1.],
        1e-12,
    );
    check(&vals(&j, "gain"), &[0., 1.5, 2., 1.], 1e-12);
    check(&vals(&j, "qini"), &[0., 1.25, 1.25, 0.], 1e-12);
    near(num(&j, "auuc"), 1.4375, 1e-12);
    near(num(&j, "qini_area"), 0.9375, 1e-12);
}
#[test]
fn causal_contracts_reject_invalid_assignment_and_unidentified_groups() {
    assert!(execute("ipw_effect", &[vec![1., 1., 0.], vec![2., 0., 0.5]], &[]).is_err());
    assert!(execute(
        "aipw_effect",
        &[vec![1., 1., 0.5, 0., 1.], vec![2., 1., 0.5, 0., 1.]],
        &[]
    )
    .is_err());
    assert!(execute(
        "propensity_match",
        &[vec![1., 1., 0.5], vec![2., 0., 0.8]],
        &[0., 1., 0.]
    )
    .is_err());
    assert!(execute(
        "causal_balance",
        &[vec![1., 1., 0.5], vec![2., 0., 0.5]],
        &[]
    )
    .is_err());
    assert!(execute("uplift_curve", &[vec![1., 2., 0.5, 0.5]], &[]).is_err());
}

#[test]
fn propensity_gradient_diagnostic_belongs_to_returned_iterate() {
    // Mutation caught: reporting the pre-update gradient at an iteration limit.
    let d = [
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 0.],
        vec![0., 1.],
        vec![1., 0.],
        vec![1., 1.],
        vec![1., 1.],
        vec![1., 1.],
    ];
    let j = execute("propensity_fit", &d, &[8., 1., 1., 1e-12]).unwrap();
    let expected = (8. * (1. / (1. + (-2_f64 / 3.).exp()) - 0.75) + 2. / 3.).abs();
    near(num(&j, "gradient_max"), expected, 1e-12);
    assert!(j.contains("\"converged\":false"));
}
#[test]
fn clipping_weights_must_not_fabricate_empirical_common_support() {
    // Mutation caught: evaluating overlap only after truncating distinct input scores to one value.
    let d = [
        vec![0., 1., 0.0001],
        vec![1., 1., 0.0002],
        vec![2., 0., 0.0003],
        vec![3., 0., 0.0004],
    ];
    let j = execute("causal_balance", &d, &[0., 0.01]).unwrap();
    assert!(j.contains("\"common_support\":null"));
    near(num(&j, "outside_common_support"), 4., 0.);
    near(num(&j, "n_clipped"), 4., 0.);
}
