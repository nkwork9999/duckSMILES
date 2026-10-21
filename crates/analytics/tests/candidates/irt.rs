#[path = "../../src/irt.rs"]
mod irt;
use irt::run as execute;
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
    let t = &j[start..];
    let end = t.find(']').unwrap();
    if end == 0 {
        vec![]
    } else {
        t[..end].split(',').map(|x| x.parse().unwrap()).collect()
    }
}
fn matrix(j: &str, key: &str) -> Vec<Vec<f64>> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let t = &j[start..];
    let end = t.find("]]").unwrap() + 1;
    t[..end]
        .split("],")
        .map(|s| {
            s.trim_matches(['[', ']'])
                .split(',')
                .map(|x| x.parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(x: f64, y: f64, t: f64) {
    assert!((x - y).abs() <= t * (1. + y.abs()), "{x} != {y}")
}
fn check(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y, tol)
    }
}

#[test]
fn three_parameter_curve_uses_guessing_and_fisher_information_not_rasch_formula() {
    let model = [1., 3., 1., 0., 61., 2., 0.5, 0.2];
    let j = execute(
        "irt_predict",
        &[vec![0.5], vec![-1e300], vec![1e300]],
        &model,
    )
    .unwrap();
    let p = matrix(&j, "probabilities");
    check(&p[0], &[0.6], 1e-13);
    check(&p[1], &[0.2], 1e-13);
    check(&p[2], &[1.], 1e-13);
    let info = matrix(&j, "item_information");
    check(&info[0], &[2. / 3.], 1e-13);
    check(&info[1], &[0.], 1e-13);
    check(&info[2], &[0.], 1e-13);
    check(&values(&j, "test_information"), &[2. / 3., 0., 0.], 1e-13);
}
#[test]
fn posterior_eap_rasch_extremes_are_finite_and_not_person_mle() {
    let model = [1., 1., 2., 0., 61., 1., 0., 0., 1., 0., 0.];
    let j = execute(
        "irt_score",
        &[vec![0., 0.], vec![0., 1.], vec![1., 1.]],
        &model,
    )
    .unwrap();
    check(
        &values(&j, "posterior_means"),
        &[-0.7042799208115517, 0., 0.7042799208115517],
        1e-9,
    );
    check(
        &values(&j, "posterior_variances"),
        &[0.7166719426236876, 0.6980147478335965, 0.7166719426236875],
        1e-9,
    );
    check(
        &values(&j, "marginal_log_probabilities"),
        &[
            -1.2262898682525063,
            -1.5768692553152777,
            -1.2262898682525063,
        ],
        1e-10,
    );
    assert!(j.contains("\"quadrature_converged\":true"), "{j}");
}
#[test]
fn omitted_missing_item_uses_prior_when_no_responses_and_is_not_an_incorrect_answer() {
    let mut model = [1., 1., 2., 2., 61., 1., 0., 0., 1., 0., 0.];
    let j = execute("irt_score", &[vec![-1., -1.], vec![1., -1.]], &model).unwrap();
    close(values(&j, "posterior_means")[0], 0., 1e-12);
    close(values(&j, "posterior_variances")[0], 1., 1e-12);
    close(values(&j, "marginal_log_probabilities")[0], 0., 1e-12);
    assert!(values(&j, "posterior_means")[1] > 0.);
    model[3] = 1.;
    let incorrect = execute("irt_score", &[vec![1., -1.]], &model).unwrap();
    close(values(&incorrect, "posterior_means")[0], 0., 1e-10);
    model[3] = 0.;
    assert!(execute("irt_score", &[vec![1., -1.]], &model).is_err());
}
#[test]
fn rasch_marginal_ml_frequency_matches_expanded_people_and_symmetric_known_likelihood() {
    let compressed = vec![
        vec![0., 0., 30.],
        vec![0., 1., 20.],
        vec![1., 0., 20.],
        vec![1., 1., 30.],
    ];
    let j = execute(
        "irt_fit",
        &compressed,
        &[1., 1., 0., 300., 1e-8, 61., 0.5, 0.],
    )
    .unwrap();
    check(&values(&j, "difficulties"), &[0., 0.], 1e-7);
    check(&values(&j, "discriminations"), &[1., 1.], 1e-12);
    check(&values(&j, "guessing"), &[0., 0.], 1e-12);
    close(
        number(&j, "log_likelihood"),
        60. * (-1.2262898682525063) + 40. * (-1.5768692553152777),
        1e-10,
    );
    let expanded = compressed
        .iter()
        .flat_map(|r| vec![r[..2].to_vec(); r[2] as usize])
        .collect::<Vec<_>>();
    let other = execute(
        "irt_fit",
        &expanded,
        &[1., 0., 0., 300., 1e-8, 61., 0.5, 0.],
    )
    .unwrap();
    close(
        number(&j, "log_likelihood"),
        number(&other, "log_likelihood"),
        1e-10,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
}
#[test]
fn irt_rejects_unidentified_items_bad_frequency_and_invalid_saved_models() {
    assert!(execute(
        "irt_fit",
        &[vec![0., 0.], vec![0., 1.], vec![0., 1.]],
        &[1.]
    )
    .is_err());
    assert!(execute("irt_fit", &[vec![0., 1.], vec![1., 0.]], &[2.]).is_err());
    assert!(execute("irt_fit", &[vec![0., 1., 0.5], vec![1., 0., 1.]], &[1., 1.]).is_err());
    for op in ["irt_fit", "irt_score", "irt_predict"] {
        assert!(execute(op, &[], &[]).is_err());
        assert!(execute(op, &[vec![f64::NAN]], &[]).is_err());
    }
    for model in [
        vec![1., 1., 1., 0., 31., 2., 0., 0.],
        vec![1., 2., 1., 0., 31., 1., 0., 0.1],
        vec![1., 3., 1., 0., 31., 0., 0., 0.2],
    ] {
        assert!(execute("irt_predict", &[vec![0.]], &model).is_err());
    }
}

#[test]
fn two_parameter_marginal_ml_matches_independent_real_line_integrated_fixture() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 0.0, 4503.0],
        vec![0.0, 0.0, 0.0, 1.0, 710.0],
        vec![0.0, 0.0, 1.0, 0.0, 847.0],
        vec![0.0, 0.0, 1.0, 1.0, 289.0],
        vec![0.0, 1.0, 0.0, 0.0, 2210.0],
        vec![0.0, 1.0, 0.0, 1.0, 651.0],
        vec![0.0, 1.0, 1.0, 0.0, 1044.0],
        vec![0.0, 1.0, 1.0, 1.0, 638.0],
        vec![1.0, 0.0, 0.0, 0.0, 4497.0],
        vec![1.0, 0.0, 0.0, 1.0, 1084.0],
        vec![1.0, 0.0, 1.0, 0.0, 1577.0],
        vec![1.0, 0.0, 1.0, 1.0, 794.0],
        vec![1.0, 1.0, 0.0, 0.0, 3654.0],
        vec![1.0, 1.0, 0.0, 1.0, 1591.0],
        vec![1.0, 1.0, 1.0, 0.0, 3096.0],
        vec![1.0, 1.0, 1.0, 1.0, 2815.0],
    ];
    let j = execute(
        "irt_fit",
        &rows,
        &[2.0, 1.0, 0.0, 1500.0, 1e-08, 61.0, 0.5, 0.0],
    )
    .unwrap();
    check(
        &values(&j, "difficulties"),
        &[
            -0.7999631379780021,
            -0.10002055708103759,
            0.49997424274335656,
            1.1000675845126098,
        ],
        3e-4,
    );
    check(
        &values(&j, "discriminations"),
        &[
            0.8001094009486159,
            1.1998611666660994,
            1.500549059633991,
            0.9999603624279422,
        ],
        3e-4,
    );
    check(&values(&j, "guessing"), &[0.0, 0.0, 0.0, 0.0], 3e-4);
    close(number(&j, "log_likelihood"), -75501.86267678159, 3e-9);
    close(
        number(&j, "penalized_log_likelihood"),
        -75501.86267678159,
        3e-9,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
}

#[test]
fn three_parameter_marginal_ml_estimates_identified_positive_guessing() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 9197.0],
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2848.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 3187.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1122.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 5419.0],
        vec![0.0, 0.0, 1.0, 0.0, 1.0, 2018.0],
        vec![0.0, 0.0, 1.0, 1.0, 0.0, 2641.0],
        vec![0.0, 0.0, 1.0, 1.0, 1.0, 1271.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 7109.0],
        vec![0.0, 1.0, 0.0, 0.0, 1.0, 2616.0],
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 3434.0],
        vec![0.0, 1.0, 0.0, 1.0, 1.0, 1581.0],
        vec![0.0, 1.0, 1.0, 0.0, 0.0, 6531.0],
        vec![0.0, 1.0, 1.0, 0.0, 1.0, 3364.0],
        vec![0.0, 1.0, 1.0, 1.0, 0.0, 4872.0],
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 3703.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 10776.0],
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 3801.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 4866.0],
        vec![1.0, 0.0, 0.0, 1.0, 1.0, 2086.0],
        vec![1.0, 0.0, 1.0, 0.0, 0.0, 8993.0],
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 4245.0],
        vec![1.0, 0.0, 1.0, 1.0, 0.0, 6081.0],
        vec![1.0, 0.0, 1.0, 1.0, 1.0, 4071.0],
        vec![1.0, 1.0, 0.0, 0.0, 0.0, 11747.0],
        vec![1.0, 1.0, 0.0, 0.0, 1.0, 5322.0],
        vec![1.0, 1.0, 0.0, 1.0, 0.0, 7616.0],
        vec![1.0, 1.0, 0.0, 1.0, 1.0, 4692.0],
        vec![1.0, 1.0, 1.0, 0.0, 0.0, 16228.0],
        vec![1.0, 1.0, 1.0, 0.0, 1.0, 11705.0],
        vec![1.0, 1.0, 1.0, 1.0, 0.0, 17004.0],
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 19856.0],
    ];
    let j = execute(
        "irt_fit",
        &rows,
        &[3.0, 1.0, 0.0, 1500.0, 1e-08, 61.0, 0.5, 0.0],
    )
    .unwrap();
    check(
        &values(&j, "difficulties"),
        &[
            -1.0008391722310146,
            -0.4000597654337309,
            0.00015946045050551974,
            0.6005339225181773,
            1.2000182178095813,
        ],
        3e-4,
    );
    check(
        &values(&j, "discriminations"),
        &[
            0.7998614011534033,
            1.0999361687909008,
            1.4002333597585983,
            1.0004007521524256,
            1.4997132319984212,
        ],
        3e-4,
    );
    check(
        &values(&j, "guessing"),
        &[
            0.07969101410691895,
            0.11996580474534009,
            0.18007496152839383,
            0.10020089987134288,
            0.19998464460571885,
        ],
        3e-4,
    );
    close(number(&j, "log_likelihood"), -644102.7326123144, 3e-9);
    close(
        number(&j, "penalized_log_likelihood"),
        -644102.7326123144,
        3e-9,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
}

#[test]
fn three_parameter_guessing_penalty_has_explicit_total_likelihood_units() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 9197.0],
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2848.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 3187.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1122.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 5419.0],
        vec![0.0, 0.0, 1.0, 0.0, 1.0, 2018.0],
        vec![0.0, 0.0, 1.0, 1.0, 0.0, 2641.0],
        vec![0.0, 0.0, 1.0, 1.0, 1.0, 1271.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 7109.0],
        vec![0.0, 1.0, 0.0, 0.0, 1.0, 2616.0],
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 3434.0],
        vec![0.0, 1.0, 0.0, 1.0, 1.0, 1581.0],
        vec![0.0, 1.0, 1.0, 0.0, 0.0, 6531.0],
        vec![0.0, 1.0, 1.0, 0.0, 1.0, 3364.0],
        vec![0.0, 1.0, 1.0, 1.0, 0.0, 4872.0],
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 3703.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 10776.0],
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 3801.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 4866.0],
        vec![1.0, 0.0, 0.0, 1.0, 1.0, 2086.0],
        vec![1.0, 0.0, 1.0, 0.0, 0.0, 8993.0],
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 4245.0],
        vec![1.0, 0.0, 1.0, 1.0, 0.0, 6081.0],
        vec![1.0, 0.0, 1.0, 1.0, 1.0, 4071.0],
        vec![1.0, 1.0, 0.0, 0.0, 0.0, 11747.0],
        vec![1.0, 1.0, 0.0, 0.0, 1.0, 5322.0],
        vec![1.0, 1.0, 0.0, 1.0, 0.0, 7616.0],
        vec![1.0, 1.0, 0.0, 1.0, 1.0, 4692.0],
        vec![1.0, 1.0, 1.0, 0.0, 0.0, 16228.0],
        vec![1.0, 1.0, 1.0, 0.0, 1.0, 11705.0],
        vec![1.0, 1.0, 1.0, 1.0, 0.0, 17004.0],
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 19856.0],
    ];
    let j = execute(
        "irt_fit",
        &rows,
        &[3.0, 1.0, 0.0, 1500.0, 1e-08, 61.0, 0.5, 100.0],
    )
    .unwrap();
    check(
        &values(&j, "difficulties"),
        &[
            -1.0866991289701133,
            -0.4289813268121235,
            -0.010022806381603369,
            0.6003882713152304,
            1.198278377600275,
        ],
        3e-4,
    );
    check(
        &values(&j, "discriminations"),
        &[
            0.7835724142283331,
            1.0852176295665414,
            1.3887043837413955,
            1.0007754812431429,
            1.5026457472174735,
        ],
        3e-4,
    );
    check(
        &values(&j, "guessing"),
        &[
            0.049077162519284345,
            0.10834540488731373,
            0.17579420270775586,
            0.10020112268446481,
            0.19987294604175002,
        ],
        3e-4,
    );
    close(number(&j, "log_likelihood"), -644102.9138489929, 3e-9);
    close(
        number(&j, "penalized_log_likelihood"),
        -644107.6658667875,
        3e-9,
    );
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(j.contains("\"information_identified\":true"), "{j}");
}

#[test]
fn separate_missing_blocks_cannot_claim_identified_two_parameter_items() {
    let rows = vec![
        vec![0., 0., -1., -1., 30.],
        vec![0., 1., -1., -1., 20.],
        vec![1., 0., -1., -1., 20.],
        vec![1., 1., -1., -1., 30.],
        vec![-1., -1., 0., 0., 30.],
        vec![-1., -1., 0., 1., 20.],
        vec![-1., -1., 1., 0., 20.],
        vec![-1., -1., 1., 1., 30.],
    ];
    let j = execute("irt_fit", &rows, &[2., 1., 2., 500., 1e-8, 61., 0.5, 0.]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"information_identified\":false"), "{j}");
}

#[test]
fn deterministic_perfect_item_agreement_hits_a_guard_instead_of_finite_discrimination_success() {
    let rows = vec![vec![0., 0., 0., 100.], vec![1., 1., 1., 100.]];
    let j = execute("irt_fit", &rows, &[2., 1., 0., 300., 1e-6, 41., 0.5, 0.]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
}

#[test]
fn finite_iteration_and_coarse_integration_limits_cannot_claim_irt_success() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 9197.0],
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 2848.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 3187.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 1122.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 5419.0],
        vec![0.0, 0.0, 1.0, 0.0, 1.0, 2018.0],
        vec![0.0, 0.0, 1.0, 1.0, 0.0, 2641.0],
        vec![0.0, 0.0, 1.0, 1.0, 1.0, 1271.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 7109.0],
        vec![0.0, 1.0, 0.0, 0.0, 1.0, 2616.0],
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 3434.0],
        vec![0.0, 1.0, 0.0, 1.0, 1.0, 1581.0],
        vec![0.0, 1.0, 1.0, 0.0, 0.0, 6531.0],
        vec![0.0, 1.0, 1.0, 0.0, 1.0, 3364.0],
        vec![0.0, 1.0, 1.0, 1.0, 0.0, 4872.0],
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 3703.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 10776.0],
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 3801.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 4866.0],
        vec![1.0, 0.0, 0.0, 1.0, 1.0, 2086.0],
        vec![1.0, 0.0, 1.0, 0.0, 0.0, 8993.0],
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 4245.0],
        vec![1.0, 0.0, 1.0, 1.0, 0.0, 6081.0],
        vec![1.0, 0.0, 1.0, 1.0, 1.0, 4071.0],
        vec![1.0, 1.0, 0.0, 0.0, 0.0, 11747.0],
        vec![1.0, 1.0, 0.0, 0.0, 1.0, 5322.0],
        vec![1.0, 1.0, 0.0, 1.0, 0.0, 7616.0],
        vec![1.0, 1.0, 0.0, 1.0, 1.0, 4692.0],
        vec![1.0, 1.0, 1.0, 0.0, 0.0, 16228.0],
        vec![1.0, 1.0, 1.0, 0.0, 1.0, 11705.0],
        vec![1.0, 1.0, 1.0, 1.0, 0.0, 17004.0],
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 19856.0],
    ];
    let short = execute("irt_fit", &rows, &[3., 1., 0., 1., 1e-8, 61., 0.5, 0.]).unwrap();
    assert!(short.contains("\"converged\":false"), "{short}");
    let coarse = execute("irt_fit", &rows, &[3., 1., 0., 1500., 1e-8, 11., 0.5, 0.]).unwrap();
    assert!(
        coarse.contains("\"quadrature_converged\":false"),
        "{coarse}"
    );
    assert!(coarse.contains("\"converged\":false"), "{coarse}");
}

#[test]
fn extremely_sharp_tail_ability_scores_report_insufficient_quadrature() {
    let mut model = vec![1., 2., 12., 0., 81.];
    for _ in 0..12 {
        model.extend([5., 8., 0.]);
    }
    let j = execute("irt_score", &[vec![1.; 12]], &model).unwrap();
    assert!(j.contains("\"quadrature_converged\":false"), "{j}");
    assert!(values(&j, "quadrature_order_deltas")[0] > 0.01);
    assert!(values(&j, "posterior_means")[0].is_finite());
}

#[test]
fn active_zero_guessing_and_small_positive_free_guessing_have_conditional_information() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 0.0, 0.0, 1924.0],
        vec![0.0, 0.0, 0.0, 0.0, 1.0, 92.0],
        vec![0.0, 0.0, 0.0, 1.0, 0.0, 408.0],
        vec![0.0, 0.0, 0.0, 1.0, 1.0, 40.0],
        vec![0.0, 0.0, 1.0, 0.0, 0.0, 583.0],
        vec![0.0, 0.0, 1.0, 0.0, 1.0, 74.0],
        vec![0.0, 0.0, 1.0, 1.0, 0.0, 241.0],
        vec![0.0, 0.0, 1.0, 1.0, 1.0, 59.0],
        vec![0.0, 1.0, 0.0, 0.0, 0.0, 1078.0],
        vec![0.0, 1.0, 0.0, 0.0, 1.0, 113.0],
        vec![0.0, 1.0, 0.0, 1.0, 0.0, 389.0],
        vec![0.0, 1.0, 0.0, 1.0, 1.0, 79.0],
        vec![0.0, 1.0, 1.0, 0.0, 0.0, 680.0],
        vec![0.0, 1.0, 1.0, 0.0, 1.0, 178.0],
        vec![0.0, 1.0, 1.0, 1.0, 0.0, 456.0],
        vec![0.0, 1.0, 1.0, 1.0, 1.0, 227.0],
        vec![1.0, 0.0, 0.0, 0.0, 0.0, 1920.0],
        vec![1.0, 0.0, 0.0, 0.0, 1.0, 163.0],
        vec![1.0, 0.0, 0.0, 1.0, 0.0, 603.0],
        vec![1.0, 0.0, 0.0, 1.0, 1.0, 100.0],
        vec![1.0, 0.0, 1.0, 0.0, 0.0, 999.0],
        vec![1.0, 0.0, 1.0, 0.0, 1.0, 215.0],
        vec![1.0, 0.0, 1.0, 1.0, 0.0, 588.0],
        vec![1.0, 0.0, 1.0, 1.0, 1.0, 241.0],
        vec![1.0, 1.0, 0.0, 0.0, 0.0, 1655.0],
        vec![1.0, 1.0, 0.0, 0.0, 1.0, 294.0],
        vec![1.0, 1.0, 0.0, 1.0, 0.0, 855.0],
        vec![1.0, 1.0, 0.0, 1.0, 1.0, 289.0],
        vec![1.0, 1.0, 1.0, 0.0, 0.0, 1716.0],
        vec![1.0, 1.0, 1.0, 0.0, 1.0, 751.0],
        vec![1.0, 1.0, 1.0, 1.0, 0.0, 1618.0],
        vec![1.0, 1.0, 1.0, 1.0, 1.0, 1372.0],
    ];
    let j = execute(
        "irt_fit",
        &rows,
        &[3., 1., 0., 2000., 1e-8, 61., 0.5, 1000.],
    )
    .unwrap();
    close(values(&j, "guessing")[0], 0., 1e-12);
    assert!(j.contains("\"information_identified\":true"), "{j}");
    assert!(j.contains("\"converged\":true"), "{j}");
    let start = j.find("\"parameter_standard_errors\":[").unwrap()
        + "\"parameter_standard_errors\":[".len();
    let errors = j[start..]
        .split(']')
        .next()
        .unwrap()
        .split(',')
        .collect::<Vec<_>>();
    assert_eq!(errors[10], "null");
    assert!(errors[11].parse::<f64>().unwrap() > 0.);
}

#[test]
fn unfinished_rasch_iterate_does_not_publish_mle_standard_errors_or_information_criteria() {
    let rows = vec![
        vec![0.0, 0.0, 0.0, 335.0],
        vec![0.0, 0.0, 1.0, 67.0],
        vec![0.0, 1.0, 0.0, 149.0],
        vec![0.0, 1.0, 1.0, 56.0],
        vec![1.0, 0.0, 0.0, 494.0],
        vec![1.0, 0.0, 1.0, 187.0],
        vec![1.0, 1.0, 0.0, 416.0],
        vec![1.0, 1.0, 1.0, 296.0],
    ];
    let j = execute("irt_fit", &rows, &[1., 1., 0., 1., 1e-8, 61., 0.5, 0.]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"parameter_standard_errors\":null"), "{j}");
    assert!(j.contains("\"parameter_covariance\":null"), "{j}");
    assert!(j.contains("\"aic\":null"), "{j}");
    assert!(j.contains("\"bic\":null"), "{j}");
    assert!(j.contains("\"inference_valid\":false"), "{j}");
}
