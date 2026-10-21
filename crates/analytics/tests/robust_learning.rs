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
    let t = &j[p..p + j[p..].find(']').unwrap()];
    if t.is_empty() {
        vec![]
    } else {
        t.split(',').map(|x| x.parse().unwrap()).collect()
    }
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
#[test]
fn linear_svm_solves_soft_margin_dual_and_scores_saved_model() {
    // Mutation caught: ridge/logistic substitute, missing dual equality, or wrong KKT sign.
    let d = [vec![-2., 0.], vec![-1., 0.], vec![1., 1.], vec![2., 1.]];
    let j = execute("svm_fit", &d, &[4., 0., 10., 1., 1e-9, 10000.]).unwrap();
    assert!(j.contains("\"converged\":true"));
    close(num(&j, "dual_objective"), 1.25, 1e-8);
    close(num(&j, "primal_objective"), 1.25, 1e-7);
    close(num(&j, "equality_residual"), 0., 1e-14);
    assert!(num(&j, "kkt_gap") <= 1e-9);
    for (a, b) in vals(&j, "dual_alphas").iter().zip([0., 1.25, 1.25, 0.]) {
        close(*a, b, 1e-7);
    }
    let q = execute(
        "svm_predict",
        &[vec![-3.], vec![0.], vec![3.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    for (a, b) in vals(&q, "decision_values").iter().zip([-3., 0., 3.]) {
        close(*a, b, 1e-7);
    }
    assert_eq!(vals(&q, "predictions"), vec![0., 1., 1.]);
}
#[test]
fn rbf_svm_separates_xor_and_reports_exhausted_optimization() {
    // Mutation caught: a linear margin labeled as RBF SVM or false convergence on limit.
    let d = [
        vec![0., 0., 0.],
        vec![0., 1., 1.],
        vec![1., 0., 1.],
        vec![1., 1., 0.],
    ];
    let j = execute("svm_fit", &d, &[4., 1., 10., 0.5, 1e-8, 10000.]).unwrap();
    assert_eq!(vals(&j, "predictions"), vec![0., 1., 1., 0.]);
    assert!(num(&j, "kkt_gap") <= 1e-8);
    close(num(&j, "duality_gap"), 0., 1e-6);
    let short = execute("svm_fit", &d, &[4., 1., 10., 0.5, 1e-10, 1.]).unwrap();
    assert!(short.contains("\"converged\":false"));
    assert!(num(&short, "kkt_gap") > 1e-10);
}
#[test]
fn svm_standardization_uses_training_prefix_and_ignores_query_labels() {
    // Mutation caught: query-driven feature centering or using suffix targets during training.
    let a = [
        vec![-2., 0.],
        vec![-1., 0.],
        vec![1., 1.],
        vec![2., 1.],
        vec![50., -100.],
    ];
    let mut b = a.clone();
    b[4][1] = 1e200;
    let j = execute("svm_fit", &a, &[4., 0., 10.]).unwrap();
    let k = execute("svm_fit", &b, &[4., 0., 10.]).unwrap();
    assert_eq!(vals(&j, "dual_alphas"), vals(&k, "dual_alphas"));
    close(vals(&j, "decision_values")[4], 50., 1e-6);
}
#[test]
fn quantile_admm_minimizes_check_loss_and_saves_offset_safe_predictions() {
    // Mutation caught: least squares, absolute loss for every quantile, or incomplete convergence checks.
    let d = [vec![0.], vec![1.], vec![2.], vec![3.], vec![100.]];
    let j = execute("quantile_regression", &d, &[0.75, 1., 50000., 1e-9]).unwrap();
    close(num(&j, "intercept"), 3., 2e-6);
    close(num(&j, "objective"), 74.25, 1e-6);
    assert!(j.contains("\"converged\":true"));
    assert!(num(&j, "primal_residual") <= num(&j, "primal_tolerance"));
    assert!(num(&j, "dual_residual") <= num(&j, "dual_tolerance"));
    let line = [
        vec![0., 0.],
        vec![1., 1.],
        vec![2., 4.],
        vec![3., 3.],
        vec![4., 4.],
    ];
    let fit = execute("quantile_regression", &line, &[0.5, 1., 50000., 1e-9]).unwrap();
    close(vals(&fit, "coefficients")[0], 1., 2e-6);
    close(num(&fit, "objective"), 1., 1e-5);
    let q = execute(
        "robust_predict",
        &[vec![2.5]],
        &vals(&fit, "predict_parameters"),
    )
    .unwrap();
    close(vals(&q, "predictions")[0], 2.5, 1e-5);
    let short = execute("quantile_regression", &line, &[0.3, 1., 1., 1e-10]).unwrap();
    assert!(short.contains("\"converged\":false"));
}
#[test]
fn huber_irls_limits_outlier_score_and_reports_fixed_scale() {
    // Mutation caught: ordinary weighted residual squares or recomputing scale silently.
    let d = [
        vec![0., 1.],
        vec![1., 3.],
        vec![2., 5.],
        vec![3., 27.],
        vec![4., 9.],
        vec![5., 11.],
        vec![6., 13.],
    ];
    let j = execute("huber_regression", &d, &[1.345, 1., 500., 1e-10]).unwrap();
    close(num(&j, "intercept"), 1. + 1.345 / 6., 1e-8);
    close(vals(&j, "coefficients")[0], 2., 1e-8);
    close(num(&j, "robust_scale"), 1., 1e-14);
    assert!(vals(&j, "weights")[3] < 0.1);
    assert!(j.contains("\"converged\":true"));
    let q = execute(
        "robust_predict",
        &[vec![10.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    close(vals(&q, "predictions")[0], 21. + 1.345 / 6., 1e-8);
}
#[test]
fn robust_contracts_reject_invalid_domains_and_unidentified_regressions() {
    assert!(execute("svm_fit", &[vec![0., 0.], vec![1., 0.]], &[2.]).is_err());
    assert!(execute("svm_fit", &[vec![0., 0.], vec![1., 1.]], &[2., 0., 0.]).is_err());
    assert!(execute("quantile_regression", &[vec![1.], vec![2.]], &[1.]).is_err());
    assert!(execute(
        "huber_regression",
        &[vec![0., 1.], vec![0., 2.], vec![0., 3.]],
        &[1.345, 1.]
    )
    .is_err());
    assert!(execute("robust_predict", &[vec![1.]], &[1., 0., 1.]).is_err());
}

#[test]
fn svm_convergence_is_checked_after_the_final_allowed_pair_update() {
    // Mutation caught: reporting false when the last permitted feasible update met KKT.
    let j = execute(
        "svm_fit",
        &[vec![-1., 0.], vec![1., 1.]],
        &[2., 0., 1., 1., 1e-10, 1.],
    )
    .unwrap();
    assert!(num(&j, "kkt_gap") <= 1e-10);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn intercept_only_quantile_saved_model_scores_dummy_rows() {
    // Mutation caught: rejecting all portable prediction calls for valid intercept-only fits.
    let j = execute(
        "quantile_regression",
        &[vec![0.], vec![1.], vec![2.], vec![3.], vec![100.]],
        &[0.5, 1., 50000., 1e-9],
    )
    .unwrap();
    let q = execute(
        "robust_predict",
        &[vec![0.], vec![100.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    for y in vals(&q, "predictions") {
        close(y, 2., 2e-6);
    }
}

#[test]
fn huber_sandwich_uncertainty_survives_large_fixed_scale_units() {
    // Mutation caught: squaring tiny standardized scores before restoring their scale.
    let d = [vec![0., 1.], vec![1., 2.], vec![2., 1.], vec![3., 4.]];
    let normal = execute("huber_regression", &d, &[1.345, 100.]).unwrap();
    let huge = execute("huber_regression", &d, &[1.345, 1e200]).unwrap();
    assert!(vals(&normal, "coefficient_standard_errors")[0] > 0.1);
    close(
        vals(&huge, "coefficient_standard_errors")[0],
        vals(&normal, "coefficient_standard_errors")[0],
        1e-12,
    );
}
