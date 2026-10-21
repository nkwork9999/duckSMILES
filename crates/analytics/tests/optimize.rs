use ducksmiles_analytics::optimize::run;

fn numbers(json: &str, key: &str) -> Vec<f64> {
    let tail = json.split_once(&format!("\"{key}\":")).unwrap().1;
    let array = tail.split_once('[').unwrap().1.split_once(']').unwrap().0;
    array
        .split(',')
        .filter(|s| !s.trim().is_empty())
        .map(|s| s.trim().parse().unwrap())
        .collect()
}

#[test]
fn expected_improvement_matches_normal_at_incumbent_and_zero_uncertainty() {
    // Dropping predictive uncertainty or treating sigma as variance breaks this fixture.
    let result = run(
        "acquisition",
        &[vec![2., 1.], vec![3., 0.], vec![1., 0.]],
        &[2., 0., 0.],
    )
    .unwrap();
    let ei = numbers(&result, "expected_improvement");
    assert!((ei[0] - 0.3989422804014327).abs() < 1e-8);
    assert_eq!(&ei[1..], &[1., 0.]);
    let result = run("acquisition", &[vec![1., 0.], vec![3., 0.]], &[2., 1., 0.]).unwrap();
    assert_eq!(numbers(&result, "expected_improvement"), vec![1., 0.]);
}

#[test]
fn gp_posterior_uses_noise_variance_and_returns_candidate_indices() {
    // One observed target centers the prior at 2. K=1, noise=1 -> latent var=1/2.
    let result = run(
        "gp_predict",
        &[vec![0., 2.], vec![0., 0.], vec![100., 0.]],
        &[1., 1., 1., 1., 0., 0.],
    )
    .unwrap();
    assert_eq!(numbers(&result, "predictions"), vec![2., 2.]);
    let variance = numbers(&result, "variances");
    assert!((variance[0] - 0.5).abs() < 1e-10);
    assert!((variance[1] - 1.).abs() < 1e-10);
    let ei = numbers(&result, "expected_improvement");
    assert!((ei[0] - 0.28209479177387814).abs() < 1e-8);
    assert_eq!(numbers(&result, "ranked_candidate_indices"), vec![1., 0.]);
}

#[test]
fn active_learning_excludes_existing_rows_and_redundant_candidates() {
    // Existing x=0: farthest first picks x=10, then x=5, never existing or duplicates.
    let result = run(
        "diversity_select",
        &[vec![0.], vec![1.], vec![5.], vec![10.], vec![10.]],
        &[2., 1.],
    )
    .unwrap();
    assert_eq!(numbers(&result, "selected_indices"), vec![3., 2.]);
}

#[test]
fn optimization_rejects_bad_domains_and_singular_noise_free_gp() {
    assert!(run("acquisition", &[vec![1., -1.]], &[0., 0., 0.]).is_err());
    assert!(run(
        "gp_predict",
        &[vec![0., 1.], vec![0., 2.], vec![1., 0.]],
        &[2., 1., 1., 0., 0., 0.]
    )
    .is_err());
    assert!(run(
        "gp_predict",
        &[vec![0., 1.], vec![1., 0.]],
        &[1., 0., 1., 0.1, 0., 0.]
    )
    .is_err());
    assert!(run("diversity_select", &[vec![0.], vec![0.]], &[2., 0.]).is_err());
    assert!(run("acquisition", &[vec![f64::NAN, 1.]], &[0., 0., 0.]).is_err());
}

#[test]
fn optimization_rejects_overflow_and_excessive_selection_work() {
    // Finite inputs must not serialize NaN/Infinity or initiate quadratic initial work.
    assert!(run("acquisition", &[vec![f64::MAX, 1.]], &[-f64::MAX, 0., 0.]).is_err());
    let many = vec![vec![0.; 100]; 10000];
    assert!(run("diversity_select", &many, &[1., 5000.]).is_err());
    let gp_rows = vec![vec![0., 0.]; 10000];
    assert!(run("gp_predict", &gp_rows, &[512., 1., 1., 0.1, 0., 0.]).is_err());
    assert!(run(
        "diversity_select",
        &[vec![0.], vec![f64::MAX], vec![-f64::MAX]],
        &[1., 1.]
    )
    .is_err());
}

#[test]
fn gp_constant_decimal_training_feature_has_unit_scale() {
    // Rounding a constant mean must not invent tiny feature variance.
    let mut data = vec![vec![0.1, 2.]; 8];
    data.push(vec![1., 0.]);
    let result = run("gp_predict", &data, &[8., 1., 1., 0.1, 0., 0.]).unwrap();
    assert_eq!(numbers(&result, "feature_scales"), vec![1.]);
}

#[test]
fn expected_improvement_resolves_small_gaussian_tail_values() {
    // Losing the lower-tail CDF to subtraction changes candidate rankings.
    // Independent scipy.stats.norm fixture: phi(-10)-10*Phi(-10).
    let result = run("acquisition", &[vec![-10., 1.]], &[0., 0., 0.]).unwrap();
    let actual = numbers(&result, "expected_improvement")[0];
    let expected = 7.474560254589328e-25;
    assert!(
        (actual / expected - 1.).abs() < 1e-10,
        "{actual} != {expected}"
    );
}
