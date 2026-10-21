use ducksmiles_analytics::learning::run;
fn array(j: &str, k: &str) -> Vec<f64> {
    j.split_once(&format!("\"{k}\":["))
        .unwrap()
        .1
        .split(']')
        .next()
        .unwrap()
        .split(',')
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().unwrap())
        .collect()
}
#[test]
fn decision_tree_uses_training_only_splits_for_regression_and_classification() {
    let j = run(
        "decision_tree",
        &[
            vec![-2., -4.],
            vec![-1., -2.],
            vec![1., 2.],
            vec![2., 4.],
            vec![-1.5, 0.],
            vec![1.5, 0.],
        ],
        &[4., 1., 1., 0.],
    )
    .unwrap();
    assert_eq!(array(&j, "predictions"), vec![-3., 3.]);
    let d = vec![
        vec![0., 0., 0.],
        vec![0., 1., 0.],
        vec![1., 0., 0.],
        vec![1., 1., 1.],
        vec![1., 1., 999.],
        vec![0., 1., 999.],
    ];
    let j = run("decision_tree", &d, &[4., 3., 1., 1.]).unwrap();
    assert_eq!(array(&j, "predictions"), vec![1., 0.]);
}
#[test]
fn random_forest_is_seeded_and_never_uses_query_labels() {
    let mut d = vec![
        vec![0., 0.],
        vec![1., 0.],
        vec![2., 1.],
        vec![3., 1.],
        vec![-10., 999.],
        vec![10., 999.],
    ];
    let a = run("random_forest", &d, &[4., 50., 3., 1., 1., 42.]).unwrap();
    d[4][1] = -1.;
    d[5][1] = -2.;
    let b = run("random_forest", &d, &[4., 50., 3., 1., 1., 42.]).unwrap();
    assert_eq!(a, b);
    assert_eq!(array(&a, "predictions"), vec![0., 1.]);
}
#[test]
fn gradient_boosting_reduces_squared_and_binomial_loss() {
    let d = vec![
        vec![0., 0.],
        vec![1., 1.],
        vec![2., 2.],
        vec![3., 3.],
        vec![0., 0.],
        vec![3., 0.],
    ];
    let j = run("gradient_boosting", &d, &[4., 60., 0.1, 2., 1., 0., 42.]).unwrap();
    let p = array(&j, "predictions");
    assert!(p[0].abs() < 0.05 && (p[1] - 3.).abs() < 0.05);
    let d = vec![
        vec![-2., 0.],
        vec![-1., 0.],
        vec![1., 1.],
        vec![2., 1.],
        vec![-1.5, 999.],
        vec![1.5, 999.],
    ];
    let j = run("gradient_boosting", &d, &[4., 100., 0.2, 1., 1., 1., 42.]).unwrap();
    let p = array(&j, "probabilities");
    assert!(p[0] < 0.2 && p[1] > 0.8);
}
#[test]
fn cpu_mlp_learns_nonlinear_xor_with_training_only_standardization() {
    let d = vec![
        vec![0., 0., 0.],
        vec![0., 1., 1.],
        vec![1., 0., 1.],
        vec![1., 1., 0.],
        vec![0., 0., 999.],
        vec![0., 1., 999.],
        vec![1., 0., 999.],
        vec![1., 1., 999.],
    ];
    let j = run("mlp", &d, &[4., 8., 2500., 0.03, 1., 42.]).unwrap();
    let p = array(&j, "probabilities");
    assert!(
        p[0] < 0.2 && p[1] > 0.8 && p[2] > 0.8 && p[3] < 0.2,
        "{p:?}"
    );
}
#[test]
fn learners_reject_unbounded_work_invalid_labels_and_empty_training() {
    assert!(run("decision_tree", &[vec![0., 0.]], &[0., 2., 1., 0.]).is_err());
    assert!(run(
        "random_forest",
        &[vec![0., 0.], vec![1., 0.]],
        &[1., 1e9, 2., 1., 0., 0.]
    )
    .is_err());
    assert!(run(
        "mlp",
        &[vec![0., 2.], vec![1., 3.]],
        &[1., 2., 10., 0.1, 1., 0.]
    )
    .is_err());
}

#[test]
fn neural_training_is_invariant_to_representable_tiny_feature_units() {
    let mut d = vec![
        vec![0., 0., 0.],
        vec![0., 1., 1.],
        vec![1., 0., 1.],
        vec![1., 1., 0.],
        vec![0., 0., 0.],
        vec![0., 1., 0.],
        vec![1., 0., 0.],
        vec![1., 1., 0.],
    ];
    let a = run("mlp", &d, &[4., 8., 1000., 0.03, 1., 42.]).unwrap();
    for r in &mut d {
        r[0] *= 1e-200;
        r[1] *= 1e-200;
    }
    let b = run("mlp", &d, &[4., 8., 1000., 0.03, 1., 42.]).unwrap();
    for (a, b) in array(&a, "probabilities")
        .iter()
        .zip(array(&b, "probabilities"))
    {
        assert!((a - b).abs() < 1e-9, "{a} != {b}");
    }
}

#[test]
fn regression_trees_and_boosting_preserve_tiny_target_units() {
    let d = vec![vec![0., 0.], vec![1., 1e-200], vec![0., 0.], vec![1., 0.]];
    let j = run("decision_tree", &d, &[2., 2., 1., 0.]).unwrap();
    let p = array(&j, "predictions");
    assert_eq!(p[0], 0.);
    assert!((p[1] / 1e-200 - 1.).abs() < 1e-12, "{p:?}");
    let j = run("gradient_boosting", &d, &[2., 50., 0.2, 1., 1., 0., 0.]).unwrap();
    let p = array(&j, "predictions");
    assert!(
        p[0].abs() < 1e-204 && (p[1] / 1e-200 - 1.).abs() < 1e-4,
        "{p:?}"
    );
}

#[test]
fn neural_prediction_restores_large_target_origin_after_centered_adjustment() {
    let mut d: Vec<_> = (0..4)
        .map(|i| vec![i as f64, 1e16 + 2. * i as f64])
        .collect();
    d.extend((0..4).map(|i| vec![i as f64, 0.]));
    let j = run("mlp", &d, &[4., 8., 4000., 0.02, 0., 42., 1e-10]).unwrap();
    assert_eq!(
        array(&j, "predictions"),
        vec![1e16, 1e16 + 2., 1e16 + 4., 1e16 + 6.]
    );
}

#[test]
fn returned_neural_training_loss_matches_returned_model_predictions() {
    let mut d: Vec<_> = (0..4).map(|i| vec![i as f64, 2. * i as f64]).collect();
    d.extend((0..4).map(|i| vec![i as f64, 0.]));
    let j = run("mlp", &d, &[4., 3., 1., 0.1, 0., 42., 1e-10]).unwrap();
    let scalar = |key: &str| -> f64 {
        j.split_once(&format!("\"{key}\":"))
            .unwrap()
            .1
            .split([',', '}', ']'])
            .next()
            .unwrap()
            .parse()
            .unwrap()
    };
    let scale = scalar("target_scale");
    let actual = array(&j, "predictions")
        .iter()
        .enumerate()
        .map(|(i, p)| 0.5 * ((p - 2. * i as f64) / scale).powi(2) / 4.)
        .sum::<f64>();
    assert!((scalar("training_loss") - actual).abs() < 1e-10);
}
