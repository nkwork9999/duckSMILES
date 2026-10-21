use ducksmiles_analytics::resampling::run;
fn value(j: &str, k: &str) -> f64 {
    j.split_once(&format!("\"{k}\":"))
        .unwrap()
        .1
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
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
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8, "{a}!={b}");
}
#[test]
fn exact_permutation_enumerates_group_assignments_and_monte_carlo_is_seeded() {
    let d = vec![vec![1., 0.], vec![2., 0.], vec![3., 1.], vec![4., 1.]];
    let j = run("permutation_test", &d, &[0., 42., 0.]).unwrap();
    close(value(&j, "p_value"), 1. / 3.);
    close(value(&j, "assignments"), 6.);
    let a = run("permutation_test", &d, &[500., 42., 0.]).unwrap();
    let b = run("permutation_test", &d, &[500., 42., 0.]).unwrap();
    assert_eq!(a, b);
    assert!((value(&a, "p_value") - 1. / 3.).abs() < 0.08);
}
#[test]
fn bootstrap_mean_is_seeded_and_preserves_a_constant_and_location_shifts() {
    let j = run("bootstrap_mean", &[vec![3.], vec![3.]], &[100., 0.95, 1.]).unwrap();
    assert_eq!(array(&j, "confidence_interval"), vec![3., 3.]);
    let a = run(
        "bootstrap_mean",
        &[vec![1.], vec![2.], vec![4.]],
        &[1000., 0.95, 42.],
    )
    .unwrap();
    let b = run(
        "bootstrap_mean",
        &[vec![101.], vec![102.], vec![104.]],
        &[1000., 0.95, 42.],
    )
    .unwrap();
    close(
        array(&b, "confidence_interval")[0] - array(&a, "confidence_interval")[0],
        100.,
    );
}
#[test]
fn row_sampling_returns_only_input_indices_and_no_replacement_duplicates() {
    let d = vec![vec![1.], vec![2.], vec![3.]];
    let j = run("sample_rows", &d, &[3., 0., 42.]).unwrap();
    let mut ix = array(&j, "indices");
    ix.sort_by(f64::total_cmp);
    assert_eq!(ix, vec![0., 1., 2.]);
    let j = run("sample_rows", &d, &[10., 1., 42.]).unwrap();
    assert_eq!(array(&j, "indices").len(), 10);
}
#[test]
fn normal_simulation_is_reproducible_and_affine() {
    let a = run("simulate_normal", &[vec![0., 1.]], &[10., 42.]).unwrap();
    let b = run("simulate_normal", &[vec![3., 2.]], &[10., 42.]).unwrap();
    let x = array(&a, "samples");
    let y = array(&b, "samples");
    for (a, b) in x.iter().zip(y) {
        close(3. + 2. * a, b);
    }
}
#[test]
fn known_sigma_power_and_margin_sample_sizes_match_reference() {
    let j = run("power_z", &[vec![0., 1.]], &[25., 0.05, 0.]).unwrap();
    close(array(&j, "power")[0], 0.05);
    let j = run("power_z", &[vec![1., 1.]], &[25., 0.05, 0.]).unwrap();
    assert!(array(&j, "power")[0] > 0.998);
    let j = run("sample_size_mean", &[vec![1.]], &[0.2, 0.95]).unwrap();
    assert_eq!(array(&j, "sample_sizes"), vec![97.]);
    let j = run("sample_size_proportion", &[vec![0.5]], &[0.05, 0.95]).unwrap();
    assert_eq!(array(&j, "sample_sizes"), vec![385.]);
    let j = run("reliability_demonstration", &[vec![0.9]], &[0.95]).unwrap();
    assert_eq!(array(&j, "sample_sizes"), vec![29.]);
}
#[test]
fn contingency_effects_match_hand_computed_odds_and_risk() {
    let d = vec![vec![20., 80.], vec![10., 90.]];
    let j = run("odds_ratio", &d, &[0.95, 0.]).unwrap();
    close(value(&j, "estimate"), 2.25);
    let j = run("relative_risk", &d, &[0.95, 0.]).unwrap();
    close(value(&j, "estimate"), 2.);
    let j = run(
        "cohen_d",
        &[
            vec![1., 0.],
            vec![2., 0.],
            vec![3., 0.],
            vec![3., 1.],
            vec![4., 1.],
            vec![5., 1.],
        ],
        &[],
    )
    .unwrap();
    close(value(&j, "cohen_d"), -2.);
    assert!(value(&j, "hedges_g").abs() < 2.);
}
#[test]
fn sampling_rejects_invalid_plans_probabilities_and_resources() {
    assert!(run("sample_rows", &[vec![1.]], &[2., 0., 1.]).is_err());
    assert!(run("power_z", &[vec![1., 0.]], &[10., 0.05, 0.]).is_err());
    assert!(run("bootstrap_mean", &[vec![1.], vec![2.]], &[1e10, 0.95, 0.]).is_err());
    assert!(run("reliability_demonstration", &[vec![1.]], &[0.95]).is_err());
}

#[test]
fn permutation_tie_policy_is_invariant_to_small_measurement_units() {
    let a = run(
        "permutation_test",
        &[vec![0., 0.], vec![1., 0.], vec![2., 1.], vec![3., 1.]],
        &[0.],
    )
    .unwrap();
    let b = run(
        "permutation_test",
        &[
            vec![0., 0.],
            vec![1e-15, 0.],
            vec![2e-15, 1.],
            vec![3e-15, 1.],
        ],
        &[0.],
    )
    .unwrap();
    close(value(&a, "p_value"), value(&b, "p_value"));
}

#[test]
fn bootstrap_standard_error_preserves_sub_ulp_resample_mean_offsets() {
    let a = run(
        "bootstrap_mean",
        &[vec![0.], vec![2.]],
        &[10000., 0.95, 42.],
    )
    .unwrap();
    let b = run(
        "bootstrap_mean",
        &[vec![1e16], vec![1e16 + 2.]],
        &[10000., 0.95, 42.],
    )
    .unwrap();
    close(
        value(&a, "bootstrap_standard_error"),
        value(&b, "bootstrap_standard_error"),
    );
}

#[test]
fn sample_size_mean_is_invariant_to_tiny_units() {
    let a = run("sample_size_mean", &[vec![1.]], &[1., 0.95]).unwrap();
    let b = run("sample_size_mean", &[vec![1e-200]], &[1e-200, 0.95]).unwrap();
    assert_eq!(array(&a, "sample_sizes"), array(&b, "sample_sizes"));
}
