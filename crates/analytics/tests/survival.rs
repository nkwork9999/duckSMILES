use ducksmiles_analytics::survival::run;

fn number(j: &str, key: &str) -> f64 {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let t = &j[start..];
    t[..t.find([',', '}', ']']).unwrap_or(t.len())]
        .parse()
        .unwrap()
}
fn array(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let t = &j[start..];
    let s = &t[..t.find(']').unwrap()];
    if s.is_empty() {
        vec![]
    } else {
        s.split(',').map(|v| v.parse().unwrap()).collect()
    }
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-8 * b.abs().max(1.), "{a} != {b}");
}
fn close_array(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        close(x, y);
    }
}
fn data() -> Vec<Vec<f64>> {
    let features = [
        [0., 1.],
        [1., 0.],
        [0., 0.],
        [2., 1.],
        [1., 2.],
        [2., 0.],
        [0., 2.],
        [1., 1.],
        [2., 2.],
        [0., 1.],
        [2., 0.],
        [1., 2.],
    ];
    let times = [1., 1., 2., 3., 3., 4., 4., 5., 6., 6., 7., 8.];
    let events = [1., 1., 0., 1., 0., 1., 0., 1., 1., 0., 1., 0.];
    (0..12)
        .map(|i| vec![features[i][0], features[i][1], times[i], events[i]])
        .collect()
}

#[test]
fn breslow_cox_matches_independent_scipy_score_root_and_observed_information() {
    // SciPy root of explicit per-event risk sets, Hessian as weighted covariance.
    // Censors at t=3 and t=6 must remain at risk for the deaths at those times.
    let j = run("cox_ph", &data(), &[0.]).unwrap();
    close_array(
        &array(&j, "coefficients"),
        &[0.22017852171500554, -0.5010300380019524],
    );
    close_array(
        &array(&j, "standard_errors"),
        &[0.5612256130912764, 0.5070171451512256],
    );
    close_array(
        &array(&j, "p_values"),
        &[0.6948237766871834, 0.3230588655938329],
    );
    close(number(&j, "partial_log_likelihood"), -11.911499078252731);
    close(
        number(&j, "null_partial_log_likelihood"),
        -12.801827480081469,
    );
    close(number(&j, "likelihood_ratio_statistic"), 1.7806568036574753);
    close(number(&j, "likelihood_ratio_p_value"), 0.41052091479388786);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn efron_cox_corrects_each_tied_death_denominator_and_defaults_to_efron() {
    let j = run("cox_ph", &data(), &[]).unwrap();
    close_array(
        &array(&j, "coefficients"),
        &[0.19974645776813277, -0.5155955813181949],
    );
    close_array(
        &array(&j, "standard_errors"),
        &[0.5579609796004918, 0.5080193998409992],
    );
    close_array(
        &array(&j, "p_values"),
        &[0.7203481057009973, 0.310147218009456],
    );
    close(number(&j, "partial_log_likelihood"), -11.814674575012589);
    close(
        number(&j, "null_partial_log_likelihood"),
        -12.71481610309184,
    );
    assert!(j.contains("\"ties\":\"efron\""));
}
#[test]
fn baseline_breslow_hazard_uses_centered_linear_predictors_at_every_death_time() {
    let j = run("cox_ph", &data(), &[0.]).unwrap();
    assert_eq!(array(&j, "baseline_times"), vec![1., 3., 4., 5., 6., 7.]);
    close_array(
        &array(&j, "baseline_cumulative_hazard"),
        &[
            0.14933530706921397,
            0.25333268292492384,
            0.38214186335934397,
            0.5737002732704883,
            0.8106480160135234,
            1.1861891921740948,
        ],
    );
    close_array(
        &array(&j, "baseline_survival"),
        &[
            0.8612802731121689,
            0.776209607181867,
            0.682398239031761,
            0.5634367068574846,
            0.44456988445556817,
            0.30538280463888534,
        ],
    );
    assert!(j.contains("\"baseline_hazard_estimator\":\"breslow\""));
    close_array(&array(&j, "feature_means"), &[1., 1.]);
}
#[test]
fn ridge_cox_penalizes_standardized_slopes_and_does_not_report_calibrated_wald_p_values() {
    let j = run("cox_ph", &data(), &[1., 0.5]).unwrap();
    close_array(
        &array(&j, "coefficients"),
        &[0.1972643906593223, -0.4792323794651989],
    );
    close_array(
        &array(&j, "standard_errors"),
        &[0.523053525700643, 0.47736478721385106],
    );
    assert!(j.contains("\"p_values\":[null,null]"));
    assert!(j.contains("\"likelihood_ratio_p_value\":null"));
    assert!(j.contains("penalized_curvature_no_calibrated_inference"));
}
#[test]
fn cox_prediction_reconstructs_training_reference_and_known_hazard_ratios() {
    let b = 2_f64.ln();
    let j = run(
        "cox_predict",
        &[vec![0.], vec![1.], vec![2.]],
        &[0.8, b, 0., 0.],
    )
    .unwrap();
    close_array(&array(&j, "relative_risks"), &[1., 2., 4.]);
    close_array(&array(&j, "survival_probabilities"), &[0.8, 0.64, 0.4096]);
    let j = run(
        "relative_risk_predict",
        &[vec![1e16 + 2.], vec![1e16 + 4.]],
        &[b, 1e16, 2.],
    )
    .unwrap();
    close_array(&array(&j, "relative_risks"), &[1., 4.]);
    for base in [0., 1.] {
        let j = run("cox_predict", &[vec![2.]], &[base, b, 0., 0.]).unwrap();
        close_array(&array(&j, "survival_probabilities"), &[base]);
    }
    let fitted = run("cox_ph", &data(), &[]).unwrap();
    let mut args = vec![array(&fitted, "baseline_survival")[2]];
    for key in ["coefficients", "feature_origins", "mean_offsets"] {
        args.extend(array(&fitted, key));
    }
    let j = run("cox_predict", &[vec![1., 1.]], &args).unwrap();
    close(array(&j, "survival_probabilities")[0], args[0]);
}
#[test]
fn tied_row_permutations_time_units_and_feature_units_preserve_the_fitted_model() {
    let baseline = run("cox_ph", &data(), &[]).unwrap();
    let coefficients = array(&baseline, "coefficients");
    let survival = array(&baseline, "baseline_survival");
    let mut seed = 0xc0a_u64;
    for _ in 0..24 {
        let mut rows = data();
        for i in (1..rows.len()).rev() {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            rows.swap(i, seed as usize % (i + 1));
        }
        for r in &mut rows {
            r[2] *= 60.;
            r[0] = 100. + 3. * r[0];
            r[1] = -20. + 2. * r[1];
        }
        let j = run("cox_ph", &rows, &[]).unwrap();
        close_array(
            &array(&j, "coefficients"),
            &[coefficients[0] / 3., coefficients[1] / 2.],
        );
        close_array(&array(&j, "baseline_survival"), &survival);
        close(
            number(&j, "partial_log_likelihood"),
            number(&baseline, "partial_log_likelihood"),
        );
    }
}
#[test]
fn co_infinite_coefficients_singular_information_and_separation_are_errors() {
    let mut duplicate = data();
    for r in &mut duplicate {
        r[1] = r[0];
    }
    assert!(run("cox_ph", &duplicate, &[]).is_err());
    let separated = (1..12)
        .map(|i| vec![i as f64, i as f64, 1.])
        .collect::<Vec<_>>();
    assert!(run("cox_ph", &separated, &[]).is_err());
    assert!(run("cox_ph", &data(), &[1., 0., 1., 1e-12]).is_err());
    let mut noevents = data();
    for r in &mut noevents {
        r[3] = 0.;
    }
    assert!(run("cox_ph", &noevents, &[]).is_err());
}
#[test]
fn invalid_cox_data_parameters_and_prediction_overflow_fail_without_partial_results() {
    for rows in [
        vec![],
        vec![vec![1., 1.]],
        vec![vec![f64::NAN, 1., 1.]],
        vec![vec![1., -1., 1.]],
        vec![vec![1., 1., 0.5]],
        vec![vec![1., 1., 1.], vec![2., 2.]],
        vec![vec![1., 1., 1.]; 100_001],
    ] {
        assert!(run("cox_ph", &rows, &[]).is_err());
    }
    for args in [
        vec![2.],
        vec![0., -1.],
        vec![0., 0., 0.],
        vec![0., 0., 1001.],
        vec![0., 0., 100., 0.],
        vec![f64::NAN],
    ] {
        assert!(run("cox_ph", &data(), &args).is_err());
    }
    assert!(run("cox_predict", &[vec![1.]], &[1.1, 1., 0., 0.]).is_err());
    assert!(run("relative_risk_predict", &[vec![1000.]], &[1., 0., 0.]).is_err());
    assert!(run("cox_predict", &[vec![1.]], &[0.5, 1.]).is_err());
    assert!(run("no_such_operation", &data(), &[]).is_err());
}
