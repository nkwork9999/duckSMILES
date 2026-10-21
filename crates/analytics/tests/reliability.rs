use ducksmiles_analytics::execute as run;

fn numbers(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let t = &j[start..];
    let end = if t.starts_with('[') {
        let mut depth = 0;
        t.char_indices()
            .find_map(|(i, c)| {
                if c == '[' {
                    depth += 1;
                }
                if c == ']' {
                    depth -= 1;
                    if depth == 0 {
                        return Some(i + 1);
                    }
                }
                None
            })
            .unwrap()
    } else {
        t.find([',', '}']).unwrap_or(t.len())
    };
    t[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-6 * b.abs().max(1.), "{a} != {b}");
}
fn near_all(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        near(x, y);
    }
}
fn data() -> Vec<Vec<f64>> {
    let x = [-1., -1., -1., -0.5, -0.5, 0., 0., 0.5, 0.5, 1., 1., 1.];
    let t = [1.2, 2.1, 1.6, 3.5, 2.7, 3.2, 5.5, 4.8, 7.1, 6.2, 8.4, 9.];
    let e = [1., 1., 0., 1., 0., 1., 1., 0., 1., 1., 0., 1.];
    (0..12).map(|i| vec![x[i], t[i], e[i]]).collect()
}

#[test]
fn censored_exponential_aft_uses_full_exposure_likelihood_and_observed_information() {
    // Independent SciPy root of complex-step scalar likelihood derivatives;
    // Hessian uses Richardson extrapolation of those derivatives.
    let j = run("aft_exponential", &data(), &[]).unwrap();
    near_all(
        &numbers(&j, "coefficients"),
        &[1.7795692316796587, 0.7636904242643491],
    );
    near_all(
        &numbers(&j, "standard_errors"),
        &[0.35355339059260316, 0.4791821847847168],
    );
    near(numbers(&j, "log_likelihood")[0], -22.23655385343727);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn censored_weibull_aft_estimates_location_effect_and_common_shape_jointly() {
    let j = run("aft_weibull", &data(), &[]).unwrap();
    near_all(
        &numbers(&j, "coefficients"),
        &[1.5422076020448903, 0.7390674142280002],
    );
    near(numbers(&j, "log_sigma")[0], -1.7040161775568994);
    near_all(
        &numbers(&j, "standard_errors"),
        &[
            0.06451602832658859,
            0.10708103929266566,
            0.27096300788586275,
        ],
    );
    near(numbers(&j, "log_likelihood")[0], -12.599177632359707);
}
#[test]
fn censored_lognormal_aft_includes_normal_survival_tail_in_the_likelihood() {
    let j = run("aft_lognormal", &data(), &[]).unwrap();
    near_all(
        &numbers(&j, "coefficients"),
        &[1.440541204073499, 0.7697624298743153],
    );
    near(numbers(&j, "log_sigma")[0], -1.3532392372547124);
    near_all(
        &numbers(&j, "standard_errors"),
        &[
            0.08402971150667052,
            0.10789552728670372,
            0.24414534906704044,
        ],
    );
    near(numbers(&j, "log_likelihood")[0], -13.428415667651329);
}
#[test]
fn intercept_only_exponential_aft_matches_exact_total_exposure_per_failure() {
    let j = run(
        "aft_exponential",
        &[vec![1., 1.], vec![2., 1.], vec![3., 0.], vec![4., 0.]],
        &[],
    )
    .unwrap();
    near(numbers(&j, "coefficients")[0], 5_f64.ln());
    near(numbers(&j, "log_likelihood")[0], -2. * 5_f64.ln() - 2.);
    near(
        numbers(&j, "standard_errors")[0],
        std::f64::consts::FRAC_1_SQRT_2,
    );
}
#[test]
fn aft_prediction_uses_saved_centering_and_returns_distribution_specific_quantiles() {
    // Model[distribution,p,log_time_reference,log_sigma,centered_intercept,origins,offsets,scales,slopes].
    let j = run(
        "aft_predict",
        &[vec![0.], vec![5.], vec![10.]],
        &[0., 0., 0., 0., 0., 5_f64.ln()],
    )
    .unwrap();
    near_all(
        &numbers(&j, "survival_probabilities"),
        &[1., (-1_f64).exp(), (-2_f64).exp()],
    );
    let j = run(
        "aft_predict",
        &[vec![0.5]],
        &[1., 0., 0., 0., 0., 5_f64.ln()],
    )
    .unwrap();
    near(numbers(&j, "time_quantiles")[0], 5. * 2_f64.ln());
    let j = run(
        "aft_predict",
        &[vec![0.5]],
        &[1., 1., 0., 0., 0.5_f64.ln(), 5_f64.ln()],
    )
    .unwrap();
    near(numbers(&j, "time_quantiles")[0], 5. * 2_f64.ln().sqrt());
    let j = run(
        "aft_predict",
        &[vec![5.]],
        &[0., 2., 0., 0., 0.4_f64.ln(), 5_f64.ln()],
    )
    .unwrap();
    near(numbers(&j, "survival_probabilities")[0], 0.5);
    let j = run(
        "aft_predict",
        &[vec![0.5]],
        &[1., 2., 0., 0., 0.4_f64.ln(), 5_f64.ln()],
    )
    .unwrap();
    near(numbers(&j, "time_quantiles")[0], 5.);
    let fitted = run("aft_weibull", &data(), &[]).unwrap();
    let mut args = vec![0.];
    args.extend(numbers(&fitted, "predict_parameters"));
    let alone = run("aft_predict", &[vec![0., 3.]], &args).unwrap();
    let batch = run("aft_predict", &[vec![0., 3.], vec![10., 3.]], &args).unwrap();
    near(
        numbers(&alone, "survival_probabilities")[0],
        numbers(&batch, "survival_probabilities")[0],
    );
}
#[test]
fn aft_estimates_and_saved_predictions_are_invariant_to_feature_units_and_time_units() {
    for op in ["aft_exponential", "aft_weibull", "aft_lognormal"] {
        let original = run(op, &data(), &[]).unwrap();
        let transformed = data()
            .iter()
            .map(|r| vec![100. + 3. * r[0], 1e50 * r[1], r[2]])
            .collect::<Vec<_>>();
        let fitted = run(op, &transformed, &[]).unwrap();
        let b = numbers(&original, "coefficients");
        let bb = numbers(&fitted, "coefficients");
        near(bb[1], b[1] / 3.);
        near(bb[0], b[0] + 50. * 10_f64.ln() - 100. * b[1] / 3.);
        near(
            numbers(&fitted, "log_likelihood")[0],
            numbers(&original, "log_likelihood")[0] - 8. * 50. * 10_f64.ln(),
        );
        let mut a = vec![0.];
        a.extend(numbers(&original, "predict_parameters"));
        let mut b = vec![0.];
        b.extend(numbers(&fitted, "predict_parameters"));
        near(
            numbers(
                &run("aft_predict", &[vec![0.5, 4.]], &a).unwrap(),
                "survival_probabilities",
            )[0],
            numbers(
                &run("aft_predict", &[vec![101.5, 4e50]], &b).unwrap(),
                "survival_probabilities",
            )[0],
        );
    }
}
#[test]
fn aft_rejects_false_convergence_and_regularized_curvature_does_not_claim_wald_calibration() {
    assert!(run("aft_weibull", &data(), &[0., 1., 1e-12]).is_err());
    let j = run("aft_lognormal", &data(), &[0.5]).unwrap();
    assert!(j.contains("\"p_values\":[null,null,null]"));
    assert!(j.contains("penalized_curvature_no_calibrated_inference"));
    assert!(run("aft_lognormal", &vec![vec![2., 1.]; 8], &[]).is_err());
}
#[test]
fn powerlaw_nhpp_matches_time_truncated_point_process_mle_and_growth_intensity() {
    let j = run(
        "nhpp_powerlaw",
        &[vec![1.], vec![2.], vec![4.], vec![7.]],
        &[10.],
    )
    .unwrap();
    near(numbers(&j, "shape")[0], 0.771457807511085);
    near(numbers(&j, "log_eta")[0], 0.5056049501753199);
    near(numbers(&j, "log_likelihood")[0], -7.518067436566543);
    near(numbers(&j, "shape_standard_error")[0], 0.3857289037555425);
    near(numbers(&j, "expected_events_at_end")[0], 4.);
    near(numbers(&j, "intensity_at_end")[0], 0.308583123004434);
    assert!(j.contains("\"growth_improving\":true"));
}
#[test]
fn nhpp_prediction_and_mle_preserve_time_rescaling() {
    let original = run(
        "nhpp_powerlaw",
        &[vec![1.], vec![2.], vec![4.], vec![7.]],
        &[10.],
    )
    .unwrap();
    let scaled = run(
        "nhpp_powerlaw",
        &[vec![1e100], vec![2e100], vec![4e100], vec![7e100]],
        &[1e101],
    )
    .unwrap();
    near(numbers(&scaled, "shape")[0], numbers(&original, "shape")[0]);
    near(
        numbers(&scaled, "log_eta")[0],
        numbers(&original, "log_eta")[0] + 100. * 10_f64.ln(),
    );
    let j = run(
        "nhpp_predict",
        &[vec![0.], vec![10.]],
        &[
            numbers(&original, "shape")[0],
            numbers(&original, "log_eta")[0],
        ],
    )
    .unwrap();
    near_all(&numbers(&j, "expected_cumulative_events"), &[0., 4.]);
}
#[test]
fn independent_component_series_parallel_and_nonidentical_k_of_n_are_exact() {
    let p = vec![vec![0.8, 0.7, 0.9]];
    near(
        numbers(
            &run("reliability_blocks", &p, &[0.]).unwrap(),
            "reliabilities",
        )[0],
        0.504,
    );
    near(
        numbers(
            &run("reliability_blocks", &p, &[1.]).unwrap(),
            "reliabilities",
        )[0],
        0.994,
    );
    near(
        numbers(
            &run("reliability_blocks", &p, &[2., 2.]).unwrap(),
            "reliabilities",
        )[0],
        0.902,
    );
    near(
        numbers(
            &run("reliability_blocks", &p, &[2., 1.]).unwrap(),
            "reliabilities",
        )[0],
        0.994,
    );
    near(
        numbers(
            &run("reliability_blocks", &p, &[2., 3.]).unwrap(),
            "reliabilities",
        )[0],
        0.504,
    );
    near(
        numbers(
            &run("reliability_blocks", &[vec![1., 0., 1.]], &[2., 2.]).unwrap(),
            "reliabilities",
        )[0],
        1.,
    );
}
#[test]
fn invalid_censoring_nonpositive_times_unidentified_aft_and_bad_blocks_are_errors() {
    for op in [
        "aft_exponential",
        "aft_weibull",
        "aft_lognormal",
        "aft_predict",
        "nhpp_powerlaw",
        "nhpp_predict",
        "reliability_blocks",
    ] {
        for d in [
            vec![],
            vec![vec![]],
            vec![vec![f64::NAN]],
            vec![vec![1., 2.], vec![3.]],
        ] {
            assert!(run(op, &d, &[]).is_err(), "{op}");
        }
    }
    for op in ["aft_exponential", "aft_weibull", "aft_lognormal"] {
        assert!(run(op, &vec![vec![1., 0.]; 8], &[]).is_err());
        assert!(run(op, &vec![vec![-1., 1.]; 8], &[]).is_err());
        assert!(run(op, &vec![vec![1., 0.5]; 8], &[]).is_err());
        let duplicate = data()
            .iter()
            .map(|r| vec![r[0], r[0], r[1], r[2]])
            .collect::<Vec<_>>();
        assert!(run(op, &duplicate, &[]).is_err());
    }
    assert!(run("nhpp_powerlaw", &[vec![10.], vec![10.]], &[10.]).is_err());
    assert!(run("nhpp_powerlaw", &[vec![11.], vec![2.]], &[10.]).is_err());
    assert!(run("reliability_blocks", &[vec![0.8, 1.1]], &[0.]).is_err());
    assert!(run("reliability_blocks", &[vec![0.8, 0.7]], &[2., 3.]).is_err());
}
