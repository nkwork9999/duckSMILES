use ducksmiles_analytics::execute as run;
fn values(j: &str, k: &str) -> Vec<f64> {
    let t = j.split_once(&format!("\"{k}\":")).expect(k).1;
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
fn near(x: f64, y: f64, t: f64) {
    assert!((x - y).abs() < t * (1. + y.abs()), "{x} != {y}");
}
fn training(y: &[f64]) -> Vec<Vec<f64>> {
    y.iter()
        .enumerate()
        .map(|(i, &value)| {
            let mut r = vec![value];
            r.extend((0..y.len()).map(|j| f64::from(i / 2 == j / 2)));
            r
        })
        .collect()
}
#[test]
fn genomic_reml_recovers_balanced_random_intercept_variances_and_blups() {
    let j = run(
        "kinship_model",
        &training(&[0., 2., 4., 6., 8., 10.]),
        &[0., 1.],
    )
    .unwrap();
    near(values(&j, "genetic_variance")[0], 15., 2e-5);
    near(values(&j, "residual_variance")[0], 2., 2e-5);
    near(values(&j, "intercept")[0], 5., 1e-10);
    for (&a, b) in values(&j, "genetic_values")
        .iter()
        .zip([-3.75, -3.75, 0., 0., 3.75, 3.75])
    {
        near(a, b, 2e-6);
    }
    let ll = -0.5
        * (5. * ((2. * std::f64::consts::PI).ln() + 1.)
            + 3. * 2_f64.ln()
            + 2. * 32_f64.ln()
            + 6_f64.ln());
    near(values(&j, "log_likelihood")[0], ll, 1e-8);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn genomic_ml_and_reml_have_different_variance_likelihood_contracts() {
    let j = run(
        "kinship_model",
        &training(&[0., 2., 4., 6., 8., 10.]),
        &[0., 0.],
    )
    .unwrap();
    near(values(&j, "genetic_variance")[0], 29. / 3., 2e-5);
    near(values(&j, "residual_variance")[0], 2., 2e-5);
    let ll = -0.5
        * (6. * ((2. * std::f64::consts::PI).ln() + 1.)
            + 3. * 2_f64.ln()
            + 3. * (64_f64 / 3.).ln());
    near(values(&j, "log_likelihood")[0], ll, 1e-8);
}
#[test]
fn saved_kinship_prediction_uses_fixed_cross_covariances_without_reestimating_variance() {
    let j = run(
        "kinship_model",
        &training(&[0., 2., 4., 6., 8., 10.]),
        &[0., 1.],
    )
    .unwrap();
    let p = values(&j, "predict_parameters");
    let j = run(
        "kinship_predict",
        &[vec![1., 1., 0., 0., 0., 0.], vec![0.; 6]],
        &p,
    )
    .unwrap();
    for (&a, b) in values(&j, "predictions").iter().zip([1.25, 5.]) {
        near(a, b, 2e-6);
    }
}
#[test]
fn genetic_zero_boundary_is_distinct_from_unidentifiable_identity_kernel() {
    let j = run(
        "kinship_model",
        &training(&[0., 2., 0., 2., 0., 2.]),
        &[0., 1.],
    )
    .unwrap();
    near(values(&j, "genetic_variance")[0], 0., 1e-10);
    near(values(&j, "residual_variance")[0], 1.2, 1e-8);
    assert!(j.contains("\"boundary\":true"));
    let rows = (0..4)
        .map(|i| {
            let mut r = vec![i as f64];
            r.extend((0..4).map(|j| f64::from(i == j)));
            r
        })
        .collect::<Vec<_>>();
    let j = run("kinship_model", &rows, &[0.]).unwrap();
    assert!(j.contains("\"variance_components_identified\":false"));
    assert!(j.contains("\"converged\":false"));
    // A flat profile has no local minima to refine; its grid already proves
    // that the two variance components cannot be distinguished here.
    assert!(values(&j, "likelihood_evaluations")[0] <= 130.);
}
#[test]
fn fixed_covariates_and_original_design_reml_determinant_are_preserved() {
    let mut rows = Vec::new();
    for (g, random) in [-3., -1., 1., 3.].iter().enumerate() {
        for (i, x) in [-1., 1.].iter().enumerate() {
            let error = if (g + i) % 2 == 0 { -1. } else { 1. };
            let mut row = vec![3. + 2. * x + random + error, *x];
            row.extend((0..8).map(|j| f64::from(j / 2 == g)));
            rows.push(row);
        }
    }
    let j = run("kinship_model", &rows, &[1.]).unwrap();
    near(values(&j, "intercept")[0], 3., 1e-10);
    near(values(&j, "coefficients")[0], 2., 1e-10);
    near(values(&j, "genetic_variance")[0], 16. / 3., 2e-5);
    near(values(&j, "residual_variance")[0], 8. / 3., 2e-5);
    for row in &mut rows {
        row[0] *= 3.;
        row[1] = row[1] * 7. + 100.;
    }
    let t = run("kinship_model", &rows, &[1.]).unwrap();
    near(
        values(&t, "log_likelihood")[0] - values(&j, "log_likelihood")[0],
        -6. * 3_f64.ln() - 7_f64.ln(),
        1e-9,
    );
}
#[test]
fn kernel_and_response_units_do_not_change_fitted_values() {
    let rows = training(&[0., 2., 4., 6., 8., 10.]);
    let j = run("kinship_model", &rows, &[0.]).unwrap();
    let mut changed = rows;
    for r in &mut changed {
        r[0] = r[0] * 1e-100 + 1e-98;
        for v in &mut r[1..] {
            *v *= 7.;
        }
    }
    let t = run("kinship_model", &changed, &[0.]).unwrap();
    near(
        values(&t, "genetic_variance")[0] / 1e-200 * 7.,
        values(&j, "genetic_variance")[0],
        2e-5,
    );
    for (&a, &b) in values(&t, "conditional_fitted")
        .iter()
        .zip(&values(&j, "conditional_fitted"))
    {
        near((a - 1e-98) / 1e-100, b, 2e-6);
    }
}
#[test]
fn subnormal_kernel_units_do_not_overflow_portable_weights_or_finite_predictions() {
    let tiny = 1e-320;
    let mut rows = training(&[0., 2., 4., 6., 8., 10.]);
    for r in &mut rows {
        r[0] *= 1e-100;
        for v in &mut r[1..] {
            *v *= tiny;
        }
    }
    let j = run("kinship_model", &rows, &[0.]).unwrap();
    near(
        values(&j, "genetic_variance")[0] / (15e-200 / tiny),
        1.,
        2e-5,
    );
    let params = values(&j, "predict_parameters");
    let j = run(
        "kinship_predict",
        &[vec![tiny, tiny, 0., 0., 0., 0.]],
        &params,
    )
    .unwrap();
    near(values(&j, "predictions")[0] / 1e-100, 1.25, 2e-6);
}
#[test]
fn profile_refines_a_true_interior_optimum_inside_the_last_grid_interval() {
    let mut rows = training(&[0., 2., 4., 6., 8., 10.]);
    for (i, row) in rows.iter_mut().enumerate() {
        row[i + 1] += 0.205;
    }
    let j = run("kinship_model", &rows, &[0., 0.]).unwrap();
    near(values(&j, "genetic_variance")[0], 29. / 3., 2e-5);
    near(values(&j, "residual_variance")[0], 11. / 600., 2e-5);
    near(
        values(&j, "normalized_genetic_fraction")[0],
        6989. / 7000.,
        1e-7,
    );
    assert!(j.contains("\"upper_bound_active\":false"));
}
#[test]
fn kinship_contract_rejects_indefinite_asymmetric_and_invalid_dimensions() {
    assert!(run("kinship_model", &training(&[0., 2., 4., 6.]), &[99.]).is_err());
    assert!(run("kinship_predict", &[vec![1.]], &[]).is_err());
}
#[test]
fn exactly_explained_response_rejects_a_spurious_finite_variance_mle() {
    let rows = (0..8)
        .map(|i| {
            let x = i as f64 - 2.;
            let z = (i % 2) as f64;
            let mut row = vec![3. + 2. * x + 4. * z, x, z];
            row.extend((0..8).map(|j| f64::from(i / 2 == j / 2) + if i == j { 0.2 } else { 0. }));
            row
        })
        .collect::<Vec<_>>();
    assert!(run("kinship_model", &rows, &[2.]).is_err());
}
#[test]
fn unrepresentable_positive_variance_does_not_appear_as_a_zero_variance_model() {
    let mut rows = training(&[0., 2., 4., 6., 8., 10.]);
    for row in &mut rows {
        row[0] *= 1e-200;
    }
    assert!(run("kinship_model", &rows, &[0.]).is_err());
}

#[test]
fn clipped_kernel_null_space_is_absent_from_blup_and_saved_training_predictions() {
    let noise = 1e-5;
    let mut rows = training(&[0., 0., 1., 1., 2., 2.]);
    let z = [
        std::f64::consts::FRAC_1_SQRT_2,
        -std::f64::consts::FRAC_1_SQRT_2,
        0.,
        0.,
        0.,
        0.,
    ];
    for (i, row) in rows.iter_mut().enumerate() {
        row[0] += noise * z[i];
        for j in 0..6 {
            row[j + 1] -= 1e-12 * z[i] * z[j];
        }
    }
    let j = run(
        "kinship_model",
        &rows,
        &[0., 1., 0.999999999999, 128., 1e-12],
    )
    .unwrap();
    assert!(j.contains("\"converged\":true"));
    let genetic = values(&j, "genetic_values");
    let null_contrast: f64 = genetic.iter().zip(z).map(|(&g, z)| g * z).sum();
    assert!(
        null_contrast.abs() < 1e-11,
        "clipped null-space BLUP {null_contrast}"
    );
    let query: Vec<_> = rows.iter().map(|r| r[1..].to_vec()).collect();
    let predicted = run("kinship_predict", &query, &values(&j, "predict_parameters")).unwrap();
    for (&a, b) in values(&predicted, "predictions")
        .iter()
        .zip(values(&j, "conditional_fitted"))
    {
        near(a, b, 1e-11);
    }
}
#[test]
fn subnormal_fixed_feature_scale_preserves_finite_zero_and_nonzero_slope_predictions() {
    for last in [1., 1. + 1e-14] {
        let mut rows = training(&[0., 2., 2., 0., 1., last]);
        for (i, row) in rows.iter_mut().enumerate() {
            row.insert(1, if i % 2 == 0 { 0. } else { 1e-320 });
        }
        let j = run("kinship_model", &rows, &[1., 1.]).unwrap();
        assert!(j.contains("\"converged\":true"));
        let expected = values(&j, "intercept")[0] + values(&j, "coefficients")[0];
        let scored = run(
            "kinship_predict",
            &[vec![1., 0., 0., 0., 0., 0., 0.]],
            &values(&j, "predict_parameters"),
        )
        .unwrap();
        near(values(&scored, "predictions")[0], expected, 2e-14);
    }
}
#[test]
fn subnormal_kernel_prediction_scales_physical_effect_before_ratio_overflow() {
    let mut rows = training(&[0., 2., 0., 2., 0., 2.]);
    for row in &mut rows {
        for v in &mut row[1..] {
            *v *= 1e-320;
        }
    }
    let j = run("kinship_model", &rows, &[0., 1.]).unwrap();
    let predicted = run(
        "kinship_predict",
        &[vec![1.; 6]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    near(values(&predicted, "predictions")[0], 1., 1e-14);
    // Valid portable layout: finite physical term despite overflow in cross/K.
    let params = [1., 4., 0., 0., 0., 1e-100, 0., 1e-320, 1e-200, 0., 0., 0.];
    let j = run("kinship_predict", &[vec![1., 0., 0., 0.]], &params).unwrap();
    let expected = 1e-100 * (1e-200 / 1e-320);
    near(values(&j, "predictions")[0], expected, 2e-14);
}
#[test]
fn physical_prediction_keeps_underflowed_products_and_cancelling_large_terms() {
    let params = [1., 4., 0., 0., 0., 1e-200, 0., 1e-100, 1e-200, 0., 0., 0.];
    let j = run("kinship_predict", &[vec![1e150, 0., 0., 0.]], &params).unwrap();
    near(values(&j, "predictions")[0] / 1e-150, 1., 2e-14);
    let params = [1., 4., 0., 0., 0., 1e-100, 0., 1e-320, 1., 1e-300, -1., 0.];
    let j = run("kinship_predict", &[vec![1., 1., 1., 0.]], &params).unwrap();
    near(
        values(&j, "predictions")[0] / (1e-100 * (1e-300 / 1e-320)),
        1.,
        2e-14,
    );
}
#[test]
fn valid_training_dimensions_exercise_psd_and_exact_symmetry_validation() {
    let mut rows = training(&[0., 2., 4., 6.]);
    rows[0][2] = 2.;
    rows[1][1] = 2.;
    assert!(run("kinship_model", &rows, &[0.])
        .unwrap_err()
        .contains("positive semidefinite"));
    let mut rows = training(&[0., 2., 4., 6.]);
    rows[0][2] += 5e-13;
    assert!(run("kinship_model", &rows, &[0.])
        .unwrap_err()
        .contains("exactly symmetric"));
}
