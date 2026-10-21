use ducksmiles_analytics::execute as run;
fn vals(s: &str, k: &str) -> Vec<f64> {
    let t = s.split_once(&format!("\"{k}\":")).expect(k).1;
    let end = if t.starts_with('[') {
        let mut d = 0;
        t.char_indices()
            .find_map(|(i, c)| {
                if c == '[' {
                    d += 1;
                }
                if c == ']' {
                    d -= 1;
                    if d == 0 {
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
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().unwrap())
        .collect()
}
fn near(x: f64, y: f64) {
    assert!((x - y).abs() < 1e-9 * (1. + y.abs()), "{x} != {y}");
}
#[test]
fn covariance_uses_four_published_correlations_and_independent_record_nugget() {
    for (kind, rho, want) in [
        (0., 2., 1.25),
        (1., 2., 4. * (-0.5_f64).exp()),
        (2., 2., 4. * (-0.25_f64).exp()),
        (3., 0.25, 1.),
    ] {
        let j = run(
            "spatial_covariance",
            &[vec![0.], vec![1.], vec![0.]],
            &[kind, 4., 1., rho],
        )
        .unwrap();
        let c = vals(&j, "covariance");
        near(c[0], 5.);
        near(c[1], want);
        near(c[2], 4.);
        near(c[4], 5.);
        near(vals(&j, "semivariance")[2], 1.);
    }
}
#[test]
fn spherical_cutoff_and_three_dimensional_euclidean_distance() {
    let j = run(
        "spatial_covariance",
        &[vec![0., 0., 0.], vec![1., 2., 2.]],
        &[0., 2., 0., 3.],
    )
    .unwrap();
    near(vals(&j, "covariance")[1], 0.);
}
#[test]
fn empirical_variogram_counts_pairs_and_uses_half_squared_differences() {
    let j = run(
        "spatial_variogram",
        &[vec![0., 0.], vec![1., 2.], vec![2., 4.]],
        &[2.],
    )
    .unwrap();
    assert_eq!(vals(&j, "pair_counts"), vec![2., 1.]);
    assert_eq!(vals(&j, "mean_distances"), vec![1., 2.]);
    assert_eq!(vals(&j, "semivariances"), vec![2., 8.]);
}
#[test]
fn empirical_variogram_keeps_empty_bins_null_and_collocated_pairs() {
    let j = run(
        "spatial_variogram",
        &[vec![0., 0.], vec![1., 2.], vec![2., 4.]],
        &[3.],
    )
    .unwrap();
    assert!(j.contains("\"semivariances\":[2,null,8]"));
    let j = run("spatial_variogram", &[vec![0., 0.], vec![0., 2.]], &[]).unwrap();
    near(vals(&j, "semivariances")[0], 2.);
    assert_eq!(vals(&j, "pair_counts"), vec![1.]);
}
#[test]
fn supplied_covariance_gls_keeps_physical_variance_and_gaussian_likelihood() {
    let j = run(
        "spatial_gls",
        &[vec![0., 1.], vec![1., 2.], vec![2., 6.]],
        &[1., 0., 0., 2., 1., 1.],
    )
    .unwrap();
    near(vals(&j, "intercept")[0], 3.);
    near(vals(&j, "coefficient_covariance")[0], 1.);
    near(
        vals(&j, "log_likelihood")[0],
        -0.5 * (3. * (2. * std::f64::consts::PI).ln() + 3. * 3_f64.ln() + 14. / 3.),
    );
}
#[test]
fn universal_kriging_keeps_nugget_out_of_cross_covariance_and_in_observation_variance() {
    let j = run(
        "spatial_gls",
        &[vec![0., 1.], vec![1., 2.], vec![2., 6.]],
        &[1., 0., 0., 2., 1., 1.],
    )
    .unwrap();
    let j = run(
        "spatial_predict",
        &[vec![0.], vec![10.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    let pred = vals(&j, "predictions");
    near(pred[0], 5. / 3.);
    near(pred[1], 3.);
    let v = vals(&j, "latent_variances");
    near(v[0], 7. / 9.);
    near(v[1], 3.);
    let v = vals(&j, "observation_variances");
    near(v[0], 16. / 9.);
    near(v[1], 4.);
}
#[test]
fn gls_trend_and_saved_prediction_account_for_fixed_effect_uncertainty() {
    let j = run(
        "spatial_gls",
        &[vec![0., -1., 1.], vec![1., 0., 4.], vec![2., 1., 5.]],
        &[1., 1., 0., 2., 1., 1.],
    )
    .unwrap();
    near(vals(&j, "intercept")[0], 10. / 3.);
    near(vals(&j, "coefficients")[0], 2.);
    let j = run(
        "spatial_predict",
        &[vec![10., 2.]],
        &vals(&j, "predict_parameters"),
    )
    .unwrap();
    near(vals(&j, "predictions")[0], 22. / 3.);
    near(vals(&j, "latent_variances")[0], 9.);
    near(vals(&j, "observation_variances")[0], 10.);
}
#[test]
fn known_covariance_accepts_exact_fixed_fit_without_inventing_residual_scale() {
    let j = run(
        "spatial_gls",
        &[vec![0., 0., 2.], vec![1., 1., 5.], vec![2., 2., 8.]],
        &[1., 1., 0., 2., 1., 1.],
    )
    .unwrap();
    near(vals(&j, "coefficients")[0], 3.);
    near(vals(&j, "mahalanobis_sse")[0], 0.);
    near(vals(&j, "coefficient_covariance")[3], 1.5);
}
#[test]
fn spatial_rejects_nonidentifiable_design_singular_kernel_and_bad_parameter_domains() {
    for a in [
        vec![0., -1., 0., 1.],
        vec![3., 1., 0., 1.],
        vec![0., 1., -1., 1.],
        vec![0., 1., 0., 0.],
    ] {
        assert!(run("spatial_covariance", &[vec![0.], vec![1.]], &a).is_err());
    }
    assert!(run(
        "spatial_gls",
        &[vec![0., 1.], vec![0., 2.], vec![1., 3.]],
        &[1., 0., 1., 1., 0., 1.]
    )
    .is_err());
    assert!(run(
        "spatial_gls",
        &[vec![0., 1., 1.], vec![1., 1., 2.], vec![2., 1., 3.]],
        &[1., 1., 1., 1., 1., 1.]
    )
    .is_err());
    assert!(run("spatial_predict", &[vec![0.]], &[]).is_err());
}

#[test]
fn malformed_saved_model_rejects_material_negative_variance_even_when_error_bound_sum_overflows() {
    let args = [
        1., 2., 1., 0., 0., 1e100, 0., 1., 0., 0., 0., 1., 1., 0., 0., 1.2e-104, 0., 1.2e-104,
        5e-209,
    ];
    assert!(run("spatial_predict", &[vec![0.]], &args).is_err());
}

#[test]
fn semivariance_preserves_independent_nugget_despite_large_process_sill() {
    let j = run(
        "spatial_covariance",
        &[vec![0.], vec![0.]],
        &[1., 1e16, 1., 1.],
    )
    .unwrap();
    near(vals(&j, "semivariance")[1], 1.);
}
#[test]
fn near_zero_distance_semivariances_do_not_subtract_rounded_correlations() {
    for (kind, d, rho, expected) in [
        (0., 1e-20, 1., 1.5e-20),
        (1., 1e-20, 1., 1e-20),
        (2., 1e-10, 1., 1e-20),
        (3., 1e-20, 0.5, 0.5_f64.ln().abs() * 1e-20),
    ] {
        let j = run(
            "spatial_covariance",
            &[vec![0.], vec![d]],
            &[kind, 1., 0., rho],
        )
        .unwrap();
        near(vals(&j, "semivariance")[1] / expected, 1.);
    }
}

#[test]
fn semivariance_keeps_finite_physical_products_when_dimensionless_distance_underflows() {
    for (kind, d, rho, sill, want) in [
        (0., 5e-324, 1e100, 1e100, 1e-323),
        (1., 5e-324, 1e100, 1e100, 5e-324),
        (2., 1e-200, 1., 1e100, 1e-300),
        (3., 5e-324, 0.5, 1e100, 5e-324 * 1e100 * 2_f64.ln()),
    ] {
        let j = run(
            "spatial_covariance",
            &[vec![0.], vec![d]],
            &[kind, sill, 0., rho],
        )
        .unwrap();
        near(vals(&j, "semivariance")[1] / want, 1.);
    }
}

#[test]
fn noiseless_training_record_has_exact_zero_prediction_variance_even_at_huge_sill() {
    for sill in [1., 1e100] {
        let j = run(
            "spatial_gls",
            &[vec![0., 1.], vec![1., 2.], vec![2., 3.]],
            &[1., 0., 1., sill, 0., 1.],
        )
        .unwrap();
        let pred = run(
            "spatial_predict",
            &[vec![0.], vec![1.], vec![2.]],
            &vals(&j, "predict_parameters"),
        )
        .unwrap();
        assert_eq!(vals(&pred, "latent_variances"), vec![0.; 3]);
        assert_eq!(vals(&pred, "observation_variances"), vec![0.; 3]);
    }
}

#[test]
fn physical_coefficient_covariance_scales_before_products_that_would_overflow() {
    let j = run(
        "spatial_gls",
        &[
            vec![0., -1e-200, -1e-200],
            vec![1., 0., 0.],
            vec![2., 1e-200, 1e-200],
        ],
        &[1., 1., 0., 1e-100, 0., 1.],
    )
    .unwrap();
    near(vals(&j, "coefficients")[0], 1.);
    let c = vals(&j, "coefficient_covariance");
    near(c[0] / (1e-100 / 3.), 1.);
    near(c[3] / 5e299, 1.);
}
