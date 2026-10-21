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
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-5 * b.abs().max(1.), "{a} != {b}");
}
fn all(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        near(x, y);
    }
}
fn times(t: &[f64]) -> Vec<Vec<f64>> {
    t.iter().map(|&v| vec![v]).collect()
}
fn slope_rows() -> Vec<Vec<f64>> {
    let mut rows = Vec::new();
    for (g, n) in [3, 4, 5, 3, 4, 5, 3, 5].into_iter().enumerate() {
        for j in 0..n {
            let u = -1.5 + 3. * j as f64 / (n - 1) as f64 + 0.07 * g as f64;
            let y = 2.3
                + 1.2 * u
                + 0.9 * (g as f64 * 1.17).sin()
                + 0.6 * (g as f64 * 0.91).cos() * u
                + 0.2 * (j as f64 * 2.31 + g as f64 * 0.57).sin();
            rows.push(vec![u, y, g as f64]);
        }
    }
    rows
}
fn gls_rows() -> Vec<Vec<f64>> {
    vec![
        vec![0., 1.1, 0., 0.],
        vec![1., 1.8, 0., 1.],
        vec![2., 2.1, 0., 2.],
        vec![0., 0.8, 1., 0.],
        vec![2., 2.4, 1., 2.],
        vec![0., 1.4, 2., 0.],
        vec![1., 1.9, 2., 1.],
        vec![2., 2.8, 2., 2.],
    ]
}

#[test]
fn nine_repeated_covariance_structures_preserve_variances_and_distinct_correlations() {
    let fixtures = [
        (
            vec![0., 1., 4., 9.],
            vec![1., 0., 0., 0., 4., 0., 0., 0., 9.],
        ),
        (
            vec![1., 1., 0.2, 4., 0.3, 0.4, 9.],
            vec![1., 0.2, 0.3, 0.2, 4., 0.4, 0.3, 0.4, 9.],
        ),
        (vec![2., 4., 0.25], vec![4., 1., 1., 1., 4., 1., 1., 1., 4.]),
        (
            vec![3., 1., 4., 9., 0.25],
            vec![1., 0.5, 0.75, 0.5, 4., 1.5, 0.75, 1.5, 9.],
        ),
        (vec![4., 4., 0.5], vec![4., 2., 1., 2., 4., 2., 1., 2., 4.]),
        (
            vec![5., 4., 0.4, 0.1],
            vec![4., 1.6, 0.4, 1.6, 4., 1.6, 0.4, 1.6, 4.],
        ),
        (
            vec![6., 1., 4., 9., 0.4, 0.1],
            vec![1., 0.8, 0.3, 0.8, 4., 2.4, 0.3, 2.4, 9.],
        ),
        (
            vec![7., 1., 4., 9., 0.5, -0.2],
            vec![1., 1., -0.3, 1., 4., -1.2, -0.3, -1.2, 9.],
        ),
        (
            vec![8., 4., 0.5, -0.2],
            vec![4., 2., -0.4, 2., 4., -0.8, -0.4, -0.8, 4.],
        ),
    ];
    for (args, expected) in fixtures {
        let j = run("repeated_covariance", &times(&[0., 1., 2.]), &args).unwrap();
        all(&numbers(&j, "covariance"), &expected);
    }
}
#[test]
fn ar1_uses_actual_nonuniform_elapsed_times_and_negative_rho_requires_integer_lags() {
    let j = run(
        "repeated_covariance",
        &times(&[0., 0.5, 2.]),
        &[4., 4., 0.25],
    )
    .unwrap();
    all(
        &numbers(&j, "covariance"),
        &[4., 2., 0.25, 2., 4., 0.5, 0.25, 0.5, 4.],
    );
    let j = run(
        "repeated_covariance",
        &times(&[0., 1., 3.]),
        &[4., 1., -0.5],
    )
    .unwrap();
    all(
        &numbers(&j, "covariance"),
        &[1., -0.5, -0.125, -0.5, 1., 0.25, -0.125, 0.25, 1.],
    );
    assert!(run(
        "repeated_covariance",
        &times(&[0., 0.5, 2.]),
        &[4., 1., -0.5]
    )
    .is_err());
    assert!(run(
        "repeated_covariance",
        &times(&[0., 0.5, 2.]),
        &[5., 1., 0.4, 0.1]
    )
    .is_err());
}
#[test]
fn known_covariance_gls_does_not_refit_scale_and_matches_dense_numpy_likelihood() {
    let j = run(
        "repeated_model",
        &gls_rows(),
        &[4., 2., 300., 1e-8, 0.7, 0.4],
    )
    .unwrap();
    all(
        &numbers(&j, "coefficients"),
        &[1.1040760869565218, 0.6666666666666667],
    );
    all(
        &numbers(&j, "coefficient_covariance"),
        &[0.226713768115942, -0.098, -0.098, 0.098],
    );
    near(numbers(&j, "log_likelihood")[0], -5.878572087085805);
    near(numbers(&j, "estimated_scale")[0], 1.);
}
#[test]
fn missing_visits_use_global_antedependent_path_instead_of_compressing_adjacency() {
    let j = run(
        "repeated_model",
        &gls_rows(),
        &[7., 2., 300., 1e-8, 1., 4., 9., 0.5, -0.2],
    )
    .unwrap();
    all(
        &numbers(&j, "covariance_by_level"),
        &[1., 1., -0.3, 1., 4., -1.2, -0.3, -1.2, 9.],
    );
}
#[test]
fn correlated_random_slopes_profile_ml_and_reml_with_unbalanced_subjects() {
    // Independent full dense covariance eigen/inverse NumPy likelihood and SciPy BFGS.
    for (method, beta, sigma, g, ll) in [
        (
            0.,
            [2.471539921441772, 1.27060673313422],
            0.03166602428341704,
            [
                0.36543285296122274,
                0.19796166464292883,
                0.19796166464292883,
                0.2124989163371142,
            ],
            -17.133021424045282,
        ),
        (
            1.,
            [2.471342534219114, 1.2706739077161133],
            0.031710741349030264,
            [
                0.4183622514310185,
                0.22608051084138833,
                0.22608051084138833,
                0.24371382904636107,
            ],
            -18.879829725553048,
        ),
    ] {
        let j = run("random_slope_model", &slope_rows(), &[method]).unwrap();
        all(&numbers(&j, "coefficients"), &beta);
        near(numbers(&j, "residual_variance")[0], sigma);
        all(&numbers(&j, "random_effect_covariance"), &g);
        near(numbers(&j, "log_likelihood")[0], ll);
        assert!(j.contains("\"converged\":true"));
    }
}
#[test]
fn random_slope_blups_predict_known_groups_and_marginal_mean_for_new_groups() {
    let data = slope_rows();
    let j = run("random_slope_model", &data, &[]).unwrap();
    let par = numbers(&j, "predict_parameters");
    let pred_rows: Vec<_> = data.iter().map(|r| vec![r[0], r[2]]).collect();
    let p = run("random_slope_predict", &pred_rows, &par).unwrap();
    all(
        &numbers(&p, "predictions"),
        &numbers(&j, "conditional_fitted_values"),
    );
    let b = numbers(&j, "coefficients");
    let p = run("random_slope_predict", &[vec![0.25, 999.]], &par).unwrap();
    near(numbers(&p, "predictions")[0], b[0] + 0.25 * b[1]);
    assert!(p.contains("\"known_subjects\":[false]"));
}
#[test]
fn random_slope_fit_and_blups_are_invariant_to_row_permutation() {
    let mut data = slope_rows();
    let a = run("random_slope_model", &data, &[]).unwrap();
    data.reverse();
    let b = run("random_slope_model", &data, &[]).unwrap();
    all(&numbers(&a, "coefficients"), &numbers(&b, "coefficients"));
    all(
        &numbers(&a, "random_effect_covariance"),
        &numbers(&b, "random_effect_covariance"),
    );
    near(
        numbers(&a, "log_likelihood")[0],
        numbers(&b, "log_likelihood")[0],
    );
}
#[test]
fn known_gls_permutation_and_feature_translation_preserve_fitted_values() {
    let mut data = gls_rows();
    let args = [4., 2., 300., 1e-8, 0.7, 0.4];
    let a = run("repeated_model", &data, &args).unwrap();
    for r in &mut data {
        r[0] += 1e8;
        r[1] = r[1] * 3. + 7.;
    }
    let b = run("repeated_model", &data, &args).unwrap();
    for (x, y) in numbers(&a, "fitted_values")
        .iter()
        .zip(numbers(&b, "fitted_values"))
    {
        near(x * 3. + 7., y);
    }
}
#[test]
fn invalid_covariance_duplicate_visits_rank_deficiency_and_work_limits_are_rejected() {
    assert!(run(
        "repeated_covariance",
        &times(&[0., 1., 2.]),
        &[5., 1., 0.99, -0.99]
    )
    .is_err());
    assert!(run("repeated_covariance", &times(&[0., 0.]), &[2., 1., 0.2]).is_err());
    assert!(run("repeated_covariance", &times(&[0., 1.]), &[0., 1., 0.]).is_err());
    let mut d = gls_rows();
    d.push(d[0].clone());
    assert!(run("repeated_model", &d, &[4.]).is_err());
    let d = vec![vec![1., 2., 0., 0.]; 257];
    assert!(run("repeated_model", &d, &[4.]).is_err());
    let d = vec![vec![1., 2., 0.]; 8];
    assert!(run("random_slope_model", &d, &[]).is_err());
    assert!(run("random_slope_predict", &[vec![0., 0.]], &[]).is_err());
}
#[test]
fn iteration_limit_is_reported_as_unconverged_instead_of_success() {
    let j = run("random_slope_model", &slope_rows(), &[1., 1., 1e-12]).unwrap();
    assert!(j.contains("\"converged\":false"));
}

fn repeated_rows() -> Vec<Vec<f64>> {
    let mut rows = Vec::new();
    for g in 0..16 {
        for t in 0..4 {
            if g % 3 == 0 && t == 1 {
                continue;
            }
            let gf = g as f64;
            let tf = t as f64;
            let y = 1.1
                + 0.7 * (g % 2) as f64
                + 0.35 * tf
                + 0.6 * (gf * 0.7).sin()
                + 0.3 * (gf * 1.3).cos() * tf
                + 0.5 * (gf * 1.7 + tf * 2.1).sin();
            rows.push(vec![(g % 2) as f64, tf, y, gf, tf]);
        }
    }
    rows
}
#[test]
fn all_nine_structures_ml_and_reml_match_independent_full_covariance_likelihood() {
    // SciPy SLSQP/Nelder-Mead optimizes every physical scale directly, unlike
    // the implementation's profiled scale and SPD shape coordinates.
    let expected = [
        -54.8065272187293,
        -58.994337440483065,
        -46.65022050193827,
        -50.633217019802494,
        -50.22869243478408,
        -53.836525259412966,
        -48.735938159961314,
        -52.45345971778257,
        -51.40043369748761,
        -54.90251958354564,
        -49.30751000998688,
        -52.96615457333559,
        -47.848561143336866,
        -51.64944381295927,
        -49.705941500618835,
        -53.360524360207215,
        -51.31519400810896,
        -54.85657646538943,
    ];
    for code in 0..9 {
        for method in 0..2 {
            let j = run(
                "repeated_model",
                &repeated_rows(),
                &[code as f64, method as f64],
            )
            .unwrap();
            near(
                numbers(&j, "log_likelihood")[0],
                expected[code * 2 + method],
            );
            assert!(
                j.contains("\"optimizer_converged\":true"),
                "{code}/{method}: {j}"
            );
            assert!(
                j.contains("\"covariance_parameters_identified\":true"),
                "{code}/{method}: {j}"
            );
        }
    }
}
#[test]
fn negative_ar1_at_large_integral_lag_does_not_truncate_to_i32() {
    let rho = -0.999999999999;
    let j = run(
        "repeated_covariance",
        &times(&[0., 3_000_000_000.]),
        &[4., 1., rho],
    )
    .unwrap();
    near(numbers(&j, "covariance")[1], (-rho).powf(3_000_000_000.));
}
#[test]
fn reported_negative_ar1_parameter_is_preserved_when_first_time_gap_is_even() {
    let mut data = gls_rows();
    for row in &mut data {
        row[3] = [0., 2., 3.][row[3] as usize];
    }
    let j = run("repeated_model", &data, &[4., 2., 300., 1e-8, 0.7, -0.4]).unwrap();
    near(numbers(&j, "covariance_parameters")[1], -0.4);
}
#[test]
fn zero_random_variance_boundary_is_explicit_instead_of_interior_convergence() {
    let mut rows = Vec::new();
    for g in 0..8 {
        for j in 0..5 {
            let u = j as f64 - 2.;
            let y = 2. + 0.7 * u + 0.2 * (j as f64 * 2.1).sin();
            rows.push(vec![u, y, g as f64]);
        }
    }
    let j = run("random_slope_model", &rows, &[]).unwrap();
    assert!(j.contains("\"variance_boundary\":true"), "{j}");
    assert!(j.contains("\"converged\":false"), "{j}");
}
#[test]
fn unobserved_unstructured_time_pair_cannot_be_reported_as_identified() {
    let mut rows = Vec::new();
    for g in 0..8 {
        for j in 0..2 {
            let t = (j + g % 2) as f64;
            rows.push(vec![1. + 0.4 * t + (g as f64 * 0.7 + t).sin(), g as f64, t]);
        }
    }
    let j = run("repeated_model", &rows, &[1.]).unwrap();
    assert!(
        j.contains("\"covariance_parameters_identified\":false"),
        "{j}"
    );
    assert!(j.contains("\"converged\":false"), "{j}");
}
#[test]
fn duplicated_fixed_feature_is_rejected_by_gls_rank_check() {
    let data: Vec<_> = gls_rows()
        .into_iter()
        .map(|r| vec![r[0], 2. * r[0], r[1], r[2], r[3]])
        .collect();
    assert!(run("repeated_model", &data, &[4., 2., 300., 1e-8, 0.7, 0.4]).is_err());
}
