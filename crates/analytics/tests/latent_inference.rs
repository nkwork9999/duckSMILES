use ducksmiles_analytics::execute;

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
    let tail = &j[start..];
    let end = tail.find(']').unwrap();
    if end == 0 {
        vec![]
    } else {
        tail[..end].split(',').map(|s| s.parse().unwrap()).collect()
    }
}
fn matrix(j: &str, key: &str) -> Vec<Vec<f64>> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &j[start..];
    if tail.starts_with(']') {
        return vec![];
    }
    let end = tail.find("]]").unwrap() + 1;
    tail[..end]
        .split("],")
        .map(|s| {
            s.trim_matches(['[', ']'])
                .split(',')
                .map(|v| v.parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1. + b.abs()), "{a} != {b}");
}
fn check(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        close(*a, *b, tol);
    }
}
fn distances(rows: &[Vec<f64>]) -> Vec<f64> {
    rows.iter()
        .flat_map(|a| {
            rows.iter()
                .map(move |b| a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum())
        })
        .collect()
}

#[test]
fn one_effect_manova_recovers_scatter_matrices_and_exact_rank_one_f() {
    // Two groups with identical square residual clouds and mean shift [2,1].
    let rows = (0..2)
        .flat_map(|group| {
            [-1., 1.].into_iter().flat_map(move |a| {
                [-1., 1.]
                    .into_iter()
                    .map(move |b| vec![group as f64, a + 2. * group as f64, b + group as f64])
            })
        })
        .collect::<Vec<_>>();
    let j = execute("manova", &rows, &[1.]).unwrap();
    check(
        &matrix(&j, "hypothesis_sscp").concat(),
        &[8., 4., 4., 2.],
        1e-10,
    );
    check(&matrix(&j, "error_sscp").concat(), &[8., 0., 0., 8.], 1e-10);
    check(&values(&j, "eigenvalues"), &[1.25, 0.], 1e-10);
    check(
        &[
            number(&j, "wilks_lambda"),
            number(&j, "pillai_trace"),
            number(&j, "hotelling_lawley_trace"),
            number(&j, "roy_largest_root"),
        ],
        &[4. / 9., 5. / 9., 1.25, 1.25],
        1e-10,
    );
    check(&values(&j, "f_statistics"), &[3.125; 4], 1e-10);
    check(&values(&j, "numerator_df"), &[2.; 4], 1e-10);
    check(&values(&j, "denominator_df"), &[5.; 4], 1e-10);
    check(&values(&j, "p_values"), &[0.1316872427983539; 4], 1e-10);
    assert!(j.contains("\"f_is_exact\":true"));
}

#[test]
fn multi_effect_manova_criteria_and_approximate_f_survive_units_and_offsets() {
    let rows = (0..16)
        .map(|i| {
            let w = (0..4)
                .map(|j| if (i >> j) & 1 == 0 { -1. } else { 1. })
                .collect::<Vec<_>>();
            vec![w[0], w[1], w[0] + w[2], 2. * w[1] + w[3]]
        })
        .collect::<Vec<_>>();
    let j = execute("manova", &rows, &[2.]).unwrap();
    check(&values(&j, "eigenvalues"), &[4., 1.], 1e-10);
    check(
        &[
            number(&j, "wilks_lambda"),
            number(&j, "pillai_trace"),
            number(&j, "hotelling_lawley_trace"),
            number(&j, "roy_largest_root"),
        ],
        &[0.1, 1.3, 5., 4.],
        1e-10,
    );
    check(
        &values(&j, "f_statistics"),
        &[12.973665961010276, 12.071428571428573, 13.75, 26.],
        1e-10,
    );
    check(
        &values(&j, "p_values"),
        &[
            9.205266807797948e-6,
            1.1176771830397329e-5,
            9.173713264964708e-6,
            2.8621670111997317e-5,
        ],
        1e-10,
    );
    assert!(j.contains("\"f_is_exact\":false"));
    assert!(j.contains("\"roy_p_value_is_lower_bound\":true"));
    let translated = rows
        .iter()
        .map(|r| {
            vec![
                1e12 + r[0] * 2.,
                -500. + r[1] * 0.5,
                10000. + r[2] * 0.25,
                -1e10 + r[3] * 4.,
            ]
        })
        .collect::<Vec<_>>();
    let k = execute("manova", &translated, &[2.]).unwrap();
    check(
        &values(&j, "eigenvalues"),
        &values(&k, "eigenvalues"),
        1e-10,
    );
    check(&values(&j, "p_values"), &values(&k, "p_values"), 1e-10);
}

#[test]
fn canonical_correlation_inference_reports_sequential_wilks_bartlett_approximations() {
    let j = execute(
        "canonical_correlation_inference",
        &[vec![0.8], vec![0.6]],
        &[100., 2., 2.],
    )
    .unwrap();
    check(&values(&j, "wilks_lambdas"), &[0.36 * 0.64, 0.64], 1e-12);
    check(
        &values(&j, "bartlett_chi_square"),
        &[141.6560507904787, -96.5_f64 * 0.64_f64.ln()],
        1e-12,
    );
    check(&values(&j, "degrees_of_freedom"), &[4., 1.], 1e-12);
    close(values(&j, "p_values")[0], 1.247594163207293e-29, 1e-35);
    assert!(j.contains("\"exact\":false"));
    assert!(j.contains("bartlett_asymptotic_chi_square"));
}

#[test]
fn canonical_inference_handles_zero_and_perfect_correlations_without_nonfinite_json() {
    let j = execute(
        "canonical_correlation_inference",
        &[vec![0.], vec![0.]],
        &[30., 2., 3.],
    )
    .unwrap();
    check(&values(&j, "p_values"), &[1., 1.], 1e-12);
    let j = execute(
        "canonical_correlation_inference",
        &[vec![1.]],
        &[30., 1., 1.],
    )
    .unwrap();
    assert!(j.contains("\"bartlett_chi_square\":[null]"));
    check(&values(&j, "p_values"), &[0.], 0.);
    for rows in [
        &[vec![0.3], vec![0.8]][..],
        &[vec![-0.2], vec![0.1]][..],
        &[vec![1.1], vec![0.1]][..],
    ] {
        assert!(execute("canonical_correlation_inference", rows, &[30., 2., 2.]).is_err());
    }
    assert!(execute(
        "canonical_correlation_inference",
        &[vec![0.8], vec![0.6]],
        &[3., 2., 2.]
    )
    .is_err());
}

#[test]
fn multiple_correspondence_builds_indicator_burt_inertia_and_adjustments() {
    let rows = vec![vec![0., 0.], vec![0., 0.], vec![1., 1.], vec![1., 1.]];
    let j = execute("multiple_correspondence_analysis", &rows, &[2.]).unwrap();
    check(&values(&j, "inertias"), &[1., 0.], 1e-10);
    close(number(&j, "total_inertia"), 1., 1e-12);
    check(&values(&j, "benzecri_inertias"), &[1., 0.], 1e-10);
    close(number(&j, "greenacre_total_inertia"), 1., 1e-10);
    let coords = matrix(&j, "row_coordinates");
    close(distances(&coords)[2], 4., 1e-10);
    close(distances(&coords)[1], 0., 1e-10);
    check(
        &matrix(&j, "burt_table").concat(),
        &[
            2., 0., 2., 0., 0., 2., 0., 2., 2., 0., 2., 0., 0., 2., 0., 2.,
        ],
        1e-12,
    );
    let independent = vec![vec![0., 0.], vec![0., 1.], vec![1., 0.], vec![1., 1.]];
    let k = execute("multiple_correspondence_analysis", &independent, &[2.]).unwrap();
    check(&values(&k, "inertias"), &[0.5, 0.5], 1e-10);
    check(&values(&k, "benzecri_inertias"), &[0., 0.], 1e-10);
    assert!(k.contains("\"benzecri_proportions\":[null,null]"));
    assert!(k.contains("\"greenacre_proportions\":[null,null]"));
}

#[test]
fn mca_geometry_is_invariant_to_category_codes_and_variable_permutation() {
    let rows = vec![
        vec![0., 1., 2.],
        vec![1., 0., 1.],
        vec![0., 0., 0.],
        vec![1., 1., 2.],
        vec![1., 0., 0.],
        vec![0., 1., 1.],
    ];
    let transformed = rows
        .iter()
        .map(|r| vec![-17. + 10. * r[2], 3. - 4. * r[0], 100. + r[1] * 7.])
        .collect::<Vec<_>>();
    let a = execute("multiple_correspondence_analysis", &rows, &[3.]).unwrap();
    let b = execute("multiple_correspondence_analysis", &transformed, &[3.]).unwrap();
    check(&values(&a, "inertias"), &values(&b, "inertias"), 1e-10);
    check(
        &distances(&matrix(&a, "row_coordinates")),
        &distances(&matrix(&b, "row_coordinates")),
        1e-9,
    );
}

fn equicorrelation(p: usize, r: f64) -> Vec<Vec<f64>> {
    (0..p)
        .map(|i| (0..p).map(|j| if i == j { 1. } else { r }).collect())
        .collect()
}
#[test]
fn maximum_likelihood_factor_recovers_exact_population_model_and_fit_test() {
    let j = execute(
        "factor_analysis_ml",
        &equicorrelation(4, 0.49),
        &[1., 0., 300., 1e-9, 1e-6, 1., 200.],
    )
    .unwrap();
    for loading in matrix(&j, "loadings").concat() {
        close(loading.abs(), 0.7, 2e-6);
    }
    check(&values(&j, "uniquenesses"), &[0.51; 4], 2e-6);
    close(number(&j, "criterion"), 0., 1e-10);
    close(number(&j, "model_df"), 2., 0.);
    close(number(&j, "model_p_value"), 1., 1e-8);
    assert!(j.contains("\"converged\":true"));
    assert!(j.contains("\"boundary\":false"));
}

#[test]
fn ml_factor_flags_iteration_limit_and_variance_boundary_and_rejects_unidentified_models() {
    let j = execute(
        "factor_analysis_ml",
        &equicorrelation(4, 0.49),
        &[1., 0., 1., 1e-12, 1e-6, 1., 200.],
    )
    .unwrap();
    assert!(j.contains("\"converged\":false"));
    assert!(j.contains("\"model_p_value\":null"));
    let j = execute(
        "factor_analysis_ml",
        &equicorrelation(4, 0.999),
        &[1., 0., 300., 1e-9, 0.01, 1., 200.],
    )
    .unwrap();
    assert!(j.contains("\"boundary\":true"));
    assert!(j.contains("\"model_p_value\":null"));
    assert!(execute(
        "factor_analysis_ml",
        &equicorrelation(4, 0.4),
        &[2., 0., 100., 1e-8, 1e-6, 1., 200.]
    )
    .is_err());
}

#[test]
fn latent_inference_rejects_invalid_data_and_domains() {
    for op in [
        "manova",
        "canonical_correlation_inference",
        "multiple_correspondence_analysis",
        "factor_analysis_ml",
    ] {
        assert!(execute(op, &[], &[]).is_err());
        assert!(execute(op, &[vec![f64::NAN]], &[]).is_err());
        assert!(execute(op, &[vec![1., 2.], vec![3.]], &[]).is_err());
    }
    assert!(execute(
        "multiple_correspondence_analysis",
        &[vec![0.5, 1.], vec![1., 2.]],
        &[]
    )
    .is_err());
    assert!(execute(
        "manova",
        &[vec![0., 1., 1.], vec![1., 2., 2.], vec![2., 4., 4.]],
        &[1.]
    )
    .is_err());
}

#[test]
fn ml_factor_optimizes_free_uniquenesses_when_other_variables_hit_bounds() {
    // Independent direct Gaussian likelihood optimized with SciPy L-BFGS-B.
    // A quasi-Newton direction must keep active constrained coordinates fixed.
    let loading = [0.95, 0.9, 0.7, 0.4, 0.2];
    let r = (0..5)
        .map(|i| {
            (0..5)
                .map(|j| if i == j { 1. } else { loading[i] * loading[j] })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let j = execute(
        "factor_analysis_ml",
        &r,
        &[1., 0., 1000., 1e-8, 0.2, 1., 200.],
    )
    .unwrap();
    assert!(j.contains("\"converged\":true"));
    assert!(j.contains("\"boundary\":true"));
    check(
        &values(&j, "uniquenesses"),
        &[
            0.2,
            0.2,
            0.49418874255701356,
            0.8357688979825946,
            0.9589766129896158,
        ],
        2e-6,
    );
    close(number(&j, "criterion"), 0.056184900001420957, 1e-9);
    assert!(number(&j, "projected_gradient_max") <= 1e-8);
    let j = execute(
        "factor_analysis_ml",
        &r,
        &[1., 0., 1000., 1e-8, 0.1, 1., 200.],
    )
    .unwrap();
    assert!(j.contains("\"converged\":true"));
    check(
        &values(&j, "uniquenesses"),
        &[
            0.1,
            0.18825439010858053,
            0.5095059156894611,
            0.8398562997339639,
            0.9599648160363297,
        ],
        2e-6,
    );
    close(number(&j, "criterion"), 2.3065227577667713e-5, 1e-10);
}
