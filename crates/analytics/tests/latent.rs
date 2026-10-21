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
        return vec![];
    }
    tail[..end].split(',').map(|s| s.parse().unwrap()).collect()
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
                .map(|x| x.parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
fn check_array(a: &[f64], b: &[f64], t: f64) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        close(*x, *y, t);
    }
}

#[test]
fn classical_mds_recovers_euclidean_distances_and_gower_eigenvalues() {
    // Mutation caught: forgetting squared-distance centering or sqrt(eigenvalue).
    let d = [vec![0., 1., 2.], vec![1., 0., 1.], vec![2., 1., 0.]];
    let j = execute("mds", &d, &[1., 1.]).unwrap();
    close(values(&j, "eigenvalues")[0], 2., 1e-10);
    let coords = matrix(&j, "coordinates");
    close((coords[0][0] - coords[1][0]).abs(), 1., 1e-10);
    close((coords[0][0] - coords[2][0]).abs(), 2., 1e-10);
    close(number(&j, "normalized_raw_stress"), 0., 1e-10);
    let non_euclidean = execute(
        "mds",
        &[vec![0., 1., 3.], vec![1., 0., 1.], vec![3., 1., 0.]],
        &[1., 1.],
    )
    .unwrap();
    assert!(number(&non_euclidean, "negative_eigenvalue_mass") > 0.);
}

#[test]
fn correspondence_analysis_obeys_chi_square_inertia_and_profile_geometry() {
    let j = execute(
        "correspondence_analysis",
        &[vec![20., 10.], vec![10., 20.]],
        &[1.],
    )
    .unwrap();
    close(values(&j, "singular_values")[0], 1. / 3., 1e-10);
    close(number(&j, "total_inertia"), 1. / 9., 1e-10);
    close(number(&j, "chi_square"), 20. / 3., 1e-10);
    check_array(&values(&j, "row_masses"), &[0.5, 0.5], 1e-12);
    let coordinates = matrix(&j, "row_coordinates");
    close(coordinates[0][0].abs(), 1. / 3., 1e-10);
    close(coordinates[0][0], -coordinates[1][0], 1e-10);
    let scaled = execute(
        "correspondence_analysis",
        &[vec![200., 100.], vec![100., 200.]],
        &[1.],
    )
    .unwrap();
    close(
        number(&scaled, "total_inertia"),
        number(&j, "total_inertia"),
        1e-12,
    );
}

fn equicorrelation() -> Vec<Vec<f64>> {
    vec![
        vec![1., 0.64, 0.64],
        vec![0.64, 1., 0.64],
        vec![0.64, 0.64, 1.],
    ]
}

#[test]
fn principal_axis_factor_analysis_recovers_known_common_factor_and_score_weights() {
    let j = execute(
        "factor_analysis",
        &equicorrelation(),
        &[1., 0., 1000., 1e-10, 1., 1.],
    )
    .unwrap();
    assert!(j.contains("\"converged\":true"));
    check_array(&values(&j, "communalities"), &[0.64; 3], 1e-8);
    check_array(&values(&j, "uniquenesses"), &[0.36; 3], 1e-8);
    for row in matrix(&j, "loadings") {
        close(row[0].abs(), 0.8, 1e-8);
    }
    let models = matrix(&j, "score_predict_parameters");
    let scored = execute("glm_predict", &[vec![1., 1., 1.]], &models[0]).unwrap();
    close(
        values(&scored, "predictions")[0].abs(),
        3. * 0.8 / 2.28,
        1e-8,
    );
}

#[test]
fn varimax_rotation_preserves_common_covariance_and_is_orthogonal() {
    let r = [
        vec![1., 0.64, 0., 0.],
        vec![0.64, 1., 0., 0.],
        vec![0., 0., 1., 0.64],
        vec![0., 0., 0.64, 1.],
    ];
    let j = execute("factor_analysis", &r, &[2., 1., 1000., 1e-10, 1., 1.]).unwrap();
    let rotation = matrix(&j, "rotation_matrix");
    for i in 0..2 {
        for k in 0..2 {
            close(
                (0..2).map(|a| rotation[a][i] * rotation[a][k]).sum(),
                if i == k { 1. } else { 0. },
                1e-9,
            );
        }
    }
    check_array(&values(&j, "communalities"), &[0.64; 4], 1e-8);
    let reproduced = matrix(&j, "reproduced_correlation");
    close(reproduced[0][1], 0.64, 1e-8);
    close(reproduced[0][2], 0., 1e-8);
}

fn canonical_data() -> Vec<Vec<f64>> {
    (0..8)
        .map(|i| {
            let x1 = if i & 1 == 0 { -1. } else { 1. };
            let x2 = if i & 2 == 0 { -1. } else { 1. };
            let e1 = if i & 4 == 0 { -1. } else { 1. };
            let e2 = x1 * x2;
            vec![x1, x2, 0.8 * x1 + 0.6 * e1, 0.6 * x2 + 0.8 * e2]
        })
        .collect()
}

#[test]
fn canonical_correlation_whitens_both_blocks_and_reuses_original_unit_scoring() {
    let data = canonical_data();
    let j = execute("canonical_correlation", &data, &[2., 2., 0.]).unwrap();
    check_array(&values(&j, "canonical_correlations"), &[0.8, 0.6], 1e-9);
    check_array(&values(&j, "score_correlations"), &[0.8, 0.6], 1e-9);
    let models = matrix(&j, "x_score_predict_parameters");
    let scored = execute(
        "glm_predict",
        &data.iter().map(|r| r[..2].to_vec()).collect::<Vec<_>>(),
        &models[0],
    )
    .unwrap();
    check_array(
        &values(&scored, "predictions"),
        &matrix(&j, "x_scores")
            .iter()
            .map(|r| r[0])
            .collect::<Vec<_>>(),
        1e-9,
    );
    let singular = [
        vec![0., 0., 1.],
        vec![1., 1., 2.],
        vec![2., 2., 1.],
        vec![3., 3., 4.],
    ];
    assert!(execute("canonical_correlation", &singular, &[2., 1., 0.]).is_err());
    assert!(execute("canonical_correlation", &singular, &[2., 1., 0.1]).is_ok());
}

#[test]
fn gaussian_mixture_one_component_matches_population_mle_and_predicts_density() {
    let j = execute(
        "gaussian_mixture",
        &[vec![-1.], vec![0.], vec![1.]],
        &[1., 100., 1e-10, 0.01, 42.],
    )
    .unwrap();
    close(matrix(&j, "component_means")[0][0], 0., 1e-10);
    close(
        matrix(&j, "component_covariances_flat")[0][0],
        2. / 3.,
        1e-10,
    );
    close(
        number(&j, "log_likelihood"),
        -1.5 * ((2. * std::f64::consts::PI * 2. / 3.).ln() + 1.),
        1e-9,
    );
    let predicted = execute(
        "gaussian_mixture_predict",
        &[vec![0.], vec![1.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert_eq!(
        matrix(&predicted, "probabilities"),
        vec![vec![1.], vec![1.]]
    );
    close(
        values(&predicted, "sample_log_likelihood")[0],
        -0.5 * (2. * std::f64::consts::PI * 2. / 3.).ln(),
        1e-9,
    );
}

#[test]
fn full_covariance_mixture_finds_separated_clusters_with_monotone_likelihood() {
    let data = [
        vec![-5.2],
        vec![-5.],
        vec![-4.8],
        vec![4.8],
        vec![5.],
        vec![5.2],
    ];
    let args = [2., 200., 1e-9, 1e-6, 42., 2.];
    let j = execute("gaussian_mixture", &data, &args).unwrap();
    assert_eq!(j, execute("gaussian_mixture", &data, &args).unwrap());
    assert!(j.contains("\"converged\":true"));
    assert_eq!(values(&j, "labels"), vec![0., 0., 0., 1., 1., 1.]);
    let means = matrix(&j, "component_means");
    close(means[0][0], -5., 1e-7);
    close(means[1][0], 5., 1e-7);
    for pair in values(&j, "log_likelihood_history").windows(2) {
        assert!(pair[1] + 1e-9 >= pair[0]);
    }
    let predicted = execute(
        "gaussian_mixture_predict",
        &[vec![-5.], vec![5.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert_eq!(values(&predicted, "labels"), vec![0., 1.]);
}

#[test]
fn batch_som_one_node_is_mean_and_two_nodes_recover_separated_codebook() {
    let j = execute(
        "self_organizing_map",
        &[vec![0.], vec![2.], vec![4.]],
        &[1., 1.],
    )
    .unwrap();
    close(matrix(&j, "codebook")[0][0], 2., 1e-10);
    close(
        number(&j, "quantization_error"),
        4. / (3. * (8_f64 / 3.).sqrt()),
        1e-9,
    );
    let two = execute(
        "self_organizing_map",
        &[vec![-2.], vec![-1.], vec![1.], vec![2.]],
        &[1., 2., 50., 0.001, 1., 1e-10],
    )
    .unwrap();
    let codebook = matrix(&two, "codebook");
    close(codebook[0][0], -1.5, 1e-9);
    close(codebook[1][0], 1.5, 1e-9);
    assert_eq!(values(&two, "labels"), vec![0., 0., 1., 1.]);
    let predicted = execute(
        "som_predict",
        &[vec![-1.5], vec![1.5]],
        &values(&two, "predict_parameters"),
    )
    .unwrap();
    assert_eq!(values(&predicted, "labels"), vec![0., 1.]);
}

#[test]
fn latent_models_keep_small_feature_spread_at_large_origins() {
    let base = [vec![0.], vec![2.], vec![4.]];
    let shifted = base.iter().map(|r| vec![r[0] + 1e16]).collect::<Vec<_>>();
    let a = execute("gaussian_mixture", &base, &[1.]).unwrap();
    let b = execute("gaussian_mixture", &shifted, &[1.]).unwrap();
    close(
        number(&a, "log_likelihood"),
        number(&b, "log_likelihood"),
        1e-9,
    );
    let p = execute(
        "gaussian_mixture_predict",
        &shifted,
        &values(&b, "predict_parameters"),
    )
    .unwrap();
    check_array(
        &values(&p, "sample_log_likelihood"),
        &values(&a, "sample_log_likelihood"),
        1e-9,
    );
    let embedded = execute("mds", &shifted, &[1., 0.]).unwrap();
    let x = matrix(&embedded, "coordinates");
    close((x[0][0] - x[2][0]).abs(), 4., 1e-9);
}

#[test]
fn latent_validation_rejects_bad_domains_singular_inputs_and_invalid_models() {
    for operation in [
        "mds",
        "correspondence_analysis",
        "factor_analysis",
        "canonical_correlation",
        "gaussian_mixture",
        "self_organizing_map",
        "gaussian_mixture_predict",
        "som_predict",
    ] {
        assert!(execute(operation, &[], &[]).is_err());
        assert!(execute(operation, &[vec![f64::NAN]], &[]).is_err());
        assert!(execute(operation, &[vec![1., 2.], vec![3.]], &[]).is_err());
    }
    assert!(execute("mds", &[vec![0., 1.], vec![2., 0.]], &[1., 1.]).is_err());
    assert!(execute(
        "correspondence_analysis",
        &[vec![1., -1.], vec![2., 3.]],
        &[1.]
    )
    .is_err());
    assert!(execute(
        "correspondence_analysis",
        &[vec![0., 0.], vec![2., 3.]],
        &[1.]
    )
    .is_err());
    assert!(execute(
        "factor_analysis",
        &equicorrelation(),
        &[3., 0., 100., 1e-8, 1., 1.]
    )
    .is_err());
    assert!(execute("gaussian_mixture", &[vec![0.], vec![1.]], &[3.]).is_err());
    assert!(execute("self_organizing_map", &[vec![0.], vec![1.]], &[32., 32.]).is_err());
    assert!(execute("gaussian_mixture_predict", &[vec![0.]], &[1., 1.]).is_err());
    assert!(execute("som_predict", &[vec![0.]], &[1., 1.]).is_err());
}

#[test]
fn gaussian_density_rejects_nonfinite_quadratic_forms_instead_of_center_density() {
    // Inverse crossproducts of extreme finite inputs can become inf-inf=NaN.
    // NaN.max(0) must never turn an extreme outlier into the modal density.
    let model = [
        2., 1., 0., 0., 0., 0., 1., 1., 1., 0., 0., 0.6, 0.4, 0.4, 0.6,
    ];
    assert!(execute("gaussian_mixture_predict", &[vec![1e308, 1e308]], &model).is_err());
}

#[test]
fn fit_rejects_models_exceeding_the_native_portable_parameter_limit() {
    // The native SQL/FFI boundary accepts at most 16384 model parameters.
    assert!(execute(
        "gaussian_mixture",
        &vec![vec![0.; 32]; 32],
        &[32., 1., 1e-6, 1e-6, 42.]
    )
    .is_err());
    assert!(execute(
        "self_organizing_map",
        &[vec![0.; 64], vec![1.; 64]],
        &[16., 16., 1.]
    )
    .is_err());
}

#[test]
fn rank_deficient_cross_covariance_keeps_orthonormal_null_canonical_scores() {
    let c = 0.23_f64.cos();
    let s = 0.23_f64.sin();
    let rows = (0..16)
        .map(|i| {
            let w = (0..4)
                .map(|j| if (i >> j) & 1 == 0 { -1. } else { 1. })
                .collect::<Vec<_>>();
            let y0 = 0.8 * w[0] + 0.6 * w[2];
            let y1 = w[3];
            vec![
                c * w[0] + s * w[1],
                -s * w[0] + c * w[1],
                c * y0 + s * y1,
                -s * y0 + c * y1,
            ]
        })
        .collect::<Vec<_>>();
    let j = execute("canonical_correlation", &rows, &[2., 2., 0.]).unwrap();
    check_array(&values(&j, "canonical_correlations"), &[0.8, 0.], 1e-12);
    for key in ["x_scores", "y_scores"] {
        let scores = matrix(&j, key);
        for a in 0..2 {
            for b in 0..2 {
                close(
                    scores.iter().map(|r| r[a] * r[b]).sum::<f64>() / 15.,
                    if a == b { 1. } else { 0. },
                    1e-10,
                );
            }
        }
    }
}

#[test]
fn saved_mixture_rejects_asymmetry_relative_to_small_covariance_scale() {
    let model = [
        2., 1., 0., 0., 0., 0., 1., 1., 1., 0., 0., 1e-20, 0., 1e-11, 1e-20,
    ];
    assert!(execute("gaussian_mixture_predict", &[vec![0., 0.]], &model).is_err());
}

#[test]
fn default_embedding_dimensions_respect_the_available_axes() {
    let j = execute("mds", &[vec![0.], vec![3.]], &[]).unwrap();
    let coordinates = matrix(&j, "coordinates");
    assert_eq!(coordinates[0].len(), 1);
    close((coordinates[0][0] - coordinates[1][0]).abs(), 3., 1e-12);
    let j = execute(
        "correspondence_analysis",
        &[vec![20., 10.], vec![10., 20.]],
        &[],
    )
    .unwrap();
    assert_eq!(matrix(&j, "row_coordinates")[0].len(), 1);
    close(values(&j, "singular_values")[0], 1. / 3., 1e-12);
}
