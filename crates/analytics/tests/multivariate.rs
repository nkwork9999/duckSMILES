use ducksmiles_analytics::multivariate::run;

fn numeric(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = &json[start..];
    let end = tail.find([',', '}', ']']).unwrap_or(tail.len());
    tail[..end].parse().expect("JSON numeric")
}
fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let end = json[start..].find(']').unwrap() + start;
    if json[start..end].is_empty() {
        return vec![];
    }
    json[start..end]
        .split(',')
        .map(|x| x.parse().unwrap())
        .collect()
}
fn matrix(json: &str, key: &str) -> Vec<Vec<f64>> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &json[start..];
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
fn close(a: f64, e: f64, t: f64) {
    assert!((a - e).abs() <= t * (1. + e.abs()), "{a} != {e}");
}

#[test]
fn covariance_sample_population_and_correlation_match_hand_computed_fixture() {
    // Mutation caught: wrong denominator or uncentered crossproducts.
    let rows = [vec![1., 2.], vec![2., 1.], vec![3., 6.]];
    let sample = run("covariance_matrix", &rows, &[]).unwrap();
    let cov = matrix(&sample, "covariance");
    close(cov[0][0], 1., 1e-12);
    close(cov[0][1], 2., 1e-12);
    close(cov[1][1], 7., 1e-12);
    let pop = run("covariance_matrix", &rows, &[0.]).unwrap();
    close(matrix(&pop, "covariance")[1][1], 14. / 3., 1e-12);
    let corr = run("correlation_matrix", &rows, &[]).unwrap();
    close(matrix(&corr, "correlation")[0][1], 2. / 7_f64.sqrt(), 1e-12);
    assert!(run("correlation_matrix", &[vec![1., 2.], vec![1., 3.]], &[]).is_err());
}

#[test]
fn covariance_preserves_small_spread_at_large_offsets() {
    let base = [vec![0., 0.], vec![2., 2.], vec![4., 6.]];
    let shifted = base
        .iter()
        .map(|r| vec![1e16 + r[0], 1e16 + r[1]])
        .collect::<Vec<_>>();
    let a = matrix(&run("covariance_matrix", &base, &[]).unwrap(), "covariance");
    let b = matrix(
        &run("covariance_matrix", &shifted, &[]).unwrap(),
        "covariance",
    );
    for (ar, br) in a.iter().zip(b) {
        for (av, bv) in ar.iter().zip(br) {
            close(*av, bv, 1e-12);
        }
    }
}

#[test]
fn hierarchical_linkage_matches_scipy_ward_fixture_and_returns_cut_labels() {
    let rows = [vec![0.], vec![1.], vec![4.], vec![5.]];
    for (method, last_distance) in [(0., 3.), (1., 5.), (2., 4.), (3., 32_f64.sqrt())] {
        let json = run("hierarchical_cluster", &rows, &[2., method]).unwrap();
        let merges = matrix(&json, "merges");
        assert_eq!(merges[0], vec![0., 1., 1., 2.]);
        assert_eq!(merges[1], vec![2., 3., 1., 2.]);
        close(merges[2][2], last_distance, 1e-12);
        assert_eq!(array(&json, "labels"), vec![0., 0., 1., 1.]);
    }
}

#[test]
fn mahalanobis_uses_only_reference_rows_and_correct_covariance_inverse() {
    let rows = [
        vec![-1., -1.],
        vec![-1., 1.],
        vec![1., -1.],
        vec![1., 1.],
        vec![2., 0.],
        vec![0., 2.],
    ];
    let json = run("mahalanobis", &rows, &[4., 0.]).unwrap();
    assert_eq!(array(&json, "reference_mean"), vec![0., 0.]);
    let covariance = matrix(&json, "covariance");
    close(covariance[0][0], 4. / 3., 1e-12);
    for d in array(&json, "squared_distances") {
        close(d, 3., 1e-12);
    }
    let changed = [
        rows[0].clone(),
        rows[1].clone(),
        rows[2].clone(),
        rows[3].clone(),
        vec![1e4, 1e4],
    ];
    close(
        matrix(
            &run("mahalanobis", &changed, &[4., 0.]).unwrap(),
            "covariance",
        )[0][0],
        4. / 3.,
        1e-12,
    );
    let singular = [vec![0., 0.], vec![1., 1.], vec![2., 2.], vec![1., 2.]];
    assert!(run("mahalanobis", &singular, &[3., 0.]).is_err());
    assert!(run("mahalanobis", &singular, &[3., 0.1]).is_ok());
}

#[test]
fn lda_and_qda_return_empirical_priors_posterior_and_heldout_predictions() {
    let rows = [
        vec![-2., 10.],
        vec![0., 10.],
        vec![2., 20.],
        vec![4., 20.],
        vec![-1., 999.5],
        vec![3., -99.5],
    ];
    for method in [0., 1.] {
        let json = run("discriminant_analysis", &rows, &[4., method, 0.]).unwrap();
        assert_eq!(array(&json, "classes"), vec![10., 20.]);
        assert_eq!(array(&json, "class_priors"), vec![0.5, 0.5]);
        assert_eq!(array(&json, "predictions"), vec![10., 20.]);
        let posterior = matrix(&json, "posterior_probabilities");
        close(posterior[0][0], 1. / (1. + (-4_f64).exp()), 1e-12);
        close(posterior[1][1], 1. / (1. + (-4_f64).exp()), 1e-12);
    }
}

#[test]
fn discriminant_analysis_rejects_unidentifiable_covariance_unless_regularized() {
    let rows = [
        vec![1., 1., 0.],
        vec![1., 1., 0.],
        vec![2., 2., 1.],
        vec![2., 2., 1.],
        vec![1., 2., 99.],
    ];
    for method in [0., 1.] {
        assert!(run("discriminant_analysis", &rows, &[4., method, 0.]).is_err());
        let fit = run("discriminant_analysis", &rows, &[4., method, 0.1]).unwrap();
        let posterior = matrix(&fit, "posterior_probabilities");
        close(posterior[0].iter().sum(), 1., 1e-12);
    }
}

#[test]
fn naive_bayes_uses_population_class_variance_and_ignores_query_labels() {
    let rows = [
        vec![-2., 0.],
        vec![0., 0.],
        vec![2., 1.],
        vec![4., 1.],
        vec![-1., 999.5],
        vec![3., -99.5],
    ];
    let json = run("gaussian_naive_bayes", &rows, &[4., 1e-9]).unwrap();
    assert_eq!(array(&json, "predictions"), vec![0., 1.]);
    let posterior = matrix(&json, "posterior_probabilities");
    close(posterior[0][0], 1. / (1. + (-8_f64).exp()), 1e-12);
    close(matrix(&json, "class_variances")[0][0], 1., 1e-12);
    let constant = [
        vec![1., 0.],
        vec![1., 0.],
        vec![1., 1.],
        vec![1., 1.],
        vec![1., -999.],
    ];
    let result = run("gaussian_naive_bayes", &constant, &[4., 0.01]).unwrap();
    close(matrix(&result, "posterior_probabilities")[0][0], 0.5, 1e-12);
}

#[test]
fn nearest_neighbors_scale_only_training_rows_and_use_deterministic_ties() {
    let rows = [
        vec![0., 0., 0.],
        vec![2., 100., 10.],
        vec![4., 0., 20.],
        vec![1., 100., 999.],
    ];
    let json = run("knn_predict", &rows, &[3., 2., 0.]).unwrap();
    // Training scales are sqrt(8/3), sqrt(20000/9); closest is row1 then row0.
    close(array(&json, "predictions")[0], 5., 1e-12);
    let tie = [vec![-1., 9.], vec![1., 3.], vec![0., 100.5]];
    assert_eq!(
        array(
            &run("knn_predict", &tie, &[2., 2., 1.]).unwrap(),
            "predictions"
        ),
        vec![3.]
    );
    assert_eq!(
        array(
            &run("knn_predict", &tie, &[2., 1., 1.]).unwrap(),
            "predictions"
        ),
        vec![9.]
    );
}

#[test]
fn dbscan_returns_density_clusters_noise_and_core_membership() {
    let rows = [
        vec![0.],
        vec![0.1],
        vec![0.2],
        vec![2.],
        vec![2.1],
        vec![2.2],
        vec![9.],
    ];
    let json = run("dbscan", &rows, &[0.21, 3.]).unwrap();
    assert_eq!(array(&json, "labels"), vec![0., 0., 0., 1., 1., 1., -1.]);
    assert_eq!(array(&json, "core_indices"), vec![0., 1., 2., 3., 4., 5.]);
    close(numeric(&json, "n_clusters"), 2., 0.);
    let border = [vec![0.], vec![0.1], vec![0.2], vec![0.4], vec![9.]];
    let fit = run("dbscan", &border, &[0.21, 3.]).unwrap();
    assert_eq!(array(&fit, "labels"), vec![0., 0., 0., 0., -1.]);
    assert!(!array(&fit, "core_indices").contains(&3.));
}

#[test]
fn multivariate_validation_rejects_nonfinite_ragged_impossible_parameters_and_work() {
    for operation in [
        "covariance_matrix",
        "correlation_matrix",
        "hierarchical_cluster",
        "mahalanobis",
        "discriminant_analysis",
        "gaussian_naive_bayes",
        "knn_predict",
        "dbscan",
    ] {
        assert!(run(operation, &[], &[]).is_err());
        assert!(run(operation, &[vec![f64::NAN]], &[]).is_err());
        assert!(run(operation, &[vec![1., 2.], vec![1.]], &[]).is_err());
    }
    assert!(run("covariance_matrix", &[vec![1.]], &[]).is_err());
    assert!(run("hierarchical_cluster", &[vec![0.], vec![1.]], &[3., 0.]).is_err());
    assert!(run("hierarchical_cluster", &[vec![0.], vec![1.]], &[1., 4.]).is_err());
    assert!(run("dbscan", &[vec![0.]], &[-1., 1.]).is_err());
    assert!(run("dbscan", &[vec![0.]], &[1., 1.5]).is_err());
    assert!(run("mahalanobis", &[vec![0.], vec![1.], vec![2.]], &[2., -1.]).is_err());
    assert!(run("knn_predict", &[vec![0., 1.], vec![1., 2.]], &[1., 2., 0.]).is_err());
    assert!(run("knn_predict", &[vec![0., 1.5], vec![1., 2.]], &[1., 1., 1.]).is_err());
    assert!(run(
        "gaussian_naive_bayes",
        &[vec![0., 0.], vec![1., 1.], vec![2., 0.]],
        &[2., 0.]
    )
    .is_err());
}

#[test]
fn unequal_size_linkage_updates_match_independent_scipy_fixture() {
    let rows = [
        vec![0., 0.],
        vec![2., 0.],
        vec![0., 1.],
        vec![3., 4.],
        vec![4., 3.],
    ];
    for (method, third, last) in [
        (0., 2., 3.605551275463989),
        (1., 2.23606797749979, 5.),
        (2., 2.118033988749895, 4.407238923866752),
        (3., 2.3804761428476167, 6.5828058860438325),
    ] {
        let fit = run("hierarchical_cluster", &rows, &[2., method]).unwrap();
        let merges = matrix(&fit, "merges");
        close(merges[2][2], third, 1e-12);
        close(merges[3][2], last, 1e-12);
        assert_eq!(array(&fit, "labels"), vec![0., 0., 0., 1., 1.]);
    }
}

#[test]
fn correlated_mahalanobis_matches_independent_numpy_inverse_fixture() {
    let rows = [
        vec![-2., 1., 0.],
        vec![-1., 0., 1.],
        vec![0., 1., 2.],
        vec![1., -1., 1.],
        vec![2., 2., 0.],
        vec![3., 0., 1.],
        vec![4., 1., 1.],
        vec![1., 3., 2.],
        vec![-1., 0.5, 1.],
        vec![2., 2., 2.],
    ];
    let fit = run("mahalanobis", &rows, &[8., 0.1]).unwrap();
    let distances = array(&fit, "squared_distances");
    close(distances[0], 1.0459972920093779, 1e-12);
    close(distances[1], 2.131102652059181, 1e-12);
}

#[test]
fn density_cluster_membership_is_permutation_equivariant_without_shared_borders() {
    let rows = [
        vec![0.],
        vec![0.1],
        vec![0.2],
        vec![2.],
        vec![2.1],
        vec![2.2],
        vec![9.],
    ];
    let order = [6, 4, 1, 3, 0, 5, 2];
    let reordered = order.iter().map(|i| rows[*i].clone()).collect::<Vec<_>>();
    let a = array(&run("dbscan", &rows, &[0.21, 3.]).unwrap(), "labels");
    let b = array(&run("dbscan", &reordered, &[0.21, 3.]).unwrap(), "labels");
    for i in 0..order.len() {
        assert_eq!(a[order[i]] < 0., b[i] < 0.);
        for j in 0..order.len() {
            assert_eq!(a[order[i]] == a[order[j]], b[i] == b[j]);
        }
    }
}

#[test]
fn classifiers_and_neighbors_preserve_small_feature_spreads_at_large_offsets() {
    let rows = [
        vec![0., 0.],
        vec![2., 0.],
        vec![6., 1.],
        vec![8., 1.],
        vec![2., -1.],
        vec![6., -1.],
    ];
    let shifted = rows
        .iter()
        .map(|r| vec![r[0] + 1e16, r[1]])
        .collect::<Vec<_>>();
    for (operation, args) in [
        ("discriminant_analysis", vec![4., 0., 0.2]),
        ("discriminant_analysis", vec![4., 1., 0.2]),
        ("gaussian_naive_bayes", vec![4., 1e-9]),
        ("knn_predict", vec![4., 1., 1.]),
    ] {
        let base = run(operation, &rows, &args).unwrap();
        let fit = run(operation, &shifted, &args).unwrap();
        assert_eq!(array(&fit, "predictions"), array(&base, "predictions"));
        if operation != "knn_predict" {
            for (a, b) in matrix(&base, "posterior_probabilities")
                .iter()
                .zip(matrix(&fit, "posterior_probabilities"))
            {
                for (x, y) in a.iter().zip(b) {
                    close(*x, y, 1e-12);
                }
            }
        }
    }
}
