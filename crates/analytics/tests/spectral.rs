use ducksmiles_analytics::execute as run;

fn numbers(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = &j[start..];
    let end = if tail.starts_with('[') {
        let mut depth = 0;
        tail.char_indices()
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
        tail.find([',', '}']).unwrap_or(tail.len())
    };
    tail[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1.), "{a} != {b}");
}
fn close_all(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        close(x, y);
    }
}

#[test]
fn snv_normalizes_each_spectrum_with_the_requested_variance_denominator() {
    let j = run(
        "spectral_snv",
        &[vec![1., 2., 3.], vec![10., 12., 14.]],
        &[],
    )
    .unwrap();
    close_all(&numbers(&j, "rows"), &[-1., 0., 1., -1., 0., 1.]);
    let j = run("spectral_snv", &[vec![1., 2., 3.]], &[0.]).unwrap();
    let q = 1.5_f64.sqrt();
    close_all(&numbers(&j, "rows"), &[-q, 0., q]);
}
#[test]
fn msc_reuses_the_supplied_reference_without_fitting_query_batch_statistics() {
    let reference = vec![1., 2., 4., 8.];
    let rows = vec![
        reference.iter().map(|v| 3. + 2. * v).collect(),
        reference.iter().map(|v| -1. + 0.5 * v).collect(),
    ];
    let j = run("spectral_msc", &rows, &reference).unwrap();
    close_all(&numbers(&j, "rows"), &[1., 2., 4., 8., 1., 2., 4., 8.]);
    close_all(&numbers(&j, "slopes"), &[2., 0.5]);
    close_all(&numbers(&j, "intercepts"), &[3., -1.]);
    let one = run("spectral_msc", &[rows[0].clone()], &reference).unwrap();
    close_all(&numbers(&one, "rows"), &reference);
    let fitted = run("spectral_msc", &rows, &[]).unwrap();
    close_all(
        &numbers(&fitted, "reference_spectrum"),
        &[2.25, 3.5, 6., 11.],
    );
}
#[test]
fn savitzky_golay_preserves_polynomials_and_analytical_derivatives_at_edges() {
    let x = vec![vec![0., 1., 4., 9., 16., 25., 36.]];
    close_all(
        &numbers(&run("savitzky_golay", &x, &[5., 2.]).unwrap(), "rows"),
        &x[0],
    );
    close_all(
        &numbers(
            &run("savitzky_golay", &x, &[5., 2., 1., 0.5]).unwrap(),
            "rows",
        ),
        &[0., 4., 8., 12., 16., 20., 24.],
    );
    close_all(
        &numbers(
            &run("savitzky_golay", &x, &[5., 2., 2., 0.5]).unwrap(),
            "rows",
        ),
        &[8.; 7],
    );
}
#[test]
fn savitzky_golay_matches_independent_scipy_interp_boundary_fits() {
    let j = run(
        "savitzky_golay",
        &[vec![2., 4., 3., 5., 7., 6., 8., 9.]],
        &[5., 2.],
    )
    .unwrap();
    close_all(
        &numbers(&j, "rows"),
        &[
            2.4285714285714297,
            2.8857142857142852,
            3.77142857142857,
            5.,
            6.228571428571425,
            6.857142857142854,
            7.828571428571427,
            8.94285714285714,
        ],
    );
}
#[test]
fn finite_derivative_uses_second_order_boundary_stencils_and_correct_spacing() {
    let x = vec![vec![0., 1., 4., 9., 16.]];
    close_all(
        &numbers(&run("curve_derivative", &x, &[1., 0.5]).unwrap(), "rows"),
        &[0., 4., 8., 12., 16.],
    );
    close_all(
        &numbers(&run("curve_derivative", &x, &[2., 0.5]).unwrap(), "rows"),
        &[8.; 5],
    );
}
#[test]
fn polynomial_detrending_matches_numpy_least_squares_in_the_documented_grid_basis() {
    let j = run(
        "spectral_detrend",
        &[vec![2., 4., 3., 5., 7., 6., 8., 9.]],
        &[1.],
    )
    .unwrap();
    close_all(&numbers(&j, "coefficients"), &[5.5, 10. / 3.]);
    close_all(
        &numbers(&j, "rows"),
        &[
            -1. / 6.,
            37. / 42.,
            -15. / 14.,
            -1. / 42.,
            43. / 42.,
            -13. / 14.,
            5. / 42.,
            1. / 6.,
        ],
    );
    let j = run("spectral_detrend", &[vec![2., 3., 6., 11., 18.]], &[2.]).unwrap();
    close_all(&numbers(&j, "rows"), &[0.; 5]);
}
#[test]
fn haar_uses_orthonormal_mallat_order_and_reverses_both_full_and_partial_levels() {
    let x = vec![vec![1., 2., 3., 4.]];
    let q = std::f64::consts::FRAC_1_SQRT_2;
    let j = run("haar_transform", &x, &[]).unwrap();
    let coeff = numbers(&j, "rows");
    close_all(&coeff, &[5., -2., -q, -q]);
    close(coeff.iter().map(|v| v * v).sum(), 30.);
    close_all(
        &numbers(&run("haar_transform", &[coeff], &[1.]).unwrap(), "rows"),
        &x[0],
    );
    let partial = numbers(&run("haar_transform", &x, &[0., 1.]).unwrap(), "rows");
    close_all(
        &numbers(
            &run("haar_transform", &[partial], &[1., 1.]).unwrap(),
            "rows",
        ),
        &x[0],
    );
}
#[test]
fn dft_fourier_coefficients_and_complex_inverse_match_numpy_fft_conventions() {
    let j = run("fourier_transform", &[vec![1., 2., 3., 4.]], &[]).unwrap();
    let coeff = numbers(&j, "rows");
    close_all(&coeff, &[10., -2., -2., -2., 0., 2., 0., -2.]);
    let inverse = run("fourier_transform", &[coeff], &[1.]).unwrap();
    close_all(&numbers(&inverse, "rows"), &[1., 2., 3., 4.]);
    close_all(&numbers(&inverse, "imaginary_rows"), &[0.; 4]);
    let complex = run("fourier_transform", &[vec![0., 0., 2., 0.]], &[1.]).unwrap();
    close_all(&numbers(&complex, "imaginary_rows"), &[1., 1.]);
}
#[test]
fn dtw_accepts_unequal_series_and_returns_a_monotone_endpoint_to_endpoint_path() {
    let data = vec![
        vec![1., 0.],
        vec![2., 0.],
        vec![3., 0.],
        vec![1., 1.],
        vec![1., 1.],
        vec![2., 1.],
        vec![3., 1.],
    ];
    let j = run("dynamic_time_warping", &data, &[]).unwrap();
    close(numbers(&j, "distance")[0], 0.);
    let path = numbers(&j, "path");
    assert_eq!(&path[..2], &[0., 0.]);
    assert_eq!(&path[path.len() - 2..], &[2., 3.]);
    for adjacent in path.chunks(2).collect::<Vec<_>>().windows(2) {
        let da = adjacent[1][0] - adjacent[0][0];
        let db = adjacent[1][1] - adjacent[0][1];
        assert!((0. ..=1.).contains(&da) && (0. ..=1.).contains(&db) && da + db >= 1.);
    }
}
#[test]
fn functional_pca_uses_trapezoid_l2_geometry_and_reconstructs_a_known_rank_one_curve_family() {
    let data = vec![vec![1., 2., 3.], vec![2., 4., 6.], vec![3., 6., 9.]];
    let j = run("functional_pca", &data, &[1.]).unwrap();
    close_all(&numbers(&j, "eigenvalues"), &[9.]);
    close_all(&numbers(&j, "quadrature_weights"), &[0.5, 1., 0.5]);
    close_all(&numbers(&j, "eigenfunctions"), &[1. / 3., 2. / 3., 1.]);
    close_all(&numbers(&j, "scores"), &[-3., 0., 3.]);
    close_all(
        &numbers(&j, "reconstructed_rows"),
        &data.into_iter().flatten().collect::<Vec<_>>(),
    );
    let scaled = run(
        "functional_pca",
        &[vec![1., 2., 3.], vec![2., 4., 6.], vec![3., 6., 9.]],
        &[1., 2.],
    )
    .unwrap();
    close(numbers(&scaled, "eigenvalues")[0], 18.);
}
#[test]
fn peak_finding_reports_global_prominence_bases_plateaus_and_distance_priority() {
    let j = run("peak_find", &[vec![0., 3., 1., 4., 0., 2., 0.]], &[]).unwrap();
    assert_eq!(numbers(&j, "peak_indices"), vec![1., 3., 5.]);
    close_all(&numbers(&j, "prominences"), &[2., 4., 2.]);
    assert_eq!(numbers(&j, "left_bases"), vec![0., 0., 4.]);
    assert_eq!(numbers(&j, "right_bases"), vec![2., 4., 6.]);
    let j = run("peak_find", &[vec![0., 2., 2., 0., 3., 0.]], &[0., 4.]).unwrap();
    assert_eq!(numbers(&j, "peak_indices"), vec![4.]);
    assert!(numbers(
        &run("peak_find", &[vec![2.; 8]], &[]).unwrap(),
        "peak_indices"
    )
    .is_empty());
}
#[test]
fn preprocessing_preserves_small_variations_at_large_origins_and_location_invariance() {
    let x = vec![vec![1e16, 1e16 + 2., 1e16 + 4., 1e16 + 6., 1e16 + 8.]];
    close_all(
        &numbers(&run("curve_derivative", &x, &[]).unwrap(), "rows"),
        &[2.; 5],
    );
    close_all(
        &numbers(&run("savitzky_golay", &x, &[5., 2., 1.]).unwrap(), "rows"),
        &[2.; 5],
    );
    let baseline = numbers(
        &run("spectral_snv", &[vec![0., 2., 4., 6., 8.]], &[]).unwrap(),
        "rows",
    );
    close_all(
        &numbers(&run("spectral_snv", &x, &[]).unwrap(), "rows"),
        &baseline,
    );
}
#[test]
fn spectral_contract_rejects_undefined_models_invalid_shapes_and_excessive_work() {
    for op in [
        "spectral_snv",
        "spectral_msc",
        "savitzky_golay",
        "curve_derivative",
        "spectral_detrend",
        "haar_transform",
        "fourier_transform",
        "dynamic_time_warping",
        "functional_pca",
        "peak_find",
    ] {
        for rows in [
            vec![],
            vec![vec![]],
            vec![vec![f64::NAN]],
            vec![vec![1., 2.], vec![3.]],
        ] {
            assert!(run(op, &rows, &[]).is_err(), "{op}");
        }
        assert!(
            run(op, &[vec![1., 2., 3.]], &[f64::INFINITY]).is_err(),
            "{op}"
        );
    }
    for (op, data, args) in [
        ("spectral_snv", vec![vec![1.; 4]], vec![]),
        ("spectral_msc", vec![vec![1.; 4]], vec![]),
        ("savitzky_golay", vec![vec![0.; 7]], vec![4., 2.]),
        ("savitzky_golay", vec![vec![0.; 7]], vec![5., 5.]),
        ("savitzky_golay", vec![vec![0.; 7]], vec![5., 2., 3.]),
        ("curve_derivative", vec![vec![1., 2.]], vec![]),
        ("spectral_detrend", vec![vec![1., 2., 3.]], vec![3.]),
        ("haar_transform", vec![vec![1., 2., 3.]], vec![]),
        ("fourier_transform", vec![vec![1., 2., 3.]], vec![1.]),
        (
            "dynamic_time_warping",
            vec![vec![1., 2.], vec![2., 0.]],
            vec![],
        ),
        ("functional_pca", vec![vec![1., 2.], vec![1., 2.]], vec![1.]),
        ("peak_find", vec![vec![0., 1., 0.]], vec![-1.]),
    ] {
        assert!(run(op, &data, &args).is_err(), "{op}");
    }
}
