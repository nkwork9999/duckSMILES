#[path = "../../src/functional_preprocess.rs"]
mod functional_preprocess;
use functional_preprocess::run as execute;
fn number(j: &str, key: &str) -> f64 {
    let s = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    j[s..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn values(j: &str, key: &str) -> Vec<f64> {
    let s = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let e = j[s..].find(']').unwrap() + s;
    if e == s {
        vec![]
    } else {
        j[s..e].split(',').map(|s| s.parse().unwrap()).collect()
    }
}
fn matrix(j: &str, key: &str) -> Vec<Vec<f64>> {
    let s = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let e = j[s..].find("]]").unwrap() + s + 1;
    j[s..e]
        .split("],")
        .map(|s| {
            s.trim_matches(['[', ']'])
                .split(',')
                .map(|s| s.parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
#[test]
fn irregular_curves_interpolate_on_training_only_common_grid() {
    let data = vec![
        vec![4., 8., 1.],
        vec![0., 0., 1.],
        vec![2., 4., 1.],
        vec![1., 3., 2.],
        vec![2., 6., 2.],
        vec![3., 9., 2.],
        vec![4., 12., 2.],
    ];
    let j = execute("functional_grid", &data, &[4.]).unwrap();
    assert_eq!(values(&j, "grid"), vec![1., 2., 3., 4.]);
    assert_eq!(
        matrix(&j, "values"),
        vec![vec![2., 4., 6., 8.], vec![3., 6., 9., 12.]]
    );
}
#[test]
fn explicit_bins_average_observations_and_keep_empty_bins_null() {
    let data = vec![
        vec![0., 1., 1.],
        vec![0.25, 3., 1.],
        vec![1., 5., 1.],
        vec![2., 9., 1.],
        vec![0., 2., 2.],
        vec![2., 10., 2.],
    ];
    let j = execute("functional_bin", &data, &[3., 0., 2.]).unwrap();
    assert!(j.contains("\"values\":[[2,5,9],[2,null,10]]"), "{j}");
    assert_eq!(
        matrix(&j, "counts"),
        vec![vec![2., 1., 1.], vec![1., 0., 1.]]
    );
}
#[test]
fn thinning_uses_sorted_per_curve_stride_and_returns_original_rows() {
    let data = (0..6)
        .rev()
        .map(|x| vec![x as f64, x as f64, 7.])
        .collect::<Vec<_>>();
    let j = execute("functional_thin", &data, &[2., 1., 0.]).unwrap();
    assert_eq!(matrix(&j, "x"), vec![vec![1., 3., 5.]]);
    assert_eq!(matrix(&j, "original_row_indices"), vec![vec![4., 2., 0.]]);
}
#[test]
fn snip_clips_peaks_in_transformed_space_and_leaves_validation_untouched() {
    let mut data = vec![];
    for (i, y) in [1., 1., 9., 1., 1.].iter().enumerate() {
        data.push(vec![i as f64, *y, 7., 0.]);
        data.push(vec![i as f64, *y, 8., 1.]);
    }
    let j = execute("functional_snip", &data, &[1., 1., 0.]).unwrap();
    let b = matrix(&j, "baselines");
    let c = matrix(&j, "corrected");
    for x in &b[0] {
        close(*x, 1., 1e-12);
    }
    close(c[0][2], 8., 1e-12);
    assert_eq!(c[1], vec![1., 1., 9., 1., 1.]);
    assert!(j.contains("\"baseline_applied\":[true,false]"));
}
#[test]
fn arwls_solves_second_difference_system_and_uses_logistic_reweighting() {
    let data = [1., 1.2, 8.8, 1.1, 1., 1.3, 1.]
        .iter()
        .enumerate()
        .map(|(i, y)| vec![i as f64, *y, 7.])
        .collect::<Vec<_>>();
    let j = execute("functional_arwls", &data, &[1., 1e-12, 2.]).unwrap();
    let want = [
        1.05390881386524,
        1.1224605428364398,
        1.137111422565701,
        1.1214999108577768,
        1.1135784906789679,
        1.1297997586496344,
        1.073041429558218,
    ];
    let b = matrix(&j, "baselines");
    for (a, b) in b[0].iter().zip(want) {
        close(*a, b, 1e-10);
    }
    assert!(j.contains("\"converged\":[false]"), "{j}");
    close(matrix(&j, "weights")[0][2], 0.0018679672837341372, 1e-10);
}
#[test]
fn supplied_baseline_is_interpolated_and_excluded_from_corrected_curves() {
    let data = vec![
        vec![0., 1., 9., 0.],
        vec![2., 3., 9., 0.],
        vec![0., 4., 1., 0.],
        vec![1., 7., 1., 0.],
        vec![2., 9., 1., 0.],
        vec![0., 8., 2., 1.],
        vec![2., 12., 2., 1.],
    ];
    let j = execute("functional_baseline", &data, &[9., 2.]).unwrap();
    assert_eq!(values(&j, "ids"), vec![1., 2.]);
    assert_eq!(
        matrix(&j, "corrected"),
        vec![vec![3., 5., 6.], vec![8., 12.]]
    );
}
#[test]
fn functional_dtw_distance_matrix_handles_unequal_lengths_and_is_symmetric() {
    let mut data = vec![];
    for (id, ys) in [
        (1., vec![1., 2., 3.]),
        (2., vec![1., 1., 2., 3.]),
        (3., vec![2., 2., 2.]),
    ] {
        for (i, y) in ys.iter().enumerate() {
            data.push(vec![i as f64, *y, id]);
        }
    }
    let j = execute("functional_dtw_matrix", &data, &[2., 0., 0.]).unwrap();
    let d = matrix(&j, "distances");
    close(d[0][1], 0., 1e-12);
    close(d[0][2], 2.0_f64.sqrt(), 1e-12);
    for (i, row) in d.iter().enumerate() {
        for (k, value) in row.iter().enumerate() {
            close(*value, d[k][i], 1e-12);
        }
    }
}
#[test]
fn validation_grid_extent_and_length_cannot_change_the_training_grid() {
    let mut data = (0..3)
        .map(|i| vec![i as f64, 2. * i as f64, 7., 0.])
        .collect::<Vec<_>>();
    data.extend((0..5).map(|i| vec![100. + i as f64, 1e100 * i as f64, 8., 1.]));
    let j = execute("functional_grid", &data, &[]).unwrap();
    assert_eq!(values(&j, "grid"), vec![0., 1., 2.]);
    assert!(j.contains("\"values\":[[0,2,4],[null,null,null]]"), "{j}");
}
#[test]
fn arwls_preserves_affine_curves_and_passes_validation_without_fitting() {
    let data = vec![
        vec![0., 2., 7., 0.],
        vec![1., 5., 7., 0.],
        vec![2., 8., 7., 0.],
        vec![0., -100., 8., 1.],
        vec![0.2, 20., 8., 1.],
    ];
    let j = execute("functional_arwls", &data, &[]).unwrap();
    assert_eq!(matrix(&j, "baselines")[0], vec![2., 5., 8.]);
    assert_eq!(matrix(&j, "corrected")[0], vec![0., 0., 0.]);
    assert_eq!(matrix(&j, "corrected")[1], vec![-100., 20.]);
    assert!(j.contains("\"converged\":[true,null]"), "{j}");
    assert!(j.contains("\"baseline_applied\":[true,false]"));
}
#[test]
fn three_point_arwls_uses_the_correct_second_difference_endpoint_penalty() {
    let data = vec![vec![0., 0., 7.], vec![1., 1., 7.], vec![2., 0., 7.]];
    let j = execute("functional_arwls", &data, &[1., 1e-12, 1.]).unwrap();
    let b = matrix(&j, "baselines");
    for (a, b) in b[0].iter().zip([2. / 7., 3. / 7., 2. / 7.]) {
        close(*a, b, 1e-12);
    }
    assert!(j.contains("\"converged\":[false]"), "{j}");
}
#[test]
fn bin_edges_must_have_representable_increasing_spacing() {
    let data = vec![vec![1e16, 1., 7.], vec![1e16 + 2., 2., 7.]];
    assert!(execute("functional_bin", &data, &[4., 1e16, 1e16 + 2.]).is_err());
}
#[test]
fn explicit_saved_grid_replays_on_only_new_validation_functions_without_refitting() {
    let train = vec![
        vec![0., 0., 7., 0.],
        vec![1., 2., 7., 0.],
        vec![2., 4., 7., 0.],
    ];
    let fitted = execute("functional_grid", &train, &[]).unwrap();
    let grid = values(&fitted, "grid");
    let validation = vec![vec![0., 1., 8., 1.], vec![2., 5., 8., 1.]];
    let j = execute(
        "functional_grid",
        &validation,
        &[grid.len() as f64, grid[0], grid[2], 2.],
    )
    .unwrap();
    assert_eq!(values(&j, "grid"), grid);
    assert_eq!(matrix(&j, "values"), vec![vec![1., 3., 5.]]);
    assert!(
        j.contains("\"grid_estimation\":\"supplied_bounds_and_size\""),
        "{j}"
    );
}
#[test]
fn dtw_preserves_identical_finite_curves_even_if_their_full_intensity_span_overflows() {
    let data = vec![
        vec![0., 1e308, 7.],
        vec![1., -1e308, 7.],
        vec![0., 1e308, 8.],
        vec![1., -1e308, 8.],
    ];
    let j = execute("functional_dtw_matrix", &data, &[]).unwrap();
    assert_eq!(matrix(&j, "distances"), vec![vec![0., 0.], vec![0., 0.]]);
    let j = execute("functional_golden_match", &data, &[7., 2.]).unwrap();
    assert_eq!(values(&j, "scores"), vec![0.]);
}
#[test]
fn population_normalization_does_not_collapse_wide_finite_curves_to_zero() {
    let a = [0., 1e308, -1e308, 1e308, -1e308];
    let b = [0., 1e308, 1e308, -1e308, -1e308];
    let mut data = vec![];
    for (i, y) in a.iter().enumerate() {
        data.push(vec![i as f64, *y, 7.]);
    }
    for (i, y) in b.iter().enumerate() {
        data.push(vec![i as f64, *y, 8.]);
    }
    let j = execute("functional_dtw_matrix", &data, &[0., 1., 0.]).unwrap();
    close(matrix(&j, "distances")[0][1], 10.0_f64.sqrt(), 1e-12);
}
#[test]
fn golden_rms_stays_finite_when_the_unnormalized_euclidean_norm_overflows() {
    let mut data = vec![];
    for i in 0..5 {
        data.push(vec![i as f64, 0., 7.]);
    }
    for (i, y) in [0., 1e308, -1e308, 1e308, -1e308].iter().enumerate() {
        data.push(vec![i as f64, *y, 8.]);
    }
    let j = execute("functional_golden_match", &data, &[7., 0.]).unwrap();
    let expected = 1e308 * 0.8_f64.sqrt();
    assert!((values(&j, "scores")[0] / expected - 1.).abs() < 1e-12);
}
#[test]
fn malformed_function_roles_arguments_and_workloads_reject_without_panics() {
    for op in [
        "functional_grid",
        "functional_bin",
        "functional_thin",
        "functional_snip",
        "functional_arwls",
        "functional_baseline",
        "functional_dtw_matrix",
        "functional_golden_match",
    ] {
        assert!(execute(op, &[], &[]).is_err());
        assert!(execute(op, &[vec![0., f64::NAN, 1.]], &[]).is_err());
        assert!(execute(op, &[vec![0., 1., 1.], vec![1., 2.]], &[]).is_err());
        assert!(execute(op, &[vec![0., 1., 1.]], &[f64::INFINITY]).is_err());
        assert!(execute(op, &[vec![0., 1., 1., 0.], vec![1., 2., 1., 1.]], &[]).is_err());
    }
    let many = (0..129)
        .map(|id| vec![0., 1., id as f64])
        .collect::<Vec<_>>();
    assert!(execute("functional_thin", &many, &[]).is_err());
    let data = vec![
        vec![0., 1., 7., 1.],
        vec![1., 2., 7., 1.],
        vec![0., 3., 8., 0.],
    ];
    assert!(execute("functional_baseline", &data, &[7.]).is_err());
    let large = (0..2200)
        .flat_map(|i| vec![vec![i as f64, 0., 7.], vec![i as f64, 1., 8.]])
        .collect::<Vec<_>>();
    assert!(execute("functional_dtw_matrix", &large, &[]).is_err());
}
#[test]
fn golden_curve_matching_excludes_the_reference_and_ranks_only_candidates() {
    let mut data = vec![];
    for (id, ys) in [
        (9., vec![1., 2., 3.]),
        (1., vec![2., 3., 4.]),
        (2., vec![1., 2., 3.]),
    ] {
        for (i, y) in ys.iter().enumerate() {
            data.push(vec![i as f64, *y, id]);
        }
    }
    let j = execute("functional_golden_match", &data, &[9., 0., -1., 0.]).unwrap();
    assert_eq!(values(&j, "candidate_ids"), vec![1., 2.]);
    let scores = values(&j, "scores");
    close(scores[0], 1., 1e-12);
    close(scores[1], 0., 1e-12);
    assert_eq!(values(&j, "ranking_ids"), vec![2., 1.]);
    assert!(j.contains("\"model_fitted\":false"));
    close(number(&j, "target_id"), 9., 0.);
}
