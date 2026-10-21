use ducksmiles_analytics::design::run;

// These tests catch missing corner/axial points, duplicates or weak information in
// designs, wrong quadratic basis ordering, and accepting unidentified models.
fn numbers(json: &str, key: &str) -> Vec<f64> {
    let needle = format!("\"{key}\":");
    let rest = json
        .split_once(&needle)
        .unwrap_or_else(|| panic!("missing {key}: {json}"))
        .1
        .trim_start();
    let end = if rest.starts_with('[') {
        let mut depth = 0;
        rest.char_indices()
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
        rest.find([',', '}']).unwrap_or(rest.len())
    };
    rest[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn factorial_covers_all_corners_and_requested_centers() {
    let output = run("factorial_design", &[vec![10., 20.], vec![-2., 2.]], &[2.]).unwrap();
    let flat = numbers(&output, "rows");
    assert_eq!(flat.len(), 12);
    for expected in [[10., -2.], [10., 2.], [20., -2.], [20., 2.]] {
        assert!(flat.chunks(2).any(|r| r == expected));
    }
    assert_eq!(flat.chunks(2).filter(|r| *r == [15., 0.]).count(), 2);
}

#[test]
fn central_composite_places_axial_runs_in_original_units() {
    let output = run(
        "central_composite_design",
        &[vec![10., 20.], vec![-2., 2.]],
        &[2., 1.],
    )
    .unwrap();
    let flat = numbers(&output, "rows");
    assert_eq!(flat.len(), 18);
    for expected in [[5., 0.], [25., 0.], [15., -4.], [15., 4.], [15., 0.]] {
        assert!(
            flat.chunks(2).any(|r| r == expected),
            "missing {expected:?}"
        );
    }
}

#[test]
fn d_optimal_selects_distinct_constrained_corner_rows() {
    let candidates = vec![
        vec![-1., -1.],
        vec![-1., 1.],
        vec![1., -1.],
        vec![1., 1.],
        vec![0., 0.],
    ];
    let output = run("d_optimal_design", &candidates, &[4., 30., 42.]).unwrap();
    let ids = numbers(&output, "indices");
    assert_eq!(ids, vec![0., 1., 2., 3.]);
    // Four balanced corners give X'X = diag(4,4,4), determinant 64.
    near(numbers(&output, "log_determinant")[0], 64_f64.ln(), 1e-10);
    assert_eq!(
        run("d_optimal_design", &candidates, &[4., 30., 42.]).unwrap(),
        output
    );
}

#[test]
fn response_surface_recovers_a_known_quadratic_and_predicts_new_points() {
    let mut data = Vec::new();
    for x in [-1., 0., 1.] {
        for y in [-1., 0., 1.] {
            data.push(vec![
                x,
                y,
                3. + 2. * x - 4. * y + 5. * x * x + 6. * y * y + 7. * x * y,
            ]);
        }
    }
    let model = run("response_surface", &data, &[0.]).unwrap();
    for (a, b) in numbers(&model, "coefficients")
        .into_iter()
        .zip([3., 2., -4., 5., 6., 7.])
    {
        near(a, b, 1e-10);
    }
    near(numbers(&model, "rmse")[0], 0., 1e-10);
    near(numbers(&model, "r2")[0], 1., 1e-10);
    let predictions = run(
        "response_surface_predict",
        &[vec![2., -1.], vec![0.5, 0.25]],
        &numbers(&model, "model_parameters"),
    )
    .unwrap();
    let predictions = numbers(&predictions, "predictions");
    near(predictions[0], 23., 1e-9);
    near(predictions[1], 5.5, 1e-9);
}

#[test]
fn response_surface_scales_offset_features_without_losing_predictions() {
    let data = vec![
        vec![1e12 - 2., 9.],
        vec![1e12 - 1., 4.],
        vec![1e12, 1.],
        vec![1e12 + 1., 0.],
        vec![1e12 + 2., 1.],
    ];
    let output = run("response_surface", &data, &[0.]).unwrap();
    near(numbers(&output, "rmse")[0], 0., 1e-9);
    let predictions = run(
        "response_surface_predict",
        &[vec![1e12 + 3.]],
        &numbers(&output, "model_parameters"),
    )
    .unwrap();
    near(numbers(&predictions, "predictions")[0], 4., 1e-8);
}

#[test]
fn ridge_identifies_an_underdetermined_quadratic() {
    let output = run("response_surface", &[vec![0., 2.], vec![1., 4.]], &[1.]).unwrap();
    assert!(numbers(&output, "coefficients")
        .iter()
        .all(|v| v.is_finite()));
}

#[test]
fn invalid_bounds_rank_deficiency_and_malformed_models_are_rejected() {
    assert!(run("factorial_design", &[vec![0., 1.]], &[0.]).is_ok());
    assert!(run("factorial_design", &[vec![1., 1.]], &[0.]).is_err());
    assert!(run("factorial_design", &[vec![0., f64::INFINITY]], &[0.]).is_err());
    assert!(run("factorial_design", &[vec![0., 1.]], &[0.5]).is_err());
    assert!(run("central_composite_design", &[vec![0., 1.]], &[0., 1.]).is_err());
    assert!(run(
        "d_optimal_design",
        &[vec![0., 0.], vec![1., 1.], vec![2., 2.]],
        &[3., 10., 1.]
    )
    .is_err());
    assert!(run("d_optimal_design", &[vec![0.], vec![0.]], &[2., 10., 1.]).is_err());
    assert!(run("response_surface", &[vec![0., 1.], vec![1., 2.]], &[0.]).is_err());
    assert!(run("response_surface", &[vec![0., 1.], vec![1., 2.]], &[-1.]).is_err());
    assert!(run(
        "response_surface_predict",
        &[vec![0.]],
        &[1., 0., 0., 1., 2., 3.]
    )
    .is_err());
    assert!(run("unknown", &[vec![0., 1.]], &[]).is_err());
}

#[test]
fn d_optimal_exchange_is_locally_optimal_on_an_irregular_candidate_set() {
    let candidates = vec![
        vec![-1., -1.],
        vec![1., -0.7],
        vec![0.8, 1.],
        vec![-0.6, 0.7],
        vec![0.1, 0.2],
        vec![-0.2, -0.9],
        vec![0.4, 0.8],
    ];
    let output = run("d_optimal_design", &candidates, &[3., 100., 7.]).unwrap();
    let selected: Vec<usize> = numbers(&output, "indices")
        .into_iter()
        .map(|x| x as usize)
        .collect();
    // For three rows [1,x,y], det(X'X) is the squared oriented triangle area.
    let determinant = |ids: &[usize]| {
        let a = &candidates[ids[0]];
        let b = &candidates[ids[1]];
        let c = &candidates[ids[2]];
        ((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])).powi(2)
    };
    let selected_det = determinant(&selected);
    for incoming in 0..candidates.len() {
        if !selected.contains(&incoming) {
            for slot in 0..3 {
                let mut swapped = selected.clone();
                swapped[slot] = incoming;
                assert!(determinant(&swapped) <= selected_det + 1e-10);
            }
        }
    }
}

#[test]
fn malformed_inputs_and_excessive_work_are_rejected_without_panics() {
    for data in [
        vec![],
        vec![vec![]],
        vec![vec![0., 1.], vec![0.]],
        vec![vec![f64::NAN, 1.]],
    ] {
        assert!(run("response_surface", &data, &[]).is_err());
    }
    assert!(run("factorial_design", &vec![vec![0., 1.]; 17], &[]).is_err());
    assert!(run("factorial_design", &[vec![0., 1.]], &[65536.]).is_err());
    assert!(run(
        "central_composite_design",
        &[vec![0., f64::MAX]],
        &[f64::MAX, 1.]
    )
    .is_err());
    assert!(run("d_optimal_design", &[vec![0.], vec![1.]], &[1., 10., 0.]).is_err());
    assert!(run(
        "response_surface_predict",
        &[vec![0.]],
        &[1., 0., 1., 1., 2., f64::NAN]
    )
    .is_err());
}
