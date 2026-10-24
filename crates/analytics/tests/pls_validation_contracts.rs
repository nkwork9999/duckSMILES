use ducksmiles_analytics::execute;
fn values(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let text = &j[start..];
    let end = if text.starts_with('[') {
        let mut depth = 0;
        text.char_indices()
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
        text.find([',', '}']).unwrap_or(text.len())
    };
    text[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * (1. + b.abs()), "{a} != {b}");
}
fn rows() -> Vec<Vec<f64>> {
    ROWS.into_iter().map(Vec::from).collect()
}
const ROWS: [[f64; 6]; 8] = [
    [0.0, 1.0, 2.0, 1.0, 4.0, 0.0],
    [1.0, 0.0, 2.0, 2.0, 3.0, 1.0],
    [2.0, 1.0, 0.0, 4.0, 1.0, 0.0],
    [3.0, 2.0, 1.0, 5.0, 2.0, 1.0],
    [4.0, 0.0, 3.0, 6.0, 6.0, 0.0],
    [5.0, 3.0, 2.0, 8.0, 3.0, 1.0],
    [6.0, 2.0, 4.0, 7.0, 9.0, 0.0],
    [7.0, 4.0, 1.0, 11.0, 4.0, 1.0],
];
#[test]
fn fold_path_matches_independent_train_only_pls_and_documented_press() {
    let result = execute("pls2_validate", &rows(), &[2., 2.]).unwrap();
    assert_eq!(values(&result, "root_mean_press").len(), 3);
    for (a, b) in values(&result, "root_mean_press").into_iter().zip([
        1.5846313486286048,
        1.5764228031048353,
        1.4531908106739075,
    ]) {
        near(a, b);
    }
    assert_eq!(values(&result, "predictions_by_components").len(), 48);
    for (a, b) in values(&result, "predictions_by_components")
        .into_iter()
        .zip([
            6.5,
            3.0,
            4.5,
            5.0,
            6.5,
            3.0,
            4.5,
            5.0,
            6.5,
            3.0,
            4.5,
            5.0,
            6.5,
            3.0,
            4.5,
            5.0,
            1.9407758569650273,
            2.498435258176537,
            3.011479739087667,
            2.7968809494735036,
            5.451656865528395,
            2.8846707384661223,
            4.167477635131034,
            4.507842535969291,
            2.7089643170202766,
            2.5829444278404114,
            6.118371206023996,
            7.395307963513838,
            5.020986012485123,
            2.8372922134274856,
            6.83626928363528,
            8.457849719010252,
            1.9178032211118818,
            2.4870217834558592,
            3.2257663712876714,
            2.7065863170913134,
            3.1366575276554323,
            1.7345116278570774,
            3.686809002383766,
            4.710383389992585,
            4.6157630263462,
            3.5302975746641314,
            4.907983901354372,
            7.905332670109656,
            8.6071844899087,
            4.6190200885856765,
            6.270648746302904,
            8.696187024337702,
        ])
    {
        near(a, b);
    }
    near(values(&result, "selected_components")[0], 2.);
    let direct = execute(
        "pls2_fit",
        &rows().iter().map(|r| r[..5].to_vec()).collect::<Vec<_>>(),
        &[2., 2.],
    )
    .unwrap();
    assert_eq!(
        values(&result, "predict_parameters"),
        values(&direct, "predict_parameters")
    );
}
#[test]
fn changing_held_out_responses_cannot_change_their_own_fixed_component_predictions() {
    for method in [0., 1.] {
        let data = rows();
        let base = execute("pls2_validate", &data, &[2., 2., method]).unwrap();
        let mut changed = data.clone();
        for (i, r) in changed.iter_mut().enumerate() {
            if i % 2 == 0 {
                r[3] += 1e6;
                r[4] -= 1e6;
            }
        }
        let altered = execute("pls2_validate", &changed, &[2., 2., method]).unwrap();
        let a = values(&base, "predictions_by_components");
        let b = values(&altered, "predictions_by_components");
        for component in 0..3 {
            for i in [0, 2, 4, 6] {
                for response in 0..2 {
                    let index = (component * 8 + i) * 2 + response;
                    near(a[index], b[index]);
                }
            }
        }
    }
}
#[test]
fn holdback_test_rows_never_enter_selection_or_training_transforms() {
    let mut data = rows();
    for (i, r) in data.iter_mut().enumerate() {
        r[5] = if i < 5 { 0. } else { 1. };
    }
    data.push(vec![9., 2., 3., 12., 8., 2.]);
    data.push(vec![10., 3., 4., 15., 9., 2.]);
    for method in [0., 1.] {
        let base = execute("pls2_validate", &data, &[2., 2., method, 1., 1., 1.]).unwrap();
        let mut changed = data.clone();
        for r in &mut changed[8..] {
            r[0] += 1e6;
            r[3] -= 1e6;
            r[4] += 1e6;
        }
        let altered = execute("pls2_validate", &changed, &[2., 2., method, 1., 1., 1.]).unwrap();
        assert_eq!(
            values(&base, "root_mean_press"),
            values(&altered, "root_mean_press")
        );
        assert_eq!(
            values(&base, "predict_parameters"),
            values(&altered, "predict_parameters")
        );
        near(values(&base, "predict_parameters")[5], 5.);
    }
}
#[test]
fn indicator_validation_preserves_training_classes_and_scores() {
    let data = vec![
        vec![-3., 10., 0.],
        vec![-2., 10., 1.],
        vec![-1., 10., 0.],
        vec![1., 30., 1.],
        vec![2., 30., 0.],
        vec![3., 30., 1.],
    ];
    for method in [0., 1.] {
        let result = execute("plsda_validate", &data, &[1., method]).unwrap();
        near(values(&result, "selected_components")[0], 1.);
        assert_eq!(values(&result, "class_labels"), vec![10., 30.]);
        assert_eq!(
            values(&result, "selected_validation_labels"),
            vec![10., 10., 10., 30., 30., 30.]
        );
        let prediction = execute(
            "plsda_predict",
            &[vec![-4.], vec![4.]],
            &values(&result, "predict_parameters"),
        )
        .unwrap();
        assert_eq!(values(&prediction, "predicted_labels"), vec![10., 30.]);
    }
}
#[test]
fn illegal_validation_roles_fold_shapes_and_undefined_variances_fail() {
    let data = rows();
    for parameters in [
        vec![2., 0.],
        vec![2., 2., 0., 1., 1., 2.],
        vec![2., 2., 0., 1., 1., 0., 1.],
    ] {
        assert!(execute("pls2_validate", &data, &parameters).is_err());
    }
    let mut bad = data.clone();
    for r in &mut bad {
        r[5] = 0.;
    }
    assert!(execute("pls2_validate", &bad, &[2., 2.]).is_err());
    bad[0][5] = 1.;
    assert!(execute("pls2_validate", &bad, &[2., 2.]).is_err());
    let constant = vec![
        vec![0., 1., 0.],
        vec![1., 1., 1.],
        vec![2., 1., 0.],
        vec![3., 1., 1.],
    ];
    assert!(execute("pls2_validate", &constant, &[1., 1.]).is_err());
    let missing = vec![
        vec![0., 10., 0.],
        vec![1., 10., 0.],
        vec![2., 30., 1.],
        vec![3., 30., 1.],
    ];
    assert!(execute("plsda_validate", &missing, &[1.]).is_err());
}

#[test]
fn normalized_press_keeps_representable_errors_when_raw_squared_loss_underflows() {
    let noise = [0., 1., -1., 2., -2., 1.];
    let ordinary = (0..6)
        .map(|i| vec![i as f64, 2. * i as f64 + noise[i] * 1e-8, (i % 2) as f64])
        .collect::<Vec<_>>();
    let tiny = ordinary
        .iter()
        .map(|r| vec![r[0], r[1] * 1e-155, r[2]])
        .collect::<Vec<_>>();
    for method in [0., 1.] {
        let a = execute("pls2_validate", &ordinary, &[1., 1., method]).unwrap();
        let b = execute("pls2_validate", &tiny, &[1., 1., method]).unwrap();
        let expected = values(&a, "root_mean_press")[1];
        let actual = values(&b, "root_mean_press")[1];
        assert!(expected > 1e-9);
        assert!(
            (actual - expected).abs() < expected * 1e-6,
            "{actual} != {expected}"
        );
    }
}

#[test]
fn balanced_three_class_intercept_ties_use_smallest_label_independent_of_row_order() {
    let labels = [10., 20., 30.];
    let data = (0..18)
        .map(|i| vec![0., labels[i % 3], (i % 2) as f64])
        .collect::<Vec<_>>();
    let result = execute("plsda_validate", &data, &[1.]).unwrap();
    near(values(&result, "selected_components")[0], 0.);
    assert_eq!(values(&result, "selected_validation_labels"), vec![10.; 18]);
    let training = data.iter().map(|r| r[..2].to_vec()).collect::<Vec<_>>();
    let fit = execute("plsda_fit", &training, &[1.]).unwrap();
    let pred = execute(
        "plsda_predict",
        &[vec![0.]],
        &values(&fit, "predict_parameters"),
    )
    .unwrap();
    assert_eq!(values(&pred, "predicted_labels"), vec![10.]);
}
