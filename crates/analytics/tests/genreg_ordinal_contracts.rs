use ducksmiles_analytics::{ds_analytics_free, ds_analytics_run, execute};

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
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
fn grouped() -> Vec<Vec<f64>> {
    let mut rows = vec![];
    for (x, counts) in [
        (-1., [2., 8., 10.]),
        (0., [5., 10., 5.]),
        (1., [10., 8., 2.]),
    ] {
        for (label, count) in [10., 20., 30.].into_iter().zip(counts) {
            rows.push(vec![x, label, count]);
        }
    }
    rows
}
fn fit(d: &[Vec<f64>]) -> String {
    execute("genreg_ordinal_fit", d, &[]).unwrap_or_else(|e| panic!("ordinal: {e}"))
}

#[test]
fn intercept_only_ordered_mle_reports_cutpoints_jmp_differences_and_information() {
    let j = fit(&[vec![30., 1.], vec![10., 1.], vec![20., 2.]]);
    assert_eq!(values(&j, "ordered_labels"), vec![10., 20., 30.]);
    let a = values(&j, "cutpoints");
    assert_eq!(a.len(), 2);
    near(a[0], -3_f64.ln(), 1e-7);
    near(a[1], 3_f64.ln(), 1e-7);
    let v = values(&j, "genreg_intercepts");
    assert_eq!(v.len(), 2);
    near(v[0], -3_f64.ln(), 1e-7);
    near(v[1], 2. * 3_f64.ln(), 1e-7);
    let cov = values(&j, "covariance");
    assert_eq!(cov.len(), 4);
    for (a, b) in cov.into_iter().zip([4. / 3., 4. / 9., 4. / 9., 4. / 3.]) {
        near(a, b, 1e-6);
    }
    let dcov = values(&j, "genreg_covariance");
    assert_eq!(dcov.len(), 4);
    for (a, b) in dcov
        .into_iter()
        .zip([4. / 3., -8. / 9., -8. / 9., 16. / 9.])
    {
        near(a, b, 1e-6);
    }
    near(values(&j, "log_likelihood")[0], -6. * 2_f64.ln(), 1e-8);
    assert!(j.contains("\"converged\":true"), "{j}");
    assert!(
        j.contains("\"covariance_method\":\"observed_information\""),
        "{j}"
    );
}

#[test]
fn common_slope_has_jmp_positive_cumulative_sign_and_portable_prediction() {
    // Cell frequencies exactly equal the probabilities for cuts +/-ln3,beta=ln3.
    let j = fit(&grouped());
    let b = values(&j, "coefficients");
    assert_eq!(b.len(), 1);
    near(b[0], 3_f64.ln(), 1e-6);
    let a = values(&j, "predict_parameters");
    assert!(a.len() < 16384);
    let pred = execute(
        "genreg_ordinal_predict",
        &[vec![-1., 0.5], vec![0., 0.5], vec![1., 0.5]],
        &a,
    )
    .unwrap();
    let probs = values(&pred, "probabilities");
    assert_eq!(probs.len(), 9);
    for (v, w) in probs
        .into_iter()
        .zip([0.1, 0.4, 0.5, 0.25, 0.5, 0.25, 0.5, 0.4, 0.1])
    {
        near(v, w, 1e-6);
    }
    assert_eq!(values(&pred, "most_likely"), vec![30., 20., 10.]);
    let one = execute("genreg_ordinal_predict", &[vec![0., 0.5]], &a).unwrap();
    let batch = execute(
        "genreg_ordinal_predict",
        &[vec![0., 0.5], vec![1e8, 0.5]],
        &a,
    )
    .unwrap();
    assert_eq!(
        &values(&one, "probabilities"),
        &values(&batch, "probabilities")[..3]
    );
}

#[test]
fn replication_frequencies_and_zero_weight_categories_preserve_likelihood_and_fit() {
    let mut rows = grouped();
    rows.push(vec![80., 999., 0.]);
    let mut expanded = vec![];
    for r in grouped() {
        for _ in 0..r[2] as usize {
            expanded.push(vec![r[0], r[1], 1.]);
        }
    }
    let a = fit(&rows);
    let b = fit(&expanded);
    for key in [
        "cutpoints",
        "coefficients",
        "covariance",
        "log_likelihood",
        "ordered_labels",
    ] {
        let av = values(&a, key);
        let bv = values(&b, key);
        assert_eq!(av.len(), bv.len());
        for (x, y) in av.into_iter().zip(bv) {
            near(x, y, 1e-6);
        }
    }
}

#[test]
fn nonconverged_fit_withholds_uncertainty_and_invalid_design_is_rejected() {
    let j = execute("genreg_ordinal_fit", &grouped(), &[0., 1e-7, 0.95]).unwrap();
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"covariance\":null"), "{j}");
    let pred = execute(
        "genreg_ordinal_predict",
        &[vec![0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert!(
        pred.contains("\"probability_standard_errors\":[[null,null,null]]"),
        "{pred}"
    );
    for data in [
        vec![vec![1., 10., 2.], vec![1., 20., 2.]],
        vec![vec![10., 1.]; 5],
        vec![vec![10., -1.], vec![20., 3.]],
    ] {
        assert!(
            execute("genreg_ordinal_fit", &data, &[]).is_err(),
            "{data:?}"
        );
    }
}

#[test]
fn exported_ordinal_c_abi_owns_finite_json_and_rejects_malformed_carrier() {
    let op = b"genreg_ordinal_fit";
    let d = [10., 1., 20., 2., 30., 1.];
    let mut len = 77;
    let ptr = unsafe {
        ds_analytics_run(
            op.as_ptr(),
            op.len(),
            d.as_ptr(),
            3,
            2,
            std::ptr::null(),
            0,
            &mut len,
        )
    };
    assert!(!ptr.is_null(), "exported ABI returned NULL");
    let j = unsafe { String::from_utf8(std::slice::from_raw_parts(ptr, len).to_vec()).unwrap() };
    unsafe { ds_analytics_free(ptr, len) };
    assert!(j.starts_with("{\"version\":1"), "{j}");
    let op = b"genreg_ordinal_predict";
    let q = [0.5];
    let a = [1., 0., 3.];
    len = 99;
    let ptr = unsafe {
        ds_analytics_run(
            op.as_ptr(),
            op.len(),
            q.as_ptr(),
            1,
            1,
            a.as_ptr(),
            a.len(),
            &mut len,
        )
    };
    assert!(ptr.is_null());
    assert_eq!(len, 0);
}

#[test]
fn certified_complete_and_quasi_separation_are_not_reported_as_finite_mles() {
    for data in [
        vec![
            vec![-2., 10., 1.],
            vec![-1., 10., 1.],
            vec![1., 20., 1.],
            vec![2., 20., 1.],
        ],
        vec![
            vec![-1., 10., 1.],
            vec![0., 10., 1.],
            vec![0., 20., 1.],
            vec![1., 20., 1.],
        ],
    ] {
        let err = execute("genreg_ordinal_fit", &data, &[]).unwrap_err();
        assert!(err.contains("separation"), "wrong reason: {err}");
    }
}

fn supplied_covariance(cov: &[f64]) -> Vec<f64> {
    let mut a = vec![1., 0., 3., 0.95, 1., 10., 20., 30., -3_f64.ln(), 3_f64.ln()];
    a.extend(cov);
    a
}

#[test]
fn predictor_units_shift_cutpoints_without_changing_probabilities_or_raw_linear_predictor() {
    let data = grouped()
        .into_iter()
        .map(|mut r| {
            r[0] = 10. + 3. * r[0];
            r
        })
        .collect::<Vec<_>>();
    let j = fit(&data);
    let beta = 3_f64.ln() / 3.;
    near(values(&j, "coefficients")[0], beta, 1e-6);
    near(values(&j, "cutpoints")[0], -3_f64.ln() - 10. * beta, 1e-6);
    near(values(&j, "genreg_intercepts")[1], 2. * 3_f64.ln(), 1e-6);
    let pred = execute(
        "genreg_ordinal_predict",
        &[vec![10., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    near(values(&pred, "linear_predictor")[0], 10. * beta, 1e-6);
    for (a, b) in values(&pred, "probabilities")
        .into_iter()
        .zip([0.25, 0.5, 0.25])
    {
        near(a, b, 1e-6);
    }
}

#[test]
fn covariance_validation_does_not_invert_and_overflow_a_valid_tiny_covariance() {
    let a = supplied_covariance(&[1e-320, 0., 0., 1e-320]);
    let j = execute("genreg_ordinal_predict", &[vec![0.5]], &a).unwrap();
    let errors = values(&j, "probability_standard_errors");
    assert_eq!(errors.len(), 3);
    assert!(errors.iter().all(|s| *s > 0. && s.is_finite()), "{j}");
}

#[test]
fn tiny_and_huge_materially_asymmetric_covariances_are_rejected() {
    for c in [
        [1e-200, 1e-201, 0., 1e-200],
        [1.7e308, 1.1e308, 1e308, 1.7e308],
    ] {
        let a = supplied_covariance(&c);
        assert!(
            execute("genreg_ordinal_predict", &[vec![0.5]], &a).is_err(),
            "accepted {c:?}"
        );
    }
}

#[test]
fn ordinal_order_ignores_numeric_label_distances_and_multifeature_separation_is_diagnosed() {
    let a = fit(&[vec![10., 1.], vec![20., 2.], vec![30., 1.]]);
    let b = fit(&[vec![-1e100, 1.], vec![0., 2.], vec![1e-100, 1.]]);
    assert_eq!(values(&a, "cutpoints"), values(&b, "cutpoints"));
    let data = [
        vec![-2., 1., 10., 1.],
        vec![1., -2., 10., 1.],
        vec![2., -1., 20., 1.],
        vec![-1., 2., 20., 1.],
    ];
    let error = execute("genreg_ordinal_fit", &data, &[]).unwrap_err();
    assert!(error.contains("separation"), "{error}");
}

#[test]
fn all_preprocessing_and_separation_searches_share_the_cumulative_work_budget() {
    let rows = (0..52000)
        .map(|i| {
            let mut r = (1..=16)
                .map(|j| ((i as f64) * 0.017 * (j as f64)).sin())
                .collect::<Vec<_>>();
            r.push((i % 16) as f64);
            r.push(1.);
            r
        })
        .collect::<Vec<_>>();
    let error = execute("genreg_ordinal_fit", &rows, &[0.]).unwrap_err();
    assert!(error.contains("work"), "wrong failure: {error}");
}

#[test]
fn tiny_predictor_units_keep_normalized_saved_inference_when_raw_covariance_overflows() {
    let data = grouped()
        .into_iter()
        .map(|mut r| {
            r[0] *= 1e-200;
            r
        })
        .collect::<Vec<_>>();
    let j = fit(&data);
    assert!(j.contains("\"information_identified\":true"), "{j}");
    assert!(j.contains("\"covariance\":null"), "{j}");
    let pred = execute(
        "genreg_ordinal_predict",
        &[vec![0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert_eq!(values(&pred, "probability_standard_errors").len(), 3);
    for (a, b) in values(&pred, "probabilities")
        .into_iter()
        .zip([0.25, 0.5, 0.25])
    {
        near(a, b, 1e-6);
    }
}
