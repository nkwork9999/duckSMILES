use ducksmiles_analytics::execute;

fn number(j: &str, key: &str) -> f64 {
    let p = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    j[p..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
fn values(j: &str, key: &str) -> Vec<f64> {
    let p = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let s = &j[p..p + j[p..].find(']').unwrap()];
    if s.is_empty() {
        vec![]
    } else {
        s.split(',').map(|x| x.parse().unwrap()).collect()
    }
}
fn matrix(j: &str, key: &str) -> Vec<Vec<f64>> {
    let p = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let s = &j[p..p + j[p..].find("]]").unwrap() + 1];
    s.split("],")
        .map(|r| {
            r.trim_matches(['[', ']'])
                .split(',')
                .map(|x| x.parse().unwrap())
                .collect()
        })
        .collect()
}
fn close(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}

#[test]
fn exponential_fit_recovers_scale_rate_and_portable_predictions() {
    // Mutation caught: optimizing logarithmic response instead of ordinary residual squares.
    let d = [
        vec![0., 2.03],
        vec![0.5, 2.64],
        vec![1., 3.7],
        vec![1.5, 4.88],
        vec![2., 6.68],
        vec![2.5, 8.94],
    ];
    let j = execute("curve_fit", &d, &[0., 500., 1e-11, 1.8, 0.4]).unwrap();
    close(values(&j, "coefficients")[0], 2.0045698299903227, 1e-7);
    close(values(&j, "coefficients")[1], 0.5987386163457715, 1e-7);
    close(number(&j, "sse"), 0.011118711912985363, 1e-9);
    assert!(j.contains("\"converged\":true"));
    let q = execute(
        "curve_predict",
        &[vec![0.], vec![1.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    close(values(&q, "predictions")[0], 2.0045698299903227, 1e-7);
    close(values(&q, "predictions")[1], 3.647959993479304, 1e-7);
    assert!(number(&j, "residual_df") == 4.);
    assert!(matrix(&j, "covariance")[0][0] > 0.);
}

#[test]
fn logistic_gompertz_power_and_michaelis_families_recover_exact_curves() {
    // Fixtures from the named formulas, independent of fitting code.
    let logistic = [
        vec![-2., 1.1897034927102672],
        vec![-1., 1.5674042596019513],
        vec![0., 2.4173747750968184],
        vec![1., 3.5826252249031816],
        vec![2., 4.43259574039805],
        vec![3., 4.810296507289733],
    ];
    let j = execute("curve_fit", &logistic, &[2., 500., 1e-11, 1., 0., 0.8, 5.2]).unwrap();
    for (a, b) in values(&j, "coefficients").iter().zip([1.2, 0.5, 1., 5.]) {
        close(*a, b, 1e-6);
    }
    let gompertz = [
        vec![-1., 1.0812764525692766],
        vec![0., 1.6946284788336938],
        vec![1., 2.821511225186546],
        vec![2., 3.80903126726258],
        vec![3., 4.412608432210386],
        vec![4., 4.724501776323641],
    ];
    let j = execute(
        "curve_fit",
        &gompertz,
        &[3., 500., 1e-11, 0.8, 5.2, 0.7, 0.5],
    )
    .unwrap();
    for (a, b) in values(&j, "coefficients").iter().zip([1., 5., 0.8, 0.7]) {
        close(*a, b, 1e-5);
    }
    let power = [
        vec![1., 5.],
        vec![2., 10.485281374238571],
        vec![3., 17.588457268119896],
        vec![4., 26.],
        vec![5., 35.54101966249685],
    ];
    let j = execute("curve_fit", &power, &[4., 500., 1e-11, 1., 2., 1.]).unwrap();
    for (a, b) in values(&j, "coefficients").iter().zip([2., 3., 1.5]) {
        close(*a, b, 1e-6);
    }
    let mm = [
        vec![1., 2.],
        vec![2., 3.2],
        vec![3., 4.],
        vec![4., 4.571428571428571],
        vec![6., 5.333333333333333],
    ];
    let j = execute("curve_fit", &mm, &[5., 500., 1e-11, 6., 2.]).unwrap();
    close(values(&j, "coefficients")[0], 8., 1e-6);
    close(values(&j, "coefficients")[1], 3., 1e-6);
}

#[test]
fn cubic_bspline_partition_endpoint_and_linear_reproduction() {
    // Mutation caught: excluding the right endpoint or using wrong knot multiplicities.
    let a = [3., 0., 0., 0., 0., 1., 2., 2., 2., 2.];
    let j = execute(
        "bspline_basis",
        &[vec![0.], vec![0.5], vec![1.], vec![2.]],
        &a,
    )
    .unwrap();
    let b = matrix(&j, "basis");
    assert_eq!(b[0], vec![1., 0., 0., 0., 0.]);
    assert_eq!(b[3], vec![0., 0., 0., 0., 1.]);
    for row in &b {
        close(row.iter().sum(), 1., 1e-14);
    }
    for (a, b) in b[1].iter().zip([0.125, 0.59375, 0.25, 0.03125, 0.]) {
        close(*a, b, 1e-14);
    }
    assert!(execute("bspline_basis", &[vec![2.1]], &a).is_err());
}

#[test]
fn cubic_regression_spline_recovers_polynomial_and_predicts_from_saved_coefficients() {
    // Mutation caught: applying penalty to an unpenalized fit or forgetting saved response offset.
    let d = (0..9)
        .map(|i| {
            let x = i as f64 / 4.;
            vec![x, 3. + 2. * x - x * x + 0.5 * x * x * x]
        })
        .collect::<Vec<_>>();
    let j = execute(
        "spline_fit",
        &d,
        &[3., 0., 0., 0., 0., 0., 1., 2., 2., 2., 2.],
    )
    .unwrap();
    close(number(&j, "sse"), 0., 1e-20);
    close(number(&j, "effective_df"), 5., 1e-11);
    let q = execute(
        "spline_predict",
        &[vec![0.5], vec![1.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    close(values(&q, "predictions")[0], 3.8125, 1e-11);
    close(values(&q, "predictions")[1], 5.4375, 1e-11);
    let penalized = execute(
        "spline_fit",
        &d,
        &[3., 100., 0., 0., 0., 0., 1., 2., 2., 2., 2.],
    )
    .unwrap();
    assert!(number(&penalized, "effective_df") < 5.);
    assert!(number(&penalized, "effective_df") > 2.);
}

#[test]
fn loess_local_quadratic_is_exact_and_reuses_saved_training_state() {
    // Mutation caught: degree-one local regression or response leakage from query rows.
    let d = (0..11)
        .map(|i| {
            let x = i as f64;
            vec![x, 2. + 3. * x + 0.5 * x * x]
        })
        .collect::<Vec<_>>();
    let j = execute("loess", &d, &[0.7, 2., 0.]).unwrap();
    close(number(&j, "sse"), 0., 1e-20);
    let q = execute(
        "loess_predict",
        &[vec![0.5], vec![5.5], vec![9.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    for (a, b) in values(&q, "predictions")
        .iter()
        .zip([3.625, 33.625, 75.625])
    {
        close(*a, b, 1e-10);
    }
}

#[test]
fn nonlinear_limits_and_rank_failures_are_explicit() {
    assert!(execute("curve_fit", &[vec![0., 1.], vec![1., 2.]], &[2.]).is_err());
    assert!(execute("curve_predict", &[vec![-1.]], &[4., 0., 1., 0.5]).is_err());
    assert!(execute("bspline_basis", &[vec![0.5]], &[3., 0., 0., 0., 1., 1., 1.]).is_err());
    assert!(execute(
        "spline_fit",
        &[vec![0., 1.], vec![0., 2.], vec![0., 3.]],
        &[1., 0., 0., 0., 1., 1.]
    )
    .is_err());
    assert!(execute(
        "loess",
        &[vec![0., 1.], vec![0., 2.], vec![0., 3.], vec![0., 4.]],
        &[1., 1., 0.]
    )
    .is_err());
    assert!(execute(
        "curve_fit",
        &[vec![0., 1.], vec![1., 2.], vec![2., 3.]],
        &[0., 1e9]
    )
    .is_err());
}

#[test]
fn robust_loess_downweights_an_isolated_response_outlier() {
    // Mutation caught: ignoring robust reweighting or saving weights from the wrong iteration.
    let mut d = (0..21)
        .map(|i| vec![i as f64, 1. + 2. * i as f64])
        .collect::<Vec<_>>();
    d[10][1] += 100.;
    let plain = execute("loess", &d, &[0.6, 1., 0.]).unwrap();
    let robust = execute("loess", &d, &[0.6, 1., 2.]).unwrap();
    assert!((values(&plain, "fitted_values")[10] - 21.).abs() > 5.);
    close(values(&robust, "fitted_values")[10], 21., 1e-8);
    close(values(&robust, "robust_weights")[10], 0., 0.);
    let saved = execute(
        "loess_predict",
        &[vec![10.]],
        &values(&robust, "predict_parameters"),
    )
    .unwrap();
    close(
        values(&saved, "predictions")[0],
        values(&robust, "fitted_values")[10],
        1e-12,
    );
}

#[test]
fn nonlinear_coefficient_uncertainty_survives_very_small_response_units() {
    // Mutation caught: forming an unscaled Jacobian Gram matrix or squaring a tiny SE first.
    let d = [
        vec![0., 2.03],
        vec![0.5, 2.64],
        vec![1., 3.7],
        vec![1.5, 4.88],
        vec![2., 6.68],
        vec![2.5, 8.94],
    ];
    let normal = execute("curve_fit", &d, &[0., 500., 1e-11, 1.8, 0.4]).unwrap();
    let tiny = d
        .iter()
        .map(|r| vec![r[0], r[1] * 1e-200])
        .collect::<Vec<_>>();
    let j = execute("curve_fit", &tiny, &[0., 500., 1e-11, 1.8e-200, 0.4]).unwrap();
    let se = values(&j, "standard_errors");
    close(se[0] / 1e-200, values(&normal, "standard_errors")[0], 1e-7);
    close(se[1], values(&normal, "standard_errors")[1], 1e-7);
}

#[test]
fn spline_rejects_a_knot_domain_whose_width_overflows() {
    // Mutation caught: accepting all-zero bases after infinite knot denominators.
    assert!(execute(
        "bspline_basis",
        &[vec![0.]],
        &[1., -1e308, -1e308, 1e308, 1e308]
    )
    .is_err());
}
