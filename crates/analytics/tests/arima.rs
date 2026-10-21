use ducksmiles_analytics::arima::run as execute;

fn number(json: &str, key: &str) -> f64 {
    let tail = &json[json.find(&format!("\"{key}\":")).unwrap() + key.len() + 3..];
    tail.split([',', '}', ']']).next().unwrap().parse().unwrap()
}

fn array(json: &str, key: &str) -> Vec<f64> {
    let tail = &json[json.find(&format!("\"{key}\":[")).unwrap() + key.len() + 4..];
    let body = tail.split(']').next().unwrap();
    if body.is_empty() {
        Vec::new()
    } else {
        body.split(',').map(|v| v.parse().unwrap()).collect()
    }
}

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} != {b}");
}
fn rows(x: &[f64]) -> Vec<Vec<f64>> {
    x.iter().map(|x| vec![*x]).collect()
}
fn fixture() -> Vec<Vec<f64>> {
    rows(&[
        0.4, -0.3, 1.2, 0.8, -0.6, -0.9, 0.2, 1.1, 0.5, -0.8, 0.1, 0.4,
    ])
}

#[test]
fn conditional_ar1_without_mean_recovers_direct_ols_and_forecasts() {
    // Full stationary likelihood or the wrong burn-in would give a different
    // AR1 coefficient from this direct conditional least-squares oracle.
    let out = execute("arima", &fixture(), &[1., 0., 0., 3., 0.]).unwrap();
    near(array(&out, "ar_coefficients")[0], 0.12212389380530976, 2e-6);
    near(number(&out, "sigma_squared"), 0.5059758648431214, 2e-8);
    near(number(&out, "effective_n"), 11., 0.);
    for (a, b) in array(&out, "forecasts").iter().zip([
        0.04884955752212391,
        0.005965698175268231,
        0.0007285542904309876,
    ]) {
        near(*a, b, 2e-6);
    }
    assert!(out.contains("conditional_sum_of_squares"));
    assert!(out.contains("\"converged\":true"));
}

#[test]
fn integrated_random_walk_forecasts_retain_last_level_and_drift_is_explicit() {
    let x = rows(&[1., 3., 2., 5.]);
    let no_drift = execute("arima", &x, &[0., 1., 0., 3., 0.]).unwrap();
    assert_eq!(array(&no_drift, "forecasts"), vec![5., 5., 5.]);
    near(number(&no_drift, "sigma_squared"), 14. / 3., 2e-12);
    let drift = execute("arima", &x, &[0., 1., 0., 3., 1.]).unwrap();
    near(number(&drift, "differenced_mean"), 4. / 3., 2e-12);
    for (a, b) in array(&drift, "forecasts")
        .iter()
        .zip([19. / 3., 23. / 3., 9.])
    {
        near(*a, b, 2e-12);
    }
    near(number(&drift, "sigma_squared"), 26. / 9., 2e-12);
}

#[test]
fn seasonal_and_ordinary_differencing_are_inverted_through_each_history_layer() {
    let x = rows(&[1., 4., 2., 5., 2., 6., 3., 8., 4., 9., 5., 11.]);
    let out = execute(
        "arima",
        &x,
        &[0., 1., 0., 5., 0., 1000., 1e-8, 0., 1., 0., 3.],
    )
    .unwrap();
    assert_eq!(array(&out, "forecasts"), vec![16., 12., 18., 23., 19.]);
    assert_eq!(array(&out, "differencing_lags"), vec![1., 3.]);
}

#[test]
fn ma1_css_matches_scipy_filter_and_independent_bounded_optimizer() {
    let out = execute("arima", &fixture(), &[0., 0., 1., 3., 0.]).unwrap();
    near(array(&out, "ma_coefficients")[0], 0.48171781511050243, 3e-6);
    near(number(&out, "sigma_squared"), 0.46041757364146496, 2e-8);
    near(array(&out, "forecasts")[0], 0.07563809521332267, 3e-6);
    near(
        number(&out, "conditional_log_likelihood"),
        -11.342405986577798,
        2e-8,
    );
    near(array(&out, "forecasts")[1], 0., 0.);
}

#[test]
fn seasonal_ar_and_ma_polynomials_compose_and_remain_stationary_invertible() {
    let x: Vec<_> = (0..90)
        .map(|i| vec![(i as f64 * 0.71).sin() + 0.25 * (i as f64 * 0.13).cos()])
        .collect();
    let out = execute("arima", &x, &[1., 0., 1., 3., 0., 2., 1e-8, 1., 0., 1., 4.]).unwrap();
    let ar = array(&out, "ar_coefficients")[0];
    let sar = array(&out, "seasonal_ar_coefficients")[0];
    let ma = array(&out, "ma_coefficients")[0];
    let sma = array(&out, "seasonal_ma_coefficients")[0];
    assert!(ar.abs() < 1. && sar.abs() < 1. && ma.abs() < 1. && sma.abs() < 1.);
    let ap = array(&out, "ar_polynomial");
    let mp = array(&out, "ma_polynomial");
    near(ap[1], -ar, 2e-12);
    near(ap[4], -sar, 2e-12);
    near(ap[5], ar * sar, 2e-12);
    near(mp[1], ma, 2e-12);
    near(mp[4], sma, 2e-12);
    near(mp[5], ma * sma, 2e-12);
    assert!(out.contains("\"converged\":false"));
}

#[test]
fn arima_rejects_invalid_orders_shapes_and_unsupported_work() {
    for args in [
        vec![],
        vec![4., 0., 0., 1.],
        vec![0., 3., 0., 1.],
        vec![0., 0., 0., -1.],
        vec![0., 0., 0., 1., 0.5],
        vec![0., 0., 0., 1., 0., 1000., 0.],
        vec![0., 0., 0., 1., 0., 1000., 1e-8, 1., 0., 0., 0.],
    ] {
        assert!(execute("arima", &fixture(), &args).is_err(), "{args:?}");
    }
    assert!(execute("arima", &[vec![f64::NAN]], &[0., 0., 0., 1.]).is_err());
    assert!(execute("arima", &[vec![1., 2.]], &[0., 0., 0., 1.]).is_err());
    assert!(execute("arima", &rows(&[2., 2., 2., 2.]), &[0., 0., 0., 1., 1.]).is_err());
    assert!(execute(
        "arima",
        &vec![vec![1.]; 20000],
        &[3., 0., 3., 1., 0., 10000., 1e-8, 1., 0., 1., 52.]
    )
    .is_err());
}

#[test]
fn conditional_mean_model_preserves_representable_variation_at_large_origin() {
    let x = rows(&[
        1., -0.25, 0.75, 0.5, -1., 0.25, 0.125, 0.875, -0.5, 1.25, -0.125, 0.625,
    ]);
    let origin = 281474976710656.0;
    let translated: Vec<_> = x.iter().map(|r| vec![r[0] + origin]).collect();
    let base = execute("arima", &x, &[1., 0., 0., 2., 1.]).unwrap();
    let shifted = execute("arima", &translated, &[1., 0., 0., 2., 1.]).unwrap();
    near(
        array(&shifted, "ar_coefficients")[0],
        array(&base, "ar_coefficients")[0],
        2e-6,
    );
    near(
        number(&shifted, "sigma_squared"),
        number(&base, "sigma_squared"),
        2e-7,
    );
    for (a, b) in array(&shifted, "forecasts")
        .iter()
        .zip(array(&base, "forecasts"))
    {
        near(*a - origin, b, 0.07);
    }
}

#[test]
fn exact_ar_decay_does_not_report_a_finite_gaussian_variance_fit() {
    let x = rows(&[1., 0.5, 0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125]);
    assert!(execute("arima", &x, &[1., 0., 0., 2., 0.]).is_err());
}
