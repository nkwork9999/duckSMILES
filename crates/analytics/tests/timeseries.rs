use ducksmiles_analytics::timeseries::run;

fn number(j: &str, key: &str) -> f64 {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = &j[start..];
    tail[..tail.find([',', '}', ']']).unwrap_or(tail.len())]
        .parse()
        .unwrap()
}
fn array(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &j[start..];
    let s = &tail[..tail.find(']').unwrap()];
    if s.is_empty() {
        vec![]
    } else {
        s.split(',').map(|s| s.parse().unwrap()).collect()
    }
}
fn nullable_array(j: &str, key: &str) -> Vec<Option<f64>> {
    let start = j.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let tail = &j[start..];
    tail[..tail.find(']').unwrap()]
        .split(',')
        .map(|s| {
            if s == "null" {
                None
            } else {
                Some(s.parse().unwrap())
            }
        })
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-9 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}
fn close_array(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&x, &y) in actual.iter().zip(expected) {
        close(x, y);
    }
}
fn x() -> Vec<Vec<f64>> {
    [2., 4., 3., 5., 7., 6., 8., 9.].map(|v| vec![v]).to_vec()
}

#[test]
fn acf_uses_biased_covariance_with_a_common_denominator() {
    // NumPy correlate(x-mean,x-mean) / dot(x-mean,x-mean).
    let j = run("acf", &x(), &[4.]).unwrap();
    close_array(
        &array(&j, "autocorrelations"),
        &[
            1.,
            0.48214285714285715,
            0.2619047619047619,
            0.05357142857142857,
            -1. / 3.,
        ],
    );
    assert_eq!(array(&j, "lags"), vec![0., 1., 2., 3., 4.]);
}
#[test]
fn pacf_matches_independent_toeplitz_solutions_for_each_lag() {
    let j = run("pacf", &x(), &[4.]).unwrap();
    close_array(
        &array(&j, "partial_autocorrelations"),
        &[
            1.,
            0.48214285714285715,
            0.038360337903337476,
            -0.1126752080124592,
            -0.43321779097253404,
        ],
    );
}
#[test]
fn difference_repeats_the_requested_lag_and_preserves_original_indices() {
    let j = run("difference", &x(), &[2., 2.]).unwrap();
    assert_eq!(array(&j, "values"), vec![3., 0., -3., 2.]);
    assert_eq!(array(&j, "original_indices"), vec![4., 5., 6., 7.]);
    assert_eq!(
        array(&run("difference", &x(), &[]).unwrap(), "values"),
        vec![2., -1., 2., 2., -1., 2., 1.]
    );
}
#[test]
fn moving_average_is_trailing_and_never_uses_future_values() {
    let j = run("moving_average", &x(), &[3.]).unwrap();
    close_array(&array(&j, "values"), &[3., 4., 5., 6., 7., 23. / 3.]);
    assert_eq!(array(&j, "original_indices"), vec![2., 3., 4., 5., 6., 7.]);
    let mut extended = x();
    extended.push(vec![1000.]);
    let more = array(&run("moving_average", &extended, &[3.]).unwrap(), "values");
    close_array(&more[..6], &array(&j, "values"));
}
#[test]
fn exponential_smoothing_initializes_at_first_value_and_forecasts_final_level() {
    let j = run("exponential_smoothing", &x(), &[0.4, 3.]).unwrap();
    close_array(
        &array(&j, "levels"),
        &[2., 2.8, 2.88, 3.728, 5.0368, 5.42208, 6.453248, 7.4719488],
    );
    close_array(&array(&j, "forecasts"), &[7.4719488; 3]);
    close_array(
        &array(&j, "fitted_values"),
        &[2., 2.8, 2.88, 3.728, 5.0368, 5.42208, 6.453248],
    );
    assert_eq!(
        array(&j, "fitted_indices"),
        vec![1., 2., 3., 4., 5., 6., 7.]
    );
    let mut extended = x();
    extended.push(vec![1000.]);
    close_array(
        &array(
            &run("exponential_smoothing", &extended, &[0.4, 1.]).unwrap(),
            "levels",
        )[..8],
        &array(&j, "levels"),
    );
}
#[test]
fn holt_damps_both_updates_and_each_future_trend_increment() {
    let j = run("holt", &x(), &[0.4, 0.3, 3., 0.9]).unwrap();
    close_array(
        &array(&j, "forecasts"),
        &[9.330366976752705, 9.917344847849883, 10.445624931837344],
    );
    close(number(&j, "final_level"), 8.678169342200285);
    close(number(&j, "final_trend"), 0.7246640383915783);
    close_array(
        &array(&j, "fitted_values"),
        &[
            3.8,
            5.5216,
            5.7180672,
            6.4378855424,
            7.6297803870208005,
            7.672196105836954,
            8.463615570333808,
        ],
    );
    let j = run(
        "holt",
        &[vec![1.], vec![2.], vec![3.], vec![4.]],
        &[0.5, 0.5, 2.],
    )
    .unwrap();
    close_array(&array(&j, "forecasts"), &[5., 6.]);
}
#[test]
fn ar_conditional_qr_fit_has_correct_lag_order_forecasts_and_gaussian_aic() {
    // NumPy lstsq([1,x[t-1],x[t-2]],x[t]); likelihood uses n-p observations,
    // variance SSE/(n-p), and counts intercept + slopes + variance in AIC/BIC.
    let j = run("ar_model", &x(), &[2., 3.]).unwrap();
    close_array(
        &array(&j, "coefficients"),
        &[0.4511494252873563, 0.617816091954023],
    );
    close(number(&j, "intercept"), 1.07183908045977);
    close_array(
        &array(&j, "residuals"),
        &[
            -1.11206896551724,
            0.103448275862069,
            1.81896551724138,
            -1.31896551724138,
            -0.103448275862069,
            0.612068965517241,
        ],
    );
    close(number(&j, "innovation_variance"), 1.1135057471264371);
    close(number(&j, "residual_variance_unbiased"), 2.2270114942528743);
    close(number(&j, "log_likelihood"), -8.83617130628482);
    close(number(&j, "aic"), 25.67234261256964);
    close(number(&j, "bic"), 24.83938048948186);
    close_array(
        &array(&j, "forecasts"),
        &[10.07471264367816, 11.177384727176637, 12.338829369418205],
    );
    assert_eq!(array(&j, "fitted_indices"), vec![2., 3., 4., 5., 6., 7.]);
}
#[test]
fn periodogram_matches_scipy_one_sided_density_and_parseval() {
    // scipy.signal.periodogram(fs=0.5,detrend='constant',scaling='density').
    let j = run("periodogram", &x(), &[2.]).unwrap();
    close_array(
        &array(&j, "frequencies"),
        &[0., 0.0625, 0.125, 0.1875, 0.25],
    );
    close_array(
        &array(&j, "power_spectral_density"),
        &[0., 49.14213562373095, 10., 20.85786437626905, 4.],
    );
    close(
        array(&j, "power_spectral_density").iter().sum::<f64>() / 16.,
        5.25,
    );
    let odd = run("periodogram", &[vec![1.], vec![0.], vec![-1.]], &[]).unwrap();
    close(
        array(&odd, "power_spectral_density").iter().sum::<f64>() / 3.,
        2. / 3.,
    );
}
#[test]
fn seasonal_decomposition_centers_even_period_filters_and_marks_unidentified_edges() {
    let season = [0., 1., -2., 1.];
    let rows = (0..12)
        .map(|i| vec![10. + 0.25 * i as f64 + season[i % 4]])
        .collect::<Vec<_>>();
    let j = run("seasonal_decompose", &rows, &[4.]).unwrap();
    let trend = nullable_array(&j, "trend");
    let residuals = nullable_array(&j, "residuals");
    assert_eq!(trend.len(), 12);
    assert!(trend[..2].iter().all(Option::is_none));
    assert!(trend[10..].iter().all(Option::is_none));
    close_array(&array(&j, "seasonal_pattern"), &season);
    for i in 2..10 {
        close(trend[i].unwrap(), 10. + 0.25 * i as f64);
        close(residuals[i].unwrap(), 0.);
    }
    let odd = (0..12)
        .map(|i| vec![10. + 0.5 * i as f64 + [-1., 0., 1.][i % 3]])
        .collect::<Vec<_>>();
    let j = run("seasonal_decompose", &odd, &[3.]).unwrap();
    close_array(&array(&j, "seasonal_pattern"), &[-1., 0., 1.]);
    assert!(nullable_array(&j, "trend")[0].is_none());
}
#[test]
fn ljung_box_uses_sample_size_correction_and_subtracts_model_df() {
    let j = run("ljung_box", &x(), &[4.]).unwrap();
    close(number(&j, "q_statistic"), 5.839434186372962);
    close(number(&j, "p_value"), 0.21146461232780495);
    close(number(&j, "degrees_of_freedom"), 4.);
    let j = run("ljung_box", &x(), &[4., 2.]).unwrap();
    close(number(&j, "degrees_of_freedom"), 2.);
    close(number(&j, "p_value"), 0.0539489476665005);
}
#[test]
fn constant_smoothing_and_transforms_work_while_unidentified_correlations_reject() {
    let rows = vec![vec![3.]; 12];
    for op in ["acf", "pacf", "ljung_box"] {
        assert!(run(op, &rows, &[3.]).is_err(), "{op}");
    }
    assert!(run("ar_model", &rows, &[1., 2.]).is_err());
    assert_eq!(
        array(&run("difference", &rows, &[]).unwrap(), "values"),
        vec![0.; 11]
    );
    close_array(
        &array(&run("moving_average", &rows, &[3.]).unwrap(), "values"),
        &[3.; 10],
    );
    close_array(
        &array(
            &run("exponential_smoothing", &rows, &[0.3, 2.]).unwrap(),
            "forecasts",
        ),
        &[3.; 2],
    );
    close_array(
        &array(&run("holt", &rows, &[0.3, 0.5, 2.]).unwrap(), "forecasts"),
        &[3.; 2],
    );
    close_array(
        &array(
            &run("seasonal_decompose", &rows, &[4.]).unwrap(),
            "seasonal_pattern",
        ),
        &[0.; 4],
    );
    close_array(
        &array(
            &run("periodogram", &rows, &[]).unwrap(),
            "power_spectral_density",
        ),
        &[0.; 7],
    );
    close_array(
        &array(
            &run("exponential_smoothing", &[vec![7.]], &[0., 2.]).unwrap(),
            "forecasts",
        ),
        &[7.; 2],
    );
}
#[test]
fn malformed_requests_and_work_explosions_reject_before_computation() {
    for op in [
        "acf",
        "pacf",
        "difference",
        "moving_average",
        "exponential_smoothing",
        "holt",
        "ar_model",
        "periodogram",
        "seasonal_decompose",
        "ljung_box",
    ] {
        for rows in [
            vec![],
            vec![vec![]],
            vec![vec![1., 2.]],
            vec![vec![f64::NAN]],
            vec![vec![1.]; 100_001],
        ] {
            assert!(run(op, &rows, &[]).is_err(), "{op}");
        }
        assert!(run(op, &x(), &[f64::NAN]).is_err(), "{op}");
    }
    for (op, args) in [
        ("acf", vec![8.]),
        ("acf", vec![1.5]),
        ("acf", vec![-1.]),
        ("pacf", vec![8.]),
        ("difference", vec![0.]),
        ("difference", vec![2., 4.]),
        ("difference", vec![1., 0.]),
        ("moving_average", vec![0.]),
        ("moving_average", vec![9.]),
        ("exponential_smoothing", vec![-0.1, 1.]),
        ("exponential_smoothing", vec![0.5, 10001.]),
        ("holt", vec![0.5, 1.1, 2.]),
        ("holt", vec![0.5, 0.3, 2., 0.]),
        ("ar_model", vec![4., 2.]),
        ("ar_model", vec![33., 2.]),
        ("periodogram", vec![0.]),
        ("seasonal_decompose", vec![1.]),
        ("seasonal_decompose", vec![5.]),
        ("ljung_box", vec![4., 4.]),
        ("ljung_box", vec![4., -1.]),
    ] {
        assert!(run(op, &x(), &args).is_err(), "{op}: {args:?}");
    }
    assert!(run("periodogram", &vec![vec![0.]; 4097], &[]).is_err());
    assert!(run(
        "acf",
        &(0..100_000).map(|i| vec![i as f64]).collect::<Vec<_>>(),
        &[1000.]
    )
    .is_err());
    assert!(run("no_such_operation", &x(), &[]).is_err());
}
#[test]
fn translation_scaling_and_time_order_properties_hold_for_nonconstant_series() {
    let mut seed = 0x51de_u64;
    for _ in 0..24 {
        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
        let shift = ((seed >> 32) % 1000) as f64 - 500.;
        let scale = ((seed >> 16) % 10 + 1) as f64;
        let transformed = x()
            .iter()
            .map(|r| vec![shift + scale * r[0]])
            .collect::<Vec<_>>();
        for (op, key, args) in [
            ("acf", "autocorrelations", vec![4.]),
            ("pacf", "partial_autocorrelations", vec![4.]),
        ] {
            close_array(
                &array(&run(op, &transformed, &args).unwrap(), key),
                &array(&run(op, &x(), &args).unwrap(), key),
            );
        }
        for (op, args) in [
            ("ar_model", vec![2., 3.]),
            ("exponential_smoothing", vec![0.4, 3.]),
            ("holt", vec![0.4, 0.3, 3., 0.9]),
        ] {
            let actual = array(&run(op, &transformed, &args).unwrap(), "forecasts");
            let expected = array(&run(op, &x(), &args).unwrap(), "forecasts")
                .iter()
                .map(|f| shift + scale * f)
                .collect::<Vec<_>>();
            close_array(&actual, &expected);
        }
    }
    let j = run(
        "exponential_smoothing",
        &x().into_iter().rev().collect::<Vec<_>>(),
        &[0.4, 1.],
    )
    .unwrap();
    assert!((array(&j, "forecasts")[0] - 7.4719488).abs() > 1.);
}

#[test]
fn ar_forecasts_keep_sub_ulp_state_at_large_observation_origins() {
    // Recursive forecasts must retain centered fractions internally even when
    // the displayed original-unit forecast rounds to a multiple of two.
    let rows = [4., 8., 6., 10., 14., 12., 16., 18.].map(|v| vec![1e16 + v]);
    let j = run("ar_model", &rows, &[2., 6.]).unwrap();
    let forecasts = array(&j, "forecasts");
    let offsets = forecasts.iter().map(|v| v - 1e16).collect::<Vec<_>>();
    assert_eq!(offsets, vec![20., 22., 24., 28., 30., 32.]);
}
