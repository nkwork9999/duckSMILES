use ducksmiles_analytics::{ds_analytics_free, ds_analytics_run, execute};

fn values(j: &str, key: &str) -> Vec<f64> {
    let s = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let t = &j[s..];
    let e = if t.starts_with('[') {
        let mut depth = 0;
        t.char_indices()
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
        t.find([',', '}']).unwrap_or(t.len())
    };
    t[..e]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64, t: f64) {
    assert!((a - b).abs() <= t * (1. + b.abs()), "{a} != {b}");
}
fn exact(y: &[f64]) -> Vec<Vec<f64>> {
    y.iter().map(|&v| vec![v, v, 0., 1.]).collect()
}
fn fit(op: &str, d: &[Vec<f64>]) -> String {
    execute(op, d, &[]).unwrap_or_else(|e| panic!("{op}: {e}"))
}

#[test]
fn censored_normal_uncensored_intercept_has_mle_scale_and_information() {
    // This catches using RMSE instead of censored-MLE scale, or missing scale information.
    let j = fit("censored_normal", &exact(&[-2., -1., 0., 1., 2.]));
    near(values(&j, "coefficients")[0], 0., 1e-6);
    near(values(&j, "sigma")[0], 2_f64.sqrt(), 1e-6);
    near(values(&j, "standard_errors")[0], (2_f64 / 5.).sqrt(), 1e-5);
    near(values(&j, "standard_errors")[1], (1_f64 / 10.).sqrt(), 1e-5);
    assert!(j.contains("\"converged\":true"), "{j}");
}
#[test]
fn exponential_right_censoring_uses_total_weighted_exposure() {
    // Failures3/exposure11 => mean11/3, sigma absent, information3.
    let j = fit(
        "censored_exponential",
        &[
            vec![1., 1., 0., 2.],
            vec![3., 3., 0., 1.],
            vec![6., 0., 2., 1.],
        ],
    );
    near(values(&j, "coefficients")[0], (11_f64 / 3.).ln(), 1e-6);
    near(values(&j, "standard_errors")[0], 1_f64 / 3_f64.sqrt(), 1e-5);
    near(values(&j, "effective_sample_size")[0], 4., 1e-12);
}
#[test]
fn all_six_families_have_exact_left_right_interval_likelihoods() {
    // All interval types must contribute, not disappear behind optional event flags.
    for (op, positive) in [
        ("censored_normal", false),
        ("censored_exponential", true),
        ("censored_gamma", true),
        ("censored_weibull", true),
        ("censored_lognormal", true),
        ("censored_beta", true),
    ] {
        let mut d = if op == "censored_beta" {
            exact(&[0.15, 0.25, 0.4, 0.55, 0.75])
        } else if positive {
            exact(&[0.4, 0.8, 1.2, 2.0, 3.0])
        } else {
            exact(&[-1., -0.2, 0.4, 1.2, 2.])
        };
        if op == "censored_beta" {
            d.extend([
                vec![0., 0.2, 1., 1.],
                vec![0.7, 0., 2., 1.],
                vec![0.3, 0.6, 3., 2.],
            ]);
        } else {
            d.extend([
                vec![0., 0.5, 1., 1.],
                vec![2.5, 0., 2., 1.],
                vec![0.6, 1.5, 3., 2.],
            ]);
        }
        let j = fit(op, &d);
        assert!(values(&j, "log_likelihood")[0].is_finite(), "{op}: {j}");
        assert!(j.contains("\"censor_counts\":[5,1,1,2]"), "{op}: {j}");
        let a = values(&j, "predict_parameters");
        assert!(a.len() < 16384);
        let v = if op == "censored_beta" { 0.45 } else { 1. };
        let p = execute("censored_predict", &[vec![v, 0.5]], &a).unwrap();
        near(values(&p, "cdf")[0] + values(&p, "survival")[0], 1., 1e-10);
        assert!(values(&p, "quantiles")[0].is_finite(), "{p}");
    }
}
#[test]
fn saved_normal_prediction_has_mean_probability_quantile_and_uncertainty() {
    let j = fit("censored_normal", &exact(&[-2., -1., 0., 1., 2.]));
    let p = execute(
        "censored_predict",
        &[vec![0., 0.5], vec![2_f64.sqrt(), 0.975]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    near(values(&p, "means")[0], 0., 1e-6);
    near(values(&p, "cdf")[0], 0.5, 1e-6);
    near(values(&p, "cdf")[1], 0.8413447460685429, 1e-6);
    near(values(&p, "quantiles")[1], 2.771807648699356, 1e-6);
    near(
        values(&p, "mean_standard_errors")[0],
        (2_f64 / 5.).sqrt(),
        1e-5,
    );
}
#[test]
fn frequency_replication_and_zero_weight_rows_preserve_estimates() {
    let a = fit(
        "censored_normal",
        &[
            vec![-2., -2., 0., 2.],
            vec![1., 1., 0., 3.],
            vec![8., 8., 0., 0.],
        ],
    );
    let b = fit("censored_normal", &exact(&[-2., -2., 1., 1., 1.]));
    for key in ["coefficients", "sigma", "standard_errors", "log_likelihood"] {
        let av = values(&a, key);
        let bv = values(&b, key);
        assert_eq!(av.len(), bv.len());
        for (x, y) in av.iter().zip(bv) {
            near(*x, y, 2e-5)
        }
    }
}
#[test]
fn malformed_censor_support_rank_and_saved_models_are_rejected() {
    for (op, d, args) in [
        ("censored_beta", exact(&[0., 0.2, 0.4, 0.6]), vec![]),
        ("censored_gamma", exact(&[-1., 1., 2., 3.]), vec![]),
        ("censored_normal", vec![vec![1., 0., 3., 1.]; 5], vec![]),
        ("censored_normal", vec![vec![1., 1., 0., -1.]; 5], vec![]),
        ("censored_normal", vec![vec![1., 1., 0., 0.]; 5], vec![]),
        ("censored_normal", vec![vec![1., 1., 1., 0., 1.]; 6], vec![]),
        ("censored_predict", vec![vec![0., 0.5]], vec![1., 0., 0.]),
    ] {
        assert!(
            execute(op, &d, &args).is_err(),
            "accepted malformed {op} {d:?}"
        );
    }
}
#[test]
fn public_c_abi_owns_censored_json_and_rejects_bad_support() {
    // A dispatch-only direct Rust API cannot satisfy the public ABI contract.
    let op = b"censored_normal";
    let data = [
        -2., -2., 0., 1., -1., -1., 0., 1., 0., 0., 0., 1., 1., 1., 0., 1., 2., 2., 0., 1.,
    ];
    unsafe {
        let mut n = 0;
        let p = ds_analytics_run(
            op.as_ptr(),
            op.len(),
            data.as_ptr(),
            5,
            4,
            std::ptr::null(),
            0,
            &mut n,
        );
        assert!(!p.is_null());
        let j = std::str::from_utf8(std::slice::from_raw_parts(p, n)).unwrap();
        near(values(j, "sigma")[0], 2_f64.sqrt(), 1e-6);
        ds_analytics_free(p, n);
        let bad_op = b"censored_beta";
        let bad = [0., 0., 0., 1., 0.2, 0.2, 0., 1., 0.4, 0.4, 0., 1.];
        n = 999;
        let bad_p = ds_analytics_run(
            bad_op.as_ptr(),
            bad_op.len(),
            bad.as_ptr(),
            3,
            4,
            std::ptr::null(),
            0,
            &mut n,
        );
        assert!(bad_p.is_null());
        assert_eq!(n, 0);
    }
}

fn saved_constant(family: f64, eta: f64, log_sigma: f64) -> Vec<f64> {
    let d = if family == 1. { 1 } else { 2 };
    let mut args = vec![1., family, 0., 0.95, 0., 1., 0., eta, log_sigma];
    args.extend(vec![0.; d * d]);
    args
}
#[test]
fn gamma_and_beta_quantiles_preserve_representable_extreme_left_tail() {
    // Gamma(shape1,scale1) has quantile -log1p(-p); uniform beta has quantile p.
    for (family, ls) in [(2., 0.), (5., -std::f64::consts::LN_2)] {
        let j = execute(
            "censored_predict",
            &[vec![0.5, 1e-100]],
            &saved_constant(family, 0., ls),
        )
        .unwrap();
        let q = values(&j, "quantiles")[0];
        assert!(
            (q - 1e-100).abs() < 1e-107,
            "lost tiny representable quantile family{family}: {q}"
        );
    }
}
#[test]
fn subnormal_left_censor_limit_has_finite_exponential_log_likelihood() {
    // P(Y<=tiny) ~ tiny/mu; two failures sum4 plus one left tail imply mu~4/3.
    let tiny = f64::from_bits(1);
    let j = fit(
        "censored_exponential",
        &[
            vec![1., 1., 0., 1.],
            vec![3., 3., 0., 1.],
            vec![0., tiny, 1., 1.],
        ],
    );
    near(values(&j, "coefficients")[0], (4_f64 / 3.).ln(), 1e-6);
    assert!(values(&j, "log_likelihood")[0].is_finite());
}
#[test]
fn original_unit_information_overflow_never_serializes_nonfinite_json() {
    // Coefficients fit in doubles while slope covariance does not; saved normalized
    // prediction remains useful and the unrepresentable original inference is null.
    let d = vec![
        vec![-2e-150, -1.8e150, -1.8e150, 0., 1.],
        vec![-1e-150, -1.1e150, -1.1e150, 0., 1.],
        vec![0., 0.1e150, 0.1e150, 0., 1.],
        vec![1e-150, 0.8e150, 0.8e150, 0., 1.],
        vec![2e-150, 2.1e150, 2.1e150, 0., 1.],
    ];
    // Data magnitude public bound is1e150; rescale responses by1/4.
    let d = d
        .into_iter()
        .map(|mut r| {
            r[1] *= 0.25;
            r[2] *= 0.25;
            r
        })
        .collect::<Vec<_>>();
    let j = fit("censored_normal", &d);
    assert!(
        !j.contains(":inf") && !j.contains("[inf") && !j.contains("NaN"),
        "{j}"
    );
    assert!(j.contains(r#""covariance":null"#), "{j}");
}

#[test]
fn mixed_censoring_predictor_fits_match_independent_scipy_likelihood_and_joint_information() {
    {
        let data = vec![
            vec![
                -1.0,
                0.2,
                -0.6600000000000001,
                -0.6600000000000001,
                0.0,
                2.0,
            ],
            vec![-0.8, -0.5, 0.85, 0.85, 0.0, 1.0],
            vec![
                -0.5,
                0.8,
                0.009999999999999953,
                0.009999999999999953,
                0.0,
                1.0,
            ],
            vec![-0.2, -0.8, 1.6400000000000001, 1.6400000000000001, 0.0, 1.0],
            vec![
                0.0,
                0.5,
                -0.050000000000000044,
                -0.050000000000000044,
                0.0,
                2.0,
            ],
            vec![0.2, -0.2, 1.16, 1.16, 0.0, 1.0],
            vec![0.4, 0.9, 0.3299999999999999, 0.3299999999999999, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 1.97, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 0.51, 1.0, 2.0],
            vec![1.0, -0.3, 1.49, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 1.0199999999999998, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 1.9619999999999997, 2.5069999999999997, 3.0, 1.0],
        ];
        let j = fit("censored_normal", &data);
        let expected = &[0.6078385544141092, 0.7067923557756827, -0.9932154872377998];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], -9.442399695974528, 2e-8);
        let se = &[
            0.10827559556464049,
            0.14102843459652833,
            0.19116724024586254,
            0.23257444892142642,
        ];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
    {
        let data = vec![
            vec![
                -1.0,
                0.2,
                0.44932896411722156,
                0.44932896411722156,
                0.0,
                2.0,
            ],
            vec![-0.8, -0.5, 1.9059870292719228, 1.9059870292719228, 0.0, 1.0],
            vec![-0.5, 0.8, 0.8394570207692074, 0.8394570207692074, 0.0, 1.0],
            vec![-0.2, -0.8, 3.781043387568781, 3.781043387568781, 0.0, 1.0],
            vec![0.0, 0.5, 0.7225273536420722, 0.7225273536420722, 0.0, 2.0],
            vec![0.2, -0.2, 2.270499837532406, 2.270499837532406, 0.0, 1.0],
            vec![0.4, 0.9, 1.0151130646157192, 1.0151130646157192, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 4.641325529953961, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 1.1107106103557052, 1.0, 2.0],
            vec![1.0, -0.3, 2.7870954605658507, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 1.7682670514337353, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 4.639652561011212, 5.928444939069881, 3.0, 1.0],
        ];
        let j = fit("censored_exponential", &data);
        let expected = &[0.39448975472281117, 0.6675444677955815, -0.8378779924293072];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], -16.078075225891638, 2e-8);
        let se = &[0.28071557435968425, 0.3769355549145644, 0.4647211222327847];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
    {
        let data = vec![
            vec![
                -1.0,
                0.2,
                0.44932896411722156,
                0.44932896411722156,
                0.0,
                2.0,
            ],
            vec![-0.8, -0.5, 1.9059870292719228, 1.9059870292719228, 0.0, 1.0],
            vec![-0.5, 0.8, 0.8394570207692074, 0.8394570207692074, 0.0, 1.0],
            vec![-0.2, -0.8, 3.781043387568781, 3.781043387568781, 0.0, 1.0],
            vec![0.0, 0.5, 0.7225273536420722, 0.7225273536420722, 0.0, 2.0],
            vec![0.2, -0.2, 2.270499837532406, 2.270499837532406, 0.0, 1.0],
            vec![0.4, 0.9, 1.0151130646157192, 1.0151130646157192, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 4.641325529953961, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 1.1107106103557052, 1.0, 2.0],
            vec![1.0, -0.3, 2.7870954605658507, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 1.7682670514337353, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 4.639652561011212, 5.928444939069881, 3.0, 1.0],
        ];
        let j = fit("censored_gamma", &data);
        let expected = &[
            0.38808168630377277,
            0.45356970676493197,
            -1.0232057603274742,
        ];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], -9.989627123291076, 2e-8);
        let se = &[
            0.09525779890624146,
            0.10568121695755883,
            0.15187232612787965,
            0.4760498640957501,
        ];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
    {
        let data = vec![
            vec![
                -1.0,
                0.2,
                0.44932896411722156,
                0.44932896411722156,
                0.0,
                2.0,
            ],
            vec![-0.8, -0.5, 1.9059870292719228, 1.9059870292719228, 0.0, 1.0],
            vec![-0.5, 0.8, 0.8394570207692074, 0.8394570207692074, 0.0, 1.0],
            vec![-0.2, -0.8, 3.781043387568781, 3.781043387568781, 0.0, 1.0],
            vec![0.0, 0.5, 0.7225273536420722, 0.7225273536420722, 0.0, 2.0],
            vec![0.2, -0.2, 2.270499837532406, 2.270499837532406, 0.0, 1.0],
            vec![0.4, 0.9, 1.0151130646157192, 1.0151130646157192, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 4.641325529953961, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 1.1107106103557052, 1.0, 2.0],
            vec![1.0, -0.3, 2.7870954605658507, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 1.7682670514337353, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 4.639652561011212, 5.928444939069881, 3.0, 1.0],
        ];
        let j = fit("censored_weibull", &data);
        let expected = &[0.47251284514853376, 0.5214893518503961, -0.8789665204703769];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], -10.093228626612703, 2e-8);
        let se = &[
            0.09534288957200678,
            0.13687772100810885,
            0.13624098002557175,
            0.25678919369673336,
        ];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
    {
        let data = vec![
            vec![
                -1.0,
                0.2,
                0.44932896411722156,
                0.44932896411722156,
                0.0,
                2.0,
            ],
            vec![-0.8, -0.5, 1.9059870292719228, 1.9059870292719228, 0.0, 1.0],
            vec![-0.5, 0.8, 0.8394570207692074, 0.8394570207692074, 0.0, 1.0],
            vec![-0.2, -0.8, 3.781043387568781, 3.781043387568781, 0.0, 1.0],
            vec![0.0, 0.5, 0.7225273536420722, 0.7225273536420722, 0.0, 2.0],
            vec![0.2, -0.2, 2.270499837532406, 2.270499837532406, 0.0, 1.0],
            vec![0.4, 0.9, 1.0151130646157192, 1.0151130646157192, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 4.641325529953961, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 1.1107106103557052, 1.0, 2.0],
            vec![1.0, -0.3, 2.7870954605658507, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 1.7682670514337353, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 4.639652561011212, 5.928444939069881, 3.0, 1.0],
        ];
        let j = fit("censored_lognormal", &data);
        let expected = &[0.30584584579937696, 0.5518078414365096, -0.9389735528965748];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], -10.555953616264873, 2e-8);
        let se = &[
            0.1071387485661854,
            0.13833129734177527,
            0.18888253455828252,
            0.23244472167024213,
        ];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
    {
        let data = vec![
            vec![-1.0, 0.2, 0.1570954688854527, 0.1570954688854527, 0.0, 2.0],
            vec![-0.8, -0.5, 0.47502081252106, 0.47502081252106, 0.0, 1.0],
            vec![
                -0.5,
                0.8,
                0.25540308417252405,
                0.25540308417252405,
                0.0,
                1.0,
            ],
            vec![-0.2, -0.8, 0.6726070170677605, 0.6726070170677605, 0.0, 1.0],
            vec![0.0, 0.5, 0.24973989440488234, 0.24973989440488234, 0.0, 2.0],
            vec![0.2, -0.2, 0.5448788923735801, 0.5448788923735801, 0.0, 1.0],
            vec![0.4, 0.9, 0.31864626621097447, 0.31864626621097447, 0.0, 1.0],
            vec![0.6, -0.9, 0.0, 0.7426905453026857, 1.0, 1.0],
            vec![0.8, 0.3, 0.0, 0.3728522336868044, 1.0, 2.0],
            vec![1.0, -0.3, 0.6271477663131956, 0.0, 2.0, 1.0],
            vec![1.2, 0.6, 0.4900013331200346, 0.0, 2.0, 1.0],
            vec![1.4, -0.6, 0.7255640142690734, 0.8255640142690734, 3.0, 1.0],
        ];
        let j = fit("censored_beta", &data);
        let expected = &[
            -0.37488813532998194,
            0.6558596289759091,
            -1.1135842191261383,
        ];
        let actual = values(&j, "coefficients");
        assert_eq!(actual.len(), expected.len());
        for (a, b) in actual.iter().zip(expected) {
            near(*a, *b, 2e-5);
        }
        near(values(&j, "log_likelihood")[0], 5.2281821492160745, 2e-8);
        let se = &[
            0.10113875034913902,
            0.1329743498644642,
            0.19075677971102553,
            0.46358494027223784,
        ];
        let actual = values(&j, "standard_errors");
        assert_eq!(actual.len(), se.len());
        for (a, b) in actual.iter().zip(se) {
            near(*a, *b, 3e-4);
        }
        assert!(j.contains("\"converged\":true"), "{j}");
    }
}

#[test]
fn unbounded_one_direction_censor_fits_are_rejected_before_claiming_inference() {
    for kind in [1., 2.] {
        let data = (1..5)
            .map(|i| vec![i as f64, i as f64, kind, 1.])
            .collect::<Vec<_>>();
        for family in [
            "normal",
            "exponential",
            "gamma",
            "weibull",
            "lognormal",
            "beta",
        ] {
            let d = if family == "beta" {
                data.iter()
                    .map(|r| vec![r[0] / 5., r[1] / 5., r[2], r[3]])
                    .collect::<Vec<_>>()
            } else {
                data.clone()
            };
            assert!(
                execute(&format!("censored_{family}"), &d, &[]).is_err(),
                "accepted unbounded {family} direction{kind}"
            );
        }
    }
}

#[test]
fn replicated_single_exponential_value_is_an_identifiable_frequency_fit() {
    let j = fit("censored_exponential", &[vec![2., 2., 0., 9.]]);
    near(values(&j, "coefficients")[0], 2_f64.ln(), 1e-7);
    near(values(&j, "standard_errors")[0], 1. / 3., 1e-5);
}
#[test]
fn common_censor_interval_without_failures_rejects_zero_dispersion_supremum() {
    for family in ["normal", "gamma", "weibull", "lognormal", "beta"] {
        let data = vec![
            vec![0.2, 0.8, 3., 2.],
            vec![0.3, 0.7, 3., 2.],
            vec![0.4, 0.6, 3., 2.],
        ];
        assert!(
            execute(&format!("censored_{family}"), &data, &[]).is_err(),
            "accepted variance-collapse {family}"
        );
    }
}

#[test]
fn prediction_budget_counts_gamma_special_function_iterations_globally() {
    // Many expensive but valid shape500000 quantiles must share one call budget.
    let model = saved_constant(2., 500000_f64.ln(), 0.);
    let data = vec![vec![500000., 0.5]; 512];
    let error = execute("censored_predict", &data, &model)
        .expect_err("unbounded repeated CDF work accepted");
    assert!(error.contains("work"), "{error}");
}
#[test]
fn nonconverged_coefficient_fit_withholds_observed_information_inference() {
    let data = vec![
        vec![-1., -1., -1., 0., 1.],
        vec![0., 1., 1., 0., 1.],
        vec![1., 1., 1., 0., 1.],
        vec![2., 0., 0., 0., 1.],
    ];
    let j = execute("censored_normal", &data, &[0., 1e-7, 0.95]).unwrap();
    assert!(j.contains(r#""converged":false"#), "{j}");
    assert!(j.contains(r#""standard_errors":null"#), "{j}");
    let pred = execute(
        "censored_predict",
        &[vec![0., 0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    assert!(pred.contains(r#""mean_standard_errors":[null]"#), "{pred}");
}

#[test]
fn adjacent_distinct_raw_endpoints_keep_positive_interval_mass_all_families() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    for family in [
        "normal",
        "exponential",
        "gamma",
        "weibull",
        "lognormal",
        "beta",
    ] {
        let mut d = if family == "beta" {
            exact(&[0.2, 0.5, 0.8])
        } else {
            exact(&[1., 2., 3.])
        };
        d.push(vec![lo, hi, 3., 1.]);
        let j = fit(&format!("censored_{family}"), &d);
        assert!(values(&j, "log_likelihood")[0].is_finite(), "{family}: {j}");
        if family == "exponential" {
            let mu = 1.525_f64;
            let expected = -4. * mu.ln() - 6.1 / mu + (hi - lo).ln();
            near(values(&j, "log_likelihood")[0], expected, 2e-8);
        }
    }
}

#[test]
fn adjacent_interval_normal_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[-2.0, -1.0, 0.0, 1.0, 2.0]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute("censored_normal", &data, &[0.0, 1e-07, 0.95, 0.0, 0.0]).unwrap();
    near(values(&j, "log_likelihood")[0], -49.33487331058497, 2e-11);
}

#[test]
fn adjacent_interval_exponential_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[0.4, 0.8, 1.2, 2., 3.]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute("censored_exponential", &data, &[0.0, 1e-07, 0.95, 0.0]).unwrap();
    near(values(&j, "log_likelihood")[0], -46.316242111356935, 2e-11);
}

#[test]
fn adjacent_interval_gamma_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[0.4, 0.8, 1.2, 2., 3.]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute("censored_gamma", &data, &[0.0, 1e-07, 0.95, 0.0, 0.0]).unwrap();
    near(values(&j, "log_likelihood")[0], -46.316242111356935, 2e-11);
}

#[test]
fn adjacent_interval_weibull_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[0.4, 0.8, 1.2, 2., 3.]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute("censored_weibull", &data, &[0.0, 1e-07, 0.95, 0.0, 0.0]).unwrap();
    near(values(&j, "log_likelihood")[0], -46.316242111356935, 2e-11);
}

#[test]
fn adjacent_interval_lognormal_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[0.4, 0.8, 1.2, 2., 3.]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute("censored_lognormal", &data, &[0.0, 1e-07, 0.95, 0.0, 0.0]).unwrap();
    near(values(&j, "log_likelihood")[0], -46.81789645297084, 2e-11);
}

#[test]
fn adjacent_interval_beta_fixed_likelihood_matches_independent_log_density_integral() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[0.15, 0.25, 0.4, 0.55, 0.75]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute(
        "censored_beta",
        &data,
        &[0.0, 1e-07, 0.95, 0.0, -std::f64::consts::LN_2],
    )
    .unwrap();
    near(values(&j, "log_likelihood")[0], -38.816242111356935, 2e-11);
}

#[test]
fn normal_centering_preserves_original_adjacent_interval_width() {
    let lo = 0.1_f64;
    let hi = f64::from_bits(lo.to_bits() + 1);
    let mut data = exact(&[1e16, 1e16 + 2., 1e16 + 4., 1e16 + 6., 1e16 + 8.]);
    data.push(vec![lo, hi, 3., 1.]);
    let j = execute(
        "censored_normal",
        &data,
        &[0., 1e-7, 0.95, 1e16, (3e15_f64).ln()],
    )
    .unwrap();
    near(values(&j, "log_likelihood")[0], -263.7097609676133, 2e-11);
}

#[test]
fn saturated_positive_beta_logit_preserves_small_shape_and_mean_derivative() {
    let model = [1., 5., 0., 0.95, 0., 1., 1., 40., 0., 1., 0., 0., 1.];
    let j = execute("censored_predict", &[vec![0.5, 1e-20]], &model).unwrap();
    for (key, want) in [
        ("cdf", 2.9447347740752107e-18),
        ("density", 8.496708510583176e-18),
        ("mean_standard_errors", 4.248354255291589e-18),
    ] {
        let got = values(&j, key)[0];
        assert!((got - want).abs() < want * 1e-8, "{key}: {got} != {want}");
    }
    near(values(&j, "quantiles")[0], 0.0023510845295342986, 1e-10);
}

#[test]
fn positive_distribution_quantiles_reject_unrepresentable_zero() {
    for (family, eta) in [(1., -740.), (3., -740.), (4., -730.)] {
        let model = saved_constant(family, eta, 0.);
        assert!(
            execute("censored_predict", &[vec![0., 1e-100]], &model).is_err(),
            "silently accepted zero quantile family{family}"
        );
    }
}
#[test]
fn subnormal_exponential_interval_uses_actual_lower_hazard() {
    let tiny = f64::from_bits(1);
    let mut data = exact(&[2. * tiny, 4. * tiny, 8. * tiny]);
    data.push(vec![tiny, 2. * tiny, 3., 1.]);
    let eta = tiny.ln();
    let expected = -3. * eta - 14. - 1. + (-(-1_f64).exp()).ln_1p();
    let j = execute("censored_exponential", &data, &[0., 1e-7, 0.95, eta]).unwrap();
    near(values(&j, "log_likelihood")[0], expected, 2e-11);
}

#[test]
fn beta_logit32_small_cdf_and_quantile_keep_complement_tail_precision() {
    let model = [1., 5., 0., 0.95, 0., 1., 1., 32., 0., 1., 0., 0., 1.];
    let j = execute("censored_predict", &[vec![0.9, 1e-16]], &model).unwrap();
    for (key, want) in [
        ("cdf", 2.916031880855248e-14),
        ("mean_standard_errors", 1.2664165549093855e-14),
    ] {
        let got = values(&j, key)[0];
        assert!((got - want).abs() < want * 1e-8, "{key}: {got} != {want}");
    }
    near(values(&j, "quantiles")[0], 0.007865202168769292, 1e-10);
}

#[test]
fn gamma_tiny_shape_upper_tail_is_computed_without_large_cdf_subtraction() {
    let model = saved_constant(2., -32., 0.);
    let data = vec![
        vec![0.1, 1. - 1e-14],
        vec![1e-100, 1. - 1e-14],
        vec![1e-300, 1. - 1e-14],
    ];
    let j = execute("censored_predict", &data, &model).unwrap();
    let expected = [
        2.308581079283308e-14,
        2.9087219261132264e-12,
        8.740785687789858e-12,
    ];
    let actual = values(&j, "survival");
    assert_eq!(actual.len(), expected.len());
    for (a, b) in actual.iter().zip(expected) {
        assert!((a - b).abs() < b * 1e-8, "{a} != {b}");
    }
    for q in values(&j, "quantiles") {
        near(q, 0.3526041408838702, 1e-9);
    }
}
#[test]
fn representable_tiny_normal_scale_and_prediction_se_do_not_underflow_when_squaring() {
    let d = exact(&[-2e-200, -1e-200, 0., 1e-200, 2e-200]);
    let j = fit("censored_normal", &d);
    let sigma = values(&j, "sigma")[0];
    assert!((sigma / (2_f64.sqrt() * 1e-200) - 1.).abs() < 1e-7, "{j}");
    assert!(
        j.contains(r#""covariance":null"#),
        "unrepresentable variance mustnotbezero precision: {j}"
    );
    let pred = execute(
        "censored_predict",
        &[vec![0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    let se = values(&pred, "mean_standard_errors")[0];
    assert!(
        (se / ((2_f64 / 5.).sqrt() * 1e-200) - 1.).abs() < 1e-5,
        "{pred}"
    );
}
#[test]
fn representable_tiny_feature_units_keep_full_rank_saved_predictions() {
    let d = vec![
        vec![-2e-200, -1.8, -1.8, 0., 1.],
        vec![-1e-200, -1.1, -1.1, 0., 1.],
        vec![0., 0.1, 0.1, 0., 1.],
        vec![1e-200, 0.8, 0.8, 0., 1.],
        vec![2e-200, 2.1, 2.1, 0., 1.],
    ];
    let j = fit("censored_normal", &d);
    let pred = execute(
        "censored_predict",
        &[vec![1e-200, 0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    near(values(&pred, "means")[0], 0.99, 1e-6);
}

#[test]
fn finite_mean_and_se_with_unrepresentable_wald_endpoints_withhold_only_interval() {
    for (model, query) in [
        (
            vec![1., 1., 0., 0.95, 0., 1., 1., 709., 0., 1.],
            vec![0., 0.1],
        ),
        (
            vec![1., 0., 0., 0.95, 1e308, 1e308, 1., 0., 0., 1., 0., 0., 1.],
            vec![0., 0.5],
        ),
    ] {
        let j = execute("censored_predict", &[query], &model).unwrap();
        assert!(values(&j, "means")[0].is_finite());
        assert!(values(&j, "mean_standard_errors")[0].is_finite());
        assert!(
            j.contains(r#""mean_confidence_intervals":[null]"#),
            "nonfiniteintervalmustnotserialize: {j}"
        );
    }
}
#[test]
fn huge_finite_asymmetric_saved_covariance_is_rejected_without_tolerance_overflow() {
    let model = [
        1., 0., 0., 0.95, 0., 1., 1., 0., 0., 1.7e308, 1.1e308, 1e308, 1.7e308,
    ];
    assert!(
        execute("censored_predict", &[vec![0., 0.5]], &model).is_err(),
        "overflowedsymmetrytoleranceacceptedmalformedcovariance"
    );
}
