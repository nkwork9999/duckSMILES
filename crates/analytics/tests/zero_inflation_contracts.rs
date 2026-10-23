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
fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * (1. + b.abs()), "{a} != {b}");
}
fn carrier(family: usize, eta: f64, dispersion: f64, pi: f64) -> Vec<f64> {
    let d = if matches!(family, 1 | 3 | 4) { 3 } else { 2 };
    let boundary = pi == 0. || (matches!(family, 1 | 3) && dispersion == 0.);
    let mut a = vec![
        1.,
        family as f64,
        0.,
        0.95,
        0.,
        0.,
        f64::from(boundary),
        eta,
        dispersion,
        pi,
    ];
    a.resize(10 + d * d, 0.);
    a
}

#[test]
fn binomial_mixture_keeps_base_zeros_and_exact_finite_support_probabilities() {
    // Rational binomial n4,p2/5 mixed with zero probability3/10.
    let a = carrier(0, (2_f64 / 3.).ln(), 0., 0.3);
    let rows = (0..5).map(|y| vec![4., y as f64, 0.5]).collect::<Vec<_>>();
    let j = execute("zi_predict", &rows, &a).unwrap();
    for (v, w) in values(&j, "query_probability_mass")
        .into_iter()
        .zip([0.39072, 0.24192, 0.24192, 0.10752, 0.01792])
    {
        near(v, w, 1e-11);
    }
    near(values(&j, "cdf")[1], 0.63264, 1e-11);
    near(values(&j, "mixture_mean")[0], 1.12, 1e-12);
    near(values(&j, "mixture_variance")[0], 1.2096, 1e-12);
    assert_eq!(values(&j, "quantile"), vec![1.; 5]);
}

#[test]
fn beta_binomial_uses_jmp_intraclass_dispersion_and_recovers_binomial_boundary() {
    // Independent exact rising-factorial fractions,denominator15625.
    let rows = (0..5).map(|y| vec![4., y as f64, 0.5]).collect::<Vec<_>>();
    let j = execute(
        "zi_predict",
        &rows,
        &carrier(1, (2_f64 / 3.).ln(), 0.2, 0.3),
    )
    .unwrap();
    for (v, w) in values(&j, "query_probability_mass")
        .into_iter()
        .zip([0.461568, 0.191488, 0.169728, 0.119808, 0.057408])
    {
        near(v, w, 1e-10);
    }
    near(values(&j, "cdf")[1], 0.653056, 1e-10);
    near(values(&j, "mixture_variance")[0], 1.6128, 1e-11);
    let base = execute("zi_predict", &rows, &carrier(0, (2_f64 / 3.).ln(), 0., 0.3)).unwrap();
    let limit = execute("zi_predict", &rows, &carrier(1, (2_f64 / 3.).ln(), 0., 0.3)).unwrap();
    assert_eq!(
        values(&base, "query_probability_mass"),
        values(&limit, "query_probability_mass")
    );
}

#[test]
fn poisson_and_negative_binomial_mixtures_have_correct_zero_mass_and_variance() {
    let rows = vec![vec![1., 0., 0.5], vec![1., 1., 0.5], vec![1., 2., 0.5]];
    let p = execute("zi_predict", &rows, &carrier(2, 2_f64.ln(), 0., 0.3)).unwrap();
    near(values(&p, "zero_probability")[0], 0.3947346982656289, 1e-12);
    near(
        values(&p, "query_probability_mass")[1],
        0.18946939653125777,
        1e-12,
    );
    near(values(&p, "mixture_mean")[0], 1.4, 1e-12);
    near(values(&p, "mixture_variance")[0], 2.24, 1e-12);
    let n = execute("zi_predict", &rows, &carrier(3, 2_f64.ln(), 0.5, 0.3)).unwrap();
    for (v, w) in values(&n, "query_probability_mass")
        .into_iter()
        .zip([0.475, 0.175, 0.13125])
    {
        near(v, w, 1e-12);
    }
    near(values(&n, "mixture_variance")[0], 3.64, 1e-12);
    near(values(&n, "cdf")[1], 0.65, 1e-12);
    let limit = execute("zi_predict", &rows, &carrier(3, 2_f64.ln(), 0., 0.3)).unwrap();
    assert_eq!(
        values(&p, "query_probability_mass"),
        values(&limit, "query_probability_mass")
    );
}

#[test]
fn gamma_has_separate_zero_atom_and_jmp_mean_over_sigma_shape() {
    let rows = vec![vec![1., 0., 0.2], vec![1., 1., 0.5], vec![1., -1., 0.9]];
    let j = execute("zi_predict", &rows, &carrier(4, 2_f64.ln(), 1., 0.3)).unwrap();
    assert_eq!(values(&j, "query_probability_mass"), vec![0.3, 0., 0.]);
    near(
        values(&j, "query_positive_density")[1],
        0.2575156088200096,
        1e-12,
    );
    near(values(&j, "cdf")[1], 0.48496878235998075, 1e-11);
    near(values(&j, "mixture_variance")[0], 2.24, 1e-12);
    assert_eq!(values(&j, "quantile")[0], 0.);
    assert_eq!(values(&j, "cdf")[2], 0.);
    assert_eq!(values(&j, "survival")[2], 1.);
}

#[test]
fn zip_intercept_mle_matches_independent_zero_truncated_count_equation() {
    let rows = vec![
        vec![0., 1., 10.],
        vec![1., 1., 3.],
        vec![2., 1., 4.],
        vec![3., 1., 2.],
        vec![4., 1., 1.],
    ];
    let j = execute("zi_poisson", &rows, &[]).unwrap();
    near(values(&j, "coefficients")[0], 0.5460013727818611, 2e-6);
    near(values(&j, "pi")[0], 0.39177549201032386, 2e-6);
    assert!(j.contains("\"converged\":true"), "{j}");
    let q = execute(
        "zi_predict",
        &[vec![1., 0., 0.5]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    near(values(&q, "mixture_mean")[0], 1.05, 2e-6);
    assert_eq!(values(&j, "covariance").len(), 4);
}

#[test]
fn finite_support_mles_recover_exact_rational_mixture_frequencies() {
    for (op, counts, disp) in [
        (
            "zi_binomial",
            vec![39072., 24192., 24192., 10752., 1792.],
            0.,
        ),
        (
            "zi_beta_binomial",
            vec![7212., 2992., 2652., 1872., 897.],
            0.2,
        ),
    ] {
        let rows = counts
            .iter()
            .enumerate()
            .map(|(y, w)| vec![y as f64, 4., *w])
            .collect::<Vec<_>>();
        let j = execute(op, &rows, &[]).unwrap();
        near(values(&j, "coefficients")[0], (2_f64 / 3.).ln(), 1e-5);
        near(values(&j, "pi")[0], 0.3, 1e-5);
        near(values(&j, "dispersion")[0], disp, 1e-5);
        assert!(j.contains("\"converged\":true"), "{op}: {j}");
    }
}

#[test]
fn gamma_and_negative_binomial_mles_match_independent_reference_fits() {
    let g = execute(
        "zi_gamma",
        &[
            vec![0., 1., 3.],
            vec![1., 1., 2.],
            vec![2., 1., 3.],
            vec![3., 1., 2.],
        ],
        &[],
    )
    .unwrap();
    near(values(&g, "coefficients")[0], 2_f64.ln(), 3e-6);
    near(values(&g, "dispersion")[0], 0.32025426596168105, 3e-6);
    near(values(&g, "pi")[0], 0.3, 3e-6);
    let counts = [
        47500., 17500., 13125., 8750., 5469., 3281., 1914., 1094., 615., 342., 188., 103., 56.,
        30., 16., 9., 5., 2., 1., 1.,
    ];
    let rows = counts
        .iter()
        .enumerate()
        .map(|(y, w)| vec![y as f64, 1., *w])
        .collect::<Vec<_>>();
    let n = execute("zi_negative_binomial", &rows, &[]).unwrap();
    near(
        values(&n, "coefficients")[0],
        2.000156667174194_f64.ln(),
        1e-5,
    );
    near(values(&n, "dispersion")[0], 0.49988967519193705, 1e-5);
    near(values(&n, "pi")[0], 0.30003183124750976, 1e-5);
}

#[test]
fn no_zero_poisson_fits_exact_pi_boundary_without_fabricated_joint_inference() {
    let j = execute(
        "zi_poisson",
        &[vec![1., 1., 1.], vec![2., 1., 1.], vec![3., 1., 1.]],
        &[],
    )
    .unwrap();
    near(values(&j, "coefficients")[0], 2_f64.ln(), 1e-6);
    assert_eq!(values(&j, "pi"), vec![0.]);
    assert!(j.contains("\"at_boundary\":true"), "{j}");
    assert!(j.contains("\"covariance\":null"), "{j}");
}

#[test]
fn unidentified_or_malformed_fits_and_carriers_are_rejected() {
    for op in [
        "zi_binomial",
        "zi_beta_binomial",
        "zi_poisson",
        "zi_negative_binomial",
        "zi_gamma",
    ] {
        assert!(
            execute(op, &[vec![0., 4., 10.]], &[]).is_err(),
            "{op}: all zero"
        );
        assert!(
            execute(op, &[vec![1., 4., -1.], vec![2., 4., 2.]], &[]).is_err(),
            "{op}: negative frequency"
        );
        assert!(
            execute(op, &[vec![1., 4., 0.5], vec![2., 4., 2.]], &[]).is_err(),
            "{op}: fractional frequency"
        );
    }
    assert!(execute("zi_binomial", &[vec![0., 1., 5.], vec![1., 1., 5.]], &[]).is_err());
    let a = carrier(2, 2_f64.ln(), 0., 0.3);
    let mut bad = a.clone();
    bad[9] = 1.;
    assert!(execute("zi_predict", &[vec![1., 1., 0.5]], &bad).is_err());
    let mut bad = a.clone();
    bad.push(0.);
    assert!(execute("zi_predict", &[vec![1., 1., 0.5]], &bad).is_err());
    assert!(execute("zi_predict", &[vec![1., 1., 0.]], &a).is_err());
    assert!(execute("zi_predict", &[vec![2., 1., 0.5]], &a).is_err());
}

#[test]
fn large_positive_logits_preserve_small_binomial_tails_and_variance() {
    // Direct mpmath80 finite-support sums at eta40. Computing1-logistic(40)
    // in binary64 loses this entire complementary probability.
    for (f, cdf, variance) in [
        (0, 8.303379160420425e-87, 2.5490125531749533e-17),
        (1, 1.9555916413246998e-19, 5.0980251063499066e-17),
    ] {
        let j = execute(
            "zi_predict",
            &[vec![6., 1., 0.5]],
            &carrier(f, 40., if f == 1 { 0.2 } else { 0. }, 0.),
        )
        .unwrap();
        let got = values(&j, "cdf")[0];
        assert!(
            (got / cdf - 1.).abs() < 1e-10,
            "family{f}:CDF {got} vs {cdf}"
        );
        let got = values(&j, "mixture_variance")[0];
        assert!(
            (got / variance - 1.).abs() < 1e-10,
            "family{f}:variance {got} vs {variance}"
        );
    }
}

#[test]
fn gamma_quantile_reports_unrepresentable_lower_tail_instead_of_iteration_truncation() {
    // Shape.01 and u1e-5 imply a positive quantile below the smallest f64.
    // A fixed240-step linear bisection returns an unrelated value near1e-73.
    let a = carrier(4, 0.01_f64.ln(), 1., 0.);
    assert!(execute("zi_predict", &[vec![1., 0.01, 1e-5]], &a).is_err());
}

#[test]
fn gamma_cdf_retains_representable_probability_when_scaled_query_underflows() {
    let a = carrier(4, 1e153_f64.ln(), 1e154, 0.);
    let j = execute("zi_predict", &[vec![1., f64::from_bits(1), 0.5]], &a).unwrap();
    let cdf = values(&j, "cdf")[0];
    assert!(cdf > 0., "scaled-positive Gamma CDF was lost: {j}");
    // Independent mpmath80 regularized integral for these binary64 inputs.
    assert!((cdf / 1.9545095271844815e-48 - 1.).abs() < 1e-10);
}

#[test]
fn gamma_quantile_preserves_tiny_representable_probabilities() {
    // Gamma(shape1,scale1) is exponential:Q(u)=-log1p(-u).
    let j = execute(
        "zi_predict",
        &[vec![1., 1., 1e-250]],
        &carrier(4, 0., 1., 0.),
    )
    .unwrap();
    let q = values(&j, "quantile")[0];
    assert!((q / 1e-250 - 1.).abs() < 1e-10, "{q}");
}

#[test]
fn mixture_quantile_just_above_zero_atom_uses_conditional_tail_without_cancellation() {
    let pi = 0.3_f64;
    let u = f64::from_bits(pi.to_bits() + 1);
    let j = execute("zi_predict", &[vec![1., 1., u]], &carrier(4, 0., 1., pi)).unwrap();
    let want = -(-((u - pi) / (1. - pi))).ln_1p();
    let got = values(&j, "quantile")[0];
    assert!((got / want - 1.).abs() < 1e-10, "{got} vs {want}");
}

#[test]
fn gamma_normalization_is_accurate_for_small_and_large_shapes() {
    // Independent mpmath75 integral/loggamma; the large-shape integral uses
    // x=a+sqrt(a)*t and integrates t from-30 to the supplied query.
    for (eta, y, cdf, sf, pdf) in [
        (
            -17.72753356339242,
            0.1,
            0.9999999635415212,
            3.6458478795160347e-8,
            1.809674773624767e-7,
        ),
        (
            -13.815510557964274,
            2.078255658098977e-44,
            0.9999,
            9.999999999998897e-5,
            4.811246374349484e37,
        ),
        (
            16.012735135300492,
            9000000.,
            0.5000443269220536,
            0.49995567307794647,
            0.00013298075890250745,
        ),
    ] {
        let j = execute(
            "zi_predict",
            &[vec![1., y, 0.999999999999]],
            &carrier(4, eta, 1., 0.),
        )
        .unwrap();
        for (key, want) in [
            ("cdf", cdf),
            ("survival", sf),
            ("query_positive_density", pdf),
        ] {
            let got = values(&j, key)[0];
            assert!(
                (got / want - 1.).abs() < 5e-10,
                "eta{eta} {key}: {got} vs {want}"
            );
        }
    }
}
