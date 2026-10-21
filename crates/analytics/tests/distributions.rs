use ducksmiles_analytics::distributions::run as execute;

fn num(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).unwrap() + key.len() + 3;
    json[start..]
        .trim_start()
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}

fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).unwrap() + key.len() + 4;
    json[start..]
        .split(']')
        .next()
        .unwrap()
        .split(',')
        .map(|v| v.parse().unwrap())
        .collect()
}

fn near(a: f64, b: f64, tol: f64) {
    assert!((a - b).abs() <= tol * b.abs().max(1.0), "{a} != {b}");
}

fn sample() -> Vec<Vec<f64>> {
    [0.5, 1.0, 1.5, 2.5, 4.0, 6.0]
        .iter()
        .map(|v| vec![*v])
        .collect()
}

#[test]
fn continuous_mle_fits_match_independent_scipy_location_zero_fixtures() {
    // A sample-SD estimator or a fitted free location would break these MLEs.
    for (op, params, ll) in [
        (
            "fit_normal",
            vec![("mean", 2.5833333333333335), ("sigma", 1.9021187017522212)],
            -12.371441425760134,
        ),
        (
            "fit_lognormal",
            vec![
                ("log_mean", 0.6344437482950532),
                ("log_sigma", 0.836813848860365),
            ],
            -11.251371872859794,
        ),
        (
            "fit_exponential",
            vec![("rate", 0.3870967741935484)],
            -11.694483328182875,
        ),
        (
            "fit_gamma",
            vec![("shape", 1.7371600335026973), ("scale", 1.4871015240458112)],
            -11.22788015721257,
        ),
        (
            "fit_weibull",
            vec![("shape", 1.3833562354097264), ("scale", 2.839505353472658)],
            -11.255887884495742,
        ),
    ] {
        let out = execute(op, &sample(), &[]).unwrap();
        near(num(&out, "version"), 1.0, 0.0);
        for (key, value) in params {
            near(num(&out, key), value, 2e-8);
        }
        near(num(&out, "log_likelihood"), ll, 2e-10);
        assert!(out.contains("\"converged\":true"));
        let count = if op == "fit_exponential" { 1.0 } else { 2.0 };
        near(num(&out, "aic"), 2.0 * count - 2.0 * ll, 2e-10);
        near(num(&out, "bic"), count * 6.0_f64.ln() - 2.0 * ll, 2e-10);
    }
}

#[test]
fn poisson_mle_includes_log_factorial_and_identifies_zero_boundary() {
    let rows: Vec<_> = [0., 1., 2., 2., 3., 7.].iter().map(|x| vec![*x]).collect();
    let out = execute("fit_poisson", &rows, &[]).unwrap();
    near(num(&out, "rate"), 2.5, 0.0);
    near(num(&out, "log_likelihood"), -12.958854213301034, 2e-12);
    let out = execute("fit_poisson", &[vec![0.], vec![0.]], &[]).unwrap();
    near(num(&out, "rate"), 0.0, 0.0);
    near(num(&out, "log_likelihood"), 0.0, 0.0);
}

#[test]
fn right_censored_lifetime_fits_include_surviving_units_exposure() {
    // Ignoring censors or using their times as exact events changes both fits.
    let rows: Vec<_> = sample()
        .iter()
        .zip([1., 0., 1., 1., 0., 1.])
        .map(|(r, e)| vec![r[0], e])
        .collect();
    let out = execute("survival_exponential", &rows, &[]).unwrap();
    near(num(&out, "rate"), 0.25806451612903225, 2e-12);
    near(num(&out, "rate_standard_error"), 0.12903225806451613, 2e-12);
    near(num(&out, "log_likelihood"), -9.41818265122124, 2e-12);
    let out = execute("survival_weibull", &rows, &[]).unwrap();
    near(num(&out, "shape"), 1.34901382272634, 2e-8);
    near(num(&out, "scale"), 3.8051173489356405, 2e-8);
    near(num(&out, "log_likelihood"), -9.168764162643598, 2e-10);
}

#[test]
fn jarque_bera_labels_asymptotic_probability_and_uses_population_moments() {
    let out = execute("normality_test", &sample(), &[]).unwrap();
    near(num(&out, "statistic"), 0.6664994693339898, 2e-12);
    near(num(&out, "p_value"), 0.7165912141396371, 2e-12);
    assert!(out.contains("jarque_bera"));
    assert!(out.contains("asymptotic"));
    assert!(out.contains("\"large_sample_recommended\":true"));
}

#[test]
fn gaussian_kde_is_average_of_normalized_kernel_densities() {
    let out = execute("kde_gaussian", &[vec![0.], vec![1.]], &[0.5, -1., 0., 1.]).unwrap();
    let density = array(&out, "densities");
    assert_eq!(density.len(), 3);
    for (a, b) in density.iter().zip([
        0.05412479673895295,
        0.452_933_246_914_620_7,
        0.452_933_246_914_620_7,
    ]) {
        near(*a, b, 2e-12);
    }
}

#[test]
fn kde_large_bandwidth_preserves_representable_subnormal_density() {
    let out = execute("kde_gaussian", &[vec![0.]], &[1e308, 0.]).unwrap();
    near(
        array(&out, "densities")[0] / 3.989422804014326e-309,
        1.0,
        2e-12,
    );
}

#[test]
fn kde_small_bandwidth_normalizes_before_large_density_multiplication() {
    let out = execute("kde_gaussian", &vec![vec![0.]; 10], &[1e-308, 0.]).unwrap();
    near(
        array(&out, "densities")[0] / 3.989422804014327e307,
        1.0,
        2e-12,
    );
}

#[test]
fn box_cox_handles_log_limit_and_estimates_profile_likelihood_lambda() {
    let out = execute("box_cox", &sample(), &[0.]).unwrap();
    for (a, b) in array(&out, "transformed").iter().zip([
        -std::f64::consts::LN_2,
        0.,
        0.4054651081081644,
        0.9162907318741551,
        1.3862943611198906,
        1.791759469228055,
    ]) {
        near(*a, b, 2e-12);
    }
    let out = execute("box_cox", &sample(), &[]).unwrap();
    near(num(&out, "lambda"), 0.12985812816832185, 2e-7);
    near(
        num(&out, "profile_log_likelihood"),
        -2.709549779362105,
        2e-10,
    );
}

#[test]
fn logrank_tied_events_use_hypergeometric_variance_with_same_time_censors_at_risk() {
    let rows = vec![
        vec![1., 1., 0.],
        vec![2., 1., 0.],
        vec![2., 0., 0.],
        vec![2., 1., 1.],
        vec![3., 1., 1.],
        vec![4., 0., 1.],
    ];
    let out = execute("logrank_test", &rows, &[]).unwrap();
    near(num(&out, "observed_group0"), 2., 0.);
    near(num(&out, "expected_group0"), 1.3, 2e-12);
    near(num(&out, "variance"), 0.61, 2e-12);
    near(num(&out, "statistic"), 0.8032786885245902, 2e-12);
    near(num(&out, "p_value"), 0.3701148999172469, 2e-12);
    let mut reversed = rows;
    reversed.reverse();
    let swapped: Vec<_> = reversed
        .iter()
        .map(|r| vec![r[0], r[1], 1. - r[2]])
        .collect();
    let again = execute("logrank_test", &swapped, &[]).unwrap();
    near(num(&again, "statistic"), num(&out, "statistic"), 2e-12);
}

#[test]
fn distribution_contract_rejects_bad_domains_and_unidentified_lifetime_models() {
    for op in [
        "fit_normal",
        "fit_lognormal",
        "fit_exponential",
        "fit_poisson",
        "fit_gamma",
        "fit_weibull",
        "normality_test",
        "box_cox",
    ] {
        assert!(execute(op, &[], &[]).is_err(), "{op}");
        assert!(execute(op, &[vec![f64::NAN]], &[]).is_err(), "{op}");
        assert!(
            execute(op, &[vec![1., 2.], vec![2., 3.]], &[]).is_err(),
            "{op}"
        );
    }
    for op in [
        "fit_normal",
        "fit_lognormal",
        "fit_gamma",
        "fit_weibull",
        "normality_test",
        "box_cox",
    ] {
        assert!(execute(op, &[vec![2.], vec![2.]], &[]).is_err(), "{op}");
    }
    for op in ["fit_lognormal", "fit_gamma", "fit_weibull", "box_cox"] {
        assert!(execute(op, &[vec![0.], vec![1.]], &[]).is_err(), "{op}");
    }
    assert!(execute("fit_poisson", &[vec![0.5], vec![2.]], &[]).is_err());
    assert!(execute("fit_exponential", &[vec![0.], vec![0.]], &[]).is_err());
    assert!(execute("kde_gaussian", &sample(), &[0., 1.]).is_err());
    for op in ["survival_exponential", "survival_weibull"] {
        assert!(
            execute(op, &[vec![1., 0.], vec![2., 0.]], &[]).is_err(),
            "{op}"
        );
        assert!(
            execute(op, &[vec![1., 2.], vec![2., 1.]], &[]).is_err(),
            "{op}"
        );
    }
    assert!(execute("logrank_test", &[vec![1., 0., 0.], vec![2., 0., 1.]], &[]).is_err());
}

#[test]
fn box_cox_estimated_lambda_is_invariant_to_measurement_units() {
    // The Jacobian must cancel the transformation's scale in profile fitting.
    let base = execute("box_cox", &sample(), &[]).unwrap();
    let large: Vec<_> = sample().iter().map(|r| vec![r[0] * 1e50]).collect();
    let scaled = execute("box_cox", &large, &[]).unwrap();
    near(num(&scaled, "lambda"), num(&base, "lambda"), 2e-6);
    near(
        num(&scaled, "profile_log_likelihood"),
        num(&base, "profile_log_likelihood") - 6.0 * 1e50_f64.ln(),
        2e-12,
    );
    let tiny: Vec<_> = sample().iter().map(|r| vec![r[0] * 1e-100]).collect();
    let scaled = execute("box_cox", &tiny, &[]).unwrap();
    near(num(&scaled, "lambda"), num(&base, "lambda"), 2e-6);
}

#[test]
fn poisson_large_equal_counts_preserve_small_log_likelihood() {
    // Stirling's convergent correction gives this independent oracle. Direct
    // subtraction of x*log(x), x and log-gamma loses these low-order digits.
    let out = execute("fit_poisson", &[vec![1e12], vec![1e12]], &[]).unwrap();
    near(num(&out, "log_likelihood"), -29.468898182338062, 2e-12);
}

#[test]
fn parametric_likelihoods_rescale_with_lifetime_units() {
    for op in ["fit_gamma", "fit_weibull"] {
        let base = execute(op, &sample(), &[]).unwrap();
        let small: Vec<_> = sample().iter().map(|r| vec![r[0] * 1e-100]).collect();
        let again = execute(op, &small, &[]).unwrap();
        near(num(&again, "shape"), num(&base, "shape"), 2e-9);
        near(num(&again, "scale") / 1e-100, num(&base, "scale"), 2e-9);
        near(
            num(&again, "log_likelihood"),
            num(&base, "log_likelihood") + 600.0 * 10.0_f64.ln(),
            2e-12,
        );
    }
}
