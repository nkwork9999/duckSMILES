use ducksmiles_analytics::execute;

fn values(s: &str, key: &str) -> Vec<f64> {
    let tail = s.split_once(&format!("\"{key}\":")).expect(key).1;
    let end = if tail.starts_with('[') {
        let mut depth = 0;
        tail.char_indices()
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
        tail.find([',', '}']).unwrap()
    };
    tail[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().unwrap())
        .collect()
}
fn near(actual: f64, expected: f64, tol: f64) {
    assert!(
        (actual - expected).abs() <= tol * (1. + expected.abs()),
        "{actual} != {expected}"
    );
}

// Wrong residual degrees of freedom or relabeling supplied-covariance GLS as an
// estimated model breaks these closed-form independent-normal likelihoods.
#[test]
fn spatial_ml_and_reml_estimate_variance_with_the_correct_residual_df() {
    let data = vec![vec![0., 0., 1.], vec![1., 0., 2.], vec![2., 0., 6.]];
    for (reml, want) in [(0., 14. / 3.), (1., 7.)] {
        let j = execute(
            "spatial_fit",
            &data,
            &[1., 0., 0., reml, 0., 100., 1e-8, 0.1, 0.5],
        )
        .unwrap();
        near(values(&j, "intercept")[0], 3., 1e-10);
        near(values(&j, "total_variance")[0], want, 1e-9);
        near(values(&j, "nugget")[0], 0., 1e-12);
        let expected = if reml == 0. {
            -0.5 * (3. * ((2. * std::f64::consts::PI).ln() + want.ln() + 1.))
        } else {
            -0.5 * (2. * ((2. * std::f64::consts::PI).ln() + want.ln() + 1.) + 3_f64.ln())
        };
        near(values(&j, "log_likelihood")[0], expected, 1e-9);
        assert!(j.contains("\"covariance_parameters_identified\":false"));
    }
}

// Dropping subject blocks would correlate observations that share coordinates.
#[test]
fn independent_subject_blocks_keep_collocated_observations_independent() {
    let d = vec![vec![0., 1., 1.], vec![0., 2., 2.], vec![0., 3., 6.]];
    let j = execute(
        "spatial_fit",
        &d,
        &[1., 0., 1., 0., 0., 100., 1e-8, 0.1, 2.],
    )
    .unwrap();
    near(values(&j, "total_variance")[0], 14. / 3., 1e-8);
    let saved = values(&j, "predict_parameters");
    let p = execute("spatial_fit_predict", &[vec![0., 4.]], &saved).unwrap();
    near(values(&p, "predictions")[0], 3., 1e-8);
    near(values(&p, "observation_variances")[0], 56. / 9., 1e-8);
}

// Rejecting an exact fitted trend is necessary: unknown variance has no positive
// finite Gaussian ML optimum, unlike existing known-covariance spatial_gls.
#[test]
fn estimated_covariance_rejects_zero_variance_aliases_and_invalid_shapes() {
    assert!(execute(
        "spatial_fit",
        &[vec![0., 0., 2.], vec![1., 0., 2.], vec![2., 0., 2.]],
        &[1., 0., 1., 1., 1.]
    )
    .is_err());
    assert!(execute(
        "spatial_fit",
        &[
            vec![0., 1., 0., 1.],
            vec![1., 1., 0., 2.],
            vec![2., 1., 0., 4.]
        ],
        &[1., 1., 1., 1., 1.]
    )
    .is_err());
    assert!(execute("spatial_fit_predict", &[vec![0., 0.]], &[]).is_err());
}

// Independent SciPy dense-covariance full three-parameter BFGS ML/REML fits.
// Fixing the range, suppressing the nugget, or not profiling the trend fails.
fn correlated_fixture() -> Vec<Vec<f64>> {
    let y = [
        3.680779827821175,
        3.49955299615519,
        4.245665084208936,
        2.985080583137176,
        3.6828878639053784,
        3.7513926741099475,
        4.823887252424179,
        2.9544118969350306,
        2.816875001593167,
        2.6808727577622538,
        1.8054882526399254,
        2.541261702417471,
        3.3318529172854947,
        4.3046048697620325,
        5.1639901527723495,
        2.80379486301126,
        3.525871037267851,
        2.163388276159442,
        4.265647365135814,
        3.596893023869619,
        2.6384275848630123,
        2.6432640650732577,
        1.9374016337868554,
        2.99800250508825,
    ];
    y.iter()
        .enumerate()
        .map(|(i, &v)| {
            let x = i as f64 / 4.;
            vec![x, x.sin(), if i < 12 { 0. } else { 1. }, v]
        })
        .collect()
}

#[test]
fn spatial_profiles_range_nugget_and_trend_against_independent_dense_likelihood() {
    let d = correlated_fixture();
    for (reml, sill, nugget, range, ll, intercept, slope) in [
        (
            0.,
            0.40382273,
            0.30558989,
            0.27121998,
            -29.3383933455387,
            3.27188052,
            0.12658942,
        ),
        (
            1.,
            0.54267705,
            0.45583567,
            1.2764756,
            -29.8049167191053,
            3.23068331,
            0.22203076,
        ),
    ] {
        let j = execute(
            "spatial_fit",
            &d,
            &[1., 1., 1., reml, 1., 1500., 1e-9, 0.01, 30.],
        )
        .unwrap();
        for (k, v) in [
            ("partial_sill", sill),
            ("nugget", nugget),
            ("range_parameter", range),
            ("log_likelihood", ll),
            ("intercept", intercept),
            ("coefficients", slope),
        ] {
            near(values(&j, k)[0], v, 3e-4);
        }
        assert!(j.contains("\"covariance_parameters_identified\":true"));
        assert_eq!(values(&j, "covariance_parameter_standard_errors").len(), 3);
        assert!(values(&j, "projected_score_norm")[0] < 1e-4);
        let saved = values(&j, "predict_parameters");
        let query = vec![vec![6., 6_f64.sin(), 1.], vec![0., 0., 0.]];
        let p = execute("spatial_fit_predict", &query, &saved).unwrap();
        assert_eq!(values(&p, "predictions").len(), 2);
        assert!(values(&p, "observation_variances").iter().all(|&v| v >= 0.));
        assert!(saved.len() <= 16384);
    }
}

// A positive Hessian at an arbitrary iteration-limited point does not justify
// coefficient-of-covariance uncertainty or a nominal confidence interval.
#[test]
fn nonconverged_spatial_search_withholds_covariance_parameter_intervals() {
    let j = execute(
        "spatial_fit",
        &correlated_fixture(),
        &[1., 1., 1., 0., 1., 1., 1e-9, 0.01, 30.],
    )
    .unwrap();
    assert!(j.contains("\"optimization_converged\":false"));
    assert!(j.contains("\"covariance_parameter_intervals\":null"));
}
