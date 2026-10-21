#[path = "../../src/sem.rs"]
mod sem;
use sem::run;

fn numbers(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let t = &j[start..];
    let end = if t.starts_with('[') {
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
    t[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-5 * b.abs().max(1.), "{a} != {b}");
}
fn all(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        near(x, y);
    }
}

fn fixture() -> Vec<Vec<f64>> {
    vec![
        vec![
            1.0,
            0.7200000000000001,
            0.5599999999999999,
            0.23799999999999996,
            0.20999999999999996,
            0.242,
            10.0,
        ],
        vec![
            0.7200000000000001,
            1.0,
            0.63,
            0.26775,
            0.23625000000000002,
            0.20475000000000002,
            -2.0,
        ],
        vec![
            0.5599999999999999,
            0.63,
            1.0,
            0.20824999999999996,
            0.18374999999999997,
            0.15924999999999997,
            0.0,
        ],
        vec![0.238, 0.26775, 0.20825, 1.0, 0.6375, 0.5525, 5.0],
        vec![
            0.20999999999999996,
            0.23624999999999996,
            0.18374999999999997,
            0.6375,
            1.0,
            0.48750000000000004,
            1.0,
        ],
        vec![
            0.242,
            0.20475,
            0.15924999999999997,
            0.5525,
            0.48750000000000004,
            1.0,
            -3.0,
        ],
    ]
}

fn controls() -> Vec<f64> {
    vec![300., 2., 1., 500., 1e-8, 0., 0., 0., 1., 1., 1.]
}
fn one_factor(load: &[f64]) -> Vec<Vec<f64>> {
    let p = load.len();
    (0..p)
        .map(|i| {
            let mut r: Vec<_> = (0..p)
                .map(|j| if i == j { 1. } else { load[i] * load[j] })
                .collect();
            r.push(i as f64 + 1.);
            r
        })
        .collect()
}
#[test]
fn cfa_true_ml_recovers_prespecified_simple_structure_and_counts_free_means() {
    let load = [0.8, 0.7, 0.9, 0.6];
    let j = run(
        "cfa_model",
        &one_factor(&load),
        &[400., 1., 0., 500., 1e-8, 0., 0., 0., 0.],
    )
    .unwrap();
    all(&numbers(&j, "loadings"), &load);
    all(
        &numbers(&j, "uniquenesses"),
        &load.iter().map(|v| 1. - v * v).collect::<Vec<_>>(),
    );
    near(numbers(&j, "n_free_parameters")[0], 12.);
    near(numbers(&j, "degrees_of_freedom")[0], 2.);
    near(numbers(&j, "chi_square")[0], 0.);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn correlated_cfa_loadings_and_observed_information_match_independent_scipy() {
    let j = run("cfa_model", &fixture(), &controls()).unwrap();
    all(
        &numbers(&j, "loadings"),
        &[
            0.8021360052186588,
            0.0,
            0.8978014708739563,
            0.0,
            0.7004290167124376,
            0.0,
            0.0,
            0.8480539927098586,
            0.0,
            0.7503326427572304,
            0.0,
            0.6522928340594946,
        ],
    );
    all(
        &numbers(&j, "uniquenesses"),
        &[
            0.356577829131852,
            0.1939525188965607,
            0.5093991925472479,
            0.2808044254488666,
            0.43700092521295053,
            0.5745140586346326,
        ],
    );
    all(
        &numbers(&j, "standard_errors"),
        &[
            0.05222877233123071,
            0.050621115392614185,
            0.05369741957663936,
            0.05559019244087833,
            0.056172406465734785,
            0.05700979101389086,
            0.04526874830157349,
            0.045791692621829395,
            0.0495134920230484,
            0.057225119073992414,
            0.05463890333676904,
            0.057146599544223874,
            0.061295083845967956,
        ],
    );
    near(numbers(&j, "chi_square")[0], 3.5494789987331465);
    near(numbers(&j, "srmr")[0], 0.010956160855163733);
    near(numbers(&j, "n_free_parameters")[0], 19.);
}
#[test]
fn saved_cfa_bartlett_scores_invert_loading_space_and_regression_reports_posterior_covariance() {
    let j = run("cfa_model", &fixture(), &controls()).unwrap();
    let l = numbers(&j, "loadings");
    let mean = numbers(&j, "means");
    let row: Vec<_> = (0..6)
        .map(|i| mean[i] + l[2 * i] - 0.5 * l[2 * i + 1])
        .collect();
    let saved = numbers(&j, "score_parameters");
    let mut args = vec![1.];
    args.extend(&saved);
    let p = run("cfa_scores", std::slice::from_ref(&row), &args).unwrap();
    all(&numbers(&p, "scores"), &[1., -0.5]);
    args[0] = 0.;
    let p = run("cfa_scores", &[row], &args).unwrap();
    assert_eq!(numbers(&p, "conditional_factor_covariance").len(), 4);
}
#[test]
fn sem_fit_indices_include_nonzero_mean_residuals_and_jmp_n_scaling() {
    let data = vec![
        vec![1., 0.3, 1.2, 0.2, 1., 1.1],
        vec![0.3, 1., 0.2, 0.8, -2., -2.2],
    ];
    let j = run("sem_fit_indices", &data, &[120., 2.]).unwrap();
    let delta = [-0.1, 0.2];
    let inverse = [[0.8 / 0.92, -0.2 / 0.92], [-0.2 / 0.92, 1.2 / 0.92]];
    let q = delta[0] * (inverse[0][0] * delta[0] + inverse[0][1] * delta[1])
        + delta[1] * (inverse[1][0] * delta[0] + inverse[1][1] * delta[1]);
    let tr = inverse[0][0] + inverse[1][1] + 0.6 * inverse[0][1];
    let chi = 120. * ((0.92_f64 / 0.91).ln() + tr + q - 2.);
    near(numbers(&j, "chi_square")[0], chi);
    near(
        numbers(&j, "srmr")[0],
        ((0.04 + 0.01 + 0.04 + 0.01 + 0.04) / 5_f64).sqrt(),
    );
}
#[test]
fn rmsea_confidence_interval_inverts_noncentral_chi_square_and_zero_df_is_undefined() {
    let c = 1.8;
    let mut data = vec![vec![0.; 8]; 3];
    for (i, r) in data.iter_mut().enumerate() {
        r[i] = 1.;
        r[3 + i] = c;
    }
    let j = run("sem_fit_indices", &data, &[200., 2.]).unwrap();
    near(
        numbers(&j, "chi_square")[0],
        600. * (1.8_f64.ln() + 1. / 1.8 - 1.),
    );
    all(
        &numbers(&j, "rmsea_confidence_interval"),
        &[0.19413133904199983, 0.28372685108686996],
    );
    let data = vec![
        vec![1., 0.2, 1., 0.2, 0., 0.],
        vec![0.2, 1., 0.2, 1., 0., 0.],
    ];
    let j = run("sem_fit_indices", &data, &[100., 5.]).unwrap();
    assert!(j.contains("\"rmsea\":null"));
    assert!(j.contains("\"chi_square_p_value\":null"));
}
#[test]
fn ram_implied_covariance_and_mean_follow_simultaneous_linear_equations() {
    let data = vec![
        vec![0., 0., 0., 1., 0., 0., 1.],
        vec![0.5, 0., 0., 0., 0.8, 0., 2.],
        vec![0.2, 0.7, 0., 0., 0., 0.6, 3.],
    ];
    let j = run("sem_ram", &data, &[]).unwrap();
    all(&numbers(&j, "means"), &[1., 2.5, 4.95]);
    all(
        &numbers(&j, "covariance"),
        &[1., 0.5, 0.55, 0.5, 1.05, 0.835, 0.55, 0.835, 1.2945],
    );
    near(numbers(&j, "indirect_effects")[6], 0.35);
    assert!(j.contains("\"dynamically_stable\":true"));
}
#[test]
fn ram_invertibility_is_separate_from_dynamic_stability() {
    let data = vec![vec![0., 2., 1., 0., 0.], vec![0.75, 0., 0., 1., 0.]];
    let j = run("sem_ram", &data, &[]).unwrap();
    assert!(j.contains("\"dynamically_stable\":false"));
    assert!(run(
        "sem_ram",
        &[vec![0., 1., 1., 0., 0.], vec![1., 0., 0., 1., 0.]],
        &[]
    )
    .is_err());
}
#[test]
fn cfa_negative_loadings_have_fixed_positive_anchor_and_rescaling_preserves_fit() {
    let load = [0.8, -0.7, 0.9, -0.6];
    let data = one_factor(&load);
    let args = [400., 1., 0., 500., 1e-8, 0., 0., 0., 0.];
    let a = run("cfa_model", &data, &args).unwrap();
    all(&numbers(&a, "loadings"), &load);
    let scales = [0.1, 2., 10., 3.];
    let mut d = data.clone();
    for i in 0..4 {
        for j in 0..4 {
            d[i][j] *= scales[i] * scales[j];
        }
        d[i][4] += 1e6;
    }
    let b = run("cfa_model", &d, &args).unwrap();
    all(
        &numbers(&b, "loadings"),
        &load
            .iter()
            .zip(scales)
            .map(|(l, s)| l * s)
            .collect::<Vec<_>>(),
    );
    near(numbers(&a, "chi_square")[0], numbers(&b, "chi_square")[0]);
}
#[test]
fn heywood_boundary_and_iteration_limit_are_not_reported_as_identified_interior_success() {
    let data = vec![
        vec![1., 0.8, 0.8, 0.],
        vec![0.8, 1., 0.4, 0.],
        vec![0.8, 0.4, 1., 0.],
    ];
    let j = run("cfa_model", &data, &[500., 1., 0., 500., 1e-8, 0., 0., 0.]).unwrap();
    assert!(j.contains("\"heywood_boundary\":true"), "{j}");
    assert!(j.contains("\"converged\":false"), "{j}");
    assert!(j.contains("\"standard_errors\":null"), "{j}");
    let mut args = controls();
    args[3] = 1.;
    let j = run("cfa_model", &fixture(), &args).unwrap();
    assert!(j.contains("\"converged\":false"));
}
#[test]
fn invalid_summary_covariance_or_unidentified_assignments_are_rejected() {
    assert!(run(
        "cfa_model",
        &fixture(),
        &[300., 3., 1., 500., 1e-8, 0., 0., 1., 1., 2., 2.]
    )
    .is_err());
    assert!(run(
        "cfa_model",
        &[vec![1., 2., 0.], vec![2., 1., 0.]],
        &[20., 1., 0., 500., 1e-8, 0., 0.]
    )
    .is_err());
    assert!(run("cfa_scores", &[vec![1., 2.]], &[]).is_err());
    assert!(run("sem_fit_indices", &[vec![1., 1., 0., 0.]], &[100., 10.]).is_err());
}
#[test]
fn nonconverged_or_boundary_cfa_does_not_report_calibrated_fit_p_or_rmsea_ci() {
    let mut args = controls();
    args[3] = 1.;
    let j = run("cfa_model", &fixture(), &args).unwrap();
    assert!(j.contains("\"chi_square_p_value\":null"), "{j}");
    assert!(j.contains("\"rmsea_confidence_interval\":null"), "{j}");
}
#[test]
fn one_manifest_variable_has_zero_independence_df_and_undefined_incremental_indices() {
    for i in 0..20 {
        let s = 10_f64.powf(i as f64 / 7.);
        let j = run("sem_fit_indices", &[vec![s, 0.9, 0.3, 0.4]], &[100., 1.]).unwrap();
        assert!(j.contains("\"cfi\":null"), "{j}");
        assert!(j.contains("\"nfi\":null"), "{j}");
    }
}

#[test]
fn exactly_identified_three_indicator_cfa_reaches_analytic_ml_without_gradient_cancellation() {
    let data = vec![
        vec![
            0.13612195085746226,
            0.10126348607558294,
            0.6209281619407684,
            -5231.405047196925,
        ],
        vec![
            0.10126348607558296,
            0.21886973850186014,
            0.7319855362324363,
            -9973.877760766825,
        ],
        vec![
            0.6209281619407684,
            0.7319855362324363,
            6.1649317037050055,
            -2940.4000207187373,
        ],
    ];
    let j = run("cfa_model", &data, &[600., 1., 1., 800., 1e-9, 0., 0., 0.]).unwrap();
    let l0 = (data[0][1] * data[0][2] / data[1][2]).sqrt();
    all(
        &numbers(&j, "loadings"),
        &[l0, data[0][1] / l0, data[0][2] / l0],
    );
    near(numbers(&j, "chi_square")[0], 0.);
    assert!(j.contains("\"converged\":true"), "{j}");
}
#[test]
fn ram_unit_circle_stability_is_unresolved_while_static_solution_exists() {
    let j = run(
        "sem_ram",
        &[vec![0., -1., 1., 0., 0.], vec![1., 0., 0., 1., 0.]],
        &[],
    )
    .unwrap();
    assert!(j.contains("\"invertible\":true"));
    assert!(j.contains("\"dynamically_stable\":null"));
}
#[test]
fn rmsea_follows_displayed_formula_without_arbitrary_clipping_and_zero_root_is_explicit() {
    let j = run("sem_fit_indices", &[vec![1., 100., 0., 0.]], &[100., 0.]).unwrap();
    assert!(numbers(&j, "rmsea")[0] > 1.);
    let j = run("sem_fit_indices", &[vec![1., 1., 0., 0.]], &[100., 0.]).unwrap();
    all(&numbers(&j, "rmsea_confidence_interval"), &[0., 0.]);
}
