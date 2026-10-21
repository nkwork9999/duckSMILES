use ducksmiles_analytics::statistics::run;

fn numeric(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = json[start..].trim_start();
    let end = tail.find([',', '}', ']']).unwrap_or(tail.len());
    tail[..end].parse().expect("JSON number")
}

fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let end = json[start..].find(']').unwrap() + start;
    json[start..end]
        .split(',')
        .map(|n| n.trim().parse().unwrap())
        .collect()
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() < 1e-9, "{actual} != {expected}");
}

#[test]
fn summaries_use_sample_variance_and_preserve_constant_columns() {
    let json = run(
        "summarize",
        &[vec![1., 4.], vec![2., 4.], vec![3., 4.]],
        &[],
    )
    .unwrap();
    assert_eq!(array(&json, "means"), vec![2., 4.]);
    assert_eq!(array(&json, "sample_variances"), vec![1., 0.]);
    assert_eq!(array(&json, "minima"), vec![1., 4.]);
    assert_eq!(array(&json, "maxima"), vec![3., 4.]);
    assert_eq!(array(&json, "standard_deviations"), vec![1., 0.]);
    let one = run("summarize", &[vec![0.1]], &[]).unwrap();
    assert!(one.contains("\"sample_variances\":[null]"));
}

#[test]
fn welch_uses_unpooled_variance_and_a_two_sided_student_probability() {
    // Group 0 mean=1,var=2,n=2; group 1 mean=0,var=0,n=2 -> t=1,df=1,p=1/2.
    let json = run(
        "welch_t_test",
        &[vec![0., 0.], vec![2., 0.], vec![0., 1.], vec![0., 1.]],
        &[],
    )
    .unwrap();
    close(numeric(&json, "t_statistic"), 1.);
    close(numeric(&json, "degrees_of_freedom"), 1.);
    close(numeric(&json, "p_value"), 0.5);
}

#[test]
fn one_way_anova_reports_mean_square_ratio_and_f_tail_probability() {
    // SSB=4,SSW=4,df=(1,2),F=2. F(1,2) CDF is sqrt(F/(F+2)).
    let json = run(
        "one_way_anova",
        &[vec![0., 0.], vec![2., 0.], vec![2., 1.], vec![4., 1.]],
        &[],
    )
    .unwrap();
    close(numeric(&json, "f_statistic"), 2.);
    close(numeric(&json, "df_between"), 1.);
    close(numeric(&json, "df_within"), 2.);
    close(numeric(&json, "p_value"), 1. - 0.5_f64.sqrt());
}

#[test]
fn holm_and_bh_adjustments_retain_original_row_order() {
    let data = vec![vec![0.01], vec![0.04], vec![0.03], vec![0.002]];
    let holm = array(
        &run("adjust_pvalues", &data, &[0.]).unwrap(),
        "adjusted_p_values",
    );
    let bh = array(
        &run("adjust_pvalues", &data, &[1.]).unwrap(),
        "adjusted_p_values",
    );
    for (actual, expected) in holm.iter().zip([0.03, 0.06, 0.06, 0.008]) {
        close(*actual, expected);
    }
    for (actual, expected) in bh.iter().zip([0.02, 0.04, 0.04, 0.008]) {
        close(*actual, expected);
    }
}

#[test]
fn capability_uses_sample_sd_and_explicit_individual_reference_limits() {
    let json = run(
        "process_capability",
        &[vec![1.], vec![2.], vec![3.]],
        &[-1., 5.],
    )
    .unwrap();
    close(numeric(&json, "pp"), 1.);
    close(numeric(&json, "ppk"), 1.);
    assert!(json.contains("\"cp\":null"));
    assert!(json.contains("\"cpk\":null"));
    close(numeric(&json, "individual_lower_reference_limit"), -1.);
    close(numeric(&json, "individual_upper_reference_limit"), 5.);
    assert!(json.contains("\"sigma_estimator\":\"overall_sample_standard_deviation\""));
}

#[test]
fn cp_and_cpk_require_a_supplied_within_process_sigma() {
    let json = run(
        "process_capability",
        &[vec![1.], vec![2.], vec![3.]],
        &[-1., 5., 0.5],
    )
    .unwrap();
    close(numeric(&json, "cp"), 2.);
    close(numeric(&json, "cpk"), 2.);
    close(numeric(&json, "pp"), 1.);
    close(numeric(&json, "ppk"), 1.);
    for sigma in [0., -1.] {
        assert!(run(
            "process_capability",
            &[vec![1.], vec![2.], vec![3.]],
            &[-1., 5., sigma]
        )
        .is_err());
    }
}

#[test]
fn welch_probability_retains_precision_for_nearly_equal_means() {
    // Variance/2=1 in group0 and zero in group1: t=1e-8,df=1.
    let json = run(
        "welch_t_test",
        &[
            vec![-1. + 1e-8, 0.],
            vec![1. + 1e-8, 0.],
            vec![0., 1.],
            vec![0., 1.],
        ],
        &[],
    )
    .unwrap();
    let p = numeric(&json, "p_value");
    assert!((p - 0.9999999936338023).abs() < 1e-12, "{p}");
}

#[test]
fn kaplan_meier_keeps_tied_censors_in_the_event_risk_set() {
    let data = vec![vec![1., 1.], vec![2., 1.], vec![2., 0.], vec![3., 1.]];
    let json = run("kaplan_meier", &data, &[]).unwrap();
    assert_eq!(array(&json, "times"), vec![1., 2., 3.]);
    assert_eq!(array(&json, "at_risk"), vec![4., 3., 1.]);
    assert_eq!(array(&json, "survival"), vec![0.75, 0.5, 0.]);
    assert_eq!(
        array(&json, "greenwood_variance"),
        vec![0.046875, 0.0625, 0.]
    );
    close(numeric(&json, "median_survival"), 2.);
}

#[test]
fn all_censored_survival_stays_one_with_no_identified_median() {
    let json = run("kaplan_meier", &[vec![1., 0.], vec![2., 0.]], &[]).unwrap();
    assert_eq!(array(&json, "survival"), vec![1., 1.]);
    assert!(json.contains("\"median_survival\":null"));
}

#[test]
fn balanced_random_intercept_separates_within_and_between_variance() {
    // MS_between=4,MS_within=2,m=2 -> between=1,within=2,ICC=1/3.
    let data = vec![vec![0., 0.], vec![2., 0.], vec![2., 1.], vec![4., 1.]];
    let json = run("random_intercept", &data, &[]).unwrap();
    close(numeric(&json, "grand_mean"), 2.);
    close(numeric(&json, "between_group_variance"), 1.);
    close(numeric(&json, "within_group_variance"), 2.);
    close(numeric(&json, "icc"), 1. / 3.);
    assert!(json.contains("balanced_random_intercept_anova"));
}

#[test]
fn statistics_reject_undefined_models_invalid_labels_and_nonfinite_inputs() {
    for op in [
        "summarize",
        "welch_t_test",
        "one_way_anova",
        "adjust_pvalues",
        "process_capability",
        "kaplan_meier",
        "random_intercept",
    ] {
        assert!(run(op, &[], &[]).is_err());
        assert!(run(op, &[vec![f64::NAN]], &[]).is_err());
        assert!(run(op, &[vec![1., 2.], vec![3.]], &[]).is_err());
    }
    assert!(run("welch_t_test", &[vec![1., 0.], vec![2., 1.]], &[]).is_err());
    assert!(run(
        "welch_t_test",
        &[vec![1., 0.], vec![1., 0.], vec![2., 1.], vec![2., 1.]],
        &[]
    )
    .is_err());
    assert!(run("one_way_anova", &[vec![1., 0.], vec![2., 0.]], &[]).is_err());
    assert!(run("adjust_pvalues", &[vec![1.1]], &[0.]).is_err());
    assert!(run("process_capability", &[vec![1.], vec![1.]], &[0., 2.]).is_err());
    assert!(run("kaplan_meier", &[vec![-1., 1.]], &[]).is_err());
    assert!(run("kaplan_meier", &[vec![1., 0.5]], &[]).is_err());
    assert!(run(
        "random_intercept",
        &[vec![1., 0.], vec![2., 0.], vec![3., 1.]],
        &[]
    )
    .is_err());
}

#[test]
fn statistical_variance_and_group_differences_preserve_small_spreads_at_large_origins() {
    let origin = 1e16;
    let data = vec![
        vec![origin],
        vec![origin + 2.],
        vec![origin + 4.],
        vec![origin + 6.],
    ];
    let json = run("summarize", &data, &[]).unwrap();
    close(array(&json, "sample_variances")[0], 20. / 3.);
    let grouped = vec![
        vec![origin, 0.],
        vec![origin + 2., 0.],
        vec![origin + 2., 1.],
        vec![origin + 4., 1.],
    ];
    let welch = run("welch_t_test", &grouped, &[]).unwrap();
    close(numeric(&welch, "mean_difference"), -2.);
    close(numeric(&welch, "p_value"), 1. - 0.5_f64.sqrt());
    let anova = run("one_way_anova", &grouped, &[]).unwrap();
    close(numeric(&anova, "f_statistic"), 2.);
    let random = run("random_intercept", &grouped, &[]).unwrap();
    close(numeric(&random, "within_group_variance"), 2.);
    close(numeric(&random, "between_group_variance"), 1.);
    let capability = run("process_capability", &data, &[origin - 4., origin + 10.]).unwrap();
    close(
        numeric(&capability, "ppk"),
        14. / (6. * (20_f64 / 3.).sqrt()),
    );
}

#[test]
fn unequal_group_tail_probabilities_match_independent_reference_fixtures() {
    // Values independently evaluated with SciPy 1.17.1 ttest_ind(equal_var=False)
    // and f_oneway, exercising noninteger Welch df and asymmetric beta shapes.
    let mut welch_rows = [1., 2., 4., 8., 16.].map(|value| vec![value, 0.]).to_vec();
    welch_rows.extend([2., 5., 11.].map(|value| vec![value, 1.]));
    let json = run("welch_t_test", &welch_rows, &[]).unwrap();
    close(numeric(&json, "degrees_of_freedom"), 5.438766354360119);
    close(numeric(&json, "p_value"), 0.9599045974132303);
    let mut anova_rows = [1., 2., 4., 8.].map(|value| vec![value, 0.]).to_vec();
    anova_rows.extend([2., 5., 11.].map(|value| vec![value, 1.]));
    anova_rows.extend([3., 6., 7., 12., 14.].map(|value| vec![value, 2.]));
    let json = run("one_way_anova", &anova_rows, &[]).unwrap();
    close(numeric(&json, "f_statistic"), 1.4304047384007896);
    close(numeric(&json, "p_value"), 0.28878644760631034);
}
