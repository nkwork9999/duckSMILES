use ducksmiles_analytics::inference::run;

// These literals were independently evaluated with SciPy 1.17.1. Tests target
// wrong distributions, wrong tails, ignored pairing/ties, and CI inversion.
fn numeric(json: &str, key: &str) -> f64 {
    let start = json.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = json[start..].trim_start();
    tail[..tail.find([',', '}', ']']).unwrap_or(tail.len())]
        .parse()
        .expect("finite JSON number")
}
fn array(json: &str, key: &str) -> Vec<f64> {
    let start = json.find(&format!("\"{key}\":[")).expect(key) + key.len() + 4;
    let end = json[start..].find(']').unwrap() + start;
    json[start..end]
        .split(',')
        .map(|s| s.parse().unwrap())
        .collect()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * expected.abs().max(1.),
        "{actual} != {expected}"
    );
}
fn a() -> Vec<Vec<f64>> {
    [1., 2., 4., 5., 8.].map(|v| vec![v]).to_vec()
}
fn paired() -> Vec<Vec<f64>> {
    [1., 2., 4., 5., 8.]
        .into_iter()
        .zip([2., 3., 3., 6., 9.])
        .map(|(x, y)| vec![x, y])
        .collect()
}
fn grouped() -> Vec<Vec<f64>> {
    let mut rows = a().into_iter().map(|r| vec![r[0], 0.]).collect::<Vec<_>>();
    rows.extend([2., 3., 3., 6., 9.].map(|v| vec![v, 1.]));
    rows
}
fn three_groups() -> Vec<Vec<f64>> {
    let mut rows = grouped();
    rows.extend([1., 1., 2., 2., 3.].map(|v| vec![v, 2.]));
    rows
}

#[test]
fn one_sample_t_inverts_student_t_for_mean_ci() {
    let json = run("one_sample_t", &a(), &[3.]).unwrap();
    close(numeric(&json, "t_statistic"), 0.8164965809277261);
    close(numeric(&json, "p_value"), 0.460050752843961);
    close(numeric(&json, "degrees_of_freedom"), 4.);
    let ci = array(&json, "confidence_interval");
    close(ci[0], 0.5995630967087218);
    close(ci[1], 7.400436903291278);
}
#[test]
fn paired_t_uses_within_pair_differences() {
    let json = run("paired_t", &paired(), &[]).unwrap();
    close(numeric(&json, "t_statistic"), -1.5);
    close(numeric(&json, "p_value"), 0.208);
    close(numeric(&json, "mean_difference"), -0.6);
    let ci = array(&json, "confidence_interval");
    close(ci[0], -1.7105780420791175);
    close(ci[1], 0.5105780420791176);
}
#[test]
fn pooled_t_uses_pooled_variance_and_combined_df() {
    let json = run("pooled_t", &grouped(), &[]).unwrap();
    close(numeric(&json, "t_statistic"), -0.33752637027780696);
    close(numeric(&json, "p_value"), 0.7444056354601131);
    close(numeric(&json, "degrees_of_freedom"), 8.);
    let ci = array(&json, "confidence_interval");
    close(ci[0], -4.699242616165666);
    close(ci[1], 3.4992426161656667);
}
#[test]
fn variance_test_uses_chi_square_equal_tail_inference() {
    let json = run("variance_test", &a(), &[2.]).unwrap();
    close(numeric(&json, "chi_square_statistic"), 7.5);
    close(numeric(&json, "p_value"), 0.22341858563208655);
    let ci = array(&json, "variance_confidence_interval");
    close(ci[0], 2.6922038880654733);
    close(ci[1], 61.929914866070064);
}
#[test]
fn variance_ratio_divides_by_null_ratio_and_inverts_f_quantiles() {
    let json = run("variance_ratio_test", &grouped(), &[]).unwrap();
    close(numeric(&json, "f_statistic"), 0.9036144578313252);
    close(numeric(&json, "p_value"), 0.9241155364991409);
    let ci = array(&json, "variance_ratio_confidence_interval");
    close(ci[0], 0.09408211215716357);
    close(ci[1], 8.67879206450861);
    let scaled = run("variance_ratio_test", &grouped(), &[2.]).unwrap();
    close(numeric(&scaled, "f_statistic"), 0.4518072289156626);
    assert_eq!(array(&scaled, "variance_ratio_confidence_interval"), ci);
}
#[test]
fn binomial_test_is_probability_ordered_exact_with_wilson_interval() {
    let rows = [1., 1., 1., 0., 0., 0., 0., 0., 0., 0.].map(|v| vec![v]);
    let json = run("proportion_test", &rows, &[0.2]).unwrap();
    close(numeric(&json, "p_value"), 0.42957465600000017);
    let ci = array(&json, "confidence_interval");
    close(ci[0], 0.10779126740630102);
    close(ci[1], 0.6032218525388546);
    close(
        numeric(
            &run("proportion_test", &vec![vec![0.]; 10], &[0.2]).unwrap(),
            "p_value",
        ),
        0.22824806400000008,
    );
    close(
        numeric(
            &run("proportion_test", &vec![vec![1.]; 10], &[0.2]).unwrap(),
            "p_value",
        ),
        1.024e-7,
    );
    close(
        numeric(
            &run("proportion_test", &vec![vec![0.]; 10], &[0.]).unwrap(),
            "p_value",
        ),
        1.,
    );
    close(
        numeric(
            &run("proportion_test", &[vec![1.]], &[0.]).unwrap(),
            "p_value",
        ),
        0.,
    );
}
#[test]
fn contingency_uses_uncorrected_pearson_and_cramers_v() {
    let json = run(
        "contingency_test",
        &[vec![12., 5., 3.], vec![4., 9., 7.]],
        &[],
    )
    .unwrap();
    close(numeric(&json, "chi_square_statistic"), 6.742857142857142);
    close(numeric(&json, "p_value"), 0.034340544364887336);
    close(numeric(&json, "degrees_of_freedom"), 2.);
    close(
        numeric(&json, "cramers_v"),
        (6.742857142857142_f64 / 40.).sqrt(),
    );
    assert!(json.contains("\"expected_counts\":[[8,7,5],[8,7,5]]"));
}
#[test]
fn fisher_sums_tables_at_most_as_probable_as_observed() {
    let json = run("fisher_exact", &[vec![1., 9.], vec![11., 3.]], &[]).unwrap();
    close(numeric(&json, "p_value"), 0.0027594561852200836);
    close(numeric(&json, "odds_ratio"), 0.030303030303030304);
    let zero = run("fisher_exact", &[vec![0., 4.], vec![2., 3.]], &[]).unwrap();
    close(numeric(&zero, "odds_ratio"), 0.);
    let infinity = run("fisher_exact", &[vec![4., 0.], vec![2., 3.]], &[]).unwrap();
    assert!(infinity.contains("\"odds_ratio\":null"));
    assert!(infinity.contains("\"odds_ratio_infinite\":true"));
}
#[test]
fn mann_whitney_uses_average_ties_and_continuity_corrected_normal_tail() {
    let json = run("mann_whitney", &grouped(), &[]).unwrap();
    close(numeric(&json, "u_statistic"), 10.5);
    close(numeric(&json, "p_value"), 0.7525671000692723);
    assert!(json.contains("tie_corrected_normal_approximation"));
    assert!(json.contains("\"continuity_correction\":true"));
}
#[test]
fn wilcoxon_drops_zero_pairs_and_corrects_tied_absolute_ranks() {
    let json = run("wilcoxon_signed_rank", &paired(), &[]).unwrap();
    close(numeric(&json, "rank_statistic"), 3.);
    close(numeric(&json, "p_value"), 0.23303798227390538);
    let mut rows = paired();
    rows.push(vec![9., 9.]);
    let same = run("wilcoxon_signed_rank", &rows, &[]).unwrap();
    close(numeric(&same, "p_value"), numeric(&json, "p_value"));
    close(numeric(&same, "n_nonzero"), 5.);
}
#[test]
fn kruskal_wallis_applies_global_tie_correction() {
    let json = run("kruskal_wallis", &three_groups(), &[]).unwrap();
    close(numeric(&json, "h_statistic"), 4.220664206642074);
    close(numeric(&json, "p_value"), 0.12119770958717949);
}
#[test]
fn levene_defaults_to_brown_forsythe_median_deviations() {
    let json = run("levene_test", &three_groups(), &[]).unwrap();
    close(numeric(&json, "f_statistic"), 1.053763440860215);
    close(numeric(&json, "p_value"), 0.37877574177837225);
    assert!(json.contains("\"center\":\"median\""));
}
#[test]
fn equivalence_uses_two_one_sided_tests_and_one_minus_two_alpha_ci() {
    let json = run("equivalence_t", &a(), &[2., 6.]).unwrap();
    close(numeric(&json, "p_lower"), 0.08890390417811062);
    close(numeric(&json, "p_upper"), 0.08890390417811062);
    close(numeric(&json, "p_value"), 0.08890390417811062);
    assert!(json.contains("\"equivalent\":false"));
    let ci = array(&json, "confidence_interval");
    close(ci[0], 1.3890315818537937);
    close(ci[1], 6.6109684181462045);
    assert!(run("equivalence_t", &a(), &[0., 8.])
        .unwrap()
        .contains("\"equivalent\":true"));
}

#[test]
fn inference_rejects_malformed_nonfinite_undefined_and_over_limit_inputs() {
    for op in [
        "one_sample_t",
        "paired_t",
        "pooled_t",
        "variance_test",
        "variance_ratio_test",
        "proportion_test",
        "contingency_test",
        "fisher_exact",
        "mann_whitney",
        "wilcoxon_signed_rank",
        "kruskal_wallis",
        "levene_test",
        "equivalence_t",
    ] {
        for rows in [
            vec![],
            vec![vec![f64::NAN]],
            vec![vec![1., 2.], vec![3.]],
            vec![vec![]],
            vec![vec![0.]; 100_001],
        ] {
            assert!(run(op, &rows, &[]).is_err(), "{op}");
        }
        assert!(run(op, &a(), &[f64::INFINITY]).is_err(), "{op}");
    }
    for (op, rows, args) in [
        ("one_sample_t", a(), vec![]),
        ("one_sample_t", a(), vec![0., 1.]),
        ("one_sample_t", vec![vec![1.]; 3], vec![0.]),
        ("one_sample_t", vec![vec![1e308], vec![-1e308]], vec![0.]),
        ("paired_t", vec![vec![1., 1.]; 3], vec![]),
        ("pooled_t", vec![vec![1., 0.], vec![2., 1.]], vec![]),
        (
            "pooled_t",
            vec![vec![1., 0.], vec![2., 0.], vec![3., 1.5], vec![4., 1.]],
            vec![],
        ),
        ("variance_test", a(), vec![0.]),
        ("variance_ratio_test", grouped(), vec![-1.]),
        ("proportion_test", vec![vec![0.5]], vec![0.5]),
        ("contingency_test", vec![vec![1., 0.], vec![3., 0.]], vec![]),
        (
            "contingency_test",
            vec![vec![1.5, 2.], vec![3., 4.]],
            vec![],
        ),
        (
            "fisher_exact",
            vec![vec![100_001., 1.], vec![1., 1.]],
            vec![],
        ),
        ("mann_whitney", vec![vec![1., 0.], vec![1., 1.]], vec![]),
        ("wilcoxon_signed_rank", vec![vec![1., 1.]; 3], vec![]),
        ("kruskal_wallis", vec![vec![1., 0.], vec![1., 1.]], vec![]),
        ("kruskal_wallis", vec![vec![1., -1.], vec![2., 1.]], vec![]),
        ("levene_test", three_groups(), vec![2.]),
        ("equivalence_t", a(), vec![6., 2.]),
        ("equivalence_t", a(), vec![2., 6., 0.5]),
    ] {
        assert!(run(op, &rows, &args).is_err(), "{op}: {args:?}");
    }
    assert!(run("not_an_operation", &a(), &[]).is_err());
}

#[test]
fn deterministic_permutations_and_translations_preserve_inferential_results() {
    let mut seed = 0x5eed_u64;
    for _ in 0..32 {
        let mut rows = grouped();
        for i in (1..rows.len()).rev() {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            rows.swap(i, seed as usize % (i + 1));
        }
        let offset = ((seed >> 32) % 1000) as f64 - 500.;
        let translated = rows
            .iter()
            .map(|r| vec![r[0] + offset, r[1]])
            .collect::<Vec<_>>();
        for op in [
            "pooled_t",
            "variance_ratio_test",
            "mann_whitney",
            "kruskal_wallis",
            "levene_test",
        ] {
            close(
                numeric(&run(op, &rows, &[]).unwrap(), "p_value"),
                numeric(&run(op, &translated, &[]).unwrap(), "p_value"),
            );
        }
        let json = run("one_sample_t", &a(), &[3.]).unwrap();
        let translated_a = a().iter().map(|r| vec![r[0] + offset]).collect::<Vec<_>>();
        close(
            numeric(&json, "p_value"),
            numeric(
                &run("one_sample_t", &translated_a, &[3. + offset]).unwrap(),
                "p_value",
            ),
        );
    }
}

#[test]
fn student_and_chi_square_confidence_limits_retain_small_tail_probabilities() {
    let json = run("one_sample_t", &a(), &[3., 0.999999]).unwrap();
    let ci = array(&json, "confidence_interval");
    close((ci[1] - 4.) / 1.5_f64.sqrt(), 49.458636757953215);
    let json = run("variance_test", &a(), &[2., 0.999999]).unwrap();
    let ci = array(&json, "variance_confidence_interval");
    close(ci[0], 30. / 34.84443002506143);
    close(ci[1], 30. / 0.002000666972381571);
}

#[test]
fn exact_test_probabilities_sum_to_one_at_large_count_modes() {
    // Every possible table/outcome is no more probable than a modal one.
    // This catches accumulated log-factorial bias in unnormalized enumeration.
    let json = run(
        "fisher_exact",
        &[vec![25000., 24999.], vec![24999., 25000.]],
        &[],
    )
    .unwrap();
    assert!((numeric(&json, "p_value") - 1.).abs() < 1e-12);
    let mut rows = vec![vec![1.]; 50_000];
    rows.extend(vec![vec![0.]; 50_000]);
    let json = run("proportion_test", &rows, &[0.5]).unwrap();
    assert!((numeric(&json, "p_value") - 1.).abs() < 1e-12);
}
