use ducksmiles_analytics::execute as run;
fn values(j: &str, k: &str) -> Vec<f64> {
    let t = j.split_once(&format!("\"{k}\":")).expect(k).1;
    let end = if t.starts_with('[') {
        t.find(']').unwrap() + 1
    } else {
        t.find([',', '}']).unwrap_or(t.len())
    };
    t[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(x: f64, y: f64) {
    assert!((x - y).abs() < 1e-10 * (1. + y.abs()), "{x} != {y}");
}
fn basket() -> Vec<Vec<f64>> {
    vec![
        vec![1., 1., 0.],
        vec![1., 1., 1.],
        vec![1., 0., 1.],
        vec![0., 1., 0.],
        vec![0., 0., 0.],
    ]
}
#[test]
fn rules_count_empty_transactions_and_allow_multivariate_consequents() {
    let j = run("association_rules", &basket(), &[0.2, 0., 0., 2., 3.]).unwrap();
    let rule = j
        .split("},{")
        .find(|r| r.contains("\"antecedent\":[0],\"consequent\":[2]"))
        .expect("A=>C");
    near(values(rule, "support")[0], 0.4);
    near(values(rule, "confidence")[0], 2. / 3.);
    near(values(rule, "lift")[0], 5. / 3.);
    near(values(rule, "leverage")[0], 0.16);
    near(values(rule, "conviction")[0], 1.8);
    assert!(j.contains("\"antecedent\":[0],\"consequent\":[1,2]"));
}
#[test]
fn rule_thresholds_filter_on_the_joint_support_confidence_and_lift() {
    let j = run("association_rules", &basket(), &[0.4, 0.6, 1.5, 1., 2.]).unwrap();
    near(values(&j, "rule_count")[0], 2.);
    assert!(j.contains("\"antecedent\":[0],\"consequent\":[2]"));
    assert!(j.contains("\"antecedent\":[2],\"consequent\":[0]"));
    let j = run(
        "association_rules",
        &[vec![1., 1.], vec![0., 0.]],
        &[0.5, 1., 0.],
    )
    .unwrap();
    assert!(j.contains("\"conviction\":null,\"conviction_infinite\":true"));
}
#[test]
fn basket_affinity_distinguishes_directional_confidence_and_symmetric_similarity() {
    let j = run("basket_affinity", &basket(), &[]).unwrap();
    let row = j
        .split("},{")
        .find(|r| r.contains("\"item0\":0,\"item1\":2"))
        .expect("pairAC");
    near(values(row, "confidence0_to1")[0], 2. / 3.);
    near(values(row, "confidence1_to0")[0], 1.);
    near(values(row, "jaccard")[0], 2. / 3.);
    near(values(row, "lift")[0], 5. / 3.);
}
#[test]
fn raw_and_standardized_item_reliability_use_their_own_covariance_scale() {
    let d = [
        vec![0., 0., 0.],
        vec![1., 2., 1.],
        vec![2., 4., 0.],
        vec![3., 6., 1.],
    ];
    let j = run("item_reliability", &d, &[]).unwrap();
    near(values(&j, "alpha")[0], 0.75);
    let deleted = values(&j, "alpha_if_deleted");
    for (&a, b) in deleted.iter().zip([0.32, 0.5, 8. / 9.]) {
        near(a, b);
    }
    for (&a, b) in values(&j, "corrected_item_total_correlations").iter().zip([
        0.9838699100999075,
        0.948683298050514,
        0.4472135954999579,
    ]) {
        near(a, b);
    }
    let j = run("item_reliability", &d, &[1.]).unwrap();
    near(values(&j, "alpha")[0], 0.8371488403210662);
    let negative = run(
        "item_reliability",
        &[vec![0., 1.], vec![1., 0.], vec![2., 0.]],
        &[],
    )
    .unwrap();
    assert!(values(&negative, "alpha")[0] < 0.);
}
#[test]
fn market_contracts_reject_nonbinary_baskets_and_report_undefined_reliability() {
    assert!(run("association_rules", &[vec![0., 2.]], &[]).is_err());
    assert!(run("association_rules", &basket(), &[0.]).is_err());
    assert!(run("basket_affinity", &basket(), &[1.]).is_err());
    let j = run("item_reliability", &[vec![0., 1.], vec![1., 0.]], &[]).unwrap();
    assert!(j.contains("\"alpha\":null"));
    assert!(run("item_reliability", &[vec![0., 1.], vec![0., 2.]], &[1.]).is_err());
}
#[test]
fn standardized_alpha_and_corrected_correlations_survive_different_item_units() {
    let d = [vec![0., 0.], vec![1e-300, 0.5e150], vec![2e-300, 1e150]];
    let j = run("item_reliability", &d, &[1.]).unwrap();
    near(values(&j, "alpha")[0], 1.);
    let d = [vec![0., 0.], vec![0.5e150, 0.5e100], vec![1e150, 1e100]];
    let j = run("item_reliability", &d, &[]).unwrap();
    let correlations = values(&j, "corrected_item_total_correlations");
    assert_eq!(correlations.len(), 2);
    for v in correlations {
        near(v, 1.);
    }
}
#[test]
fn deleting_a_dominant_item_preserves_the_variance_of_remaining_small_items() {
    let d = [
        vec![0., 0., 0.],
        vec![0.5e150, 0.5e100, 0.5e100],
        vec![1e150, 1e100, 1e100],
    ];
    let j = run("item_reliability", &d, &[]).unwrap();
    near(values(&j, "alpha_if_deleted")[0], 1.);
}
