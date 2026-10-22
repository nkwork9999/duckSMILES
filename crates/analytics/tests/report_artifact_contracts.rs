use ducksmiles_analytics::execute;

fn rows(values: &[f64]) -> Vec<Vec<f64>> {
    values
        .iter()
        .enumerate()
        .map(|(i, &v)| vec![v, 0., 1., 0., 0., 1., i as f64])
        .collect()
}
fn nums(json: &str, key: &str) -> Vec<f64> {
    let s = json.split_once(&format!("\"{key}\":")).expect(key).1;
    let s = s
        .strip_prefix('[')
        .expect("numeric array")
        .split(']')
        .next()
        .unwrap();
    if s.is_empty() {
        return vec![];
    }
    s.split(',').map(|v| v.parse().unwrap()).collect()
}
fn hist(data: &[Vec<f64>], scale: f64, orientation: f64, edges: &[f64]) -> String {
    let mut a = vec![1., 900., 600., scale, orientation, edges.len() as f64];
    a.extend(edges);
    execute("report_histogram_svg", data, &a).unwrap()
}
fn bx(data: &[Vec<f64>], kind: f64, orientation: f64, visible: f64, labels: f64) -> String {
    execute(
        "report_box_plot_svg",
        data,
        &[1., 900., 650., kind, orientation, visible, labels],
    )
    .unwrap()
}
#[test]
fn explicit_histogram_edges_include_final_edge_and_preserve_density_units() {
    let j = hist(
        &rows(&[-2., -1., -0.5, 0., 0.5, 1., 2.]),
        0.,
        0.,
        &[-2., -1., 0., 1., 2.],
    );
    assert_eq!(nums(&j, "counts"), vec![1., 2., 2., 2.]);
    let d = hist(&rows(&[0., 0.5, 1., 2., 3.]), 2., 1., &[0., 1., 3.]);
    assert_eq!(nums(&d, "counts"), vec![2., 3.]);
    assert_eq!(nums(&d, "densities"), vec![0.4, 0.3]);
    assert!(j.contains("\"schema\":\"ducksmiles-cpu-plot-v1\""));
    assert!(j.contains("http://www.w3.org/2000/svg"));
}
#[test]
fn box_quartiles_use_type6_and_observed_tukey_whiskers() {
    let j = bx(&rows(&[1., 2., 3., 4.]), 0., 0., 1., 1.);
    assert_eq!(nums(&j, "five_number"), vec![1., 1.25, 2.5, 3.75, 4.]);
    assert_eq!(nums(&j, "fences"), vec![-2.5, 7.5]);
    let j = bx(&rows(&[0., 1., 2., 3., 4., 5., 100.]), 0., 1., 1., 0.);
    assert_eq!(nums(&j, "whiskers"), vec![0., 5.]);
    assert!(j.contains("\"row_id\":6,\"value\":100,\"frequency\":1,\"hidden\":false"));
    let j = bx(&rows(&[-5., 1., 2., 3., 4., 5., 11.]), 0., 0., 1., 0.);
    assert_eq!(nums(&j, "whiskers"), vec![-5., 11.]);
    assert!(j.contains("\"outliers\":[]"));
}
#[test]
fn frequency_flags_and_empty_groups_keep_source_rows_distinct_from_n() {
    let mut d = rows(&[1., 2., 10., 999., 999., 999.]);
    d[1][2] = 2.;
    d[2][4] = 1.;
    d[3][3] = 1.;
    d[4][5] = 0.;
    d[5][2] = 0.;
    d[5][1] = 3.;
    let j = bx(&d, 0., 0., 1., 0.);
    assert_eq!(nums(&j, "five_number"), vec![1., 1.25, 2., 8., 10.]);
    assert!(j.contains("\"source_rows\":6,\"selected_rows\":3,\"N\":4"));
    assert!(j.contains("\"group_id\":3,\"status\":\"empty\""));
    assert!(j.contains("\"five_number\":null"));
    let h = hist(&d, 0., 0., &[0., 2., 10.]);
    assert_eq!(nums(&h, "counts"), vec![1., 3.]);
}
#[test]
fn hidden_or_disabled_outliers_keep_statistics_and_provenance_but_hide_marks() {
    let mut d = rows(&[0., 0., 0., 0., 0., 0., 10.]);
    d[6][4] = 1.;
    d[6][2] = 3.;
    // Twelve zero replications keep Q3 at zero with a 3-count outlier.
    d[0][2] = 7.;
    let j = bx(&d, 0., 0., 1., 0.);
    assert!(j.contains("\"row_id\":6,\"value\":10,\"frequency\":3,\"hidden\":true"));
    assert!(!j.contains("data-role=\\\"outlier\\\""));
    d[6][4] = 0.;
    let j = bx(&d, 0., 0., 0., 0.);
    assert!(!j.contains("data-role=\\\"outlier\\\""));
    assert_eq!(nums(&j, "whiskers"), vec![0., 0.]);
}
#[test]
fn quantile_box_has_full_range_and_four_additional_type6_ticks() {
    let j = bx(
        &rows(&(0..=100).map(|x| x as f64).collect::<Vec<_>>()),
        1.,
        0.,
        1.,
        0.,
    );
    assert_eq!(nums(&j, "whiskers"), vec![0., 100.]);
    let q = nums(&j, "quantile_ticks");
    for (a, b) in q.iter().zip([1.55, 9.2, 90.8, 98.45]) {
        assert!((a - b).abs() < 1e-12);
    }
}
#[test]
fn singleton_and_permuted_rows_have_deterministic_finite_artifacts() {
    let d = rows(&[42.]);
    let j = bx(&d, 0., 0., 1., 1.);
    assert_eq!(nums(&j, "five_number"), vec![42.; 5]);
    assert_eq!(j, bx(&d, 0., 0., 1., 1.));
    let d = rows(&[3., 1., 2., 100., 4., 5., 6.]);
    let mut rev = d.clone();
    rev.reverse();
    assert_eq!(bx(&d, 0., 1., 1., 0.), bx(&rev, 0., 1., 1., 0.));
}
#[test]
fn reports_reject_invalid_states_unrepresentable_ranges_and_dense_layouts() {
    let a = [1., 900., 600., 0., 0., 1., 0.];
    for (column, value) in [(1, 0.5), (2, 0.5), (3, 2.), (4, -1.), (5, 0.5), (6, -1.)] {
        let mut d = rows(&[1., 2.]);
        d[0][column] = value;
        assert!(execute("report_box_plot_svg", &d, &a).is_err());
    }
    let mut d = rows(&[1., 2.]);
    d[1][6] = 0.;
    assert!(execute("report_box_plot_svg", &d, &a).is_err());
    assert!(execute("report_box_plot_svg", &rows(&[-f64::MAX, f64::MAX]), &a).is_err());
    assert!(execute(
        "report_histogram_svg",
        &rows(&[0.]),
        &[1., 320., 240., 0., 0., 3., 0., 1e-300, 1.]
    )
    .is_err());
    let edges: Vec<f64> = (0..=512).map(|v| v as f64).collect();
    let mut args = vec![1., 320., 240., 0., 0., 513.];
    args.extend(edges);
    assert!(execute("report_histogram_svg", &rows(&[1.]), &args).is_err());
    assert!(execute(
        "report_histogram_svg",
        &rows(&[3.]),
        &[1., 900., 600., 0., 0., 3., 0., 1., 2.]
    )
    .is_err());
    let mut d = rows(&[1.]);
    d[0][2] = 1_000_001.;
    assert!(execute("report_box_plot_svg", &d, &a).is_err());
}

// A finite numerical domain still must have distinguishable displayed ticks.
#[test]
fn narrow_value_domains_keep_distinguishable_axis_labels() {
    let d = rows(&[1., 1. + 1e-8, 1. + 2e-8]);
    let a = [1., 900., 600., 0., 0., 3., 1., 1. + 1e-8, 1. + 2e-8];
    let actual = execute("report_histogram_svg", &d, &a).unwrap();
    assert!(
        actual.contains(">1.00000002</text>"),
        "distinct endpoint must have readable precise numeric label"
    );
    let d = rows(&[1e12, 1e12 + 0.01, 1e12 + 0.02]);
    let actual = execute("report_box_plot_svg", &d, &[1., 900., 650., 0., 0., 1., 0.]).unwrap();
    assert!(actual.contains(">1000000000000.02</text>"));
}

#[test]
fn long_group_identifiers_reject_overlapping_labels() {
    let mut d = rows(&[0., 1., 2., 3.]);
    for (i, r) in d.iter_mut().enumerate() {
        r[1] = 9_007_199_254_740_990. - i as f64;
    }
    assert!(execute("report_box_plot_svg", &d, &[1., 400., 600., 0., 0., 1., 0.]).is_err());
    assert!(execute("report_box_plot_svg", &d, &[1., 400., 600., 0., 1., 1., 0.]).is_ok());
}

#[test]
fn exact_and_adjacent_bin_edges_and_million_frequency_remain_bounded() {
    let below_one = f64::from_bits(1f64.to_bits() - 1);
    let above_one = f64::from_bits(1f64.to_bits() + 1);
    let below_two = f64::from_bits(2f64.to_bits() - 1);
    let j = hist(
        &rows(&[0., below_one, 1., above_one, below_two, 2.]),
        1.,
        0.,
        &[0., 1., 2.],
    );
    assert_eq!(nums(&j, "counts"), vec![2., 4.]);
    let mut d = rows(&[42.]);
    d[0][2] = 1_000_000.;
    let j = bx(&d, 0., 1., 1., 0.);
    assert!(j.contains("\"source_rows\":1,\"selected_rows\":1,\"N\":1000000"));
    assert_eq!(nums(&j, "five_number"), vec![42.; 5]);
}

#[test]
fn stable_geometry_origin_is_explicit() {
    let d = rows(&[1e12, 1e12 + 0.01, 1e12 + 0.02]);
    let j = bx(&d, 0., 0., 1., 0.);
    assert!(j.contains("\"value_origin\":1000000000000"));
    assert!(j.contains("\"svg_origin\":"));
}
#[test]
fn accessible_group_statistics_are_explicit_without_visible_table() {
    let j = bx(&rows(&[1., 2., 3., 4.]), 0., 0., 1., 0.);
    assert!(j.contains("Minimum 1; Q1 1.25; median 2.5; Q3 3.75; maximum 4"));
}
#[test]
fn empty_group_statistics_are_explicit_null() {
    let mut d = rows(&[1., 2.]);
    d[1][1] = 2.;
    d[1][3] = 1.;
    let j = bx(&d, 0., 0., 1., 0.);
    assert!(j.contains("\"whiskers\":null"));
    assert!(j.contains("\"quantile_ticks\":null"));
}

#[test]
fn visible_five_number_table_preserves_distinct_large_offset_values() {
    let values: Vec<f64> = [0., 0.005, 0.01, 0.015, 0.02]
        .iter()
        .map(|v| 1e12 + v)
        .collect();
    let j = execute(
        "report_box_plot_svg",
        &rows(&values),
        &[1., 1000., 600., 0., 0., 1., 1.],
    )
    .unwrap();
    let after_header = j.split_once(">Maximum</text>").unwrap().1;
    let numbers: Vec<f64> = after_header
        .split("<text ")
        .skip(1)
        .filter_map(|s| {
            s.split_once('>')
                .unwrap()
                .1
                .split("</text>")
                .next()
                .unwrap()
                .parse()
                .ok()
        })
        .collect();
    let visible = &numbers[numbers.len() - 5..];
    assert!(
        visible.windows(2).all(|p| p[0] < p[1]),
        "five different values collapsed in visible table: {visible:?}"
    );
    assert_eq!(visible, nums(&j, "five_number"));
}
