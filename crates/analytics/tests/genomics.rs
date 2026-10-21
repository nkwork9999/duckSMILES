use ducksmiles_analytics::execute as run;
fn values(j: &str, k: &str) -> Vec<f64> {
    let t = j.split_once(&format!("\"{k}\":")).expect(k).1;
    let end = if t.starts_with('[') {
        let mut d = 0;
        t.char_indices()
            .find_map(|(i, c)| {
                if c == '[' {
                    d += 1;
                }
                if c == ']' {
                    d -= 1;
                    if d == 0 {
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
        .filter(|v| !v.is_empty())
        .map(|v| v.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1.), "{a} != {b}");
}
fn vector(j: &str, k: &str, w: &[f64]) {
    let a = values(j, k);
    assert_eq!(a.len(), w.len());
    for (&x, &y) in a.iter().zip(w) {
        near(x, y);
    }
}
#[test]
fn marker_statistics_recovers_diploid_counts_diversity_and_pearson_hwe() {
    let j = run(
        "marker_statistics",
        &[vec![0., 0.], vec![1., 1.], vec![1., 2.], vec![2., -1.]],
        &[],
    )
    .unwrap();
    vector(&j, "called_count", &[4., 3.]);
    vector(&j, "genotype_counts", &[1., 2., 1., 1., 1., 1.]);
    vector(&j, "allele_frequency", &[0.5, 0.5]);
    vector(&j, "minor_allele_frequency", &[0.5, 0.5]);
    vector(&j, "observed_heterozygosity", &[0.5, 1. / 3.]);
    vector(&j, "polymorphism_information_content", &[0.375, 0.375]);
    vector(&j, "hwe_chi_square", &[0., 1. / 3.]);
    vector(&j, "hwe_p_value", &[1., 0.5637028616507731]);
}
#[test]
fn mean_and_mode_imputation_save_training_values_for_new_rows() {
    let data = vec![vec![0., -1.], vec![1., 1.], vec![-1., 2.]];
    let mean = run("marker_impute", &data, &[0.]).unwrap();
    vector(&mean, "rows", &[0., 1.5, 1., 1., 0.5, 2.]);
    vector(&mean, "predict_parameters", &[0.5, 1.5]);
    let pred = run(
        "marker_impute_predict",
        &[vec![-1., -1.], vec![2., 0.]],
        &values(&mean, "predict_parameters"),
    )
    .unwrap();
    vector(&pred, "rows", &[0.5, 1.5, 2., 0.]);
    let mode = run("marker_impute", &data, &[1.]).unwrap();
    vector(&mode, "rows", &[0., 1., 1., 1., 0., 2.]);
}
#[test]
fn ld_knn_uses_correlated_markers_and_nearest_observed_sample() {
    let j = run(
        "marker_impute",
        &[vec![0., 0.], vec![1., 1.], vec![2., 2.], vec![-1., 2.]],
        &[2., 1., 1.],
    )
    .unwrap();
    vector(&j, "rows", &[0., 0., 1., 1., 2., 2., 2., 2.]);
}
#[test]
fn additive_and_dominance_relationships_have_known_vanraden_kernels() {
    let data = vec![vec![0., 0., 2.], vec![1., 1., 2.], vec![2., 2., 2.]];
    for method in [0., 1.] {
        let j = run("marker_relatedness", &data, &[method]).unwrap();
        vector(&j, "matrix", &[2., 0., -2., 0., 0., 0., -2., 0., 2.]);
    }
    let d = run("marker_relatedness", &data, &[2.]).unwrap();
    vector(&d, "matrix", &[1., -1., 1., -1., 1., -1., 1., -1., 1.]);
}
#[test]
fn distance_matrix_implements_numeric_euclidean_jaccard_bray_gower_and_hamming() {
    let data = vec![vec![0., 1.], vec![2., 1.]];
    for (method, want) in [(0., 2.), (1., 0.5), (2., 0.5), (3., 0.5), (4., 0.5)] {
        let j = run("distance_matrix", &data, &[method]).unwrap();
        vector(&j, "matrix", &[0., want, want, 0.]);
    }
    let j = run(
        "distance_matrix",
        &[vec![-1., 2.], vec![1., 0.]],
        &[0., -1.],
    )
    .unwrap();
    vector(&j, "matrix", &[0., 2., 2., 0.]);
}
#[test]
fn seeded_genotype_simulation_obeys_allele_frequency_and_inbreeding() {
    let a = run(
        "marker_simulate",
        &[vec![0.25], vec![0.5]],
        &[10000., 42., 0.5],
    )
    .unwrap();
    assert_eq!(
        a,
        run(
            "marker_simulate",
            &[vec![0.25], vec![0.5]],
            &[10000., 42., 0.5]
        )
        .unwrap()
    );
    let rows = values(&a, "rows");
    assert_eq!(rows.len(), 20000);
    for (j, p) in [(0, 0.25), (1, 0.5)] {
        let x: Vec<_> = rows.chunks(2).map(|r| r[j]).collect();
        assert!(x.iter().all(|v| [0., 1., 2.].contains(v)));
        assert!((x.iter().sum::<f64>() / 10000. - 2. * p).abs() < 0.04);
        assert!(
            (x.iter().filter(|&&v| v == 1.).count() as f64 / 10000. - p * (1. - p)).abs() < 0.03
        );
    }
}
#[test]
fn unlinked_parental_cross_preserves_mendelian_fixed_parent_genotypes() {
    let j = run(
        "marker_cross",
        &[vec![0., 2., 1.], vec![2., 0., 1.]],
        &[0., 1., 200., 42.],
    )
    .unwrap();
    let x = values(&j, "rows");
    for r in x.chunks(3) {
        assert_eq!(r[0], 1.);
        assert_eq!(r[1], 1.);
        assert!([0., 1., 2.].contains(&r[2]));
    }
}
#[test]
fn genomic_blup_known_variance_ratio_matches_ridge_and_reuses_training_preprocessing() {
    let mut d = vec![vec![0., 1.], vec![1., 2.], vec![2., 5.], vec![1., 999.]];
    let j = run("genomic_blup", &d, &[3., 1.]).unwrap();
    vector(&j, "marker_effects", &[1.6]);
    vector(&j, "predictions", &[8. / 3.]);
    vector(&j, "fitted", &[8. / 3. - 1.6, 8. / 3., 8. / 3. + 1.6]);
    d[3][1] = -999.;
    assert_eq!(j, run("genomic_blup", &d, &[3., 1.]).unwrap());
    let p = run(
        "genomic_blup_predict",
        &[vec![0.], vec![-1.], vec![2.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    vector(&p, "predictions", &[8. / 3. - 1.6, 8. / 3., 8. / 3. + 1.6]);
}
#[test]
fn genetics_rejects_invalid_calls_all_missing_loci_and_unbounded_pair_work() {
    assert!(run("marker_statistics", &[vec![3.]], &[]).is_err());
    assert!(run("marker_impute", &[vec![-1.], vec![-1.]], &[0.]).is_err());
    assert!(run("marker_relatedness", &[vec![2.], vec![2.]], &[]).is_err());
    assert!(run("marker_simulate", &[vec![1.5]], &[3., 42.]).is_err());
    assert!(run("genomic_blup", &[vec![0., 1.], vec![2., 2.]], &[2., 0.]).is_err());
    assert!(run("distance_matrix", &vec![vec![0.; 256]; 1000], &[]).is_err());
}

#[test]
fn stochastic_imputation_respects_hwe_observed_uniform_and_specified_policies() {
    let mut data = vec![vec![0.], vec![0.], vec![2.], vec![2.]];
    data.extend(vec![vec![-1.]; 20000]);
    for (method, want) in [
        (3., [0.25, 0.5, 0.25]),
        (4., [0.5, 0., 0.5]),
        (5., [1. / 3.; 3]),
        (6., [0., 1., 0.]),
    ] {
        let j = run("marker_impute", &data, &[method, 10., 5., 42., 1.]).unwrap();
        let rows = values(&j, "rows");
        assert_eq!(&rows[..4], &[0., 0., 2., 2.]);
        for (g, expected) in want.iter().enumerate() {
            let proportion = rows[4..].iter().filter(|&&v| v == g as f64).count() as f64 / 20000.;
            assert!(
                (proportion - expected).abs() < 0.02,
                "method{method} g{g}: {proportion} != {expected}"
            );
        }
        assert_eq!(
            j,
            run("marker_impute", &data, &[method, 10., 5., 42., 1.]).unwrap()
        );
    }
}

#[test]
fn simulation_keeps_certain_alleles_even_when_the_random_midpoint_rounds_to_one() {
    // Independently inverted SplitMix64 output with its upper53bits all one.
    // This valid exact seed exposed the f64 midpoint rounding to 1.0.
    let seed = 608_688_947_055_533.;
    let j = run("marker_simulate", &[vec![0.]], &[1., seed, 0.]).unwrap();
    vector(&j, "rows", &[0.]);
    let j = run("marker_cross", &[vec![2.]], &[0., 0., 1., seed]).unwrap();
    vector(&j, "rows", &[2.]);
}
