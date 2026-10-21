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
fn near(x: f64, y: f64) {
    assert!((x - y).abs() < 2e-7 * y.abs().max(1.), "{x} != {y}");
}
fn vector(j: &str, k: &str, w: &[f64]) {
    let got = values(j, k);
    assert_eq!(got.len(), w.len());
    for (&a, &b) in got.iter().zip(w) {
        near(a, b);
    }
}

#[test]
fn count_normalization_distinguishes_cpm_rpkm_tpm_and_relative_abundance() {
    let rows = vec![vec![1., 3.], vec![2., 6.]];
    for (method, args, want) in [
        (1., vec![], vec![250000., 750000.]),
        (2., vec![1000., 2000.], vec![250000., 375000.]),
        (3., vec![1000., 2000.], vec![400000., 600000.]),
        (5., vec![], vec![25., 75.]),
        (6., vec![], vec![0.25, 0.75]),
    ] {
        let mut p = vec![method];
        p.extend(args);
        let j = run("assay_normalize", &rows, &p).unwrap();
        vector(
            &j,
            "rows",
            &want.iter().chain(&want).copied().collect::<Vec<_>>(),
        );
    }
}
#[test]
fn row_standardization_and_clr_preserve_mean_centering_and_explicit_pseudocount() {
    let j = run("assay_normalize", &[vec![1., 3.], vec![2., 6.]], &[0.]).unwrap();
    let z = 1. / 2_f64.sqrt();
    vector(&j, "rows", &[-z, z, -z, z]);
    let j = run("assay_normalize", &[vec![1., 4.]], &[4.]).unwrap();
    vector(&j, "rows", &[-2_f64.ln(), 2_f64.ln()]);
    let j = run("assay_normalize", &[vec![0., 3.]], &[4., 1.]).unwrap();
    vector(&j, "rows", &[-2_f64.ln(), 2_f64.ln()]);
}
#[test]
fn row_standardization_is_invariant_when_measurement_units_are_subnormal() {
    // Centering or forming the SD in physical units loses fractional subnormals.
    let tiny = f64::from_bits(1);
    let j = run("assay_normalize", &[vec![0., 0., 0., tiny]], &[0.]).unwrap();
    vector(&j, "rows", &[-0.5, -0.5, -0.5, 1.5]);
    let j = run("assay_normalize", &[vec![0., 0., 0., tiny]], &[0., 0.]).unwrap();
    let s = 3_f64.sqrt();
    vector(&j, "rows", &[-1. / s, -1. / s, -1. / s, s]);
}
#[test]
fn quantile_normalization_uses_equal_rank_means_and_tie_averaging() {
    let j = run(
        "quantile_normalize",
        &[vec![5., 1., 3.], vec![4., 2., 6.]],
        &[],
    )
    .unwrap();
    vector(&j, "reference_quantiles", &[1.5, 3.5, 5.5]);
    vector(&j, "rows", &[5.5, 1.5, 3.5, 3.5, 1.5, 5.5]);
    let j = run(
        "quantile_normalize",
        &[vec![1., 1., 5.], vec![2., 4., 6.]],
        &[],
    )
    .unwrap();
    vector(&j, "rows", &[2., 2., 5.5, 1.5, 2.5, 5.5]);
    let p = run(
        "quantile_normalize_predict",
        &[vec![4., 4., 0.]],
        &values(&j, "reference_quantiles"),
    )
    .unwrap();
    vector(&p, "rows", &[4., 4., 1.5]);
}
#[test]
fn single_ancestral_population_has_closed_form_allele_mle_and_full_binomial_likelihood() {
    let j = run(
        "marker_admixture",
        &[vec![0., 0.], vec![1., 2.], vec![2., -1.]],
        &[1., 1000., 1e-9, 42., 1.],
    )
    .unwrap();
    vector(&j, "ancestry", &[1., 1., 1.]);
    vector(&j, "allele_frequencies", &[0.5, 0.5]);
    near(values(&j, "log_likelihood")[0], -9. * 2_f64.ln());
    let p = run(
        "admixture_predict",
        &[vec![0., 2.], vec![1., 1.]],
        &values(&j, "predict_parameters"),
    )
    .unwrap();
    vector(&p, "ancestry", &[1., 1.]);
}
#[test]
fn fixed_ancestral_frequencies_infer_new_rows_without_refitting_allele_frequencies() {
    // Portable [K,p,maxiter,tolerance,floor,initial_Q(K),row_major_F(K*p)].
    let params = [2., 2., 10000., 1e-10, 1e-8, 0.5, 0.5, 0.1, 0.9, 0.9, 0.1];
    let j = run(
        "admixture_predict",
        &[vec![0., 2.], vec![1., 1.], vec![2., 0.]],
        &params,
    )
    .unwrap();
    vector(&j, "ancestry", &[1., 0., 0.5, 0.5, 0., 1.]);
}
#[test]
fn saved_ancestry_initialization_cannot_lock_a_population_out_of_inference() {
    let params = [2., 2., 10000., 1e-10, 1e-8, 1., 0., 0.1, 0.9, 0.9, 0.1];
    // A zero EM starting weight never revives, even for a query strongly
    // favoring the omitted population. Reject such a malformed saved model.
    assert!(run("admixture_predict", &[vec![2., 0.]], &params).is_err());
}
#[test]
fn admixture_iteration_limit_reports_nonconvergence_and_scoring_is_batch_independent() {
    let rows = vec![
        vec![0., 2.],
        vec![0., 2.],
        vec![1., 1.],
        vec![2., 0.],
        vec![2., 0.],
    ];
    let j = run("marker_admixture", &rows, &[2., 1., 1e-12, 42., 1.]).unwrap();
    assert!(j.contains("\"converged\":false"));
    vector(&j, "iterations", &[1.]);
    let params = [2., 2., 10000., 1e-10, 1e-8, 0.5, 0.5, 0.1, 0.9, 0.9, 0.1];
    let batch = run("admixture_predict", &[vec![0., 2.], vec![2., 0.]], &params).unwrap();
    let first = run("admixture_predict", &[vec![0., 2.]], &params).unwrap();
    vector(&first, "ancestry", &values(&batch, "ancestry")[..2]);
}
#[test]
fn count_and_log_ratio_normalization_remain_finite_for_extreme_units() {
    let huge = f64::MAX;
    let j = run("assay_normalize", &[vec![huge, huge]], &[1.]).unwrap();
    vector(&j, "rows", &[500000., 500000.]);
    let j = run("assay_normalize", &[vec![huge, huge]], &[4., huge]).unwrap();
    vector(&j, "rows", &[0., 0.]);
    let j = run("assay_normalize", &[vec![huge, huge]], &[3., huge, huge]).unwrap();
    vector(&j, "rows", &[500000., 500000.]);
}
#[test]
fn admixture_estimation_is_seeded_and_reports_monotone_finite_likelihood() {
    let rows = vec![
        vec![0., 2.],
        vec![0., 2.],
        vec![1., 1.],
        vec![2., 0.],
        vec![2., 0.],
    ];
    let params = [2., 1000., 1e-8, 42., 3.];
    let j = run("marker_admixture", &rows, &params).unwrap();
    assert_eq!(j, run("marker_admixture", &rows, &params).unwrap());
    let q = values(&j, "ancestry");
    for r in q.chunks(2) {
        near(r.iter().sum(), 1.);
        assert!(r.iter().all(|&v| (0. ..=1.).contains(&v)));
    }
    let history = values(&j, "log_likelihood_history");
    assert!(history.len() > 1);
    assert!(history.windows(2).all(|w| w[1] >= w[0] - 1e-9));
    // At separated ancestry maximum: eight homozygous calls have probability~1,
    // two heterozygous calls have probability.5.
    near(values(&j, "log_likelihood")[0], -2. * 2_f64.ln());
}
#[test]
fn omics_rejects_undefined_normalization_and_unidentified_or_unbounded_models() {
    assert!(run("assay_normalize", &[vec![0., 0.]], &[1.]).is_err());
    assert!(run("assay_normalize", &[vec![0., 1.]], &[4.]).is_err());
    assert!(run("assay_normalize", &[vec![1., 2.]], &[2., 0., 1.]).is_err());
    assert!(run("quantile_normalize_predict", &[vec![1., 2.]], &[2., 1.]).is_err());
    assert!(run("marker_admixture", &[vec![-1., 0.], vec![-1., 2.]], &[2.]).is_err());
    assert!(run("marker_admixture", &[vec![0.], vec![2.]], &[1000.]).is_err());
}
