#[path = "../../src/text_analysis.rs"]
mod text_analysis;
use text_analysis::run;
fn values(j: &str, k: &str) -> Vec<f64> {
    let t = j.split_once(&format!("\"{k}\":")).expect(k).1;
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
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(x: f64, y: f64, t: f64) {
    assert!((x - y).abs() < t * (1. + y.abs()), "{x} != {y}");
}
fn counts() -> Vec<Vec<f64>> {
    vec![vec![2., 1., 0.], vec![0., 1., 3.], vec![0., 0., 0.]]
}
#[test]
fn stacked_pretokenized_counts_preserve_empty_documents_and_repeated_terms() {
    let j = run(
        "text_dtm",
        &[vec![0., 2., 1.], vec![0., 2., 2.], vec![2., 0., 4.]],
        &[4., 3.],
    )
    .unwrap();
    assert_eq!(
        values(&j, "counts"),
        vec![0., 0., 3., 0., 0., 0., 4., 0., 0., 0., 0., 0.]
    );
    assert_eq!(values(&j, "document_frequency"), vec![1., 0., 1.]);
}
#[test]
fn jmp_tfidf_is_unsmoothed_base_ten_with_frequency_not_document_normalization() {
    let j = run("text_weight", &counts(), &[4.]).unwrap();
    let m = values(&j, "matrix");
    near(m[1], 2. * 3_f64.log10(), 1e-14);
    near(m[3], 3. * 3_f64.log10(), 1e-14);
    assert_eq!(values(&j, "selected_terms"), vec![2., 0., 1.]);
    let binary = run("text_weight", &counts(), &[0.]).unwrap();
    assert_eq!(
        values(&binary, "matrix"),
        vec![0., 1., 1., 1., 0., 1., 0., 0., 0.]
    );
}
#[test]
fn saved_weighting_uses_training_idf_without_fitting_query_documents() {
    let j = run("text_weight", &counts(), &[4.]).unwrap();
    let q = run(
        "text_weight_predict",
        &[vec![0., 4., 2.]],
        &values(&j, "weight_parameters"),
    )
    .unwrap();
    let m = values(&q, "matrix");
    near(m[0], 2. * 3_f64.log10(), 1e-14);
    near(m[2], 4. * 1.5_f64.log10(), 1e-14);
}
fn diagonal() -> Vec<Vec<f64>> {
    vec![vec![3., 0.], vec![0., 2.], vec![0., 0.]]
}
#[test]
fn lsa_returns_normalized_svd_document_scores_and_reconstruction() {
    let j = run("text_lsa", &diagonal(), &[2., 0.]).unwrap();
    let s = values(&j, "singular_values");
    near(s[0], 1., 1e-13);
    near(s[1], 2. / 3., 1e-13);
    near(values(&j, "relative_reconstruction_error")[0], 0., 1e-13);
}
#[test]
fn saved_projection_reuses_training_center_scale_and_corpus_divisor() {
    let j = run("text_lsa", &diagonal(), &[2., 0.]).unwrap();
    let q = run(
        "text_projection",
        &diagonal(),
        &values(&j, "projection_parameters"),
    )
    .unwrap();
    for (a, b) in values(&q, "scores")
        .into_iter()
        .zip(values(&j, "document_scores"))
    {
        near(a, b, 1e-13);
    }
}
#[test]
fn rotated_svd_topics_preserve_rank_reconstruction_and_return_true_rotation() {
    let j = run("text_topic_svd", &diagonal(), &[2., 0., 100., 1e-10, 1.]).unwrap();
    assert!(j.contains("\"rotation_converged\":true"));
    let r = values(&j, "rotation");
    near(r[0] * r[0] + r[2] * r[2], 1., 1e-13);
    near(values(&j, "relative_reconstruction_error")[0], 0., 1e-13);
    let q = run(
        "text_projection",
        &diagonal(),
        &values(&j, "projection_parameters"),
    )
    .unwrap();
    for (a, b) in values(&q, "scores")
        .into_iter()
        .zip(values(&j, "document_scores"))
    {
        near(a, b, 1e-13);
    }
}
#[test]
fn one_class_bernoulli_model_is_empirical_term_presence_and_true_likelihood() {
    let j = run("text_latent_class", &counts(), &[1., 0., 200., 1e-10, 1e-8]).unwrap();
    let p = values(&j, "term_probabilities");
    near(p[0], 1. / 3., 1e-13);
    near(p[1], 2. / 3., 1e-13);
    near(p[2], 1. / 3., 1e-13);
    let ll = 3. * (1. / 3_f64).ln() + 6. * (2. / 3_f64).ln();
    near(values(&j, "log_likelihood")[0], ll, 1e-13);
    assert!(j.contains("\"converged\":true"));
}
#[test]
fn saved_latent_class_posterior_includes_both_present_and_absent_terms() {
    let params = [1., 2., 2., 1e-8, 0.5, 0.5, 0.9, 0.1, 0.1, 0.9];
    let j = run(
        "text_latent_class_predict",
        &[vec![1., 0.], vec![0., 1.], vec![0., 0.]],
        &params,
    )
    .unwrap();
    let p = values(&j, "probabilities");
    near(p[0], 81. / 82., 1e-13);
    near(p[3], 81. / 82., 1e-13);
    near(p[4], 0.5, 1e-13);
}
#[test]
fn numeric_token_and_dense_matrix_contracts_reject_malformed_inputs() {
    for (op, d, a) in [
        ("text_dtm", vec![vec![0., 3., 1.]], vec![2., 3.]),
        ("text_weight", vec![vec![-1.]], vec![]),
        ("text_lsa", vec![vec![0., 0.], vec![0., 0.]], vec![1., 0.]),
        ("text_projection", vec![vec![1.]], vec![]),
    ] {
        assert!(run(op, &d, &a).is_err());
    }
}
#[test]
fn weighting_saved_document_count_supports_every_accepted_training_corpus_size() {
    let data = vec![vec![1.]; 4097];
    let j = run("text_weight", &data, &[2.]).unwrap();
    let q = run(
        "text_weight_predict",
        &[vec![2.]],
        &values(&j, "weight_parameters"),
    )
    .unwrap();
    near(values(&q, "matrix")[0], 2., 1e-14);
}
#[test]
fn saved_models_bound_expanded_output_cells_independently_of_input_width() {
    let data = vec![vec![0.]; 65537];
    assert!(run("text_projection", &data, &[1., 1., 2., 1., 0., 1., 1., 1.]).is_err());
    assert!(run(
        "text_latent_class_predict",
        &data,
        &[1., 2., 1., 1e-8, 0.5, 0.5, 0.5, 0.5]
    )
    .is_err());
}
#[test]
fn centered_scaled_lsa_handles_constant_columns_rank_loss_and_tiny_units() {
    let data = vec![
        vec![1., 5., 4.],
        vec![2., 5., 3.],
        vec![3., 5., 2.],
        vec![4., 5., 1.],
    ];
    let j = run("text_lsa", &data, &[1., 2.]).unwrap();
    assert_eq!(values(&j, "numerical_rank"), vec![1.]);
    assert!(run("text_lsa", &data, &[2., 2.]).is_err());
    let tiny = data
        .iter()
        .map(|r| r.iter().map(|v| v * 1e-300).collect())
        .collect::<Vec<Vec<f64>>>();
    let t = run("text_lsa", &tiny, &[1., 2.]).unwrap();
    near(
        values(&j, "singular_values")[0],
        values(&t, "singular_values")[0],
        1e-13,
    );
    let scored = run(
        "text_projection",
        &tiny,
        &values(&t, "projection_parameters"),
    )
    .unwrap();
    for (a, b) in values(&scored, "scores")
        .into_iter()
        .zip(values(&t, "document_scores"))
    {
        near(a, b, 1e-13);
    }
}
#[test]
fn ubiquitous_tfidf_terms_have_zero_weight_and_do_not_create_false_singular_vectors() {
    let j = run("text_weight", &[vec![1., 2.], vec![3., 1.]], &[4.]).unwrap();
    assert!(values(&j, "matrix").iter().all(|&v| v == 0.));
    assert!(run("text_lsa", &[vec![0., 0.], vec![0., 0.]], &[1., 0.]).is_err());
}
#[test]
fn collapsed_and_unfinished_mixtures_do_not_claim_convergence_or_model_selection() {
    let same = vec![vec![1., 0., 1.]; 8];
    let j = run("text_latent_class", &same, &[2., 3., 100., 1e-10]).unwrap();
    assert!(j.contains("\"components_distinct\":false"));
    assert!(j.contains("\"converged\":false"));
    assert!(j.contains("\"bic\":null"));
    let j = run("text_latent_class", &counts(), &[1., 0., 1.]).unwrap();
    assert!(j.contains("\"optimizer_converged\":false"));
    assert!(j.contains("\"aic\":null"));
}
#[test]
fn full_rank_lsa_reconstruction_error_measures_residual_instead_of_cancelled_energies() {
    let data = vec![
        vec![4., 2., 3.],
        vec![4., 1., 2.],
        vec![3., 0., 1.],
        vec![3., 3., 0.],
        vec![7., 4., 1.],
        vec![0., 3., 5.],
        vec![3., 3., 4.],
    ];
    let j = run("text_lsa", &data, &[3., 0.]).unwrap();
    assert!(values(&j, "relative_reconstruction_error")[0] < 1e-12);
}
#[test]
fn varimax_rotation_is_invariant_to_representable_uniform_dtm_units() {
    let data = vec![
        vec![4., 2., 3.],
        vec![4., 1., 2.],
        vec![3., 0., 1.],
        vec![3., 3., 0.],
        vec![7., 4., 1.],
        vec![0., 3., 5.],
        vec![3., 3., 4.],
    ];
    let j = run("text_topic_svd", &data, &[2., 0., 1000., 1e-12, 0.]).unwrap();
    let tiny = data
        .iter()
        .map(|r| r.iter().map(|v| v * 1e-100).collect())
        .collect::<Vec<Vec<f64>>>();
    let t = run("text_topic_svd", &tiny, &[2., 0., 1000., 1e-12, 0.]).unwrap();
    for (a, b) in values(&j, "rotation")
        .into_iter()
        .zip(values(&t, "rotation"))
    {
        near(a, b, 1e-11);
    }
}
