use ducksmiles_analytics::execute as run;

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
fn treatments(n: usize) -> Vec<Vec<f64>> {
    (0..n).map(|i| vec![i as f64]).collect()
}
fn check_blocks(j: &str, v: usize, k: usize, r: usize, lambda: usize) {
    let treatments = numbers(j, "treatment_ids");
    let blocks = numbers(j, "block_ids");
    let mut incidence = std::collections::BTreeMap::<usize, Vec<usize>>::new();
    for (t, b) in treatments.iter().zip(&blocks) {
        incidence.entry(*b as usize).or_default().push(*t as usize);
    }
    let mut counts = vec![0; v];
    let mut pair = vec![vec![0; v]; v];
    for b in incidence.values() {
        assert_eq!(b.len(), k);
        for (i, &t) in b.iter().enumerate() {
            counts[t] += 1;
            for &u in &b[..i] {
                assert_ne!(t, u);
                pair[t][u] += 1;
                pair[u][t] += 1;
            }
        }
    }
    assert!(counts.iter().all(|&x| x == r));
    for (i, row) in pair.iter().enumerate() {
        for (j, &x) in row.iter().enumerate() {
            if i != j {
                assert_eq!(x, lambda);
            }
        }
    }
}
#[test]
fn randomized_complete_blocks_contain_each_treatment_once_with_contiguous_blocks() {
    let j = run("rcbd_design", &treatments(5), &[3., 42.]).unwrap();
    check_blocks(&j, 5, 5, 3, 3);
    let b = numbers(&j, "block_ids");
    for chunk in b.chunks(5) {
        assert!(chunk.iter().all(|x| *x == chunk[0]));
    }
    all(
        &numbers(&j, "run_order"),
        &(1..=15).map(|i| i as f64).collect::<Vec<_>>(),
    );
}
#[test]
fn fano_incomplete_blocks_have_exact_pair_concurrence_and_replication() {
    let j = run("balanced_incomplete_block", &treatments(7), &[3., 1., 91.]).unwrap();
    check_blocks(&j, 7, 3, 3, 1);
    near(numbers(&j, "n_blocks")[0], 7.);
    assert!(j.contains("\"search_complete\":true"));
}
#[test]
fn higher_lambda_bibd_preserves_each_pair_and_all_treatments() {
    let j = run("balanced_incomplete_block", &treatments(6), &[3., 2., 5.]).unwrap();
    check_blocks(&j, 6, 3, 5, 2);
}
#[test]
fn resolvable_affine_blocks_partition_treatments_in_each_parallel_class() {
    let j = run("resolvable_incomplete_block", &treatments(9), &[13.]).unwrap();
    check_blocks(&j, 9, 3, 4, 1);
    let cl = numbers(&j, "parallel_class_ids");
    let t = numbers(&j, "treatment_ids");
    let mut groups = std::collections::BTreeMap::<usize, Vec<usize>>::new();
    for (c, t) in cl.iter().zip(t) {
        groups.entry(*c as usize).or_default().push(t as usize);
    }
    assert_eq!(groups.len(), 4);
    for ids in groups.values_mut() {
        ids.sort();
        assert_eq!(ids, &(0..9).collect::<Vec<_>>());
    }
}
#[test]
fn split_plot_whole_factors_stay_constant_and_subplots_randomize_within_whole_plot() {
    let bounds = vec![vec![-2., 2., 0.], vec![10., 20., 1.], vec![0., 4., 1.]];
    let j = run("split_plot_design", &bounds, &[3., 73.]).unwrap();
    let d = numbers(&j, "design");
    let wp = numbers(&j, "whole_plot_ids");
    let rep = numbers(&j, "replicate_ids");
    assert_eq!(wp.len(), 24);
    for c in 0..6 {
        let ix = c * 4;
        assert!(wp[ix..ix + 4].iter().all(|x| *x == wp[ix]));
        assert!(d[ix * 3..(ix + 4) * 3].chunks(3).all(|r| r[0] == d[ix * 3]));
    }
    for r in 0..3 {
        let mut combos = std::collections::BTreeSet::new();
        for (i, rr) in rep.iter().enumerate() {
            if *rr == r as f64 {
                combos.insert(
                    d[i * 3..i * 3 + 3]
                        .iter()
                        .map(|x| *x as i64)
                        .collect::<Vec<_>>(),
                );
            }
        }
        assert_eq!(combos.len(), 8);
    }
}
#[test]
fn strip_plot_has_crossed_unit_ids_and_two_valid_stage_orders() {
    let j = run(
        "strip_plot_design",
        &[vec![0., 1., 0.], vec![2., 3., 1.], vec![4., 5., 1.]],
        &[2., 29.],
    )
    .unwrap();
    let rows = numbers(&j, "row_strip_ids");
    let cols = numbers(&j, "column_strip_ids");
    let stage = numbers(&j, "stage2_order");
    let mut sorted = stage.clone();
    sorted.sort_by(f64::total_cmp);
    all(&sorted, &(0..16).map(|i| i as f64).collect::<Vec<_>>());
    let mut cells = std::collections::BTreeSet::new();
    for (&r, &c) in rows.iter().zip(&cols) {
        assert!(cells.insert((r as usize, c as usize)));
    }
    for group in stage.chunks(2) {
        let first = cols[group[0] as usize];
        assert!(group.iter().all(|i| cols[*i as usize] == first));
    }
    assert_eq!(cells.len(), 16);
}
#[test]
fn mixed_level_l18_is_pairwise_orthogonal_and_l4_reports_interaction_aliasing() {
    let levels = vec![
        vec![2.],
        vec![3.],
        vec![3.],
        vec![3.],
        vec![3.],
        vec![3.],
        vec![3.],
        vec![3.],
    ];
    let j = run("taguchi_orthogonal_array", &levels, &[18., 11.]).unwrap();
    let d = numbers(&j, "design");
    assert_eq!(d.len(), 144);
    for a in 0..8 {
        for b in 0..a {
            let la = levels[a][0] as usize;
            let lb = levels[b][0] as usize;
            let mut pairs = vec![0; la * lb];
            for row in d.chunks(8) {
                pairs[row[a] as usize * lb + row[b] as usize] += 1;
            }
            assert!(pairs.iter().all(|x| *x == 18 / (la * lb)));
        }
    }
    let gram = numbers(&j, "main_effect_gram");
    for i in 0..16 {
        for k in 0..16 {
            near(gram[i * 16 + k], if i == k { 1. } else { 0. });
        }
    }
    let j = run(
        "taguchi_orthogonal_array",
        &[vec![2.], vec![2.], vec![2.]],
        &[4., 0.],
    )
    .unwrap();
    near(numbers(&j, "max_main_interaction_alias")[0], 1.);
}
#[test]
fn seeds_are_reproducible_and_invalid_or_unsupported_designs_are_rejected() {
    let data = treatments(7);
    let a = run("balanced_incomplete_block", &data, &[3., 1., 42.]).unwrap();
    assert_eq!(
        a,
        run("balanced_incomplete_block", &data, &[3., 1., 42.]).unwrap()
    );
    assert_ne!(
        a,
        run("balanced_incomplete_block", &data, &[3., 1., 43.]).unwrap()
    );
    assert!(run("balanced_incomplete_block", &treatments(6), &[3., 1.]).is_err());
    assert!(run("resolvable_incomplete_block", &treatments(16), &[]).is_err());
    assert!(run("split_plot_design", &[vec![0., 1., 0.]], &[]).is_err());
    assert!(run("taguchi_orthogonal_array", &[vec![2.], vec![2.]], &[18.]).is_err());
    assert!(run("rcbd_design", &treatments(2), &[100_000.]).is_err());
}
#[test]
fn repeated_block_bibd_is_valid_when_one_candidate_must_cover_pair_more_than_once() {
    let j = run("balanced_incomplete_block", &treatments(3), &[2., 3., 4.]).unwrap();
    check_blocks(&j, 3, 2, 6, 3);
}
#[test]
fn unreplicated_split_and_strip_plans_report_zero_error_df() {
    let data = [vec![0., 1., 0.], vec![0., 1., 1.]];
    for op in ["split_plot_design", "strip_plot_design"] {
        let j = run(op, &data, &[]).unwrap();
        assert!(j.contains("\"inference_ready\":false"));
        assert!(numbers(&j, "replication_error_df").iter().all(|x| *x == 0.));
    }
}
#[test]
fn affine_eight_and_cyclic_nine_treatment_bibds_do_not_need_exponential_search() {
    for v in [8, 9] {
        let j = run("balanced_incomplete_block", &treatments(v), &[4., 3., 37.]).unwrap();
        check_blocks(&j, v, 4, v - 1, 3);
    }
}
