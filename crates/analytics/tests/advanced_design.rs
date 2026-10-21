use ducksmiles_analytics::execute as run;

fn values(j: &str, key: &str) -> Vec<f64> {
    let start = j.find(&format!("\"{key}\":")).expect(key) + key.len() + 3;
    let tail = &j[start..];
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
        tail.find([',', '}']).unwrap_or(tail.len())
    };
    tail[..end]
        .split(|c: char| !c.is_ascii_digit() && !matches!(c, '-' | '+' | '.' | 'e' | 'E'))
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().unwrap())
        .collect()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 2e-9 * b.abs().max(1.), "{a} != {b}");
}
fn near_vec(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        near(x, y);
    }
}
fn bounds(p: usize) -> Vec<Vec<f64>> {
    vec![vec![-1., 1.]; p]
}

#[test]
fn box_behnken_enumerates_every_two_factor_edge_and_requested_centers() {
    let limits = vec![vec![10., 20.], vec![-2., 2.], vec![100., 104.]];
    let j = run("box_behnken", &limits, &[2.]).unwrap();
    assert_eq!(values(&j, "run_count"), vec![14.]);
    let coded = values(&j, "coded_rows");
    let rows = values(&j, "rows");
    let mut edges = std::collections::BTreeSet::new();
    let mut centers = 0;
    for (z, x) in coded.chunks(3).zip(rows.chunks(3)) {
        near_vec(x, &[15. + 5. * z[0], 2. * z[1], 102. + 2. * z[2]]);
        if z == [0.; 3] {
            centers += 1;
        } else {
            assert_eq!(z.iter().filter(|&&v| v == 0.).count(), 1);
            assert!(z.iter().all(|v| [-1., 0., 1.].contains(v)));
            edges.insert(z.iter().map(|&v| v as i32).collect::<Vec<_>>());
        }
    }
    assert_eq!(edges.len(), 12);
    assert_eq!(centers, 2);
    let default = run("box_behnken", &bounds(4), &[]).unwrap();
    assert_eq!(values(&default, "run_count"), vec![25.]);
}
#[test]
fn latin_hypercube_stratifies_each_factor_and_uses_reproducible_independent_permutations() {
    let limits = vec![vec![10., 20.], vec![-2., 2.], vec![0., 1.]];
    let j = run("latin_hypercube", &limits, &[17., 42.]).unwrap();
    assert_eq!(j, run("latin_hypercube", &limits, &[17., 42.]).unwrap());
    assert_ne!(j, run("latin_hypercube", &limits, &[17., 43.]).unwrap());
    let flat = values(&j, "rows");
    assert_eq!(flat.len(), 51);
    for factor in 0..3 {
        let mut strata = flat
            .chunks(3)
            .map(|r| {
                assert!(r[factor] >= limits[factor][0] && r[factor] < limits[factor][1]);
                (((r[factor] - limits[factor][0]) / (limits[factor][1] - limits[factor][0])) * 17.)
                    .floor() as usize
            })
            .collect::<Vec<_>>();
        strata.sort_unstable();
        assert_eq!(strata, (0..17).collect::<Vec<_>>());
    }
}
#[test]
fn mixture_lattice_enumerates_simplex_compositions_and_filters_upper_constraints() {
    let j = run(
        "mixture_lattice",
        &[vec![0., 0.8], vec![0., 1.], vec![0., 1.]],
        &[2.],
    )
    .unwrap();
    let rows = values(&j, "rows");
    assert_eq!(
        rows,
        vec![0., 0., 1., 0., 0.5, 0.5, 0., 1., 0., 0.5, 0., 0.5, 0.5, 0.5, 0.]
    );
    assert_eq!(values(&j, "enumerated_compositions"), vec![6.]);
    assert_eq!(values(&j, "run_count"), vec![5.]);
    for r in rows.chunks(3) {
        near(r.iter().sum(), 1.);
        assert!(r[0] <= 0.8);
    }
}
#[test]
fn mixture_lattice_offsets_lower_bounds_before_distributing_remaining_total() {
    let limits = vec![vec![0.1, 1.], vec![0.2, 1.], vec![0.1, 1.]];
    let j = run("mixture_lattice", &limits, &[2., 1.]).unwrap();
    let rows = values(&j, "rows");
    assert_eq!(rows.len(), 18);
    for r in rows.chunks(3) {
        near(r.iter().sum(), 1.);
        for j in 0..3 {
            assert!(r[j] >= limits[j][0]);
            near(
                ((r[j] - limits[j][0]) / 0.3).round(),
                (r[j] - limits[j][0]) / 0.3,
            );
        }
    }
    let fixed = run(
        "mixture_lattice",
        &[vec![0.4, 0.4], vec![0.6, 0.6]],
        &[1000.],
    )
    .unwrap();
    near_vec(&values(&fixed, "rows"), &[0.4, 0.6]);
}
#[test]
fn fractional_factorial_generators_are_products_of_selected_base_columns() {
    let j = run("fractional_factorial", &bounds(4), &[3., 1., 7.]).unwrap();
    let coded = values(&j, "coded_rows");
    assert_eq!(values(&j, "run_count"), vec![9.]);
    for r in coded.chunks(4).take(8) {
        near(r[3], r[0] * r[1] * r[2]);
        assert!(r.iter().all(|v| v.abs() == 1.));
    }
    assert_eq!(&coded[32..], &[0.; 4]);
    for j in 0..4 {
        near(coded.chunks(4).take(8).map(|r| r[j]).sum(), 0.);
        for k in 0..j {
            near(coded.chunks(4).take(8).map(|r| r[j] * r[k]).sum(), 0.);
        }
    }
}
#[test]
fn orthogonal_screening_uses_nonconstant_sylvester_hadamard_columns() {
    let j = run("orthogonal_screening", &bounds(5), &[]).unwrap();
    let coded = values(&j, "coded_rows");
    assert_eq!(values(&j, "run_count"), vec![8.]);
    assert!(j.contains("sylvester_hadamard_main_effect_screening"));
    for a in 0..5 {
        near(coded.chunks(5).map(|r| r[a]).sum(), 0.);
        for b in 0..5 {
            near(
                coded.chunks(5).map(|r| r[a] * r[b]).sum(),
                if a == b { 8. } else { 0. },
            );
        }
    }
}
#[test]
fn randomized_replicates_are_lossless_and_each_replicate_has_balanced_blocks() {
    let data = (0..6)
        .map(|i| vec![i as f64, 100. + i as f64])
        .collect::<Vec<_>>();
    let j = run("randomize_design", &data, &[3., 42., 2.]).unwrap();
    assert_eq!(j, run("randomize_design", &data, &[3., 42., 2.]).unwrap());
    assert_ne!(j, run("randomize_design", &data, &[3., 43., 2.]).unwrap());
    let sources = values(&j, "source_indices");
    let replicate = values(&j, "replicate_labels");
    let blocks = values(&j, "block_labels");
    let rows = values(&j, "rows");
    assert_eq!(sources.len(), 18);
    for r in 0..3 {
        let mut ids = sources[6 * r..6 * r + 6]
            .iter()
            .map(|&v| v as usize)
            .collect::<Vec<_>>();
        ids.sort_unstable();
        assert_eq!(ids, (0..6).collect::<Vec<_>>());
        assert_eq!(&replicate[6 * r..6 * r + 6], &[r as f64; 6]);
        assert_eq!(&blocks[6 * r..6 * r + 6], &[0., 0., 0., 1., 1., 1.]);
    }
    for (i, r) in rows.chunks(2).enumerate() {
        assert_eq!(r, data[sources[i] as usize]);
    }
}
fn evaluation_rows() -> Vec<Vec<f64>> {
    vec![
        vec![-1., 0.],
        vec![0., 1.],
        vec![1., 0.],
        vec![2., 1.],
        vec![-2., -1.],
    ]
}
#[test]
fn design_evaluation_matches_independent_numpy_information_inverse_singular_values_and_leverage() {
    let j = run("design_evaluate", &evaluation_rows(), &[]).unwrap();
    assert_eq!(values(&j, "rank"), vec![3.]);
    assert_eq!(
        values(&j, "information_matrix"),
        vec![5., 0., 1., 0., 10., 4., 1., 4., 3.]
    );
    near_vec(
        &values(&j, "inverse_information_matrix"),
        &[
            7. / 30.,
            1. / 15.,
            -1. / 6.,
            1. / 15.,
            7. / 30.,
            -1. / 3.,
            -1. / 6.,
            -1. / 3.,
            5. / 6.,
        ],
    );
    near(values(&j, "log_determinant")[0], 4.0943445622221);
    near(values(&j, "coded_log_determinant")[0], 3.9810158769150976);
    near(values(&j, "condition_number_coded")[0], 2.6822257700311996);
    near_vec(
        &values(&j, "leverages"),
        &[1. / 3., 11. / 15., 0.6, 0.6, 11. / 15.],
    );
    near(values(&j, "leverages").iter().sum(), 3.);
}
#[test]
fn design_evaluation_retains_identifiability_and_leverage_under_feature_unit_changes() {
    let baseline = run("design_evaluate", &evaluation_rows(), &[]).unwrap();
    let transformed = evaluation_rows()
        .iter()
        .map(|r| vec![100. + 3. * r[0], -20. + 2. * r[1]])
        .collect::<Vec<_>>();
    let j = run("design_evaluate", &transformed, &[]).unwrap();
    near_vec(&values(&j, "leverages"), &values(&baseline, "leverages"));
    near(
        values(&j, "coded_log_determinant")[0],
        values(&baseline, "coded_log_determinant")[0],
    );
    near(
        values(&j, "log_determinant")[0],
        values(&baseline, "log_determinant")[0] + 2. * 6_f64.ln(),
    );
    near(
        values(&j, "condition_number_coded")[0],
        values(&baseline, "condition_number_coded")[0],
    );
}
#[test]
fn advanced_design_rejects_malformed_generators_infeasible_mixtures_and_work_explosions() {
    for op in [
        "box_behnken",
        "latin_hypercube",
        "mixture_lattice",
        "fractional_factorial",
        "orthogonal_screening",
        "randomize_design",
        "design_evaluate",
    ] {
        for rows in [
            vec![],
            vec![vec![]],
            vec![vec![f64::NAN]],
            vec![vec![1., 2.], vec![3.]],
        ] {
            assert!(run(op, &rows, &[]).is_err(), "{op}");
        }
        assert!(run(op, &bounds(3), &[f64::INFINITY]).is_err(), "{op}");
    }
    for (op, rows, args) in [
        ("box_behnken", bounds(2), vec![]),
        ("box_behnken", bounds(3), vec![-1.]),
        ("latin_hypercube", bounds(3), vec![100_001., 0.]),
        ("latin_hypercube", bounds(3), vec![2., 0.5]),
        (
            "mixture_lattice",
            vec![vec![0.6, 1.], vec![0.6, 1.]],
            vec![2.],
        ),
        (
            "mixture_lattice",
            vec![vec![0., 0.1], vec![0., 0.1]],
            vec![2.],
        ),
        ("mixture_lattice", vec![vec![0., 1.]; 64], vec![1000.]),
        ("fractional_factorial", bounds(4), vec![3., 0., 1.]),
        ("fractional_factorial", bounds(4), vec![3., 0., 8.]),
        ("fractional_factorial", bounds(5), vec![3., 0., 3., 3.]),
        ("fractional_factorial", bounds(4), vec![3., 0.]),
        ("orthogonal_screening", bounds(3), vec![1.]),
        ("randomize_design", vec![vec![1.]; 5], vec![1., 0., 2.]),
        ("randomize_design", vec![vec![1.]; 1000], vec![101., 0.]),
        (
            "design_evaluate",
            vec![vec![0., 0.], vec![1., 1.], vec![2., 2.]],
            vec![],
        ),
        ("design_evaluate", vec![vec![1.], vec![1.]], vec![]),
        (
            "design_evaluate",
            vec![vec![1e308], vec![-1e308], vec![0.]],
            vec![],
        ),
    ] {
        assert!(run(op, &rows, &args).is_err(), "{op}: {args:?}");
    }
}

#[test]
fn mixture_sum_tolerance_respects_small_measurement_units() {
    // An absolute tolerance based on max(total,1) accepts this impossible
    // mixture, whose fixed lower-bound sum is 85 orders above its total.
    assert!(run(
        "mixture_lattice",
        &[vec![1e-15, 1e-15], vec![1e-15, 1e-15]],
        &[2., 1e-100]
    )
    .is_err());
    let total = 1e-100;
    let j = run("mixture_lattice", &vec![vec![0., total]; 3], &[2., total]).unwrap();
    let normalized = values(&j, "rows")
        .iter()
        .map(|v| v / total)
        .collect::<Vec<_>>();
    near_vec(
        &normalized,
        &[
            0., 0., 1., 0., 0.5, 0.5, 0., 1., 0., 0.5, 0., 0.5, 0.5, 0.5, 0., 1., 0., 0.,
        ],
    );
    assert!(values(&j, "sum_tolerance")[0] < 1e-110);
}
