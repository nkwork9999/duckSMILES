use ducksmiles_smiles::{descriptor_profile::SPECS, workbench::analyze};

#[test]
fn non_tetrahedral_original_inputs_match_all_registered_reference_values() {
    let mut differences = Vec::new();
    let mut count = 0;
    for line in include_str!("fixtures/descriptor_non_tetra.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        assert_eq!(f.len(), SPECS.len() + 2);
        count += 1;
        let Some(actual) = analyze(f[0]) else {
            differences.push(format!("{}: unsupported", f[0]));
            continue;
        };
        if actual.formula != f[1] {
            differences.push(format!("{} formula: {} vs {}", f[0], actual.formula, f[1]));
        }
        for ((spec, actual), expected) in SPECS.iter().zip(actual.descriptor_values).zip(&f[2..]) {
            let expected = if *expected == "null" {
                None
            } else {
                Some(expected.parse::<f64>().unwrap())
            };
            let same = match (actual, expected) {
                (None, None) => true,
                (Some(a), Some(e)) => {
                    (a - e).abs() <= spec.tolerance.max(spec.relative_tolerance * e.abs())
                }
                _ => false,
            };
            if !same {
                differences.push(format!(
                    "{} {}: {:?} vs {:?}",
                    f[0], spec.id, actual, expected
                ));
            }
        }
    }
    assert_eq!(count, 132);
    assert!(
        differences.is_empty(),
        "{} differences; first25:\n{}",
        differences.len(),
        differences
            .iter()
            .take(25)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn invalid_non_tetrahedral_permutations_are_not_accepted_as_prefixes() {
    let mut count = 0;
    for source in include_str!("fixtures/non_tetra_rejected.txt").lines() {
        assert!(analyze(source).is_none(), "pinned parser rejects {source}");
        count += 1;
    }
    assert_eq!(count, 6);
}

#[test]
fn extended_prepared_geometry_does_not_change_the_raw_graph_grammar() {
    for source in [
        "F[Pt@SP1](Cl)(Br)I",
        "F[Fe@TB3](Cl)(Br)(I)N",
        "F[Co@OH1](Cl)(Br)(I)(N)O",
    ] {
        // SAFETY: readable source, no output allocation or dereference.
        let n = unsafe {
            ducksmiles_smiles::ds_mol_graph_json(
                source.as_ptr(),
                source.len(),
                std::ptr::null_mut(),
                0,
            )
        };
        assert!(
            n <= 0,
            "raw graph must not advertise geometry it does not model: {source}"
        );
    }
}
