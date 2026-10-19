use ducksmiles_smiles::workbench::analyze;

#[test]
fn coordinate_free_cx_wedges_match_original_annotated_references() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut differences = Vec::new();
    let mut count = 0;
    for line in include_str!("fixtures/descriptor_cx_wedge.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), SPECS.len() + 2);
        count += 1;
        let Some(actual) = analyze(fields[0]) else {
            differences.push(format!("{}: unsupported", fields[0]));
            continue;
        };
        if actual.formula != fields[1] {
            differences.push(format!(
                "{}: formula {} vs {}",
                fields[0], actual.formula, fields[1]
            ));
        }
        for ((spec, actual), expected) in
            SPECS.iter().zip(actual.descriptor_values).zip(&fields[2..])
        {
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
                    fields[0], spec.id, actual, expected
                ));
            }
        }
    }
    assert_eq!(
        count, 482,
        "original population plus independent synthetic controls"
    );
    assert!(
        differences.is_empty(),
        "{} differences, first20:\n{}",
        differences.len(),
        differences
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn cx_wedges_reject_malformed_or_chemically_unhandled_extensions() {
    for input in [
        "CC |wU:0.0",
        "CC|wU:0.0|",
        "CC |wU:0.0| trailing",
        "CC |wU:0.0| |wD:0.0|",
        "CC ||",
        "CC |wU:|",
        "CC |wX:0.0|",
        "CC |wU:0.0,|",
        "CC |wU:0.0,wD:1.0|",
        "CC |wU:0.1|",
        "CC |wU:9.0|",
        "CCC |wU:0.1|",
        "CC |wU:-1.0|",
        "CC |wU:0.18446744073709551616|",
        "CC |wU:0.0,ctu:0|",
        "CC |w:0.0|",
        "CC |(0,0,;1,0,),wU:0.0|",
        "CC |^1:0,wU:0.0|",
        "CO |C:0.0,wU:0.0|",
        "CC |wU:0.0,$foo;bar$|",
        "CC |wU:٠.0|",
        "N->[Pt] |wU:1.0|",
    ] {
        assert!(
            analyze(input).is_none(),
            "must not silently discard unsupported meaning: {input}"
        );
    }
    let oversized = format!("CC {}|wU:0.0|", " ".repeat(4096));
    assert!(analyze(&oversized).is_none());
}

#[test]
fn registered_information_descriptors_match_exact_graph_definitions() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let ipc = SPECS
        .iter()
        .position(|s| s.id == "Ipc")
        .expect("Ipc is missing");
    let average = SPECS
        .iter()
        .position(|s| s.id == "AvgIpc")
        .expect("AvgIpc is missing");
    let mut count = 0;
    for line in include_str!("fixtures/descriptor_ipc_smiles.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let actual = analyze(fields[0]).expect(fields[0]);
        assert_eq!(actual.formula, fields[1], "{} formula", fields[0]);
        for (index, column) in [(ipc, 2), (average, 3)] {
            let expected: f64 = fields[column].parse().unwrap();
            let value = actual.descriptor_values[index].expect("information value must be defined");
            assert!(
                (value - expected).abs() <= (1e-12 * expected.abs()).max(1e-12),
                "{} {} {value} vs {expected}",
                fields[0],
                SPECS[index].id
            );
        }
        count += 1;
    }
    assert_eq!(count, 34);
}

#[test]
fn attached_hydrogens_on_hydrogen_vertices_are_included_in_profiles() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let ids = [
        "ds_ExplicitHydrogens",
        "ds_ImplicitHydrogens",
        "ds_TotalHydrogens",
    ];
    let indices = ids.map(|id| SPECS.iter().position(|s| s.id == id).unwrap());
    let mut differences = Vec::new();
    let mut count = 0;
    for line in include_str!("fixtures/descriptor_hydrogen_profile.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let actual = analyze(f[0]).expect(f[0]);
        assert_eq!(actual.formula, f[1], "{} formula", f[0]);
        for (column, (&index, id)) in indices.iter().zip(ids).enumerate() {
            let expected = f[column + 2].parse::<f64>().unwrap();
            if actual.descriptor_values[index] != Some(expected) {
                differences.push(format!(
                    "{} {id}: {:?} vs {expected}",
                    f[0], actual.descriptor_values[index]
                ));
            }
        }
        count += 1;
    }
    assert_eq!(count, 10);
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}

#[test]
fn charged_oxygen_and_sulfur_acceptors_match_pinned_smarts() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let index = SPECS.iter().position(|s| s.id == "NumHAcceptors").unwrap();
    let mut differences = Vec::new();
    let mut count = 0;
    for line in include_str!("fixtures/descriptor_hba_charge.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let expected: usize = f[1].parse().unwrap();
        let actual = analyze(f[0]).expect(f[0]);
        if actual.hba != expected || actual.descriptor_values[index] != Some(expected as f64) {
            differences.push(format!(
                "{}: {} / {:?} vs {expected}",
                f[0], actual.hba, actual.descriptor_values[index]
            ));
        }
        count += 1;
    }
    assert_eq!(count, 36);
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}

#[test]
fn registered_sps_matches_independent_modern_stereo_pipeline() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let index = SPECS
        .iter()
        .position(|s| s.id == "SPS")
        .expect("SPS standard descriptor missing");
    let mut differences = Vec::new();
    let mut sources = 0;
    let mut undefined = 0;
    for line in include_str!("fixtures/descriptor_sps.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let actual = analyze(f[0]).expect(f[0]).descriptor_values[index];
        let expected = if f[1] == "null" {
            undefined += 1;
            None
        } else {
            Some(f[1].parse::<f64>().unwrap())
        };
        let same = match (actual, expected) {
            (Some(a), Some(e)) => (a - e).abs() <= 1e-9,
            (None, None) => true,
            _ => false,
        };
        if !same {
            differences.push(format!("{}: {:?} vs {:?}", f[0], actual, expected));
        }
        sources += 1;
    }
    assert!(
        sources >= 650 && undefined >= 3,
        "reference coverage unexpectedly reduced"
    );
    assert!(
        differences.is_empty(),
        "{} SPS differences, first30:\n{}",
        differences.len(),
        differences
            .iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn rejected_carbon_stereo_tags_reclassify_hydrogens_without_losing_radicals() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_stereo_h.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let actual = analyze(f[0]).expect(f[0]);
        if f[1] == "@formula" {
            if actual.formula != f[2] {
                differences.push(format!("{} formula {} vs {}", f[0], actual.formula, f[2]));
            }
        } else {
            let index = SPECS.iter().position(|s| s.id == f[1]).unwrap();
            let expected: f64 = f[2].parse().unwrap();
            if !actual.descriptor_values[index]
                .is_some_and(|v| (v - expected).abs() <= SPECS[index].tolerance)
            {
                differences.push(format!(
                    "{} {} {:?} vs {expected}",
                    f[0], f[1], actual.descriptor_values[index]
                ));
            }
        }
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}

#[test]
fn registered_legacy_stereocounts_match_original_input_reference() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    for id in ["NumAtomStereoCenters", "NumUnspecifiedAtomStereoCenters"] {
        assert!(
            SPECS.iter().any(|s| s.id == id),
            "standard descriptor missing: {id}"
        );
    }
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    let mut sources = 0;
    let mut cells = 0;
    for line in include_str!("fixtures/descriptor_stereo_counts.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        if f[0] != previous {
            actual = Some(analyze(f[0]).expect(f[0]));
            previous = f[0];
            sources += 1;
        }
        let index = SPECS.iter().position(|s| s.id == f[1]).unwrap();
        let expected: f64 = f[2].parse().unwrap();
        let observed = actual.as_ref().unwrap().descriptor_values[index];
        if observed != Some(expected) {
            differences.push(format!("{} {}: {observed:?} vs {expected}", f[0], f[1]));
        }
        cells += 1;
    }
    assert_eq!((sources, cells), (658, 1316));
    assert!(
        differences.is_empty(),
        "{} stereo differences, first30:\n{}",
        differences.len(),
        differences
            .iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn additional_fifty_descriptors_match_independent_dense_graph_reference() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let lines: Vec<_> = include_str!("fixtures/descriptor_extra50.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .collect();
    let ids: std::collections::BTreeSet<_> = lines.iter().map(|f| f[1]).collect();
    assert_eq!(ids.len(), 50);
    assert_eq!(lines.len(), 24_000);
    assert_eq!(
        lines
            .iter()
            .map(|f| f[0])
            .collect::<std::collections::BTreeSet<_>>()
            .len(),
        480
    );
    let positive: std::collections::BTreeSet<_> = lines
        .iter()
        .filter(|f| f[2].parse::<f64>().is_ok_and(|v| v > 0.))
        .map(|f| f[1])
        .collect();
    assert_eq!(
        positive, ids,
        "every added descriptor needs a positive independent control"
    );
    let missing: Vec<_> = ids
        .iter()
        .filter(|&&id| !SPECS.iter().any(|s| s.id == id))
        .collect();
    assert!(
        missing.is_empty(),
        "additional descriptors not registered: {missing:?}"
    );
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for fields in &lines {
        if fields[0] != previous {
            actual = Some(analyze(fields[0]).expect(fields[0]));
            previous = fields[0];
        }
        let index = SPECS.iter().position(|s| s.id == fields[1]).unwrap();
        let observed = actual.as_ref().unwrap().descriptor_values[index];
        let matches = if fields[2] == "null" {
            observed.is_none()
        } else {
            let expected: f64 = fields[2].parse().unwrap();
            observed.is_some_and(|v| (v - expected).abs() <= SPECS[index].tolerance)
        };
        if !matches {
            differences.push(format!(
                "{} {}: {observed:?} vs {}",
                fields[0], fields[1], fields[2]
            ));
        }
    }
    assert!(
        differences.is_empty(),
        "{} extra50 differences, first20:\n{}",
        differences.len(),
        differences
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn partial_aromatic_bond_hydrogens_and_retained_formula_residuals_match_reference() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut sources = 0;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_valence_h.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        if fields[0] != previous {
            actual = Some(analyze(fields[0]).expect(fields[0]));
            previous = fields[0];
            sources += 1;
        }
        let row = actual.as_ref().unwrap();
        if fields[1] == "@formula" {
            if row.formula != fields[2] {
                differences.push(format!(
                    "{} formula: {} vs {}",
                    fields[0], row.formula, fields[2]
                ));
            }
        } else {
            let index = SPECS.iter().position(|s| s.id == fields[1]).unwrap();
            let expected: f64 = fields[2].parse().unwrap();
            let matches = if expected.is_finite() {
                row.descriptor_values[index]
                    .is_some_and(|v| (v - expected).abs() <= SPECS[index].tolerance)
            } else {
                row.descriptor_values[index].is_none()
            };
            if !matches {
                differences.push(format!(
                    "{} {}: {:?} vs {expected}",
                    fields[0], fields[1], row.descriptor_values[index]
                ));
            }
        }
    }
    assert_eq!(sources, 41);
    assert!(
        differences.is_empty(),
        "{} differences, first20:\n{}",
        differences.len(),
        differences
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn registered_qed_matches_independent_scores_with_original_inputs() {
    let index = ducksmiles_smiles::descriptor_profile::SPECS
        .iter()
        .position(|s| s.id == "qed")
        .expect("QED must be registered in the workbench catalog");
    let mut differences = Vec::new();
    let mut checked = 0;
    for line in include_str!("fixtures/descriptor_qed.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let expected: f64 = fields[1].parse().unwrap();
        let row = analyze(fields[0]).expect(fields[0]);
        if !row.descriptor_values[index].is_some_and(|v| (v - expected).abs() <= 1e-9) {
            differences.push(format!(
                "{}: {:?} vs {expected}",
                fields[0], row.descriptor_values[index]
            ));
        }
        checked += 1;
    }
    assert_eq!(checked, 413);
    assert!(
        differences.is_empty(),
        "{} QED differences:\n{}",
        differences.len(),
        differences
            .iter()
            .take(30)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn aromatic_radical_and_default_valence_preparation_matches_all_catalog_controls() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut errors = Vec::new();
    for line in include_str!("fixtures/descriptor_aromatic_candidate.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        if previous != f[0] {
            actual = Some(analyze(f[0]).expect(f[0]));
            previous = f[0];
        }
        let row = actual.as_ref().unwrap();
        if f[1] == "@formula" {
            if row.formula != f[2] {
                errors.push(format!("{} formula: {} vs {}", f[0], row.formula, f[2]));
            }
        } else {
            let i = SPECS.iter().position(|s| s.id == f[1]).expect(f[1]);
            let expected: f64 = f[2].parse().unwrap();
            let value = row.descriptor_values[i];
            let matches = if expected.is_finite() {
                value.is_some_and(|v| (v - expected).abs() <= SPECS[i].tolerance)
            } else {
                value.is_none()
            };
            if !matches {
                errors.push(format!("{} {}: {value:?} vs {expected}", f[0], f[1]));
            }
        }
    }
    assert!(
        errors.is_empty(),
        "{} differences:\n{}",
        errors.len(),
        errors.into_iter().take(20).collect::<Vec<_>>().join("\n")
    );
}

#[test]
fn registered_morgan_density_uses_unfolded_sparse_keys_and_heavy_denominator() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let indices: Vec<_> = ["FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3"]
        .iter()
        .map(|id| SPECS.iter().position(|s| s.id == *id).expect(id))
        .collect();
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_morgan_density.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let row = analyze(fields[0]).expect(fields[0]);
        for (radius, &index) in indices.iter().enumerate() {
            let expected: f64 = fields[radius + 1].parse().unwrap();
            let actual = row.descriptor_values[index];
            if !actual.is_some_and(|v| (v - expected).abs() <= SPECS[index].tolerance) {
                differences.push(format!(
                    "{} radius{}: {actual:?} vs {expected}",
                    fields[0],
                    radius + 1
                ));
            }
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences
            .into_iter()
            .take(25)
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn original_and_permuted_clusters_match_all_descriptors_and_formula() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut failures = Vec::new();
    for line in include_str!("fixtures/descriptor_clusters.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        assert_eq!(fields.len(), SPECS.len() + 2);
        let Some(row) = analyze(fields[0]) else {
            failures.push(format!("{} unsupported", fields[0]));
            continue;
        };
        if row.formula != fields[1] {
            failures.push(format!(
                "{} formula {} != {}",
                fields[0], row.formula, fields[1]
            ));
        }
        for (i, spec) in SPECS.iter().enumerate() {
            let expected = fields[i + 2];
            let actual = row.descriptor_values[i];
            let good = if expected == "null" {
                actual.is_none()
            } else {
                let expected: f64 = expected.parse().unwrap();
                actual.is_some_and(|v| {
                    (v - expected).abs()
                        <= spec.tolerance.max(spec.relative_tolerance * expected.abs())
                })
            };
            if !good {
                failures.push(format!(
                    "{} {}: {actual:?} != {expected}",
                    fields[0], spec.id
                ));
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{} differences:\n{}",
        failures.len(),
        failures
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn directed_coordination_matches_independent_all_catalog_and_formula_controls() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_coordination.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        if previous != fields[0] {
            actual = Some(analyze(fields[0]).expect(fields[0]));
            previous = fields[0];
        }
        let row = actual.as_ref().unwrap();
        if fields[1] == "@formula" {
            if row.formula != fields[2] {
                differences.push(format!(
                    "{} formula: {} vs {}",
                    fields[0], row.formula, fields[2]
                ));
            }
            continue;
        }
        let i = SPECS
            .iter()
            .position(|s| s.id == fields[1])
            .expect(fields[1]);
        let expected: f64 = fields[2].parse().unwrap();
        let value = row.descriptor_values[i];
        let matches = if expected.is_finite() {
            value.is_some_and(|x| (x - expected).abs() <= SPECS[i].tolerance)
        } else {
            value.is_none()
        };
        if !matches {
            differences.push(format!(
                "{} {}: {value:?} vs {expected}",
                fields[0], fields[1]
            ));
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences
            .iter()
            .take(50)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn strict_rotors_keep_explicit_vertex_degree_and_total_h_query_semantics() {
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_rotors.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let expected: usize = fields[1].parse().unwrap();
        let actual = analyze(fields[0]).expect(fields[0]).rotatable_bonds;
        if actual != expected {
            differences.push(format!("{}: {actual} vs {expected}", fields[0]));
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn native_workbench_cleanup_matches_pinned_all_catalog_and_formula_controls() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    // Catches omitted/wrong N/P/halogen charge and bond-order cleanup,
    // including explicit H, isotope position and aromaticity after cleanup.
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_cleanup.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        if previous != fields[0] {
            actual = Some(analyze(fields[0]).expect(fields[0]));
            previous = fields[0];
        }
        let row = actual.as_ref().unwrap();
        if fields[1] == "@formula" {
            if row.formula != fields[2] {
                differences.push(format!(
                    "{} formula: {} vs {}",
                    fields[0], row.formula, fields[2]
                ));
            }
            continue;
        }
        let i = SPECS
            .iter()
            .position(|s| s.id == fields[1])
            .expect(fields[1]);
        let expected: f64 = fields[2].parse().unwrap();
        let value = row.descriptor_values[i];
        let matches = if expected.is_finite() {
            value.is_some_and(|x| (x - expected).abs() <= SPECS[i].tolerance)
        } else {
            value.is_none()
        };
        if !matches {
            differences.push(format!(
                "{} {}: {value:?} vs {expected}",
                fields[0], fields[1]
            ));
        }
    }
    assert!(
        differences.is_empty(),
        "{} mismatches:\n{}",
        differences.len(),
        differences
            .iter()
            .take(40)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

#[test]
fn registered_structural_counts_match_independent_source_references() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    // Catches wrong amide O-ring/bond rules and pairwise-ring intersection
    // counting (including deduplication), using independent registered outputs.
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_structural.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        if previous != fields[0] {
            actual = Some(analyze(fields[0]).expect(fields[0]));
            previous = fields[0];
        }
        let index = SPECS
            .iter()
            .position(|s| s.id == fields[1])
            .unwrap_or_else(|| panic!("missing registered structural descriptor {}", fields[1]));
        let expected: f64 = fields[2].parse().unwrap();
        let value = actual.as_ref().unwrap().descriptor_values[index];
        if value != Some(expected) {
            differences.push(format!(
                "{} {}: {value:?}, expected {expected}",
                fields[0], fields[1]
            ));
        }
    }
    assert!(
        differences.is_empty(),
        "{} mismatches:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn ethanol_has_physically_expected_descriptors() {
    let value = analyze("CCO").expect("valid ethanol");
    assert!((value.mw - 46.069).abs() < 0.01);
    assert!((value.exact_mass - 46.0418648).abs() < 0.00001);
    assert_eq!(value.formula, "C2H6O");
    assert_eq!(value.heavy_atoms, 3);
    assert_eq!(value.bonds, 2);
    assert_eq!(value.hbd, 1);
    assert_eq!(value.hba, 1);
    assert!((value.tpsa - 20.23).abs() < 0.01);
    assert_eq!(value.fragments, 1);
}

#[test]
fn salt_keeps_charge_and_disconnected_components() {
    let value = analyze("[NH4+].[Cl-]").unwrap();
    assert_eq!(value.fragments, 2);
    assert_eq!(value.formal_charge, 0);
    assert_eq!(value.formula, "H4ClN");
}

#[test]
fn malformed_and_unbounded_inputs_do_not_produce_results() {
    for input in ["", "C(", "not-a-molecule", "C=", "C..O"] {
        assert!(analyze(input).is_none(), "{input}");
    }
    assert!(analyze(&"C".repeat(641)).is_none());
}

#[test]
fn native_batch_preserves_source_positions_and_obeys_buffer_contract() {
    use ducksmiles_smiles::workbench::ds_workbench_analyze;
    let input = b"CCO\ninvalid\n[Na+].[Cl-]";
    let required =
        unsafe { ds_workbench_analyze(input.as_ptr(), input.len(), std::ptr::null_mut(), 0) };
    assert!(required > 0, "query required size");
    let mut tiny = [0x55u8; 2];
    let result = unsafe { ds_workbench_analyze(input.as_ptr(), input.len(), tiny.as_mut_ptr(), 2) };
    assert_eq!(result, required);
    assert_eq!(tiny, [0x55; 2], "undersized destination stays untouched");
    let mut output = vec![0; required as usize];
    assert_eq!(
        unsafe {
            ds_workbench_analyze(
                input.as_ptr(),
                input.len(),
                output.as_mut_ptr(),
                output.len(),
            )
        },
        required
    );
    let text = String::from_utf8(output).unwrap();
    assert!(text.contains("\"index\":1,\"status\":\"invalid_or_unsupported\""));
    assert!(text.contains("\"index\":2,\"status\":\"ok\""));
    assert_eq!(
        unsafe { ds_workbench_analyze([0xff].as_ptr(), 1, std::ptr::null_mut(), 0) },
        -1
    );
}

#[test]
fn mass_uses_isotopes_charge_and_complete_element_data() {
    // RDKit 2026.03.6 literals, independent of the native mass helpers.
    for (smiles, mw, exact) in [
        ("CCS", 62.137, 62.019021192),
        ("[13CH4]", 17.03535484, 17.034654968),
        ("[NH4+]", 18.039, 18.03382554809),
        ("[Cl-]", 35.453, 34.96940125991),
        ("[Hg]", 200.59, 201.970643),
        ("[Th]", 232.038, 232.0380553),
    ] {
        let d = analyze(smiles).expect(smiles);
        assert!((d.mw - mw).abs() < 1e-6, "MW {smiles}: {}", d.mw);
        assert!(
            (d.exact_mass - exact).abs() < 1e-7,
            "exact {smiles}: {}",
            d.exact_mass
        );
    }
    assert!(
        analyze("[999C]").is_none(),
        "unknown isotope must not silently use carbon-12"
    );
}

#[test]
fn cache_revision_obeys_sized_output_contract() {
    use ducksmiles_smiles::workbench::ds_workbench_revision;
    unsafe {
        let required = ds_workbench_revision(std::ptr::null_mut(), 0);
        assert!(required > 0);
        let mut short = [0x55u8; 2];
        assert_eq!(
            ds_workbench_revision(short.as_mut_ptr(), short.len()),
            required
        );
        assert_eq!(short, [0x55; 2]);
        let mut out = vec![0u8; required as usize];
        assert_eq!(ds_workbench_revision(out.as_mut_ptr(), out.len()), required);
        assert!(std::str::from_utf8(&out).is_ok());
        assert_eq!(ds_workbench_revision(std::ptr::null_mut(), 4), -1);
    }
}

#[test]
fn phosphorus_implicit_hydrogen_and_carbon_free_formula_match_rdkit() {
    let p = analyze("CCOP(=O)OCC").unwrap();
    assert_eq!(p.formula, "C4H11O3P");
    assert!((p.mw - 138.103).abs() < 1e-6);
    assert!((p.exact_mass - 138.044580842).abs() < 1e-7);
    assert_eq!(analyze("OP(O)(F)=O").unwrap().formula, "H2FO3P");
    assert_eq!(analyze("[NH4+].[Cl-]").unwrap().formula, "H4ClN");
}

#[test]
fn extended_batch_preserves_legacy_values_and_exposes_descriptor_vector() {
    let input = b"CCO\ninvalid";
    let mut output = vec![0u8; 65536];
    let n = unsafe {
        ducksmiles_smiles::workbench::ds_workbench_analyze(
            input.as_ptr(),
            input.len(),
            output.as_mut_ptr(),
            output.len(),
        )
    };
    assert!(n > 0 && (n as usize) <= output.len());
    let json = std::str::from_utf8(&output[..n as usize]).unwrap();
    assert!(
        json.contains("\"descriptor_values\":["),
        "extended descriptors missing from valid row: {json}"
    );
    assert!(json.contains("\"mw\":46.069"));
    assert!(json.contains("\"index\":1,\"status\":\"invalid_or_unsupported\""));
}

#[test]
fn catalog_values_match_pinned_independent_rdkit_literals() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let fixture = include_str!("fixtures/descriptor_phase1.tsv");
    let mut previous = String::new();
    let mut actual = None;
    let mut seen = std::collections::HashSet::new();
    let mut differences = Vec::new();
    for spec in SPECS {
        assert!(seen.insert(spec.id), "duplicate descriptor ID");
    }
    for line in fixture.lines().filter(|line| !line.starts_with('#')) {
        let mut fields = line.split('\t');
        let smiles = fields.next().unwrap();
        let id = fields.next().unwrap();
        let expected: f64 = fields.next().unwrap().parse().unwrap();
        if previous != smiles {
            actual = Some(analyze(smiles).expect(smiles));
            previous = smiles.to_owned();
        }
        let index = SPECS.iter().position(|s| s.id == id).expect(id);
        let row = actual.as_ref().unwrap();
        assert_eq!(row.descriptor_values.len(), SPECS.len());
        let value = row.descriptor_values[index].expect(id);
        if (value - expected).abs() > SPECS[index].tolerance {
            differences.push(format!("{smiles} {id}: native {value}, RDKit {expected}"));
        }
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}

#[test]
fn explicit_ring_bonds_are_accepted_without_moving_the_bond() {
    for (smiles, formula) in [
        ("C=1CCCCC1", "C6H10"),
        ("C1CCCCC=1", "C6H10"),
        ("C1=CCCCC=1", "C6H8"),
        ("C=1CCCCC=1", "C6H10"),
        ("C%10CCCCC=%10", "C6H10"),
        ("C1:C:C:C:C:C:1", "C6H6"),
        ("C1:C:N:S:N:1", "C2H2N2S"),
    ] {
        let value = analyze(smiles).unwrap_or_else(|| panic!("valid explicit ring bond: {smiles}"));
        assert_eq!(value.formula, formula, "{smiles}");
    }
    let input = b"C=1CNCCC1";
    let mut out = vec![0u8; 65536];
    let count = unsafe {
        ducksmiles_smiles::ds_mol_graph_json(
            input.as_ptr(),
            input.len(),
            out.as_mut_ptr(),
            out.len(),
        )
    };
    assert!(count > 0 && (count as usize) <= out.len());
    let graph = std::str::from_utf8(&out[..count as usize]).unwrap();
    assert!(
        graph.contains("\"index\":0,\"begin\":0,\"end\":1,\"bond_type\":\"single\""),
        "opening annotation must not leak onto the following bond"
    );
    assert!(
        graph.contains("\"index\":5,\"begin\":0,\"end\":5,\"bond_type\":\"double\""),
        "opening annotation belongs to ring closure"
    );
}

#[test]
fn registered_functional_groups_match_independent_rdkit_counts() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    let mut positive_ids = std::collections::HashSet::new();
    for line in include_str!("fixtures/descriptor_fragments.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let (smiles, id) = (fields[0], fields[1]);
        let expected: f64 = fields[2].parse().unwrap();
        if expected > 0.0 {
            positive_ids.insert(id);
        }
        if previous != smiles {
            actual = Some(analyze(smiles).unwrap_or_else(|| panic!("fixture rejected: {smiles}")));
            previous = smiles;
        }
        let index = SPECS
            .iter()
            .position(|s| s.id == id)
            .unwrap_or_else(|| panic!("missing registered descriptor: {id}"));
        let value = actual.as_ref().unwrap().descriptor_values[index];
        if value != Some(expected) {
            differences.push(format!("{smiles} {id}: native {value:?}, RDKit {expected}"));
        }
    }
    assert_eq!(
        positive_ids.len(),
        85,
        "every registered fragment must have a nonzero independent example"
    );
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn registered_surface_descriptors_match_independent_rdkit_values() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_surface.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let fields: Vec<_> = line.split('\t').collect();
        let (smiles, id) = (fields[0], fields[1]);
        let expected: f64 = fields[2].parse().unwrap();
        if previous != smiles {
            actual = Some(analyze(smiles).unwrap_or_else(|| panic!("fixture rejected: {smiles}")));
            previous = smiles;
        }
        let index = SPECS
            .iter()
            .position(|s| s.id == id)
            .unwrap_or_else(|| panic!("missing registered descriptor: {id}"));
        let value = actual.as_ref().unwrap().descriptor_values[index];
        let agrees = if expected.is_finite() {
            value.is_some_and(|v| (v - expected).abs() <= 1e-9)
        } else {
            value.is_none()
        };
        if !agrees {
            differences.push(format!("{smiles} {id}: native {value:?}, RDKit {expected}"));
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn partial_explicit_hydrogens_keep_expanded_totals_and_original_atom_surfaces() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    for smiles in ["[2H]O", "[H]O"] {
        let row = analyze(smiles).unwrap();
        assert!((row.logp - (-0.8247)).abs() <= 1e-9);
        for (id, expected) in [
            ("MolMR", 3.6138),
            ("SMR_VSA1", 6.909672289380502),
            ("SMR_VSA2", 0.0),
            ("SlogP_VSA2", 5.478472632147868),
            ("SlogP_VSA5", 1.4311996572326342),
        ] {
            let i = SPECS.iter().position(|s| s.id == id).unwrap();
            let value = row.descriptor_values[i].unwrap();
            assert!(
                (value - expected).abs() <= 1e-9,
                "{smiles} {id}: {value} vs {expected}"
            );
        }
    }
}

#[test]
fn registered_graph_descriptors_match_end_to_end_reference_inputs() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    // End-to-end cases fixed before Graph work. The full 177-case graph test
    // uses independent sanitized atom/bond contexts, including parser gaps;
    // compare_catalog.py separately reports parser and accepted-value gaps.
    let inputs: std::collections::HashSet<_> = include_str!("fixtures/descriptor_phase1.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    for line in include_str!("fixtures/descriptor_graph.tsv")
        .lines()
        .filter(|s| !s.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let (smiles, id) = (f[0], f[1]);
        if !inputs.contains(smiles) {
            continue;
        }
        let expected: f64 = f[2].parse().unwrap();
        if previous != smiles {
            actual = Some(analyze(smiles).unwrap_or_else(|| panic!("fixture rejected: {smiles}")));
            previous = smiles;
        }
        let index = SPECS
            .iter()
            .position(|s| s.id == id)
            .unwrap_or_else(|| panic!("missing registered descriptor: {id}"));
        let value = actual.as_ref().unwrap().descriptor_values[index];
        if !(if expected.is_finite() {
            value.is_some_and(|v| (v - expected).abs() <= 1e-9)
        } else {
            value.is_none()
        }) {
            differences.push(format!("{smiles} {id}: native {value:?}, RDKit {expected}"));
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn fraction_csp3_requires_four_carbon_neighbors_including_attached_h() {
    // RDKit CalcFractionCSP3 counts four-coordinate carbons, not its SP3
    // hybridization enum: isolated [CH3]/[CH2]/[C] are SP3 but excluded.
    for (smiles, expected) in [
        ("[CH3]", 0.),
        ("[CH2]", 0.),
        ("[C]", 0.),
        ("[CH3-]", 0.),
        ("C", 1.),
        ("[H]C([H])([H])[H]", 1.),
        ("[CH2]C", 0.5),
    ] {
        let actual = analyze(smiles).unwrap().fraction_csp3;
        assert_eq!(actual, expected, "{smiles}");
    }
}

#[test]
fn registered_charge_descriptors_match_end_to_end_reference_inputs() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let inputs: std::collections::HashSet<_> = include_str!("fixtures/descriptor_phase1.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    let mut previous = "";
    let mut actual = None;
    let mut differences = Vec::new();
    for line in include_str!("fixtures/descriptor_charge.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = line.split('\t').collect();
        let (smiles, id) = (f[0], f[1]);
        if !inputs.contains(smiles) {
            continue;
        }
        let expected: f64 = f[2].parse().unwrap();
        if previous != smiles {
            actual = Some(analyze(smiles).unwrap());
            previous = smiles;
        }
        let i = SPECS
            .iter()
            .position(|s| s.id == id)
            .unwrap_or_else(|| panic!("missing registered descriptor: {id}"));
        let value = actual.as_ref().unwrap().descriptor_values[i];
        if !(if expected.is_finite() {
            value.is_some_and(|v| {
                (v - expected).abs()
                    <= if id == "NumRadicalElectrons" {
                        0.
                    } else {
                        1e-9
                    }
            })
        } else {
            value.is_none()
        }) {
            differences.push(format!("{smiles} {id}: {value:?} vs {expected}"));
        }
    }
    assert!(
        differences.is_empty(),
        "{} differences:\n{}",
        differences.len(),
        differences.join("\n")
    );
}

#[test]
fn registered_bcut_descriptors_match_end_to_end_reference_inputs() {
    use ducksmiles_smiles::descriptor_profile::SPECS;
    let inputs: std::collections::HashSet<_> = include_str!("fixtures/descriptor_phase1.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    let mut previous = "";
    let mut row = None;
    let mut errors = Vec::new();
    for l in include_str!("fixtures/descriptor_bcut.tsv")
        .lines()
        .filter(|l| !l.starts_with('#'))
    {
        let f: Vec<_> = l.split('\t').collect();
        if !inputs.contains(f[0]) {
            continue;
        }
        if previous != f[0] {
            row = Some(analyze(f[0]).unwrap());
            previous = f[0];
        }
        let i = SPECS
            .iter()
            .position(|s| s.id == f[1])
            .unwrap_or_else(|| panic!("missing registered descriptor: {}", f[1]));
        let expected: f64 = f[2].parse().unwrap();
        let value = row.as_ref().unwrap().descriptor_values[i];
        if !(if expected.is_finite() {
            value.is_some_and(|v| (v - expected).abs() <= 1e-9)
        } else {
            value.is_none()
        }) {
            errors.push(format!("{} {}: {value:?} vs {expected}", f[0], f[1]));
        }
    }
    assert!(errors.is_empty(), "{}", errors.join("\n"));
}
