mod common;
use ducksmiles_smiles::*;

#[test]
fn fixture_is_complete_and_contains_regressions() {
    let rows = common::tsv("rdkit_profile.tsv");
    assert_eq!(rows.len(), 4510);
    let inputs: std::collections::HashSet<_> = rows.iter().map(|r| r["smiles"].as_str()).collect();
    assert_eq!(inputs.len(), rows.len());
    for s in [
        "C1=CNC=C1",
        "C1=C[N]C=C1",
        "C1=CC=CC=CC=CC=C1",
        "O=C1C=CC(=O)C2=C1OC=CO2",
        "[13CH3][18OH]",
        "[235U]",
        "[H]O[H]",
    ] {
        assert!(inputs.contains(s), "missing regression: {s}");
    }
    for r in &rows {
        assert_eq!(r.len(), 22);
    }
}

#[test]
fn profile_hydrogen_partition_and_invalid_inputs() {
    for s in ["CCO", "[H]O[H]", "C1=CNC=C1", "c1cc[nH]c1", "[NH4+]"] {
        let (p, n) = (s.as_ptr(), s.len());
        assert_eq!(
            ds_mol_num_total_h(p, n),
            ds_mol_num_explicit_h(p, n) + ds_mol_num_implicit_h(p, n)
        );
    }
    let bad_utf8 = [0xffu8];
    for (p, n) in [
        (std::ptr::null(), 0),
        (bad_utf8.as_ptr(), 1),
        (b"invalid".as_ptr(), 7),
    ] {
        assert_eq!(ds_mol_formal_charge(p, n), i32::MIN);
        assert_eq!(ds_mol_num_total_h(p, n), -1);
        assert!(ds_mol_heavy_atom_mass(p, n).is_nan());
        assert!(ds_mol_mean_degree(p, n).is_nan());
    }
}
#[test]
fn mol_num_fragments() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_fragments(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_fragments"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_fragments", &bad, rows.len())
    );
}
#[test]
fn mol_formal_charge() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_formal_charge(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_formal_charge"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_formal_charge", &bad, rows.len())
    );
}
#[test]
fn mol_num_explicit_h() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_explicit_h(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_explicit_h"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_explicit_h", &bad, rows.len())
    );
}
#[test]
fn mol_num_implicit_h() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_implicit_h(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_implicit_h"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_implicit_h", &bad, rows.len())
    );
}
#[test]
fn mol_num_total_h() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_total_h(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_total_h"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_total_h", &bad, rows.len())
    );
}
#[test]
fn mol_num_single_bonds() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_single_bonds(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_single_bonds"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_single_bonds", &bad, rows.len())
    );
}
#[test]
fn mol_num_double_bonds() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_double_bonds(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_double_bonds"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_double_bonds", &bad, rows.len())
    );
}
#[test]
fn mol_num_triple_bonds() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_triple_bonds(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_triple_bonds"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_triple_bonds", &bad, rows.len())
    );
}
#[test]
fn mol_num_aromatic_bonds() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_aromatic_bonds(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_aromatic_bonds"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_aromatic_bonds", &bad, rows.len())
    );
}
#[test]
fn mol_num_ring_atoms() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_ring_atoms(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_ring_atoms"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_ring_atoms", &bad, rows.len())
    );
}
#[test]
fn mol_num_ring_bonds() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_ring_bonds(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_ring_bonds"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_ring_bonds", &bad, rows.len())
    );
}
#[test]
fn mol_largest_ring_size() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_largest_ring_size(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_largest_ring_size"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_largest_ring_size", &bad, rows.len())
    );
}
#[test]
fn mol_num_aromatic_atoms() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_aromatic_atoms(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_aromatic_atoms"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_aromatic_atoms", &bad, rows.len())
    );
}
#[test]
fn mol_num_carbons() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_carbons(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_carbons"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_carbons", &bad, rows.len())
    );
}
#[test]
fn mol_num_nitrogens() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_nitrogens(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_nitrogens"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_nitrogens", &bad, rows.len())
    );
}
#[test]
fn mol_num_oxygens() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_oxygens(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_oxygens"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_oxygens", &bad, rows.len())
    );
}
#[test]
fn mol_num_halogens() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_num_halogens(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_num_halogens"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 0.0 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_num_halogens", &bad, rows.len())
    );
}
#[test]
fn mol_heteroatom_fraction() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_heteroatom_fraction(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_heteroatom_fraction"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 1e-9 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_heteroatom_fraction", &bad, rows.len())
    );
}
#[test]
fn mol_aromatic_fraction() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_aromatic_fraction(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_aromatic_fraction"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 1e-9 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_aromatic_fraction", &bad, rows.len())
    );
}
#[test]
fn mol_heavy_atom_mass() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_heavy_atom_mass(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_heavy_atom_mass"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 1e-9 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_heavy_atom_mass", &bad, rows.len())
    );
}
#[test]
fn mol_mean_degree() {
    let rows = common::tsv("rdkit_profile.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ds_mol_mean_degree(s.as_ptr(), s.len()) as f64;
        let want: f64 = r["mol_mean_degree"].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > 1e-9 {
            bad.push(format!("{s}: ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        common::mismatch_report("mol_mean_degree", &bad, rows.len())
    );
}
