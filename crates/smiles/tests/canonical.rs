mod common;

use common::{ffi_string, tsv};
use ducksmiles_smiles::*;
use std::collections::{BTreeMap, HashMap};

fn canonical(s: &str) -> String {
    ffi_string(ds_canonical_smiles, s).unwrap()
}

fn groups() -> BTreeMap<String, Vec<String>> {
    let mut groups = BTreeMap::new();
    for r in tsv("smiles_variants.tsv") {
        groups
            .entry(r["id"].clone())
            .or_insert_with(Vec::new)
            .push(r["smiles"].clone());
    }
    groups
}

#[test]
fn order_invariance() {
    let mut failures = Vec::new();
    for (id, spellings) in groups() {
        let values: Vec<_> = spellings.iter().map(|s| (s, canonical(s))).collect();
        if values.iter().any(|x| x.1 != values[0].1) {
            failures.push(format!("id={id}: {values:?}"));
        }
    }
    assert!(
        failures.is_empty(),
        "{} invariant groups failed\n{}",
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
fn idempotence() {
    for r in tsv("rdkit_descriptors.tsv") {
        let a = canonical(&r["smiles"]);
        assert_eq!(canonical(&a), a, "{}", r["smiles"]);
    }
}

fn ints(s: &str) -> [i32; 11] {
    [
        ds_mol_num_atoms(s.as_ptr(), s.len()),
        ds_mol_num_bonds(s.as_ptr(), s.len()),
        ds_ring_count(s.as_ptr(), s.len()),
        ds_num_aromatic_rings(s.as_ptr(), s.len()),
        ds_num_aliphatic_rings(s.as_ptr(), s.len()),
        ds_num_saturated_rings(s.as_ptr(), s.len()),
        ds_num_aromatic_carbocycles(s.as_ptr(), s.len()),
        ds_num_aromatic_heterocycles(s.as_ptr(), s.len()),
        ds_num_saturated_carbocycles(s.as_ptr(), s.len()),
        ds_num_saturated_heterocycles(s.as_ptr(), s.len()),
        ds_num_aliphatic_heterocycles(s.as_ptr(), s.len()),
    ]
}

#[test]
fn round_trip_preserves_structure() {
    for r in tsv("rdkit_descriptors.tsv") {
        let s = &r["smiles"];
        let c = canonical(s);
        assert!(
            ds_mol_num_atoms(c.as_ptr(), c.len()) >= 0,
            "canonical did not parse: {s} -> {c}"
        );
        assert_eq!(ints(s), ints(&c), "ring/count round trip: {s} -> {c}");
        assert_eq!(
            ffi_string(ds_mol_formula, s).unwrap(),
            ffi_string(ds_mol_formula, &c).unwrap(),
            "formula: {s}"
        );
        assert!(
            (ds_mol_weight(s.as_ptr(), s.len()) - ds_mol_weight(c.as_ptr(), c.len())).abs() < 1e-9,
            "weight: {s}"
        );
    }
}

#[test]
fn distinct_molecules_have_distinct_canonical_forms() {
    let mut seen: HashMap<String, String> = HashMap::new();
    for r in tsv("rdkit_descriptors.tsv") {
        let s = &r["smiles"];
        let c = canonical(s);
        if let Some(old) = seen.insert(c.clone(), s.clone()) {
            panic!("distinct molecules collide: {old} and {s} -> {c}");
        }
    }
}

#[test]
fn stereo_is_invariant_and_enantiomers_differ() {
    let rows = tsv("rdkit_descriptors.tsv");
    let pairs = [
        ("C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"),
        ("N[C@@H](Cc1ccccc1)C(=O)O", "N[C@H](Cc1ccccc1)C(=O)O"),
    ];
    for (a, b) in pairs {
        assert_ne!(
            canonical(a),
            canonical(b),
            "enantiomers collapsed: {a} / {b}"
        );
    }
    let stereo_ids: Vec<String> = rows
        .iter()
        .enumerate()
        .filter(|(_, r)| r["smiles"].contains('@'))
        .map(|(i, _)| i.to_string())
        .collect();
    let all = groups();
    for id in stereo_ids {
        let spellings = &all[&id];
        let first = canonical(&spellings[0]);
        for s in spellings {
            assert_eq!(canonical(s), first, "stereo spelling {s}");
        }
    }
}
