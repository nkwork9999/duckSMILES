mod common;

use common::{ffi_string, tsv};
use ducksmiles_smiles::*;

fn scaffold(s: &str) -> String {
    ffi_string(ds_murcko_scaffold, s).unwrap()
}

#[test]
fn outputs_reparse() {
    for r in tsv("rdkit_descriptors.tsv") {
        let x = scaffold(&r["smiles"]);
        if !x.is_empty() {
            assert!(
                ds_mol_num_atoms(x.as_ptr(), x.len()) >= 0,
                "{} -> {x}",
                r["smiles"]
            );
        }
    }
}

#[test]
fn exocyclic_double_bonded_framework_atoms_are_retained() {
    for s in [
        "O=C1C=CC(=O)C=C1",
        "O=C1CCCCC1",
        "c1ccc(S(=O)(=O)N2CCCCC2)cc1",
    ] {
        let x = scaffold(s);
        assert!(x.contains('O'), "exocyclic oxygen removed: {s} -> {x}");
    }
    let x = scaffold("c1ccc(S(=O)(=O)N2CCCCC2)cc1");
    assert!(x.contains('S'), "sulfonamide linker removed: {x}");
}

#[test]
fn acyclic_molecules_have_empty_scaffolds() {
    for s in ["C", "CCO", "CC(=O)O", "NCC(=O)NCC(=O)O"] {
        assert_eq!(scaffold(s), "", "{s}");
    }
}

#[test]
fn scaffold_is_idempotent() {
    for r in tsv("rdkit_descriptors.tsv") {
        let x = scaffold(&r["smiles"]);
        if !x.is_empty() {
            assert_eq!(scaffold(&x), x, "{}", r["smiles"]);
        }
    }
}

#[test]
fn scaffold_is_order_invariant() {
    let mut groups = std::collections::BTreeMap::<String, Vec<String>>::new();
    for r in tsv("smiles_variants.tsv") {
        groups
            .entry(r["id"].clone())
            .or_default()
            .push(r["smiles"].clone());
    }
    for (id, v) in groups {
        let first = scaffold(&v[0]);
        for s in &v[1..] {
            assert_eq!(scaffold(s), first, "id={id}, {s}");
        }
    }
}
