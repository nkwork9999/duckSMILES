mod common;

use common::{ffi_string, mismatch_report, tsv};
use ducksmiles_smiles::*;

fn check_i32(name: &str, col: &str, f: extern "C" fn(*const u8, usize) -> i32) {
    let rows = tsv("rdkit_descriptors.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = f(s.as_ptr(), s.len());
        let want: i32 = r[col].parse().unwrap();
        if got != want {
            bad.push(format!("{s}, ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        mismatch_report(name, &bad, rows.len())
    );
}

fn check_f64(name: &str, col: &str, tolerance: f64, f: extern "C" fn(*const u8, usize) -> f64) {
    let rows = tsv("rdkit_descriptors.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = f(s.as_ptr(), s.len());
        let want: f64 = r[col].parse().unwrap();
        if !got.is_finite() || (got - want).abs() > tolerance {
            bad.push(format!("{s}, ducksmiles={got:.6}, rdkit={want:.6}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        mismatch_report(name, &bad, rows.len())
    );
}

/// Exact string comparison, for outputs with one correct spelling.
fn check_string(name: &str, col: &str, f: common::StringFn) {
    let rows = tsv("rdkit_descriptors.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ffi_string(f, s).unwrap_or_else(|e| format!("ERR({e})"));
        let want = &r[col];
        if &got != want {
            bad.push(format!("{s}, ducksmiles={got}, rdkit={want}"));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        mismatch_report(name, &bad, rows.len())
    );
}

/// Compare a SMILES-returning function against RDKit *structurally*.
///
/// Two canonicalisers are free to pick different spellings of the same
/// molecule, so `C(C)C` and `CCC` are both correct answers for propane and a
/// byte comparison against RDKit's string would only be measuring whose
/// traversal order we copied. What has to hold is that RDKit's spelling and
/// ours describe the same graph: feeding RDKit's answer back through
/// `ds_canonical_smiles` must land on exactly the string we produced.
fn check_smiles_structure(name: &str, col: &str, f: common::StringFn) {
    let rows = tsv("rdkit_descriptors.tsv");
    let mut bad = Vec::new();
    for r in &rows {
        let s = &r["smiles"];
        let got = ffi_string(f, s).unwrap_or_else(|e| format!("ERR({e})"));
        let want = &r[col];
        if want.is_empty() || got.is_empty() {
            if got != *want {
                bad.push(format!("{s}, ducksmiles={got:?}, rdkit={want:?}"));
            }
            continue;
        }
        let ours = ffi_string(ds_canonical_smiles, &got).unwrap_or_else(|e| format!("ERR({e})"));
        let theirs = ffi_string(ds_canonical_smiles, want).unwrap_or_else(|e| format!("ERR({e})"));
        if ours != theirs {
            bad.push(format!(
                "{s}, ducksmiles={got} (canonical {ours}), rdkit={want} (canonical {theirs})"
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "{}",
        mismatch_report(name, &bad, rows.len())
    );
}

macro_rules! itest {
    ($n:ident,$c:literal,$f:ident) => {
        #[test]
        fn $n() {
            check_i32(stringify!($n), $c, $f)
        }
    };
}
macro_rules! ftest {
    ($n:ident,$c:literal,$t:expr,$f:ident) => {
        #[test]
        fn $n() {
            check_f64(stringify!($n), $c, $t, $f)
        }
    };
}

itest!(num_atoms, "num_atoms", ds_mol_num_atoms);
itest!(num_bonds, "num_bonds", ds_mol_num_bonds);
ftest!(mol_weight, "mol_weight", 0.05, ds_mol_weight);
ftest!(exact_mass, "exact_mass", 0.01, ds_mol_exact_mass);
ftest!(tpsa, "tpsa", 0.1, ds_tpsa);
ftest!(logp, "logp", 0.05, ds_logp_crippen);
ftest!(mol_mr, "mol_mr", 0.5, ds_mol_mr);
ftest!(qed, "qed", 0.02, ds_qed);
itest!(h_donors, "hbd", ds_num_h_donors);
itest!(h_acceptors, "hba", ds_num_h_acceptors);
itest!(rotatable_bonds, "rotb", ds_num_rotatable_bonds);
itest!(heteroatoms, "heteroatoms", ds_num_heteroatoms);
ftest!(fraction_csp3, "fcsp3", 1e-3, ds_fraction_csp3);
itest!(ring_count, "ring_count", ds_ring_count);
itest!(aromatic_rings, "arom_rings", ds_num_aromatic_rings);
itest!(aliphatic_rings, "aliph_rings", ds_num_aliphatic_rings);
itest!(saturated_rings, "satur_rings", ds_num_saturated_rings);
itest!(
    aromatic_carbocycles,
    "arom_carbo",
    ds_num_aromatic_carbocycles
);
itest!(
    aromatic_heterocycles,
    "arom_hetero",
    ds_num_aromatic_heterocycles
);
itest!(
    saturated_carbocycles,
    "satur_carbo",
    ds_num_saturated_carbocycles
);
itest!(
    saturated_heterocycles,
    "satur_hetero",
    ds_num_saturated_heterocycles
);
itest!(
    aliphatic_carbocycles,
    "aliph_carbo",
    ds_num_aliphatic_carbocycles
);
itest!(
    aliphatic_heterocycles,
    "aliph_hetero",
    ds_num_aliphatic_heterocycles
);
itest!(lipinski_violations, "lipinski_viol", ds_lipinski_violations);
itest!(
    structural_alert_count,
    "structural_alerts",
    ds_structural_alert_count
);

#[test]
fn formula() {
    check_string("formula", "formula", ds_mol_formula)
}
#[test]
fn canonical_smiles() {
    check_smiles_structure("canonical_smiles", "canonical", ds_canonical_smiles)
}
#[test]
fn murcko_scaffold() {
    check_smiles_structure("murcko_scaffold", "murcko", ds_murcko_scaffold)
}
