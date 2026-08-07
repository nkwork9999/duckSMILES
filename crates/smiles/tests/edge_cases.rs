mod common;

use common::StringFn;
use ducksmiles_smiles::*;

#[test]
fn invalid_inputs_are_rejected() {
    for s in [
        "",
        "   ",
        "not smiles",
        "C(",
        "C1CC",
        "C=1CCCCC",
        "C(C))",
        "[C",
    ] {
        assert_eq!(
            ds_mol_num_atoms(s.as_ptr(), s.len()),
            -1,
            "accepted invalid input {s:?}"
        );
        assert!(
            ds_mol_weight(s.as_ptr(), s.len()).is_nan(),
            "weight accepted {s:?}"
        );
        assert!(
            ds_canonical_smiles(s.as_ptr(), s.len(), std::ptr::null_mut(), 0) < 0,
            "canonical accepted {s:?}"
        );
    }
}

#[test]
fn very_long_chain_parses() {
    let s = "C".repeat(4096);
    assert_eq!(ds_mol_num_atoms(s.as_ptr(), s.len()), 4096);
    assert_eq!(ds_mol_num_bonds(s.as_ptr(), s.len()), 4095);
}

#[test]
fn single_atoms_and_disconnected_fragments() {
    for s in ["C", "N", "O", "F", "Cl", "Br", "[Na+]", "[Cl-]"] {
        assert_eq!(ds_mol_num_atoms(s.as_ptr(), s.len()), 1, "{s}");
        assert_eq!(ds_mol_num_bonds(s.as_ptr(), s.len()), 0, "{s}");
    }
    let s = "CC.[Na+].[Cl-]";
    assert_eq!(ds_mol_num_atoms(s.as_ptr(), s.len()), 4);
    assert_eq!(ds_mol_num_bonds(s.as_ptr(), s.len()), 1);
}

fn sizing_error(name: &str, f: StringFn) -> Option<String> {
    let s = "CC(=O)Oc1ccccc1C(=O)O.[Na+]";
    let required = f(s.as_ptr(), s.len(), std::ptr::null_mut(), 0);
    if required < 0 {
        return Some(format!("{name}: sizing call returned {required}"));
    }
    let mut b = vec![0u8; required as usize];
    let written = f(s.as_ptr(), s.len(), b.as_mut_ptr(), b.len());
    (written != required).then(|| format!("{name}: required={required}, written={written}"))
}

#[test]
fn every_string_function_obeys_small_buffer_contract() {
    let funcs: [(&str, StringFn); 12] = [
        ("formula", ds_mol_formula),
        ("canonical", ds_canonical_smiles),
        ("murcko", ds_murcko_scaffold),
        ("generic_scaffold", ds_generic_scaffold),
        ("largest_fragment", ds_largest_fragment),
        ("strip_salts", ds_strip_salts),
        ("neutralize", ds_neutralize_charges),
        ("normalize", ds_normalize_smiles),
        ("fragment_parent", ds_fragment_parent),
        ("admet_json", ds_admet_json),
        ("alerts_json", ds_structural_alerts_json),
        ("ring_systems_json", ds_ring_systems_json),
    ];
    let failures: Vec<_> = funcs
        .into_iter()
        .filter_map(|(name, f)| sizing_error(name, f))
        .collect();
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn undersized_non_null_buffer_returns_required_length() {
    let s = "c1ccccc1";
    let funcs: [(&str, StringFn); 12] = [
        ("formula", ds_mol_formula),
        ("canonical", ds_canonical_smiles),
        ("murcko", ds_murcko_scaffold),
        ("generic_scaffold", ds_generic_scaffold),
        ("largest_fragment", ds_largest_fragment),
        ("strip_salts", ds_strip_salts),
        ("neutralize", ds_neutralize_charges),
        ("normalize", ds_normalize_smiles),
        ("fragment_parent", ds_fragment_parent),
        ("admet_json", ds_admet_json),
        ("alerts_json", ds_structural_alerts_json),
        ("ring_systems_json", ds_ring_systems_json),
    ];
    let mut failures = Vec::new();
    for (name, f) in funcs {
        let need = f(s.as_ptr(), s.len(), std::ptr::null_mut(), 0);
        let mut byte = 0u8;
        let got = f(s.as_ptr(), s.len(), &mut byte, 1);
        if got != need {
            failures.push(format!("{name}: sizing={need}, undersized={got}"));
        }
    }
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}
