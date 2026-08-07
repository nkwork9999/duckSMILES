mod common;

use common::{ffi_string, on_bits, tsv};
use ducksmiles_smiles::*;

fn mac(s: &str) -> std::collections::BTreeSet<usize> {
    let mut b = [0u8; 21];
    assert_eq!(
        ds_maccs_keys(s.as_ptr(), s.len(), b.as_mut_ptr(), b.len()),
        21
    );
    on_bits(&b, 1, 167)
}
fn counts(s: &str) -> [i32; 10] {
    [
        ds_ring_count(s.as_ptr(), s.len()),
        ds_num_aromatic_rings(s.as_ptr(), s.len()),
        ds_num_aliphatic_rings(s.as_ptr(), s.len()),
        ds_num_saturated_rings(s.as_ptr(), s.len()),
        ds_num_aromatic_carbocycles(s.as_ptr(), s.len()),
        ds_num_aromatic_heterocycles(s.as_ptr(), s.len()),
        ds_num_saturated_carbocycles(s.as_ptr(), s.len()),
        ds_num_saturated_heterocycles(s.as_ptr(), s.len()),
        ds_num_aliphatic_carbocycles(s.as_ptr(), s.len()),
        ds_num_aliphatic_heterocycles(s.as_ptr(), s.len()),
    ]
}

#[test]
fn aromatic_and_kekule_spellings_agree() {
    let desc = tsv("rdkit_descriptors.tsv");
    let vars = tsv("smiles_variants.tsv");
    let mut by_id = std::collections::BTreeMap::<String, Vec<String>>::new();
    for r in vars {
        by_id
            .entry(r["id"].clone())
            .or_default()
            .push(r["smiles"].clone());
    }
    for (id, v) in by_id {
        if !v[0].contains(|c: char| matches!(c, 'c' | 'n' | 'o' | 's')) || v.len() < 2 {
            continue;
        }
        let a = &v[0];
        let k = &v[1];
        assert_eq!(counts(a), counts(k), "ring flavours: {a} / {k}");
        assert_eq!(
            counts(a)[1],
            desc[id.parse::<usize>().unwrap()]["arom_rings"]
                .parse()
                .unwrap(),
            "absolute aromatic rings: {a}"
        );
        assert_eq!(
            counts(k)[1],
            desc[id.parse::<usize>().unwrap()]["arom_rings"]
                .parse()
                .unwrap(),
            "absolute aromatic rings: {k}"
        );
        for (name, x, y, tol) in [
            (
                "tpsa",
                ds_tpsa(a.as_ptr(), a.len()),
                ds_tpsa(k.as_ptr(), k.len()),
                1e-9,
            ),
            (
                "logp",
                ds_logp_crippen(a.as_ptr(), a.len()),
                ds_logp_crippen(k.as_ptr(), k.len()),
                1e-9,
            ),
            (
                "mr",
                ds_mol_mr(a.as_ptr(), a.len()),
                ds_mol_mr(k.as_ptr(), k.len()),
                1e-9,
            ),
            (
                "qed",
                ds_qed(a.as_ptr(), a.len()),
                ds_qed(k.as_ptr(), k.len()),
                1e-9,
            ),
        ] {
            assert!((x - y).abs() <= tol, "{name}: {a} / {k}: {x} != {y}");
        }
        assert_eq!(
            ds_num_h_acceptors(a.as_ptr(), a.len()),
            ds_num_h_acceptors(k.as_ptr(), k.len()),
            "HBA"
        );
        assert_eq!(
            ds_num_h_donors(a.as_ptr(), a.len()),
            ds_num_h_donors(k.as_ptr(), k.len()),
            "HBD"
        );
        assert_eq!(mac(a), mac(k), "MACCS: {a} / {k}");
        assert_eq!(
            ffi_string(ds_murcko_scaffold, a),
            ffi_string(ds_murcko_scaffold, k),
            "Murcko: {a} / {k}"
        );
    }
}
