use ducksmiles_smiles::{ds_mol_has_substructure, ds_mol_substructure_count, ds_tanimoto_bit};

#[test]
fn word_popcount_agrees_with_byte_reference() {
    let left: [u8; 11] = [0, 1, 2, 5, 9, 10, 11, 80, 128, 255, 254];
    let right: [u8; 11] = [0, 2, 4, 5, 13, 1, 13, 8, 127, 254, 255];
    let zero: [u8; 11] = [0; 11];
    for len in 0..10 {
        for (a, b) in [
            (&left[1..1 + len], &right[1..1 + len]),
            (&zero[1..1 + len], &zero[1..1 + len]),
        ] {
            let intersection: u32 = a.iter().zip(b).map(|(x, y)| (x & y).count_ones()).sum();
            let union: u32 = a.iter().zip(b).map(|(x, y)| (x | y).count_ones()).sum();
            let expected = if union == 0 {
                0.0
            } else {
                f64::from(intersection) / f64::from(union)
            };
            assert_eq!(
                ds_tanimoto_bit(a.as_ptr(), a.len(), b.as_ptr(), b.len()),
                expected
            );
        }
    }
}

#[test]
fn existence_agrees_with_nonempty_matches_on_symmetric_inputs() {
    for (molecule, pattern) in [("c1ccccc1", "cc"), ("C1CCCCC1", "CCC"), ("C1CCCCC1", "N")] {
        let has = ds_mol_has_substructure(
            molecule.as_ptr(),
            molecule.len(),
            pattern.as_ptr(),
            pattern.len(),
        );
        let count = ds_mol_substructure_count(
            molecule.as_ptr(),
            molecule.len(),
            pattern.as_ptr(),
            pattern.len(),
        );
        assert!(count >= 0);
        assert_eq!(has, i32::from(count > 0));
    }
}
