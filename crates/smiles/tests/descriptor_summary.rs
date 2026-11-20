use ducksmiles_smiles::{ds_descriptor_summary, ds_mol_num_atoms, ds_mol_weight, ds_ring_count};

#[test]
fn bundle_agrees_with_existing_descriptor_functions() {
    for input in ["O", "CCO", "c1ccccc1", "[13CH4]", "invalid", ""] {
        let (ptr, len) = (input.as_ptr(), input.len());
        let summary = ds_descriptor_summary(ptr, len);
        let expected_weight = ds_mol_weight(ptr, len);
        if expected_weight.is_nan() {
            assert!(summary.weight.is_nan(), "{input}");
        } else {
            assert_eq!(summary.weight, expected_weight, "{input}");
        }
        assert_eq!(summary.atom_count, ds_mol_num_atoms(ptr, len), "{input}");
        assert_eq!(summary.ring_count, ds_ring_count(ptr, len), "{input}");
    }
}
