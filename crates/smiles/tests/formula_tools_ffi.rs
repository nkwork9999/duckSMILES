use ducksmiles_smiles::mass_context::{
    ds_formula_composition_json, ds_formula_dbe, ds_formula_weight, ds_mass_error_ppm,
};

#[test]
fn numeric_formula_abi_returns_values_and_nan_sentinels() {
    // SAFETY: every non-null pointer addresses its declared readable byte slice.
    unsafe {
        assert!((ds_formula_weight(b"C2H6O".as_ptr(), 5) - 46.069).abs() < 1e-9);
        assert_eq!(ds_formula_dbe(b"C6H6".as_ptr(), 4), 4.0);
        for function in [ds_formula_weight, ds_formula_dbe] {
            assert!(function(std::ptr::null(), 1).is_nan());
            assert!(function(b"C".as_ptr(), usize::MAX).is_nan());
            assert!(function([0xffu8].as_ptr(), 1).is_nan());
        }
    }
    assert!((ds_mass_error_ppm(100.001, 100.0) - 10.0).abs() < 1e-9);
    assert!(ds_mass_error_ppm(1.0, 0.0).is_nan());
}

#[test]
fn composition_abi_reports_required_capacity_without_partial_writes() {
    let formula = b"CH3COOH";
    let expected = b"[{\"element\":\"C\",\"isotope\":null,\"count\":2},{\"element\":\"H\",\"isotope\":null,\"count\":4},{\"element\":\"O\",\"isotope\":null,\"count\":2}]";
    // SAFETY: input/output allocations are readable/writable for declared sizes.
    unsafe {
        let needed =
            ds_formula_composition_json(formula.as_ptr(), formula.len(), std::ptr::null_mut(), 0);
        assert_eq!(needed, expected.len() as i32);
        let mut output = vec![0xa5; expected.len()];
        assert_eq!(
            ds_formula_composition_json(
                formula.as_ptr(),
                formula.len(),
                output.as_mut_ptr(),
                output.len() - 1
            ),
            needed
        );
        assert!(output.iter().all(|value| *value == 0xa5));
        assert_eq!(
            ds_formula_composition_json(
                formula.as_ptr(),
                formula.len(),
                output.as_mut_ptr(),
                output.len()
            ),
            needed
        );
        assert_eq!(output, expected);
        assert_eq!(
            ds_formula_composition_json([0xffu8].as_ptr(), 1, std::ptr::null_mut(), 0),
            -1
        );
        assert_eq!(
            ds_formula_composition_json(std::ptr::null(), 1, std::ptr::null_mut(), 0),
            -1
        );
    }
}
