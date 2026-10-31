use ducksmiles_smiles::mass_context::{
    formula_composition_json, formula_dbe, formula_weight, mass_error_ppm,
};

#[test]
fn average_formula_weight_uses_isotope_mass_only_for_labelled_atoms() {
    for (formula, expected) in [
        ("C2H6O", 46.069),
        ("[13C]2H6O", 48.05370968),
        ("Li", 6.941),
        ("CH3COOH", 60.052),
    ] {
        assert!(
            (formula_weight(formula).unwrap() - expected).abs() < 1e-9,
            "{formula}"
        );
    }
    for formula in ["Xx", "[999C]", "C0", "C1000001", "NH4+"] {
        assert!(formula_weight(formula).is_none());
    }
}

#[test]
fn formula_composition_aggregates_repeated_tokens_and_keeps_isotopes_distinct() {
    assert_eq!(formula_composition_json("[13C]CH3[13C]OH"), Some("[{\"element\":\"C\",\"isotope\":null,\"count\":1},{\"element\":\"C\",\"isotope\":13,\"count\":2},{\"element\":\"H\",\"isotope\":null,\"count\":4},{\"element\":\"O\",\"isotope\":null,\"count\":1}]".to_owned()));
    for formula in ["", "[13CH4]", "[0C]", "C H", "C999999[2H]2"] {
        assert!(formula_composition_json(formula).is_none());
    }
}

#[test]
fn dbe_counts_isotope_elements_and_rejects_unsupported_chemistry() {
    for (formula, expected) in [
        ("C6H6", 4.0),
        ("C8H10N4O2", 6.0),
        ("C6H5Cl", 4.0),
        ("[13C]6[2H]6", 4.0),
        ("C2H6O", 0.0),
        ("C6H12O6", 1.0),
    ] {
        assert_eq!(formula_dbe(formula), Some(expected), "{formula}");
    }
    for formula in ["C2H7", "CH6", "NaCl", "Fe", "C2H5P", "C2H6S", "NH4+"] {
        assert!(formula_dbe(formula).is_none(), "{formula}");
    }
}

#[test]
fn ppm_error_is_signed_reference_relative_and_retains_small_differences() {
    assert!((mass_error_ppm(100.001, 100.0).unwrap() - 10.0).abs() < 1e-9);
    assert!((mass_error_ppm(99.999, 100.0).unwrap() + 10.0).abs() < 1e-9);
    assert_eq!(mass_error_ppm(1e300, 1e300), Some(0.0));
    assert_eq!(
        mass_error_ppm(f64::from_bits(1), f64::from_bits(1)),
        Some(0.0)
    );
    assert!((mass_error_ppm(1.000000000001, 1.0).unwrap() - 0.000001000088900582341).abs() < 1e-18);
    for (observed, reference) in [
        (0.0, 100.0),
        (100.0, 0.0),
        (-1.0, 1.0),
        (f64::NAN, 1.0),
        (1.0, f64::INFINITY),
        (f64::MAX, 1e-300),
    ] {
        assert!(mass_error_ppm(observed, reference).is_none());
    }
}
