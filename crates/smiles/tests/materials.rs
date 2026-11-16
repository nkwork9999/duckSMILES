use ducksmiles_smiles::materials::*;
#[test]
fn fractional_and_grouped_compositions_preserve_stoichiometry() {
    let c = composition("Cs0.5Rb0.5Pb(I0.8Br0.2)3").unwrap();
    assert_eq!(c["Cs"], 0.5);
    assert_eq!(c["Pb"], 1.);
    assert!((c["I"] - 2.4).abs() < 1e-12);
    assert!((c["Br"] - 0.6).abs() < 1e-12);
    assert_eq!(
        composition("Ca(OH)2").unwrap(),
        composition("CaO2H2").unwrap()
    );
    let j = composition_json("CsPbI3").unwrap();
    assert!(j.contains("\"formula_amount\":5"));
    assert!(j.contains("\"atomic_fraction\":0.6"));
}
#[test]
fn composition_rejects_partial_unknown_overflow_and_ambiguous_inputs() {
    for s in [
        "",
        "()",
        "Cs0PbI3",
        "Cs-1PbI3",
        "Cs1e3",
        "XxO",
        "Cs(PbI3",
        "CsPbI3)",
        "Cs PbI3",
        "[13C]",
        "NH4+",
        "CuSO4·5H2O",
        "C1000001",
        "C.0",
        "C1.",
    ] {
        assert!(composition(s).is_none(), "{s}");
    }
    assert!(composition(&format!("{}H{}", "(".repeat(18), ")".repeat(18))).is_none());
}
#[test]
fn solvent_profile_is_bounded_and_not_a_bulk_property_prediction() {
    for s in ["CS(C)=O", "CN(C)C=O", "CCO", "O", "CC#N", "[NH3+]C"] {
        let j = solvent_json(s).expect(s);
        assert!(j.contains("\"MolWt\":"));
        assert!(j.contains("\"MolMR\":"));
        assert!(!j.contains("NaN"));
    }
    for s in ["", "C1CC", "CCO.O", "[Na+].[Cl-]", "[H]"] {
        assert!(solvent_json(s).is_none(), "{s}");
    }
}
#[test]
fn public_abi_capacity_null_and_utf8_contracts() {
    for f in [ds_composition_descriptors_json, ds_solvent_descriptors_json] {
        let s = b"CCO"; // composition C,C,O and molecular ethanol are both valid.
        unsafe {
            let n = f(s.as_ptr(), s.len(), std::ptr::null_mut(), 0);
            assert!(n > 0);
            let mut out = vec![0xa5; n as usize];
            assert_eq!(f(s.as_ptr(), s.len(), out.as_mut_ptr(), out.len() - 1), n);
            assert!(out.iter().all(|v| *v == 0xa5));
            assert_eq!(f(s.as_ptr(), s.len(), out.as_mut_ptr(), out.len()), n);
            assert!(std::str::from_utf8(&out).unwrap().starts_with('{'));
            assert_eq!(f(std::ptr::null(), 1, std::ptr::null_mut(), 0), -1);
            assert_eq!(f([0xff].as_ptr(), 1, std::ptr::null_mut(), 0), -1);
            assert_eq!(f(s.as_ptr(), usize::MAX, std::ptr::null_mut(), 0), -1);
            assert_eq!(f(s.as_ptr(), s.len(), std::ptr::null_mut(), 1), -1);
        }
    }
}

#[test]
fn composition_underflow_never_emits_nonfinite_json_or_zero_amounts() {
    let tiny = format!("0.{}5", "0".repeat(323));
    assert!(composition_json(&format!("C{tiny}H999999")).is_none());
    assert!(composition(&format!("(C{tiny}){tiny}")).is_none());
}
