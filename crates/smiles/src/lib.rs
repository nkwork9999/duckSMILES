mod logp_crippen;
mod morgan;
mod parser;
mod smarts;
mod tanimoto;
mod weights;

use logp_crippen::calc_logp;
use morgan::morgan_bits;
use parser::parse;
use tanimoto::tanimoto_bit;

// =============================================================================
// C FFI exports
// =============================================================================

/// Returns 1 if valid SMILES, 0 otherwise
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_is_valid(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    if parse(s).is_some() { 1 } else { 0 }
}

/// Returns heavy atom count, or -1 on invalid
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_num_atoms(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => mol.heavy_atom_count() as i32,
        None => -1,
    }
}

/// Returns bond count, or -1 on invalid
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_num_bonds(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => mol.bond_count as i32,
        None => -1,
    }
}

/// Writes molecular formula to buffer. Returns length written, or -1 on invalid.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_formula(
    ptr: *const u8, len: usize,
    out: *mut u8, out_cap: usize,
) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let formula = mol.formula();
            let bytes = formula.as_bytes();
            let n = bytes.len().min(out_cap);
            unsafe { std::ptr::copy_nonoverlapping(bytes.as_ptr(), out, n); }
            n as i32
        }
        None => -1,
    }
}

/// Returns molecular weight (average), or NaN on invalid
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_weight(ptr: *const u8, len: usize) -> f64 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => mol.molecular_weight(),
        None => f64::NAN,
    }
}

/// Returns monoisotopic mass, or NaN on invalid
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_exact_mass(ptr: *const u8, len: usize) -> f64 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => mol.exact_mass(),
        None => f64::NAN,
    }
}

/// Returns Wildman-Crippen LogP (RDKit-compatible), or NaN on invalid SMILES.
#[unsafe(no_mangle)]
pub extern "C" fn ds_logp_crippen(ptr: *const u8, len: usize) -> f64 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => calc_logp(&mol),
        None => f64::NAN,
    }
}

/// Writes SMILES with explicit H atoms to buffer. Returns length written, or -1 on invalid.
/// Result is a verbose bracket-form SMILES that round-trips through parse.
#[unsafe(no_mangle)]
pub extern "C" fn ds_add_hydrogens(
    ptr: *const u8, len: usize,
    out: *mut u8, out_cap: usize,
) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let smi = mol.with_explicit_hydrogens().to_smiles_verbose();
            let bytes = smi.as_bytes();
            let n = bytes.len().min(out_cap);
            unsafe { std::ptr::copy_nonoverlapping(bytes.as_ptr(), out, n); }
            n as i32
        }
        None => -1,
    }
}

/// Tanimoto similarity between two equal-length fingerprint BLOBs.
///
/// Returns:
///   - `popcount(a & b) / popcount(a | b)` for normal inputs (range `[0.0, 1.0]`).
///   - `0.0` if both blobs are all-zero (matches RDKit `CalcBitmapTanimoto`).
///   - `NaN` if the lengths differ — the DuckDB extension surfaces this as
///     an `InvalidInputException` so the user sees a clear error rather than
///     a silent NULL.
#[unsafe(no_mangle)]
pub extern "C" fn ds_tanimoto_bit(
    a_ptr: *const u8, a_len: usize,
    b_ptr: *const u8, b_len: usize,
) -> f64 {
    if a_len != b_len {
        return f64::NAN;
    }
    if a_len == 0 {
        // Both empty → degenerate match RDKit's "both zero" path.
        return 0.0;
    }
    let a = unsafe { std::slice::from_raw_parts(a_ptr, a_len) };
    let b = unsafe { std::slice::from_raw_parts(b_ptr, b_len) };
    tanimoto_bit(a, b)
}

/// Writes Morgan/ECFP fingerprint bits to buffer. Returns bytes written
/// (= ceil(n_bits/8)), or -1 on invalid SMILES / buffer too small.
///
/// `radius`: ECFPn corresponds to radius = n/2 (so ECFP4 → radius=2).
/// `n_bits`: target bit-vector width (e.g. 2048). Capped to `out_cap * 8`.
#[unsafe(no_mangle)]
pub extern "C" fn ds_morgan_fp_bits(
    ptr: *const u8, len: usize,
    radius: u32, n_bits: u32,
    out: *mut u8, out_cap: usize,
) -> i32 {
    if n_bits == 0 { return -1; }
    let n_bytes = ((n_bits as usize) + 7) / 8;
    if n_bytes > out_cap { return -1; }
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let bits = morgan_bits(&mol, radius, n_bits);
            unsafe { std::ptr::copy_nonoverlapping(bits.as_ptr(), out, n_bytes); }
            n_bytes as i32
        }
        None => -1,
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_water() {
        let mol = parse("O").unwrap();
        assert_eq!(mol.formula(), "H2O");
        assert_eq!(mol.heavy_atom_count(), 1);
        assert!((mol.molecular_weight() - 18.015).abs() < 0.01);
    }

    #[test]
    fn test_ethanol() {
        let mol = parse("CCO").unwrap();
        assert_eq!(mol.formula(), "C2H6O");
        assert_eq!(mol.heavy_atom_count(), 3);
        assert_eq!(mol.bond_count, 2);
    }

    #[test]
    fn test_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(mol.formula(), "C6H6");
        assert_eq!(mol.heavy_atom_count(), 6);
        assert_eq!(mol.bond_count, 6);
    }

    #[test]
    fn test_acetic_acid() {
        let mol = parse("CC(=O)O").unwrap();
        assert_eq!(mol.formula(), "C2H4O2");
    }

    #[test]
    fn test_salt() {
        let mol = parse("[Na+].[Cl-]").unwrap();
        assert_eq!(mol.formula(), "ClNa");
    }

    #[test]
    fn test_aspirin() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(mol.formula(), "C9H8O4");
    }

    #[test]
    fn test_invalid() {
        assert!(parse("").is_none());
        assert!(parse("not_a_molecule").is_none());
        assert!(parse("C(C").is_none());
    }
}
