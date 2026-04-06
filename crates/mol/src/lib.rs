pub mod parser;

// =============================================================================
// C FFI exports
// =============================================================================

fn as_str(ptr: *const u8, len: usize) -> &'static str {
    unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) }
}

fn write_buf(src: &[u8], out: *mut u8, cap: usize) -> i32 {
    let n = src.len().min(cap);
    unsafe { std::ptr::copy_nonoverlapping(src.as_ptr(), out, n); }
    n as i32
}

/// Count molecules in SDF text
#[unsafe(no_mangle)]
pub extern "C" fn ds_sdf_count(data: *const u8, len: usize) -> i32 {
    if data.is_null() || len == 0 { return 0; }
    parser::parse_sdf(as_str(data, len)).len() as i32
}

/// Parse first MOL block, write formula to out. Returns length or -1.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_formula(data: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => write_buf(mol.formula().as_bytes(), out, cap),
        None => -1,
    }
}

/// Parse first MOL block, return molecular weight. NaN on invalid.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_weight(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => mol.weight(),
        None => f64::NAN,
    }
}

/// Parse first MOL block, return atom count. -1 on invalid.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_num_atoms(data: *const u8, len: usize) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => mol.atoms.len() as i32,
        None => -1,
    }
}

/// Parse first MOL block, return bond count. -1 on invalid.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_num_bonds(data: *const u8, len: usize) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => mol.bonds.len() as i32,
        None => -1,
    }
}

/// Parse first MOL block, write molecule name to out. Returns length or -1.
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_name(data: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => write_buf(mol.name.as_bytes(), out, cap),
        None => -1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const MOL: &str = "ethanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0\n    2.3100    1.3300    0.0000 O   0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\nM  END\n";

    #[test]
    fn test_ffi_formula() {
        let mut buf = [0u8; 64];
        let len = ds_mol_block_formula(MOL.as_ptr(), MOL.len(), buf.as_mut_ptr(), buf.len());
        assert!(len > 0);
        assert_eq!(&buf[..len as usize], b"C2O");
    }

    #[test]
    fn test_ffi_weight() {
        let w = ds_mol_block_weight(MOL.as_ptr(), MOL.len());
        assert!((w - (12.011 * 2.0 + 15.999)).abs() < 0.01);
    }
}
