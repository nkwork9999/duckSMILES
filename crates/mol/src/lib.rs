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

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_has_3d(data: *const u8, len: usize) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => if mol.has_3d() { 1 } else { 0 },
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_centroid_x(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.centroid().map(|c| c.0))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_centroid_y(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.centroid().map(|c| c.1))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_centroid_z(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.centroid().map(|c| c.2))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_radius_of_gyration(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.radius_of_gyration())
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_min_x(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[0].0))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_max_x(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[0].1))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_min_y(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[1].0))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_max_y(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[1].1))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_min_z(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[2].0))
        .unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_max_z(data: *const u8, len: usize) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    parser::parse_mol(as_str(data, len))
        .and_then(|mol| mol.coordinate_bounds().map(|b| b[2].1))
        .unwrap_or(f64::NAN)
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

    #[test]
    fn test_ffi_3d_coordinates() {
        assert_eq!(ds_mol_block_has_3d(MOL.as_ptr(), MOL.len()), 0);
        assert!((ds_mol_block_centroid_x(MOL.as_ptr(), MOL.len()) - 1.283333).abs() < 0.01);
        assert!((ds_mol_block_centroid_y(MOL.as_ptr(), MOL.len()) - 0.443333).abs() < 0.01);
        assert_eq!(ds_mol_block_centroid_z(MOL.as_ptr(), MOL.len()), 0.0);
        assert!((ds_mol_block_radius_of_gyration(MOL.as_ptr(), MOL.len()) - 1.147).abs() < 0.01);
        assert_eq!(ds_mol_block_min_x(MOL.as_ptr(), MOL.len()), 0.0);
        assert!((ds_mol_block_max_x(MOL.as_ptr(), MOL.len()) - 2.31).abs() < 0.01);
        assert_eq!(ds_mol_block_min_y(MOL.as_ptr(), MOL.len()), 0.0);
        assert!((ds_mol_block_max_y(MOL.as_ptr(), MOL.len()) - 1.33).abs() < 0.01);
    }
}
