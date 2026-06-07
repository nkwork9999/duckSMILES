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

fn write_buf_required(src: &[u8], out: *mut u8, cap: usize) -> i32 {
    if out.is_null() || cap < src.len() {
        return src.len() as i32;
    }
    unsafe { std::ptr::copy_nonoverlapping(src.as_ptr(), out, src.len()); }
    src.len() as i32
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
pub extern "C" fn ds_mol_block_property(
    data: *const u8,
    len: usize,
    key: *const u8,
    key_len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    if data.is_null() || key.is_null() || len == 0 {
        return -1;
    }
    let key = as_str(key, key_len);
    match parser::parse_mol(as_str(data, len)).and_then(|mol| mol.property(key).map(str::to_owned))
    {
        Some(value) => write_buf_required(value.as_bytes(), out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_block_properties_json(
    data: *const u8,
    len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    if data.is_null() || len == 0 {
        return -1;
    }
    match parser::parse_mol(as_str(data, len)) {
        Some(mol) => write_buf_required(mol.properties_json().as_bytes(), out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_sdf_property(
    data: *const u8,
    len: usize,
    record_index: i32,
    key: *const u8,
    key_len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    if data.is_null() || key.is_null() || len == 0 || record_index <= 0 {
        return -1;
    }
    let key = as_str(key, key_len);
    let mols = parser::parse_sdf(as_str(data, len));
    let idx = (record_index - 1) as usize;
    match mols.get(idx).and_then(|mol| mol.property(key).map(str::to_owned)) {
        Some(value) => write_buf_required(value.as_bytes(), out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_sdf_properties_json(
    data: *const u8,
    len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    if data.is_null() || len == 0 {
        return -1;
    }
    let mols = parser::parse_sdf(as_str(data, len));
    write_buf_required(parser::sdf_properties_json(&mols).as_bytes(), out, cap)
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
    const MOL_V3000: &str = "ethanol_v3000\n  test\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.0000 0.0000 0.0000 0\nM  V30 2 C 1.5400 0.0000 0.0000 0 CFG=0\nM  V30 3 O 2.3100 1.3300 1.0000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=0\nM  V30 2 2 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";

    #[test]
    fn test_ffi_formula() {
        let mut buf = [0u8; 64];
        let len = ds_mol_block_formula(MOL.as_ptr(), MOL.len(), buf.as_mut_ptr(), buf.len());
        assert!(len > 0);
        assert_eq!(&buf[..len as usize], b"C2O");

        let len = ds_mol_block_formula(MOL_V3000.as_ptr(), MOL_V3000.len(), buf.as_mut_ptr(), buf.len());
        assert!(len > 0);
        assert_eq!(&buf[..len as usize], b"C2O");
    }

    #[test]
    fn test_ffi_weight() {
        let w = ds_mol_block_weight(MOL.as_ptr(), MOL.len());
        assert!((w - (12.011 * 2.0 + 15.999)).abs() < 0.01);
        let w_v3000 = ds_mol_block_weight(MOL_V3000.as_ptr(), MOL_V3000.len());
        assert!((w_v3000 - (12.011 * 2.0 + 15.999)).abs() < 0.01);
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

        assert_eq!(ds_mol_block_has_3d(MOL_V3000.as_ptr(), MOL_V3000.len()), 1);
        assert_eq!(ds_mol_block_num_atoms(MOL_V3000.as_ptr(), MOL_V3000.len()), 3);
        assert_eq!(ds_mol_block_num_bonds(MOL_V3000.as_ptr(), MOL_V3000.len()), 2);
        assert!((ds_mol_block_centroid_x(MOL_V3000.as_ptr(), MOL_V3000.len()) - 1.283333).abs() < 0.01);
        assert!((ds_mol_block_centroid_y(MOL_V3000.as_ptr(), MOL_V3000.len()) - 0.443333).abs() < 0.01);
        assert!((ds_mol_block_centroid_z(MOL_V3000.as_ptr(), MOL_V3000.len()) - 0.333333).abs() < 0.01);
        assert_eq!(ds_mol_block_min_z(MOL_V3000.as_ptr(), MOL_V3000.len()), 0.0);
        assert_eq!(ds_mol_block_max_z(MOL_V3000.as_ptr(), MOL_V3000.len()), 1.0);
    }

    fn read_ffi_string<F: FnMut(*mut u8, usize) -> i32>(mut f: F) -> Option<String> {
        let needed = f(std::ptr::null_mut(), 0);
        if needed < 0 {
            return None;
        }
        let mut buf = vec![0u8; needed as usize];
        let n = f(buf.as_mut_ptr(), buf.len());
        assert_eq!(n, needed);
        Some(String::from_utf8(buf).unwrap())
    }

    #[test]
    fn test_ffi_mol_block_property_and_json() {
        let mol = format!("{}> <ID>\n123\n\n> <NOTE>\nline one\nline two\n\n", MOL);
        let id = "ID";
        let note = "NOTE";

        let id_value = read_ffi_string(|out, cap| {
            ds_mol_block_property(
                mol.as_ptr(), mol.len(),
                id.as_ptr(), id.len(),
                out, cap,
            )
        }).unwrap();
        assert_eq!(id_value, "123");

        let note_value = read_ffi_string(|out, cap| {
            ds_mol_block_property(
                mol.as_ptr(), mol.len(),
                note.as_ptr(), note.len(),
                out, cap,
            )
        }).unwrap();
        assert_eq!(note_value, "line one\nline two");

        let json = read_ffi_string(|out, cap| {
            ds_mol_block_properties_json(mol.as_ptr(), mol.len(), out, cap)
        }).unwrap();
        assert_eq!(
            json,
            "[{\"name\":\"ID\",\"value\":\"123\"},{\"name\":\"NOTE\",\"value\":\"line one\\nline two\"}]"
        );

        let mol_v3000 = format!("{}> <ID>\nV3000-1\n\n", MOL_V3000);
        let id_value_v3000 = read_ffi_string(|out, cap| {
            ds_mol_block_property(
                mol_v3000.as_ptr(), mol_v3000.len(),
                id.as_ptr(), id.len(),
                out, cap,
            )
        }).unwrap();
        assert_eq!(id_value_v3000, "V3000-1");
    }

    #[test]
    fn test_ffi_sdf_property_and_json() {
        let sdf = format!("{}> <ID>\n1\n\n$$$$\n{}> <ID>\n2\n\n$$$$\n", MOL, MOL_V3000);
        let key = "ID";

        let second = read_ffi_string(|out, cap| {
            ds_sdf_property(
                sdf.as_ptr(), sdf.len(), 2,
                key.as_ptr(), key.len(),
                out, cap,
            )
        }).unwrap();
        assert_eq!(second, "2");

        let json = read_ffi_string(|out, cap| {
            ds_sdf_properties_json(sdf.as_ptr(), sdf.len(), out, cap)
        }).unwrap();
        assert_eq!(
            json,
            "[{\"record\":1,\"name\":\"ethanol\",\"properties\":[{\"name\":\"ID\",\"value\":\"1\"}]},{\"record\":2,\"name\":\"ethanol_v3000\",\"properties\":[{\"name\":\"ID\",\"value\":\"2\"}]}]"
        );
    }
}
