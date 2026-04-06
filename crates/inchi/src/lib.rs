pub mod parser;

use std::slice;

fn as_bytes(ptr: *const u8, len: usize) -> &'static [u8] {
    unsafe { slice::from_raw_parts(ptr, len) }
}

fn write_buf(src: &[u8], out: *mut u8, cap: usize) -> i32 {
    let n = src.len().min(cap);
    unsafe { std::ptr::copy_nonoverlapping(src.as_ptr(), out, n); }
    n as i32
}

// =============================================================================
// C FFI exports
// =============================================================================

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_is_valid(ptr: *const u8, len: usize) -> i32 {
    if parser::is_valid(as_bytes(ptr, len)) { 1 } else { 0 }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_is_standard(ptr: *const u8, len: usize) -> i32 {
    if !parser::is_valid(as_bytes(ptr, len)) { return -1; }
    if parser::is_standard(as_bytes(ptr, len)) { 1 } else { 0 }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_version(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::version(as_bytes(ptr, len)) {
        Some(v) => write_buf(v, out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_formula(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::formula(as_bytes(ptr, len)) {
        Some(v) => write_buf(v, out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_connections(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::layer(as_bytes(ptr, len), b'c') {
        Some(v) => write_buf(v, out, cap),
        None => 0, // absent layer = empty, not error
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_hydrogens(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::layer(as_bytes(ptr, len), b'h') {
        Some(v) => write_buf(v, out, cap),
        None => 0,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_charge(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::layer(as_bytes(ptr, len), b'q') {
        Some(v) => write_buf(v, out, cap),
        None => 0,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_stereo_bond(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::layer(as_bytes(ptr, len), b'b') {
        Some(v) => write_buf(v, out, cap),
        None => 0,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_stereo_tetrahedral(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::layer(as_bytes(ptr, len), b't') {
        Some(v) => write_buf(v, out, cap),
        None => 0,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_has_stereo(ptr: *const u8, len: usize) -> i32 {
    let s = as_bytes(ptr, len);
    if !parser::is_valid(s) { return -1; }
    if parser::has_stereo(s) { 1 } else { 0 }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_num_stereo_centers(ptr: *const u8, len: usize) -> i32 {
    let s = as_bytes(ptr, len);
    if !parser::is_valid(s) { return -1; }
    parser::num_stereo_centers(s)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchi_skeleton_match(
    a_ptr: *const u8, a_len: usize,
    b_ptr: *const u8, b_len: usize,
) -> i32 {
    let a = as_bytes(a_ptr, a_len);
    let b = as_bytes(b_ptr, b_len);
    if !parser::is_valid(a) || !parser::is_valid(b) { return -1; }
    if parser::skeleton_match(a, b) { 1 } else { 0 }
}

// ---- InChIKey ----

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchikey_is_valid(ptr: *const u8, len: usize) -> i32 {
    if parser::inchikey_is_valid(as_bytes(ptr, len)) { 1 } else { 0 }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchikey_connectivity(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::inchikey_connectivity(as_bytes(ptr, len)) {
        Some(v) => write_buf(v, out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchikey_stereo(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::inchikey_stereo(as_bytes(ptr, len)) {
        Some(v) => write_buf(v, out, cap),
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_inchikey_protonation(ptr: *const u8, len: usize, out: *mut u8, cap: usize) -> i32 {
    match parser::inchikey_protonation(as_bytes(ptr, len)) {
        Some(v) => write_buf(v, out, cap),
        None => -1,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ACETIC: &[u8] = b"InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)";
    const KEY: &[u8] = b"QTBSBXVTEAMEQO-UHFFFAOYSA-N";

    #[test]
    fn test_ffi_valid() {
        assert_eq!(ds_inchi_is_valid(ACETIC.as_ptr(), ACETIC.len()), 1);
        let bad = b"garbage";
        assert_eq!(ds_inchi_is_valid(bad.as_ptr(), bad.len()), 0);
    }

    #[test]
    fn test_ffi_formula() {
        let mut buf = [0u8; 64];
        let len = ds_inchi_formula(ACETIC.as_ptr(), ACETIC.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(&buf[..len as usize], b"C2H4O2");
    }

    #[test]
    fn test_ffi_inchikey() {
        assert_eq!(ds_inchikey_is_valid(KEY.as_ptr(), KEY.len()), 1);
        let mut buf = [0u8; 64];
        let len = ds_inchikey_connectivity(KEY.as_ptr(), KEY.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(&buf[..len as usize], b"QTBSBXVTEAMEQO");
    }
}
