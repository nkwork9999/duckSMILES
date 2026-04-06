pub mod encoder;
pub mod decoder;

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

/// SMILES → SELFIES. Returns length written, or -1 on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ds_smiles_to_selfies(
    ptr: *const u8, len: usize,
    out: *mut u8, cap: usize,
) -> i32 {
    let smiles = as_str(ptr, len);
    match encoder::smiles_to_selfies(smiles) {
        Some(s) => write_buf(s.as_bytes(), out, cap),
        None => -1,
    }
}

/// SELFIES → SMILES. Returns length written, or -1 on failure.
#[unsafe(no_mangle)]
pub extern "C" fn ds_selfies_to_smiles(
    ptr: *const u8, len: usize,
    out: *mut u8, cap: usize,
) -> i32 {
    let selfies = as_str(ptr, len);
    match decoder::selfies_to_smiles(selfies) {
        Some(s) => write_buf(s.as_bytes(), out, cap),
        None => -1,
    }
}

/// Validate SELFIES (any non-empty SELFIES with brackets is valid by design)
#[unsafe(no_mangle)]
pub extern "C" fn ds_selfies_is_valid(ptr: *const u8, len: usize) -> i32 {
    let s = as_str(ptr, len);
    if s.is_empty() { return 0; }
    // SELFIES is valid if it contains at least one token
    if s.contains('[') && s.contains(']') { 1 } else { 0 }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_roundtrip_simple() {
        let smiles = "CCO";
        let selfies = encoder::smiles_to_selfies(smiles).unwrap();
        let back = decoder::selfies_to_smiles(&selfies).unwrap();
        assert_eq!(back, "CCO");
    }

    #[test]
    fn test_roundtrip_double_bond() {
        let smiles = "C=O";
        let selfies = encoder::smiles_to_selfies(smiles).unwrap();
        let back = decoder::selfies_to_smiles(&selfies).unwrap();
        assert_eq!(back, "C=O");
    }

    #[test]
    fn test_roundtrip_triple() {
        let smiles = "C#N";
        let selfies = encoder::smiles_to_selfies(smiles).unwrap();
        let back = decoder::selfies_to_smiles(&selfies).unwrap();
        assert_eq!(back, "C#N");
    }

    #[test]
    fn test_ffi_encode() {
        let smiles = b"CCO";
        let mut buf = [0u8; 256];
        let len = ds_smiles_to_selfies(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert!(len > 0);
        assert_eq!(&buf[..len as usize], b"[C][C][O]");
    }

    #[test]
    fn test_ffi_decode() {
        let selfies = b"[C][C][O]";
        let mut buf = [0u8; 256];
        let len = ds_selfies_to_smiles(selfies.as_ptr(), selfies.len(), buf.as_mut_ptr(), buf.len());
        assert!(len > 0);
        assert_eq!(&buf[..len as usize], b"CCO");
    }

    #[test]
    fn test_ffi_valid() {
        assert_eq!(ds_selfies_is_valid(b"[C][O]".as_ptr(), 6), 1);
        assert_eq!(ds_selfies_is_valid(b"".as_ptr(), 0), 0);
    }
}
