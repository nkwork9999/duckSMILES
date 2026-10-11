//! Formula and adduct arithmetic for a neutral molecule M.
//! No structure identity, adduct inference, or spectrum modification is performed.

pub const SPEC_VERSION: &str = "duck-formula-adduct-v1";
pub const MAX_FORMULA_BYTES: usize = 4096;
const MAX_ELEMENT_COUNT: u32 = 1_000_000;
const ELECTRON: f64 = 0.00054857990;
const H: f64 = 1.00782503207;
const N: f64 = 14.0030740048;
const O: f64 = 15.9949146196;
const NA: f64 = 22.9897692809;
const K: f64 = 38.96370668;
const WATER: f64 = 2.0 * H + O;

// Versioned numerical table from the public CASMI AFIX formula-mass routine.
// Keep this separate from the legacy SMILES mass API, whose table differs.
fn atom_mass(symbol: &str) -> Option<f64> {
    Some(match symbol {
        "C" => 12.0,
        "H" => H,
        "N" => N,
        "O" => O,
        "S" => 31.97207100,
        "P" => 30.97376163,
        "Cl" => 34.96885268,
        "Br" => 78.9183371,
        "F" => 18.99840322,
        "I" => 126.904473,
        "Si" => 27.9769265325,
        "Se" => 79.9165213,
        "Na" => NA,
        "K" => K,
        "B" => 11.0093054,
        "As" => 74.9215965,
        "Fe" => 55.9349375,
        "D" => 2.0141017778,
        _ => return None,
    })
}

/// Sum neutral atomic monoisotopic masses from a complete, flat formula.
///
/// Tokens are `[A-Z][a-z]?[1-9][0-9]*` with an optional count (default 1).
/// Repeated symbols are summed. Parentheses, charges, isotopic brackets,
/// whitespace, mixtures and unknown elements are rejected, not partly parsed.
pub fn formula_exact_mass(formula: &str) -> Option<f64> {
    if formula.is_empty() || formula.len() > MAX_FORMULA_BYTES || !formula.is_ascii() {
        return None;
    }
    let b = formula.as_bytes();
    let mut pos = 0;
    let mut total = 0.0;
    let mut atoms = 0u32;
    while pos < b.len() {
        let start = pos;
        if !b[pos].is_ascii_uppercase() {
            return None;
        }
        pos += 1;
        if pos < b.len() && b[pos].is_ascii_lowercase() {
            pos += 1;
        }
        let mass = atom_mass(&formula[start..pos])?;
        let mut count = 0u32;
        if pos < b.len() && b[pos].is_ascii_digit() {
            if b[pos] == b'0' {
                return None;
            }
            while pos < b.len() && b[pos].is_ascii_digit() {
                count = count
                    .checked_mul(10)?
                    .checked_add(u32::from(b[pos] - b'0'))?;
                pos += 1;
            }
        } else {
            count = 1;
        }
        atoms = atoms.checked_add(count)?;
        if atoms > MAX_ELEMENT_COUNT {
            return None;
        }
        total += mass * f64::from(count);
    }
    Some(total)
}

/// Explicitly supported labels. No aliases or missing-charge guesses are made.
pub const SUPPORTED_ADDUCTS: &[&str] = &[
    "[M+H]+",
    "[M+NH4]+",
    "[M+Na]+",
    "[M+K]+",
    "[M-H2O+H]+",
    "[M-2H2O+H]+",
    "[M+2H]2+",
    "[M]+",
    "[M-H2O]+",
    "[M+CH3OH+H]+",
    "[M+CH3CN+H]+",
    "[M-H]-",
    "[M-H2O-H]-",
    "[M+CH2O2-H]-",
    "[M+C2H4O2-H]-",
    "[M+Cl]-",
    "[M]-",
    "[M-2H]2-",
    "[M+Na-2H]-",
    "[2M+H]+",
    "[2M+Na]+",
    "[2M+NH4]+",
    "[2M+K]+",
    "[2M-H]-",
    "[2M+CH2O2-H]-",
    "[2M+C2H4O2-H]-",
    "[2M+Na-2H]-",
    "[3M+H]+",
    "[3M-H]-",
];

// (molecule count, signed charge, sum of neutral atomic additions/losses).
fn adduct_parameters(adduct: &str) -> Option<(f64, f64, f64)> {
    let (n, z, atoms) = match adduct {
        "[M+H]+" => (1.0, 1.0, H),
        "[M+NH4]+" => (1.0, 1.0, N + 4.0 * H),
        "[M+Na]+" => (1.0, 1.0, NA),
        "[M+K]+" => (1.0, 1.0, K),
        "[M-H2O+H]+" => (1.0, 1.0, H - WATER),
        "[M-2H2O+H]+" => (1.0, 1.0, H - 2.0 * WATER),
        "[M+2H]2+" => (1.0, 2.0, 2.0 * H),
        "[M]+" => (1.0, 1.0, 0.0),
        "[M-H2O]+" => (1.0, 1.0, -WATER),
        "[M+CH3OH+H]+" => (1.0, 1.0, 12.0 + 5.0 * H + O),
        "[M+CH3CN+H]+" => (1.0, 1.0, 24.0 + 4.0 * H + N),
        "[M-H]-" => (1.0, -1.0, -H),
        "[M-H2O-H]-" => (1.0, -1.0, -H - WATER),
        "[M+CH2O2-H]-" => (1.0, -1.0, 12.0 + H + 2.0 * O),
        "[M+C2H4O2-H]-" => (1.0, -1.0, 24.0 + 3.0 * H + 2.0 * O),
        "[M+Cl]-" => (1.0, -1.0, 34.96885268),
        "[M]-" => (1.0, -1.0, 0.0),
        "[M-2H]2-" => (1.0, -2.0, -2.0 * H),
        "[M+Na-2H]-" => (1.0, -1.0, NA - 2.0 * H),
        "[2M+H]+" => (2.0, 1.0, H),
        "[2M+Na]+" => (2.0, 1.0, NA),
        "[2M+NH4]+" => (2.0, 1.0, N + 4.0 * H),
        "[2M+K]+" => (2.0, 1.0, K),
        "[2M-H]-" => (2.0, -1.0, -H),
        "[2M+CH2O2-H]-" => (2.0, -1.0, 12.0 + H + 2.0 * O),
        "[2M+C2H4O2-H]-" => (2.0, -1.0, 24.0 + 3.0 * H + 2.0 * O),
        "[2M+Na-2H]-" => (2.0, -1.0, NA - 2.0 * H),
        "[3M+H]+" => (3.0, 1.0, H),
        "[3M-H]-" => (3.0, -1.0, -H),
        _ => return None,
    };
    Some((n, z, atoms - z * ELECTRON))
}

fn positive_finite(x: f64) -> Option<f64> {
    (x.is_finite() && x > 0.0).then_some(x)
}

/// Invert a supported adduct: M = (m/z * abs(z) - delta) / n.
pub fn adduct_neutral_mass(mz: f64, adduct: &str) -> Option<f64> {
    positive_finite(mz)?;
    let (n, z, delta) = adduct_parameters(adduct)?;
    positive_finite((mz * z.abs() - delta) / n)
}

/// Apply a supported adduct: m/z = (n * M + delta) / abs(z).
/// This is arithmetic, not evidence that this ion exists for the molecule.
pub fn adduct_precursor_mz(neutral_mass: f64, adduct: &str) -> Option<f64> {
    positive_finite(neutral_mass)?;
    let (n, z, delta) = adduct_parameters(adduct)?;
    positive_finite((n * neutral_mass + delta) / z.abs())
}

/// Caller must provide a readable allocation of `len` bytes for a non-null ptr.
unsafe fn input_text<'a>(ptr: *const u8, len: usize) -> Option<&'a str> {
    if ptr.is_null() || len == 0 || len > MAX_FORMULA_BYTES {
        return None;
    }
    // SAFETY: the FFI caller guarantees readability; length is bounded above.
    std::str::from_utf8(unsafe { std::slice::from_raw_parts(ptr, len) }).ok()
}

/// Formula mass or NaN for invalid/unsupported input.
/// # Safety
/// A non-null `ptr` must address at least `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_formula_exact_mass(ptr: *const u8, len: usize) -> f64 {
    // SAFETY: inherited from this function's caller contract.
    unsafe { input_text(ptr, len) }
        .and_then(formula_exact_mass)
        .unwrap_or(f64::NAN)
}

/// Neutral mass or NaN for invalid/unsupported input.
/// # Safety
/// A non-null `ptr` must address at least `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_adduct_neutral_mass(mz: f64, ptr: *const u8, len: usize) -> f64 {
    // SAFETY: inherited from this function's caller contract.
    unsafe { input_text(ptr, len) }
        .and_then(|a| adduct_neutral_mass(mz, a))
        .unwrap_or(f64::NAN)
}

/// Precursor m/z or NaN for invalid/unsupported input.
/// # Safety
/// A non-null `ptr` must address at least `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_adduct_precursor_mz(mass: f64, ptr: *const u8, len: usize) -> f64 {
    // SAFETY: inherited from this function's caller contract.
    unsafe { input_text(ptr, len) }
        .and_then(|a| adduct_precursor_mz(mass, a))
        .unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formula_examples_and_nitrogen() {
        assert!((formula_exact_mass("C8H10N4O2").unwrap() - 194.0803755791).abs() < 1e-10);
        assert!((formula_exact_mass("C5H5N").unwrap() - 79.04219916515).abs() < 1e-10);
        assert_eq!(formula_exact_mass("CH3COOH"), formula_exact_mass("C2H4O2"));
        assert!((formula_exact_mass("D2O").unwrap() - 20.0231181752).abs() < 1e-10);
        assert_eq!(formula_exact_mass("C1000000"), Some(12_000_000.0));
    }

    #[test]
    fn reject_partial_formula_or_overflow() {
        for s in [
            "",
            " ",
            "C0",
            "C01",
            "C-1",
            "C1.5",
            "(CH3)2",
            "CuSO4.5H2O",
            "NH4+",
            "[13C]H4",
            "13CH4",
            "C H",
            "C\0H",
            "cH",
            "Cé",
            "Xx",
            "He",
            "C1000001",
            "C9999999999999999999",
            "C999999H2",
        ] {
            assert!(formula_exact_mass(s).is_none(), "{s:?}");
        }
        assert!(formula_exact_mass(&"C".repeat(MAX_FORMULA_BYTES + 1)).is_none());
    }

    #[test]
    fn every_supported_adduct_roundtrips() {
        for a in SUPPORTED_ADDUCTS {
            for m in [50.0, 194.0803755791, 1000.0, 100_000.0] {
                let mz = adduct_precursor_mz(m, a).unwrap();
                assert!(
                    (adduct_neutral_mass(mz, a).unwrap() - m).abs() < 1e-9,
                    "{a}"
                );
            }
        }
    }

    #[test]
    fn charge_and_multimer_mass_balance() {
        let m = 194.0803755791;
        assert!((adduct_precursor_mz(m, "[M+H]+").unwrap() - 195.08765203127).abs() < 1e-10);
        assert!(
            (adduct_precursor_mz(m, "[M+2H]2+").unwrap() - (m / 2.0 + H - ELECTRON)).abs() < 1e-10
        );
        assert!(
            (adduct_precursor_mz(m, "[M-2H]2-").unwrap() - (m / 2.0 - H + ELECTRON)).abs() < 1e-10
        );
        assert!(
            (adduct_precursor_mz(m, "[2M+Na-2H]-").unwrap() - (2.0 * m + NA - 2.0 * H + ELECTRON))
                .abs()
                < 1e-10
        );
        assert!(
            (adduct_precursor_mz(m, "[3M-H]-").unwrap() - (3.0 * m - H + ELECTRON)).abs() < 1e-10
        );
    }

    #[test]
    fn invalid_and_impossible_mass() {
        for x in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY, 0.0, -1.0] {
            assert!(adduct_neutral_mass(x, "[M+H]+").is_none());
            assert!(adduct_precursor_mz(x, "[M+H]+").is_none());
        }
        assert!(adduct_precursor_mz(f64::MAX, "[2M+H]+").is_none());
        assert!(adduct_neutral_mass(f64::MAX, "[M+2H]2+").is_none());
        assert!(adduct_neutral_mass(1.0, "[M+Na]+").is_none());
        assert!(adduct_precursor_mz(1.0, "[M-H]-").is_none());
        for a in [
            "", "[M-2H]-", "[M+H]", "M+H", "[M+H]+ ", "[M+H]0+", "[M+Li]+",
        ] {
            assert!(adduct_precursor_mz(200.0, a).is_none());
            assert!(adduct_neutral_mass(200.0, a).is_none());
        }
    }

    #[test]
    fn ffi_valid_null_utf8_and_length() {
        // SAFETY: every non-null input below is backed by a readable byte slice.
        unsafe {
            assert!(ds_formula_exact_mass(std::ptr::null(), 0).is_nan());
            assert!(ds_formula_exact_mass(std::ptr::null(), 10).is_nan());
            assert!(ds_formula_exact_mass(b"C".as_ptr(), usize::MAX).is_nan());
            assert!(ds_formula_exact_mass([0xffu8].as_ptr(), 1).is_nan());
            assert_eq!(ds_formula_exact_mass(b"C".as_ptr(), 1), 12.0);
            let a = "[M+H]+";
            let mz = ds_adduct_precursor_mz(100.0, a.as_ptr(), a.len());
            assert!((ds_adduct_neutral_mass(mz, a.as_ptr(), a.len()) - 100.0).abs() < 1e-12);
            assert!(ds_adduct_neutral_mass(100.0, std::ptr::null(), 1).is_nan());
            assert!(ds_adduct_precursor_mz(100.0, [0xffu8].as_ptr(), 1).is_nan());
        }
    }
}
