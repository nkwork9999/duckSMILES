//! Composition and solvent feature profiles; no learned property predictions.
use std::collections::BTreeMap;

/// Positive fractional stoichiometry, element symbols and parenthesized groups.
/// No isotope, charge, hydrate-dot, variable occupancy or phase inference.
pub fn composition(formula: &str) -> Option<BTreeMap<String, f64>> {
    if formula.is_empty() || formula.len() > 4096 || !formula.is_ascii() {
        return None;
    }
    fn amount(b: &[u8], i: &mut usize) -> Option<f64> {
        let start = *i;
        while *i < b.len() && b[*i].is_ascii_digit() {
            *i += 1;
        }
        if *i < b.len() && b[*i] == b'.' {
            *i += 1;
            let after = *i;
            while *i < b.len() && b[*i].is_ascii_digit() {
                *i += 1;
            }
            if after == *i {
                return None;
            }
        }
        if start == *i {
            return Some(1.);
        }
        let n: f64 = std::str::from_utf8(&b[start..*i]).ok()?.parse().ok()?;
        (n.is_finite() && n > 0. && n <= 1e6).then_some(n)
    }
    fn group(b: &[u8], i: &mut usize, depth: usize) -> Option<BTreeMap<String, f64>> {
        if depth > 16 {
            return None;
        }
        let mut out = BTreeMap::new();
        while *i < b.len() && b[*i] != b')' {
            let terms = if b[*i] == b'(' {
                *i += 1;
                let terms = group(b, i, depth + 1)?;
                if *i >= b.len() || b[*i] != b')' {
                    return None;
                }
                *i += 1;
                terms
            } else {
                let start = *i;
                if !b[*i].is_ascii_uppercase() {
                    return None;
                }
                *i += 1;
                if *i < b.len() && b[*i].is_ascii_lowercase() {
                    *i += 1;
                }
                let symbol = std::str::from_utf8(&b[start..*i]).ok()?;
                if !crate::periodic_data::ELEMENTS.iter().any(|e| e.0 == symbol) {
                    return None;
                }
                BTreeMap::from([(symbol.to_owned(), 1.)])
            };
            let scale = amount(b, i)?;
            for (e, n) in terms {
                *out.entry(e).or_insert(0.) += n * scale;
            }
            if out.values().sum::<f64>() > 1e6 {
                return None;
            }
        }
        (!out.is_empty()).then_some(out)
    }
    let mut i = 0;
    let out = group(formula.as_bytes(), &mut i, 0)?;
    (i == formula.len()).then_some(out)
}

pub fn composition_json(formula: &str) -> Option<String> {
    let counts = composition(formula)?;
    let total = counts.values().sum::<f64>();
    let mut weight = 0.;
    let mut zmean = 0.;
    let mut entropy = 0.;
    let mut rows = Vec::new();
    for (e, n) in &counts {
        let z = crate::periodic_data::ELEMENTS
            .iter()
            .position(|v| v.0 == e)?
            + 1;
        let mass = crate::weights::atomic_weight(e)?;
        let fraction = n / total;
        weight += n * mass;
        zmean += fraction * z as f64;
        entropy -= fraction * fraction.ln();
        rows.push(format!(
            "{{\"element\":\"{e}\",\"amount\":{n},\"atomic_fraction\":{fraction}}}"
        ));
    }
    let zstd = counts
        .iter()
        .map(|(e, n)| {
            let z = crate::periodic_data::ELEMENTS
                .iter()
                .position(|v| v.0 == e)
                .unwrap()
                + 1;
            n / total * (z as f64 - zmean).powi(2)
        })
        .sum::<f64>()
        .sqrt();
    let norms: Vec<String> = [2, 3, 5, 7, 10]
        .iter()
        .map(|p| {
            counts
                .values()
                .map(|n| (n / total).powi(*p))
                .sum::<f64>()
                .powf(1. / *p as f64)
                .to_string()
        })
        .collect();
    Some(format!("{{\"schema\":1,\"method\":\"composition-fractions-v1\",\"weight_reference\":\"RDKit-2026.03.6\",\"elements\":[{}],\"n_elements\":{},\"formula_amount\":{total},\"formula_weight_g_mol\":{weight},\"mean_atomic_weight\":{},\"mean_atomic_number\":{zmean},\"std_atomic_number\":{zstd},\"mixing_entropy_over_R\":{entropy},\"stoichiometry_p\":[2,3,5,7,10],\"stoichiometry_norms\":[{}]}}",rows.join(","),counts.len(),weight/total,norms.join(",")))
}

pub const SOLVENT_IDS: &[&str] = &[
    "MolWt",
    "MolLogP",
    "MolMR",
    "TPSA",
    "NumHDonors",
    "NumHAcceptors",
    "NumRotatableBonds",
    "HeavyAtomCount",
    "FractionCSP3",
    "MaxPartialCharge",
    "MinPartialCharge",
];
pub fn solvent_json(smiles: &str) -> Option<String> {
    // Finite molecular features, not dielectric constant, donor number or viscosity.
    let d = crate::workbench::analyze(smiles)?;
    if d.heavy_atoms == 0 || d.fragments != 1 {
        return None;
    }
    let mut features = Vec::new();
    for id in SOLVENT_IDS {
        let index = crate::descriptor_catalog::SPECS
            .iter()
            .position(|s| s.id == *id)?;
        let v = d.descriptor_values[index]
            .filter(|v| v.is_finite())
            .map_or_else(|| "null".into(), |v| v.to_string());
        features.push(format!("\"{id}\":{v}"));
    }
    Some(format!("{{\"schema\":1,\"method\":\"solvent-molecular-profile-v1\",\"descriptor_reference\":\"RDKit-2026.03.6\",\"formal_charge\":{},\"descriptors\":{{{}}}}}",d.formal_charge,features.join(",")))
}

unsafe fn json_ffi(
    ptr: *const u8,
    len: usize,
    out: *mut u8,
    cap: usize,
    f: fn(&str) -> Option<String>,
) -> i32 {
    std::panic::catch_unwind(|| {
        if ptr.is_null() || len == 0 || len > 16384 || (out.is_null() && cap != 0) {
            return -1;
        }
        // SAFETY: caller guarantees readable input and disjoint writable output.
        let Ok(s) = std::str::from_utf8(unsafe { std::slice::from_raw_parts(ptr, len) }) else {
            return -1;
        };
        let Some(json) = f(s) else {
            return -1;
        };
        let Ok(n) = i32::try_from(json.len()) else {
            return -1;
        };
        if !out.is_null() && cap >= json.len() {
            // SAFETY: capacity was checked; input and output do not overlap.
            unsafe {
                std::ptr::copy_nonoverlapping(json.as_ptr(), out, json.len());
            }
        }
        n
    })
    .unwrap_or(-1)
}
/// # Safety
/// Input/output pointers must cover their declared lengths, and must not overlap.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_composition_descriptors_json(
    ptr: *const u8,
    len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    unsafe { json_ffi(ptr, len, out, cap, composition_json) }
}
/// # Safety
/// Input/output pointers must cover their declared lengths, and must not overlap.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_solvent_descriptors_json(
    ptr: *const u8,
    len: usize,
    out: *mut u8,
    cap: usize,
) -> i32 {
    unsafe { json_ffi(ptr, len, out, cap, solvent_json) }
}
