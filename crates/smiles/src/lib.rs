mod logp_crippen;
mod maccs;
mod morgan;
mod parser;
mod smarts;
mod tanimoto;
#[cfg(test)]
pub(crate) mod test_fixtures;
mod tpsa;
mod weights;

use logp_crippen::calc_logp;
use morgan::morgan_bits;
use parser::{parse, BondOrder, Molecule};
use smarts::{count_unique, matches_mol, parse_smarts};
use tanimoto::tanimoto_bit;
use tpsa::calc_tpsa;

/// Re-exports used by the `maccs_verify` example (Python-RDKit cross-check).
/// Not part of the C ABI; kept minimal.
pub mod verify {
    pub use crate::maccs::{maccs_bits, on_bits, MACCS_N_BYTES};
    pub use crate::parser::parse;
}

fn is_hetero_atom(symbol: &str) -> bool {
    symbol != "C" && symbol != "H"
}

fn is_h_acceptor_atom(mol: &Molecule, idx: usize) -> bool {
    let atom = &mol.atoms[idx];
    if atom.charge > 0 {
        return false;
    }
    match atom.symbol.as_str() {
        "O" | "S" => true,
        "N" => {
            if atom.aromatic && atom.hydrogen > 0 {
                return false;
            }
            // Amide-like nitrogens are poor acceptors:
            // N-C(=O), where the carbonyl carbon has a double-bonded oxygen.
            for (nbr, order) in mol.neighbors(idx) {
                if order != BondOrder::Single || mol.atoms[nbr].symbol != "C" {
                    continue;
                }
                let has_carbonyl_o = mol.neighbors(nbr).into_iter().any(|(nbr2, order2)| {
                    nbr2 != idx && order2 == BondOrder::Double && mol.atoms[nbr2].symbol == "O"
                });
                if has_carbonyl_o {
                    return false;
                }
            }
            true
        }
        _ => false,
    }
}

fn num_h_acceptors(mol: &Molecule) -> usize {
    (0..mol.atoms.len())
        .filter(|&idx| is_h_acceptor_atom(mol, idx))
        .count()
}

fn num_h_donors(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|atom| matches!(atom.symbol.as_str(), "N" | "O" | "S") && atom.hydrogen > 0)
        .count()
}

fn num_heteroatoms(mol: &Molecule) -> usize {
    mol.atoms
        .iter()
        .filter(|atom| is_hetero_atom(&atom.symbol))
        .count()
}

fn ring_count(mol: &Molecule) -> usize {
    mol.ring_info().rings.len()
}

fn num_aromatic_rings(mol: &Molecule) -> usize {
    let ring_info = mol.ring_info();
    ring_info
        .rings
        .iter()
        .filter(|ring| {
            ring.iter().all(|&bond_idx| {
                mol.bonds
                    .get(bond_idx)
                    .map(|bond| bond.order == BondOrder::Aromatic)
                    .unwrap_or(false)
            })
        })
        .count()
}

fn is_terminal_heavy_atom(mol: &Molecule, idx: usize) -> bool {
    mol.neighbors(idx)
        .iter()
        .filter(|(nbr, _)| mol.atoms[*nbr].symbol != "H")
        .count()
        <= 1
}

fn is_amide_like_bond(mol: &Molecule, a: usize, b: usize) -> bool {
    let (n_idx, c_idx) = match (mol.atoms[a].symbol.as_str(), mol.atoms[b].symbol.as_str()) {
        ("N", "C") => (a, b),
        ("C", "N") => (b, a),
        _ => return false,
    };
    mol.neighbors(c_idx).into_iter().any(|(nbr, order)| {
        nbr != n_idx && order == BondOrder::Double && mol.atoms[nbr].symbol == "O"
    })
}

fn num_rotatable_bonds(mol: &Molecule) -> usize {
    let ring_info = mol.ring_info();
    mol.bonds
        .iter()
        .enumerate()
        .filter(|(idx, bond)| {
            bond.order == BondOrder::Single
                && !ring_info.bond_in_ring.get(*idx).copied().unwrap_or(false)
                && mol.atoms[bond.a].symbol != "H"
                && mol.atoms[bond.b].symbol != "H"
                && !is_terminal_heavy_atom(mol, bond.a)
                && !is_terminal_heavy_atom(mol, bond.b)
                && !is_amide_like_bond(mol, bond.a, bond.b)
        })
        .count()
}

fn fraction_csp3(mol: &Molecule) -> f64 {
    let carbon_indices = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, atom)| atom.symbol == "C")
        .map(|(idx, _)| idx)
        .collect::<Vec<_>>();
    if carbon_indices.is_empty() {
        return 0.0;
    }
    let sp3 = carbon_indices
        .iter()
        .filter(|&&idx| {
            let atom = &mol.atoms[idx];
            !atom.aromatic
                && mol
                    .neighbors(idx)
                    .iter()
                    .all(|(_, order)| *order == BondOrder::Single)
        })
        .count();
    sp3 as f64 / carbon_indices.len() as f64
}

// =============================================================================
// C FFI exports
// =============================================================================

/// Returns 1 if valid SMILES, 0 otherwise
#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_is_valid(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    if parse(s).is_some() {
        1
    } else {
        0
    }
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
pub extern "C" fn ds_mol_formula(ptr: *const u8, len: usize, out: *mut u8, out_cap: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let formula = mol.formula();
            let bytes = formula.as_bytes();
            let n = bytes.len().min(out_cap);
            unsafe {
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), out, n);
            }
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

/// Returns RDKit-default TPSA (N/O only), or NaN on invalid SMILES.
#[unsafe(no_mangle)]
pub extern "C" fn ds_tpsa(ptr: *const u8, len: usize) -> f64 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => calc_tpsa(&mol),
        None => f64::NAN,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_canonical_smiles(
    ptr: *const u8,
    len: usize,
    out: *mut u8,
    out_cap: usize,
) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let canonical = mol.canonical_smiles();
            let bytes = canonical.as_bytes();
            if bytes.len() > out_cap {
                return -1;
            }
            unsafe {
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), out, bytes.len());
            }
            bytes.len() as i32
        }
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_num_h_acceptors(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => num_h_acceptors(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_num_h_donors(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => num_h_donors(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_num_rotatable_bonds(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => num_rotatable_bonds(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_ring_count(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => ring_count(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_num_aromatic_rings(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => num_aromatic_rings(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_num_heteroatoms(ptr: *const u8, len: usize) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => num_heteroatoms(&mol) as i32,
        None => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_fraction_csp3(ptr: *const u8, len: usize) -> f64 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => fraction_csp3(&mol),
        None => f64::NAN,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_has_substructure(
    smiles_ptr: *const u8,
    smiles_len: usize,
    smarts_ptr: *const u8,
    smarts_len: usize,
) -> i32 {
    let smiles = unsafe {
        std::str::from_utf8_unchecked(std::slice::from_raw_parts(smiles_ptr, smiles_len))
    };
    let smarts = unsafe {
        std::str::from_utf8_unchecked(std::slice::from_raw_parts(smarts_ptr, smarts_len))
    };
    match (parse(smiles), parse_smarts(smarts)) {
        (Some(mol), Some(pattern)) => {
            if matches_mol(&pattern, &mol) {
                1
            } else {
                0
            }
        }
        _ => -1,
    }
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_mol_substructure_count(
    smiles_ptr: *const u8,
    smiles_len: usize,
    smarts_ptr: *const u8,
    smarts_len: usize,
) -> i32 {
    let smiles = unsafe {
        std::str::from_utf8_unchecked(std::slice::from_raw_parts(smiles_ptr, smiles_len))
    };
    let smarts = unsafe {
        std::str::from_utf8_unchecked(std::slice::from_raw_parts(smarts_ptr, smarts_len))
    };
    match (parse(smiles), parse_smarts(smarts)) {
        (Some(mol), Some(pattern)) => count_unique(&pattern, &mol) as i32,
        _ => -1,
    }
}

/// Writes SMILES with explicit H atoms to buffer. Returns length written, or -1 on invalid.
/// Result is a verbose bracket-form SMILES that round-trips through parse.
#[unsafe(no_mangle)]
pub extern "C" fn ds_add_hydrogens(
    ptr: *const u8,
    len: usize,
    out: *mut u8,
    out_cap: usize,
) -> i32 {
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let smi = mol.with_explicit_hydrogens().to_smiles_verbose();
            let bytes = smi.as_bytes();
            let n = bytes.len().min(out_cap);
            unsafe {
                std::ptr::copy_nonoverlapping(bytes.as_ptr(), out, n);
            }
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
    a_ptr: *const u8,
    a_len: usize,
    b_ptr: *const u8,
    b_len: usize,
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
    ptr: *const u8,
    len: usize,
    radius: u32,
    n_bits: u32,
    out: *mut u8,
    out_cap: usize,
) -> i32 {
    if n_bits == 0 {
        return -1;
    }
    let n_bytes = ((n_bits as usize) + 7) / 8;
    if n_bytes > out_cap {
        return -1;
    }
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let bits = morgan_bits(&mol, radius, n_bits);
            unsafe {
                std::ptr::copy_nonoverlapping(bits.as_ptr(), out, n_bytes);
            }
            n_bytes as i32
        }
        None => -1,
    }
}

/// Writes the 166 MACCS keys to `out` as a fixed 21-byte (167-bit) buffer.
/// Returns bytes written (always 21), or -1 on invalid SMILES / buffer too small.
///
/// Bit `n` (1..=166) is stored at `byte = n / 8`, `offset = n % 8`. Bit 0 and
/// bit 1 (isotope) are always 0, matching RDKit's `ExplicitBitVect(167)`.
#[unsafe(no_mangle)]
pub extern "C" fn ds_maccs_keys(ptr: *const u8, len: usize, out: *mut u8, out_cap: usize) -> i32 {
    if out_cap < maccs::MACCS_N_BYTES {
        return -1;
    }
    let s = unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) };
    match parse(s) {
        Some(mol) => {
            let bits = maccs::maccs_bits(&mol);
            unsafe {
                std::ptr::copy_nonoverlapping(bits.as_ptr(), out, maccs::MACCS_N_BYTES);
            }
            maccs::MACCS_N_BYTES as i32
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

    fn assert_close(label: &str, got: f64, expected: f64, tolerance: f64) {
        assert!(
            (got - expected).abs() <= tolerance,
            "{label}: got {got:.6}, expected {expected:.6} ± {tolerance}"
        );
    }

    fn logp_ffi(smiles: &str) -> f64 {
        ds_logp_crippen(smiles.as_ptr(), smiles.len())
    }

    fn tpsa_ffi(smiles: &str) -> f64 {
        ds_tpsa(smiles.as_ptr(), smiles.len())
    }

    fn canonical_ffi(smiles: &str) -> String {
        let mut buf = [0u8; 1024];
        let n = ds_canonical_smiles(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert!(n >= 0, "canonical_smiles failed for {smiles}");
        String::from_utf8_lossy(&buf[..n as usize]).to_string()
    }

    fn maccs_ffi(smiles: &str) -> [u8; 21] {
        let mut buf = [0u8; 21];
        let n = ds_maccs_keys(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(n, 21, "MACCS generation failed for {smiles}");
        buf
    }

    fn on_bits(blob: &[u8]) -> Vec<usize> {
        let mut bits = Vec::new();
        for bit in 0..blob.len() * 8 {
            if (blob[bit / 8] & (1 << (bit % 8))) != 0 {
                bits.push(bit);
            }
        }
        bits
    }

    fn tanimoto_ffi(a: &[u8], b: &[u8]) -> f64 {
        ds_tanimoto_bit(a.as_ptr(), a.len(), b.as_ptr(), b.len())
    }

    #[test]
    fn lipinski_descriptor_ffi_basics() {
        let ethanol = "CCO";
        assert_eq!(ds_num_h_acceptors(ethanol.as_ptr(), ethanol.len()), 1);
        assert_eq!(ds_num_h_donors(ethanol.as_ptr(), ethanol.len()), 1);
        assert_eq!(ds_num_rotatable_bonds(ethanol.as_ptr(), ethanol.len()), 0);
        assert_eq!(ds_ring_count(ethanol.as_ptr(), ethanol.len()), 0);
        assert_eq!(ds_num_aromatic_rings(ethanol.as_ptr(), ethanol.len()), 0);
        assert_eq!(ds_num_heteroatoms(ethanol.as_ptr(), ethanol.len()), 1);
        assert_close(
            "ethanol fraction_csp3",
            ds_fraction_csp3(ethanol.as_ptr(), ethanol.len()),
            1.0,
            1e-12,
        );

        let benzene = "c1ccccc1";
        assert_eq!(ds_ring_count(benzene.as_ptr(), benzene.len()), 1);
        assert_eq!(ds_num_aromatic_rings(benzene.as_ptr(), benzene.len()), 1);
        assert_close(
            "benzene fraction_csp3",
            ds_fraction_csp3(benzene.as_ptr(), benzene.len()),
            0.0,
            1e-12,
        );

        let butane = "CCCC";
        assert_eq!(ds_num_rotatable_bonds(butane.as_ptr(), butane.len()), 1);
    }

    #[test]
    fn substructure_ffi_basics() {
        let aspirin = "CC(=O)Oc1ccccc1C(=O)O";
        let carbonyl = "[#6]=[#8]";
        let aromatic_carbon = "c";
        assert_eq!(
            ds_mol_has_substructure(
                aspirin.as_ptr(),
                aspirin.len(),
                carbonyl.as_ptr(),
                carbonyl.len()
            ),
            1
        );
        assert_eq!(
            ds_mol_substructure_count(
                aspirin.as_ptr(),
                aspirin.len(),
                carbonyl.as_ptr(),
                carbonyl.len()
            ),
            2
        );
        assert_eq!(
            ds_mol_substructure_count(
                aspirin.as_ptr(),
                aspirin.len(),
                aromatic_carbon.as_ptr(),
                aromatic_carbon.len()
            ),
            6
        );
    }

    #[test]
    fn canonical_smiles_aromatizes_kekule_forms() {
        assert_eq!(canonical_ffi("C1=CC=CC=C1"), "c1ccccc1");
        assert_eq!(canonical_ffi("c1ccccc1"), "c1ccccc1");
        assert_eq!(canonical_ffi("C1CCCCC1"), "C1CCCCC1");
    }

    #[test]
    fn rdkit_compat_logp_exact_small_molecules() {
        // Reference values captured from RDKit Descriptors.MolLogP.
        // Larger drug-like molecules are tracked separately because the current
        // parser/aromaticity model intentionally does not claim full RDKit parity.
        let cases = [
            ("methane", "C", 0.6361),
            ("water", "O", -0.8247),
            ("ethanol", "CCO", -0.0014),
            ("methanol", "CO", -0.3915),
            ("benzene", "c1ccccc1", 1.6866),
        ];

        for (name, smiles, expected) in cases {
            assert_close(name, logp_ffi(smiles), expected, 0.0001);
        }
    }

    #[test]
    fn rdkit_compat_tpsa_reference_values() {
        // Reference values captured from RDKit rdMolDescriptors.CalcTPSA with
        // the default N/O-only scope.
        let cases = [
            ("methane", "C", 0.00),
            ("benzene", "c1ccccc1", 0.00),
            ("water", "O", 31.50),
            ("ethanol", "CCO", 20.23),
            ("phenol", "c1ccccc1O", 20.23),
            ("aspirin", "CC(=O)Oc1ccccc1C(=O)O", 63.60),
            ("salicylic acid", "OC(=O)c1ccccc1O", 57.53),
            ("paracetamol", "CC(=O)Nc1ccc(O)cc1", 49.33),
            ("ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1", 37.30),
            (
                "cholesterol",
                "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C",
                20.23,
            ),
        ];

        for (name, smiles, expected) in cases {
            assert_close(name, tpsa_ffi(smiles), expected, 0.01);
        }
    }

    #[test]
    fn rdkit_compat_maccs_ffi_exact_on_bits() {
        // Exact on-bit sets captured from RDKit MACCSkeys.GenMACCSKeys.
        // Use lowercase-aromatic SMILES because Kekule aromaticity perception is
        // outside the current parser's claimed compatibility surface.
        let cases: &[(&str, &[usize])] = &[
            ("CCO", &[82, 109, 114, 139, 153, 155, 157, 160, 164]),
            ("c1ccccc1", &[162, 163, 165]),
            (
                "CC(=O)Oc1ccccc1C(=O)O",
                &[
                    89, 113, 123, 126, 127, 136, 139, 140, 143, 144, 146, 150, 152, 154, 157, 159,
                    160, 162, 163, 164, 165,
                ],
            ),
            (
                "CC(=O)Nc1ccc(O)cc1",
                &[
                    92, 110, 113, 117, 127, 131, 133, 135, 139, 143, 151, 152, 154, 156, 157, 158,
                    159, 160, 161, 162, 163, 164, 165,
                ],
            ),
            ("c1ccc2ccccc2c1", &[101, 105, 125, 145, 162, 163, 165]),
            ("[Na+].[Cl-]", &[35, 44, 49, 103, 134, 166]),
        ];

        for (smiles, expected) in cases {
            let got = on_bits(&maccs_ffi(smiles));
            assert_eq!(
                got,
                expected.to_vec(),
                "MACCS on-bits mismatch for {smiles}"
            );
        }
    }

    #[test]
    fn rdkit_compat_tanimoto_on_rdkit_maccs_bits() {
        // Since the MACCS bitsets above are RDKit-exact, these Tanimoto values
        // match RDKit DataStructs.TanimotoSimilarity for the same fingerprints.
        let aspirin = maccs_ffi("CC(=O)Oc1ccccc1C(=O)O");
        let paracetamol = maccs_ffi("CC(=O)Nc1ccc(O)cc1");
        let benzene = maccs_ffi("c1ccccc1");
        let naphthalene = maccs_ffi("c1ccc2ccccc2c1");
        let ethanol = maccs_ffi("CCO");

        assert_close(
            "aspirin/paracetamol MACCS Tanimoto",
            tanimoto_ffi(&aspirin, &paracetamol),
            13.0 / 31.0,
            1e-12,
        );
        assert_close(
            "aspirin/benzene MACCS Tanimoto",
            tanimoto_ffi(&aspirin, &benzene),
            3.0 / 21.0,
            1e-12,
        );
        assert_close(
            "benzene/naphthalene MACCS Tanimoto",
            tanimoto_ffi(&benzene, &naphthalene),
            3.0 / 7.0,
            1e-12,
        );
        assert_close(
            "ethanol/benzene MACCS Tanimoto",
            tanimoto_ffi(&ethanol, &benzene),
            0.0,
            1e-12,
        );
    }

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
    fn test_tpsa_ffi() {
        let smiles = b"CC(=O)Oc1ccccc1C(=O)O";
        let got = ds_tpsa(smiles.as_ptr(), smiles.len());
        assert!((got - 63.60).abs() < 0.01);
    }

    #[test]
    fn test_tpsa_ffi_invalid_is_nan() {
        let smiles = b"not_a_molecule";
        assert!(ds_tpsa(smiles.as_ptr(), smiles.len()).is_nan());
    }

    #[test]
    fn test_maccs_ffi() {
        let smiles = b"CCO";
        let mut buf = [0u8; 21];
        let n = ds_maccs_keys(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(n, 21);
        // bit 164 (any oxygen) → byte 20, offset 4
        assert_ne!(buf[164 / 8] & (1 << (164 % 8)), 0);
        // bit 0 and 1 never set
        assert_eq!(buf[0] & 0b11, 0);
    }

    #[test]
    fn test_maccs_ffi_invalid_is_neg1() {
        let smiles = b"not_a_molecule";
        let mut buf = [0u8; 21];
        let n = ds_maccs_keys(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(n, -1);
    }

    #[test]
    fn test_maccs_ffi_small_buffer_is_neg1() {
        let smiles = b"CCO";
        let mut buf = [0u8; 20]; // too small (need 21)
        let n = ds_maccs_keys(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
        assert_eq!(n, -1);
    }

    #[test]
    fn test_invalid() {
        assert!(parse("").is_none());
        assert!(parse("not_a_molecule").is_none());
        assert!(parse("C(C").is_none());
    }
}
