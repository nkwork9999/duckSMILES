//! MACCS keys fingerprint, ported from RDKit `MACCS.cpp` (BSD 3-Clause).
//!
//! RDKit defines 166 public MACCS keys stored in a 167-bit vector: bit 0 is
//! unused and bit 1 (ISOTOPE) is left permanently 0, so the meaningful keys are
//! bits 2..=166. Each bit is set by one of three mechanisms:
//!
//!   1. **Presence** — a SMARTS substructure match exists (most bits).
//!   2. **Count threshold** — a SMARTS pattern matches more than N times
//!      (RDKit's `SubstructMatch(..., uniquify=true)` count). Several bits share
//!      one count (e.g. the oxygen-count cascade 140/146/159/164).
//!   3. **Special** — non-SMARTS C++ logic: a per-atom atomic-number scan
//!      (metals, halogens, S, P, …), the "≥2 aromatic rings" rule (bit 125),
//!      and the multi-fragment rule (bit 166).
//!
//! The output is a 21-byte BLOB (`ceil(167/8) = 21`). Bit `n` lives in
//! `byte = n / 8`, `offset = n % 8` (little-endian within each byte), matching
//! the convention used by ducksmiles' Morgan BLOB and the article's `set_bit`.
//!
//! Reference: RDKit `Code/GraphMol/Fingerprints/MACCS.cpp`.

use crate::parser::Molecule;
use crate::smarts::{count_unique, matches_mol, parse_smarts, Pattern};
use std::sync::OnceLock;

pub const MACCS_N_BITS: usize = 167;
pub const MACCS_N_BYTES: usize = 21; // ceil(167 / 8)

/// One presence key: set `bit` if `smarts` matches anywhere.
struct PresenceKey {
    bit: usize,
    smarts: &'static str,
}

/// One count key: compute the unique-match count of `smarts` once, then set
/// each `(threshold, bit)` whose `count > threshold`.
struct CountKey {
    smarts: &'static str,
    thresholds: &'static [(usize, usize)], // (count_strictly_greater_than, bit)
}

// -----------------------------------------------------------------------------
// Presence keys (set bit if SMARTS matches at least once)
// -----------------------------------------------------------------------------
const PRESENCE_KEYS: &[PresenceKey] = &[
    PresenceKey { bit: 8,  smarts: "[!#6!#1]1~*~*~*~1" },
    PresenceKey { bit: 11, smarts: "*1~*~*~*~1" },
    PresenceKey { bit: 13, smarts: "[#8]~[#7](~[#6])~[#6]" },
    PresenceKey { bit: 14, smarts: "[#16]-[#16]" },
    PresenceKey { bit: 15, smarts: "[#8]~[#6](~[#8])~[#8]" },
    PresenceKey { bit: 16, smarts: "[!#6!#1]1~*~*~1" },
    PresenceKey { bit: 17, smarts: "[#6]#[#6]" },
    PresenceKey { bit: 19, smarts: "*1~*~*~*~*~*~*~1" },
    PresenceKey { bit: 20, smarts: "[#14]" },
    PresenceKey { bit: 21, smarts: "[#6]=[#6](~[!#6!#1])~[!#6!#1]" },
    PresenceKey { bit: 22, smarts: "*1~*~*~1" },
    PresenceKey { bit: 23, smarts: "[#7]~[#6](~[#8])~[#8]" },
    PresenceKey { bit: 24, smarts: "[#7]-[#8]" },
    PresenceKey { bit: 25, smarts: "[#7]~[#6](~[#7])~[#7]" },
    PresenceKey { bit: 26, smarts: "[#6]=@[#6](@*)@*" },
    PresenceKey { bit: 28, smarts: "[!#6!#1]~[CH2]~[!#6!#1]" },
    PresenceKey { bit: 30, smarts: "[#6]~[!#6!#1](~[#6])(~[#6])~*" },
    PresenceKey { bit: 31, smarts: "[!#6!#1]~[F,Cl,Br,I]" },
    PresenceKey { bit: 32, smarts: "[#6]~[#16]~[#7]" },
    PresenceKey { bit: 33, smarts: "[#7]~[#16]" },
    PresenceKey { bit: 34, smarts: "[CH2]=*" },
    PresenceKey { bit: 36, smarts: "[#16R]" },
    PresenceKey { bit: 37, smarts: "[#7]~[#6](~[#8])~[#7]" },
    PresenceKey { bit: 38, smarts: "[#7]~[#6](~[#6])~[#7]" },
    PresenceKey { bit: 39, smarts: "[#8]~[#16](~[#8])~[#8]" },
    PresenceKey { bit: 40, smarts: "[#16]-[#8]" },
    PresenceKey { bit: 41, smarts: "[#6]#[#7]" },
    PresenceKey { bit: 43, smarts: "[!#6!#1!H0]~*~[!#6!#1!H0]" },
    PresenceKey { bit: 44, smarts: "[!#1;!#6;!#7;!#8;!#9;!#14;!#15;!#16;!#17;!#35;!#53]" },
    PresenceKey { bit: 45, smarts: "[#6]=[#6]~[#7]" },
    PresenceKey { bit: 47, smarts: "[#16]~*~[#7]" },
    PresenceKey { bit: 48, smarts: "[#8]~[!#6!#1](~[#8])~[#8]" },
    PresenceKey { bit: 49, smarts: "[!+0]" },
    PresenceKey { bit: 50, smarts: "[#6]=[#6](~[#6])~[#6]" },
    PresenceKey { bit: 51, smarts: "[#6]~[#16]~[#8]" },
    PresenceKey { bit: 52, smarts: "[#7]~[#7]" },
    PresenceKey { bit: 53, smarts: "[!#6!#1!H0]~*~*~*~[!#6!#1!H0]" },
    PresenceKey { bit: 54, smarts: "[!#6!#1!H0]~*~*~[!#6!#1!H0]" },
    PresenceKey { bit: 55, smarts: "[#8]~[#16]~[#8]" },
    PresenceKey { bit: 56, smarts: "[#8]~[#7](~[#8])~[#6]" },
    PresenceKey { bit: 57, smarts: "[#8R]" },
    PresenceKey { bit: 58, smarts: "[!#6!#1]~[#16]~[!#6!#1]" },
    PresenceKey { bit: 59, smarts: "[#16]!:*:*" },
    PresenceKey { bit: 60, smarts: "[#16]=[#8]" },
    PresenceKey { bit: 61, smarts: "*~[#16](~*)~*" },
    PresenceKey { bit: 62, smarts: "*@*!@*@*" },
    PresenceKey { bit: 63, smarts: "[#7]=[#8]" },
    PresenceKey { bit: 64, smarts: "*@*!@[#16]" },
    PresenceKey { bit: 65, smarts: "c:n" },
    PresenceKey { bit: 66, smarts: "[#6]~[#6](~[#6])(~[#6])~*" },
    PresenceKey { bit: 67, smarts: "[!#6!#1]~[#16]" },
    PresenceKey { bit: 68, smarts: "[!#6!#1!H0]~[!#6!#1!H0]" },
    PresenceKey { bit: 69, smarts: "[!#6!#1]~[!#6!#1!H0]" },
    PresenceKey { bit: 70, smarts: "[!#6!#1]~[#7]~[!#6!#1]" },
    PresenceKey { bit: 71, smarts: "[#7]~[#8]" },
    PresenceKey { bit: 72, smarts: "[#8]~*~*~[#8]" },
    PresenceKey { bit: 73, smarts: "[#16]=*" },
    PresenceKey { bit: 74, smarts: "[CH3]~*~[CH3]" },
    PresenceKey { bit: 75, smarts: "*!@[#7]@*" },
    PresenceKey { bit: 76, smarts: "[#6]=[#6](~*)~*" },
    PresenceKey { bit: 77, smarts: "[#7]~*~[#7]" },
    PresenceKey { bit: 78, smarts: "[#6]=[#7]" },
    PresenceKey { bit: 79, smarts: "[#7]~*~*~[#7]" },
    PresenceKey { bit: 80, smarts: "[#7]~*~*~*~[#7]" },
    PresenceKey { bit: 81, smarts: "[#16]~*(~*)~*" },
    PresenceKey { bit: 82, smarts: "*~[CH2]~[!#6!#1!H0]" },
    PresenceKey { bit: 83, smarts: "[!#6!#1]1~*~*~*~*~1" },
    PresenceKey { bit: 84, smarts: "[NH2]" },
    PresenceKey { bit: 85, smarts: "[#6]~[#7](~[#6])~[#6]" },
    PresenceKey { bit: 86, smarts: "[C;H2,H3][!#6!#1][C;H2,H3]" },
    PresenceKey { bit: 87, smarts: "[F,Cl,Br,I]!@*@*" },
    PresenceKey { bit: 89, smarts: "[#8]~*~*~*~[#8]" },
    PresenceKey { bit: 90, smarts: "[$([!#6!#1!H0]~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[CH2R]1)]" },
    PresenceKey { bit: 91, smarts: "[$([!#6!#1!H0]~*~*~*~[CH2]~*),$([!#6!#1!H0R]1@[R]@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~[R]1@[R]@[R]@[CH2R]1),$([!#6!#1!H0]~*~[R]1@[R]@[CH2R]1)]" },
    PresenceKey { bit: 92, smarts: "[#8]~[#6](~[#7])~[#6]" },
    PresenceKey { bit: 93, smarts: "[!#6!#1]~[CH3]" },
    PresenceKey { bit: 94, smarts: "[!#6!#1]~[#7]" },
    PresenceKey { bit: 95, smarts: "[#7]~*~*~[#8]" },
    PresenceKey { bit: 96, smarts: "*1~*~*~*~*~1" },
    PresenceKey { bit: 97, smarts: "[#7]~*~*~*~[#8]" },
    PresenceKey { bit: 98, smarts: "[!#6!#1]1~*~*~*~*~*~1" },
    PresenceKey { bit: 99, smarts: "[#6]=[#6]" },
    PresenceKey { bit: 100, smarts: "*~[CH2]~[#7]" },
    PresenceKey { bit: 101, smarts: "[$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1),$([R]1@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@[R]@1)]" },
    PresenceKey { bit: 102, smarts: "[!#6!#1]~[#8]" },
    PresenceKey { bit: 104, smarts: "[!#6!#1!H0]~*~[CH2]~*" },
    PresenceKey { bit: 105, smarts: "*@*(@*)@*" },
    PresenceKey { bit: 106, smarts: "[!#6!#1]~*(~[!#6!#1])~[!#6!#1]" },
    PresenceKey { bit: 107, smarts: "[F,Cl,Br,I]~*(~*)~*" },
    PresenceKey { bit: 108, smarts: "[CH3]~*~*~*~[CH2]~*" },
    PresenceKey { bit: 109, smarts: "*~[CH2]~[#8]" },
    PresenceKey { bit: 110, smarts: "[#7]~[#6]~[#8]" },
    PresenceKey { bit: 111, smarts: "[#7]~*~[CH2]~*" },
    PresenceKey { bit: 112, smarts: "*~*(~*)(~*)~*" },
    PresenceKey { bit: 113, smarts: "[#8]!:*:*" },
    PresenceKey { bit: 114, smarts: "[CH3]~[CH2]~*" },
    PresenceKey { bit: 115, smarts: "[CH3]~*~[CH2]~*" },
    PresenceKey { bit: 116, smarts: "[$([CH3]~*~*~[CH2]~*),$([CH3]~*1~*~[CH2]1)]" },
    PresenceKey { bit: 117, smarts: "[#7]~*~[#8]" },
    PresenceKey { bit: 119, smarts: "[#7]=*" },
    PresenceKey { bit: 121, smarts: "[#7R]" },
    PresenceKey { bit: 122, smarts: "*~[#7](~*)~*" },
    PresenceKey { bit: 123, smarts: "[#8]~[#6]~[#8]" },
    PresenceKey { bit: 126, smarts: "*!@[#8]!@*" },
    PresenceKey { bit: 128, smarts: "[$(*~[CH2]~*~*~*~[CH2]~*),$([R]1@[CH2R]@[R]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[R]@[CH2R]1),$(*~[CH2]~*~[R]1@[R]@[CH2R]1)]" },
    PresenceKey { bit: 129, smarts: "[$(*~[CH2]~*~*~[CH2]~*),$([R]1@[CH2]@[R]@[R]@[CH2R]1),$(*~[CH2]~[R]1@[R]@[CH2R]1)]" },
    PresenceKey { bit: 132, smarts: "[#8]~*~[CH2]~*" },
    PresenceKey { bit: 133, smarts: "*@*!@[#7]" },
    PresenceKey { bit: 135, smarts: "[#7]!:*:*" },
    PresenceKey { bit: 137, smarts: "[!C!cR]" },
    PresenceKey { bit: 139, smarts: "[O!H0]" },
    PresenceKey { bit: 144, smarts: "*!:*:*!:*" },
    PresenceKey { bit: 147, smarts: "[$(*~[CH2]~[CH2]~*),$([R]1@[CH2R]@[CH2R]1)]" },
    PresenceKey { bit: 148, smarts: "*~[!#6!#1](~*)~*" },
    PresenceKey { bit: 150, smarts: "*!@*@*!@*" },
    PresenceKey { bit: 151, smarts: "[#7!H0]" },
    PresenceKey { bit: 152, smarts: "[#8]~[#6](~[#6])~[#6]" },
    PresenceKey { bit: 154, smarts: "[#6]=[#8]" },
    PresenceKey { bit: 155, smarts: "*!@[CH2]!@*" },
    PresenceKey { bit: 156, smarts: "[#7]~*(~*)~*" },
    PresenceKey { bit: 157, smarts: "[#6]-[#8]" },
    PresenceKey { bit: 158, smarts: "[#6]-[#7]" },
    PresenceKey { bit: 162, smarts: "a" },
    PresenceKey { bit: 165, smarts: "[R]" },
];

// -----------------------------------------------------------------------------
// Count keys (one SMARTS, possibly multiple count thresholds → multiple bits)
// -----------------------------------------------------------------------------
const COUNT_KEYS: &[CountKey] = &[
    CountKey { smarts: "[$(*~[CH2]~[CH2]~*),$(*1~[CH2]~[CH2]1)]", thresholds: &[(1, 118)] },
    CountKey { smarts: "[!#6R]", thresholds: &[(1, 120)] },
    CountKey { smarts: "[!#6!#1]~[!#6!#1]", thresholds: &[(0, 124), (1, 130)] },
    CountKey { smarts: "*@*!@[#8]", thresholds: &[(1, 127), (0, 143)] },
    CountKey { smarts: "[!#6!#1!H0]", thresholds: &[(1, 131)] },
    CountKey { smarts: "[#8]=*", thresholds: &[(1, 136)] },
    CountKey { smarts: "[!#6!#1]~[CH2]~*", thresholds: &[(1, 138), (0, 153)] },
    CountKey { smarts: "[#8]", thresholds: &[(3, 140), (2, 146), (1, 159), (0, 164)] },
    CountKey { smarts: "[CH3]", thresholds: &[(2, 141)] },
    CountKey { smarts: "[#7]", thresholds: &[(1, 142), (0, 161)] },
    CountKey { smarts: "*1~*~*~*~*~*~1", thresholds: &[(1, 145), (0, 163)] },
    CountKey { smarts: "[C;H3,H4]", thresholds: &[(1, 149), (0, 160)] },
];

struct Compiled {
    presence: Vec<(usize, Pattern)>,
    counts: Vec<(Pattern, &'static [(usize, usize)])>,
}

fn compiled() -> &'static Compiled {
    static CACHE: OnceLock<Compiled> = OnceLock::new();
    CACHE.get_or_init(|| {
        let presence = PRESENCE_KEYS
            .iter()
            .map(|k| {
                let pat = parse_smarts(k.smarts)
                    .unwrap_or_else(|| panic!("MACCS presence SMARTS failed to parse: {}", k.smarts));
                (k.bit, pat)
            })
            .collect();
        let counts = COUNT_KEYS
            .iter()
            .map(|k| {
                let pat = parse_smarts(k.smarts)
                    .unwrap_or_else(|| panic!("MACCS count SMARTS failed to parse: {}", k.smarts));
                (pat, k.thresholds)
            })
            .collect();
        Compiled { presence, counts }
    })
}

#[inline]
fn set_bit(buf: &mut [u8; MACCS_N_BYTES], bit: usize) {
    if bit >= MACCS_N_BITS {
        return;
    }
    buf[bit / 8] |= 1 << (bit % 8);
}

/// Compute the 166 MACCS keys for a molecule, returning a 21-byte bit buffer.
///
/// Bit 0 and bit 1 are never set (RDKit leaves bit 1, the isotope key, off).
pub fn maccs_bits(mol: &Molecule) -> [u8; MACCS_N_BYTES] {
    let mut fp = [0u8; MACCS_N_BYTES];
    if mol.atoms.is_empty() {
        return fp;
    }

    let comp = compiled();

    // --- Presence keys ---
    for (bit, pat) in &comp.presence {
        if matches_mol(pat, mol) {
            set_bit(&mut fp, *bit);
        }
    }

    // --- Count keys (shared counts computed once each) ---
    for (pat, thresholds) in &comp.counts {
        let count = count_unique(pat, mol);
        for &(thr, bit) in thresholds.iter() {
            if count > thr {
                set_bit(&mut fp, bit);
            }
        }
    }

    // --- Special: per-atom atomic-number scan (metals, halogens, S, P, …) ---
    for atom in &mol.atoms {
        match atomic_num(&atom.symbol) {
            3 | 11 | 19 | 37 | 55 | 87 => set_bit(&mut fp, 35),
            4 | 12 | 20 | 38 | 56 | 88 => set_bit(&mut fp, 10),
            5 | 13 | 31 | 49 | 81 => set_bit(&mut fp, 18),
            9 => { set_bit(&mut fp, 42); set_bit(&mut fp, 134); }
            15 => set_bit(&mut fp, 29),
            16 => set_bit(&mut fp, 88),
            17 => { set_bit(&mut fp, 103); set_bit(&mut fp, 134); }
            21 | 22 | 39 | 40 | 72 => set_bit(&mut fp, 5),
            23 | 24 | 25 | 41 | 42 | 43 | 73 | 74 | 75 => set_bit(&mut fp, 7),
            26 | 27 | 28 | 44 | 45 | 46 | 76 | 77 | 78 => set_bit(&mut fp, 9),
            29 | 30 | 47 | 48 | 79 | 80 => set_bit(&mut fp, 12),
            32 | 33 | 34 | 50 | 51 | 52 | 82 | 83 | 84 => set_bit(&mut fp, 3),
            35 => { set_bit(&mut fp, 46); set_bit(&mut fp, 134); }
            53 => { set_bit(&mut fp, 27); set_bit(&mut fp, 134); }
            57..=71 => set_bit(&mut fp, 6),
            89..=103 => set_bit(&mut fp, 4),
            104 => set_bit(&mut fp, 2),
            _ => {}
        }
    }

    // --- Special bit 125: two or more fully-aromatic rings ---
    {
        let ri = mol.ring_info();
        let mut n_arom = 0usize;
        for ring_bonds in &ri.rings {
            let is_arom = ring_bonds
                .iter()
                .all(|&bi| mol.bonds[bi].order == crate::parser::BondOrder::Aromatic);
            if is_arom {
                if n_arom >= 1 {
                    set_bit(&mut fp, 125);
                    break;
                }
                n_arom += 1;
            }
        }
    }

    // --- Special bit 166: more than one fragment ---
    if mol.fragment_count() > 1 {
        set_bit(&mut fp, 166);
    }

    fp
}

/// Atomic number for an element symbol (0 if unknown). Mirrors the table in
/// `smarts.rs`; duplicated here to keep MACCS self-contained.
fn atomic_num(sym: &str) -> u8 {
    match sym {
        "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7,
        "O" => 8, "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14,
        "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20,
        "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25, "Fe" => 26,
        "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30, "Ga" => 31, "Ge" => 32,
        "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37, "Sr" => 38,
        "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44,
        "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50,
        "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56,
        "La" => 57, "Ce" => 58, "Pr" => 59, "Nd" => 60, "Pm" => 61, "Sm" => 62,
        "Eu" => 63, "Gd" => 64, "Tb" => 65, "Dy" => 66, "Ho" => 67, "Er" => 68,
        "Tm" => 69, "Yb" => 70, "Lu" => 71,
        "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77,
        "Pt" => 78, "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83,
        "Po" => 84, "At" => 85, "Rn" => 86, "Fr" => 87, "Ra" => 88,
        "Ac" => 89, "Th" => 90, "Pa" => 91, "U" => 92, "Np" => 93, "Pu" => 94,
        "Am" => 95, "Cm" => 96, "Bk" => 97, "Cf" => 98, "Es" => 99, "Fm" => 100,
        "Md" => 101, "No" => 102, "Lr" => 103, "Rf" => 104,
        _ => 0,
    }
}

/// Convenience: list the on-bit indices (1..=166) for testing/verification.
pub fn on_bits(fp: &[u8; MACCS_N_BYTES]) -> Vec<usize> {
    let mut bits = Vec::new();
    for bit in 0..MACCS_N_BITS {
        if fp[bit / 8] & (1 << (bit % 8)) != 0 {
            bits.push(bit);
        }
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn maccs_on_bits(smiles: &str) -> Vec<usize> {
        let mol = parse(smiles).unwrap();
        on_bits(&maccs_bits(&mol))
    }

    #[test]
    fn all_smarts_parse() {
        // compiled() panics on any unparseable SMARTS — calling it is the test.
        let c = compiled();
        assert_eq!(c.presence.len(), PRESENCE_KEYS.len());
        assert_eq!(c.counts.len(), COUNT_KEYS.len());
    }

    #[test]
    fn empty_molecule_all_zero() {
        let mol = parse("").unwrap_or(crate::parser::Molecule {
            atoms: vec![],
            bonds: vec![],
            bond_count: 0,
        });
        let fp = maccs_bits(&mol);
        assert!(fp.iter().all(|&b| b == 0));
    }

    #[test]
    fn byte_length_is_21() {
        let fp = maccs_bits(&parse("CCO").unwrap());
        assert_eq!(fp.len(), 21);
    }

    #[test]
    fn ethanol_reasonable_bits() {
        // CCO: should at least include oxygen (164), methyl (160), C-O (157).
        let bits = maccs_on_bits("CCO");
        assert!(bits.contains(&164), "ethanol has oxygen (bit 164)");
        assert!(bits.contains(&157), "ethanol has C-O (bit 157)");
        assert!(bits.contains(&139), "ethanol has O with H (bit 139)");
        assert!(!bits.contains(&0));
        assert!(!bits.contains(&1));
    }

    #[test]
    fn benzene_is_aromatic_ring() {
        let bits = maccs_on_bits("c1ccccc1");
        assert!(bits.contains(&162), "benzene aromatic (bit 162)");
        assert!(bits.contains(&165), "benzene ring atoms (bit 165)");
        // single aromatic ring → bit 125 (>=2 aromatic rings) NOT set
        assert!(!bits.contains(&125));
    }

    #[test]
    fn naphthalene_two_aromatic_rings() {
        let bits = maccs_on_bits("c1ccc2ccccc2c1");
        assert!(bits.contains(&125), "naphthalene has >=2 aromatic rings (bit 125)");
    }

    #[test]
    fn salt_sets_fragment_bit() {
        let bits = maccs_on_bits("[Na+].[Cl-]");
        assert!(bits.contains(&166), "salt has >1 fragment (bit 166)");
    }

    #[test]
    fn single_fragment_no_166() {
        let bits = maccs_on_bits("CCO");
        assert!(!bits.contains(&166));
    }

    /// Exact-match regression against Python RDKit `MACCSkeys.GenMACCSKeys`.
    /// These on-bit sets were captured from RDKit and must match bit-for-bit.
    /// All use lowercase-aromatic SMILES (the parser's supported aromatic form);
    /// Kekulé-written aromatic rings depend on aromaticity perception, which is
    /// tracked separately and intentionally excluded here.
    #[test]
    fn rdkit_exact_match_regression() {
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
            (
                "Cn1cnc2c1c(=O)n(C)c(=O)n2C", // caffeine, aromatic form
                &[
                    37, 38, 65, 75, 77, 79, 80, 83, 85, 89, 92, 93, 95, 96, 97, 98, 101, 105, 106,
                    110, 113, 117, 120, 121, 122, 125, 127, 136, 137, 141, 142, 143, 144, 148, 149,
                    150, 154, 156, 158, 159, 160, 161, 162, 163, 164, 165,
                ],
            ),
        ];
        for (smi, expected) in cases {
            let got = maccs_on_bits(smi);
            assert_eq!(
                got,
                expected.to_vec(),
                "MACCS on-bits mismatch for {}\n got: {:?}\n want: {:?}",
                smi,
                got,
                expected
            );
        }
    }
}
