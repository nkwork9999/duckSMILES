//! Versioned native graph sidecar, not a pretrained ms-pred tensor format.
//! Atom indices follow the input parse order. No canonicalization, H expansion,
//! neutralization, CIP/EZ assignment or RDKit sanitization is performed.
use crate::parser::{parse_with_ring_info, BondOrder, Molecule, RingInfo};
use std::collections::HashSet;
use std::fmt::Write;

pub const SCHEMA: &str = "ducksmiles-native-graph-v1";
pub const MAX_INPUT_BYTES: usize = 4096;
pub const MAX_ATOMS: usize = 256;
pub const MAX_BONDS: usize = 512;

const ELEMENTS: [&str; 118] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
];

fn atomic_number(symbol: &str) -> Option<usize> {
    ELEMENTS.iter().position(|&s| s == symbol).map(|i| i + 1)
}

// This API-specific conservative grammar guard avoids permissive native-parser
// cases (dangling bonds, empty branches, overflowing isotope/map annotations).
// It does not change the parser used by any existing SQL function.
fn preflight(smiles: &str, allow_dative: bool, max_atoms: usize) -> bool {
    if smiles.is_empty() || smiles.len() > MAX_INPUT_BYTES || !smiles.is_ascii() {
        return false;
    }
    let b = smiles.as_bytes();
    let mut i = 0;
    let mut state = b'.'; // A atom/ring/closed branch, B bond, ( open, . start/dot
    let mut branches = Vec::new();
    let mut atoms = 0;
    while i < b.len() {
        match b[i] {
            b'[' => {
                let start = i + 1;
                let Some(relative_end) = b[start..].iter().position(|&x| x == b']') else {
                    return false;
                };
                let end = start + relative_end;
                let mut p = start;
                while p < end && b[p].is_ascii_digit() {
                    p += 1;
                }
                if p > start && smiles[start..p].parse::<u16>().is_err() {
                    return false;
                }
                if let Some(map_start) = b[start..end].iter().position(|&x| x == b':') {
                    if smiles[start + map_start + 1..end].parse::<u32>().is_err() {
                        return false;
                    }
                }
                atoms += 1;
                state = b'A';
                i = end + 1;
            }
            b'A'..=b'Z' | b'b' | b'c' | b'n' | b'o' | b'p' | b's' => {
                if i + 1 < b.len()
                    && ((b[i] == b'C' && b[i + 1] == b'l') || (b[i] == b'B' && b[i + 1] == b'r'))
                {
                    i += 1;
                }
                atoms += 1;
                state = b'A';
                i += 1;
            }
            b'(' => {
                if state != b'A' {
                    return false;
                }
                branches.push(atoms);
                state = b'(';
                i += 1;
            }
            b')' => {
                if state != b'A' || branches.pop().is_none_or(|n| n == atoms) {
                    return false;
                }
                i += 1;
            }
            b'.' => {
                if state != b'A' || !branches.is_empty() {
                    return false;
                }
                state = b'.';
                i += 1;
            }
            b'<' if allow_dative && b.get(i + 1) == Some(&b'-') => {
                if state != b'A' && state != b'(' {
                    return false;
                }
                state = b'B';
                i += 2;
            }
            b'-' | b'=' | b'#' | b':' | b'/' | b'\\' => {
                if state != b'A' && state != b'(' {
                    return false;
                }
                state = b'B';
                i += if allow_dative && b[i] == b'-' && b.get(i + 1) == Some(&b'>') {
                    2
                } else {
                    1
                };
            }
            b'0'..=b'9' | b'%' => {
                if state != b'A' && state != b'B' {
                    return false;
                }
                if b[i] == b'%' {
                    if i + 2 >= b.len() || !b[i + 1].is_ascii_digit() || !b[i + 2].is_ascii_digit()
                    {
                        return false;
                    }
                    i += 3;
                } else {
                    i += 1;
                }
                state = b'A';
            }
            _ => return false,
        }
        if atoms > max_atoms {
            return false;
        }
    }
    state == b'A' && branches.is_empty()
}

pub(crate) fn graph(smiles: &str) -> Option<Molecule> {
    graph_with_ring_info(smiles).map(|(mol, _)| mol)
}

pub(crate) fn graph_with_ring_info(smiles: &str) -> Option<(Molecule, RingInfo)> {
    checked_graph(smiles, false)
}

#[cfg(test)]
pub(crate) fn graph_as_written_with_ring_info(smiles: &str) -> Option<(Molecule, RingInfo)> {
    checked_graph(smiles, true)
}

pub(crate) fn graph_as_written_bounded(
    smiles: &str,
    max_atoms: usize,
    max_bonds: usize,
) -> Option<Molecule> {
    if !preflight(smiles, true, max_atoms) {
        return None;
    }
    let mol = crate::parser::parse_as_written_graph(smiles)?;
    supported_graph(&mol, max_atoms, max_bonds).then_some(mol)
}

fn checked_graph(smiles: &str, as_written: bool) -> Option<(Molecule, RingInfo)> {
    if !preflight(smiles, as_written, MAX_ATOMS) {
        return None;
    }
    let (mol, rings) = if as_written {
        crate::parser::parse_as_written_with_ring_info(smiles)?
    } else {
        parse_with_ring_info(smiles)?
    };
    supported_graph(&mol, MAX_ATOMS, MAX_BONDS).then_some((mol, rings))
}

fn supported_graph(mol: &Molecule, max_atoms: usize, max_bonds: usize) -> bool {
    if mol.atoms.len() > max_atoms || mol.bonds.len() > max_bonds {
        return false;
    }
    if mol
        .atoms
        .iter()
        .any(|atom| atomic_number(&atom.symbol).is_none())
    {
        return false;
    }
    let mut edges = HashSet::new();
    if mol
        .bonds
        .iter()
        .any(|b| b.a == b.b || !edges.insert((b.a.min(b.b), b.a.max(b.b))))
    {
        return false;
    }
    true
}

fn json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn number_or_null<T: std::fmt::Display>(value: Option<T>) -> String {
    value
        .map(|n| n.to_string())
        .unwrap_or_else(|| "null".to_string())
}

/// Deterministic JSON for this exact input and backend/schema version.
/// None means invalid, unsupported or resource-limited input, not a zero graph.
pub fn calculate(smiles: &str) -> Option<String> {
    let mol = graph(smiles)?;
    let mut out = format!(
        "{{\"schema\":\"{SCHEMA}\",\"backend\":\"ducksmiles-native\",\"backend_version\":\"{}\",\"input_smiles\":{},\"atom_order\":\"input_parse_order\",\"hydrogen_policy\":\"preserve_explicit_vertices_no_expansion\",\"pretrained_tensor_compatible\":false,\"computed\":{{\"native_aromaticity\":true,\"cip\":false,\"ez\":false,\"gasteiger\":false,\"hybridization\":false,\"ring_basis\":false,\"radical_electrons\":false}},\"atom_count\":{},\"bond_count\":{},\"atoms\":[",
        env!("CARGO_PKG_VERSION"), json_string(smiles), mol.atoms.len(), mol.bonds.len()
    );
    for (i, atom) in mol.atoms.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let explicit_neighbors = mol
            .neighbors(i)
            .iter()
            .filter(|(j, _)| mol.atoms[*j].symbol == "H")
            .count();
        let bracket_h = if atom.in_bracket { atom.hydrogen } else { 0 };
        let implicit_h = if atom.in_bracket { 0 } else { atom.hydrogen };
        let chirality = atom
            .chirality
            .as_ref()
            .map(|s| json_string(s))
            .unwrap_or_else(|| "null".into());
        let neighbors = if atom.chirality.is_some() {
            atom.nbr_order
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let _ = write!(out,
            "{{\"index\":{i},\"symbol\":{},\"atomic_number\":{},\"isotope\":{},\"formal_charge\":{},\"aromatic\":{},\"in_bracket\":{},\"bracket_h_count\":{bracket_h},\"implicit_h_count\":{implicit_h},\"explicit_h_neighbor_count\":{explicit_neighbors},\"atom_map\":{},\"chirality_token\":{chirality},\"stereo_neighbor_order\":[{neighbors}]}}",
            json_string(&atom.symbol), atomic_number(&atom.symbol)?, number_or_null(atom.isotope),
            atom.charge, atom.aromatic, atom.in_bracket, number_or_null(atom.atom_map));
    }
    out.push_str("],\"bonds\":[");
    for (i, bond) in mol.bonds.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let kind = match bond.order {
            BondOrder::Single => "single",
            BondOrder::Double => "double",
            BondOrder::Triple => "triple",
            BondOrder::Aromatic => "aromatic",
            BondOrder::Dative => "dative",
        };
        let _ = write!(
            out,
            "{{\"index\":{i},\"begin\":{},\"end\":{},\"bond_type\":\"{kind}\",\"direction\":{}}}",
            bond.a, bond.b, bond.direction
        );
    }
    out.push_str("]}");
    Some(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn input_order_and_hydrogen_semantics() {
        let mol = graph("[H][13C@@H:7](N)C(=O)O").unwrap();
        assert_eq!(mol.atoms[0].symbol, "H");
        assert_eq!(mol.atoms[1].isotope, Some(13));
        assert_eq!(mol.atoms[1].atom_map, Some(7));
        assert_eq!(mol.atoms[1].hydrogen, 1);
        assert_eq!(mol.atoms[1].chirality.as_deref(), Some("@@"));
        assert_eq!(mol.atoms[1].nbr_order, [0, -1, 2, 3]);
        let json = calculate("[H][13C@@H:7](N)C(=O)O").unwrap();
        assert!(json.contains(
            "\"bracket_h_count\":1,\"implicit_h_count\":0,\"explicit_h_neighbor_count\":1"
        ));
        assert!(json.contains("\"stereo_neighbor_order\":[0,-1,2,3]"));
        assert_eq!(graph("[H]O[H]").unwrap().atoms.len(), 3);
        assert_eq!(graph("O").unwrap().atoms.len(), 1);
    }

    #[test]
    fn charges_aromaticity_bond_direction_and_components() {
        let mol = graph("[NH4+].[Cl-]").unwrap();
        assert_eq!(mol.atoms[0].charge, 1);
        assert_eq!(mol.atoms[1].charge, -1);
        assert!(mol.bonds.is_empty());
        assert!(graph("C1=CC=CC=C1")
            .unwrap()
            .atoms
            .iter()
            .all(|a| a.aromatic));
        let mol = graph("F/C=C\\F").unwrap();
        assert_eq!(mol.bonds[0].direction, 1);
        assert_eq!(mol.bonds[1].order, BondOrder::Double);
        assert_eq!(mol.bonds[2].direction, -1);
        assert!(calculate("F/C=C\\F").unwrap().contains("F/C=C\\\\F"));
    }

    #[test]
    fn reject_malformed_unsupported_and_lossy_inputs() {
        for s in [
            "",
            "C(",
            "C1CC",
            "C=",
            "=C",
            "C()",
            ".C",
            "C.",
            "C..N",
            "C=(N)",
            "C(-)N",
            "C(C.)N",
            "C11",
            "C1C1",
            "[Xx]",
            "[*]",
            "*",
            "[99999C]",
            "[C:4294967296]",
            "C C",
            "C\0C",
            "炭素",
        ] {
            assert!(calculate(s).is_none(), "{s:?}");
        }
        assert!(calculate(&"C".repeat(MAX_ATOMS + 1)).is_none());
        assert!(calculate(&"C".repeat(MAX_INPUT_BYTES + 1)).is_none());
        assert!(calculate("C(=O)N").is_some());
        assert!(calculate("C%10CCCCC%10").is_some());
        assert!(calculate("C=1CCCCC1").is_some());
    }

    #[test]
    fn version_determinism_and_no_pretrained_compatibility_claim() {
        let text = calculate("CN").unwrap();
        assert_eq!(Some(text.clone()), calculate("CN"));
        assert_ne!(Some(text.clone()), calculate("NC"));
        assert!(text.contains(SCHEMA));
        assert!(text.contains("\"pretrained_tensor_compatible\":false"));
        assert!(text.contains("\"gasteiger\":false"));
        assert!(text.contains("\"radical_electrons\":false"));
        assert_eq!(atomic_number("N"), Some(7));
        assert_eq!(atomic_number("Og"), Some(118));
        assert_eq!(atomic_number("D"), None);
    }

    #[test]
    fn ffi_validates_input_and_obeys_exact_buffer_contract() {
        use crate::ds_mol_graph_json;
        use std::ptr::null_mut;
        // SAFETY: normal calls use owned readable/writable disjoint buffers;
        // null and oversized inputs exercise guards before any dereference.
        unsafe {
            assert_eq!(ds_mol_graph_json(std::ptr::null(), 1, null_mut(), 0), -1);
            assert_eq!(ds_mol_graph_json([0xff].as_ptr(), 1, null_mut(), 0), -1);
            let smiles = "N[C@@H](C)C(=O)O";
            // Oversized lengths are rejected before dereferencing the input.
            assert_eq!(
                ds_mol_graph_json(smiles.as_ptr(), usize::MAX, null_mut(), 0),
                -1
            );
            let needed = ds_mol_graph_json(smiles.as_ptr(), smiles.len(), null_mut(), 0);
            assert!(needed > 0);
            let mut short = [0x55; 4];
            assert_eq!(
                ds_mol_graph_json(
                    smiles.as_ptr(),
                    smiles.len(),
                    short.as_mut_ptr(),
                    short.len()
                ),
                needed
            );
            assert_eq!(short, [0x55; 4]);
            let mut result = vec![0x55; needed as usize + 1];
            assert_eq!(
                ds_mol_graph_json(
                    smiles.as_ptr(),
                    smiles.len(),
                    result.as_mut_ptr(),
                    needed as usize
                ),
                needed
            );
            assert_eq!(result[needed as usize], 0x55);
            assert_eq!(
                std::str::from_utf8(&result[..needed as usize]).unwrap(),
                calculate(smiles).unwrap()
            );
        }
    }
}
