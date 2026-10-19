//! Coordinate-free CX wedge projection for the current 291 descriptors only.
//!
//! `wU`/`wD` without coordinates encode wedge/axial configuration; the current
//! descriptor definitions do not measure the handedness of axial stereoisomers.
//! The caller keeps the original input as the identity/export/cache key. This
//! projection does not provide an axial stereochemical graph or depiction.
//! Coordinates, unknown/wiggly stereo, ring stereo, charges/radicals, labels and
//! every other CX block remain unsupported instead of being silently discarded.
//! A future configuration-sensitive descriptor must extend this representation.
use crate::parser::{BondOrder, Molecule};

pub(crate) fn graph(input: &str) -> Option<Molecule> {
    let Some((head, extension)) = input.split_once('|') else {
        return crate::native_graph::graph_as_written_bounded(
            input,
            crate::workbench::MAX_ATOMS,
            crate::workbench::MAX_BONDS,
        );
    };
    if input.len() > crate::native_graph::MAX_INPUT_BYTES || !input.is_ascii() {
        return None;
    }
    // A whitespace separator is required. Names, extra blocks and unclosed
    // annotations must not become successfully parsed prefixes.
    if !head.ends_with([' ', '\t']) {
        return None;
    }
    let base = head.trim_end_matches([' ', '\t']);
    let extension = extension.strip_suffix('|')?;
    if extension.is_empty() || extension.contains('|') {
        return None;
    }
    let mol = crate::native_graph::graph_as_written_bounded(
        base,
        crate::workbench::MAX_ATOMS,
        crate::workbench::MAX_BONDS,
    )?;
    let mut has_direction = false;
    let mut seen_bonds = vec![false; mol.bonds.len()];
    for token in extension.split(',') {
        let pair = if let Some(pair) = token
            .strip_prefix("wU:")
            .or_else(|| token.strip_prefix("wD:"))
        {
            has_direction = true;
            pair
        } else if has_direction {
            token
        } else {
            return None;
        };
        let (atom, bond) = pair.split_once('.')?;
        if atom.is_empty()
            || bond.is_empty()
            || !atom.bytes().all(|b| b.is_ascii_digit())
            || !bond.bytes().all(|b| b.is_ascii_digit())
        {
            return None;
        }
        let atom: usize = atom.parse().ok()?;
        let bond: usize = bond.parse().ok()?;
        // Native bond order is the CX parse order, with each ring closure at
        // its textual closing digit (not RDKit's ring-bonds-last storage order).
        let edge = mol.bonds.get(bond)?;
        if atom >= mol.atoms.len()
            || (edge.a != atom && edge.b != atom)
            || seen_bonds[bond]
            || edge.order == BondOrder::Dative
        {
            return None;
        }
        // Reversing a dative edge changes chemistry; it is intentionally excluded.
        // Ordinary wedge orientation has no value in this descriptor projection.
        seen_bonds[bond] = true;
    }
    Some(mol)
}
