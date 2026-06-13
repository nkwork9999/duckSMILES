pub(crate) mod params;
pub(crate) mod bounds;
pub(crate) mod smoothing;
pub(crate) mod embed;
pub(crate) mod minimize;

use crate::parser::parse;
use bounds::build_bounds;
use embed::embed;
use minimize::minimize;
use smoothing::smooth;

// ── generate_conformer ────────────────────────────────────────────────────────
// Full pipeline: SMILES → bounds → smooth → embed → minimize → coords.
// Returns None if the SMILES is invalid or the embedding fails.

pub fn generate_conformer(smiles: &str, seed: u64) -> Option<Vec<f64>> {
    let mol_bare = parse(smiles)?;
    let mol = mol_bare.with_explicit_hydrogens();
    let mut bm = build_bounds(&mol);
    smooth(&mut bm);
    let mut coords = embed(&bm, seed)?;
    minimize(&mol, &mut coords);
    Some(coords)
}

// ── atom_count_for_smiles ─────────────────────────────────────────────────────
// Returns the number of atoms (including explicit H) for the molecule.
// Used by the C FFI to allocate the right buffer size.
pub fn atom_count(smiles: &str) -> Option<usize> {
    let mol = parse(smiles)?;
    Some(mol.with_explicit_hydrogens().atoms.len())
}
