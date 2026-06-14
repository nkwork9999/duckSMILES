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
// Full pipeline: SMILES → bounds → smooth → embed (×K) → minimize → best coords.
//
// ETKDG-lite ensemble: distance-geometry embedding is stochastic, so a single
// shot often lands in a poor local minimum (puckered ring, bad clash). We embed
// and minimise `N_TRIES` independent conformers from derived seeds and keep the
// lowest force-field energy — the conformer-ensemble idea RDKit's EmbedMultiple
// uses. Planarity of sp2/aromatic systems is enforced inside `minimize` via the
// out-of-plane term, so the lowest-energy member is also the flattest.
// Returns None if the SMILES is invalid or every embedding fails.

const N_TRIES: u64 = 8;

pub fn generate_conformer(smiles: &str, seed: u64) -> Option<Vec<f64>> {
    let mol_bare = parse(smiles)?;
    let mol = mol_bare.with_explicit_hydrogens();
    let mut bm = build_bounds(&mol);
    smooth(&mut bm);

    let mut best: Option<(f64, Vec<f64>)> = None;
    for k in 0..N_TRIES {
        // Derive distinct seeds deterministically from the caller's seed.
        let s = seed.wrapping_mul(0x9E3779B97F4A7C15).wrapping_add(k.wrapping_add(1));
        let Some(mut coords) = embed(&bm, s) else { continue };
        let energy = minimize(&mol, &mut coords);
        match &best {
            Some((be, _)) if *be <= energy => {}
            _ => best = Some((energy, coords)),
        }
    }

    best.map(|(_, coords)| coords)
}

// ── atom_count_for_smiles ─────────────────────────────────────────────────────
// Returns the number of atoms (including explicit H) for the molecule.
// Used by the C FFI to allocate the right buffer size.
pub fn atom_count(smiles: &str) -> Option<usize> {
    let mol = parse(smiles)?;
    Some(mol.with_explicit_hydrogens().atoms.len())
}
