//! Morgan / ECFP fingerprint.
//!
//! Faithful port of the algorithm in `RDKit/Code/GraphMol/Fingerprints/MorganGenerator.cpp`:
//!
//!   1. Compute round-0 atom invariants (Daylight ECFP style: heavy-degree,
//!      atomic_num, total H count, formal charge, isotope, in_ring).
//!   2. For each round r = 1..=radius:
//!        - For each atom, gather (bond_invariant, neighbor_invariant) pairs,
//!          sort them, then hash_combine(layer, current_inv, [sorted pairs])
//!          to produce that atom's new invariant.
//!        - Drop atoms whose growing neighborhood (set of covered bonds)
//!          duplicates an already-seen one; otherwise emit its invariant.
//!   3. Fold all emitted invariants into a fixed-width bit vector
//!      (`invariant % n_bits` -> set bit).
//!
//! Hash function: we use a boost-style `hash_combine` over u32 identity hashes
//! so the algorithm is deterministic across platforms. This is **not**
//! bit-exact RDKit compatible (RDKit's bits depend on boost::hash internals);
//! it is an ECFP-equivalent fingerprint with the same algorithmic structure.

use crate::parser::{BondOrder, Molecule};

const HASH_SEED_LAYER: u32 = 0x9e37_79b9;

/// Boost-style hash_combine for u32.
#[inline]
fn hash_combine(seed: u32, v: u32) -> u32 {
    seed ^ (v
        .wrapping_add(HASH_SEED_LAYER)
        .wrapping_add(seed << 6)
        .wrapping_add(seed >> 2))
}

/// Periodic table lookup. Covers the elements our SMILES parser emits.
fn atomic_num(sym: &str) -> u32 {
    match sym {
        "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5,
        "C" => 6, "N" => 7, "O" => 8, "F" => 9, "Ne" => 10,
        "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14, "P" => 15,
        "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20,
        "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25,
        "Fe" => 26, "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30,
        "Ga" => 31, "Ge" => 32, "As" => 33, "Se" => 34, "Br" => 35,
        "Kr" => 36, "Rb" => 37, "Sr" => 38, "Ag" => 47, "I" => 53,
        "Cs" => 55, "Ba" => 56, "Pt" => 78, "Au" => 79, "Hg" => 80,
        "Pb" => 82,
        _ => 0,
    }
}

/// For each atom, true if it belongs to at least one cycle.
/// Plain DFS — every non-tree edge marks every atom on the back-edge path.
fn compute_in_ring(mol: &Molecule) -> Vec<bool> {
    let n = mol.atoms.len();
    let mut in_ring = vec![false; n];
    let mut parent: Vec<isize> = vec![-1; n];
    let mut depth = vec![0i32; n];
    let mut visited = vec![false; n];

    // Build adjacency list once
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for b in &mol.bonds {
        adj[b.a].push(b.b);
        adj[b.b].push(b.a);
    }

    for start in 0..n {
        if visited[start] { continue; }
        let mut stack: Vec<(usize, usize)> = vec![(start, 0)]; // (node, next_neighbor_idx)
        visited[start] = true;
        while let Some(&(u, ni)) = stack.last() {
            if ni < adj[u].len() {
                let v = adj[u][ni];
                stack.last_mut().unwrap().1 += 1;
                if !visited[v] {
                    visited[v] = true;
                    parent[v] = u as isize;
                    depth[v] = depth[u] + 1;
                    stack.push((v, 0));
                } else if parent[u] != v as isize {
                    // back-edge u→v: every node on the path from u up to v is in a ring
                    let mut cur = u;
                    while cur != v {
                        in_ring[cur] = true;
                        if parent[cur] < 0 { break; }
                        cur = parent[cur] as usize;
                    }
                    in_ring[v] = true;
                }
            } else {
                stack.pop();
            }
        }
    }
    in_ring
}

/// Heavy-atom degree of an atom (excludes hydrogens that the parser already split out).
fn heavy_degree(mol: &Molecule, atom_idx: usize) -> u32 {
    mol.neighbors(atom_idx)
        .iter()
        .filter(|(n, _)| mol.atoms[*n].symbol != "H")
        .count() as u32
}

/// Initial (round-0) ECFP-style atom invariant.
fn connectivity_invariant(mol: &Molecule, in_ring: &[bool], atom_idx: usize) -> u32 {
    let a = &mol.atoms[atom_idx];
    let an = atomic_num(&a.symbol);
    let deg = heavy_degree(mol, atom_idx);
    let h = a.hydrogen.max(0) as u32;
    let chg = a.charge as i32 as u32; // wrap negatives
    let isotope: u32 = 0;
    let ring: u32 = if in_ring[atom_idx] { 1 } else { 0 };

    let mut seed: u32 = 0;
    for v in [deg, an, h, chg, isotope, ring] {
        seed = hash_combine(seed, v);
    }
    seed
}

fn bond_invariant(order: BondOrder) -> u32 {
    match order {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}

/// Compute the set of Morgan bit invariants for `mol` at `radius`.
/// Returns the raw u32 invariants emitted across all rounds (one per atom-layer
/// that produced a previously-unseen environment).
pub fn morgan_invariants(mol: &Molecule, radius: u32) -> Vec<u32> {
    let n = mol.atoms.len();
    if n == 0 {
        return Vec::new();
    }

    let in_ring = compute_in_ring(mol);
    let mut current: Vec<u32> = (0..n)
        .map(|i| connectivity_invariant(mol, &in_ring, i))
        .collect();

    // Round 0: emit each atom's initial invariant.
    let mut result: Vec<u32> = current.clone();

    // For each atom, the set of bonds covered by its current neighborhood.
    let mut atom_nbhd: Vec<Vec<bool>> = vec![vec![false; mol.bonds.len()]; n];
    let mut seen_nbhds: std::collections::HashSet<Vec<bool>> = std::collections::HashSet::new();
    seen_nbhds.reserve(n * (radius as usize + 1));
    let mut dead = vec![false; n];

    // Precompute bond indices per atom
    let mut bonds_of: Vec<Vec<(usize, usize, BondOrder)>> = vec![Vec::new(); n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        bonds_of[b.a].push((bi, b.b, b.order));
        bonds_of[b.b].push((bi, b.a, b.order));
    }

    for layer in 0..radius {
        let mut next = current.clone();
        // Tuple per atom: (round_nbhd, new_invariant, atom_idx)
        let mut this_round: Vec<(Vec<bool>, u32, usize)> = Vec::with_capacity(n);

        for atom_idx in 0..n {
            if dead[atom_idx] { continue; }
            if bonds_of[atom_idx].is_empty() {
                dead[atom_idx] = true;
                continue;
            }

            // Grow neighborhood: previous atom_nbhd ∪ each bond ∪ neighbor's previous nbhd
            let mut round_nbhd = atom_nbhd[atom_idx].clone();
            let mut neighbor_pairs: Vec<(u32, u32)> = Vec::with_capacity(bonds_of[atom_idx].len());

            for &(bi, other, order) in &bonds_of[atom_idx] {
                round_nbhd[bi] = true;
                for b in 0..mol.bonds.len() {
                    if atom_nbhd[other][b] {
                        round_nbhd[b] = true;
                    }
                }
                neighbor_pairs.push((bond_invariant(order), current[other]));
            }

            neighbor_pairs.sort();
            let mut invar = layer;
            invar = hash_combine(invar, current[atom_idx]);
            for (bi, ni) in &neighbor_pairs {
                invar = hash_combine(invar, *bi);
                invar = hash_combine(invar, *ni);
            }

            next[atom_idx] = invar;
            this_round.push((round_nbhd, invar, atom_idx));
        }

        // Sort for deterministic dedup order, then emit those with unseen environments.
        this_round.sort_by(|a, b| a.1.cmp(&b.1).then(a.2.cmp(&b.2)));
        for (nbhd, invar, atom_idx) in &this_round {
            if seen_nbhds.insert(nbhd.clone()) {
                result.push(*invar);
            } else {
                dead[*atom_idx] = true;
            }
        }

        // Promote round neighborhoods and invariants.
        for (nbhd, _, atom_idx) in this_round {
            atom_nbhd[atom_idx] = nbhd;
        }
        current = next;
    }

    result
}

/// Pack Morgan invariants into a fixed-width bit vector of `n_bits` bits
/// (rounded up to the nearest byte; little-endian bit order within each byte).
pub fn morgan_bits(mol: &Molecule, radius: u32, n_bits: u32) -> Vec<u8> {
    let n_bytes = ((n_bits + 7) / 8) as usize;
    let mut bits = vec![0u8; n_bytes];
    if n_bits == 0 {
        return bits;
    }
    for inv in morgan_invariants(mol, radius) {
        let bit = (inv % n_bits) as usize;
        bits[bit / 8] |= 1 << (bit % 8);
    }
    bits
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn methane_round0_only() {
        let mol = parse("C").unwrap();
        let invs = morgan_invariants(&mol, 0);
        assert_eq!(invs.len(), 1, "1 heavy atom, radius 0 → 1 invariant");
    }

    #[test]
    fn benzene_has_ring_invariants() {
        let mol = parse("c1ccccc1").unwrap();
        let in_ring = compute_in_ring(&mol);
        assert!(in_ring.iter().all(|&x| x), "all benzene atoms in ring");
    }

    #[test]
    fn ethanol_chain_no_ring() {
        let mol = parse("CCO").unwrap();
        let in_ring = compute_in_ring(&mol);
        assert!(in_ring.iter().all(|&x| !x), "no ring in CCO");
    }

    #[test]
    fn cyclohexane_in_ring() {
        let mol = parse("C1CCCCC1").unwrap();
        let in_ring = compute_in_ring(&mol);
        assert_eq!(in_ring.iter().filter(|&&x| x).count(), 6);
    }

    #[test]
    fn same_smiles_same_bits() {
        let a = morgan_bits(&parse("CCO").unwrap(), 2, 2048);
        let b = morgan_bits(&parse("CCO").unwrap(), 2, 2048);
        assert_eq!(a, b);
    }

    #[test]
    fn different_smiles_different_bits() {
        let a = morgan_bits(&parse("CCO").unwrap(), 2, 2048);
        let b = morgan_bits(&parse("c1ccccc1").unwrap(), 2, 2048);
        assert_ne!(a, b);
    }

    #[test]
    fn bit_count_grows_with_radius() {
        let mol = parse("c1ccccc1O").unwrap(); // phenol
        let r0_pop = morgan_bits(&mol, 0, 2048).iter().map(|b| b.count_ones()).sum::<u32>();
        let r2_pop = morgan_bits(&mol, 2, 2048).iter().map(|b| b.count_ones()).sum::<u32>();
        assert!(r2_pop >= r0_pop, "radius-2 should set at least as many bits as radius-0");
    }

    #[test]
    fn bits_are_2048() {
        let bits = morgan_bits(&parse("CCO").unwrap(), 2, 2048);
        assert_eq!(bits.len(), 256);
    }

    #[test]
    fn aspirin_radius2_sets_many_bits() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        let bits = morgan_bits(&mol, 2, 2048);
        let pop: u32 = bits.iter().map(|b| b.count_ones()).sum();
        assert!(pop >= 10, "aspirin ECFP4 should set many bits, got {}", pop);
    }

    #[test]
    fn empty_radius_still_returns_correct_size() {
        let bits = morgan_bits(&parse("C").unwrap(), 0, 512);
        assert_eq!(bits.len(), 64);
        assert!(bits.iter().any(|b| *b != 0));
    }

    #[test]
    fn caffeine_radius2_distinct_from_benzene() {
        let caffeine = morgan_bits(&parse("Cn1c(=O)c2c(ncn2C)n(C)c1=O").unwrap(), 2, 2048);
        let benzene  = morgan_bits(&parse("c1ccccc1").unwrap(), 2, 2048);
        assert_ne!(caffeine, benzene);
    }

    #[test]
    fn popcount_report() {
        // Run with: cargo test -p ducksmiles_smiles morgan::tests::popcount_report -- --nocapture
        let mols = [
            ("methane",     "C"),
            ("water",       "O"),
            ("ethanol",     "CCO"),
            ("benzene",     "c1ccccc1"),
            ("phenol",      "c1ccccc1O"),
            ("aspirin",     "CC(=O)Oc1ccccc1C(=O)O"),
            ("caffeine",    "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
            ("cholesterol", "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"),
        ];
        eprintln!("\n=== ECFP4 (radius=2, 2048 bits) popcounts ===");
        for (name, smi) in mols {
            let mol = parse(smi).unwrap_or_else(|| panic!("parse failed for {}", smi));
            let bits = morgan_bits(&mol, 2, 2048);
            let pop: u32 = bits.iter().map(|b| b.count_ones()).sum();
            eprintln!("  {:<12} ({:<2} atoms) → {} bits", name, mol.atoms.len(), pop);
        }
        eprintln!("\n=== ECFP6 (radius=3, 4096 bits) popcounts ===");
        for (name, smi) in mols {
            let mol = parse(smi).unwrap();
            let bits = morgan_bits(&mol, 3, 4096);
            let pop: u32 = bits.iter().map(|b| b.count_ones()).sum();
            eprintln!("  {:<12} → {} bits", name, pop);
        }
        // First 8 bytes of the BLOB for the article's "fp" column display
        eprintln!("\n=== First bytes of ECFP4/2048 BLOB ===");
        for (name, smi) in mols {
            let mol = parse(smi).unwrap();
            let bits = morgan_bits(&mol, 2, 2048);
            let hex: String = bits.iter().take(8).map(|b| format!("\\x{:02x}", b)).collect();
            eprintln!("  {:<12} → {}...", name, hex);
        }
    }
}
