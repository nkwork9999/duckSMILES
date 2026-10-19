//! Pinned registered Morgan densities; independent from legacy fingerprints.
//! MorganGenerator.cpp / FingerprintUtil.cpp, Release_2026_03_6 (BSD notice).
use crate::parser::{BondOrder, Molecule, RingInfo};
use std::collections::{BTreeSet, HashSet};

pub(crate) fn calculate(mol: &Molecule, rings: &RingInfo) -> Option<[f64; 3]> {
    let heavy = mol.atoms.iter().filter(|a| a.symbol != "H").count();
    if heavy == 0 {
        return Some([0.; 3]);
    }
    let initial = connectivity_invariants(mol, rings)?;
    let keys = sparse_keys(mol, &initial);
    Some(std::array::from_fn(|i| {
        keys[i + 1].len() as f64 / heavy as f64
    }))
}

// RDKit's portable hash_result_t is uint32_t on every platform.
fn combine(seed: &mut u32, value: u32) {
    *seed ^= value
        .wrapping_add(0x9e3779b9)
        .wrapping_add(*seed << 6)
        .wrapping_add(*seed >> 2);
}

fn connectivity_invariants(mol: &Molecule, rings: &RingInfo) -> Option<Vec<u32>> {
    let z = crate::electronic::atomic_numbers(mol)?;
    let adj = mol.adjacency();
    mol.atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let average = crate::weights::atomic_weight(&a.symbol)?;
            let mass = match a.isotope.filter(|&v| v != 0) {
                Some(v) => crate::weights::isotope_mass(&a.symbol, v)?,
                None => average,
            };
            // The pinned invariant truncates this mass difference toward zero.
            // In particular, it does not round a 13C-average(C) difference to1.
            let mut seed = 0;
            let h = a.hydrogen.max(0) as u32;
            let all_h = h + adj[i]
                .iter()
                .filter(|&&(j, _)| mol.atoms[j].symbol == "H")
                .count() as u32;
            for value in [
                z[i] as u32,
                adj[i].len() as u32 + h,
                all_h,
                a.charge as u32,
                (mass - average) as i32 as u32,
            ] {
                combine(&mut seed, value);
            }
            // Absent ring membership is an absent vector element, not a zero.
            if rings.atom_in_ring[i] {
                combine(&mut seed, 1);
            }
            Some(seed)
        })
        .collect()
}

/// Cumulative unfolded key sets at radius0..3. Chirality is disabled in the
/// registered descriptors; duplicates are bond-set environments, not hashes.
fn sparse_keys(mol: &Molecule, invariants: &[u32]) -> [Vec<u32>; 4] {
    let n = mol.atoms.len();
    let adj = mol.adjacency();
    let words = mol.bonds.len().div_ceil(64);
    let mut previous = vec![vec![0u64; words]; n];
    let mut current = invariants.to_vec();
    let mut emitted = HashSet::<Vec<u64>>::new();
    let mut unique: BTreeSet<u32> = invariants.iter().copied().collect();
    let mut dead = vec![false; n];
    let mut output = std::array::from_fn(|_| Vec::new());
    output[0] = unique.iter().copied().collect();
    for layer in 0..3 {
        let mut round = previous.clone();
        // Dead atoms contribute zero in the following round upstream.
        let mut next = vec![0u32; n];
        let mut candidates = Vec::with_capacity(n);
        for i in 0..n {
            if dead[i] {
                continue;
            }
            if adj[i].is_empty() {
                dead[i] = true;
                continue;
            }
            let mut neighbors = Vec::with_capacity(adj[i].len());
            for &(j, bi) in &adj[i] {
                round[i][bi / 64] |= 1u64 << (bi % 64);
                for (value, &old) in round[i].iter_mut().zip(&previous[j]) {
                    *value |= old;
                }
                let order: u32 = match mol.bonds[bi].order {
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Aromatic => 12,
                    BondOrder::Dative => 17,
                };
                neighbors.push((order, current[j]));
            }
            neighbors.sort_unstable();
            let mut key = layer as u32;
            combine(&mut key, current[i]);
            for (order, other) in neighbors {
                let mut pair = 0;
                combine(&mut pair, order);
                combine(&mut pair, other);
                combine(&mut key, pair);
            }
            next[i] = key;
            candidates.push((round[i].clone(), key, i));
        }
        // dynamic_bitset compares its most significant block first. This
        // order selects which atom survives an identical bond-set environment.
        candidates.sort_unstable_by(|a, b| {
            a.0.iter()
                .rev()
                .cmp(b.0.iter().rev())
                .then_with(|| a.1.cmp(&b.1))
                .then_with(|| a.2.cmp(&b.2))
        });
        for (environment, key, i) in candidates {
            if emitted.insert(environment) {
                unique.insert(key);
            } else {
                dead[i] = true;
            }
        }
        current = next;
        previous = round;
        output[layer + 1] = unique.iter().copied().collect();
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{Atom, Bond, BondOrder};
    use std::collections::BTreeMap;

    fn reference_contexts() -> BTreeMap<&'static str, (Molecule, RingInfo)> {
        include_str!("../tests/fixtures/morgan_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|line| {
                let f: Vec<_> = line.split('\t').collect();
                let atoms = f[1]
                    .split(';')
                    .map(|s| {
                        let a: Vec<_> = s.split(',').collect();
                        Atom {
                            symbol: a[0].into(),
                            hydrogen: a[1].parse().unwrap(),
                            charge: a[2].parse().unwrap(),
                            aromatic: a[3] == "1",
                            in_bracket: a[4] == "1",
                            isotope: match a[5] {
                                "0" => None,
                                s => Some(s.parse().unwrap()),
                            },
                            ..Default::default()
                        }
                    })
                    .collect::<Vec<_>>();
                let bonds = f[2]
                    .split(';')
                    .filter(|s| !s.is_empty())
                    .map(|s| {
                        let b: Vec<_> = s.split(',').collect();
                        let order = match b[2] {
                            "1" => BondOrder::Single,
                            "2" => BondOrder::Double,
                            "3" => BondOrder::Triple,
                            "12" => BondOrder::Aromatic,
                            "17" => BondOrder::Dative,
                            s => panic!("unhandled reference type {s}"),
                        };
                        Bond::new(b[0].parse().unwrap(), b[1].parse().unwrap(), order)
                    })
                    .collect::<Vec<_>>();
                let mut rings = RingInfo {
                    atom_in_ring: vec![false; atoms.len()],
                    bond_in_ring: vec![false; bonds.len()],
                    rings: Vec::new(),
                };
                for r in f[3].split(';').filter(|s| !s.is_empty()) {
                    let indices: Vec<usize> = r.split(',').map(|v| v.parse().unwrap()).collect();
                    for &bi in &indices {
                        rings.bond_in_ring[bi] = true;
                        rings.atom_in_ring[bonds[bi].a] = true;
                        rings.atom_in_ring[bonds[bi].b] = true;
                    }
                    rings.rings.push(indices);
                }
                (
                    f[0],
                    (
                        Molecule {
                            atoms,
                            bond_count: bonds.len() as i32,
                            bonds,
                        },
                        rings,
                    ),
                )
            })
            .collect()
    }

    #[test]
    fn atom_invariants_and_all_sparse_keys_match_independent_reference_graphs() {
        let contexts = reference_contexts();
        let mut errors = Vec::new();
        for line in include_str!("../tests/fixtures/morgan_sparse_keys.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let (mol, rings) = &contexts[f[0]];
            let expected: Vec<u32> = f[1].split(',').map(|v| v.parse().unwrap()).collect();
            let actual = connectivity_invariants(mol, rings)
                .expect("registered connectivity invariants missing");
            if actual != expected {
                errors.push(format!(
                    "{} atom invariants {actual:?} vs {expected:?}",
                    f[0]
                ));
            }
            let keys = sparse_keys(mol, &actual);
            for r in 0..4 {
                let expected: Vec<u32> = f[r + 2]
                    .split(',')
                    .filter(|v| !v.is_empty())
                    .map(|v| v.parse().unwrap())
                    .collect();
                if keys[r] != expected {
                    errors.push(format!(
                        "{} radius{r} keys {:?} vs {expected:?}",
                        f[0], keys[r]
                    ));
                }
            }
        }
        assert!(
            errors.is_empty(),
            "{} errors:\n{}",
            errors.len(),
            errors.into_iter().take(20).collect::<Vec<_>>().join("\n")
        );
    }
}
