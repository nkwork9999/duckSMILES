//! Pinned legacy atomic stereocenter flags (RDKit2026.03.6 Chirality.cpp).
//! Counts use _ChiralityPossible, not @ token count or modern FindStereo.
//! BSD-3-Clause notice: RDKIT-DATA-LICENSE.txt.
use crate::{
    electronic::{Context, Hybridization},
    parser::{BondOrder, Molecule, RingInfo},
};
mod sps;
pub(crate) use sps::calculate as spacial_score;

pub(crate) struct Legacy {
    possible: Vec<bool>,
    codes: Vec<Option<u8>>, // 1=S,2=R, used only for legacy rank refinement.
    ranks: Vec<usize>,
    non_tetrahedral: Vec<bool>,
}
impl Legacy {
    /// Legacy cleanup reclassifies the H on a discarded neutral, acyclic
    /// single-bonded carbon stereo tag. Its total H and cached radicals do not
    /// change. Preserve raw atom syntax; do not apply this bounded rule to
    /// ring stereo, other elements or partial aromatic bond valence.
    pub(crate) fn implicit_hydrogen_transfers(&self, mol: &Molecule, rings: &RingInfo) -> usize {
        mol.atoms
            .iter()
            .enumerate()
            .filter(|&(i, a)| {
                a.symbol == "C"
                    && a.in_bracket
                    && a.hydrogen == 1
                    && a.charge == 0
                    && !a.aromatic
                    && !rings.atom_in_ring[i]
                    && a.is_tetrahedral()
                    && !self.possible[i]
                    && mol
                        .bonds
                        .iter()
                        .filter(|b| b.a == i || b.b == i)
                        .all(|b| b.order == BondOrder::Single)
            })
            .count()
    }

    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        Some(match id {
            "NumAtomStereoCenters" => self.possible.iter().filter(|&&b| b).count(),
            "NumUnspecifiedAtomStereoCenters" => self
                .possible
                .iter()
                .enumerate()
                .filter(|(i, p)| **p && self.codes[*i].is_none() && !self.non_tetrahedral[*i])
                .count(),
            _ => return None,
        } as f64)
    }
}

fn dense<K: Ord>(keys: &[K]) -> Vec<usize> {
    let mut order: Vec<_> = (0..keys.len()).collect();
    order.sort_unstable_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut result = vec![0; keys.len()];
    let mut rank = 0;
    for i in 0..order.len() {
        if i > 0 && keys[order[i]] != keys[order[i - 1]] {
            rank += 1;
        }
        result[order[i]] = rank;
    }
    result
}
fn refine(mol: &Molecule, adj: &[Vec<(usize, usize)>], mut ranks: Vec<usize>) -> Vec<usize> {
    for _ in 0..mol.atoms.len() / 2 + 1 {
        let keys: Vec<_> = adj
            .iter()
            .enumerate()
            .map(|(i, neighbors)| {
                let mut nbrs: Vec<_> = neighbors
                    .iter()
                    .map(|&(j, b)| {
                        let weight = if mol.bonds[b].order == BondOrder::Double
                            && mol.atoms[j].symbol == "P"
                            && matches!(adj[j].len(), 3 | 4)
                        {
                            1
                        } else {
                            match mol.bonds[b].order {
                                BondOrder::Single | BondOrder::Dative => 2,
                                BondOrder::Double => 4,
                                BondOrder::Triple => 6,
                                BondOrder::Aromatic => 3,
                            }
                        };
                        (ranks[j] + 1, weight)
                    })
                    .collect();
                nbrs.sort_unstable_by(|a, b| b.0.cmp(&a.0));
                let mut entries = Vec::new();
                for (r, n) in nbrs {
                    entries.extend(std::iter::repeat_n(r, n));
                }
                entries.extend(std::iter::repeat_n(
                    0,
                    mol.atoms[i].hydrogen.max(0) as usize,
                ));
                (ranks[i], entries)
            })
            .collect();
        let next = dense(&keys);
        let unchanged = next.iter().max() == ranks.iter().max();
        ranks = next;
        if unchanged {
            break;
        }
    }
    ranks
}
fn initial_ranks(mol: &Molecule, context: &Context, adj: &[Vec<(usize, usize)>]) -> Vec<usize> {
    let keys: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mut mass = a.isotope.filter(|&v| v != 0).map_or(0, |v| {
                let delta =
                    v as i32 - crate::periodic_data::ELEMENTS[context.z[i] - 1].2.round() as i32;
                if delta >= 0 {
                    delta + 1
                } else {
                    delta
                }
            }) + 512;
            mass = if mass < 0 { 0 } else { mass % 1024 };
            let map = a.atom_map.map_or(0, |v| (u64::from(v) + 1) % 1024);
            (((context.z[i] as u64 % 128) << 10 | mass as u64) << 10) | map
        })
        .collect();
    refine(mol, adj, dense(&keys))
}
fn parity(reference: &[i32], target: &[i32]) -> Option<bool> {
    if reference.len() != target.len() {
        return None;
    }
    let order: Option<Vec<_>> = target
        .iter()
        .map(|v| reference.iter().position(|x| x == v))
        .collect();
    let order = order?;
    Some(
        (0..order.len())
            .map(|i| {
                (i + 1..order.len())
                    .filter(|&j| order[i] > order[j])
                    .count()
            })
            .sum::<usize>()
            % 2
            == 1,
    )
}
fn assign_code(
    mol: &Molecule,
    i: usize,
    adj: &[Vec<(usize, usize)>],
    ranks: &[usize],
) -> Option<u8> {
    let a = &mol.atoms[i];
    let tag = match a.chirality.as_deref()? {
        "@" => 1,
        "@@" => 2,
        _ => return None,
    };
    let include_h = adj[i].len() == 3 && a.hydrogen == 1;
    let reference: Vec<_> = a
        .nbr_order
        .iter()
        .copied()
        .filter(|&v| v >= 0 || (include_h && v == -1))
        .collect();
    let mut neighbors: Vec<_> = adj[i].iter().map(|&(v, _)| v).collect();
    neighbors.sort_unstable_by_key(|&j| ranks[j]);
    let mut target: Vec<i32> = neighbors.into_iter().map(|j| j as i32).collect();
    if include_h {
        target.insert(0, -1);
    }
    Some(if parity(&reference, &target)? {
        3 - tag
    } else {
        tag
    })
}
fn bond_stereo(
    mol: &Molecule,
    bi: usize,
    adj: &[Vec<(usize, usize)>],
    ranks: &[usize],
    rings: &RingInfo,
) -> Option<u8> {
    let bond = &mol.bonds[bi];
    if bond.order != BondOrder::Double || rings.rings.iter().any(|r| r.len() < 8 && r.contains(&bi))
    {
        return None;
    }
    let direction = |i: usize| {
        if !matches!(adj[i].len(), 2 | 3) {
            return None;
        }
        let mut other: Vec<_> = adj[i]
            .iter()
            .filter(|&&(_, b)| b != bi)
            .map(|&(j, b)| {
                let edge = &mol.bonds[b];
                (
                    j,
                    if edge.a == i {
                        edge.direction
                    } else {
                        -edge.direction
                    },
                )
            })
            .collect();
        if other.iter().all(|&(_, d)| d == 0) {
            return None;
        }
        if other.len() == 2 {
            if ranks[other[0].0] == ranks[other[1].0] || other[0].1 == other[1].1 {
                return None;
            }
            if other[0].1 == 0 {
                other[0].1 = -other[1].1;
            }
            if other[1].1 == 0 {
                other[1].1 = -other[0].1;
            }
        }
        other
            .iter()
            .max_by_key(|&&(j, _)| ranks[j])
            .map(|&(_, d)| d)
    };
    Some(if direction(bond.a)? == direction(bond.b)? {
        2
    } else {
        1
    }) // Z=2,E=1.
}

pub(crate) fn legacy(
    mol: &Molecule,
    rings: &RingInfo,
    context: &Context,
    bridgehead: &[bool],
) -> Legacy {
    let n = mol.atoms.len();
    let adj = mol.adjacency();
    let eligible: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let degree = adj[i]
                .iter()
                .filter(|&&(_, b)| mol.bonds[b].valence_contribution(i) > 0.)
                .count();
            let total = degree as i32 + a.hydrogen;
            if !(3..=4).contains(&total) || (degree < 3 && !matches!(a.symbol.as_str(), "P" | "As"))
            {
                return false;
            }
            if degree != 3 {
                return true;
            }
            if a.hydrogen == 1 {
                return !adj[i].iter().any(|&(j, _)| {
                    mol.atoms[j].symbol == "H" && mol.atoms[j].isotope.unwrap_or(0) <= 1
                });
            }
            match a.symbol.as_str() {
                "P" | "As" => true,
                "S" | "Se" => {
                    context.states[i].valence == 4
                        || (context.states[i].valence == 3 && a.charge == 1)
                }
                "N" => {
                    context.hybridization[i] == Hybridization::Sp3
                        && !adj[i].iter().any(|&(_, b)| context.conjugated[b])
                        && (bridgehead[i]
                            || rings.rings.iter().any(|r| {
                                r.len() == 3
                                    && r.iter()
                                        .any(|&b| mol.bonds[b].a == i || mol.bonds[b].b == i)
                            }))
                }
                _ => false,
            }
        })
        .collect();
    let mut result = Legacy {
        possible: vec![false; n],
        codes: vec![None; n],
        ranks: initial_ranks(mol, context, &adj),
        non_tetrahedral: mol
            .atoms
            .iter()
            .map(|a| a.non_tetrahedral_geometry().is_some())
            .collect(),
    };
    if !eligible.iter().any(|&b| b) {
        return result;
    }
    let mut ranks = result.ranks.clone();
    let mut bonds = vec![None; mol.bonds.len()];
    for _ in 0..=n + mol.bonds.len() {
        let mut changed = false;
        let mut pending = false;
        for i in 0..n {
            if !eligible[i] || result.codes[i].is_some() {
                continue;
            }
            let mut seen = Vec::new();
            let unique = adj[i].iter().all(|&(j, b)| {
                if mol.bonds[b].valence_contribution(i) == 0. {
                    return true;
                }
                if seen.contains(&ranks[j]) {
                    false
                } else {
                    seen.push(ranks[j]);
                    true
                }
            });
            if unique {
                result.possible[i] = true;
                result.codes[i] = assign_code(mol, i, &adj, &ranks);
                changed |= result.codes[i].is_some();
            }
            pending |= result.codes[i].is_none();
        }
        for (bi, known) in bonds.iter_mut().enumerate() {
            if known.is_some() {
                continue;
            }
            let b = &mol.bonds[bi];
            if b.order != BondOrder::Double
                || !matches!(adj[b.a].len(), 2 | 3)
                || !matches!(adj[b.b].len(), 2 | 3)
                || rings.rings.iter().any(|r| r.len() < 8 && r.contains(&bi))
            {
                continue;
            }
            *known = bond_stereo(mol, bi, &adj, &ranks, rings);
            changed |= known.is_some();
            pending |= known.is_none();
        }
        if !changed || !pending {
            break;
        }
        let mut factor = 100;
        while factor < n {
            factor *= 10;
        }
        let seeds: Vec<_> = ranks
            .iter()
            .enumerate()
            .map(|(i, &r)| {
                r * factor
                    + result.codes[i].unwrap_or(0) as usize * 10
                    + adj[i]
                        .iter()
                        .map(|&(_, b)| bonds[b].unwrap_or(0) as usize)
                        .sum::<usize>()
            })
            .collect();
        ranks = refine(mol, &adj, dense(&seeds));
    }
    result.ranks = ranks;
    result
}

#[cfg(test)]
mod tests {
    #[test]
    fn per_atom_legacy_possible_flags_and_cip_refinement_match_reference() {
        let mut differences = Vec::new();
        let mut sources = 0;
        for line in include_str!("../tests/fixtures/stereo_legacy_atoms.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let (m, rings) = crate::workbench::prepare(f[0]).expect(f[0]);
            let context = crate::electronic::Context::with_ring_info(&m, &rings).unwrap();
            let structural = crate::structural_descriptors::calculate(&m, &rings);
            let result = super::legacy(&m, &rings, &context, &structural.bridgehead_atoms);
            let possible = result
                .possible
                .iter()
                .map(|&b| if b { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(",");
            let codes = result
                .codes
                .iter()
                .map(|c| match c {
                    Some(1) => "S",
                    Some(2) => "R",
                    None => "-",
                    _ => unreachable!(),
                })
                .collect::<Vec<_>>()
                .join(",");
            if possible != f[1] {
                differences.push(format!("{} atom flags: {possible} vs {}", f[0], f[1]));
            }
            if codes != f[2] {
                differences.push(format!(
                    "{} rank-refinement codes: {codes} vs {}",
                    f[0], f[2]
                ));
            }
            sources += 1;
        }
        assert_eq!(sources, 658);
        assert!(
            differences.is_empty(),
            "{} context differences, first30:\n{}",
            differences.len(),
            differences
                .iter()
                .take(30)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
