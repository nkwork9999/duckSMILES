//! RDKit2026.03.6 SpacialScore/FindStereo pipeline. BSD notice applies.
//! Modern potential atom centers differ from the public legacy counts.
//! Local symmetry classes never individualize ordinary equivalent atoms.
use super::{dense, Legacy};
use crate::{
    electronic::{Context, Hybridization},
    parser::{BondOrder, Molecule, RingInfo},
};

struct Perception {
    atoms: Vec<bool>,
    double_stereo: Vec<bool>,
}

fn eligible(
    mol: &Molecule,
    rings: &RingInfo,
    ctx: &Context,
    bridgehead: &[bool],
    i: usize,
    adj: &[Vec<(usize, usize)>],
) -> bool {
    let a = &mol.atoms[i];
    // FindStereo considers non-tetrahedral geometry for Be and elements >= Mg.
    // Lighter atoms still follow its tetrahedral path even with an extended tag.
    if a.non_tetrahedral_geometry().is_some() && (ctx.z[i] >= 12 || ctx.z[i] == 4) {
        return false;
    }
    let degree = adj[i]
        .iter()
        .filter(|&&(_, b)| mol.bonds[b].valence_contribution(i) > 0.)
        .count();
    // FindStereo's non-tetrahedral fallback reports an unspecified total-
    // degree-four atom as tetrahedral, including e.g. [PH3][Cu+].
    if (ctx.z[i] >= 12 || ctx.z[i] == 4)
        && mol.atoms[i].chirality.is_none()
        && degree as i32 + a.hydrogen == 4
    {
        return true;
    }
    if degree as i32 + a.hydrogen > 4 {
        return false;
    }
    if degree == 4 {
        return true;
    }
    if degree <= 1 {
        return false;
    }
    if matches!(a.symbol.as_str(), "P" | "As") {
        return true;
    }
    if degree != 3 {
        return false;
    }
    if a.hydrogen == 1 {
        return !adj[i]
            .iter()
            .any(|&(j, _)| mol.atoms[j].symbol == "H" && mol.atoms[j].isotope.unwrap_or(0) <= 1);
    }
    match a.symbol.as_str() {
        "S" | "Se" => ctx.states[i].valence == 4 || (ctx.states[i].valence == 3 && a.charge == 1),
        "N" => {
            ctx.hybridization[i] == Hybridization::Sp3
                && !adj[i].iter().any(|&(_, b)| ctx.conjugated[b])
                && (bridgehead[i]
                    || rings.rings.iter().any(|r| {
                        r.len() == 3
                            && r.iter()
                                .any(|&b| mol.bonds[b].a == i || mol.bonds[b].b == i)
                    }))
        }
        _ => false,
    }
}

fn bond_symbol(mol: &Molecule, ctx: &Context, b: usize) -> &'static str {
    if ctx.aromatic_bonds[b] {
        return ":";
    }
    match mol.bonds[b].order {
        BondOrder::Single => "-",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
        BondOrder::Aromatic => ":",
        BondOrder::Dative => "?",
    }
}

fn symmetry(
    mol: &Molecule,
    adj: &[Vec<(usize, usize)>],
    atom_symbols: &[String],
    bond_symbols: &[String],
    ctx: &Context,
) -> Vec<usize> {
    // new_canon AtomCompareFunctor: degree precedes supplied atom symbol;
    // supplied symbol replaces atomic-number/H/charge comparisons.
    let keys: Vec<_> = (0..mol.atoms.len())
        .map(|i| (adj[i].len(), &atom_symbols[i]))
        .collect();
    let mut ranks = dense(&keys);
    for _ in 0..mol.atoms.len() {
        let keys: Vec<_> = (0..mol.atoms.len())
            .map(|i| {
                let mut neighbors: Vec<_> = adj[i]
                    .iter()
                    .map(|&(j, b)| {
                        let kind = if ctx.aromatic_bonds[b] {
                            12
                        } else {
                            match mol.bonds[b].order {
                                BondOrder::Single => 1,
                                BondOrder::Double => 2,
                                BondOrder::Triple => 3,
                                BondOrder::Aromatic => 12,
                                BondOrder::Dative => 17,
                            }
                        };
                        (&bond_symbols[b], kind, ranks[j])
                    })
                    .collect();
                neighbors.sort_unstable_by(|a, b| b.cmp(a));
                (ranks[i], neighbors)
            })
            .collect();
        let next = dense(&keys);
        let stable = next.iter().max() == ranks.iter().max();
        ranks = next;
        if stable {
            break;
        }
    }
    ranks
}

fn ring_atoms(mol: &Molecule, rings: &RingInfo) -> Vec<Vec<usize>> {
    rings
        .rings
        .iter()
        .map(|ring| {
            let first = &mol.bonds[ring[0]];
            let mut result = vec![first.a, first.b];
            while result.len() < ring.len() {
                let i = *result.last().unwrap();
                let previous = result[result.len() - 2];
                let next = ring
                    .iter()
                    .find_map(|&b| {
                        let b = &mol.bonds[b];
                        if b.a == i && b.b != previous {
                            Some(b.b)
                        } else if b.b == i && b.a != previous {
                            Some(b.a)
                        } else {
                            None
                        }
                    })
                    .expect("cycle bond set");
                result.push(next);
            }
            result
        })
        .collect()
}

fn flag_ring_stereo(
    mol: &Molecule,
    rings: &RingInfo,
    cycles: &[Vec<usize>],
    active_a: &[bool],
    active_b: &[bool],
    adj: &[Vec<(usize, usize)>],
    ring_a: &mut [usize],
    ring_b: &mut [usize],
) {
    let mut memberships = vec![0; mol.bonds.len()];
    for ring in &rings.rings {
        for &b in ring {
            memberships[b] += 1;
        }
    }
    let atom_members: Vec<_> = (0..mol.atoms.len())
        .map(|i| cycles.iter().filter(|r| r.contains(&i)).count())
        .collect();
    for (ri, cycle) in cycles.iter().enumerate() {
        let size = cycle.len();
        let bring = &rings.rings[ri];
        let mut found = vec![false; mol.atoms.len()];
        let mut n_here = 0;
        for (ai, &a) in cycle.iter().enumerate() {
            if !active_a[a] {
                continue;
            }
            for divisor in [2, 3] {
                if size % divisor != 0 {
                    continue;
                }
                let step = size / divisor;
                let mut bonds = 0;
                let mut atoms = 0;
                for offset in (step..size).step_by(step) {
                    let j = cycle[(ai + offset) % size];
                    if adj[j]
                        .iter()
                        .any(|&(_, b)| active_b[b] && !bring.contains(&b))
                    {
                        bonds += 1;
                    }
                    if bonds == 0 && active_a[j] {
                        atoms += 1;
                    }
                }
                if bonds == divisor - 1 || atoms == divisor - 1 {
                    n_here += 1 + bonds;
                    for offset in (0..size).step_by(step) {
                        found[cycle[(ai + offset) % size]] = true;
                    }
                }
            }
            if atom_members[a] > 1 {
                let mut previous = a;
                for step in 1..=size.div_ceil(2) {
                    let j = cycle[(ai + step) % size];
                    let b = adj[previous].iter().find(|&&(k, _)| k == j).unwrap().1;
                    if memberships[b] < 2 {
                        break;
                    }
                    if active_a[j] {
                        n_here += 2;
                        found[a] = true;
                        found[j] = true;
                        break;
                    }
                    previous = j;
                }
            }
        }
        if n_here > 1 {
            for &a in cycle {
                if found[a] {
                    ring_a[a] += 1;
                }
            }
            for &b in bring {
                ring_b[b] += 1;
            }
        }
    }
}

fn perceive(
    mol: &Molecule,
    rings: &RingInfo,
    ctx: &Context,
    bridgehead: &[bool],
    legacy: &Legacy,
) -> Perception {
    let n = mol.atoms.len();
    let m = mol.bonds.len();
    let adj = mol.adjacency();
    let cycles = ring_atoms(mol, rings);
    let base_a: Vec<_> = mol
        .atoms
        .iter()
        .map(|a| format!("{}{}{}", a.isotope.unwrap_or(0), a.symbol, a.charge))
        .collect();
    let base_b: Vec<_> = (0..m)
        .map(|b| bond_symbol(mol, ctx, b).to_owned())
        .collect();
    let potential_a: Vec<_> = (0..n)
        .map(|i| eligible(mol, rings, ctx, bridgehead, i, &adj))
        .collect();
    // Legacy FindPotentialStereoBonds ignores ring bonds unless parsed E/Z
    // stereo was already present. STEREOANY contributes to the SPS factor.
    let double_stereo: Vec<_> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(bi, b)| {
            if b.order != BondOrder::Double {
                return false;
            }
            if super::bond_stereo(mol, bi, &adj, &legacy.ranks, rings).is_some() {
                return true;
            }
            !rings.bond_in_ring[bi]
                && [b.a, b.b].iter().all(|&i| {
                    if !matches!(adj[i].len(), 2 | 3) {
                        return false;
                    }
                    let neighbors: Vec<_> = adj[i]
                        .iter()
                        .filter(|&&(_, b)| b != bi)
                        .filter(|&&(_, b)| {
                            matches!(mol.bonds[b].order, BondOrder::Single | BondOrder::Aromatic)
                        })
                        .map(|&(j, _)| j)
                        .collect();
                    !neighbors.is_empty()
                        && (neighbors.len() != 2
                            || legacy.ranks[neighbors[0]] != legacy.ranks[neighbors[1]])
                })
        })
        .collect();
    let potential_b: Vec<_> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(bi, b)| {
            b.order == BondOrder::Double
                && [b.a, b.b].iter().all(|&i| {
                    let total = adj[i].len() as i32 + mol.atoms[i].hydrogen;
                    let h = mol.atoms[i].hydrogen
                        + adj[i]
                            .iter()
                            .filter(|&&(j, _)| mol.atoms[j].symbol == "H")
                            .count() as i32;
                    total > 1 && total < 4 && h < 2
                })
                && !rings.rings.iter().any(|r| r.len() < 8 && r.contains(&bi))
        })
        .collect();
    let mut known_a: Vec<_> = (0..n)
        .map(|i| {
            potential_a[i]
                && mol.atoms[i].is_tetrahedral()
                && ctx.hybridization[i] == Hybridization::Sp3
                && (legacy.possible[i] || rings.atom_in_ring[i])
        })
        .collect();
    let mut possible_a: Vec<_> = (0..n).map(|i| potential_a[i] && !known_a[i]).collect();
    let known_b: Vec<_> = (0..m).map(|b| potential_b[b] && double_stereo[b]).collect();
    let mut possible_b: Vec<_> = (0..m).map(|b| potential_b[b] && !known_b[b]).collect();
    let original_a = possible_a.clone();
    let original_b = possible_b.clone();
    let mut atom_symbols = base_a.clone();
    let mut bond_symbols = base_b.clone();
    for i in 0..n {
        if possible_a[i] {
            atom_symbols[i] = format!("{}_{}", base_a[i], i);
        } else if known_a[i] {
            let ranks: Vec<_> = (0..n).collect();
            atom_symbols[i] = format!(
                "{}_{:?}",
                base_a[i],
                super::assign_code(mol, i, &adj, &ranks)
            );
        }
    }
    for b in 0..m {
        if let Some(code) = super::bond_stereo(mol, b, &adj, &legacy.ranks, rings) {
            bond_symbols[b] = format!("{}_{code}", base_b[b]);
        } else if known_b[b] || possible_b[b] {
            bond_symbols[b] = format!("{}_{}", base_b[b], b);
        }
    }
    let mut ring_a = vec![0; n];
    let mut ring_b = vec![0; m];
    let mut fixed = vec![false; n];
    let mut output = vec![false; n];
    for pass in 0..2 {
        let active_a: Vec<_> = (0..n).map(|i| known_a[i] || possible_a[i]).collect();
        let active_b: Vec<_> = (0..m).map(|b| known_b[b] || possible_b[b]).collect();
        flag_ring_stereo(
            mol,
            rings,
            &cycles,
            &active_a,
            &active_b,
            &adj,
            &mut ring_a,
            &mut ring_b,
        );
        for _ in 0..=n + m + 1 {
            let ranks = symmetry(mol, &adj, &atom_symbols, &bond_symbols, ctx);
            let mut changed = false;
            output.fill(false);
            for i in 0..n {
                if !known_a[i] && !possible_a[i] {
                    continue;
                }
                if fixed[i] {
                    output[i] = true;
                    continue;
                }
                let mut seen = Vec::new();
                let unique = adj[i].iter().all(|&(j, b)| {
                    if seen.contains(&ranks[j]) {
                        ring_a[i] > 0 && ring_b[b] > 0
                    } else {
                        seen.push(ranks[j]);
                        true
                    }
                });
                if unique {
                    output[i] = true;
                    if !possible_a[i] {
                        let symbol = format!(
                            "{}_{:?}",
                            base_a[i],
                            super::assign_code(mol, i, &adj, &ranks)
                        );
                        changed |= symbol != atom_symbols[i];
                        atom_symbols[i] = symbol;
                        fixed[i] = true;
                    }
                } else {
                    changed |= possible_a[i];
                    possible_a[i] = false;
                    atom_symbols[i] = base_a[i].clone();
                    if ring_a[i] > 0 {
                        ring_a[i] = 0;
                        changed = true;
                        for (ri, ring) in cycles.iter().enumerate() {
                            for &a in ring {
                                fixed[a] = false;
                            }
                            let count = ring.iter().filter(|&&a| ring_a[a] > 0).count();
                            if count <= 1 {
                                if count == 1 {
                                    let a = *ring.iter().find(|&&a| ring_a[a] > 0).unwrap();
                                    ring_a[a] -= 1;
                                }
                                for &b in &rings.rings[ri] {
                                    ring_b[b] = ring_b[b].saturating_sub(1);
                                }
                            }
                        }
                    }
                }
            }
            for bi in 0..m {
                if !known_b[bi] && !possible_b[bi] {
                    continue;
                }
                let b = &mol.bonds[bi];
                let dupe = [b.a, b.b].iter().any(|&i| {
                    let other: Vec<_> = adj[i]
                        .iter()
                        .filter(|&&(_, b)| b != bi)
                        .map(|&(j, _)| j)
                        .collect();
                    other.len() == 2
                        && ranks[other[0]] == ranks[other[1]]
                        && !cycles.iter().any(|ring| {
                            if ring.len() % 2 != 0
                                || !ring.contains(&other[0])
                                || !ring.contains(&other[1])
                            {
                                return false;
                            }
                            [b.a, b.b].iter().any(|end| {
                                ring.iter().position(|a| a == end).is_some_and(|pos| {
                                    let opposite = ring[(pos + ring.len() / 2) % ring.len()];
                                    known_a[opposite]
                                        || possible_a[opposite]
                                        || (adj[opposite].len() == 3
                                            && adj[opposite].iter().any(|&(j, edge)| {
                                                !ring.contains(&j)
                                                    && (known_b[edge] || possible_b[edge])
                                            }))
                                })
                            })
                        })
                });
                if dupe && possible_b[bi] {
                    possible_b[bi] = false;
                    bond_symbols[bi] = base_b[bi].clone();
                    changed = true;
                }
            }
            if !changed {
                break;
            }
        }
        if pass == 1 || (possible_a == original_a && possible_b == original_b) {
            break;
        }
        possible_a = original_a.clone();
        possible_b = original_b.clone();
        for i in 0..n {
            if !fixed[i] && known_a[i] {
                possible_a[i] = true;
                known_a[i] = false;
            }
            if possible_a[i] {
                atom_symbols[i] = format!("{}_{}", atom_symbols[i], i);
            }
        }
        for b in 0..m {
            if possible_b[b] {
                bond_symbols[b] = format!("{}_{}", bond_symbols[b], b);
            }
        }
    }
    Perception {
        atoms: output,
        double_stereo,
    }
}

pub(crate) fn calculate(
    mol: &Molecule,
    rings: &RingInfo,
    ctx: &Context,
    bridgehead: &[bool],
    legacy: &Legacy,
) -> Option<f64> {
    let heavy = ctx.z.iter().filter(|&&z| z > 1).count();
    if heavy == 0 {
        return None;
    }
    let result = perceive(mol, rings, ctx, bridgehead, legacy);
    let adj = mol.adjacency();
    let total: usize = (0..mol.atoms.len())
        .map(|i| {
            let hybrid = match ctx.hybridization[i] {
                Hybridization::Sp => 1,
                Hybridization::Sp2 => 2,
                Hybridization::Sp3 => 3,
                _ => 4,
            };
            let stereo = if result.atoms[i] || adj[i].iter().any(|&(_, b)| result.double_stereo[b])
            {
                2
            } else {
                1
            };
            let ring = if !mol.atoms[i].aromatic && rings.atom_in_ring[i] {
                2
            } else {
                1
            };
            hybrid * stereo * ring * adj[i].len().pow(2)
        })
        .sum();
    Some(total as f64 / heavy as f64)
}

#[cfg(test)]
mod tests {
    #[test]
    fn modern_atom_and_double_bond_masks_match_pinned_pipeline() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../../tests/fixtures/sps_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let (mol, rings) = crate::workbench::prepare(f[0]).expect(f[0]);
            let ctx = crate::electronic::Context::with_ring_info(&mol, &rings).unwrap();
            let structural = crate::structural_descriptors::calculate(&mol, &rings);
            let legacy = super::super::legacy(&mol, &rings, &ctx, &structural.bridgehead_atoms);
            let result = super::perceive(&mol, &rings, &ctx, &structural.bridgehead_atoms, &legacy);
            let adj = mol.adjacency();
            let actual = [
                result
                    .atoms
                    .iter()
                    .map(|&b| usize::from(b))
                    .collect::<Vec<_>>(),
                {
                    let mut bonds: Vec<_> = mol
                        .bonds
                        .iter()
                        .enumerate()
                        .map(|(b, e)| {
                            (
                                (e.a.min(e.b), e.a.max(e.b)),
                                usize::from(result.double_stereo[b]),
                            )
                        })
                        .collect();
                    bonds.sort_unstable_by_key(|&(edge, _)| edge);
                    bonds.into_iter().map(|(_, value)| value).collect()
                },
                ctx.hybridization
                    .iter()
                    .map(|h| match h {
                        crate::electronic::Hybridization::Sp => 1,
                        crate::electronic::Hybridization::Sp2 => 2,
                        crate::electronic::Hybridization::Sp3 => 3,
                        _ => 4,
                    })
                    .collect(),
                mol.atoms
                    .iter()
                    .enumerate()
                    .map(|(i, a)| {
                        if !a.aromatic && rings.atom_in_ring[i] {
                            2
                        } else {
                            1
                        }
                    })
                    .collect(),
                adj.iter().map(Vec::len).collect(),
            ];
            for (i, v) in actual.iter().enumerate() {
                let text = v.iter().map(usize::to_string).collect::<Vec<_>>().join(",");
                if text != f[i + 1] {
                    differences.push(format!("{} factor{} {} vs {}", f[0], i, text, f[i + 1]));
                }
            }
            count += 1;
        }
        assert_eq!(count, 919);
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
