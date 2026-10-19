//! Shared valence, radicals, conjugation and hybridization for 2D descriptors.
//! Pinned RDKit Atom.cpp/MolOps.cpp/ConjugHybrid.cpp; BSD-3-Clause notice applies.
use crate::parser::{BondOrder, Molecule, RingInfo};
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) enum Hybridization {
    Unspecified,
    S,
    Sp,
    Sp2,
    Sp3,
    Sp2d,
    Sp3d,
    Sp3d2,
}
pub(crate) struct Context {
    pub(crate) z: Vec<usize>,
    pub(crate) states: Vec<AtomState>,
    pub(crate) hybridization: Vec<Hybridization>,
    pub(crate) conjugated: Vec<bool>,
    pub(crate) aromatic_bonds: Vec<bool>,
    pub(crate) alpha: Vec<f64>,
}
impl Context {
    pub(crate) fn new(mol: &Molecule) -> Option<Self> {
        Self::with_ring_info(mol, &mol.ring_info())
    }

    pub(crate) fn with_ring_info(mol: &Molecule, rings: &RingInfo) -> Option<Self> {
        let z = atomic_numbers(mol)?;
        let states = atom_states(mol, &z);
        let aromatic_bonds = mol.aromatic_bond_flags(rings);
        let (hybridization, conjugated) = hybridizations(mol, &z, &states, aromatic_bonds.clone());
        let alpha = alpha_for_hybridizations(&z, &hybridization);
        Some(Self {
            z,
            states,
            hybridization,
            conjugated,
            aromatic_bonds,
            alpha,
        })
    }
}
pub(crate) fn atomic_numbers(mol: &Molecule) -> Option<Vec<usize>> {
    mol.atoms
        .iter()
        .map(|a| {
            crate::periodic_data::ELEMENTS
                .iter()
                .position(|e| e.0 == a.symbol)
                .map(|z| z + 1)
        })
        .collect()
}

pub(crate) struct AtomState {
    pub(crate) valence: i32,
    pub(crate) radicals: i32,
    orbitals: i32,
}

/// Atom.cpp integral valence. Aromaticity here includes incident aromatic
/// bonds even when the atom's aromatic flag is false (e.g. acyclic C:O).
pub(crate) fn integral_valence(z: usize, charge: i32, aromatic: bool, mut accum: f64) -> i32 {
    let original = crate::graph_data::VALENCES[z - 1];
    let effective = if original == [-1] {
        z
    } else {
        (z as i32 - charge).clamp(1, 118) as usize
    };
    let allowed = crate::graph_data::VALENCES[effective - 1];
    let default = allowed[0];
    if aromatic && default >= 0 && accum > default as f64 {
        let mut lower = default;
        for &value in allowed {
            if value < 0 || value as f64 > accum {
                break;
            }
            lower = value;
        }
        if accum - lower as f64 <= 1.5 {
            accum = lower as f64;
        }
    }
    (accum + 0.1).round() as i32
}

pub(crate) fn total_valence(mol: &Molecule, index: usize) -> i32 {
    let atom = &mol.atoms[index];
    let mut aromatic = atom.aromatic;
    let mut accum = atom.hydrogen.max(0) as f64;
    for bond in &mol.bonds {
        if bond.a == index || bond.b == index {
            aromatic |= bond.order == BondOrder::Aromatic;
            accum += bond.valence_contribution(index);
        }
    }
    match crate::periodic_data::ELEMENTS
        .iter()
        .position(|e| e.0 == atom.symbol)
    {
        Some(z) => integral_valence(z + 1, atom.charge, aromatic, accum),
        None => -1,
    }
}

pub(crate) fn atom_states(mol: &Molecule, z: &[usize]) -> Vec<AtomState> {
    let adj = mol.adjacency();
    let mut states = Vec::with_capacity(z.len());
    for (i, a) in mol.atoms.iter().enumerate() {
        let outer = crate::descriptor_catalog::OUTER_ELECTRONS[z[i] - 1];
        let vals = crate::graph_data::VALENCES[z[i] - 1];
        let accum = a.hydrogen as f64
            + adj[i]
                .iter()
                .map(|(_, o)| mol.bonds[*o].valence_contribution(i))
                .sum::<f64>();
        let aromatic = a.aromatic
            || adj[i]
                .iter()
                .any(|&(_, b)| mol.bonds[b].order == BondOrder::Aromatic);
        let valence = integral_valence(z[i], a.charge, aromatic, accum);
        let radicals = if !a.in_bracket {
            0
        } else if vals == [-1] {
            if adj[i].is_empty() {
                (outer - a.charge).max(0) % 2
            } else {
                0
            }
        } else {
            // assignRadicals sees the explicit Kekulé valence before aromatic
            // bonds are restored. A two-electron aromatic heteroatom has two
            // single ring bonds, not the later 1.5+1.5 bond sum. The adjusted
            // aromatic valence above retains that integral valence.
            let raw = valence;
            let mut r = if z[i] <= 2 { 2 } else { 8 } - outer - raw + a.charge;
            if r < 0 {
                r = 0;
                if vals.len() > 1 {
                    if let Some(v) = vals.iter().find(|&&v| v - raw + a.charge >= 0) {
                        r = v - raw + a.charge;
                    }
                }
            }
            let r2 = outer - raw - a.charge;
            if r2 >= 0 {
                r = r.min(r2);
            }
            r
        };
        let degree = adj[i]
            .iter()
            .filter(|&&(_, b)| mol.bonds[b].order != BondOrder::Dative || mol.bonds[b].b == i)
            .count() as i32
            + a.hydrogen;
        let free = outer - (valence + a.charge);
        let orbitals = if z[i] <= 1 || z[i] >= 89 {
            degree
        } else if valence + outer - a.charge < 8 {
            degree + (free - radicals) / 2 + radicals
        } else {
            degree + free / 2
        };
        states.push(AtomState {
            valence,
            radicals,
            orbitals,
        });
    }
    states
}
fn hybridizations(
    mol: &Molecule,
    z: &[usize],
    states: &[AtomState],
    mut conjugated: Vec<bool>,
) -> (Vec<Hybridization>, Vec<bool>) {
    let adj = mol.adjacency();
    let degree: Vec<i32> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| adj[i].len() as i32 + a.hydrogen)
        .collect();
    let candidates: Vec<bool> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let vals = crate::graph_data::VALENCES[z[i] - 1];
            let outer = crate::descriptor_catalog::OUTER_ELECTRONS[z[i] - 1];
            if a.charge == 0 && vals[0] >= 0 && states[i].valence > vals[0] {
                return false;
            }
            if !(z[i] <= 10 || (outer != 5 && outer != 6) || (outer == 6 && degree[i] < 2)) {
                return false;
            }
            let dv = vals[0];
            if dv <= 1 || degree[i] > 3 {
                return false;
            }
            (dv - degree[i]) + (outer - dv - a.charge).max(0) - states[i].radicals > 0
        })
        .collect();
    // Only aromatic flags seed conjugation. Acyclic AROMATIC bond order alone
    // does not establish it; the candidate/donor rules below still apply.
    for h in 0..z.len() {
        if !candidates[h] || !(2..=3).contains(&degree[h]) {
            continue;
        }
        for (i, b1) in mol.bonds.iter().enumerate() {
            let v = if b1.a == h {
                b1.b
            } else if b1.b == h {
                b1.a
            } else {
                continue;
            };
            if b1.valence_contribution(h) < 1.5 || !candidates[v] {
                continue;
            }
            for (j, b2) in mol.bonds.iter().enumerate() {
                if i == j {
                    continue;
                }
                let v = if b2.a == h {
                    b2.b
                } else if b2.b == h {
                    b2.a
                } else {
                    continue;
                };
                if candidates[v] && degree[v] <= 3 {
                    conjugated[i] = true;
                    conjugated[j] = true;
                }
            }
        }
    }
    let mut has_conj = vec![false; z.len()];
    for (b, &c) in mol.bonds.iter().zip(&conjugated) {
        if c {
            has_conj[b.a] = true;
            has_conj[b.b] = true;
        }
    }
    let hybridizations = (0..z.len())
        .map(|i| {
            use crate::parser::NonTetrahedralGeometry as G;
            if let Some(geometry) = mol.atoms[i].non_tetrahedral_geometry() {
                if (2..=geometry.maximum_degree()).contains(&degree[i]) {
                    return match geometry {
                        G::SquarePlanar => Hybridization::Sp2d,
                        G::TrigonalBipyramidal => Hybridization::Sp3d,
                        G::Octahedral => Hybridization::Sp3d2,
                    };
                }
            }
            if mol.atoms[i].is_tetrahedral() && degree[i] == 4 {
                return Hybridization::Sp3;
            }
            match states[i].orbitals {
                0 | 1 => Hybridization::S,
                2 => Hybridization::Sp,
                3 => Hybridization::Sp2,
                4 => {
                    if degree[i] <= 3 && has_conj[i] {
                        Hybridization::Sp2
                    } else {
                        Hybridization::Sp3
                    }
                }
                5 => Hybridization::Sp3d,
                6 => Hybridization::Sp3d2,
                _ => Hybridization::Unspecified,
            }
        })
        .collect();
    (hybridizations, conjugated)
}

fn alpha_for_hybridizations(z: &[usize], hybrid: &[Hybridization]) -> Vec<f64> {
    z.iter()
        .zip(hybrid)
        .map(|(&z, &h)| {
            let sp = h == Hybridization::Sp;
            let sp2 = h == Hybridization::Sp2;
            match z {
                1 => 0.,
                6 => {
                    if sp {
                        -0.22
                    } else if sp2 {
                        -0.13
                    } else {
                        0.
                    }
                }
                7 => {
                    if sp {
                        -0.29
                    } else if sp2 {
                        -0.20
                    } else {
                        -0.04
                    }
                }
                8 => {
                    if sp2 {
                        -0.20
                    } else {
                        -0.04
                    }
                }
                9 => -0.07,
                15 => {
                    if sp2 {
                        0.30
                    } else {
                        0.43
                    }
                }
                16 => {
                    if sp2 {
                        0.22
                    } else {
                        0.35
                    }
                }
                17 => 0.29,
                35 => 0.48,
                53 => 0.73,
                _ => crate::surface_data::RADII[z - 1] / crate::surface_data::RADII[5] - 1.,
            }
        })
        .collect()
}
#[cfg(test)]
pub(crate) fn alpha_contributions(mol: &Molecule, z: &[usize], states: &[AtomState]) -> Vec<f64> {
    let (hybrid, _) = hybridizations(mol, z, states, mol.aromatic_bond_flags(&mol.ring_info()));
    alpha_for_hybridizations(z, &hybrid)
}
