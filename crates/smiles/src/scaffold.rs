use crate::parser::{Atom, Bond, BondOrder, Molecule};
use std::collections::HashSet;

#[derive(Debug, Clone)]
struct RingSummary {
    index: usize,
    atoms: Vec<usize>,
    bonds: Vec<usize>,
    aromatic: bool,
}

#[derive(Debug, Clone)]
struct RingSystem {
    atoms: Vec<usize>,
    bonds: Vec<usize>,
    rings: Vec<RingSummary>,
}

pub fn murcko_scaffold_mol(mol: &Molecule) -> Option<Molecule> {
    let ring_info = mol.ring_info();
    if !ring_info.atom_in_ring.iter().any(|&in_ring| in_ring) {
        return None;
    }

    let mut keep = mol
        .atoms
        .iter()
        .map(|atom| atom.symbol != "H")
        .collect::<Vec<_>>();

    loop {
        let mut degree = vec![0usize; mol.atoms.len()];
        for bond in &mol.bonds {
            if keep.get(bond.a).copied().unwrap_or(false)
                && keep.get(bond.b).copied().unwrap_or(false)
            {
                degree[bond.a] += 1;
                degree[bond.b] += 1;
            }
        }

        let mut changed = false;
        for idx in 0..keep.len() {
            if keep[idx] && !ring_info.atom_in_ring[idx] && degree[idx] <= 1 {
                keep[idx] = false;
                changed = true;
            }
        }
        if !changed {
            break;
        }
    }

    // A Murcko scaffold keeps atoms that hang off the framework by a double or
    // triple bond, so a cyclohexanone stays a ketone and a sulfonamide linker
    // keeps its oxygens. These are added back only after pruning has settled —
    // adding them earlier would keep whole ester and acetyl side chains alive.
    let framework = keep.clone();
    for bond in &mol.bonds {
        if !matches!(bond.order, BondOrder::Double | BondOrder::Triple) {
            continue;
        }
        for (kept, other) in [(bond.a, bond.b), (bond.b, bond.a)] {
            if framework.get(kept).copied().unwrap_or(false)
                && mol.atoms[other].symbol != "H"
                && !framework[other]
            {
                keep[other] = true;
            }
        }
    }

    subgraph(mol, &keep)
}

pub fn murcko_scaffold_smiles(mol: &Molecule) -> String {
    murcko_scaffold_mol(mol)
        .map(|scaffold| scaffold.canonical_smiles())
        .unwrap_or_default()
}

pub fn generic_scaffold_smiles(mol: &Molecule) -> String {
    let Some(mut scaffold) = murcko_scaffold_mol(mol) else {
        return String::new();
    };
    for atom in &mut scaffold.atoms {
        *atom = Atom {
            symbol: "C".to_string(),
            hydrogen: 0,
            charge: 0,
            aromatic: false,
            in_bracket: false,
            ..Default::default()
        };
    }
    for bond in &mut scaffold.bonds {
        bond.order = BondOrder::Single;
    }
    // Every atom just became a plain carbon, so its hydrogen count has to be
    // re-derived or the writer would have to bracket all of them.
    scaffold.recompute_implicit_hydrogens();
    scaffold.canonical_smiles()
}

pub fn ring_systems_json(mol: &Molecule) -> String {
    let systems = ring_systems(mol);
    let mut out = String::from("[");
    for (i, system) in systems.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        let aromatic = system
            .bonds
            .iter()
            .all(|&bond_idx| mol.bonds[bond_idx].order == BondOrder::Aromatic);

        out.push_str("{\"index\":");
        out.push_str(&(i + 1).to_string());
        out.push_str(",\"num_atoms\":");
        out.push_str(&system.atoms.len().to_string());
        out.push_str(",\"num_bonds\":");
        out.push_str(&system.bonds.len().to_string());
        out.push_str(",\"num_rings\":");
        out.push_str(&system.rings.len().to_string());
        out.push_str(",\"aromatic\":");
        out.push_str(if aromatic { "true" } else { "false" });
        out.push_str(",\"atoms\":");
        push_usize_array_1based(&mut out, &system.atoms);
        out.push_str(",\"bonds\":");
        push_usize_array_1based(&mut out, &system.bonds);
        out.push_str(",\"rings\":[");
        for (j, ring) in system.rings.iter().enumerate() {
            if j > 0 {
                out.push(',');
            }
            out.push_str("{\"index\":");
            out.push_str(&ring.index.to_string());
            out.push_str(",\"size\":");
            out.push_str(&ring.atoms.len().to_string());
            out.push_str(",\"aromatic\":");
            out.push_str(if ring.aromatic { "true" } else { "false" });
            out.push_str(",\"atoms\":");
            push_usize_array_1based(&mut out, &ring.atoms);
            out.push_str(",\"bonds\":");
            push_usize_array_1based(&mut out, &ring.bonds);
            out.push('}');
        }
        out.push_str("]}");
    }
    out.push(']');
    out
}

fn subgraph(mol: &Molecule, keep: &[bool]) -> Option<Molecule> {
    let mut index_map = vec![usize::MAX; mol.atoms.len()];
    let mut next = 0;
    for idx in 0..mol.atoms.len() {
        if keep.get(idx).copied().unwrap_or(false) {
            index_map[idx] = next;
            next += 1;
        }
    }
    if next == 0 {
        return None;
    }

    let mut atoms = Vec::with_capacity(next);
    for (idx, source) in mol.atoms.iter().enumerate() {
        if index_map[idx] == usize::MAX {
            continue;
        }
        let mut atom = source.clone();
        // Every bond that was cut away is replaced by hydrogens, which is what
        // turns a stripped N-methyl nitrogen back into `[nH]` and a stripped
        // sugar carbon back into a plain `C`.
        let lost: i32 = mol
            .bonds
            .iter()
            .filter(|bond| {
                (bond.a == idx && index_map[bond.b] == usize::MAX)
                    || (bond.b == idx && index_map[bond.a] == usize::MAX)
            })
            .map(|bond| match bond.order {
                BondOrder::Single | BondOrder::Aromatic => 1,
                BondOrder::Double => 2,
                BondOrder::Triple => 3,
            })
            .sum();
        if lost > 0 {
            atom.hydrogen += lost;
            // The neighbour set changed, so the recorded tetrahedral order no
            // longer describes this atom and the stereo tag has to go.
            atom.chirality = None;
            atom.nbr_order.clear();
        } else {
            atom.nbr_order = atom
                .nbr_order
                .iter()
                .map(|&slot| {
                    if slot < 0 {
                        slot
                    } else {
                        index_map[slot as usize] as i32
                    }
                })
                .collect();
        }
        atoms.push(atom);
    }

    let mut bonds = Vec::new();
    for bond in &mol.bonds {
        let a = index_map.get(bond.a).copied().unwrap_or(usize::MAX);
        let b = index_map.get(bond.b).copied().unwrap_or(usize::MAX);
        if a != usize::MAX && b != usize::MAX {
            bonds.push(Bond {
                a,
                b,
                order: bond.order,
                direction: 0,
            });
        }
    }

    let bond_count = bonds.len() as i32;
    Some(Molecule {
        atoms,
        bonds,
        bond_count,
    })
}

fn ring_systems(mol: &Molecule) -> Vec<RingSystem> {
    let ring_info = mol.ring_info();
    let mut systems = Vec::new();
    let mut seen = vec![false; mol.atoms.len()];
    let mut ring_adj = vec![Vec::new(); mol.atoms.len()];

    for (bond_idx, bond) in mol.bonds.iter().enumerate() {
        if ring_info
            .bond_in_ring
            .get(bond_idx)
            .copied()
            .unwrap_or(false)
        {
            ring_adj[bond.a].push((bond.b, bond_idx));
            ring_adj[bond.b].push((bond.a, bond_idx));
        }
    }

    for start in 0..mol.atoms.len() {
        if seen[start] || !ring_info.atom_in_ring[start] {
            continue;
        }
        let mut atoms = Vec::new();
        let mut bond_set = HashSet::new();
        let mut stack = vec![start];
        seen[start] = true;
        while let Some(atom_idx) = stack.pop() {
            atoms.push(atom_idx);
            for &(next, bond_idx) in &ring_adj[atom_idx] {
                bond_set.insert(bond_idx);
                if !seen[next] {
                    seen[next] = true;
                    stack.push(next);
                }
            }
        }
        atoms.sort_unstable();
        let mut bonds = bond_set.into_iter().collect::<Vec<_>>();
        bonds.sort_unstable();
        systems.push(RingSystem {
            atoms,
            bonds,
            rings: Vec::new(),
        });
    }

    for (ring_idx, ring_bonds) in ring_info.rings.iter().enumerate() {
        let ring_bond_set = ring_bonds.iter().copied().collect::<HashSet<_>>();
        let mut ring_atoms = HashSet::new();
        for &bond_idx in ring_bonds {
            if let Some(bond) = mol.bonds.get(bond_idx) {
                ring_atoms.insert(bond.a);
                ring_atoms.insert(bond.b);
            }
        }
        let mut atoms = ring_atoms.into_iter().collect::<Vec<_>>();
        atoms.sort_unstable();
        let mut bonds = ring_bonds.clone();
        bonds.sort_unstable();
        let aromatic = bonds
            .iter()
            .all(|&bond_idx| mol.bonds[bond_idx].order == BondOrder::Aromatic);
        let summary = RingSummary {
            index: ring_idx + 1,
            atoms,
            bonds,
            aromatic,
        };
        if let Some(system) = systems.iter_mut().find(|system| {
            ring_bond_set
                .iter()
                .all(|bond_idx| system.bonds.contains(bond_idx))
        }) {
            system.rings.push(summary);
        }
    }

    systems.sort_by(|a, b| a.atoms.cmp(&b.atoms));
    systems
}

fn push_usize_array_1based(out: &mut String, values: &[usize]) {
    out.push('[');
    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&(value + 1).to_string());
    }
    out.push(']');
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn murcko_removes_side_chains() {
        let mol = parse("Cc1ccccc1").unwrap();
        assert_eq!(murcko_scaffold_smiles(&mol), "c1ccccc1");

        let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        assert_eq!(murcko_scaffold_smiles(&aspirin), "c1ccccc1");
    }

    #[test]
    fn murcko_keeps_linkers_between_ring_systems() {
        let mol = parse("c1ccccc1CCc2ccccc2").unwrap();
        let scaffold = murcko_scaffold_mol(&mol).unwrap();
        assert_eq!(scaffold.atoms.len(), 14);
        assert_eq!(scaffold.bonds.len(), 15);
        assert!(!murcko_scaffold_smiles(&mol).is_empty());
    }

    #[test]
    fn acyclic_molecules_have_empty_scaffold() {
        let mol = parse("CCCC").unwrap();
        assert_eq!(murcko_scaffold_smiles(&mol), "");
        assert_eq!(generic_scaffold_smiles(&mol), "");
    }

    #[test]
    fn generic_scaffold_converts_atoms_and_bonds() {
        let benzene = parse("c1ccccc1").unwrap();
        assert_eq!(generic_scaffold_smiles(&benzene), "C1CCCCC1");

        let pyridine = parse("c1ccccn1").unwrap();
        assert_eq!(generic_scaffold_smiles(&pyridine), "C1CCCCC1");
    }

    #[test]
    fn ring_systems_json_reports_rings() {
        let benzene = parse("c1ccccc1").unwrap();
        assert_eq!(
            ring_systems_json(&benzene),
            "[{\"index\":1,\"num_atoms\":6,\"num_bonds\":6,\"num_rings\":1,\"aromatic\":true,\"atoms\":[1,2,3,4,5,6],\"bonds\":[1,2,3,4,5,6],\"rings\":[{\"index\":1,\"size\":6,\"aromatic\":true,\"atoms\":[1,2,3,4,5,6],\"bonds\":[1,2,3,4,5,6]}]}]"
        );

        let naphthalene = parse("c1ccc2ccccc2c1").unwrap();
        let json = ring_systems_json(&naphthalene);
        assert!(json.contains("\"num_atoms\":10"));
        assert!(json.contains("\"num_rings\":2"));
        assert!(json.contains("\"aromatic\":true"));
    }
}
