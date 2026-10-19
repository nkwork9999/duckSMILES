//! Registered structural counts over the shared symmetric SSSR topology.
use crate::parser::{BondOrder, Molecule, RingInfo};

pub(crate) struct Values {
    amides: usize,
    spiro: usize,
    bridgeheads: usize,
    /// QueryOps::queryIsAtomBridgehead, used by stereochemical eligibility.
    /// This is deliberately distinct from CalcNumBridgeheadAtoms' endpoints.
    pub(crate) bridgehead_atoms: Vec<bool>,
}

impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        Some(match id {
            "NumAmideBonds" => self.amides,
            "NumSpiroAtoms" => self.spiro,
            "NumBridgeheadAtoms" => self.bridgeheads,
            _ => return None,
        } as f64)
    }
}

pub(crate) fn calculate(mol: &Molecule, rings: &RingInfo) -> Values {
    let adj = mol.adjacency();
    // The registered query is C(=[O;!R])N. Ring carbonyl O is excluded;
    // the N connection uses the SMARTS default single-or-aromatic bond.
    let amides = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| a.symbol == "C" && !a.aromatic)
        .map(|(i, _)| {
            let oxygen = adj[i]
                .iter()
                .filter(|&&(j, b)| {
                    let a = &mol.atoms[j];
                    a.symbol == "O"
                        && !a.aromatic
                        && !rings.atom_in_ring[j]
                        && mol.bonds[b].order == BondOrder::Double
                })
                .count();
            let nitrogen = adj[i]
                .iter()
                .filter(|&&(j, b)| {
                    let a = &mol.atoms[j];
                    a.symbol == "N"
                        && !a.aromatic
                        && matches!(mol.bonds[b].order, BondOrder::Single | BondOrder::Aromatic)
                })
                .count();
            oxygen * nitrogen
        })
        .sum();

    let bonds: Vec<Vec<usize>> = rings
        .rings
        .iter()
        .map(|r| {
            let mut r = r.clone();
            r.sort_unstable();
            r.dedup();
            r
        })
        .collect();
    let atoms: Vec<Vec<usize>> = bonds
        .iter()
        .map(|ring| {
            let mut a: Vec<_> = ring
                .iter()
                .flat_map(|&b| [mol.bonds[b].a, mol.bonds[b].b])
                .collect();
            a.sort_unstable();
            a.dedup();
            a
        })
        .collect();
    let mut spiro = vec![false; mol.atoms.len()];
    let mut bridgeheads = vec![false; mol.atoms.len()];
    let mut overlaps = vec![vec![false; bonds.len()]; bonds.len()];
    for i in 0..bonds.len() {
        for j in i + 1..bonds.len() {
            let shared_atoms = intersection(&atoms[i], &atoms[j]);
            if shared_atoms.len() == 1 {
                spiro[shared_atoms[0]] = true;
            }
            let shared_bonds = intersection(&bonds[i], &bonds[j]);
            if shared_bonds.len() > 1 {
                overlaps[i][j] = true;
                overlaps[j][i] = true;
                // Fused rings sharing only one bond do not count. Endpoints
                // are defined within each shared bond subgraph, then deduped
                // across pairs; total molecular degree is not a substitute.
                let mut degrees = vec![0usize; mol.atoms.len()];
                for b in shared_bonds {
                    degrees[mol.bonds[b].a] += 1;
                    degrees[mol.bonds[b].b] += 1;
                }
                for (a, degree) in degrees.into_iter().enumerate() {
                    if degree == 1 {
                        bridgeheads[a] = true;
                    }
                }
            }
        }
    }
    // Pinned QueryOps.cpp: at least three incident ring bonds, and EVERY
    // ring containing this atom must overlap another containing ring by two
    // or more bonds. The descriptor above only counts shared-path endpoints;
    // substituting those endpoint flags admits some invertible ring nitrogens.
    let stereo_bridgeheads = (0..mol.atoms.len())
        .map(|a| {
            if adj[a]
                .iter()
                .filter(|&&(_, b)| rings.bond_in_ring[b])
                .count()
                < 3
            {
                return false;
            }
            let memberships: Vec<_> = atoms
                .iter()
                .enumerate()
                .filter_map(|(i, r)| r.contains(&a).then_some(i))
                .collect();
            !memberships.is_empty()
                && memberships
                    .iter()
                    .all(|&i| memberships.iter().any(|&j| overlaps[i][j]))
        })
        .collect();
    Values {
        amides,
        spiro: spiro.into_iter().filter(|x| *x).count(),
        bridgeheads: bridgeheads.iter().filter(|&&x| x).count(),
        bridgehead_atoms: stereo_bridgeheads,
    }
}

fn intersection(a: &[usize], b: &[usize]) -> Vec<usize> {
    let (mut i, mut j) = (0, 0);
    let mut out = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Less => i += 1,
            std::cmp::Ordering::Greater => j += 1,
            std::cmp::Ordering::Equal => {
                out.push(a[i]);
                i += 1;
                j += 1;
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    #[test]
    fn stereochemical_bridgehead_query_matches_independent_atom_flags() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/stereo_bridgehead_query.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let (mol, rings) = crate::workbench::prepare(f[0]).expect(f[0]);
            let actual = super::calculate(&mol, &rings)
                .bridgehead_atoms
                .iter()
                .map(|&v| if v { "1" } else { "0" })
                .collect::<Vec<_>>()
                .join(",");
            if actual != f[1] {
                differences.push(format!("{}: {actual} vs {}", f[0], f[1]));
            }
            count += 1;
        }
        assert_eq!(count, 658);
        assert!(
            differences.is_empty(),
            "{} query differences, first20:\n{}",
            differences.len(),
            differences
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
