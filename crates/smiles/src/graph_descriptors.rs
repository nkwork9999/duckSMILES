//! Pinned RDKit 2026.03.6 connectivity and GraphDescriptors definitions.
//! BSD-3-Clause upstream algorithms: see RDKIT-DATA-LICENSE.txt.
use crate::parser::{BondOrder, Molecule};
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap};

pub(crate) struct Values(BTreeMap<&'static str, f64>);
impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        self.0.get(id).copied().filter(|v| v.is_finite())
    }
}
fn order(b: BondOrder) -> f64 {
    match b {
        BondOrder::Single => 1.,
        BondOrder::Dative => 17.,
        BondOrder::Double => 2.,
        BondOrder::Triple => 3.,
        BondOrder::Aromatic => 1.5,
    }
}
// Each bond set is one path. Its final atom may repeat any earlier atom;
// traversal stops at that closure. H vertices are excluded (useHs=False).
fn paths(mol: &Molecule) -> [Vec<Vec<usize>>; 4] {
    let mut adj = vec![Vec::new(); mol.atoms.len()];
    for (i, b) in mol.bonds.iter().enumerate() {
        if mol.atoms[b.a].symbol != "H" && mol.atoms[b.b].symbol != "H" {
            adj[b.a].push((b.b, i));
            adj[b.b].push((b.a, i));
        }
    }
    for a in &mut adj {
        a.sort_unstable();
    }
    fn visit(
        adj: &[Vec<(usize, usize)>],
        atoms: &mut Vec<usize>,
        bonds: &mut Vec<usize>,
        seen: &mut [BTreeSet<Vec<usize>>; 4],
        out: &mut [Vec<Vec<usize>>; 4],
    ) {
        let u = *atoms.last().unwrap();
        for &(v, bi) in &adj[u] {
            if bonds.contains(&bi) {
                continue;
            }
            let closed = atoms.contains(&v);
            bonds.push(bi);
            atoms.push(v);
            let len = bonds.len();
            let mut key = bonds.clone();
            key.sort_unstable();
            if seen[len - 1].insert(key) {
                out[len - 1].push(atoms.clone());
            }
            if len < 4 && !closed {
                visit(adj, atoms, bonds, seen, out);
            }
            atoms.pop();
            bonds.pop();
        }
    }
    let mut out = std::array::from_fn(|_| Vec::new());
    let mut seen = std::array::from_fn(|_| BTreeSet::new());
    for i in 0..mol.atoms.len() {
        visit(&adj, &mut vec![i], &mut Vec::new(), &mut seen, &mut out);
    }
    out
}

// Hall/Kier needs sanitized conjugation/hybridization, not the conformer helper
// which only inspects incident multiple bonds. Keep this context local for now.
fn distances(mol: &Molecule, aromatic: &[bool]) -> Vec<Vec<f64>> {
    let adj = mol.adjacency();
    let n = adj.len();
    let mut matrix = vec![vec![1e8; n]; n];
    for (start, row) in matrix.iter_mut().enumerate() {
        let mut d = vec![usize::MAX; n];
        d[start] = 0;
        let mut heap = BinaryHeap::from([Reverse((0usize, start))]);
        while let Some(Reverse((dist, u))) = heap.pop() {
            if dist != d[u] {
                continue;
            }
            for &(v, o) in &adj[u] {
                let w = if aromatic[o] {
                    4
                } else {
                    match mol.bonds[o].order {
                        BondOrder::Single | BondOrder::Dative => 6,
                        BondOrder::Double => 3,
                        BondOrder::Triple => 2,
                        BondOrder::Aromatic => 4,
                    }
                };
                if dist + w < d[v] {
                    d[v] = dist + w;
                    heap.push(Reverse((d[v], v)));
                }
            }
        }
        for (value, d) in row.iter_mut().zip(d) {
            if d != usize::MAX {
                *value = d as f64 / 6.;
            }
        }
    }
    matrix
}
fn entropy(values: &[f64]) -> f64 {
    let total: f64 = values.iter().sum();
    if total == 0. {
        return 0.;
    }
    -values
        .iter()
        .filter(|&&x| x != 0.)
        .map(|&x| {
            let p = x / total;
            p * p.ln()
        })
        .sum::<f64>()
        / std::f64::consts::LN_2
}
fn bertz(mol: &Molecule, z: &[usize], dist: &[Vec<f64>], aromatic: &[bool]) -> f64 {
    let n = z.len();
    if n < 2 {
        return 0.;
    }
    let mut keys = BTreeMap::new();
    let mut classes = Vec::with_capacity(n);
    for row in dist {
        let mut sorted = row.clone();
        sorted.sort_by(f64::total_cmp);
        let key: String = sorted
            .iter()
            .take(100)
            .map(|x| format!("{x:0.4},"))
            .collect();
        let next = keys.len() + 1;
        classes.push(*keys.entry(key).or_insert(next));
    }
    let mut types = BTreeMap::new();
    for &z in z {
        *types.entry(z).or_insert(0.) += 1.;
    }
    let mut connections: BTreeMap<Vec<usize>, f64> = BTreeMap::new();
    let mut adj = mol.adjacency();
    for a in &mut adj {
        a.sort_by_key(|a| a.0);
    }
    for h in 0..n {
        for (i, &(v, o)) in adj[h].iter().enumerate() {
            let b = if aromatic[o] {
                1.5
            } else {
                order(mol.bonds[o].order)
            };
            let (c1, c2) = (classes[h], classes[v]);
            if b > 1. && v > h {
                *connections.entry(vec![c1.min(c2), c1.max(c2)]).or_default() += b * (b - 1.) / 2.;
            }
            for &(w, o2) in &adj[h][i + 1..] {
                let c3 = classes[w];
                *connections
                    .entry(vec![c2.min(c3), c1, c2.max(c3)])
                    .or_default() += b * if aromatic[o2] {
                    1.5
                } else {
                    order(mol.bonds[o2].order)
                };
            }
        }
    }
    if connections.is_empty() {
        connections.insert(vec![0], 1.);
    }
    let c: Vec<_> = connections.into_values().collect();
    let total: f64 = c.iter().sum();
    n as f64 * entropy(&types.into_values().collect::<Vec<_>>())
        + total * (entropy(&c) + total.ln() / std::f64::consts::LN_2)
}
pub(crate) fn calculate(mol: &Molecule, context: &crate::electronic::Context) -> Option<Values> {
    let z = &context.z;
    let adj = mol.adjacency();
    let n = z.len();
    let paths = paths(mol);
    let mut values = BTreeMap::new();
    let delta = |x: f64| if x == 0. { 0. } else { 1. / x.sqrt() };
    let dn: Vec<f64> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            // RDKit's GetTotalNumHs() is unsigned. Its subtraction promotes
            // the outer-electron count to uint32, including underflow for BH4−.
            delta(
                (crate::descriptor_catalog::OUTER_ELECTRONS[z[i] - 1] as u32)
                    .wrapping_sub(a.hydrogen as u32) as f64,
            )
        })
        .collect();
    let dv: Vec<f64> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            if z[i] <= 1 {
                0.
            } else {
                let outer = crate::descriptor_catalog::OUTER_ELECTRONS[z[i] - 1];
                let raw = (outer as u32).wrapping_sub(a.hydrogen as u32) as f64;
                delta(if z[i] <= 10 {
                    raw
                } else {
                    raw / (z[i] as u32).wrapping_sub(outer as u32).wrapping_sub(1) as f64
                })
            }
        })
        .collect();
    let degree: Vec<f64> = adj.iter().map(|v| v.len() as f64).collect();
    values.insert(
        "Chi0",
        degree
            .iter()
            .filter(|&&d| d != 0.)
            .map(|d| 1. / d.sqrt())
            .sum(),
    );
    values.insert(
        "Chi1",
        mol.bonds
            .iter()
            .map(|b| 1. / (degree[b.a] * degree[b.b]).sqrt())
            .sum(),
    );
    for (ds, ids) in [
        (&dn, ["Chi0n", "Chi1n", "Chi2n", "Chi3n", "Chi4n"]),
        (&dv, ["Chi0v", "Chi1v", "Chi2v", "Chi3v", "Chi4v"]),
    ] {
        values.insert(ids[0], ds.iter().sum());
        values.insert(ids[1], mol.bonds.iter().map(|b| ds[b.a] * ds[b.b]).sum());
        for k in 2..=4 {
            values.insert(
                ids[k],
                paths[k - 1]
                    .iter()
                    .map(|p| {
                        let len = if p.first() == p.last() {
                            p.len() - 1
                        } else {
                            p.len()
                        };
                        p[..len].iter().fold(1., |a, &i| a * ds[i])
                    })
                    .sum(),
            );
        }
    }
    let alpha: f64 = context.alpha.iter().sum();
    values.insert("HallKierAlpha", alpha);
    let a = z.iter().filter(|&&z| z > 1).count() as f64;
    let denom1 = mol.bonds.len() as f64 + alpha;
    let k1 = if denom1 == 0. {
        0.
    } else {
        (a + alpha) * (a + alpha - 1.) * (a + alpha - 1.) / (denom1 * denom1)
    };
    let denom2 = (paths[1].len() as f64 + alpha) * (paths[1].len() as f64 + alpha);
    let k2 = if denom2 == 0. {
        0.
    } else {
        (a + alpha - 1.) * (a + alpha - 2.) * (a + alpha - 2.) / denom2
    };
    let denom3 = (paths[2].len() as f64 + alpha) * (paths[2].len() as f64 + alpha);
    let k3 = if denom3 == 0. {
        0.
    } else {
        let start = if a as usize % 2 == 1 {
            a + alpha - 1.
        } else {
            a + alpha - 2.
        };
        start * (a + alpha - 3.) * (a + alpha - 3.) / denom3
    };
    values.insert("Kappa1", k1);
    values.insert("Kappa2", k2);
    values.insert("Kappa3", k3);
    values.insert("Phi", if a == 0. { 0. } else { k1 * k2 / a });
    let dist = distances(mol, &context.aromatic_bonds);
    let sums: Vec<f64> = dist.iter().map(|r| r.iter().sum()).collect();
    let q = mol.bonds.len() as f64;
    let mu = q - n as f64 + 1.;
    let balaban = if mu + 1. == 0. {
        0.
    } else {
        q / (mu + 1.)
            * mol
                .bonds
                .iter()
                .map(|b| 1. / (sums[b.a] * sums[b.b]).sqrt())
                .sum::<f64>()
    };
    values.insert("BalabanJ", balaban);
    values.insert("BertzCT", bertz(mol, &z, &dist, &context.aromatic_bonds));
    Some(Values(values))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::electronic::{alpha_contributions, atom_states, atomic_numbers, Context};
    fn reference_contexts() -> BTreeMap<&'static str, Molecule> {
        use crate::parser::{Atom, Bond};
        include_str!("../tests/fixtures/graph_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<_> = l.split('\t').collect();
                let atoms = f[1]
                    .split(';')
                    .map(|a| {
                        let a: Vec<_> = a.split(',').collect();
                        Atom {
                            symbol: a[0].into(),
                            hydrogen: a[1].parse().unwrap(),
                            charge: a[2].parse().unwrap(),
                            aromatic: a[3] == "1",
                            in_bracket: a[4] == "1",
                            ..Default::default()
                        }
                    })
                    .collect();
                let bonds: Vec<_> = f[2]
                    .split(';')
                    .filter(|s| !s.is_empty())
                    .map(|b| {
                        let b: Vec<_> = b.split(',').collect();
                        let order = match b[2] {
                            "1.0" => BondOrder::Single,
                            "2.0" => BondOrder::Double,
                            "3.0" => BondOrder::Triple,
                            "1.5" => BondOrder::Aromatic,
                            _ => panic!("unexpected bond"),
                        };
                        Bond::new(b[0].parse().unwrap(), b[1].parse().unwrap(), order)
                    })
                    .collect();
                (
                    f[0],
                    Molecule {
                        atoms,
                        bond_count: bonds.len() as i32,
                        bonds,
                    },
                )
            })
            .collect()
    }
    #[test]
    fn all_graph_algorithms_match_independent_sanitized_contexts() {
        let contexts = reference_contexts();
        let mut previous = "";
        let mut actual = None;
        let mut differences = Vec::new();
        for l in include_str!("../tests/fixtures/descriptor_graph.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = l.split('\t').collect();
            let expected: f64 = f[2].parse().unwrap();
            if previous != f[0] {
                actual = calculate(&contexts[f[0]], &Context::new(&contexts[f[0]]).unwrap());
                previous = f[0];
            }
            let value = actual.as_ref().unwrap().value(f[1]);
            let agrees = if expected.is_finite() {
                value.is_some_and(|v| (v - expected).abs() <= 1e-9)
            } else {
                value.is_none()
            };
            if !agrees {
                differences.push(format!("{} {}: {value:?} vs {expected}", f[0], f[1]));
            }
        }
        assert!(
            differences.is_empty(),
            "{} differences:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }
    #[test]
    fn atomic_hall_kier_context_matches_independent_reference() {
        let mut previous = "";
        let mut states = Vec::new();
        let mut alpha = Vec::new();
        let mut differences = Vec::new();
        for line in include_str!("../tests/fixtures/graph_atomic.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let s = f[0];
            let i: usize = f[1].parse().unwrap();
            if previous != s {
                let mol = crate::parser::parse(s).expect(s);
                let z = atomic_numbers(&mol).unwrap();
                states = atom_states(&mol, &z);
                alpha = alpha_contributions(&mol, &z, &states);
                previous = s;
            }
            let valence: i32 = f[3].parse().unwrap();
            let radicals: i32 = f[4].parse().unwrap();
            let expected: f64 = f[5].parse().unwrap();
            if states[i].valence != valence
                || states[i].radicals != radicals
                || (alpha[i] - expected).abs() > 1e-12
            {
                differences.push(format!("{s} atom{i}: valence {} vs {valence}, radicals {} vs {radicals}, alpha {} vs {expected}",states[i].valence,states[i].radicals,alpha[i]));
            }
        }
        assert!(
            differences.is_empty(),
            "{} differences:\n{}",
            differences.len(),
            differences.join("\n")
        );
    }
    #[test]
    fn short_path_counts_match_bond_sets_with_closed_rings_and_h_policy() {
        let mut previous = "";
        let mut actual = std::array::from_fn::<_, 4, _>(|_| Vec::new());
        let mut differences = Vec::new();
        for line in include_str!("../tests/fixtures/graph_paths.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let s = f[0];
            let k: usize = f[1].parse().unwrap();
            let expected: usize = f[2].parse().unwrap();
            if previous != s {
                actual = paths(&crate::parser::parse(s).expect(s));
                previous = s;
            }
            if actual[k - 1].len() != expected {
                differences.push(format!(
                    "{s} length{k}: {} vs {expected}",
                    actual[k - 1].len()
                ));
            }
        }
        assert!(differences.is_empty(), "{}", differences.join("\n"));
    }
}
