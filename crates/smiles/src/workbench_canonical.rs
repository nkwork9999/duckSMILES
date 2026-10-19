//! Pinned new_canon.cpp/.h partition refinement for a restricted raw graph.
//! Copyright2014 Greg Landrum, adapted from Roger Sayle pseudocode.
//! Fast ring DFS: FindRings.cpp; BSD3 notice in RDKIT-DATA-LICENSE.txt.
//! Single bonds and no stereo are required for full ranks. Other molecules
//! retain the existing unsupported result instead of an input-order guess.
use crate::parser::{BondOrder, Molecule};
use std::collections::BTreeSet;

pub(crate) type InitialKey = (u32, usize, usize, u16, i32, u32);

/// Prefix compared before any chirality or neighbor refinement. Distinct
/// initial classes never cross during the reference's partition refinement.
pub(crate) fn initial_keys(mol: &Molecule) -> Option<Vec<InitialKey>> {
    let n = mol.atoms.len();
    if n == 0
        || n > 64
        || mol.atoms.iter().any(|a| {
            !a.in_bracket
                || a.prepared_h_is_implicit
                || a.aromatic
                || a.hydrogen < 0
                || a.atom_map.is_some_and(|m| m > i32::MAX as u32)
        })
        || mol
            .bonds
            .iter()
            .any(|b| b.order != BondOrder::Single || b.direction != 0)
    {
        return None;
    }
    let z = crate::electronic::atomic_numbers(mol)?;
    let adj = mol.adjacency();
    Some(
        mol.atoms
            .iter()
            .enumerate()
            .map(|(i, a)| {
                (
                    a.atom_map.unwrap_or(0),
                    adj[i].len(),
                    z[i],
                    a.isotope.unwrap_or(0),
                    a.hydrogen,
                    a.charge as u32,
                )
            })
            .collect(),
    )
}

fn adjacency(mol: &Molecule) -> Vec<Vec<usize>> {
    // Match CloseMolRings' appended closure order without mutating native bonds.
    let mut order: Vec<_> = (0..mol.bonds.len()).collect();
    order.sort_by_key(|&i| match mol.bonds[i].smiles_ring {
        None => (false, 0, i),
        Some((label, _)) => (true, label, i),
    });
    let mut adj = vec![Vec::new(); mol.atoms.len()];
    for i in order {
        let b = &mol.bonds[i];
        adj[b.a].push(b.b);
        adj[b.b].push(b.a);
    }
    adj
}

fn fast_ring_counts(adj: &[Vec<usize>]) -> Vec<usize> {
    fn visit(
        v: usize,
        parent: Option<usize>,
        adj: &[Vec<usize>],
        colors: &mut [u8],
        path: &mut Vec<usize>,
        rings: &mut [usize],
    ) {
        colors[v] = 1;
        path.push(v);
        for &j in &adj[v] {
            if colors[j] == 0 {
                if adj[j].len() < 2 {
                    colors[j] = 2;
                } else {
                    visit(j, Some(v), adj, colors, path, rings);
                }
            } else if colors[j] == 1 && parent.is_some_and(|p| j != p) {
                for &a in path.iter().rev() {
                    rings[a] += 1;
                    if a == j {
                        break;
                    }
                }
            }
        }
        colors[v] = 2;
        path.pop();
    }
    let mut colors = vec![0; adj.len()];
    let mut rings = vec![0; adj.len()];
    for i in 0..adj.len() {
        if colors[i] != 0 {
            continue;
        }
        if adj[i].len() < 2 {
            colors[i] = 2;
        } else {
            visit(i, None, adj, &mut colors, &mut Vec::new(), &mut rings);
        }
    }
    rings
}

#[derive(PartialEq, Eq, PartialOrd, Ord)]
enum Key {
    Standard(usize, InitialKey, Vec<usize>),
    Special(Vec<i32>, Vec<i32>, Vec<usize>),
}

struct Ranking {
    adj: Vec<Vec<usize>>,
    base: Vec<InitialKey>,
    ring_counts: Vec<usize>,
    rank: Vec<usize>,
    order: Vec<usize>,
    count: Vec<usize>,
    next: Vec<i32>,
    active: i32,
    profiles: Vec<(Vec<i32>, Vec<i32>)>,
    remaining: usize,
}

impl Ranking {
    fn key(&self, v: usize, special: bool) -> Key {
        let mut neighbors: Vec<_> = self.adj[v].iter().map(|&j| self.rank[j]).collect();
        neighbors.sort_by(|a, b| b.cmp(a));
        if special {
            Key::Special(
                self.profiles[v].0.clone(),
                self.profiles[v].1.clone(),
                neighbors,
            )
        } else {
            Key::Standard(self.rank[v], self.base[v], neighbors)
        }
    }

    fn activate(&mut self) {
        self.active = -1;
        self.next.fill(-2);
        let mut position = 0;
        while position < self.order.len() {
            let head = self.order[position];
            debug_assert!(self.count[head] > 0);
            if self.count[head] > 1 {
                self.next[head] = self.active;
                self.active = head as i32;
                position += self.count[head];
            } else {
                position += 1;
            }
        }
    }

    fn touch(&mut self, vertices: &[usize]) {
        let partitions: BTreeSet<_> = vertices
            .iter()
            .flat_map(|&v| &self.adj[v])
            .map(|&j| self.rank[j])
            .collect();
        for position in partitions {
            let head = self.order[position];
            if self.count[head] > 1 && self.next[head] == -2 {
                self.next[head] = self.active;
                self.active = head as i32;
            }
        }
    }

    fn refine(&mut self, special: bool) -> Option<()> {
        while self.active != -1 {
            self.remaining = self.remaining.checked_sub(1)?;
            let head = self.active as usize;
            self.active = self.next[head];
            self.next[head] = -2;
            let (length, offset) = (self.count[head], self.rank[head]);
            let mut keyed: Vec<_> = self.order[offset..offset + length]
                .iter()
                .map(|&v| (v, self.key(v, special)))
                .collect();
            // Stable equivalence-group ordering matches hanoiSort. Recompute
            // keys rather than relying on its changed-neighbor optimization.
            keyed.sort_by(|a, b| a.1.cmp(&b.1));
            let members: Vec<_> = keyed.iter().map(|x| x.0).collect();
            self.order[offset..offset + length].copy_from_slice(&members);
            for &v in &members {
                self.count[v] = 0;
            }
            let mut group = 0;
            for i in 1..=length {
                if i == length || keyed[i].1 != keyed[group].1 {
                    self.count[members[group]] = i - group;
                    group = i;
                }
            }
            let first = self.count[members[0]];
            let mut current = offset;
            for i in first..length {
                let v = members[i];
                if self.count[v] != 0 {
                    current = offset + i;
                }
                self.rank[v] = current;
            }
            self.touch(&members[first..]);
        }
        Some(())
    }

    fn ring_profiles(&mut self) {
        let n = self.order.len();
        for seed in 0..n {
            if self.ring_counts[seed] == 0 {
                continue;
            }
            let mut visited = vec![false; n];
            let mut last = vec![false; n];
            let mut level = vec![seed];
            let mut numbers = Vec::new();
            let mut revisited = Vec::new();
            while !level.is_empty() {
                let mut current = vec![false; n];
                let mut following = Vec::new();
                let mut counts = vec![0; n];
                for v in level {
                    if self.ring_counts[v] == 0 {
                        continue;
                    }
                    last[v] = true;
                    visited[v] = true;
                    for &j in &self.adj[v] {
                        if !visited[j] {
                            visited[j] = true;
                            current[j] = true;
                            following.push(j);
                        }
                    }
                }
                for &v in &following {
                    for &j in &self.adj[v] {
                        if current[j] || last[j] {
                            counts[j] += 1;
                        }
                    }
                }
                let mut nonzero: Vec<_> = counts.into_iter().filter(|&c| c > 0).collect();
                nonzero.sort();
                revisited.extend(nonzero);
                revisited.push(-1);
                numbers.extend([following.len() as i32, -1]);
                last = current;
                level = following;
            }
            self.profiles[seed] = (numbers, revisited);
        }
    }

    fn run(mut self) -> Option<Vec<usize>> {
        self.activate();
        self.refine(false)?;
        // Special chirality refinement is inert in this non-stereo subset.
        let sym: usize = self
            .order
            .iter()
            .filter(|&&v| self.ring_counts[v] > 0 && self.count[v] > 2)
            .map(|&v| self.count[v])
            .sum();
        let ring_atoms = self.ring_counts.iter().filter(|&&c| c > 0).count();
        let branching = (0..self.order.len()).any(|v| self.ring_counts[v] > 1 && self.count[v] > 1);
        if self.count.contains(&0) && ring_atoms > 0 && 2 * sym > ring_atoms && branching {
            self.ring_profiles();
            self.activate();
            self.refine(true)?;
        }
        let mut position = 0;
        while position < self.order.len() {
            let head = self.order[position];
            let old = self.rank[head];
            while self.count[head] > 1 {
                let end = self.rank[head] + self.count[head] - 1;
                let vertex = self.order[end];
                self.rank[vertex] = end;
                self.count[head] -= 1;
                self.count[vertex] = 1;
                if self.adj[vertex].is_empty() {
                    continue;
                }
                self.touch(&[vertex]);
                self.refine(false)?;
            }
            if self.rank[head] == old {
                position += 1;
            }
        }
        Some(self.rank)
    }
}

pub(crate) fn ranks(mol: &Molecule) -> Option<Vec<usize>> {
    let base = initial_keys(mol)?;
    if mol.atoms.iter().any(|a| a.chirality.is_some()) {
        return None;
    }
    let n = base.len();
    let adj = adjacency(mol);
    let ring_counts = fast_ring_counts(&adj);
    let mut count = vec![0; n];
    count[0] = n;
    Ranking {
        adj,
        base,
        ring_counts,
        rank: vec![0; n],
        order: (0..n).collect(),
        count,
        next: vec![-2; n],
        active: -1,
        profiles: vec![(Vec::new(), Vec::new()); n],
        remaining: 8 * n * n + 32,
    }
    .run()
}

#[cfg(test)]
mod tests {
    #[test]
    fn pinned_cluster_and_random_graph_ranks() {
        let mut failures = Vec::new();
        for line in include_str!("../tests/fixtures/cluster_ranks.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let (source, expected) = line.split_once('\t').unwrap();
            let expected: Vec<usize> = expected.split(',').map(|v| v.parse().unwrap()).collect();
            let mol = crate::parser::parse_as_written_graph(source).expect(source);
            let actual = super::ranks(&mol);
            if actual.as_ref() != Some(&expected) {
                failures.push(format!("{source}: {actual:?} != {expected:?}"));
            }
        }
        assert!(
            failures.is_empty(),
            "{} rank differences:\n{}",
            failures.len(),
            failures
                .iter()
                .take(12)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
