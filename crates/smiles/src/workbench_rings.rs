//! Prepared ring discovery following RDKit2026.03.6 FindRings.cpp.
//! BSD-3-Clause attribution: RDKIT-DATA-LICENSE.txt. Raw parser API is separate.
use crate::parser::{BondOrder, Molecule, RingInfo};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

type Cycle = Vec<usize>;
type Invariant = Vec<u64>;
const UNSET: usize = usize::MAX;
const MAX_BFS_QUEUE: usize = 200_000;

struct Graph<'a> {
    mol: &'a Molecule,
    adj: Vec<Vec<(usize, usize)>>,
    order: Vec<usize>,
}

impl<'a> Graph<'a> {
    fn new(mol: &'a Molecule) -> Self {
        // CloseMolRings appends closures in numeric-label order, preserving
        // reuse order within a label. Do not reorder the native bond array.
        let mut order: Vec<_> = (0..mol.bonds.len()).collect();
        order.sort_by_key(|&i| match mol.bonds[i].smiles_ring {
            None => (false, 0, i),
            Some((label, _)) => (true, label, i),
        });
        let mut adj = vec![Vec::new(); mol.atoms.len()];
        for &i in &order {
            let b = &mol.bonds[i];
            adj[b.a].push((b.b, i));
            adj[b.b].push((b.a, i));
        }
        Self { mol, adj, order }
    }

    fn invariant(&self, ring: &[usize]) -> Invariant {
        // boost::dynamic_bitset compares the most significant word first.
        let words = self.adj.len().div_ceil(64);
        let mut key = vec![0; words];
        for &a in ring {
            key[words - 1 - a / 64] |= 1u64 << (a % 64);
        }
        key
    }

    fn bonds(&self, ring: &[usize]) -> Option<Vec<usize>> {
        if ring.len() < 3 {
            return None;
        }
        (0..ring.len())
            .map(|i| {
                self.adj[ring[i]]
                    .iter()
                    .find_map(|&(j, b)| (j == ring[(i + 1) % ring.len()]).then_some(b))
            })
            .collect()
    }

    fn reference_endpoints(&self, bond: usize) -> (usize, usize) {
        let b = &self.mol.bonds[bond];
        if b.order != BondOrder::Dative && b.smiles_ring.is_some_and(|(_, reverse)| reverse) {
            (b.b, b.a)
        } else {
            (b.a, b.b)
        }
    }

    fn trim(
        &self,
        a: usize,
        degree: &mut [usize],
        active: &mut [bool],
        changed: &mut VecDeque<usize>,
    ) {
        for &(j, b) in &self.adj[a] {
            if !active[b] {
                continue;
            }
            if degree[j] <= 2 {
                changed.push_back(j);
            }
            active[b] = false;
            degree[j] -= 1;
            degree[a] -= 1;
        }
    }

    fn smallest(&self, root: usize, active: &[bool], forbidden: &[usize]) -> Option<Vec<Cycle>> {
        let n = self.adj.len();
        let mut colors = vec![0u8; n];
        for &a in forbidden {
            colors[a] = 2;
        }
        let mut parents = vec![UNSET; n];
        let mut depths = vec![0; n];
        let mut queue = VecDeque::from([root]);
        let mut size = UNSET;
        let mut result = Vec::new();
        while let Some(a) = queue.pop_front() {
            if queue.len() >= MAX_BFS_QUEUE {
                return None;
            }
            colors[a] = 2;
            let depth = depths[a] + 1;
            if depth > size {
                break;
            }
            for &(j, b) in &self.adj[a] {
                if !active[b] || colors[j] == 2 || parents[a] == j {
                    continue;
                }
                if colors[j] == 0 {
                    parents[j] = a;
                    colors[j] = 1;
                    depths[j] = depth;
                    queue.push_back(j);
                } else {
                    let mut ring = vec![j];
                    let mut p = parents[j];
                    while p != UNSET && p != root {
                        ring.push(p);
                        p = parents[p];
                    }
                    ring.insert(0, a);
                    p = parents[a];
                    while p != UNSET {
                        if ring.contains(&p) {
                            ring.clear();
                            break;
                        }
                        ring.insert(0, p);
                        p = parents[p];
                    }
                    if ring.len() > 1 {
                        if ring.len() > size {
                            return Some(result);
                        }
                        size = ring.len();
                        result.push(ring);
                    }
                }
            }
        }
        Some(result)
    }

    fn d2_roots(&self, fragment: &[usize], degree: &[usize], active: &[bool]) -> Vec<usize> {
        let mut forbidden = vec![false; self.adj.len()];
        let mut result = Vec::new();
        for &a in fragment {
            if degree[a] != 2 || forbidden[a] {
                continue;
            }
            forbidden[a] = true;
            result.push(a);
            let mut stack = vec![a];
            while let Some(i) = stack.pop() {
                for &(j, b) in &self.adj[i] {
                    if active[b] && degree[j] == 2 && !forbidden[j] {
                        forbidden[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        result
    }

    fn add_marked(
        &self,
        ring: Cycle,
        result: &mut Vec<Cycle>,
        invariants: &mut BTreeSet<Invariant>,
        ring_atoms: &mut [bool],
        ring_bonds: &mut [bool],
    ) -> Option<()> {
        if invariants.insert(self.invariant(&ring)) {
            for &a in &ring {
                ring_atoms[a] = true;
            }
            for b in self.bonds(&ring)? {
                ring_bonds[b] = true;
            }
            result.push(ring);
        }
        Some(())
    }

    fn d2(
        &self,
        roots: &[usize],
        result: &mut Vec<Cycle>,
        invariants: &mut BTreeSet<Invariant>,
        degree: &mut [usize],
        active: &mut [bool],
        ring_atoms: &mut [bool],
        ring_bonds: &mut [bool],
    ) -> Option<()> {
        let mut duplicates: BTreeMap<Invariant, Vec<usize>> = BTreeMap::new();
        let mut peers: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for &root in roots {
            let small = self.smallest(root, active, &[])?;
            for ring in &small {
                let key = self.invariant(ring);
                let other_roots = duplicates.entry(key.clone()).or_default();
                if invariants.contains(&key) {
                    for &other in other_roots.iter() {
                        peers.entry(root).or_default().push(other);
                        peers.entry(other).or_default().push(root);
                    }
                } else {
                    self.add_marked(ring.clone(), result, invariants, ring_atoms, ring_bonds)?;
                }
                other_roots.push(root);
            }
            if small.is_empty() {
                let mut changed = VecDeque::from([root]);
                while let Some(a) = changed.pop_front() {
                    self.trim(a, degree, active, &mut changed);
                }
            }
        }
        for duplicate_roots in duplicates.values().filter(|v| v.len() > 1) {
            let mut candidates = Vec::new();
            let mut min_size = UNSET;
            for &root in duplicate_roots {
                let mut degrees_copy = degree.to_vec();
                let mut active_copy = active.to_vec();
                let mut unused = VecDeque::new();
                for &other in peers.get(&root)? {
                    self.trim(other, &mut degrees_copy, &mut active_copy, &mut unused);
                }
                for ring in self.smallest(root, &active_copy, &[])? {
                    min_size = min_size.min(ring.len());
                    candidates.push(ring);
                }
            }
            for ring in candidates.into_iter().filter(|r| r.len() == min_size) {
                if invariants.insert(self.invariant(&ring)) {
                    result.push(ring);
                }
            }
        }
        Some(())
    }

    fn d3(
        &self,
        root: usize,
        result: &mut Vec<Cycle>,
        invariants: &mut BTreeSet<Invariant>,
        active: &[bool],
    ) -> Option<()> {
        let small = self.smallest(root, active, &[])?;
        for r in &small {
            if invariants.insert(self.invariant(r)) {
                result.push(r.clone());
            }
        }
        if small.len() >= 3 {
            return Some(());
        }
        let neighbors: Vec<_> = self.adj[root]
            .iter()
            .filter(|&&(_, b)| active[b])
            .map(|&(j, _)| j)
            .collect();
        if neighbors.len() != 3 {
            return None;
        }
        let mut forbidden = Vec::new();
        if small.len() == 2 {
            forbidden.push(
                *neighbors
                    .iter()
                    .find(|&&j| small[0].contains(&j) && small[1].contains(&j))?,
            );
        } else if small.len() == 1 {
            let missing = neighbors.iter().position(|j| !small[0].contains(j))?;
            let inside: Vec<_> = neighbors
                .iter()
                .enumerate()
                .filter_map(|(i, &j)| (i != missing).then_some(j))
                .collect();
            forbidden.extend([inside[1], inside[0]]);
        }
        for f in forbidden {
            for r in self.smallest(root, active, &[f])? {
                if invariants.insert(self.invariant(&r)) {
                    result.push(r);
                }
            }
        }
        Some(())
    }

    fn connecting(
        &self,
        bond: usize,
        ring_atoms: &[bool],
        invariants: &BTreeSet<Invariant>,
    ) -> Option<Option<Cycle>> {
        let (start, end) = self.reference_endpoints(bond);
        let mut queue = VecDeque::from([vec![start]]);
        while let Some(path) = queue.pop_front() {
            if queue.len() >= MAX_BFS_QUEUE {
                return None;
            }
            let a = *path.last()?;
            for &(j, _) in &self.adj[a] {
                if j == end {
                    if a != start {
                        let mut ring = path.clone();
                        ring.push(j);
                        if !invariants.contains(&self.invariant(&ring)) {
                            return Some(Some(ring));
                        }
                    }
                } else if ring_atoms[j] && !path.contains(&j) {
                    let mut next = path.clone();
                    next.push(j);
                    queue.push_back(next);
                }
            }
        }
        Some(None)
    }

    fn prune(&self, result: &mut Vec<Cycle>, extras: &mut Vec<Cycle>) -> Option<()> {
        result.sort_by_key(Vec::len);
        let bonds: Vec<BTreeSet<_>> = result
            .iter()
            .map(|r| self.bonds(r).map(|b| b.into_iter().collect()))
            .collect::<Option<_>>()?;
        let mut available = vec![true; bonds.len()];
        let mut keep = vec![false; bonds.len()];
        let mut union = BTreeSet::new();
        for i in 0..bonds.len() {
            if bonds[i].is_subset(&union) {
                available[i] = false;
            }
            if !available[i] {
                continue;
            }
            union.extend(&bonds[i]);
            keep[i] = true;
            let mut consider: Vec<_> = (i + 1..bonds.len())
                .filter(|&j| available[j] && bonds[j].len() == bonds[i].len())
                .collect();
            while !consider.is_empty() {
                let best = *consider.iter().max_by_key(|&&j| {
                    (bonds[j].intersection(&union).count(), std::cmp::Reverse(j))
                })?;
                consider.retain(|&j| j != best);
                available[best] = false;
                if !bonds[best].is_subset(&union) {
                    keep[best] = true;
                    union.extend(&bonds[best]);
                }
            }
        }
        let mut kept = Vec::new();
        for (i, ring) in result.drain(..).enumerate() {
            if keep[i] {
                kept.push(ring);
            } else {
                extras.push(ring);
            }
        }
        *result = kept;
        Some(())
    }

    fn fast_dfs(&self) -> Vec<Cycle> {
        fn visit(
            g: &Graph<'_>,
            a: usize,
            from: usize,
            colors: &mut [u8],
            path: &mut Vec<usize>,
            result: &mut Vec<Cycle>,
        ) {
            colors[a] = 1;
            path.push(a);
            for &(j, _) in &g.adj[a] {
                if colors[j] == 0 {
                    if g.adj[j].len() < 2 {
                        colors[j] = 2;
                    } else {
                        visit(g, j, a, colors, path, result);
                    }
                } else if colors[j] == 1 && from != UNSET && j != from {
                    if let Some(begin) = path.iter().position(|&v| v == j) {
                        result.push(path[begin..].iter().rev().copied().collect());
                    }
                }
            }
            colors[a] = 2;
            path.pop();
        }
        let mut colors = vec![0; self.adj.len()];
        let mut result = Vec::new();
        for a in 0..self.adj.len() {
            if colors[a] != 0 {
                continue;
            }
            if self.adj[a].len() < 2 {
                colors[a] = 2;
            } else {
                visit(self, a, UNSET, &mut colors, &mut Vec::new(), &mut result);
            }
        }
        result
    }

    fn discover(&self) -> Option<Vec<Cycle>> {
        let n = self.adj.len();
        let m = self.order.len();
        let mut active: Vec<_> = self
            .mol
            .bonds
            .iter()
            .map(|b| b.order != BondOrder::Dative)
            .collect();
        let mut degree: Vec<_> = self
            .adj
            .iter()
            .map(|r| r.iter().filter(|&&(_, b)| active[b]).count())
            .collect();
        let mut invariants = BTreeSet::new();
        let mut ring_atoms = vec![false; n];
        let mut ring_bonds = vec![false; m];
        let mut result = Vec::new();
        let mut extras = Vec::new();
        for mut fragment in self.mol.components() {
            fragment.sort_unstable();
            if fragment.len() < 3 {
                continue;
            }
            let all_edges = fragment.iter().map(|&a| self.adj[a].len()).sum::<usize>() / 2;
            if all_edges < fragment.len() {
                continue;
            }
            let edge_count = fragment.iter().map(|&a| degree[a]).sum::<usize>() / 2;
            let expected = edge_count as isize - fragment.len() as isize + 1;
            let mut changed: VecDeque<_> = fragment
                .iter()
                .copied()
                .filter(|&a| degree[a] < 2)
                .collect();
            let mut done = vec![false; n];
            let mut n_done = 0;
            let mut local = Vec::new();
            while n_done <= fragment.len() - 3 {
                while let Some(a) = changed.pop_front() {
                    if !done[a] {
                        done[a] = true;
                        n_done += 1;
                        self.trim(a, &mut degree, &mut active, &mut changed);
                    }
                }
                let roots = self.d2_roots(&fragment, &degree, &active);
                if !roots.is_empty() {
                    self.d2(
                        &roots,
                        &mut local,
                        &mut invariants,
                        &mut degree,
                        &mut active,
                        &mut ring_atoms,
                        &mut ring_bonds,
                    )?;
                    for a in roots {
                        done[a] = true;
                        n_done += 1;
                        self.trim(a, &mut degree, &mut active, &mut changed);
                    }
                } else if n_done <= fragment.len() - 3 {
                    let Some(a) = fragment.iter().copied().find(|&a| degree[a] == 3) else {
                        break;
                    };
                    self.d3(a, &mut local, &mut invariants, &active)?;
                    done[a] = true;
                    n_done += 1;
                    self.trim(a, &mut degree, &mut active, &mut changed);
                }
            }
            if (local.len() as isize) < expected {
                let mut dead = vec![false; m];
                loop {
                    // Preserve the pinned global-bond-prefix rescue scan.
                    let Some(b) = self.order.iter().take(edge_count).copied().find(|&b| {
                        let edge = &self.mol.bonds[b];
                        !ring_bonds[b] && !dead[b] && ring_atoms[edge.a] && ring_atoms[edge.b]
                    }) else {
                        break;
                    };
                    if let Some(ring) = self.connecting(b, &ring_atoms, &invariants)? {
                        self.add_marked(
                            ring,
                            &mut local,
                            &mut invariants,
                            &mut ring_atoms,
                            &mut ring_bonds,
                        )?;
                    } else {
                        dead[b] = true;
                    }
                }
                if (local.len() as isize) < expected {
                    return Some(self.fast_dfs());
                }
            }
            if (local.len() as isize) > expected {
                self.prune(&mut local, &mut extras)?;
            }
            result.extend(local);
        }
        let base: Vec<BTreeSet<_>> = result
            .iter()
            .map(|r| self.bonds(r).map(|v| v.into_iter().collect()))
            .collect::<Option<_>>()?;
        let mut counts = vec![0; m];
        for ring in &base {
            for &b in ring {
                counts[b] += 1;
            }
        }
        for extra in extras {
            let bonds: BTreeSet<_> = self.bonds(&extra)?.into_iter().collect();
            if base.iter().any(|ring| {
                ring.len() == bonds.len()
                    && !ring.is_disjoint(&bonds)
                    && ring.iter().all(|&b| counts[b] != 1 || bonds.contains(&b))
            }) {
                result.push(extra);
            }
        }
        Some(result)
    }
}

pub(crate) fn calculate(mol: &Molecule) -> Option<RingInfo> {
    let graph = Graph::new(mol);
    let cycles = graph.discover()?;
    let mut info = RingInfo {
        atom_in_ring: vec![false; mol.atoms.len()],
        bond_in_ring: vec![false; mol.bonds.len()],
        rings: Vec::new(),
    };
    for ring in cycles {
        for &a in &ring {
            info.atom_in_ring[a] = true;
        }
        let bonds = graph.bonds(&ring)?;
        for &b in &bonds {
            info.bond_in_ring[b] = true;
        }
        info.rings.push(bonds);
    }
    Some(info)
}

#[cfg(test)]
mod tests {
    #[test]
    fn logical_bond_iteration_preserves_smiles_closure_provenance() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/smiles_bond_iteration.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let (source, expected) = (fields[0], fields[1]);
            let (mol, _) =
                crate::native_graph::graph_as_written_with_ring_info(source).expect(source);
            let graph = super::Graph::new(&mol);
            let actual = graph
                .order
                .iter()
                .map(|&b| {
                    let (a, b) = graph.reference_endpoints(b);
                    format!("{a}-{b}")
                })
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(actual, expected, "{source}");
            let types = graph
                .order
                .iter()
                .map(|&b| format!("{:?}", mol.bonds[b].order).to_uppercase())
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(types, fields[2], "{source}: bond types");
            count += 1;
        }
        assert_eq!(count, 28);
    }
}
