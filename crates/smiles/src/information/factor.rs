//! Exact characteristic coefficients via components and articulation vertices.
//! For components Ci of G-v, qi=P(Ci), pi=P(Ci+v):
//! P(G)=sum_i pi*product(j!=i,qj) -(k-1)*x*product_i qi.
//! Small or biconnected leaves use the independent dense integer recurrence.
use super::Integer;
use std::collections::BTreeMap;

fn one() -> Integer {
    Integer {
        negative: false,
        digits: vec![1],
    }
}

fn product(a: &[Integer], b: &[Integer]) -> Vec<Integer> {
    let mut result = vec![Integer::default(); a.len() + b.len() - 1];
    for (i, left) in a.iter().enumerate().filter(|(_, c)| !c.digits.is_empty()) {
        for (j, right) in b.iter().enumerate().filter(|(_, c)| !c.digits.is_empty()) {
            result[i + j].add_assign(&left.multiply(right));
        }
    }
    result
}

struct Solver<'a> {
    adjacency: &'a [Vec<usize>],
    cutoff: usize,
    cache: BTreeMap<Vec<usize>, Vec<Integer>>,
    remaining_visits: usize,
}

impl Solver<'_> {
    fn components(&self, vertices: &[usize], excluded: Option<usize>) -> Vec<Vec<usize>> {
        let mut unseen = vec![false; self.adjacency.len()];
        for &v in vertices {
            unseen[v] = Some(v) != excluded;
        }
        let mut result = Vec::new();
        for &seed in vertices {
            if !unseen[seed] {
                continue;
            }
            unseen[seed] = false;
            let mut queue = vec![seed];
            let mut head = 0;
            while head < queue.len() {
                let vertex = queue[head];
                head += 1;
                for &neighbor in &self.adjacency[vertex] {
                    if unseen[neighbor] {
                        unseen[neighbor] = false;
                        queue.push(neighbor);
                    }
                }
            }
            queue.sort_unstable();
            result.push(queue);
        }
        result
    }

    fn solve(&mut self, vertices: &[usize]) -> Option<Vec<Integer>> {
        // Count cache hits too: bound all recursive work, including branches
        // that have not yet produced an entry. The caller falls back to dense.
        if self.remaining_visits == 0 {
            return None;
        }
        self.remaining_visits -= 1;
        if let Some(result) = self.cache.get(vertices) {
            return Some(result.clone());
        }
        let n = vertices.len();
        if n == 0 {
            return Some(vec![one()]);
        }
        let parts = self.components(vertices, None);
        let result = if parts.len() > 1 {
            let mut value = vec![one()];
            for part in parts {
                value = product(&value, &self.solve(&part)?);
            }
            value
        } else if let Some(cut) = (n > self.cutoff)
            .then(|| balanced_cut(self.adjacency, vertices))
            .flatten()
        {
            let parts = self.components(vertices, Some(cut));
            let mut deleted = Vec::with_capacity(parts.len());
            let mut arms = Vec::with_capacity(parts.len());
            for part in &parts {
                deleted.push(self.solve(part)?);
                let mut arm = part.clone();
                arm.push(cut);
                arm.sort_unstable();
                arms.push(self.solve(&arm)?);
            }
            let mut prefix = vec![vec![one()]];
            for p in &deleted {
                prefix.push(product(prefix.last()?, p));
            }
            let mut suffix = vec![vec![one()]; parts.len() + 1];
            for i in (0..parts.len()).rev() {
                suffix[i] = product(&deleted[i], &suffix[i + 1]);
            }
            let scale = Integer {
                negative: true,
                digits: vec![(parts.len() - 1) as u32],
            };
            let mut value: Vec<_> = prefix.last()?.iter().map(|c| c.multiply(&scale)).collect();
            value.push(Integer::default()); // multiplication by x, descending coefficients
            for (i, arm) in arms.iter().enumerate() {
                let term = product(arm, &product(&prefix[i], &suffix[i + 1]));
                if term.len() != value.len() {
                    return None;
                }
                for (a, b) in value.iter_mut().zip(term) {
                    a.add_assign(&b);
                }
            }
            value
        } else {
            let mut remap = vec![usize::MAX; self.adjacency.len()];
            for (i, &v) in vertices.iter().enumerate() {
                remap[v] = i;
            }
            let graph: Vec<Vec<_>> = vertices
                .iter()
                .map(|&v| {
                    self.adjacency[v]
                        .iter()
                        .filter_map(|&j| (remap[j] != usize::MAX).then_some(remap[j]))
                        .collect()
                })
                .collect();
            super::dense_coefficients(&graph)?
        };
        if result.len() != n + 1 {
            return None;
        }
        self.cache.insert(vertices.to_vec(), result.clone());
        Some(result)
    }
}

fn balanced_cut(adjacency: &[Vec<usize>], vertices: &[usize]) -> Option<usize> {
    struct Search<'a> {
        adjacency: &'a [Vec<usize>],
        allowed: Vec<bool>,
        discovery: Vec<usize>,
        low: Vec<usize>,
        sizes: Vec<usize>,
        time: usize,
        total: usize,
        best: Option<(usize, usize)>,
    }
    impl Search<'_> {
        fn visit(&mut self, vertex: usize, parent: Option<usize>) {
            self.discovery[vertex] = self.time;
            self.low[vertex] = self.time;
            self.time += 1;
            self.sizes[vertex] = 1;
            let mut separated = 0;
            let mut count = 0;
            let mut maximum = 0;
            for index in 0..self.adjacency[vertex].len() {
                let neighbor = self.adjacency[vertex][index];
                if !self.allowed[neighbor] || Some(neighbor) == parent {
                    continue;
                }
                if self.discovery[neighbor] == usize::MAX {
                    self.visit(neighbor, Some(vertex));
                    self.sizes[vertex] += self.sizes[neighbor];
                    self.low[vertex] = self.low[vertex].min(self.low[neighbor]);
                    if self.low[neighbor] >= self.discovery[vertex] {
                        separated += self.sizes[neighbor];
                        count += 1;
                        maximum = maximum.max(self.sizes[neighbor]);
                    }
                } else {
                    self.low[vertex] = self.low[vertex].min(self.discovery[neighbor]);
                }
            }
            let remaining = self.total - 1 - separated;
            if remaining > 0 {
                count += 1;
                maximum = maximum.max(remaining);
            }
            if count > 1 {
                let candidate = (maximum, vertex);
                if self.best.is_none_or(|best| candidate < best) {
                    self.best = Some(candidate);
                }
            }
        }
    }
    let n = adjacency.len();
    let mut search = Search {
        adjacency,
        allowed: vec![false; n],
        discovery: vec![usize::MAX; n],
        low: vec![0; n],
        sizes: vec![0; n],
        time: 0,
        total: vertices.len(),
        best: None,
    };
    for &v in vertices {
        search.allowed[v] = true;
    }
    search.visit(vertices[0], None);
    // Peeling leaves off a large cyclic core duplicates nearly the whole
    // problem. Require a useful size reduction before applying the identity.
    search
        .best
        .filter(|&(largest, _)| largest * 4 <= vertices.len() * 3)
        .map(|(_, v)| v)
}

/// Caller validates the simple undirected adjacency and workbench budgets.
pub(super) fn coefficients(adjacency: &[Vec<usize>], cutoff: usize) -> Option<Vec<Integer>> {
    coefficients_limited(adjacency, cutoff, adjacency.len() * 4 + 64)
}

pub(super) fn coefficients_limited(
    adjacency: &[Vec<usize>],
    cutoff: usize,
    budget: usize,
) -> Option<Vec<Integer>> {
    let mut solver = Solver {
        adjacency,
        cutoff,
        cache: BTreeMap::new(),
        remaining_visits: budget,
    };
    solver.solve(&(0..adjacency.len()).collect::<Vec<_>>())
}

#[cfg(test)]
mod tests {
    #[test]
    fn unbalanced_ring_with_leaves_declines_an_unproductive_cut() {
        let mut graph = vec![Vec::new(); 96];
        for i in 0..48 {
            let j = (i + 1) % 48;
            graph[i].push(j);
            graph[j].push(i);
            graph[i].push(i + 48);
            graph[i + 48].push(i);
        }
        assert!(
            super::balanced_cut(&graph, &(0..96).collect::<Vec<_>>()).is_none(),
            "peeling pendant atoms from a large cyclic core causes redundant expansions"
        );
        let actual = super::coefficients(&graph, 32).unwrap();
        let expected = super::super::dense_coefficients(&graph).unwrap();
        assert_eq!(
            actual
                .iter()
                .map(super::Integer::decimal)
                .collect::<Vec<_>>(),
            expected
                .iter()
                .map(super::Integer::decimal)
                .collect::<Vec<_>>()
        );
    }

    #[test]
    fn explicit_work_limit_stops_recursive_expansion() {
        let graph: Vec<Vec<usize>> = (0..65usize)
            .map(|i| {
                [i.checked_sub(1), (i + 1 < 65).then_some(i + 1)]
                    .into_iter()
                    .flatten()
                    .collect()
            })
            .collect();
        assert!(super::coefficients_limited(&graph, 3, 1).is_none());
        assert!(super::coefficients_limited(&graph, 3, 1024).is_some());
    }
}
