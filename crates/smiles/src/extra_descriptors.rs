//! Fifty published scalar descriptors absent from RDKit's pinned descriptor API.
//! Definitions/domain: tasks/descriptors/EXTRA50-DESIGN.md. No trained model/3D.
//! Sparse walk DP and BFS share one hydrogen-suppressed, unweighted graph.
use crate::parser::Molecule;

pub(crate) const IDS: [&str; 50] = [
    "ds_MDEC11",
    "ds_MDEC12",
    "ds_MDEC13",
    "ds_MDEC14",
    "ds_MDEC22",
    "ds_MDEC23",
    "ds_MDEC24",
    "ds_MDEC33",
    "ds_MDEC34",
    "ds_MDEC44",
    "ds_MDEN11",
    "ds_MDEN12",
    "ds_MDEN13",
    "ds_MDEN22",
    "ds_MDEN23",
    "ds_MDEN33",
    "ds_MDEO11",
    "ds_MDEO12",
    "ds_MDEO22",
    "ds_GGI1",
    "ds_GGI2",
    "ds_GGI3",
    "ds_GGI4",
    "ds_GGI5",
    "ds_GGI6",
    "ds_GGI7",
    "ds_GGI8",
    "ds_GGI9",
    "ds_GGI10",
    "ds_MWC3",
    "ds_MWC4",
    "ds_MWC5",
    "ds_MWC6",
    "ds_MWC7",
    "ds_MWC8",
    "ds_MWC9",
    "ds_MWC10",
    "ds_SRW3",
    "ds_SRW4",
    "ds_SRW5",
    "ds_SRW6",
    "ds_SRW7",
    "ds_SRW8",
    "ds_SRW9",
    "ds_SRW10",
    "ds_Wiener",
    "ds_Harary",
    "ds_EccentricConnectivity",
    "ds_DegreeDistance",
    "ds_Gutman",
];
pub(crate) struct Values([Option<f64>; 50]);
impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        if !id.starts_with("ds_") {
            return None;
        }
        IDS.iter()
            .position(|&key| key == id)
            .and_then(|i| self.0[i])
    }
}

pub(crate) fn calculate(mol: &Molecule) -> Values {
    let mut remap = vec![usize::MAX; mol.atoms.len()];
    let mut elements = Vec::new();
    for (i, a) in mol.atoms.iter().enumerate() {
        if a.symbol != "H" {
            remap[i] = elements.len();
            elements.push(a.symbol.as_str());
        }
    }
    let n = elements.len();
    if n == 0 {
        return Values([None; 50]);
    }
    let mut adj = vec![Vec::new(); n];
    for bond in &mol.bonds {
        let (a, b) = (remap[bond.a], remap[bond.b]);
        if a != usize::MAX && b != usize::MAX {
            // Dative connectivity is deliberately undirected in these indices.
            adj[a].push(b);
            adj[b].push(a);
        }
    }
    let degree: Vec<_> = adj.iter().map(Vec::len).collect();
    let mut distances = vec![vec![usize::MAX; n]; n];
    for (i, row) in distances.iter_mut().enumerate() {
        row[i] = 0;
        let mut queue = Vec::with_capacity(n);
        queue.push(i);
        let mut head = 0;
        while head < queue.len() {
            let u = queue[head];
            head += 1;
            for &v in &adj[u] {
                if row[v] == usize::MAX {
                    row[v] = row[u] + 1;
                    queue.push(v);
                }
            }
        }
    }
    let mut result = [Some(0.); 50];
    let mut sums = [0.; 19];
    let mut compensation = [0.; 19];
    let mut counts = [0usize; 19];
    let reciprocal_square = |d: usize| {
        if d == 0 || d == usize::MAX {
            0.
        } else {
            1. / (d as f64).powi(2)
        }
    };
    let mut wiener = 0.;
    let mut distance_counts = vec![0usize; n];
    let mut eccentric = 0.;
    let mut degree_distance = 0.;
    let mut gutman = 0.;
    for i in 0..n {
        let ecc = distances[i]
            .iter()
            .copied()
            .filter(|&d| d != usize::MAX)
            .max()
            .unwrap_or(0);
        eccentric += (degree[i] * ecc) as f64;
        for j in i + 1..n {
            let d = distances[i][j];
            if d == usize::MAX {
                continue;
            }
            let df = d as f64;
            wiener += df;
            distance_counts[d] += 1;
            degree_distance += (degree[i] + degree[j]) as f64 * df;
            gutman += (degree[i] * degree[j]) as f64 * df;
            if elements[i] == elements[j] {
                let (base, maximum) = match elements[i] {
                    "C" => (0, 4),
                    "N" => (10, 3),
                    "O" => (16, 2),
                    _ => (0, 0),
                };
                let (a, b) = (degree[i].min(degree[j]), degree[i].max(degree[j]));
                if a > 0 && b <= maximum {
                    let index = base + (1..a).map(|x| maximum + 1 - x).sum::<usize>() + b - a;
                    // Compensated log sum avoids the product overflow of older
                    // MDE implementations and reduces atom-order sensitivity.
                    let y = df.ln() - compensation[index];
                    let t = sums[index] + y;
                    compensation[index] = (t - sums[index]) - y;
                    sums[index] = t;
                    counts[index] += 1;
                }
            }
            if d <= 10 {
                let left: f64 = adj[i]
                    .iter()
                    .map(|&v| reciprocal_square(distances[v][j]))
                    .sum();
                let right: f64 = adj[j]
                    .iter()
                    .map(|&v| reciprocal_square(distances[v][i]))
                    .sum();
                *result[19 + d - 1].as_mut().unwrap() += (left - right).abs();
            }
        }
    }
    for i in 0..19 {
        result[i] = (counts[i] > 0).then(|| counts[i] as f64 * (-sums[i] / counts[i] as f64).exp());
    }
    // A^k from each start vertex. Walks may repeat edges/vertices. Reuse two
    // flat buffers; no dense matrix multiplication or eigensolver is added.
    let mut walks = vec![0.; n * n];
    let mut next = vec![0.; n * n];
    for i in 0..n {
        walks[i * n + i] = 1.;
    }
    for k in 1..=10 {
        next.fill(0.);
        for start in 0..n {
            for (u, neighbors) in adj.iter().enumerate() {
                let count = walks[start * n + u];
                for &v in neighbors {
                    next[start * n + v] += count;
                }
            }
        }
        std::mem::swap(&mut walks, &mut next);
        if k >= 3 {
            let total: f64 = walks.iter().sum();
            let returning: f64 = (0..n).map(|i| walks[i * n + i]).sum();
            result[29 + k - 3] = Some(total.ln_1p());
            result[37 + k - 3] = Some(returning.ln_1p());
        }
    }
    // Count pairs at each distance exactly, then reduce in a fixed order.
    // This avoids accumulating hundreds of thousands of rounded reciprocals
    // and makes Harary independent of the atom enumeration order.
    let mut harary = 0.;
    let mut harary_compensation = 0.;
    for (distance, count) in distance_counts.into_iter().enumerate().skip(1) {
        let y = count as f64 / distance as f64 - harary_compensation;
        let next = harary + y;
        harary_compensation = (next - harary) - y;
        harary = next;
    }
    for (i, v) in [wiener, harary, eccentric, degree_distance, gutman]
        .into_iter()
        .enumerate()
    {
        result[45 + i] = Some(v);
    }
    for value in &mut result {
        *value = value.filter(|v| v.is_finite());
    }
    Values(result)
}

#[cfg(test)]
mod tests {
    #[test]
    fn large_harary_matches_exact_rational_distance_sums() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/large_harary.tsv")
            .lines()
            .filter(|line| !line.starts_with('#'))
        {
            let (source, reference) = line.split_once('\t').unwrap();
            let (mol, _) = crate::workbench::prepare(source).unwrap();
            let expected: f64 = reference.parse().unwrap();
            let actual = super::calculate(&mol).value("ds_Harary").unwrap();
            if (actual - expected).abs() > 1e-9 {
                differences.push(format!("case{count}: {actual} vs {expected}"));
            }
            count += 1;
        }
        assert_eq!(count, 7);
        assert!(differences.is_empty(), "{}", differences.join("\n"));
    }
}
