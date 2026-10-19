// LeVerrier-Faddeev recurrence used by RDKit Graphs.CharacteristicPolynomial.
// A is sparse and 0/1. Accumulate neighbours in atom-index order; compute both
// IPC variants once. No eigenvalue approximation or per-molecule reference data.
fn ipc(mol: &Molecule) -> (f64, f64) {
    let n = mol.atoms.len();
    let mut adj = mol.adjacency();
    for row in &mut adj {
        row.sort_by_key(|x| x.0);
    }
    let mut an = vec![0.; n * n];
    for i in 0..n {
        for &(j, _) in &adj[i] {
            an[i * n + j] = 1.;
        }
    }
    let mut coefficients = vec![1.];
    let mut next = vec![0.; n * n];
    for k in 1..=n {
        let c = (1. / k as f64) * (0..n).map(|i| an[i * n + i]).sum::<f64>();
        coefficients.push(c.abs());
        for i in 0..n {
            an[i * n + i] -= c;
        }
        next.fill(0.);
        for i in 0..n {
            for &(u, _) in &adj[i] {
                for j in 0..n {
                    next[i * n + j] += an[u * n + j];
                }
            }
        }
        std::mem::swap(&mut an, &mut next);
    }
    let avg = entropy(&coefficients);
    (coefficients.iter().sum::<f64>() * avg, avg)
}

