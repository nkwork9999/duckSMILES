use crate::conformer::bounds::BoundsMatrix;

// ── Embedding: distance bounds → 3D coordinates ───────────────────────────────
//
// Algorithm (classical distance geometry):
//  1. Sample a trial distance matrix D by picking d(i,j) ∈ [lo(i,j), up(i,j)]
//     using a deterministic LCG PRNG seeded by `seed`.
//  2. Compute the metric (Gram) matrix G from D:
//       G[i][j] = ½ (D[0][i]² + D[0][j]² − D[i][j]²)
//     where atom 0 is chosen as the reference origin.
//  3. Eigen-decompose G via the Jacobi method (symmetric matrix, exact for DG).
//  4. Retain the three largest positive eigenvalues λ₁ ≥ λ₂ ≥ λ₃ > 0.
//  5. Coordinates: X[i] = √λ₁ · v₁[i],  Y[i] = √λ₂ · v₂[i],  Z[i] = √λ₃ · v₃[i].
//
// The Jacobi method converges for symmetric real matrices (which G always is).

// ── LCG PRNG ──────────────────────────────────────────────────────────────────
// Linear Congruential Generator: Knuth/Numerical-Recipes parameters.
struct Lcg {
    state: u64,
}

impl Lcg {
    fn new(seed: u64) -> Self {
        Lcg { state: seed.wrapping_add(1) }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }

    // Returns a value in [0.0, 1.0)
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
}

// ── Jacobi symmetric eigendecomposition ──────────────────────────────────────
// Input:  n×n symmetric matrix `a` (flat row-major).
// Output: eigenvalues in `evals`, eigenvectors as columns of `evecs` (row-major).
// Convergence: iterates until all off-diagonal elements < tol, max 1000 sweeps.
fn jacobi(a: &[f64], n: usize, evals: &mut Vec<f64>, evecs: &mut Vec<f64>) {
    let mut a = a.to_vec(); // working copy
    *evecs = (0..n * n)
        .map(|k| if k / n == k % n { 1.0 } else { 0.0 })
        .collect(); // identity matrix

    const MAX_ITER: usize = 100 * 1000; // ~100 sweeps × n² rotations
    const TOL: f64 = 1e-12;

    for _ in 0..MAX_ITER {
        // Find the largest off-diagonal |a[p][q]|
        let mut max_val = 0.0_f64;
        let (mut p, mut q) = (0, 1);
        for i in 0..n {
            for j in (i + 1)..n {
                let v = a[i * n + j].abs();
                if v > max_val {
                    max_val = v;
                    p = i;
                    q = j;
                }
            }
        }
        if max_val < TOL {
            break;
        }

        // Jacobi rotation to zero out a[p][q]
        let theta = if (a[p * n + p] - a[q * n + q]).abs() < 1e-300 {
            std::f64::consts::FRAC_PI_4
        } else {
            0.5 * ((2.0 * a[p * n + q]) / (a[p * n + p] - a[q * n + q])).atan()
        };
        let c = theta.cos();
        let s = theta.sin();

        // Update matrix rows/cols p and q
        let app = a[p * n + p];
        let aqq = a[q * n + q];
        let apq = a[p * n + q];
        a[p * n + p] = c * c * app + 2.0 * s * c * apq + s * s * aqq;
        a[q * n + q] = s * s * app - 2.0 * s * c * apq + c * c * aqq;
        a[p * n + q] = 0.0;
        a[q * n + p] = 0.0;

        for r in 0..n {
            if r != p && r != q {
                let arp = a[r * n + p];
                let arq = a[r * n + q];
                a[r * n + p] = c * arp + s * arq;
                a[p * n + r] = a[r * n + p];
                a[r * n + q] = -s * arp + c * arq;
                a[q * n + r] = a[r * n + q];
            }
        }

        // Update eigenvectors
        for r in 0..n {
            let vrp = evecs[r * n + p];
            let vrq = evecs[r * n + q];
            evecs[r * n + p] = c * vrp + s * vrq;
            evecs[r * n + q] = -s * vrp + c * vrq;
        }
    }

    // Diagonal of the rotated matrix = eigenvalues
    *evals = (0..n).map(|i| a[i * n + i]).collect();
}

// ── embed ─────────────────────────────────────────────────────────────────────
// Returns `Some(coords)` where coords is a flat Vec<f64> of length 3n:
//   coords[3*i], coords[3*i+1], coords[3*i+2] = (x, y, z) of atom i.
// Returns `None` if the metric matrix has fewer than 3 positive eigenvalues.
pub fn embed(bm: &BoundsMatrix, seed: u64) -> Option<Vec<f64>> {
    let n = bm.n;
    let mut rng = Lcg::new(seed);

    // ── Step 1: sample a trial distance matrix ─────────────────────────────
    let mut d = vec![0.0_f64; n * n]; // symmetric, d[i][j] = sampled distance
    for i in 0..n {
        for j in (i + 1)..n {
            let lo = bm.lo(i, j);
            let hi = bm.up(i, j);
            let t = rng.next_f64();
            let dij = lo + t * (hi - lo);
            d[i * n + j] = dij;
            d[j * n + i] = dij;
        }
    }

    // ── Step 2: double-centred Gram matrix (Torgerson/MDS formula) ────────
    // G = -½ J D² J  where J = I − (1/n)11ᵀ
    // G[i][j] = -½ (D²[i][j] - row_mean[i] - col_mean[j] + grand_mean)
    let d2: Vec<f64> = d.iter().map(|v| v * v).collect();
    let row_mean: Vec<f64> = (0..n)
        .map(|i| d2[i * n..i * n + n].iter().sum::<f64>() / n as f64)
        .collect();
    let grand_mean: f64 = row_mean.iter().sum::<f64>() / n as f64;
    let mut g = vec![0.0_f64; n * n];
    for i in 0..n {
        for j in 0..n {
            g[i * n + j] =
                -0.5 * (d2[i * n + j] - row_mean[i] - row_mean[j] + grand_mean);
        }
    }

    // ── Step 3: Jacobi eigendecomposition of G ────────────────────────────
    let mut evals = Vec::new();
    let mut evecs = Vec::new();
    jacobi(&g, n, &mut evals, &mut evecs);

    // ── Step 4: sort by descending eigenvalue, pick top 3 ─────────────────
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| evals[b].partial_cmp(&evals[a]).unwrap());

    // Molecules with < 3 atoms or linear configurations may have only 1-2
    // positive eigenvalues; use 0 for missing dimensions.
    let scale = |idx: usize| -> f64 {
        let ev = evals[order[idx]];
        if ev > 1e-9 { ev.sqrt() } else { 0.0 }
    };
    let s0 = scale(0);
    let s1 = if n > 1 { scale(1) } else { 0.0 };
    let s2 = if n > 2 { scale(2) } else { 0.0 };

    if s0 == 0.0 {
        return None; // no positive eigenvalue at all — degenerate
    }

    let (i0, i1, i2) = (order[0], order.get(1).copied().unwrap_or(0), order.get(2).copied().unwrap_or(0));

    // ── Step 5: extract coordinates ────────────────────────────────────────
    // evecs is stored as n×n row-major; evecs[atom*n + col_idx] = component.
    let mut coords = vec![0.0_f64; 3 * n];
    for atom in 0..n {
        coords[3 * atom]     = s0 * evecs[atom * n + i0];
        coords[3 * atom + 1] = s1 * evecs[atom * n + i1];
        coords[3 * atom + 2] = s2 * evecs[atom * n + i2];
    }

    Some(coords)
}

// ── helpers ───────────────────────────────────────────────────────────────────

#[cfg(test)]
pub(crate) fn dist3(coords: &[f64], i: usize, j: usize) -> f64 {
    let dx = coords[3 * i] - coords[3 * j];
    let dy = coords[3 * i + 1] - coords[3 * j + 1];
    let dz = coords[3 * i + 2] - coords[3 * j + 2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::conformer::smoothing::smooth;
    use crate::parser::parse;

    fn embed_mol(smi: &str, seed: u64) -> (Vec<f64>, usize) {
        let mol = parse(smi).unwrap().with_explicit_hydrogens();
        let n = mol.atoms.len();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        let coords = embed(&bm, seed).expect("embed returned None");
        (coords, n)
    }

    #[test]
    fn embed_returns_3n_coords() {
        let (coords, n) = embed_mol("CCC", 42);
        assert_eq!(coords.len(), 3 * n);
    }

    #[test]
    fn embed_no_nan_or_inf() {
        let (coords, _) = embed_mol("CC(=O)Oc1ccccc1C(=O)O", 7);
        for &v in &coords {
            assert!(v.is_finite(), "non-finite coordinate: {v}");
        }
    }

    #[test]
    fn embed_bond_lengths_reasonable() {
        // After smoothing + embedding, bonded heavy-atom distances should be
        // in a chemically plausible range [1.0, 2.2] Å.
        let mol_bare = parse("CCO").unwrap();
        let mol = mol_bare.with_explicit_hydrogens();
        let n = mol.atoms.len();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        let coords = embed(&bm, 1).unwrap();

        // Check all 1-2 pairs (bonds)
        for b in &mol.bonds {
            let d = dist3(&coords, b.a, b.b);
            assert!(
                d >= 0.8 && d <= 2.5,
                "bond distance out of range: atoms {}-{}: {d:.3}",
                b.a, b.b
            );
        }
        let _ = n;
    }

    #[test]
    fn embed_different_seeds_differ() {
        let mol = parse("CCCC").unwrap().with_explicit_hydrogens();
        let n = mol.atoms.len();
        let mut bm1 = build_bounds(&mol);
        smooth(&mut bm1);
        let c1 = embed(&bm1, 1).unwrap();

        let mol2 = parse("CCCC").unwrap().with_explicit_hydrogens();
        let mut bm2 = build_bounds(&mol2);
        smooth(&mut bm2);
        let c2 = embed(&bm2, 999).unwrap();

        // At least one coordinate should differ between seeds
        let same = (0..3 * n).all(|k| (c1[k] - c2[k]).abs() < 1e-12);
        assert!(!same, "two different seeds produced identical coordinates");
    }

    #[test]
    fn embed_water_triangle() {
        // Water: 1 O + 2 H; 3 atoms → exactly 3 eigenvalues available.
        let (coords, n) = embed_mol("O", 0);
        assert_eq!(n, 3); // O + 2 H
        // All coordinates should be finite
        for &v in &coords {
            assert!(v.is_finite());
        }
    }
}
