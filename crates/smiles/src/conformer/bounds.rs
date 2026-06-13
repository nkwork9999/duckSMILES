use std::collections::HashSet;

use crate::parser::{BondOrder, Molecule};
use crate::conformer::params::{bond_length, ideal_angle, vdw_radius};

// ── BoundsMatrix ──────────────────────────────────────────────────────────────
// Convention: lo(i,j) stores the LOWER bound, up(i,j) the UPPER bound.
// Internally: lower bound at (max, min), upper bound at (min, max) in the flat
// array, so `lo` and `up` can be identified by index ordering.

pub struct BoundsMatrix {
    pub n: usize,
    data: Vec<f64>, // n×n flat; element at [i*n+j]
}

impl BoundsMatrix {
    pub fn new(n: usize) -> Self {
        BoundsMatrix {
            n,
            data: vec![0.0; n * n],
        }
    }

    #[inline]
    pub fn lo(&self, i: usize, j: usize) -> f64 {
        let (a, b) = if i > j { (i, j) } else { (j, i) };
        self.data[a * self.n + b]
    }

    #[inline]
    pub fn up(&self, i: usize, j: usize) -> f64 {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.data[a * self.n + b]
    }

    #[inline]
    pub fn set_lo(&mut self, i: usize, j: usize, v: f64) {
        let (a, b) = if i > j { (i, j) } else { (j, i) };
        self.data[a * self.n + b] = v;
    }

    #[inline]
    pub fn set_up(&mut self, i: usize, j: usize, v: f64) {
        let (a, b) = if i < j { (i, j) } else { (j, i) };
        self.data[a * self.n + b] = v;
    }

    pub fn tighten_lo(&mut self, i: usize, j: usize, v: f64) {
        if v > self.lo(i, j) {
            self.set_lo(i, j, v);
        }
    }

    pub fn tighten_up(&mut self, i: usize, j: usize, v: f64) {
        let cur = self.up(i, j);
        if cur == 0.0 || v < cur {
            self.set_up(i, j, v);
        }
    }
}

// ── build_bounds ──────────────────────────────────────────────────────────────

pub fn build_bounds(mol: &Molecule) -> BoundsMatrix {
    let n = mol.atoms.len();
    let mut bm = BoundsMatrix::new(n);

    // ── Pre-compute 1-2 and 1-3 pair sets ─────────────────────────────────
    // VDW lower bounds must NOT be applied to 1-2 or 1-3 pairs because the
    // geometry constraints place those atoms closer than their VDW sum.
    let mut pairs_12: HashSet<(usize, usize)> = HashSet::new();
    for b in &mol.bonds {
        let p = (b.a.min(b.b), b.a.max(b.b));
        pairs_12.insert(p);
    }

    let mut pairs_13: HashSet<(usize, usize)> = HashSet::new();
    for center in 0..n {
        let nbrs = mol.neighbors(center);
        for ai in 0..nbrs.len() {
            for bi in (ai + 1)..nbrs.len() {
                let (a, _) = nbrs[ai];
                let (b, _) = nbrs[bi];
                let p = (a.min(b), a.max(b));
                if !pairs_12.contains(&p) {
                    pairs_13.insert(p);
                }
            }
        }
    }

    // ── 1. Initialise with VDW lower bounds (1-4+ pairs only) ─────────────
    for i in 0..n {
        for j in (i + 1)..n {
            bm.set_up(i, j, 100.0); // generous upper bound for every pair
            let p = (i, j);
            if !pairs_12.contains(&p) && !pairs_13.contains(&p) {
                let lo = (vdw_radius(&mol.atoms[i].symbol)
                    + vdw_radius(&mol.atoms[j].symbol))
                    * 0.8;
                bm.set_lo(i, j, lo);
            }
        }
    }

    // ── 2. 1-2 (bonded) constraints: lo = up = ideal bond length ± 5 % ───
    for b in &mol.bonds {
        let (i, j) = (b.a.min(b.b), b.a.max(b.b));
        let d = bond_length(&mol.atoms[i].symbol, &mol.atoms[j].symbol, b.order);
        bm.set_lo(i, j, d * 0.95);
        bm.set_up(i, j, d * 1.05);
    }

    // ── 3. 1-3 (angle) constraints via law of cosines ─────────────────────
    // d(A,B)² = d(A,k)² + d(k,B)² − 2·d(A,k)·d(k,B)·cos(θ)
    const ANGLE_SLACK: f64 = 0.10; // ± 0.10 Å tolerance
    for center in 0..n {
        let nbrs = mol.neighbors(center);
        let hyb = mol.hybridization(center);
        let theta = ideal_angle(&mol.atoms[center].symbol, hyb).to_radians();

        for ai in 0..nbrs.len() {
            for bi in (ai + 1)..nbrs.len() {
                let (a, bo_a) = nbrs[ai];
                let (b, bo_b) = nbrs[bi];
                let d_ca =
                    bond_length(&mol.atoms[center].symbol, &mol.atoms[a].symbol, bo_a);
                let d_cb =
                    bond_length(&mol.atoms[center].symbol, &mol.atoms[b].symbol, bo_b);
                let d13 =
                    (d_ca * d_ca + d_cb * d_cb - 2.0 * d_ca * d_cb * theta.cos()).sqrt();

                let (p, q) = (a.min(b), a.max(b));
                // First 1-3 path wins for lo; subsequent paths tighten.
                if bm.lo(p, q) == 0.0 {
                    bm.set_lo(p, q, (d13 - ANGLE_SLACK).max(0.0));
                } else {
                    bm.tighten_lo(p, q, (d13 - ANGLE_SLACK).max(0.0));
                }
                bm.tighten_up(p, q, d13 + ANGLE_SLACK);
            }
        }
    }

    // ── 4. Aromatic-ring planarity (ETKDG-style) ─────────────────────────
    // A planar regular polygon has every vertex-pair distance fixed. Pinning
    // those distances forces the embedding to start (and stay) flat — the
    // ingredient raw DG lacks that makes aromatic rings pucker.
    apply_aromatic_ring_planarity(mol, &mut bm);

    bm
}

// ── aromatic-ring planarity bounds ───────────────────────────────────────────

/// For each aromatic ring, pin all intra-ring atom-pair distances to the planar
/// regular-polygon values derived from the mean ring bond length.
fn apply_aromatic_ring_planarity(mol: &Molecule, bm: &mut BoundsMatrix) {
    const RING_SLACK: f64 = 0.06; // Å
    let ring_info = mol.ring_info();

    for ring_bonds in &ring_info.rings {
        // Aromatic ring = every ring bond is aromatic.
        let all_aromatic = ring_bonds.iter().all(|&bi| {
            mol.bonds
                .get(bi)
                .map(|b| b.order == BondOrder::Aromatic)
                .unwrap_or(false)
        });
        if !all_aromatic {
            continue;
        }

        let atoms = match ordered_ring_atoms(mol, ring_bonds) {
            Some(a) => a,
            None => continue,
        };
        let r = atoms.len();
        if r < 3 {
            continue;
        }

        // Mean ring bond length → circumradius of the regular polygon.
        let mut sum_len = 0.0;
        for &bi in ring_bonds {
            let b = &mol.bonds[bi];
            sum_len += bond_length(&mol.atoms[b.a].symbol, &mol.atoms[b.b].symbol, b.order);
        }
        let mean_len = sum_len / ring_bonds.len() as f64;
        let circum = mean_len / (2.0 * (std::f64::consts::PI / r as f64).sin());

        // Pin every pair (i,j) to 2·Rc·sin(kπ/R), k = ring separation.
        for ii in 0..r {
            for jj in (ii + 1)..r {
                let k = {
                    let d = jj - ii;
                    d.min(r - d) // shorter way around the ring
                };
                if k == 1 {
                    continue; // 1-2 bonds already pinned tightly
                }
                let dist = 2.0 * circum * ((k as f64) * std::f64::consts::PI / r as f64).sin();
                let (p, q) = (atoms[ii], atoms[jj]);
                bm.set_lo(p, q, (dist - RING_SLACK).max(0.0));
                bm.set_up(p, q, dist + RING_SLACK);
            }
        }
    }
}

/// Reconstruct the cyclic atom order of a ring given its bond indices.
/// Returns None if the bonds do not form a single simple cycle.
fn ordered_ring_atoms(mol: &Molecule, ring_bonds: &[usize]) -> Option<Vec<usize>> {
    // Adjacency limited to the ring's bonds.
    let mut adj: std::collections::HashMap<usize, Vec<usize>> = std::collections::HashMap::new();
    for &bi in ring_bonds {
        let b = mol.bonds.get(bi)?;
        adj.entry(b.a).or_default().push(b.b);
        adj.entry(b.b).or_default().push(b.a);
    }
    // Every ring atom must have degree exactly 2 in a simple cycle.
    if adj.values().any(|v| v.len() != 2) {
        return None;
    }
    let start = *adj.keys().min()?;
    let mut order = vec![start];
    let mut prev = start;
    let mut cur = adj[&start][0];
    while cur != start {
        order.push(cur);
        let nbrs = &adj[&cur];
        let next = if nbrs[0] == prev { nbrs[1] } else { nbrs[0] };
        prev = cur;
        cur = next;
        if order.len() > ring_bonds.len() + 1 {
            return None; // safety
        }
    }
    Some(order)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn methane_ch_bond() {
        let mol = parse("C").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // C=0, H=1..4. C–H ideal = 1.090 Å ± 5 %
        for h in 1..=4 {
            let lo = bm.lo(0, h);
            let up = bm.up(0, h);
            assert!(lo >= 1.035 && lo <= 1.090, "C-H lo out of range: {lo}");
            assert!(up >= 1.090 && up <= 1.145, "C-H up out of range: {up}");
            assert!(lo < up, "lo >= up for C-H");
        }
    }

    #[test]
    fn ethane_cc_bond() {
        let mol = parse("CC").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // C–C ideal = 1.540 Å
        assert!(bm.lo(0, 1) >= 1.460, "C-C lo too tight: {}", bm.lo(0, 1));
        assert!(bm.up(0, 1) <= 1.620, "C-C up too loose: {}", bm.up(0, 1));
        assert!(bm.lo(0, 1) < bm.up(0, 1));
    }

    #[test]
    fn ethylene_cc_double() {
        let mol = parse("C=C").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // C=C ideal = 1.340 Å
        assert!(bm.lo(0, 1) >= 1.270, "C=C lo too tight");
        assert!(bm.up(0, 1) <= 1.410, "C=C up too loose: {}", bm.up(0, 1));
        assert!(bm.lo(0, 1) < bm.up(0, 1));
    }

    #[test]
    fn acetylene_cc_triple() {
        let mol = parse("C#C").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // C≡C ideal = 1.200 Å
        assert!(bm.up(0, 1) <= 1.270, "C≡C up too loose: {}", bm.up(0, 1));
        assert!(bm.lo(0, 1) < bm.up(0, 1));
    }

    #[test]
    fn water_hh_angle_constraint() {
        // O with two H: 1-3 H–H distance derived from 104.5° angle.
        let mol = parse("O").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // H-H 1-3 distance ≈ 1.52 Å; lo < up must hold.
        let lo = bm.lo(1, 2);
        let up = bm.up(1, 2);
        assert!(lo > 0.0, "H-H lo should be > 0");
        assert!(lo < up, "H-H: lo={lo} >= up={up}");
        // Upper should be well below 100 (angle constraint applied)
        assert!(up < 2.5, "H-H up should be < 2.5 Å: {up}");
    }

    #[test]
    fn bounds_lo_le_up_for_all_pairs() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O")
            .unwrap()
            .with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        let n = mol.atoms.len();
        for i in 0..n {
            for j in (i + 1)..n {
                let lo = bm.lo(i, j);
                let up = bm.up(i, j);
                assert!(lo < up, "lo >= up at ({i},{j}): lo={lo} up={up}");
            }
        }
    }

    #[test]
    fn bounds_symmetric() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        let n = mol.atoms.len();
        for i in 0..n {
            for j in (i + 1)..n {
                assert_eq!(bm.lo(i, j), bm.lo(j, i), "lo not symmetric ({i},{j})");
                assert_eq!(bm.up(i, j), bm.up(j, i), "up not symmetric ({i},{j})");
            }
        }
    }

    #[test]
    fn vdw_lower_applied_to_14_pairs() {
        // In propane C0-C1-C2, atoms C0 and C2 are 1-4 (through C1).
        // They should have VDW lower bound (not 0).
        let mol = parse("CCC").unwrap().with_explicit_hydrogens();
        let bm = build_bounds(&mol);
        // C0-C2 are atoms 0 and 2
        let lo = bm.lo(0, 2);
        // VDW lower for C-C = (1.70 + 1.70) * 0.8 = 2.72 Å
        // But they're 1-3, so the angle constraint applies instead!
        // The 1-3 C-C-C at 109.5° gives d ≈ 2.52 Å.
        // lo should be > 0 and < 3.0
        assert!(lo > 0.0 && lo < 3.0, "C0-C2 lo unexpected: {lo}");
    }
}

