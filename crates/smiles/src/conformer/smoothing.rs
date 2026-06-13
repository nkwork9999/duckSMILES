use crate::conformer::bounds::BoundsMatrix;

// ── Triangle Inequality Smoothing (Floyd-Warshall) ────────────────────────────
//
// For every triple (i, j, k) update:
//   up(i,j)  = min(up(i,j), up(i,k) + up(k,j))   — upper can only shrink
//   lo(i,j)  = max(lo(i,j), lo(i,k) − up(k,j))   — lower can only grow
//              = max(lo(i,j), lo(k,j) − up(i,k))
//
// Convergence is guaranteed in O(n³) passes of a single triple loop.
// We run the standard Floyd-Warshall order (k outer, i/j inner) which
// provides exact convergence for upper bounds; for lower bounds we run
// a second pass to propagate the tightened lowers.

pub fn smooth(bm: &mut BoundsMatrix) {
    let n = bm.n;

    // Two passes ensure both upper and lower bounds converge.
    for _ in 0..2 {
        for k in 0..n {
            for i in 0..n {
                if i == k {
                    continue;
                }
                for j in (i + 1)..n {
                    if j == k {
                        continue;
                    }

                    // Upper bound: triangle inequality
                    let new_up = bm.up(i, k) + bm.up(k, j);
                    bm.tighten_up(i, j, new_up);

                    // Lower bound: reverse triangle inequality
                    let new_lo_a = bm.lo(i, k) - bm.up(k, j);
                    let new_lo_b = bm.lo(k, j) - bm.up(i, k);
                    let new_lo = new_lo_a.max(new_lo_b).max(0.0);
                    bm.tighten_lo(i, j, new_lo);
                }
            }
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::parser::parse;

    fn check_triangle(bm: &BoundsMatrix) {
        let n = bm.n;
        for i in 0..n {
            for j in (i + 1)..n {
                assert!(
                    bm.lo(i, j) <= bm.up(i, j),
                    "lo > up at ({i},{j}): lo={} up={}",
                    bm.lo(i, j),
                    bm.up(i, j)
                );
                for k in 0..n {
                    if k == i || k == j {
                        continue;
                    }
                    // upper triangle inequality: up(i,j) ≤ up(i,k) + up(k,j)
                    let sum_up = bm.up(i, k) + bm.up(k, j);
                    assert!(
                        bm.up(i, j) <= sum_up + 1e-9,
                        "upper TI violated at ({i},{j},{k}): up={} > {}",
                        bm.up(i, j),
                        sum_up
                    );
                }
            }
        }
    }

    #[test]
    fn smooth_ethane_triangle_ok() {
        let mol = parse("CC").unwrap().with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        check_triangle(&bm);
    }

    #[test]
    fn smooth_propane_triangle_ok() {
        let mol = parse("CCC").unwrap().with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        check_triangle(&bm);
    }

    #[test]
    fn smooth_aspirin_triangle_ok() {
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O")
            .unwrap()
            .with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        check_triangle(&bm);
    }

    #[test]
    fn smooth_tightens_upper() {
        // In propane, after smoothing the C0-C2 upper bound should be tighter
        // than 100 Å (the loose initial value for non-bonded pairs).
        let mol = parse("CCC").unwrap().with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        let up_before = bm.up(0, 2);
        smooth(&mut bm);
        let up_after = bm.up(0, 2);
        assert!(
            up_after <= up_before,
            "upper should not grow after smoothing: {up_before} -> {up_after}"
        );
        // C0-C2 should be < 5 Å after smoothing (it's a 1-3 pair ≈ 2.5 Å)
        assert!(up_after < 5.0, "C0-C2 upper too loose after smoothing: {up_after}");
    }

    #[test]
    fn smooth_benzene_triangle_ok() {
        let mol = parse("c1ccccc1").unwrap().with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        check_triangle(&bm);
    }
}
