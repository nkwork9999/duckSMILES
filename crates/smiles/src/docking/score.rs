use crate::docking::atomtype::VinaType;

// ── Vina VDW radii (Å) ───────────────────────────────────────────────────────
// From Vina source: atom_constants.h / model.cpp
pub fn vdw(t: VinaType) -> f64 {
    match t {
        VinaType::C | VinaType::A => 2.00,
        VinaType::N | VinaType::NA | VinaType::NS => 1.75,
        VinaType::O | VinaType::OA => 1.60,
        VinaType::S | VinaType::SA => 2.00,
        VinaType::P => 2.10,
        VinaType::F => 1.54,
        VinaType::CL => 1.80,
        VinaType::BR => 2.00,
        VinaType::I => 2.20,
        VinaType::H | VinaType::HD => 1.00,
        VinaType::MG | VinaType::FE | VinaType::MN | VinaType::CU => 0.65,
        VinaType::ZN => 0.74,
        VinaType::CA => 0.99,
        VinaType::Unknown => 1.70,
    }
}

// ── Interaction classes ───────────────────────────────────────────────────────

pub fn is_hydrophobic(t: VinaType) -> bool {
    matches!(t, VinaType::C | VinaType::A | VinaType::S | VinaType::SA
              | VinaType::F | VinaType::CL | VinaType::BR | VinaType::I)
}

pub fn is_hbond_acceptor(t: VinaType) -> bool {
    matches!(t, VinaType::NA | VinaType::OA | VinaType::SA | VinaType::NS)
}

pub fn is_hbond_donor(t: VinaType) -> bool {
    t == VinaType::HD
}

// ── Scoring terms (d = surface distance = r − r_vdw_i − r_vdw_j) ─────────────

#[inline]
pub fn gauss1(d: f64) -> f64 {
    (-(d / 0.5).powi(2)).exp()
}

#[inline]
pub fn gauss2(d: f64) -> f64 {
    (-((d - 3.0) / 2.0).powi(2)).exp()
}

#[inline]
pub fn repulsion(d: f64) -> f64 {
    if d < 0.0 { d * d } else { 0.0 }
}

#[inline]
pub fn hydrophobic(d: f64) -> f64 {
    if d < 0.5 { 1.0 } else if d > 1.5 { 0.0 } else { 1.5 - d }
}

#[inline]
pub fn hbond(d: f64) -> f64 {
    if d < -0.7 { 1.0 } else if d > 0.0 { 0.0 } else { -d / 0.7 }
}

// ── Weights ────────────────────────────────────────────────────────────────────

pub const W_GAUSS1: f64 = -0.035579;
pub const W_GAUSS2: f64 = -0.005156;
pub const W_REPULSION: f64 = 0.840245;
pub const W_HYDROPHOBIC: f64 = -0.035069;
pub const W_HBOND: f64 = -0.587439;
pub const W_ROT: f64 = 0.058459; // per rotatable bond

// ── Pairwise interaction energy ───────────────────────────────────────────────
// r: actual atom-atom distance (Å)
// Cutoff distance in terms of surface distance: 8 Å (Vina default).

const D_CUTOFF: f64 = 8.0;

pub fn pairwise(ti: VinaType, tj: VinaType, r: f64) -> f64 {
    let d = r - vdw(ti) - vdw(tj);
    if d > D_CUTOFF {
        return 0.0;
    }
    let mut e = W_GAUSS1 * gauss1(d) + W_GAUSS2 * gauss2(d) + W_REPULSION * repulsion(d);

    // Hydrophobic term — only if both atoms are hydrophobic
    if is_hydrophobic(ti) && is_hydrophobic(tj) {
        e += W_HYDROPHOBIC * hydrophobic(d);
    }

    // H-bond term — donor + acceptor pair (either order)
    let hb = (is_hbond_donor(ti) && is_hbond_acceptor(tj))
          || (is_hbond_donor(tj) && is_hbond_acceptor(ti));
    if hb {
        e += W_HBOND * hbond(d);
    }

    e
}

// ── Total score ────────────────────────────────────────────────────────────────
// Inter-molecular score between ligand atoms (with coords + types) and
// protein atoms (with coords + types).  Intra-molecular is excluded here;
// the search engine handles that separately.
//
// n_rot: number of rotatable bonds in the ligand (for the flexibility penalty)

#[allow(dead_code)]
pub fn total_score(
    lig_coords: &[[f64; 3]],
    lig_types: &[VinaType],
    prot_coords: &[[f64; 3]],
    prot_types: &[VinaType],
    n_rot: u32,
) -> f64 {
    let mut e = 0.0;
    for (lc, lt) in lig_coords.iter().zip(lig_types.iter()) {
        for (pc, pt) in prot_coords.iter().zip(prot_types.iter()) {
            let dx = lc[0] - pc[0];
            let dy = lc[1] - pc[1];
            let dz = lc[2] - pc[2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            e += pairwise(*lt, *pt, r);
        }
    }
    // Flexibility penalty
    e / (1.0 + W_ROT * n_rot as f64)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gauss1_at_zero() {
        assert!((gauss1(0.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn gauss1_at_half_sigma() {
        // d = 0.5: gauss1 = e^(-1) ≈ 0.3679
        assert!((gauss1(0.5) - (-1.0_f64).exp()).abs() < 1e-9);
    }

    #[test]
    fn gauss2_at_three() {
        assert!((gauss2(3.0) - 1.0).abs() < 1e-9);
    }

    #[test]
    fn gauss2_at_one() {
        // d = 1: gauss2 = e^(-((1-3)/2)^2) = e^(-1)
        assert!((gauss2(1.0) - (-1.0_f64).exp()).abs() < 1e-9);
    }

    #[test]
    fn repulsion_negative() {
        assert!((repulsion(-1.0) - 1.0).abs() < 1e-9);
        assert!((repulsion(-0.5) - 0.25).abs() < 1e-9);
    }

    #[test]
    fn repulsion_positive() {
        assert_eq!(repulsion(0.1), 0.0);
        assert_eq!(repulsion(0.0), 0.0);
    }

    #[test]
    fn hydrophobic_boundaries() {
        assert!((hydrophobic(0.0) - 1.0).abs() < 1e-9);
        assert!((hydrophobic(0.5) - 1.0).abs() < 1e-9);
        assert!((hydrophobic(1.0) - 0.5).abs() < 1e-9);
        assert!((hydrophobic(1.5) - 0.0).abs() < 1e-9);
        assert_eq!(hydrophobic(2.0), 0.0);
    }

    #[test]
    fn hbond_boundaries() {
        assert!((hbond(-1.0) - 1.0).abs() < 1e-9);
        assert!((hbond(-0.7) - 1.0).abs() < 1e-9);
        assert!((hbond(-0.35) - 0.5).abs() < 1e-6);
        assert!((hbond(0.0) - 0.0).abs() < 1e-9);
        assert_eq!(hbond(0.1), 0.0);
    }

    #[test]
    fn is_hydrophobic_classification() {
        assert!(is_hydrophobic(VinaType::C));
        assert!(is_hydrophobic(VinaType::A));
        assert!(is_hydrophobic(VinaType::CL));
        assert!(!is_hydrophobic(VinaType::N));
        assert!(!is_hydrophobic(VinaType::OA));
        assert!(!is_hydrophobic(VinaType::HD));
    }

    #[test]
    fn hbond_pair_detection() {
        // HD + OA should trigger H-bond term
        let r = vdw(VinaType::HD) + vdw(VinaType::OA) + (-0.3); // surface d = -0.3
        let e_hd_oa = pairwise(VinaType::HD, VinaType::OA, r);
        // OA + OA does NOT trigger H-bond
        let e_oa_oa = pairwise(VinaType::OA, VinaType::OA, r);
        // HD-OA should have more negative energy (better) due to H-bond
        assert!(e_hd_oa < e_oa_oa, "HD-OA should be more attractive: {e_hd_oa} vs {e_oa_oa}");
    }

    #[test]
    fn c_c_hydrophobic_pair() {
        // C-C at short surface distance should give more negative score
        // than N-N at same distance (N-N has no hydrophobic term)
        let r_cc = vdw(VinaType::C) + vdw(VinaType::C) + 0.8; // d = 0.8 (in hydrophobic range)
        let r_nn = vdw(VinaType::N) + vdw(VinaType::N) + 0.8;
        let e_cc = pairwise(VinaType::C, VinaType::C, r_cc);
        let e_nn = pairwise(VinaType::N, VinaType::N, r_nn);
        assert!(e_cc < e_nn, "C-C should be more attractive due to hydrophobic: {e_cc} vs {e_nn}");
    }

    #[test]
    fn pairwise_beyond_cutoff() {
        // At d >> cutoff, score should be 0
        let r = 100.0;
        assert_eq!(pairwise(VinaType::C, VinaType::C, r), 0.0);
    }

    #[test]
    fn total_score_decreases_with_rotation() {
        // At r=4.0 Å, C-C interaction is attractive (d=0 Å, gauss+hydrophobic dominate).
        // More rotatable bonds divide by (1 + W_ROT*n_rot), making the score less negative.
        let lc = vec![[0.0, 0.0, 0.0]];
        let lt = vec![VinaType::C];
        let pc = vec![[4.0, 0.0, 0.0]]; // d = 4 - 2 - 2 = 0 Å (contact)
        let pt = vec![VinaType::C];
        let s0 = total_score(&lc, &lt, &pc, &pt, 0);
        let s5 = total_score(&lc, &lt, &pc, &pt, 5);
        // s0 should be negative (attractive) and s0 < s5 (s5 less negative)
        assert!(s0 < 0.0, "C-C at contact should be attractive: s0={s0}");
        assert!(s0 < s5, "more rot bonds should make score less negative: s0={s0} s5={s5}");
    }
}
