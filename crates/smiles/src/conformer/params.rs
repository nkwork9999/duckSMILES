use crate::parser::BondOrder;

// ── Hybridization ─────────────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
}

// ── VDW radii (Å) ─────────────────────────────────────────────────────────────
// Bondi 1964 values used by RDKit / AutoDock

pub fn vdw_radius(sym: &str) -> f64 {
    match sym {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "S" => 1.80,
        "P" => 1.80,
        "F" => 1.47,
        "Cl" => 1.75,
        "Br" => 1.85,
        "I" => 1.98,
        "Si" => 2.10,
        "Se" => 1.90,
        _ => 1.70, // fallback to carbon
    }
}

// ── Covalent radii (Å) ────────────────────────────────────────────────────────
// Used to estimate VDW lower bounds for non-bonded pairs.

pub fn covalent_radius(sym: &str) -> f64 {
    match sym {
        "H" => 0.31,
        "C" => 0.76,
        "N" => 0.71,
        "O" => 0.66,
        "S" => 1.05,
        "P" => 1.07,
        "F" => 0.57,
        "Cl" => 1.02,
        "Br" => 1.20,
        "I" => 1.39,
        "Si" => 1.11,
        "Se" => 1.20,
        _ => 0.76,
    }
}

// ── Bond lengths (Å) ──────────────────────────────────────────────────────────
// Returns ideal bond length for a pair of elements at the given bond order.
// Ordered so (a, b) and (b, a) both match.

pub fn bond_length(sym_a: &str, sym_b: &str, order: BondOrder) -> f64 {
    // Canonicalise order of elements so the match arms need only one direction.
    let (a, b) = canonical_pair(sym_a, sym_b);

    match (a, b, order) {
        // ── C–* ──────────────────────────────────────────────────────────
        ("C", "C", BondOrder::Single) => 1.540,
        ("C", "C", BondOrder::Double) => 1.340,
        ("C", "C", BondOrder::Triple) => 1.200,
        ("C", "C", BondOrder::Aromatic) => 1.400,

        ("C", "H", BondOrder::Single | BondOrder::Aromatic) => 1.090,

        ("C", "N", BondOrder::Single) => 1.470,
        ("C", "N", BondOrder::Double) => 1.270,
        ("C", "N", BondOrder::Triple) => 1.150,
        ("C", "N", BondOrder::Aromatic) => 1.340,

        ("C", "O", BondOrder::Single) => 1.430,
        ("C", "O", BondOrder::Double) => 1.200,
        ("C", "O", BondOrder::Aromatic) => 1.360,

        ("C", "S", BondOrder::Single) => 1.820,
        ("C", "S", BondOrder::Double) => 1.600,
        ("C", "S", BondOrder::Aromatic) => 1.760,

        ("C", "P", BondOrder::Single) => 1.840,
        ("C", "P", BondOrder::Double) => 1.660,

        ("C", "F", BondOrder::Single | BondOrder::Aromatic) => 1.350,
        ("C", "Cl", BondOrder::Single | BondOrder::Aromatic) => 1.770,
        ("C", "Br", BondOrder::Single | BondOrder::Aromatic) => 1.940,
        ("C", "I", BondOrder::Single | BondOrder::Aromatic) => 2.140,
        ("C", "Si", BondOrder::Single) => 1.870,
        ("C", "Se", BondOrder::Single) => 1.970,

        // ── N–* ──────────────────────────────────────────────────────────
        ("H", "N", BondOrder::Single | BondOrder::Aromatic) => 1.010,
        ("N", "N", BondOrder::Single) => 1.450,
        ("N", "N", BondOrder::Double) => 1.250,
        ("N", "N", BondOrder::Aromatic) => 1.350,
        ("N", "O", BondOrder::Single) => 1.400,
        ("N", "O", BondOrder::Double) => 1.210,

        // ── O–* ──────────────────────────────────────────────────────────
        ("H", "O", BondOrder::Single | BondOrder::Aromatic) => 0.960,
        ("O", "O", BondOrder::Single) => 1.480,
        ("O", "O", BondOrder::Double) => 1.210,
        ("O", "P", BondOrder::Single) => 1.610,
        ("O", "P", BondOrder::Double) => 1.480,
        ("O", "S", BondOrder::Single) => 1.650,
        ("O", "S", BondOrder::Double) => 1.430,

        // ── S–* ──────────────────────────────────────────────────────────
        ("H", "S", BondOrder::Single) => 1.340,
        ("S", "S", BondOrder::Single) => 2.050,
        ("S", "S", BondOrder::Double) => 1.890,

        // ── P–* ──────────────────────────────────────────────────────────
        ("H", "P", BondOrder::Single) => 1.420,

        // ── fallback: sum of covalent radii ──────────────────────────────
        _ => covalent_radius(sym_a) + covalent_radius(sym_b),
    }
}

// ── Bond angles (°) ───────────────────────────────────────────────────────────
// Ideal valence angle at `center` with the given hybridization.

pub fn ideal_angle(center: &str, hyb: Hybridization) -> f64 {
    match (center, hyb) {
        // Special cases override hybridization rule
        ("O", Hybridization::SP3) => 104.5, // water / ether
        ("N", Hybridization::SP3) => 107.0, // ammonia-like
        ("S", Hybridization::SP3) => 100.0,
        // General
        (_, Hybridization::SP) => 180.0,
        (_, Hybridization::SP2) => 120.0,
        (_, Hybridization::SP3) => 109.5,
    }
}

// ── helpers ───────────────────────────────────────────────────────────────────

fn canonical_pair<'a>(a: &'a str, b: &'a str) -> (&'a str, &'a str) {
    if a <= b { (a, b) } else { (b, a) }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bond_length_cc_single() {
        assert!((bond_length("C", "C", BondOrder::Single) - 1.540).abs() < 1e-9);
    }

    #[test]
    fn bond_length_cc_double() {
        assert!((bond_length("C", "C", BondOrder::Double) - 1.340).abs() < 1e-9);
    }

    #[test]
    fn bond_length_cc_triple() {
        assert!((bond_length("C", "C", BondOrder::Triple) - 1.200).abs() < 1e-9);
    }

    #[test]
    fn bond_length_cc_aromatic() {
        assert!((bond_length("C", "C", BondOrder::Aromatic) - 1.400).abs() < 1e-9);
    }

    #[test]
    fn bond_length_cn_single() {
        assert!((bond_length("C", "N", BondOrder::Single) - 1.470).abs() < 1e-9);
    }

    #[test]
    fn bond_length_co_double() {
        assert!((bond_length("C", "O", BondOrder::Double) - 1.200).abs() < 1e-9);
    }

    #[test]
    fn bond_length_ch() {
        assert!((bond_length("C", "H", BondOrder::Single) - 1.090).abs() < 1e-9);
    }

    #[test]
    fn bond_length_oh() {
        assert!((bond_length("O", "H", BondOrder::Single) - 0.960).abs() < 1e-9);
    }

    #[test]
    fn bond_length_symmetric() {
        // order of elements must not matter
        assert_eq!(
            bond_length("N", "C", BondOrder::Double),
            bond_length("C", "N", BondOrder::Double)
        );
        assert_eq!(
            bond_length("O", "C", BondOrder::Single),
            bond_length("C", "O", BondOrder::Single)
        );
    }

    #[test]
    fn bond_length_fallback_uses_covalent_radii() {
        // Fe-C: not in the table, should return sum of cov radii > 0
        let d = bond_length("Fe", "C", BondOrder::Single);
        assert!(d > 1.0 && d < 3.0);
    }

    #[test]
    fn vdw_radius_known() {
        assert!((vdw_radius("C") - 1.70).abs() < 1e-9);
        assert!((vdw_radius("N") - 1.55).abs() < 1e-9);
        assert!((vdw_radius("O") - 1.52).abs() < 1e-9);
        assert!((vdw_radius("H") - 1.20).abs() < 1e-9);
        assert!((vdw_radius("Cl") - 1.75).abs() < 1e-9);
    }

    #[test]
    fn ideal_angle_sp3_carbon() {
        assert!((ideal_angle("C", Hybridization::SP3) - 109.5).abs() < 1e-9);
    }

    #[test]
    fn ideal_angle_sp2_carbon() {
        assert!((ideal_angle("C", Hybridization::SP2) - 120.0).abs() < 1e-9);
    }

    #[test]
    fn ideal_angle_sp_carbon() {
        assert!((ideal_angle("C", Hybridization::SP) - 180.0).abs() < 1e-9);
    }

    #[test]
    fn ideal_angle_water_oxygen() {
        assert!((ideal_angle("O", Hybridization::SP3) - 104.5).abs() < 1e-9);
    }

    #[test]
    fn ideal_angle_nitrogen_sp3() {
        assert!((ideal_angle("N", Hybridization::SP3) - 107.0).abs() < 1e-9);
    }
}
