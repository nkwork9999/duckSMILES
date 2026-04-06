/// Average atomic weight (IUPAC 2021)
pub fn atomic_weight(symbol: &str) -> Option<f64> {
    match symbol {
        "H" => Some(1.008),
        "He" => Some(4.0026),
        "Li" => Some(6.941),
        "Be" => Some(9.0122),
        "B" => Some(10.81),
        "C" => Some(12.011),
        "N" => Some(14.007),
        "O" => Some(15.999),
        "F" => Some(18.998),
        "Ne" => Some(20.180),
        "Na" => Some(22.990),
        "Mg" => Some(24.305),
        "Al" => Some(26.982),
        "Si" => Some(28.085),
        "P" => Some(30.974),
        "S" => Some(32.06),
        "Cl" => Some(35.45),
        "Ar" => Some(39.948),
        "K" => Some(39.098),
        "Ca" => Some(40.078),
        "Ti" => Some(47.867),
        "V" => Some(50.942),
        "Cr" => Some(51.996),
        "Mn" => Some(54.938),
        "Fe" => Some(55.845),
        "Co" => Some(58.933),
        "Ni" => Some(58.693),
        "Cu" => Some(63.546),
        "Zn" => Some(65.38),
        "Ga" => Some(69.723),
        "Ge" => Some(72.630),
        "As" => Some(74.922),
        "Se" => Some(78.971),
        "Br" => Some(79.904),
        "Kr" => Some(83.798),
        "Rb" => Some(85.468),
        "Sr" => Some(87.62),
        "Zr" => Some(91.224),
        "Mo" => Some(95.95),
        "Ru" => Some(101.07),
        "Rh" => Some(102.91),
        "Pd" => Some(106.42),
        "Ag" => Some(107.87),
        "Cd" => Some(112.41),
        "In" => Some(114.82),
        "Sn" => Some(118.71),
        "Sb" => Some(121.76),
        "Te" => Some(127.60),
        "I" => Some(126.90),
        "Xe" => Some(131.29),
        "Cs" => Some(132.91),
        "Ba" => Some(137.33),
        "Pt" => Some(195.08),
        "Au" => Some(196.97),
        "Hg" => Some(200.59),
        "Pb" => Some(207.2),
        "Bi" => Some(208.98),
        _ => None,
    }
}

/// Monoisotopic mass (NIST)
pub fn monoisotopic_mass(symbol: &str) -> Option<f64> {
    match symbol {
        "H" => Some(1.00782503),
        "B" => Some(11.00930536),
        "C" => Some(12.0),
        "N" => Some(14.00307401),
        "O" => Some(15.99491462),
        "F" => Some(18.99840322),
        "Na" => Some(22.98976928),
        "Mg" => Some(23.98504170),
        "Al" => Some(26.98153863),
        "Si" => Some(27.97692653),
        "P" => Some(30.97376163),
        "S" => Some(31.97207100),
        "Cl" => Some(34.96885268),
        "K" => Some(38.96370668),
        "Ca" => Some(39.96259098),
        "Fe" => Some(55.93493633),
        "Cu" => Some(62.92959772),
        "Zn" => Some(63.92914201),
        "Br" => Some(78.91833710),
        "Se" => Some(79.91652130),
        "I" => Some(126.90447190),
        "Ag" => Some(106.90509300),
        "Au" => Some(196.96656870),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_carbon_weight() {
        assert!((atomic_weight("C").unwrap() - 12.011).abs() < 0.001);
    }

    #[test]
    fn test_carbon_mono() {
        assert!((monoisotopic_mass("C").unwrap() - 12.0).abs() < 0.0001);
    }

    #[test]
    fn test_unknown() {
        assert!(atomic_weight("Xx").is_none());
        assert!(monoisotopic_mass("Xx").is_none());
    }
}
