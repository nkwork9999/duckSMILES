//! RDKit 2026.03.6-compatible element/isotope masses, generated reference data.
fn element(symbol: &str) -> Option<(usize, f64, f64)> {
    crate::periodic_data::ELEMENTS
        .iter()
        .enumerate()
        .find(|(_, (s, _, _))| *s == symbol)
        .map(|(index, (_, average, exact))| (index + 1, *average, *exact))
}

pub fn atomic_weight(symbol: &str) -> Option<f64> {
    element(symbol).map(|(_, average, _)| average)
}

pub fn monoisotopic_mass(symbol: &str) -> Option<f64> {
    element(symbol)
        .map(|(_, _, exact)| exact)
        .filter(|m| *m > 0.0)
}

pub fn isotope_mass(symbol: &str, isotope: u16) -> Option<f64> {
    let (z, _, _) = element(symbol)?;
    let values = crate::periodic_data::ISOTOPES;
    let index = values
        .binary_search_by_key(&(z as u8, isotope), |(z, i, _)| (*z, *i))
        .ok()?;
    Some(values[index].2)
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
