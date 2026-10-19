//! RDKit registered fragment counts from pinned SMARTS; no embedded Python.
use crate::fragment_patterns::PATTERNS;
use crate::parser::{Molecule, RingInfo};
use crate::smarts::{parse_smarts, Pattern, PreparedMatcher};
use std::sync::OnceLock;

fn patterns() -> &'static [Option<Pattern>] {
    static COMPILED: OnceLock<Vec<Option<Pattern>>> = OnceLock::new();
    COMPILED.get_or_init(|| PATTERNS.iter().map(|(_, s)| parse_smarts(s)).collect())
}

pub(crate) fn calculate(mol: &Molecule, ring: &RingInfo) -> Vec<Option<f64>> {
    let matcher = PreparedMatcher::with_ring_info(mol, ring);
    patterns()
        .iter()
        .map(|pattern| {
            // A missing pattern is an implementation gap, never zero occurrences.
            pattern
                .as_ref()
                .map(|p| matcher.count_unique_capped(p, 1000) as f64)
        })
        .collect()
}

pub(crate) fn index(id: &str) -> Option<usize> {
    PATTERNS.binary_search_by_key(&id, |p| p.0).ok()
}
