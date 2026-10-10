//! MinHash compression from Probst & Reymond (2018), Eq. 2,
//! DOI 10.1186/s13321-018-0321-8; also MAP4 Eq. 1 (2020).
//! The token source is explicit: use molecular shingles for MHFP/MAP4, or
//! unfolded Morgan invariants for MHECFP. Coefficients use our versioned SplitMix64
//! stream, NOT NumPy's author-code RNG; sketches are not bit-compatible with it.
//! Compare matching coordinates, never cosine/distance between integer magnitudes.
//! Inputs must already be well-distributed 32-bit hashes (as in the paper's SHA-1
//! shingles); consecutive integer IDs can bias affine-family Jaccard estimates.
use std::collections::BTreeSet;
const PRIME: u64 = (1u64 << 61) - 1;
const MAX: u64 = u32::MAX as u64;
pub const SPEC: &str = "paper-minhash-splitmix64-v1";
fn next(s: &mut u64) -> u64 {
    *s = s.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *s;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}
#[derive(Clone)]
pub struct MinHash {
    coefficients: Vec<(u32, u32)>,
}
impl MinHash {
    pub fn new(dimensions: usize, seed: u64) -> Result<Self, &'static str> {
        if !(1..=4096).contains(&dimensions) {
            return Err("dimensions outside 1..4096");
        }
        let mut state = seed;
        let mut seen_a = BTreeSet::new();
        let mut seen_b = BTreeSet::new();
        let mut coefficients = Vec::new();
        for _ in 0..dimensions {
            let mut a = (next(&mut state) % (MAX - 1) + 1) as u32;
            while !seen_a.insert(a) {
                a = (next(&mut state) % (MAX - 1) + 1) as u32;
            }
            let mut b = (next(&mut state) % MAX) as u32;
            while !seen_b.insert(b) {
                b = (next(&mut state) % MAX) as u32;
            }
            coefficients.push((a, b));
        }
        Ok(Self { coefficients })
    }
    pub fn coefficients(&self) -> &[(u32, u32)] {
        &self.coefficients
    }
    pub fn encode(&self, tokens: &[u32]) -> Vec<u32> {
        self.coefficients
            .iter()
            .map(|&(a, b)| {
                tokens
                    .iter()
                    .map(|&x| {
                        // Maximum product plus offset fits u64 for these bounded inputs.
                        (((u64::from(a) * u64::from(x) + u64::from(b)) % PRIME) % MAX) as u32
                    })
                    .min()
                    .unwrap_or(u32::MAX)
            })
            .collect()
    }
}
pub fn similarity(a: &[u32], b: &[u32]) -> Result<f64, &'static str> {
    if a.len() != b.len() || a.is_empty() {
        return Err("nonempty equal dimensions required");
    }
    Ok(a.iter().zip(b).filter(|(x, y)| x == y).count() as f64 / a.len() as f64)
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn equation_boundaries_and_set_invariance() {
        let m = MinHash::new(64, 42).unwrap();
        let t = [0, 1, u32::MAX - 1, u32::MAX];
        let fp = m.encode(&t);
        for (i, &(a, b)) in m.coefficients().iter().enumerate() {
            let expected = t
                .iter()
                .map(|&x| {
                    (((a as u128 * x as u128 + b as u128) % PRIME as u128) % MAX as u128) as u32
                })
                .min()
                .unwrap();
            assert_eq!(fp[i], expected);
        }
        assert_eq!(m.encode(&[1, 2, 3]), m.encode(&[3, 2, 1, 1]));
        assert_eq!(similarity(&fp, &fp), Ok(1.));
        assert!(similarity(&[], &[]).is_err());
        assert!(MinHash::new(0, 42).is_err());
    }
    #[test]
    fn approximates_jaccard_not_numeric_distance() {
        let m = MinHash::new(2048, 42).unwrap();
        let mut state = 12345;
        let tokens: Vec<u32> = (0..150).map(|_| next(&mut state) as u32).collect();
        let a = m.encode(&tokens[..100]);
        let b = m.encode(&tokens[50..]);
        assert!((similarity(&a, &b).unwrap() - 1. / 3.).abs() < 0.05);
        assert_eq!(similarity(&m.encode(&[]), &a), Ok(0.));
    }
}
