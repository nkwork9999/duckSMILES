//! Bit-fingerprint similarity metrics beyond Tanimoto.
//!
//! Faithful ports of the symmetric similarity family in
//! `RDKit/Code/DataStructs/BitOps.cpp`. Each metric is expressed in terms of
//! three popcounts over the raw BLOB bytes:
//!
//! ```text
//! common = popcount(a AND b)   // NumOnBitsInCommon
//! na     = popcount(a)         // bv1.getNumOnBits()
//! nb     = popcount(b)         // bv2.getNumOnBits()
//! ```
//!
//! `tanimoto_bit` lives in its own module (it predates this one); the metrics
//! here reuse the same chunk-at-a-time popcount strategy and the same
//! degenerate-input conventions RDKit uses (return 0.0 rather than dividing by
//! zero, so `ORDER BY sim DESC` keeps working on all-zero fingerprints).
//!
//! The caller must ensure `a.len() == b.len()`; the FFI wrappers enforce this
//! and surface a clear error from the DuckDB extension on mismatch.

/// On-bit counts for a pair of equal-length bit vectors, computed in a single
/// pass: `(common, na, nb)` = `(popcount(a&b), popcount(a), popcount(b))`.
fn bit_counts(a: &[u8], b: &[u8]) -> (u64, u64, u64) {
    debug_assert_eq!(a.len(), b.len());
    let (mut common, mut na, mut nb) = (0u64, 0u64, 0u64);

    let mut iter_a = a.chunks_exact(8);
    let mut iter_b = b.chunks_exact(8);
    for (ca, cb) in (&mut iter_a).zip(&mut iter_b) {
        let x = u64::from_le_bytes(ca.try_into().unwrap());
        let y = u64::from_le_bytes(cb.try_into().unwrap());
        common += (x & y).count_ones() as u64;
        na += x.count_ones() as u64;
        nb += y.count_ones() as u64;
    }
    for (xa, xb) in iter_a.remainder().iter().zip(iter_b.remainder()) {
        common += (xa & xb).count_ones() as u64;
        na += xa.count_ones() as u64;
        nb += xb.count_ones() as u64;
    }
    (common, na, nb)
}

/// Dice: `2*common / (na + nb)`. Returns 0.0 if both vectors are empty.
/// (`BitOps.cpp` `DiceSimilarity`.)
pub fn dice_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    if y + z > 0 {
        (2 * x) as f64 / (y + z) as f64
    } else {
        0.0
    }
}

/// Cosine (Ochiai): `common / sqrt(na * nb)`. Returns 0.0 if either is empty.
/// (`BitOps.cpp` `CosineSimilarity`.)
pub fn cosine_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    if y * z > 0 {
        x as f64 / ((y * z) as f64).sqrt()
    } else {
        0.0
    }
}

/// Kulczynski: `common * (na + nb) / (2 * na * nb)`. Returns 0.0 if either is
/// empty. (`BitOps.cpp` `KulczynskiSimilarity`.)
pub fn kulczynski_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    if y * z > 0 {
        (x * (y + z)) as f64 / (2 * y * z) as f64
    } else {
        0.0
    }
}

/// Sokal: `common / (2*na + 2*nb - 3*common)`. Returns 0.0 if either is empty.
/// (`BitOps.cpp` `SokalSimilarity`.)
pub fn sokal_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    if y == 0 || z == 0 {
        return 0.0;
    }
    x as f64 / (2.0 * y as f64 + 2.0 * z as f64 - 3.0 * x as f64)
}

/// McConnaughey: `(common*(na+nb) - na*nb) / (na*nb)`, range `[-1, 1]`.
/// Returns 0.0 if either is empty. (`BitOps.cpp` `McConnaugheySimilarity`.)
pub fn mcconnaughey_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    if y * z > 0 {
        (x as f64 * (y + z) as f64 - (y * z) as f64) / (y * z) as f64
    } else {
        0.0
    }
}

/// Asymmetric: `common / min(na, nb)`. Returns 0.0 if either is empty.
/// (`BitOps.cpp` `AsymmetricSimilarity`.)
pub fn asymmetric_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    let min = y.min(z);
    if min > 0 {
        x as f64 / min as f64
    } else {
        0.0
    }
}

/// Braun-Blanquet: `common / max(na, nb)`. Returns 0.0 if both are empty.
/// (`BitOps.cpp` `BraunBlanquetSimilarity`.)
pub fn braun_blanquet_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, y, z) = bit_counts(a, b);
    let max = y.max(z);
    if max > 0 {
        x as f64 / max as f64
    } else {
        0.0
    }
}

/// Russel: `common / total_bits` (denominator is the full vector width, not the
/// on-bit count). (`BitOps.cpp` `RusselSimilarity`.)
pub fn russel_bit(a: &[u8], b: &[u8]) -> f64 {
    let (x, _, _) = bit_counts(a, b);
    let total_bits = (a.len() * 8) as f64;
    if total_bits > 0.0 {
        x as f64 / total_bits
    } else {
        0.0
    }
}

/// Tversky with asymmetric weights `alpha`, `beta` (both in `[0, 1]`):
/// `common / (alpha*na + beta*nb + (1-alpha-beta)*common)`.
///
/// Matches `BitOps.cpp` `TverskySimilarity`: returns 0.0 if either vector is
/// empty, and 1.0 if the denominator evaluates to exactly zero. Returns `NaN`
/// if `alpha` or `beta` is outside `[0, 1]` (the FFI surfaces this as an error;
/// RDKit throws a `RANGE_CHECK` here).
pub fn tversky_bit(a: &[u8], b: &[u8], alpha: f64, beta: f64) -> f64 {
    if !(0.0..=1.0).contains(&alpha) || !(0.0..=1.0).contains(&beta) {
        return f64::NAN;
    }
    let (x, y, z) = bit_counts(a, b);
    if y == 0 || z == 0 {
        return 0.0;
    }
    let denom = alpha * y as f64 + beta * z as f64 + (1.0 - alpha - beta) * x as f64;
    if denom == 0.0 {
        1.0
    } else {
        x as f64 / denom
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // a = 0b1111_0000 (na=4), b = 0b0011_1100 (nb=4), common = 0b0011_0000 (2).
    const A: &[u8] = &[0b1111_0000];
    const B: &[u8] = &[0b0011_1100];

    fn approx(x: f64, y: f64) {
        assert!((x - y).abs() < 1e-12, "{x} != {y}");
    }

    #[test]
    fn dice_half_overlap() {
        // 2*2 / (4+4) = 0.5
        approx(dice_bit(A, B), 0.5);
    }

    #[test]
    fn cosine_half_overlap() {
        // 2 / sqrt(4*4) = 2/4 = 0.5
        approx(cosine_bit(A, B), 0.5);
    }

    #[test]
    fn kulczynski_half_overlap() {
        // 2*(4+4) / (2*4*4) = 16/32 = 0.5
        approx(kulczynski_bit(A, B), 0.5);
    }

    #[test]
    fn sokal_half_overlap() {
        // 2 / (8 + 8 - 6) = 2/10 = 0.2
        approx(sokal_bit(A, B), 0.2);
    }

    #[test]
    fn mcconnaughey_half_overlap() {
        // (2*8 - 16) / 16 = 0/16 = 0.0
        approx(mcconnaughey_bit(A, B), 0.0);
    }

    #[test]
    fn asymmetric_and_braun_blanquet() {
        // common/min = 2/4 = 0.5; common/max = 2/4 = 0.5 (na == nb here)
        approx(asymmetric_bit(A, B), 0.5);
        approx(braun_blanquet_bit(A, B), 0.5);
    }

    #[test]
    fn russel_uses_total_width() {
        // common=2, total bits = 8 → 0.25
        approx(russel_bit(A, B), 0.25);
    }

    #[test]
    fn tversky_reduces_to_tanimoto_at_one_one() {
        // alpha=beta=1: 2 / (4 + 4 - 2) = 2/6 = 1/3 (== Tanimoto)
        approx(tversky_bit(A, B, 1.0, 1.0), 1.0 / 3.0);
    }

    #[test]
    fn tversky_reduces_to_dice_at_half_half() {
        // alpha=beta=0.5: 2 / (2 + 2 + 0) = 0.5 (== Dice)
        approx(tversky_bit(A, B, 0.5, 0.5), 0.5);
    }

    #[test]
    fn tversky_substructure_weighting() {
        // alpha=0, beta=1: 2 / (0 + 4 + (0)*2) = 2/4 = 0.5
        approx(tversky_bit(A, B, 0.0, 1.0), 0.5);
    }

    #[test]
    fn tversky_out_of_range_is_nan() {
        assert!(tversky_bit(A, B, 1.5, 0.0).is_nan());
        assert!(tversky_bit(A, B, 0.0, -0.1).is_nan());
    }

    #[test]
    fn empty_vectors_degenerate_to_zero() {
        let z = &[0u8; 4][..];
        approx(dice_bit(z, z), 0.0);
        approx(cosine_bit(z, z), 0.0);
        approx(kulczynski_bit(z, z), 0.0);
        approx(sokal_bit(z, z), 0.0);
        approx(mcconnaughey_bit(z, z), 0.0);
        approx(asymmetric_bit(z, z), 0.0);
        approx(braun_blanquet_bit(z, z), 0.0);
        approx(russel_bit(z, z), 0.0);
        approx(tversky_bit(z, z, 1.0, 1.0), 0.0);
    }

    #[test]
    fn identical_vectors() {
        let v = &[0b1010_1010u8; 32][..];
        approx(dice_bit(v, v), 1.0);
        approx(cosine_bit(v, v), 1.0);
        approx(sokal_bit(v, v), 1.0);
        approx(asymmetric_bit(v, v), 1.0);
        approx(braun_blanquet_bit(v, v), 1.0);
        approx(tversky_bit(v, v, 0.3, 0.7), 1.0);
    }

    #[test]
    fn tail_bytes_handled() {
        // 9-byte vectors exercise the chunks_exact remainder path.
        let a = vec![0xFFu8; 9];
        let b = vec![0x0Fu8; 9];
        // na=72, nb=36, common=36 → dice = 72/108 = 2/3
        approx(dice_bit(&a, &b), 2.0 / 3.0);
    }
}
