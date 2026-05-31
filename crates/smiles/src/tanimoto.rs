//! Tanimoto similarity over raw bit-fingerprint BLOBs.
//!
//! Faithful port of the algorithm in `RDKit/Code/DataStructs/BitOps.cpp`
//! (`CalcBitmapTanimoto`):
//!
//! ```text
//! Tanimoto(a, b) = popcount(a AND b) / popcount(a OR b)
//! ```
//!
//! Two popcounts are computed in a single loop (AND and OR over the same
//! bytes), matching RDKit's optimization.
//!
//! RDKit dispatches between a 256-entry lookup table (`byte_popcounts[]`)
//! and `__builtin_popcountll` via `#ifdef RDK_OPTIMIZE_POPCNT`. In Rust,
//! `u64::count_ones()` already lowers to POPCNT (x86_64) or CNT (aarch64)
//! when the target supports it, with a SWAR fallback otherwise — so we
//! don't carry our own lookup table.

/// Tanimoto similarity between two equal-length byte slices interpreted
/// as bit vectors.
///
/// - Returns `0.0` if both vectors are all-zeros (RDKit convention — avoids
///   the `0 / 0` case and lets downstream `ORDER BY sim DESC` keep going).
/// - Caller must ensure `a.len() == b.len()`; the FFI wrapper enforces this
///   and the DuckDB extension surfaces a clear error to the user.
pub fn tanimoto_bit(a: &[u8], b: &[u8]) -> f64 {
    debug_assert_eq!(a.len(), b.len());

    let mut intersect: u64 = 0;
    let mut union_: u64 = 0;

    // Process 8 bytes at a time as u64 — one POPCNT per chunk.
    let mut iter_a = a.chunks_exact(8);
    let mut iter_b = b.chunks_exact(8);
    for (ca, cb) in (&mut iter_a).zip(&mut iter_b) {
        let x = u64::from_le_bytes(ca.try_into().unwrap());
        let y = u64::from_le_bytes(cb.try_into().unwrap());
        intersect += (x & y).count_ones() as u64;
        union_ += (x | y).count_ones() as u64;
    }

    // Tail bytes for fingerprint widths that aren't multiples of 64 bits
    // (256 bytes / 2048 bits has no tail; this is a safety net for
    // 1024-bit or other variants we may add later).
    for (xa, xb) in iter_a.remainder().iter().zip(iter_b.remainder()) {
        intersect += (xa & xb).count_ones() as u64;
        union_ += (xa | xb).count_ones() as u64;
    }

    if union_ == 0 {
        0.0
    } else {
        intersect as f64 / union_ as f64
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn identical_blobs_give_one() {
        let a = vec![0b10101010_u8; 256];
        assert_eq!(tanimoto_bit(&a, &a), 1.0);
    }

    #[test]
    fn empty_blobs_give_zero() {
        let a = vec![0u8; 256];
        let b = vec![0u8; 256];
        assert_eq!(tanimoto_bit(&a, &b), 0.0);
    }

    #[test]
    fn disjoint_bits_give_zero() {
        // a has only even-index bits, b has only odd-index bits → no overlap
        let a = vec![0b01010101_u8; 256];
        let b = vec![0b10101010_u8; 256];
        assert_eq!(tanimoto_bit(&a, &b), 0.0);
    }

    #[test]
    fn half_overlap_gives_one_third() {
        // a = 0b1111_0000 (4 bits), b = 0b0011_1100 (4 bits)
        //   AND = 0b0011_0000 (2 bits), OR = 0b1111_1100 (6 bits)
        //   → 2/6 = 1/3
        let a = vec![0b1111_0000_u8];
        let b = vec![0b0011_1100_u8];
        let t = tanimoto_bit(&a, &b);
        assert!((t - 1.0 / 3.0).abs() < 1e-12);
    }

    #[test]
    fn symmetric() {
        let a = vec![0b11001100_u8; 32];
        let b = vec![0b10101010_u8; 32];
        assert_eq!(tanimoto_bit(&a, &b), tanimoto_bit(&b, &a));
    }

    #[test]
    fn in_range_zero_to_one() {
        for seed in 0..16u8 {
            let a: Vec<u8> = (0..256)
                .map(|i| (i as u8).wrapping_mul(seed).wrapping_add(seed))
                .collect();
            let b: Vec<u8> = (0..256)
                .map(|i| {
                    (i as u8)
                        .wrapping_mul(seed.wrapping_add(7))
                        .wrapping_add(13)
                })
                .collect();
            let t = tanimoto_bit(&a, &b);
            assert!((0.0..=1.0).contains(&t), "out of range: {t}");
        }
    }

    #[test]
    fn tail_bytes_handled() {
        // 9-byte blobs (not multiple of 8) — exercise the remainder loop.
        let a = vec![0xFF_u8; 9];
        let b = vec![0xFF_u8; 9];
        assert_eq!(tanimoto_bit(&a, &b), 1.0);
    }
}
