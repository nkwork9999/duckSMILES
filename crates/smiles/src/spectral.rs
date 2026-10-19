//! Dense symmetric extrema: scaled Householder reduction and Sturm bisection.
//! No eigenvectors/full spectrum are computed. The recurrence is the standard
//! symmetric tridiagonal LDL/Sturm count (see LAPACK DLAEBZ documentation).

pub(crate) fn extrema(matrix: &[f64], n: usize) -> Option<(f64, f64)> {
    if n == 0 || matrix.len() != n.checked_mul(n)? || matrix.iter().any(|x| !x.is_finite()) {
        return None;
    }
    for i in 0..n {
        for j in 0..i {
            if matrix[i * n + j] != matrix[j * n + i] {
                return None;
            }
        }
    }
    let scale = matrix.iter().map(|x| x.abs()).fold(0., f64::max);
    if scale == 0. {
        return Some((0., 0.));
    }
    if n == 1 {
        return Some((matrix[0], matrix[0]));
    }
    let mut a: Vec<f64> = matrix.iter().map(|x| x / scale).collect();
    let mut e = vec![0.; n - 1];
    for k in 0..n - 1 {
        let start = k + 1;
        let size = n - start;
        let mut v: Vec<f64> = (start..n).map(|i| a[i * n + k]).collect();
        let norm = v.iter().fold(0f64, |sum, &x| sum.hypot(x));
        if norm == 0. {
            continue;
        }
        let sign = if v[0] < 0. { -1. } else { 1. };
        e[k] = -sign * norm;
        for x in &mut v {
            *x /= norm;
        }
        v[0] += sign;
        let vn = v.iter().fold(0f64, |sum, &x| sum.hypot(x));
        for x in &mut v {
            *x /= vn;
        }
        let mut w = vec![0.; size];
        for i in 0..size {
            // A is exactly symmetric. Read each off-diagonal element once,
            // sharing its factor of two between the two row products.
            // Each w row still receives columns in increasing order, so this
            // preserves the original floating-point accumulation order.
            let offset = (start + i) * n + start;
            let (previous, current) = w.split_at_mut(i);
            let mut sum = 0.;
            for ((&value, &vj), wj) in a[offset..offset + i].iter().zip(&v[..i]).zip(previous) {
                let twice = 2. * value;
                sum += twice * vj;
                *wj += twice * v[i];
            }
            current[0] = sum + 2. * a[offset + i] * v[i];
        }
        let projection = v.iter().zip(&w).map(|(v, w)| v * w).sum::<f64>();
        for i in 0..size {
            w[i] -= projection * v[i];
        }
        for i in 0..size {
            // All subsequent reduction reads use the lower triangle. Avoid
            // writing a mirrored upper half which is never read again.
            let offset = (start + i) * n + start;
            for ((value, &vj), &wj) in a[offset..offset + i + 1]
                .iter_mut()
                .zip(&v[..=i])
                .zip(&w[..=i])
            {
                *value = *value - v[i] * wj - w[i] * vj;
            }
        }
    }
    let d: Vec<f64> = (0..n).map(|i| a[i * n + i]).collect();
    if d.iter().chain(&e).any(|x| !x.is_finite()) {
        return None;
    }
    let mut low = f64::INFINITY;
    let mut high = f64::NEG_INFINITY;
    for i in 0..n {
        let radius =
            if i > 0 { e[i - 1].abs() } else { 0. } + if i + 1 < n { e[i].abs() } else { 0. };
        low = low.min(d[i] - radius);
        high = high.max(d[i] + radius);
    }
    let margin = (low.abs() + high.abs() + 1.) * n as f64 * f64::EPSILON * 4.;
    low -= margin;
    high += margin;
    let squares: Vec<f64> = e.iter().map(|x| x * x).collect();
    let pivot_min = f64::MIN_POSITIVE * squares.iter().copied().fold(1., f64::max);
    let count = |x: f64| {
        let mut pivot = d[0] - x;
        if pivot.abs() < pivot_min {
            pivot = -pivot_min;
        }
        let mut count = usize::from(pivot < 0.);
        for i in 1..n {
            pivot = d[i] - x - squares[i - 1] / pivot;
            if pivot.abs() < pivot_min {
                pivot = -pivot_min;
            }
            count += usize::from(pivot < 0.);
        }
        count
    };
    let find = |target: usize| -> Option<f64> {
        let (mut lo, mut hi) = (low, high);
        for _ in 0..128 {
            let mid = lo + (hi - lo) * 0.5;
            if hi - lo <= 8. * f64::EPSILON * lo.abs().max(hi.abs()).max(1.)
                || mid == lo
                || mid == hi
            {
                return Some(mid * scale);
            }
            if count(mid) >= target {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        None
    };
    let lo = find(1)?;
    let hi = find(n)?;
    (lo.is_finite() && hi.is_finite()).then_some((lo, hi))
}
#[cfg(test)]
mod tests {
    use super::*;
    fn input(family: &str, n: usize, seed: usize) -> Vec<f64> {
        (0..n)
            .flat_map(|i| {
                (0..n).map(move |j| match family {
                    "diagonal" => {
                        if i == j {
                            ((i * 13 + seed) % 17) as f64 - 8.
                        } else {
                            0.
                        }
                    }
                    "rankone" => {
                        if i == j {
                            5.
                        } else {
                            0.001
                        }
                    }
                    "path" => {
                        if i == j {
                            12.011
                        } else if i.abs_diff(j) == 1 {
                            1.
                        } else {
                            0.
                        }
                    }
                    "star" => {
                        if i == j {
                            0.
                        } else if i.min(j) == 0 {
                            1.
                        } else {
                            0.
                        }
                    }
                    "clustered" => {
                        if i == j {
                            1. + (i % 3) as f64 * 1e-10
                        } else {
                            1e-12
                        }
                    }
                    "random" => {
                        if i == j {
                            ((i * 13 + seed) % 50 + 1) as f64
                        } else {
                            ((i.min(j) * 131 + i.max(j) * 17 + seed * 31) % 101) as f64 / 250. - 0.2
                        }
                    }
                    _ => panic!("unknown independent matrix family"),
                })
            })
            .collect()
    }
    #[test]
    fn extrema_match_numpy_for_dense_sparse_degenerate_and_clustered_matrices() {
        for line in include_str!("../tests/fixtures/symmetric_extrema.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let n: usize = f[1].parse().unwrap();
            let seed: usize = f[2].parse().unwrap();
            let expected = (f[3].parse::<f64>().unwrap(), f[4].parse::<f64>().unwrap());
            let matrix = input(f[0], n, seed);
            let actual =
                extrema(&matrix, n).unwrap_or_else(|| panic!("no extrema for {} n{n}", f[0]));
            assert!(
                (actual.0 - expected.0).abs() <= 1e-9 && (actual.1 - expected.1).abs() <= 1e-9,
                "{} n{n}: {actual:?} vs {expected:?}",
                f[0]
            );
        }
    }
    #[test]
    fn invalid_shapes_nonfinite_and_asymmetric_inputs_have_no_result() {
        assert!(extrema(&[], 0).is_none());
        assert!(extrema(&[1., 2.], 2).is_none());
        assert!(extrema(&[f64::NAN], 1).is_none());
        assert!(extrema(&[f64::INFINITY], 1).is_none());
        assert!(extrema(&[1., 2., 3., 4.], 2).is_none());
    }
    #[test]
    fn scalar_shift_scale_permutation_and_zero_invariants() {
        let original = input("random", 17, 23);
        let reference = extrema(&original, 17).unwrap();
        for scale in [1e-6, 1., 100.] {
            let m: Vec<_> = original.iter().map(|x| x * scale).collect();
            let actual = extrema(&m, 17).unwrap();
            assert!(
                (actual.0 - reference.0 * scale).abs() <= 1e-9
                    && (actual.1 - reference.1 * scale).abs() <= 1e-9
            );
        }
        let mut shifted = original.clone();
        for i in 0..17 {
            shifted[i * 17 + i] += 100.;
        }
        let actual = extrema(&shifted, 17).unwrap();
        assert!(
            (actual.0 - reference.0 - 100.).abs() <= 1e-9
                && (actual.1 - reference.1 - 100.).abs() <= 1e-9
        );
        let mut reversed = original.clone();
        for i in 0..17 {
            for j in 0..17 {
                reversed[i * 17 + j] = original[(16 - i) * 17 + (16 - j)];
            }
        }
        let actual = extrema(&reversed, 17).unwrap();
        assert!((actual.0 - reference.0).abs() <= 1e-9 && (actual.1 - reference.1).abs() <= 1e-9);
        assert_eq!(extrema(&vec![0.; 64], 8), Some((0., 0.)));
    }
}
