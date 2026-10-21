//! Additional bounded experimental designs and linear-model design diagnostics.
//!
//! Bounds are one [low, high] row per factor. Box–Behnken uses all four signs for
//! each factor pair, with other factors centered. Fractional factorial generator
//! masks multiply selected base columns; no resolution or alias-freedom claim
//! is inferred. Orthogonal screening is a Sylvester Hadamard main-effect design,
//! not general Plackett–Burman or definitive screening.
//!
//! Design evaluation adds an intercept to the supplied model terms. Raw-unit
//! information/inverse/log determinant are accompanied by centered, RMS-scaled
//! information and a coded-basis spectral condition number. Rank failure and
//! unrepresentable diagnostics return errors rather than a pseudoinverse.

use std::collections::BTreeSet;

const MAX_ROWS: usize = 100_000;
const MAX_CELLS: usize = 1_000_000;
const MAX_FACTORS: usize = 64;
const MAX_WORK: usize = 50_000_000;
const MAX_SEED: f64 = 9_007_199_254_740_991.;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > MAX_ROWS
        || width == 0
        || width > 1024
        || data.len() > MAX_CELLS / width
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("design input must be finite rectangular data within 100000 rows, 1024 columns and 1000000 cells".into());
    }
    match operation {
        "box_behnken" => box_behnken(data, args),
        "latin_hypercube" => latin_hypercube(data, args),
        "mixture_lattice" => mixture_lattice(data, args),
        "fractional_factorial" => fractional_factorial(data, args),
        "orthogonal_screening" => orthogonal_screening(data, args),
        "randomize_design" => randomize(data, args),
        "design_evaluate" => evaluate(data, args),
        _ => Err(format!("unknown advanced design operation: {operation}")),
    }
}
fn arguments(args: &[f64], min: usize, max: usize) -> Result<(), String> {
    if (min..=max).contains(&args.len()) {
        Ok(())
    } else {
        Err(format!("operation requires {min} to {max} arguments"))
    }
}
fn integer(x: f64, min: usize, max: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < min as f64 || x > max as f64 {
        Err(format!("{name} must be an integer in [{min},{max}]"))
    } else {
        Ok(x as usize)
    }
}
fn seed(x: f64) -> Result<u64, String> {
    if !(0. ..=MAX_SEED).contains(&x) || x.fract() != 0. {
        Err("seed must be an exact nonnegative integer no greater than 2^53-1".into())
    } else {
        Ok(x as u64)
    }
}
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn vector(x: &[f64]) -> Result<String, String> {
    if x.iter().any(|v| !v.is_finite()) {
        return Err("design output exceeds finite numerical range".into());
    }
    Ok(format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    ))
}
fn matrix(x: &[Vec<f64>]) -> Result<String, String> {
    Ok(format!(
        "[{}]",
        x.iter()
            .map(|r| vector(r))
            .collect::<Result<Vec<_>, _>>()?
            .join(",")
    ))
}
fn index_vector(x: &[usize]) -> String {
    format!(
        "[{}]",
        x.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
    )
}
fn limit_output(n: usize, p: usize) -> Result<(), String> {
    if n == 0 || n > MAX_ROWS || n > MAX_CELLS / p {
        Err("generated design exceeds 100000 runs or 1000000 cells".into())
    } else {
        Ok(())
    }
}
fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.checked_mul(b)
        .and_then(|v| v.checked_mul(c))
        .is_none_or(|v| v > MAX_WORK)
    {
        Err("design computation exceeds 50000000 work units".into())
    } else {
        Ok(())
    }
}
fn strict_bounds(bounds: &[Vec<f64>]) -> Result<(Vec<f64>, Vec<f64>), String> {
    if bounds.len() > MAX_FACTORS || bounds[0].len() != 2 || bounds.iter().any(|r| r[0] >= r[1]) {
        return Err("design bounds require 1 to 64 [low,high] rows with low < high".into());
    }
    let mid = bounds
        .iter()
        .map(|r| r[0] * 0.5 + r[1] * 0.5)
        .collect::<Vec<_>>();
    let half = bounds
        .iter()
        .map(|r| r[1] * 0.5 - r[0] * 0.5)
        .collect::<Vec<_>>();
    if mid
        .iter()
        .zip(&half)
        .zip(bounds)
        .any(|((&m, &h), b)| !m.is_finite() || !h.is_finite() || h <= 0. || m <= b[0] || m >= b[1])
    {
        return Err("factor bounds cannot represent distinct low, center and high levels".into());
    }
    Ok((mid, half))
}
fn physical(bounds: &[Vec<f64>], coded: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let (mid, half) = strict_bounds(bounds)?;
    coded
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .map(|(j, &z)| {
                    if z == -1. {
                        Ok(bounds[j][0])
                    } else if z == 1. {
                        Ok(bounds[j][1])
                    } else {
                        finite(mid[j] + half[j] * z, "factor level")
                    }
                })
                .collect()
        })
        .collect()
}
fn design_result(
    name: &str,
    bounds: &[Vec<f64>],
    coded: &[Vec<f64>],
    extra: &str,
) -> Result<String, String> {
    let rows = physical(bounds, coded)?;
    Ok(format!("{{\"version\":1,\"operation\":\"{name}\",\"factors\":{},\"run_count\":{},\"rows\":{},\"coded_rows\":{}{extra}}}",bounds.len(),rows.len(),matrix(&rows)?,matrix(coded)?))
}
fn box_behnken(bounds: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 1)?;
    strict_bounds(bounds)?;
    let p = bounds.len();
    if p < 3 {
        return Err("Box-Behnken requires at least three factors".into());
    }
    let centers = integer(
        args.first().copied().unwrap_or(1.),
        0,
        MAX_ROWS,
        "center points",
    )?;
    let count = 2 * p * (p - 1) + centers;
    limit_output(count, p)?;
    let mut coded = Vec::with_capacity(count);
    for a in 0..p {
        for b in a + 1..p {
            for u in [-1., 1.] {
                for v in [-1., 1.] {
                    let mut row = vec![0.; p];
                    row[a] = u;
                    row[b] = v;
                    coded.push(row);
                }
            }
        }
    }
    coded.extend((0..centers).map(|_| vec![0.; p]));
    design_result(
        "box_behnken",
        bounds,
        &coded,
        &format!(",\"center_points\":{centers},\"construction\":\"all_two_factor_edges\""),
    )
}
struct Random(u64);
impl Random {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn unit(&mut self) -> f64 {
        (self.next() >> 11) as f64 / 9_007_199_254_740_992.
    }
    fn index(&mut self, n: usize) -> usize {
        let n = n as u64;
        let threshold = n.wrapping_neg() % n;
        loop {
            let v = self.next();
            if v >= threshold {
                return (v % n) as usize;
            }
        }
    }
    fn shuffle(&mut self, values: &mut [usize]) {
        for i in (1..values.len()).rev() {
            let j = self.index(i + 1);
            values.swap(i, j);
        }
    }
}
fn latin_hypercube(bounds: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 2, 2)?;
    strict_bounds(bounds)?;
    let n = integer(args[0], 1, MAX_ROWS, "Latin hypercube runs")?;
    let initial = seed(args[1])?;
    let p = bounds.len();
    limit_output(n, p)?;
    let mut rng = Random(initial);
    let mut unit = vec![vec![0.; p]; n];
    let mut strata = vec![vec![0.; p]; n];
    for j in 0..p {
        let mut permutation = (0..n).collect::<Vec<_>>();
        rng.shuffle(&mut permutation);
        for i in 0..n {
            unit[i][j] = (permutation[i] as f64 + rng.unit()) / n as f64;
            strata[i][j] = permutation[i] as f64;
        }
    }
    let coded = unit
        .iter()
        .map(|r| r.iter().map(|v| 2. * v - 1.).collect())
        .collect::<Vec<_>>();
    design_result(
        "latin_hypercube",
        bounds,
        &coded,
        &format!(
            ",\"seed\":{initial},\"jitter\":true,\"unit_rows\":{},\"stratum_indices\":{}",
            matrix(&unit)?,
            matrix(&strata)?
        ),
    )
}
fn combination_count(n: usize, k: usize, cap: usize) -> Result<usize, String> {
    let k = k.min(n - k);
    let mut value = 1usize;
    for i in 0..k {
        value = value
            .checked_mul(n - i)
            .ok_or("mixture composition count overflow")?
            / (i + 1);
        if value > cap {
            return Err("mixture enumeration exceeds 1000000 compositions".into());
        }
    }
    Ok(value)
}
fn mixture_lattice(bounds: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 1, 2)?;
    let p = bounds.len();
    if p > MAX_FACTORS || bounds[0].len() != 2 || bounds.iter().any(|r| r[0] < 0. || r[0] > r[1]) {
        return Err(
            "mixture bounds require 1 to 64 [nonnegative lower, upper] rows with lower <= upper"
                .into(),
        );
    }
    let degree = integer(args[0], 1, 1000, "mixture lattice degree")?;
    let total = args.get(1).copied().unwrap_or(1.);
    if total <= 0. {
        return Err("mixture total must be positive".into());
    }
    let lower = finite(
        bounds.iter().map(|r| r[0]).sum(),
        "sum of mixture lower bounds",
    )?;
    let upper = finite(
        bounds.iter().map(|r| r[1]).sum(),
        "sum of mixture upper bounds",
    )?;
    let tolerance = (64. * f64::EPSILON * p as f64) * total;
    if lower > total + tolerance || upper < total - tolerance {
        return Err("mixture lower/upper totals make the simplex infeasible".into());
    }
    let remaining = (total - lower).max(0.);
    let mut rows = Vec::new();
    let mut compositions = Vec::new();
    let mut visited = 0usize;
    let mut seen = BTreeSet::new();
    if remaining == 0. {
        let row = bounds.iter().map(|r| r[0]).collect::<Vec<_>>();
        if (row.iter().sum::<f64>() - total).abs() > tolerance {
            return Err("fixed mixture does not sum to total".into());
        }
        let mut composition = vec![0; p];
        composition[p - 1] = degree;
        rows.push(row);
        compositions.push(composition);
        visited = 1;
    } else {
        combination_count(degree + p - 1, p - 1, 1_000_000)?;
        struct Grid<'a> {
            bounds: &'a [Vec<f64>],
            degree: usize,
            total: f64,
            remaining: f64,
            tolerance: f64,
            rows: &'a mut Vec<Vec<f64>>,
            compositions: &'a mut Vec<Vec<usize>>,
            visited: &'a mut usize,
            seen: &'a mut BTreeSet<Vec<u64>>,
        }
        impl Grid<'_> {
            fn visit(
                &mut self,
                counts: &mut [usize],
                position: usize,
                left: usize,
            ) -> Result<(), String> {
                if position + 1 < counts.len() {
                    for k in 0..=left {
                        counts[position] = k;
                        self.visit(counts, position + 1, left - k)?;
                    }
                    return Ok(());
                }
                counts[position] = left;
                *self.visited += 1;
                let p = counts.len();
                let mut row = counts
                    .iter()
                    .zip(self.bounds)
                    .map(|(&k, b)| b[0] + self.remaining * (k as f64 / self.degree as f64))
                    .collect::<Vec<_>>();
                row[p - 1] = self.total - row[..p - 1].iter().sum::<f64>();
                if row.iter().zip(self.bounds).any(|(&v, b)| {
                    !v.is_finite() || v < b[0] - self.tolerance || v > b[1] + self.tolerance
                }) {
                    return Ok(());
                }
                for (v, b) in row.iter_mut().zip(self.bounds) {
                    *v = v.clamp(b[0], b[1]);
                }
                if (row.iter().sum::<f64>() - self.total).abs() > self.tolerance {
                    return Err("mixture sum lost numerical precision".into());
                }
                let bits = row.iter().map(|v| v.to_bits()).collect::<Vec<_>>();
                if !self.seen.insert(bits) {
                    return Err(
                        "mixture lattice points collapse at the supplied numerical scale".into(),
                    );
                }
                limit_output(self.rows.len() + 1, p)?;
                self.rows.push(row);
                self.compositions.push(counts.to_vec());
                Ok(())
            }
        }
        let mut grid = Grid {
            bounds,
            degree,
            total,
            remaining,
            tolerance,
            rows: &mut rows,
            compositions: &mut compositions,
            visited: &mut visited,
            seen: &mut seen,
        };
        grid.visit(&mut vec![0; p], 0, degree)?;
    }
    if rows.is_empty() {
        return Err("no simplex lattice points satisfy the supplied bounds".into());
    }
    let max_error = rows
        .iter()
        .map(|r| (r.iter().sum::<f64>() - total).abs())
        .fold(0., f64::max);
    let counts = compositions
        .iter()
        .map(|r| r.iter().map(|&v| v as f64).collect())
        .collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"operation\":\"mixture_lattice\",\"construction\":\"lower_bound_offset_simplex_lattice_with_upper_filter\",\"factors\":{p},\"degree\":{degree},\"total\":{total},\"remaining_total\":{remaining},\"enumerated_compositions\":{visited},\"run_count\":{},\"rows\":{},\"integer_compositions\":{},\"sum_tolerance\":{tolerance},\"maximum_sum_error\":{max_error}}}",rows.len(),matrix(&rows)?,matrix(&counts)?))
}
fn fractional_factorial(bounds: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 2, MAX_FACTORS + 2)?;
    strict_bounds(bounds)?;
    let p = bounds.len();
    let base = integer(args[0], 1, p.min(16), "base factors")?;
    let centers = integer(args[1], 0, MAX_ROWS, "center points")?;
    if args.len() != 2 + p - base {
        return Err(
            "fractional factorial requires exactly one generator mask per nonbase factor".into(),
        );
    }
    let n = 1usize << base;
    if p > n - 1 {
        return Err("not enough distinct nonconstant product columns for requested factors".into());
    }
    limit_output(n + centers, p)?;
    let mut masks = (0..base).map(|j| 1usize << j).collect::<Vec<_>>();
    let mut seen = masks.iter().copied().collect::<BTreeSet<_>>();
    for &v in &args[2..] {
        let mask = integer(v, 1, n - 1, "generator mask")?;
        if mask.count_ones() < 2 || !seen.insert(mask) {
            return Err(
                "generator masks must be distinct products of at least two base factors".into(),
            );
        }
        masks.push(mask);
    }
    let mut coded = (0..n)
        .map(|run| {
            masks
                .iter()
                .map(|&mask| {
                    if (mask.count_ones() - (run & mask).count_ones()).is_multiple_of(2) {
                        1.
                    } else {
                        -1.
                    }
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    coded.extend((0..centers).map(|_| vec![0.; p]));
    design_result("fractional_factorial",bounds,&coded,&format!(",\"base_factors\":{base},\"center_points\":{centers},\"generator_masks\":{},\"construction\":\"positive_product_generators\",\"factor_order\":\"base_columns_then_generators\",\"alias_analysis_performed\":false",index_vector(&masks[base..])))
}
fn orthogonal_screening(bounds: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 0)?;
    strict_bounds(bounds)?;
    let p = bounds.len();
    let n = (p + 1).next_power_of_two();
    limit_output(n, p)?;
    let coded = (0..n)
        .map(|i| {
            (1..=p)
                .map(|j| {
                    if (i & j).count_ones().is_multiple_of(2) {
                        1.
                    } else {
                        -1.
                    }
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    design_result("orthogonal_screening",bounds,&coded,",\"construction\":\"sylvester_hadamard_main_effect_screening\",\"interactions_may_be_aliased\":true")
}
fn randomize(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 2, 3)?;
    let replicates = integer(args[0], 1, 1000, "replicates")?;
    let initial = seed(args[1])?;
    let blocks = integer(
        args.get(2).copied().unwrap_or(1.),
        1,
        data.len(),
        "blocks per replicate",
    )?;
    if !data.len().is_multiple_of(blocks) {
        return Err("input run count must divide evenly into blocks for each replicate".into());
    }
    let count = data
        .len()
        .checked_mul(replicates)
        .ok_or("replicate count overflow")?;
    limit_output(count, data[0].len())?;
    let mut rng = Random(initial);
    let mut rows = Vec::with_capacity(count);
    let mut source = Vec::with_capacity(count);
    let mut replicate = Vec::with_capacity(count);
    let mut block = Vec::with_capacity(count);
    for r in 0..replicates {
        let mut permutation = (0..data.len()).collect::<Vec<_>>();
        rng.shuffle(&mut permutation);
        for (position, i) in permutation.into_iter().enumerate() {
            rows.push(data[i].clone());
            source.push(i);
            replicate.push(r);
            block.push(position / (data.len() / blocks));
        }
    }
    Ok(format!("{{\"version\":1,\"operation\":\"randomize_design\",\"randomization\":\"independent_complete_shuffle_per_replicate_then_equal_size_blocks\",\"replicates\":{replicates},\"blocks_per_replicate\":{blocks},\"seed\":{initial},\"run_count\":{count},\"rows\":{},\"source_indices\":{},\"replicate_labels\":{},\"block_labels\":{}}}",matrix(&rows)?,index_vector(&source),index_vector(&replicate),index_vector(&block)))
}

fn inner(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(
        a.iter().zip(b).map(|(a, b)| a * b).sum(),
        "information inner product",
    )
}
fn centered_column(data: &[Vec<f64>], j: usize) -> Result<(f64, f64, f64, Vec<f64>), String> {
    let origin = data[0][j];
    let mut offset = 0.;
    for (i, r) in data.iter().enumerate() {
        let value = finite(r[j] - origin, "origin-relative model term")?;
        offset = finite(
            offset + (value - offset) / (i + 1) as f64,
            "model-term mean",
        )?;
    }
    let mut column = data
        .iter()
        .map(|r| finite((r[j] - origin) - offset, "centered model term"))
        .collect::<Result<Vec<_>, _>>()?;
    let largest = column.iter().map(|v| v.abs()).fold(0., f64::max);
    if largest == 0. {
        return Err("design matrix is rank deficient: a supplied term is constant".into());
    }
    let scale = finite(
        largest
            * (column.iter().map(|v| (v / largest).powi(2)).sum::<f64>() / data.len() as f64)
                .sqrt(),
        "model-term RMS scale",
    )?;
    if scale == 0. {
        return Err("model-term scale underflow".into());
    }
    for v in &mut column {
        *v /= scale;
    }
    Ok((origin, offset, scale, column))
}
fn spectral_condition(gram: &[Vec<f64>]) -> Result<f64, String> {
    let mut a = gram.to_vec();
    let p = a.len();
    let maximum = (0..p).map(|j| a[j][j]).fold(0., f64::max);
    let threshold = maximum * 1e-14;
    let mut converged = false;
    for _ in 0..64 {
        let mut largest = 0_f64;
        for i in 0..p {
            for j in i + 1..p {
                let off = a[i][j];
                largest = largest.max(off.abs());
                if off.abs() <= threshold {
                    continue;
                }
                let tau = (a[j][j] - a[i][i]) / (2. * off);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / ((-tau) + tau.hypot(1.))
                };
                let c = 1. / (1. + t * t).sqrt();
                let s = t * c;
                a[i][i] -= t * off;
                a[j][j] += t * off;
                a[i][j] = 0.;
                a[j][i] = 0.;
                let updates = a
                    .iter()
                    .enumerate()
                    .filter(|(k, _)| *k != i && *k != j)
                    .map(|(k, row)| (k, c * row[i] - s * row[j], s * row[i] + c * row[j]))
                    .collect::<Vec<_>>();
                for (k, first, second) in updates {
                    a[k][i] = first;
                    a[i][k] = first;
                    a[k][j] = second;
                    a[j][k] = second;
                }
            }
        }
        if largest <= threshold {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("coded information spectrum did not converge".into());
    }
    let lo = (0..p).map(|j| a[j][j]).fold(f64::INFINITY, f64::min);
    let hi = (0..p).map(|j| a[j][j]).fold(0., f64::max);
    if lo <= hi * 1e-14 {
        return Err(
            "coded information is too ill-conditioned for reliable spectral diagnostics".into(),
        );
    }
    finite((hi / lo).sqrt(), "coded spectral condition number")
}
fn evaluate(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 0)?;
    let n = data.len();
    let p = data[0].len();
    let m = p + 1;
    if p > MAX_FACTORS || n < m {
        return Err(
            "design evaluation requires at most 64 model terms and at least terms+1 runs".into(),
        );
    }
    work(n, m, 4 * m)?;
    let mut columns = vec![vec![1.; n]];
    let mut origins = Vec::new();
    let mut offsets = Vec::new();
    let mut means = Vec::new();
    let mut scales = Vec::new();
    for j in 0..p {
        let (origin, offset, scale, column) = centered_column(data, j)?;
        origins.push(origin);
        offsets.push(offset);
        means.push(finite(origin + offset, "model-term mean")?);
        scales.push(scale);
        columns.push(column);
    }
    let mut q = Vec::<Vec<f64>>::new();
    let mut r = vec![vec![0.; m]; m];
    for j in 0..m {
        let mut v = columns[j].clone();
        let original = inner(&v, &v)?.sqrt();
        for _ in 0..2 {
            for i in 0..j {
                let amount = inner(&q[i], &v)?;
                r[i][j] += amount;
                for (v, q) in v.iter_mut().zip(&q[i]) {
                    *v -= amount * q;
                }
            }
        }
        let length = inner(&v, &v)?.sqrt();
        if length <= original * 1e-10 {
            return Err("design model matrix is numerically rank deficient".into());
        }
        r[j][j] = length;
        for v in &mut v {
            *v /= length;
        }
        q.push(v);
    }
    let coded_logdet = finite(
        2. * (0..m).map(|j| r[j][j].ln()).sum::<f64>(),
        "coded log determinant",
    )?;
    let raw_logdet = finite(
        coded_logdet + 2. * scales.iter().map(|v| v.ln()).sum::<f64>(),
        "raw log determinant",
    )?;
    let mut coded_gram = vec![vec![0.; m]; m];
    for i in 0..m {
        for j in 0..m {
            coded_gram[i][j] = inner(&columns[i], &columns[j])?;
        }
    }
    let condition = spectral_condition(&coded_gram)?;
    let inverse_columns = (0..m)
        .map(|j| {
            let mut column = vec![0.; m];
            for i in (0..m).rev() {
                let rhs = f64::from(i == j);
                column[i] = finite(
                    (rhs - (i + 1..m).map(|k| r[i][k] * column[k]).sum::<f64>()) / r[i][i],
                    "QR inverse",
                )?;
            }
            Ok(column)
        })
        .collect::<Result<Vec<_>, String>>()?;
    let rinv = (0..m)
        .map(|i| inverse_columns.iter().map(|c| c[i]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut coded_inverse = vec![vec![0.; m]; m];
    for i in 0..m {
        for j in 0..m {
            coded_inverse[i][j] = inner(&rinv[i], &rinv[j])?;
        }
    }
    let mut transform = vec![vec![0.; m]; m];
    transform[0][0] = 1.;
    for j in 0..p {
        transform[0][j + 1] = finite(-means[j] / scales[j], "raw model mean-to-scale ratio")?;
        transform[j + 1][j + 1] = finite(1. / scales[j], "inverse model scale")?;
    }
    let mut temp = vec![vec![0.; m]; m];
    for i in 0..m {
        for j in 0..m {
            temp[i][j] = finite(
                (0..m).map(|k| transform[i][k] * coded_inverse[k][j]).sum(),
                "raw inverse transform",
            )?;
        }
    }
    let mut raw_inverse = vec![vec![0.; m]; m];
    for i in 0..m {
        for j in 0..m {
            raw_inverse[i][j] = inner(&temp[i], &transform[j])?;
        }
    }
    let raw_columns = (0..m)
        .map(|j| {
            if j == 0 {
                vec![1.; n]
            } else {
                data.iter().map(|r| r[j - 1]).collect()
            }
        })
        .collect::<Vec<Vec<f64>>>();
    let mut raw_gram = vec![vec![0.; m]; m];
    for i in 0..m {
        for j in 0..m {
            raw_gram[i][j] = inner(&raw_columns[i], &raw_columns[j])?;
        }
    }
    let leverages = (0..n)
        .map(|i| {
            finite(
                (0..m).map(|j| q[j][i].powi(2)).sum::<f64>(),
                "design leverage",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let coded_rows = (0..n)
        .map(|i| (0..m).map(|j| columns[j][i]).collect())
        .collect::<Vec<Vec<f64>>>();
    let efficiency = (coded_logdet / m as f64 - (n as f64).ln()).exp();
    Ok(format!("{{\"version\":1,\"operation\":\"design_evaluate\",\"rank_method\":\"twice_reorthogonalized_qr_on_centered_rms_scaled_terms\",\"n_runs\":{n},\"model_terms_excluding_intercept\":{p},\"parameter_count\":{m},\"rank\":{m},\"information_basis\":\"supplied_model_terms_with_intercept_original_units\",\"information_matrix\":{},\"inverse_information_matrix\":{},\"log_determinant\":{raw_logdet},\"coding\":\"intercept_then_origin_relative_mean_centered_population_rms_scaled_terms\",\"feature_origins\":{},\"mean_offsets\":{},\"feature_means\":{},\"feature_scales\":{},\"coded_rows\":{},\"coded_information_matrix\":{},\"coded_inverse_information_matrix\":{},\"coded_log_determinant\":{coded_logdet},\"condition_number_coded\":{condition},\"condition_norm\":\"spectral_2_norm_of_coded_design_matrix\",\"d_efficiency_coded\":{efficiency},\"leverages\":{}}}",matrix(&raw_gram)?,matrix(&raw_inverse)?,vector(&origins)?,vector(&offsets)?,vector(&means)?,vector(&scales)?,matrix(&coded_rows)?,matrix(&coded_gram)?,matrix(&coded_inverse)?,vector(&leverages)?))
}
