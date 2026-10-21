//! Deterministic, dependency-free exploration of descriptor matrices.
//!
//! Rows are observations and columns are descriptors. PCA centers every column
//! and uses sample covariance (`n - 1`); optional scaling uses sample standard
//! deviations. K-means uses squared Euclidean distance on the supplied columns,
//! so callers should scale descriptors beforehand when their units differ.

use std::fmt::Write;

const MAX_ROWS: usize = 100_000;
const MAX_COLUMNS: usize = 256;
const MAX_CELLS: usize = 1_000_000;
const MAX_WORK: usize = 500_000_000;
const MAX_ABSOLUTE_VALUE: f64 = 1e150;
const MAX_SEED: f64 = 9_007_199_254_740_991.0;

/// Run `pca` or `kmeans` and return a version-1 JSON result.
///
/// PCA parameters are `[components, standardize]`, defaulting to
/// `[min(columns, rows - 1), 0]`. `standardize` must be 0 or 1. K-means
/// parameters are `[clusters, max_iterations, seed]`, defaulting to
/// `[min(2, rows), 300, 0]`. Integer parameters must be exact nonnegative
/// integers; the seed is limited to the exactly represented JSON integer range.
///
/// Bounds are 100,000 rows, 256 columns and 1,000,000 cells. Excessive estimated
/// work, nonfinite inputs, and magnitudes above 1e150 return an error before
/// expensive work. K-means additionally caps clusters at 256 and iterations at
/// 10,000. These limits bound allocation, serialization, and numerical overflow.
/// PCA rejects offsets so large that the represented mean cannot center the
/// observed variation reliably; subtract a reference value beforehand.
pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    if operation != "pca" && operation != "kmeans" {
        return Err(format!("unknown exploration operation: {operation}"));
    }
    validate_data(data)?;
    match operation {
        "pca" => pca(data, parameters),
        "kmeans" => kmeans(data, parameters),
        _ => unreachable!(),
    }
}

fn validate_data(data: &[Vec<f64>]) -> Result<(), String> {
    let columns = data.first().map_or(0, Vec::len);
    if data.is_empty() || columns == 0 {
        return Err("data must be a nonempty matrix with at least one column".into());
    }
    if data.len() > MAX_ROWS || columns > MAX_COLUMNS || data.len() * columns > MAX_CELLS {
        return Err("data exceeds the row, column, or cell limit".into());
    }
    for row in data {
        if row.len() != columns {
            return Err("data must be rectangular".into());
        }
        if row
            .iter()
            .any(|value| !value.is_finite() || value.abs() > MAX_ABSOLUTE_VALUE)
        {
            return Err("data must be finite with absolute values at most 1e150".into());
        }
    }
    Ok(())
}

fn integer_parameter(value: f64, name: &str, min: usize, max: usize) -> Result<usize, String> {
    if !value.is_finite() || value.fract() != 0.0 || value < min as f64 || value > max as f64 {
        return Err(format!("{name} must be an integer in [{min}, {max}]"));
    }
    Ok(value as usize)
}

fn append_vector(output: &mut String, values: &[f64]) {
    output.push('[');
    for (index, value) in values.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        // Scientific notation bounds serialized size even for extreme scales.
        let _ = write!(output, "{value:.17e}");
    }
    output.push(']');
}

fn append_matrix(output: &mut String, matrix: &[Vec<f64>]) {
    output.push('[');
    for (index, row) in matrix.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        append_vector(output, row);
    }
    output.push(']');
}

fn pca(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let rows = data.len();
    let columns = data[0].len();
    if rows < 2 {
        return Err("PCA requires at least two observations".into());
    }
    if parameters.len() > 2 {
        return Err("PCA accepts [components, standardize]".into());
    }
    let limit = columns.min(rows - 1);
    let count = integer_parameter(
        parameters.first().copied().unwrap_or(limit as f64),
        "components",
        1,
        limit,
    )?;
    let standardize = integer_parameter(
        parameters.get(1).copied().unwrap_or(0.0),
        "standardize",
        0,
        1,
    )? == 1;
    if rows * columns * columns > MAX_WORK {
        return Err("PCA covariance computation exceeds the work limit".into());
    }

    let mut mean = vec![0.0; columns];
    let mut correction = vec![0.0; columns];
    for row in data {
        for column in 0..columns {
            let value = row[column];
            let sum = mean[column] + value;
            correction[column] += if mean[column].abs() >= value.abs() {
                (mean[column] - sum) + value
            } else {
                (value - sum) + mean[column]
            };
            mean[column] = sum;
        }
    }
    for (average, correction) in mean.iter_mut().zip(correction) {
        *average = (*average + correction) / rows as f64;
    }
    for column in 0..columns {
        if data.iter().all(|row| row[column] == data[0][column]) {
            mean[column] = data[0][column];
        }
    }
    let mut centered: Vec<Vec<f64>> = data
        .iter()
        .map(|row| {
            row.iter()
                .zip(&mean)
                .map(|(x, average)| x - average)
                .collect()
        })
        .collect();
    for column in 0..columns {
        let maximum = centered
            .iter()
            .map(|row| row[column].abs())
            .fold(0.0_f64, f64::max);
        if maximum != 0.0 {
            let mut residual = 0.0_f64;
            let mut correction = 0.0_f64;
            for row in &centered {
                let term = (row[column] / maximum) / rows as f64 - correction;
                let sum = residual + term;
                correction = (sum - residual) - term;
                residual = sum;
            }
            // A mean error below 1e-8 of the maximum deviation contributes at
            // most about machine precision to the covariance's squared scale.
            if residual.abs() > 1e-8 {
                return Err("PCA cannot reliably center this feature offset; subtract a reference value before analysis".into());
            }
        }
    }
    let mut scale = vec![1.0; columns];
    if standardize {
        for column in 0..columns {
            let maximum = centered
                .iter()
                .map(|row| row[column].abs())
                .fold(0.0_f64, f64::max);
            if maximum != 0.0 {
                let squares: f64 = centered
                    .iter()
                    .map(|row| (row[column] / maximum).powi(2))
                    .sum();
                scale[column] = maximum * (squares / (rows - 1) as f64).sqrt();
                for row in &mut centered {
                    row[column] /= scale[column];
                }
            }
        }
    }

    // Normalize covariance construction to preserve relative precision for
    // tiny descriptor units and avoid overflow for very large descriptor units.
    let magnitude = centered
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let divisor = if magnitude == 0.0 { 1.0 } else { magnitude };
    let mut covariance = vec![vec![0.0; columns]; columns];
    for row in &centered {
        for a in 0..columns {
            for b in a..columns {
                covariance[a][b] += (row[a] / divisor) * (row[b] / divisor);
            }
        }
    }
    for row in &mut covariance {
        for value in row {
            *value /= (rows - 1) as f64;
        }
    }
    for a in 0..columns {
        let (preceding, following) = covariance.split_at_mut(a);
        for (b, row) in preceding.iter().enumerate() {
            following[0][b] = row[a];
        }
    }
    let (eigenvalues, vectors) = symmetric_eigen(covariance)?;
    let mut order: Vec<usize> = (0..columns).collect();
    order.sort_by(|&a, &b| eigenvalues[b].total_cmp(&eigenvalues[a]).then(a.cmp(&b)));
    let total: f64 = eigenvalues.iter().sum();
    let mut components = Vec::with_capacity(count);
    let mut variance = Vec::with_capacity(count);
    let mut ratios = Vec::with_capacity(count);
    for &index in order.iter().take(count) {
        let mut component: Vec<f64> = vectors.iter().map(|row| row[index]).collect();
        // Eigenvector signs have no statistical meaning; fix the largest
        // loading positive so repeated results are easy to compare.
        let pivot = component
            .iter()
            .enumerate()
            .max_by(|(i, a), (j, b)| a.abs().total_cmp(&b.abs()).then(j.cmp(i)))
            .map_or(0, |(index, _)| index);
        if component[pivot] < 0.0 {
            for value in &mut component {
                *value = -*value;
            }
        }
        components.push(component);
        variance.push((eigenvalues[index] * divisor) * divisor);
        ratios.push(if total == 0.0 {
            0.0
        } else {
            eigenvalues[index] / total
        });
    }
    let scores: Vec<Vec<f64>> = centered
        .iter()
        .map(|row| {
            components
                .iter()
                .map(|component| {
                    row.iter()
                        .zip(component)
                        .map(|(x, loading)| x * loading)
                        .sum()
                })
                .collect()
        })
        .collect();
    if scale
        .iter()
        .chain(&variance)
        .chain(&ratios)
        .chain(scores.iter().flatten())
        .any(|value| !value.is_finite())
    {
        return Err("PCA exceeded the finite numerical range".into());
    }
    let mut output = String::from("{\"version\":1,\"mean\":");
    append_vector(&mut output, &mean);
    output.push_str(",\"scale\":");
    append_vector(&mut output, &scale);
    output.push_str(",\"components\":");
    append_matrix(&mut output, &components);
    output.push_str(",\"scores\":");
    append_matrix(&mut output, &scores);
    output.push_str(",\"explained_variance\":");
    append_vector(&mut output, &variance);
    output.push_str(",\"explained_variance_ratio\":");
    append_vector(&mut output, &ratios);
    output.push('}');
    Ok(output)
}

/// Cyclic Jacobi rotations for a real symmetric matrix. This returns vectors
/// as columns and uses a relative residual to distinguish convergence from a
/// merely exhausted iteration budget.
fn symmetric_eigen(mut matrix: Vec<Vec<f64>>) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let dimension = matrix.len();
    let mut vectors = vec![vec![0.0; dimension]; dimension];
    for (index, row) in vectors.iter_mut().enumerate() {
        row[index] = 1.0;
    }
    let norm = matrix
        .iter()
        .flatten()
        .map(|value| value.abs())
        .fold(0.0_f64, f64::max);
    let tolerance = norm * f64::EPSILON * 16.0;
    let mut converged = norm == 0.0;
    for _ in 0..100 {
        let mut largest = 0.0_f64;
        for p in 0..dimension {
            for q in p + 1..dimension {
                largest = largest.max(matrix[p][q].abs());
                if matrix[p][q].abs() <= tolerance {
                    continue;
                }
                let off = matrix[p][q];
                let theta = (matrix[q][q] - matrix[p][p]) / (2.0 * off);
                let t = 1.0_f64.copysign(theta) / (theta.abs() + theta.hypot(1.0));
                let cosine = 1.0 / (1.0 + t * t).sqrt();
                let sine = t * cosine;
                matrix[p][p] -= t * off;
                matrix[q][q] += t * off;
                matrix[p][q] = 0.0;
                matrix[q][p] = 0.0;
                for r in 0..dimension {
                    if r != p && r != q {
                        let a = matrix[r][p];
                        let b = matrix[r][q];
                        matrix[r][p] = cosine * a - sine * b;
                        matrix[p][r] = matrix[r][p];
                        matrix[r][q] = sine * a + cosine * b;
                        matrix[q][r] = matrix[r][q];
                    }
                    let a = vectors[r][p];
                    let b = vectors[r][q];
                    vectors[r][p] = cosine * a - sine * b;
                    vectors[r][q] = sine * a + cosine * b;
                }
            }
        }
        if largest <= tolerance {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("PCA symmetric eigensolver did not converge".into());
    }
    let mut eigenvalues = Vec::with_capacity(dimension);
    for (index, row) in matrix.iter().enumerate() {
        if row[index] < -tolerance * dimension as f64 {
            return Err("PCA covariance has a materially negative eigenvalue".into());
        }
        eigenvalues.push(row[index].max(0.0));
    }
    Ok((eigenvalues, vectors))
}

fn squared_distance(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| (x - y).powi(2)).sum()
}

// SplitMix64 has a fully specified sequence and needs no platform RNG or seed
// state. Only initialization uses randomness; ties always use stable row order.
fn random(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut value = *state;
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

fn kmeans(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    if parameters.len() > 3 {
        return Err("k-means accepts [clusters, max_iterations, seed]".into());
    }
    let rows = data.len();
    let columns = data[0].len();
    let clusters = integer_parameter(
        parameters.first().copied().unwrap_or(rows.min(2) as f64),
        "clusters",
        1,
        rows.min(256),
    )?;
    let max_iterations = integer_parameter(
        parameters.get(1).copied().unwrap_or(300.0),
        "max_iterations",
        1,
        10_000,
    )?;
    let seed = parameters.get(2).copied().unwrap_or(0.0);
    if !seed.is_finite() || seed.fract() != 0.0 || !(0.0..=MAX_SEED).contains(&seed) {
        return Err("seed must be an exact integer in [0, 9007199254740991]".into());
    }
    let work = rows
        .checked_mul(columns)
        .and_then(|n| n.checked_mul(clusters))
        .and_then(|n| n.checked_mul(max_iterations));
    if work.is_none_or(|work| work > MAX_WORK) {
        return Err("k-means exceeds the estimated work limit".into());
    }

    // K-means++ sampling spreads initial centers while preserving a specified
    // deterministic seed. Scale weights before summing to avoid overflow.
    let mut state = seed as u64;
    let first = (random(&mut state) % rows as u64) as usize;
    let mut selected = vec![false; rows];
    selected[first] = true;
    let mut centers = vec![data[first].clone()];
    let mut closest = vec![f64::INFINITY; rows];
    while centers.len() < clusters {
        let last = centers.last().ok_or("k-means initialization failed")?;
        for (index, row) in data.iter().enumerate() {
            closest[index] = closest[index].min(squared_distance(row, last));
        }
        let maximum = closest.iter().copied().fold(0.0_f64, f64::max);
        let next = if maximum == 0.0 {
            selected
                .iter()
                .position(|&used| !used)
                .ok_or("k-means initialization failed")?
        } else {
            let total: f64 = closest.iter().map(|distance| distance / maximum).sum();
            let uniform = (random(&mut state) >> 11) as f64 / 9_007_199_254_740_992.0;
            let threshold = uniform * total;
            let mut cumulative = 0.0;
            let mut candidate = None;
            for index in 0..rows {
                cumulative += closest[index] / maximum;
                if !selected[index] && closest[index] > 0.0 && cumulative > threshold {
                    candidate = Some(index);
                    break;
                }
            }
            candidate
                .or_else(|| {
                    (0..rows)
                        .rev()
                        .find(|&index| !selected[index] && closest[index] > 0.0)
                })
                .ok_or("k-means initialization failed")?
        };
        selected[next] = true;
        centers.push(data[next].clone());
    }

    let mut labels = vec![usize::MAX; rows];
    let mut iterations = 0;
    let mut converged = false;
    for iteration in 1..=max_iterations {
        let previous = labels.clone();
        let mut counts = vec![0_usize; clusters];
        for (index, row) in data.iter().enumerate() {
            // Retain the previous cluster on an exact tie. This avoids
            // oscillation after duplicate points repair an empty cluster.
            let mut best = if labels[index] < clusters {
                labels[index]
            } else {
                0
            };
            let mut distance = squared_distance(row, &centers[best]);
            for (cluster, center) in centers.iter().enumerate() {
                let candidate = squared_distance(row, center);
                if candidate < distance {
                    distance = candidate;
                    best = cluster;
                }
            }
            labels[index] = best;
            counts[best] += 1;
        }
        // Give an empty cluster the worst represented point from a cluster
        // containing at least two points. Every returned cluster is nonempty,
        // including the valid case with fewer distinct points than clusters.
        for empty in 0..clusters {
            if counts[empty] != 0 {
                continue;
            }
            let mut donor = None;
            let mut farthest = -1.0;
            for (index, row) in data.iter().enumerate() {
                if counts[labels[index]] > 1 {
                    let distance = squared_distance(row, &centers[labels[index]]);
                    if distance > farthest {
                        farthest = distance;
                        donor = Some(index);
                    }
                }
            }
            let donor = donor.ok_or("k-means could not repair an empty cluster")?;
            counts[labels[donor]] -= 1;
            labels[donor] = empty;
            counts[empty] += 1;
        }
        centers = vec![vec![0.0; columns]; clusters];
        let mut corrections = vec![vec![0.0; columns]; clusters];
        let mut origins = vec![usize::MAX; clusters];
        for (index, row) in data.iter().enumerate() {
            let cluster = labels[index];
            if origins[cluster] == usize::MAX {
                origins[cluster] = index;
            }
            for column in 0..columns {
                let offset = row[column] - data[origins[cluster]][column];
                let old = centers[cluster][column];
                let sum = old + offset;
                corrections[cluster][column] += if old.abs() >= offset.abs() {
                    (old - sum) + offset
                } else {
                    (offset - sum) + old
                };
                centers[cluster][column] = sum;
            }
        }
        for cluster in 0..clusters {
            for column in 0..columns {
                centers[cluster][column] = data[origins[cluster]][column]
                    + (centers[cluster][column] + corrections[cluster][column])
                        / counts[cluster] as f64;
            }
        }
        iterations = iteration;
        if labels == previous {
            converged = true;
            break;
        }
    }
    // Even when the iteration limit is reached, returned centers are exactly
    // the means of returned labels, and inertia describes that same partition.
    let inertia: f64 = data
        .iter()
        .zip(&labels)
        .map(|(row, &label)| squared_distance(row, &centers[label]))
        .sum();
    if !inertia.is_finite() || centers.iter().flatten().any(|value| !value.is_finite()) {
        return Err("k-means exceeded the finite numerical range".into());
    }
    let mut output = String::from("{\"version\":1,\"labels\":[");
    for (index, label) in labels.iter().enumerate() {
        if index != 0 {
            output.push(',');
        }
        let _ = write!(output, "{label}");
    }
    output.push_str("],\"centers\":");
    append_matrix(&mut output, &centers);
    let _ = write!(
        output,
        ",\"inertia\":{inertia:.17e},\"iterations\":{iterations},\"converged\":{converged}}}"
    );
    Ok(output)
}
