//! Deterministic experimental designs and full quadratic response surfaces.
//!
//! `factorial_design` consumes one `[low, high]` row per factor and optional
//! `[center_points=0]`. `central_composite_design` uses the same bounds and
//! `[alpha=2^(p/4), center_points=1]`. Axial points with alpha > 1 extend outside
//! the specified factorial bounds; these bounds are not hard constraints.
//!
//! `d_optimal_design` selects distinct rows of a candidate feature matrix using
//! `[n_runs, max_iterations=50, seed=0]`. An intercept is added automatically.
//! Its determinant uses mean-centered, half-range-scaled features; only supplied candidates can
//! be selected. Greedy initialization plus exchange finds a local optimum,
//! without claiming a globally optimal design. Indices are zero based.
//!
//! `response_surface` consumes `[features..., target]`, `[lambda=0]`, and fits
//! `b0 + sum(b_i*z_i) + sum(b_ii*z_i^2) + sum(b_ij*z_i*z_j)`, where
//! `z_i=(x_i-means_i)/scales_i`. Terms are ordered intercept, linear, squares,
//! then interactions (i < j, lexicographic). Ridge penalizes all non-intercept
//! coefficients in this normalized basis. `model_parameters` is the portable
//! model `[p, means..., scales..., coefficients...]` accepted unchanged by
//! `response_surface_predict` on raw feature rows. No train/test separation or
//! predictive uncertainty is implied by the in-sample fit statistics.

const MAX_CELLS: usize = 1_000_000;
const MAX_FACTORS: usize = 16;

pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    validate_matrix(data)?;
    if parameters.iter().any(|v| !v.is_finite()) {
        return Err("parameters must be finite".into());
    }
    match operation {
        "factorial_design" => factorial(data, parameters, false),
        "central_composite_design" => factorial(data, parameters, true),
        "d_optimal_design" => d_optimal(data, parameters),
        "response_surface" => response_surface(data, parameters),
        "response_surface_predict" => predict(data, parameters),
        _ => Err(format!("unknown design operation: {operation}")),
    }
}

fn validate_matrix(data: &[Vec<f64>]) -> Result<(), String> {
    let columns = data.first().map_or(0, Vec::len);
    if data.is_empty() || columns == 0 || data.len() > 65_536 || columns > 32 {
        return Err("expected a nonempty matrix with at most 65536 rows and 32 columns".into());
    }
    if data.len().saturating_mul(columns) > MAX_CELLS {
        return Err("matrix exceeds the 1000000-cell limit".into());
    }
    if data
        .iter()
        .any(|r| r.len() != columns || r.iter().any(|x| !x.is_finite()))
    {
        return Err("matrix must be rectangular and finite".into());
    }
    Ok(())
}

fn integer(value: f64, name: &str, low: usize, high: usize) -> Result<usize, String> {
    if value.fract() != 0. || value < low as f64 || value > high as f64 {
        Err(format!("{name} must be an integer in [{low},{high}]"))
    } else {
        Ok(value as usize)
    }
}

fn vector(values: &[f64]) -> Result<String, String> {
    if values.iter().any(|v| !v.is_finite()) {
        return Err("numerical overflow in result".into());
    }
    Ok(format!(
        "[{}]",
        values
            .iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(",")
    ))
}

fn matrix(rows: &[Vec<f64>]) -> Result<String, String> {
    Ok(format!(
        "[{}]",
        rows.iter()
            .map(|r| vector(r))
            .collect::<Result<Vec<_>, _>>()?
            .join(",")
    ))
}

fn factorial(bounds: &[Vec<f64>], parameters: &[f64], composite: bool) -> Result<String, String> {
    let p = bounds.len();
    if p > MAX_FACTORS || bounds[0].len() != 2 || bounds.iter().any(|r| r[0] >= r[1]) {
        return Err("bounds require 1 to 16 [low,high] rows with low < high".into());
    }
    if parameters.len() > if composite { 2 } else { 1 } {
        return Err("too many design parameters".into());
    }
    let alpha = if composite {
        parameters
            .first()
            .copied()
            .unwrap_or(2_f64.powf(p as f64 / 4.))
    } else {
        1.
    };
    if alpha <= 0. {
        return Err("alpha must be positive".into());
    }
    let centers = integer(
        parameters
            .get(usize::from(composite))
            .copied()
            .unwrap_or(f64::from(composite)),
        "center_points",
        0,
        65_536,
    )?;
    let corners = 1_usize << p;
    let runs = corners + if composite { 2 * p } else { 0 } + centers;
    if runs > 65_536 || runs.saturating_mul(p) > MAX_CELLS {
        return Err("design exceeds the run or cell limit".into());
    }
    let mid: Vec<f64> = bounds.iter().map(|r| r[0] * 0.5 + r[1] * 0.5).collect();
    let half: Vec<f64> = bounds.iter().map(|r| r[1] * 0.5 - r[0] * 0.5).collect();
    if half.iter().any(|v| *v <= 0.) {
        return Err("bounds are too close to represent their midpoint".into());
    }
    let mut coded: Vec<Vec<f64>> = (0..corners)
        .map(|mask| {
            (0..p)
                .map(|j| if mask & (1 << j) == 0 { -1. } else { 1. })
                .collect()
        })
        .collect();
    if composite {
        for j in 0..p {
            for sign in [-1., 1.] {
                let mut row = vec![0.; p];
                row[j] = sign * alpha;
                coded.push(row);
            }
        }
    }
    coded.extend((0..centers).map(|_| vec![0.; p]));
    let rows: Vec<Vec<f64>> = coded
        .iter()
        .map(|r| {
            r.iter()
                .enumerate()
                .map(|(j, z)| {
                    if *z == -1. {
                        bounds[j][0]
                    } else if *z == 1. {
                        bounds[j][1]
                    } else {
                        mid[j] + half[j] * z
                    }
                })
                .collect()
        })
        .collect();
    Ok(format!("{{\"version\":1,\"operation\":\"{}\",\"factors\":{p},\"run_count\":{runs},\"center_points\":{centers},\"alpha\":{alpha},\"rows\":{},\"coded_rows\":{}}}",if composite { "central_composite_design" } else { "factorial_design" },matrix(&rows)?,matrix(&coded)?))
}

fn norm(values: impl IntoIterator<Item = f64>) -> f64 {
    values.into_iter().fold(0., f64::hypot)
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}

// Subtraction may overflow for finite values of opposite signs. Halving first
// preserves their normalized ratio without overflowing that subtraction.
fn scaled(x: f64, center: f64, scale: f64) -> f64 {
    let difference = x - center;
    if difference.is_finite() {
        difference / scale
    } else {
        (x * 0.5 - center * 0.5) / (scale * 0.5)
    }
}

fn normalization(data: &[Vec<f64>], p: usize) -> (Vec<f64>, Vec<f64>) {
    let mut means = Vec::with_capacity(p);
    let mut scales = Vec::with_capacity(p);
    for j in 0..p {
        let lo = data.iter().map(|r| r[j]).fold(f64::INFINITY, f64::min);
        let hi = data.iter().map(|r| r[j]).fold(f64::NEG_INFINITY, f64::max);
        let midpoint = lo * 0.5 + hi * 0.5;
        let half = hi * 0.5 - lo * 0.5;
        if half == 0. {
            means.push(lo);
            scales.push(1.);
            continue;
        }
        let normalized_mean: f64 = data
            .iter()
            .map(|r| scaled(r[j], midpoint, half) / data.len() as f64)
            .sum();
        let mean = midpoint + half * normalized_mean;
        means.push(mean.clamp(lo, hi));
        scales.push(half);
    }
    (means, scales)
}

fn cholesky_inverse(gram: &[Vec<f64>]) -> Result<(Vec<Vec<f64>>, f64), String> {
    let m = gram.len();
    let mut l = vec![vec![0.; m]; m];
    let threshold = gram
        .iter()
        .enumerate()
        .map(|(i, r)| r[i])
        .fold(0., f64::max)
        * 1e-12;
    let mut logdet = 0.;
    for i in 0..m {
        for j in 0..=i {
            let s = gram[i][j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if !s.is_finite() || s <= threshold {
                    return Err("candidate design is numerically rank deficient".into());
                }
                l[i][j] = s.sqrt();
                logdet += s.ln();
            } else {
                l[i][j] = s / l[j][j];
            }
        }
    }
    let mut inverse_columns = Vec::with_capacity(m);
    for k in 0..m {
        let mut y = vec![0.; m];
        for i in 0..m {
            y[i] = (if i == k { 1. } else { 0. } - (0..i).map(|j| l[i][j] * y[j]).sum::<f64>())
                / l[i][i];
        }
        let mut column = vec![0.; m];
        for i in (0..m).rev() {
            column[i] = (y[i] - (i + 1..m).map(|j| l[j][i] * column[j]).sum::<f64>()) / l[i][i];
        }
        inverse_columns.push(column);
    }
    let inverse = (0..m)
        .map(|i| inverse_columns.iter().map(|column| column[i]).collect())
        .collect();
    Ok((inverse, logdet))
}

fn information(features: &[Vec<f64>], selected: &[usize]) -> Vec<Vec<f64>> {
    let m = features[0].len();
    let mut gram = vec![vec![0.; m]; m];
    for &i in selected {
        for j in 0..m {
            for k in 0..m {
                gram[j][k] += features[i][j] * features[i][k];
            }
        }
    }
    gram
}

fn product(matrix: &[Vec<f64>], row: &[f64]) -> Vec<f64> {
    matrix.iter().map(|r| dot(r, row)).collect()
}

fn d_optimal(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let p = data[0].len();
    let m = p + 1;
    if parameters.is_empty() || parameters.len() > 3 || p > MAX_FACTORS || data.len() > 2048 {
        return Err("D-optimal design needs [n_runs, max_iterations=50, seed=0], <=2048 candidates and <=16 features".into());
    }
    let n = integer(parameters[0], "n_runs", m, 256.min(data.len()))?;
    let iterations = integer(
        parameters.get(1).copied().unwrap_or(50.),
        "max_iterations",
        0,
        1000,
    )?;
    let seed_value = parameters.get(2).copied().unwrap_or(0.);
    if seed_value.fract() != 0. || !(0. ..=9_007_199_254_740_991.).contains(&seed_value) {
        return Err("seed must be an integer in [0,9007199254740991]".into());
    }
    let seed = seed_value as u64;
    if data
        .len()
        .saturating_mul(n)
        .saturating_mul(m * m)
        .saturating_mul(iterations.max(1))
        > 500_000_000
    {
        return Err("D-optimal search exceeds the 500000000-operation work budget; reduce iterations, candidates, factors or runs".into());
    }
    let mut unique = Vec::new();
    for i in 0..data.len() {
        if !unique.iter().any(|&j| data[j] == data[i]) {
            unique.push(i);
        }
    }
    if unique.len() < n {
        return Err("not enough distinct candidate rows".into());
    }
    let (means, scales) = normalization(data, p);
    let features: Vec<Vec<f64>> = unique
        .iter()
        .map(|&i| {
            let mut r = vec![1.];
            r.extend((0..p).map(|j| scaled(data[i][j], means[j], scales[j])));
            r
        })
        .collect();
    let mut order: Vec<usize> = (0..unique.len()).collect();
    let mut state = seed ^ 0x9e3779b97f4a7c15;
    for i in (1..order.len()).rev() {
        state ^= state << 13;
        state ^= state >> 7;
        state ^= state << 17;
        order.swap(i, (state % (i as u64 + 1)) as usize);
    }
    let mut selected = Vec::with_capacity(n);
    let mut basis: Vec<Vec<f64>> = Vec::with_capacity(m);
    for _ in 0..m {
        let mut best = None;
        let mut best_norm = 0.;
        let mut best_residual = Vec::new();
        for &i in &order {
            if selected.contains(&i) {
                continue;
            }
            let mut residual = features[i].clone();
            for _ in 0..2 {
                for q in &basis {
                    let projection = dot(&residual, q);
                    for j in 0..m {
                        residual[j] -= projection * q[j];
                    }
                }
            }
            let length = norm(residual.iter().copied());
            if length > best_norm {
                best = Some(i);
                best_norm = length;
                best_residual = residual;
            }
        }
        if best_norm <= 1e-10 {
            return Err("candidate design is numerically rank deficient".into());
        }
        selected.push(best.ok_or("candidate design is rank deficient")?);
        basis.push(best_residual.into_iter().map(|v| v / best_norm).collect());
    }
    while selected.len() < n {
        let (inverse, _) = cholesky_inverse(&information(&features, &selected))?;
        let best = order
            .iter()
            .copied()
            .filter(|i| !selected.contains(i))
            .max_by(|&a, &b| {
                dot(&features[a], &product(&inverse, &features[a]))
                    .total_cmp(&dot(&features[b], &product(&inverse, &features[b])))
            })
            .ok_or("not enough candidates")?;
        selected.push(best);
    }
    let mut exchanges = 0;
    let mut converged = false;
    for _ in 0..iterations {
        let (inverse, _) = cholesky_inverse(&information(&features, &selected))?;
        let mut best_ratio = 1. + 1e-10;
        let mut best = None;
        for &incoming in &order {
            if selected.contains(&incoming) {
                continue;
            }
            let inv_new = product(&inverse, &features[incoming]);
            let h_new = dot(&features[incoming], &inv_new);
            for (slot, &outgoing) in selected.iter().enumerate() {
                let h_old = dot(&features[outgoing], &product(&inverse, &features[outgoing]));
                let cross = dot(&features[outgoing], &inv_new);
                let ratio = (1. - h_old) * (1. + h_new) + cross * cross;
                if ratio > best_ratio {
                    best_ratio = ratio;
                    best = Some((slot, incoming));
                }
            }
        }
        if let Some((slot, incoming)) = best {
            selected[slot] = incoming;
            exchanges += 1;
        } else {
            converged = true;
            break;
        }
    }
    let (_, logdet) = cholesky_inverse(&information(&features, &selected))?;
    let mut indices: Vec<usize> = selected.iter().map(|&i| unique[i]).collect();
    indices.sort_unstable();
    let rows: Vec<Vec<f64>> = indices.iter().map(|&i| data[i].clone()).collect();
    Ok(format!("{{\"version\":1,\"operation\":\"d_optimal_design\",\"indices\":{:?},\"rows\":{},\"log_determinant\":{logdet},\"intercept\":true,\"means\":{},\"scales\":{},\"exchanges\":{exchanges},\"converged\":{converged},\"optimality\":\"local_exchange\"}}",indices,matrix(&rows)?,vector(&means)?,vector(&scales)?))
}

fn terms(p: usize) -> usize {
    1 + 2 * p + p * (p - 1) / 2
}

fn quadratic(row: &[f64], means: &[f64], scales: &[f64]) -> Vec<f64> {
    let z: Vec<f64> = row
        .iter()
        .zip(means)
        .zip(scales)
        .map(|((&x, &mean), &scale)| scaled(x, mean, scale))
        .collect();
    let mut x = Vec::with_capacity(terms(z.len()));
    x.push(1.);
    x.extend(z.iter().copied());
    x.extend(z.iter().map(|v| v * v));
    for i in 0..z.len() {
        for j in i + 1..z.len() {
            x.push(z[i] * z[j]);
        }
    }
    x
}

// Column-equilibrated, column-pivoted Householder QR avoids normal equations
// and allows ridge rows to identify otherwise underdetermined models.
fn least_squares(mut a: Vec<Vec<f64>>, mut y: Vec<f64>) -> Result<Vec<f64>, String> {
    let n = a.len();
    let m = a[0].len();
    if n < m {
        return Err("quadratic model is underdetermined; add runs or a positive lambda".into());
    }
    let mut lengths = Vec::with_capacity(m);
    for j in 0..m {
        let length = norm(a.iter().map(|r| r[j]));
        if !length.is_finite() || length == 0. {
            return Err(
                "quadratic model is rank deficient; add independent runs or a positive lambda"
                    .into(),
            );
        }
        lengths.push(length);
        for row in &mut a {
            row[j] /= length;
        }
    }
    let mut permutation: Vec<usize> = (0..m).collect();
    for k in 0..m {
        let pivot = (k..m)
            .max_by(|&i, &j| {
                norm(a[k..].iter().map(|r| r[i])).total_cmp(&norm(a[k..].iter().map(|r| r[j])))
            })
            .unwrap_or(k);
        for row in &mut a {
            row.swap(k, pivot);
        }
        permutation.swap(k, pivot);
        let length = norm(a[k..].iter().map(|r| r[k]));
        if length <= 1e-10 {
            return Err("quadratic model is numerically rank deficient; add independent runs or a positive lambda".into());
        }
        let alpha = if a[k][k] >= 0. { -length } else { length };
        let mut v: Vec<f64> = a[k..].iter().map(|r| r[k]).collect();
        v[0] -= alpha;
        let vnorm = norm(v.iter().copied());
        for value in &mut v {
            *value /= vnorm;
        }
        let projections: Vec<f64> = (k..m)
            .map(|j| 2. * (k..n).map(|i| v[i - k] * a[i][j]).sum::<f64>())
            .collect();
        for i in k..n {
            for (entry, projection) in a[i][k..].iter_mut().zip(&projections) {
                *entry -= projection * v[i - k];
            }
        }
        let projection = 2. * dot(&v, &y[k..]);
        for i in k..n {
            y[i] -= projection * v[i - k];
        }
        a[k][k] = alpha;
    }
    let mut coefficients = vec![0.; m];
    for i in (0..m).rev() {
        coefficients[i] =
            (y[i] - (i + 1..m).map(|j| a[i][j] * coefficients[j]).sum::<f64>()) / a[i][i];
    }
    let mut ordered = vec![0.; m];
    for j in 0..m {
        ordered[permutation[j]] = coefficients[j] / lengths[permutation[j]];
    }
    Ok(ordered)
}

fn response_surface(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let p = data[0].len() - 1;
    if p == 0 || p > MAX_FACTORS || parameters.len() > 1 {
        return Err("response_surface needs 1 to 16 features plus target and [lambda=0]".into());
    }
    let lambda = parameters.first().copied().unwrap_or(0.);
    if lambda < 0. {
        return Err("lambda must be nonnegative".into());
    }
    let m = terms(p);
    if data.len().saturating_mul(m * m) > 100_000_000
        || (data.len() + m).saturating_mul(m) > MAX_CELLS
    {
        return Err("quadratic fit exceeds the cell or work limit; reduce rows or features".into());
    }
    let (means, scales) = normalization(data, p);
    let features: Vec<Vec<f64>> = data
        .iter()
        .map(|r| quadratic(&r[..p], &means, &scales))
        .collect();
    let mut a = features.clone();
    let targets: Vec<Vec<f64>> = data.iter().map(|r| vec![r[p]]).collect();
    let (target_means, target_scales) = normalization(&targets, 1);
    let ymean = target_means[0];
    let yscale = target_scales[0];
    let mut y: Vec<f64> = data.iter().map(|r| scaled(r[p], ymean, yscale)).collect();
    if lambda > 0. {
        for j in 1..m {
            let mut row = vec![0.; m];
            row[j] = lambda.sqrt();
            a.push(row);
            y.push(0.);
        }
    }
    let mut coefficients = least_squares(a, y)?;
    for coefficient in &mut coefficients {
        *coefficient *= yscale;
    }
    coefficients[0] += ymean;
    let fitted: Vec<f64> = features.iter().map(|r| dot(r, &coefficients)).collect();
    let rmse = norm(data.iter().zip(&fitted).map(|(r, &f)| r[p] - f)) / (data.len() as f64).sqrt();
    let residual_scaled = norm(
        data.iter()
            .zip(&fitted)
            .map(|(r, &f)| scaled(r[p], f, yscale)),
    );
    let target_scaled = norm(data.iter().map(|r| scaled(r[p], ymean, yscale)));
    let r2 = if target_scaled == 0. {
        if residual_scaled == 0. {
            1.
        } else {
            0.
        }
    } else {
        1. - (residual_scaled / target_scaled).powi(2)
    };
    if !rmse.is_finite() || !r2.is_finite() {
        return Err("fit statistics overflowed".into());
    }
    let mut model = vec![p as f64];
    model.extend_from_slice(&means);
    model.extend_from_slice(&scales);
    model.extend_from_slice(&coefficients);
    let mut names = vec!["\"intercept\"".to_string()];
    names.extend((0..p).map(|j| format!("\"z{j}\"")));
    names.extend((0..p).map(|j| format!("\"z{j}^2\"")));
    for i in 0..p {
        for j in i + 1..p {
            names.push(format!("\"z{i}*z{j}\""));
        }
    }
    Ok(format!("{{\"version\":1,\"operation\":\"response_surface\",\"feature_count\":{p},\"basis\":\"z=(x-means)/scales\",\"terms\":[{}],\"means\":{},\"scales\":{},\"coefficients\":{},\"model_parameters\":{},\"fitted\":{},\"rmse\":{rmse},\"r2\":{r2},\"lambda\":{lambda},\"statistics\":\"in_sample\"}}",names.join(","),vector(&means)?,vector(&scales)?,vector(&coefficients)?,vector(&model)?,vector(&fitted)?))
}

fn predict(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let p = integer(
        *parameters.first().ok_or("missing response surface model")?,
        "model feature count",
        1,
        MAX_FACTORS,
    )?;
    if parameters.len() != 1 + 2 * p + terms(p) || data[0].len() != p {
        return Err("model or prediction feature dimensions do not match".into());
    }
    let means = &parameters[1..1 + p];
    let scales = &parameters[1 + p..1 + 2 * p];
    let coefficients = &parameters[1 + 2 * p..];
    if scales.iter().any(|s| *s <= 0.) {
        return Err("model scales must be positive".into());
    }
    let predictions: Vec<f64> = data
        .iter()
        .map(|r| dot(&quadratic(r, means, scales), coefficients))
        .collect();
    Ok(format!(
        "{{\"version\":1,\"operation\":\"response_surface_predict\",\"predictions\":{}}}",
        vector(&predictions)?
    ))
}
