//! Regression, validation, and deterministic group partitioning.
//!
//! Ridge minimizes sum((y - intercept - X_scaled * beta)^2) + lambda * ||beta||^2,
//! using training-only population standard deviations and an unpenalized intercept.
//! An augmented Householder QR solve avoids squaring the design condition number.
//! Zero-penalty rank-deficient designs are rejected; PLS1 stops at predictive rank.

use std::collections::BTreeMap;

const MAX_CELLS: usize = 1_000_000;
const MAX_FEATURES: usize = 256;
const MAX_WORK: usize = 100_000_000;
const MAX_EXACT_INTEGER: f64 = 9_007_199_254_740_991.;

#[derive(Debug)]
struct Model {
    intercept: f64,
    coefficients: Vec<f64>,
    means: Vec<f64>,
    scales: Vec<f64>,
    origins: Vec<f64>,
    mean_offsets: Vec<f64>,
    response_origin: f64,
    response_mean_offset: f64,
    components_used: usize,
}

/// Runs a numeric operation and returns a versioned JSON result.
///
/// Fit rows are `[features..., target]`; validation rows append a nonnegative fold
/// id. Prediction parameters are `[intercept, coefficients...]`, or the returned
/// `stable_predict_parameters` vector for offset-safe model portability. Group
/// ids and seeds must be exact nonnegative integers representable in a double.
pub fn run(operation: &str, data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    let width = validate_matrix(data)?;
    if parameters.iter().any(|x| !x.is_finite()) {
        return Err("parameters must be finite".into());
    }
    match operation {
        "ridge" | "pls" => {
            if parameters.len() != 1 {
                return Err("fit requires one hyperparameter".into());
            }
            let model = fit(operation, data, parameters[0])?;
            let predictions = data
                .iter()
                .map(|row| predict(&model, &row[..width - 1]))
                .collect::<Result<Vec<_>, _>>()?;
            let actual = data.iter().map(|row| row[width - 1]).collect::<Vec<_>>();
            let (rmse, r2) = metrics(&actual, &predictions)?;
            let mut stable_parameters = vec![model.response_origin, model.response_mean_offset];
            stable_parameters.extend_from_slice(&model.origins);
            stable_parameters.extend_from_slice(&model.mean_offsets);
            stable_parameters.extend_from_slice(&model.coefficients);
            Ok(format!("{{\"version\":1,\"algorithm\":\"{operation}\",\"n_samples\":{},\"n_features\":{},\"intercept\":{},\"coefficients\":{},\"feature_means\":{},\"feature_scales\":{},\"stable_predict_parameters\":{},\"standardization\":\"training_population_standard_deviation\",\"hyperparameter\":{},\"components_used\":{},\"metrics\":{{\"rmse\":{rmse},\"r2\":{r2}}}}}", data.len(), width - 1, model.intercept, numbers(&model.coefficients), numbers(&model.means), numbers(&model.scales), numbers(&stable_parameters), parameters[0], model.components_used))
        }
        "predict" => {
            let stable = parameters.len() == 3 * width + 2;
            if parameters.len() != width + 1 && !stable {
                return Err(
                    "prediction requires [intercept, coefficients...] or stable_predict_parameters"
                        .into(),
                );
            }
            let model = Model {
                intercept: if stable { 0. } else { parameters[0] },
                coefficients: if stable {
                    parameters[2 + 2 * width..].to_vec()
                } else {
                    parameters[1..].to_vec()
                },
                means: Vec::new(),
                scales: Vec::new(),
                origins: if stable {
                    parameters[2..2 + width].to_vec()
                } else {
                    Vec::new()
                },
                mean_offsets: if stable {
                    parameters[2 + width..2 + 2 * width].to_vec()
                } else {
                    Vec::new()
                },
                response_origin: if stable { parameters[0] } else { 0. },
                response_mean_offset: if stable { parameters[1] } else { 0. },
                components_used: 0,
            };
            let predictions = data
                .iter()
                .map(|row| predict(&model, row))
                .collect::<Result<Vec<_>, _>>()?;
            Ok(format!(
                "{{\"version\":1,\"algorithm\":\"linear_predict\",\"predictions\":{}}}",
                numbers(&predictions)
            ))
        }
        "cross_validate" => cross_validate(data, parameters),
        "group_folds" => group_folds(data, parameters),
        _ => Err(format!("unknown regression operation: {operation}")),
    }
}

fn validate_matrix(data: &[Vec<f64>]) -> Result<usize, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0 || data.len() > MAX_CELLS / width {
        return Err("data must be a nonempty matrix with at most 1000000 cells".into());
    }
    if data
        .iter()
        .any(|row| row.len() != width || row.iter().any(|x| !x.is_finite()))
    {
        return Err("data must be rectangular and finite".into());
    }
    Ok(width)
}

fn integer(value: f64, name: &str) -> Result<u64, String> {
    if !value.is_finite() || !(0. ..=MAX_EXACT_INTEGER).contains(&value) || value.fract() != 0. {
        return Err(format!("{name} must be an exact nonnegative integer"));
    }
    Ok(value as u64)
}

fn mean(values: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut result = 0.;
    for (i, value) in values.enumerate() {
        // Preserve identical inputs exactly and avoid overflowing a raw sum.
        let difference = value - result;
        result = if difference.is_finite() {
            result + difference / (i + 1) as f64
        } else {
            let weight = 1. / (i + 1) as f64;
            result * (1. - weight) + value * weight
        };
    }
    finite(result, "mean")
}

fn finite(value: f64, name: &str) -> Result<f64, String> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(format!("{name} exceeds finite numeric range"))
    }
}

fn norm(values: impl Iterator<Item = f64>) -> f64 {
    values.fold(0_f64, f64::hypot)
}

fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

fn fit(method: &str, data: &[Vec<f64>], hyperparameter: f64) -> Result<Model, String> {
    let width = validate_matrix(data)?;
    let n = data.len();
    if n < 2 || width < 2 || width - 1 > MAX_FEATURES {
        return Err("fit requires at least two rows and 1 to 256 features".into());
    }
    let p = width - 1;
    if n.saturating_mul(p).saturating_mul(p) > MAX_WORK {
        return Err("fit exceeds numerical work limit".into());
    }
    let components = if method == "pls" {
        let count = integer(hyperparameter, "PLS components")?;
        if count == 0 || count > p as u64 {
            return Err("PLS components must be between 1 and feature count".into());
        }
        count as usize
    } else {
        if !hyperparameter.is_finite() || hyperparameter < 0. {
            return Err("ridge lambda must be finite and nonnegative".into());
        }
        0
    };
    // Center small differences before adding any large origin back. A rounded
    // absolute mean can otherwise erase a significant part of the spread.
    let origins = data[0][..p].to_vec();
    let mean_offsets = (0..p)
        .map(|j| mean(data.iter().map(|row| row[j] - origins[j])))
        .collect::<Result<Vec<_>, _>>()?;
    let means = origins
        .iter()
        .zip(&mean_offsets)
        .map(|(origin, offset)| finite(origin + offset, "feature mean"))
        .collect::<Result<Vec<_>, _>>()?;
    let response_origin = data[0][p];
    let response_mean_offset = mean(data.iter().map(|row| row[p] - response_origin))?;
    let scales = (0..p)
        .map(|j| {
            let deviation = finite(
                norm(
                    data.iter()
                        .map(|row| (row[j] - origins[j]) - mean_offsets[j]),
                ) / (n as f64).sqrt(),
                "feature scale",
            )?;
            Ok(if deviation == 0. { 1. } else { deviation })
        })
        .collect::<Result<Vec<_>, String>>()?;
    let x = data
        .iter()
        .map(|row| {
            (0..p)
                .map(|j| {
                    finite(
                        ((row[j] - origins[j]) - mean_offsets[j]) / scales[j],
                        "scaled feature",
                    )
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let y = data
        .iter()
        .map(|row| {
            finite(
                (row[p] - response_origin) - response_mean_offset,
                "centered target",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let (scaled_coefficients, components_used) = if method == "pls" {
        pls(x, y, components)?
    } else {
        (ridge(x, y, hyperparameter)?, 0)
    };
    let coefficients = scaled_coefficients
        .iter()
        .zip(&scales)
        .map(|(coefficient, scale)| finite(coefficient / scale, "coefficient"))
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        (response_origin - dot(&origins, &coefficients))
            + (response_mean_offset - dot(&mean_offsets, &coefficients)),
        "intercept",
    )?;
    Ok(Model {
        intercept,
        coefficients,
        means,
        scales,
        origins,
        mean_offsets,
        response_origin,
        response_mean_offset,
        components_used,
    })
}

// Householder updates traverse columns across separate rows.
#[allow(clippy::needless_range_loop)]
fn ridge(mut x: Vec<Vec<f64>>, mut y: Vec<f64>, lambda: f64) -> Result<Vec<f64>, String> {
    let p = x[0].len();
    if lambda > 0. {
        let penalty = lambda.sqrt();
        for j in 0..p {
            let mut row = vec![0.; p];
            row[j] = penalty;
            x.push(row);
            y.push(0.);
        }
    }
    if x.len() < p {
        return Err("rank-deficient design; use positive ridge lambda or PLS".into());
    }
    let tolerance =
        f64::EPSILON * 64. * (x.len().max(p) as f64).sqrt() * norm(x.iter().flatten().copied());
    for j in 0..p {
        let column_norm = finite(norm(x[j..].iter().map(|row| row[j])), "QR norm")?;
        if column_norm <= tolerance {
            return Err(
                "rank-deficient or numerically singular design; increase ridge lambda or use PLS"
                    .into(),
            );
        }
        let diagonal = -column_norm.copysign(x[j][j]);
        let mut reflector = x[j..]
            .iter()
            .map(|row| row[j] / column_norm)
            .collect::<Vec<_>>();
        reflector[0] -= diagonal / column_norm;
        let length = norm(reflector.iter().copied());
        for value in &mut reflector {
            *value /= length;
        }
        for k in j..p {
            let projection = 2.
                * reflector
                    .iter()
                    .enumerate()
                    .map(|(i, v)| v * x[j + i][k])
                    .sum::<f64>();
            for (i, v) in reflector.iter().enumerate() {
                x[j + i][k] -= projection * v;
            }
        }
        let projection = 2. * dot(&reflector, &y[j..]);
        for (i, v) in reflector.iter().enumerate() {
            y[j + i] -= projection * v;
        }
        x[j][j] = diagonal;
    }
    let mut coefficients = vec![0.; p];
    for j in (0..p).rev() {
        coefficients[j] = finite(
            (y[j] - dot(&x[j][j + 1..], &coefficients[j + 1..])) / x[j][j],
            "QR coefficient",
        )?;
    }
    Ok(coefficients)
}

fn pls(
    mut x: Vec<Vec<f64>>,
    mut y: Vec<f64>,
    components: usize,
) -> Result<(Vec<f64>, usize), String> {
    let p = x[0].len();
    let mut coefficients = vec![0.; p];
    let mut rotations: Vec<Vec<f64>> = Vec::new();
    let mut loadings: Vec<Vec<f64>> = Vec::new();
    let mut initial_covariance = 0.;
    for component in 0..components.min(x.len() - 1) {
        let mut weights = (0..p)
            .map(|j| {
                x.iter()
                    .zip(&y)
                    .map(|(row, target)| row[j] * target)
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        let covariance = finite(norm(weights.iter().copied()), "PLS covariance")?;
        if component == 0 {
            initial_covariance = covariance;
        }
        if covariance == 0. || covariance <= initial_covariance * 1e-12 {
            break;
        }
        for weight in &mut weights {
            *weight /= covariance;
        }
        let scores = x.iter().map(|row| dot(row, &weights)).collect::<Vec<_>>();
        let score_squared = finite(dot(&scores, &scores), "PLS score norm")?;
        if score_squared <= (x.len() as f64) * 1e-24 {
            break;
        }
        let loading = (0..p)
            .map(|j| {
                x.iter()
                    .zip(&scores)
                    .map(|(row, score)| row[j] * score)
                    .sum::<f64>()
                    / score_squared
            })
            .collect::<Vec<_>>();
        let response_loading = finite(dot(&y, &scores) / score_squared, "PLS response loading")?;
        let mut rotation = weights.clone();
        for (previous_rotation, previous_loading) in rotations.iter().zip(&loadings) {
            let projection = dot(previous_loading, &weights);
            for j in 0..p {
                rotation[j] -= previous_rotation[j] * projection;
            }
        }
        for j in 0..p {
            coefficients[j] = finite(
                coefficients[j] + rotation[j] * response_loading,
                "PLS coefficient",
            )?;
        }
        for (i, row) in x.iter_mut().enumerate() {
            for j in 0..p {
                row[j] -= scores[i] * loading[j];
            }
            y[i] -= scores[i] * response_loading;
        }
        rotations.push(rotation);
        loadings.push(loading);
    }
    Ok((coefficients, rotations.len()))
}

fn predict(model: &Model, row: &[f64]) -> Result<f64, String> {
    if model.origins.is_empty() {
        finite(
            model.intercept + dot(&model.coefficients, row),
            "prediction",
        )
    } else {
        let adjustment = model
            .coefficients
            .iter()
            .enumerate()
            .map(|(j, coefficient)| {
                coefficient * ((row[j] - model.origins[j]) - model.mean_offsets[j])
            })
            .sum::<f64>();
        finite(
            model.response_origin + (model.response_mean_offset + adjustment),
            "prediction",
        )
    }
}

fn metrics(actual: &[f64], predictions: &[f64]) -> Result<(f64, String), String> {
    let residual_norm = finite(
        norm(
            actual
                .iter()
                .zip(predictions)
                .map(|(y, prediction)| y - prediction),
        ),
        "residual norm",
    )?;
    let rmse = finite(residual_norm / (actual.len() as f64).sqrt(), "RMSE")?;
    let target_origin = actual[0];
    let target_mean_offset = mean(actual.iter().map(|y| y - target_origin))?;
    let target_norm = finite(
        norm(
            actual
                .iter()
                .map(|y| (y - target_origin) - target_mean_offset),
        ),
        "target norm",
    )?;
    let r2 = if target_norm == 0. {
        "null".to_string()
    } else {
        finite(1. - (residual_norm / target_norm).powi(2), "R squared")?.to_string()
    };
    Ok((rmse, r2))
}

fn cross_validate(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    if parameters.len() != 2 || data[0].len() < 3 {
        return Err("cross_validate requires [method (0 ridge, 1 PLS), hyperparameter] and rows [features..., target, fold_id]".into());
    }
    let method = match parameters[0] {
        0. => "ridge",
        1. => "pls",
        _ => return Err("method must be 0 (ridge) or 1 (PLS)".into()),
    };
    let width = data[0].len();
    let mut folds = BTreeMap::<u64, Vec<usize>>::new();
    for (i, row) in data.iter().enumerate() {
        folds
            .entry(integer(row[width - 1], "fold id")?)
            .or_default()
            .push(i);
    }
    if folds.len() < 2 || folds.len() > 1000 {
        return Err("cross validation requires between 2 and 1000 nonempty folds".into());
    }
    if folds
        .len()
        .saturating_mul(data.len())
        .saturating_mul(width - 2)
        .saturating_mul(width - 2)
        > MAX_WORK
    {
        return Err("cross validation exceeds numerical work limit".into());
    }
    let mut predictions = vec![0.; data.len()];
    let mut fold_metrics = Vec::new();
    for (&fold, indices) in &folds {
        let training = data
            .iter()
            .filter(|row| row[width - 1] != fold as f64)
            .map(|row| row[..width - 1].to_vec())
            .collect::<Vec<_>>();
        let model = fit(method, &training, parameters[1])
            .map_err(|error| format!("fold {fold}: {error}"))?;
        let mut actual = Vec::with_capacity(indices.len());
        let mut fold_predictions = Vec::with_capacity(indices.len());
        for &i in indices {
            predictions[i] = predict(&model, &data[i][..width - 2])?;
            actual.push(data[i][width - 2]);
            fold_predictions.push(predictions[i]);
        }
        let (rmse, r2) = metrics(&actual, &fold_predictions)?;
        fold_metrics.push(format!(
            "{{\"fold_id\":{fold},\"n_train\":{},\"n_test\":{},\"rmse\":{rmse},\"r2\":{r2}}}",
            training.len(),
            indices.len()
        ));
    }
    let actual = data.iter().map(|row| row[width - 2]).collect::<Vec<_>>();
    let (rmse, r2) = metrics(&actual, &predictions)?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"{method}_cross_validation\",\"n_samples\":{},\"n_features\":{},\"n_folds\":{},\"preprocessing\":\"training_fold_only\",\"out_of_fold_predictions\":{},\"rmse\":{rmse},\"r2\":{r2},\"fold_metrics\":[{}]}}", data.len(), width - 2, folds.len(), numbers(&predictions), fold_metrics.join(",")))
}

fn mix(mut value: u64) -> u64 {
    value = (value ^ (value >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    value = (value ^ (value >> 27)).wrapping_mul(0x94d049bb133111eb);
    value ^ (value >> 31)
}

fn group_folds(data: &[Vec<f64>], parameters: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || parameters.len() != 2 {
        return Err("group_folds requires one group id per row and [fold_count, seed]".into());
    }
    let count = integer(parameters[0], "fold count")?;
    let seed = integer(parameters[1], "seed")?;
    if !(2..=1000).contains(&count) {
        return Err("fold count must be between 2 and 1000".into());
    }
    let count = count as usize;
    let mut groups = BTreeMap::<u64, usize>::new();
    for row in data {
        *groups.entry(integer(row[0], "group id")?).or_default() += 1;
    }
    if groups.len() < count {
        return Err("number of distinct groups must be at least fold count".into());
    }
    let mut ordered = groups.into_iter().collect::<Vec<_>>();
    // Largest groups first; seeded tie-breaking depends only on identity, never row order.
    ordered.sort_by_key(|&(id, size)| (std::cmp::Reverse(size), mix(id ^ mix(seed)), id));
    let mut sizes = vec![0; count];
    let mut assignment = BTreeMap::<u64, usize>::new();
    for (id, size) in ordered {
        let fold = (0..count)
            .min_by_key(|&fold| (sizes[fold], fold))
            .unwrap_or(0);
        sizes[fold] += size;
        assignment.insert(id, fold);
    }
    let folds = data
        .iter()
        .map(|row| assignment[&(row[0] as u64)].to_string())
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!("{{\"version\":1,\"algorithm\":\"balanced_group_folds\",\"seed\":{seed},\"n_folds\":{count},\"n_groups\":{},\"folds\":[{folds}],\"fold_sizes\":[{}]}}", assignment.len(), sizes.iter().map(usize::to_string).collect::<Vec<_>>().join(",")))
}

fn numbers(values: &[f64]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(f64::to_string)
            .collect::<Vec<_>>()
            .join(",")
    )
}
