//! Explicitly bounded MANOVA, sequential CCA inference, MCA and ML factoring.
//!
//! MANOVA tests the whole centered numeric design against intercept-only;
//! arbitrary effect contrasts/repeated-response designs are not supplied.
//! CCA inference consumes all unregularized fitted correlations and uses an
//! asymptotic Bartlett chi-square approximation. MCA uses complete integer-coded
//! categories and reports indicator/Burt eigenvalue conventions explicitly.
//! Gaussian ML factor analysis bounds unique variances and reports convergence
//! and boundary diagnostics; no oblique rotation or missing-data ML is implied.

const MAX_CELLS: usize = 1_000_000;
const MAX_WORK: usize = 200_000_000;
const MAX_FEATURES: usize = 64;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 10000
        || width > 128
        || data.len() > MAX_CELLS / width
        || args.len() > 16384
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("latent inference requires finite rectangular data, <=10000 rows,128columns,1000000cells".into());
    }
    match operation {
        "manova" => manova(data, args),
        "canonical_correlation_inference" => canonical_inference(data, args),
        "multiple_correspondence_analysis" => mca(data, args),
        "factor_analysis_ml" => ml_factor(data, args),
        _ => Err("unknown latent inference operation".into()),
    }
}
fn manova(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() != 1 || data[0].len() < 2 || data[0].len() > 64 {
        return Err("manova requires numeric [X...,Y...] and [x_columns], <=64totalcolumns".into());
    }
    let n = data.len();
    let xcols = integer(args[0], 1, data[0].len() - 1, "X columns")?;
    let p = data[0].len() - xcols;
    if xcols > 32 || p > 32 || n <= xcols + 1 {
        return Err("MANOVA requires <=32predictors/32responses and positive residual df".into());
    }
    work(n, xcols * xcols + p * p, 1)?;
    let xf = Frame::new(
        &data.iter().map(|r| r[..xcols].to_vec()).collect::<Vec<_>>(),
        true,
    )?;
    let yf = Frame::new(
        &data.iter().map(|r| r[xcols..].to_vec()).collect::<Vec<_>>(),
        true,
    )?;
    let mut q: Vec<Vec<f64>> = Vec::with_capacity(xcols);
    for j in 0..xcols {
        let mut column = xf.x.iter().map(|r| r[j]).collect::<Vec<_>>();
        for _ in 0..2 {
            for previous in &q {
                let projection = dot(&column, previous);
                for (value, basis) in column.iter_mut().zip(previous) {
                    *value -= projection * basis;
                }
            }
        }
        let length = norm(column.iter().copied());
        if length <= 1e-10 * (n as f64).sqrt() {
            return Err("MANOVA centered design is rank deficient".into());
        }
        q.push(column.iter().map(|v| v / length).collect());
    }
    let mut fitted = vec![vec![0.; p]; n];
    for basis in &q {
        for j in 0..p {
            let coefficient = dot(basis, &yf.x.iter().map(|r| r[j]).collect::<Vec<_>>());
            for i in 0..n {
                fitted[i][j] += basis[i] * coefficient;
            }
        }
    }
    let residuals =
        yf.x.iter()
            .zip(&fitted)
            .map(|(y, m)| y.iter().zip(m).map(|(a, b)| a - b).collect())
            .collect::<Vec<Vec<f64>>>();
    let h = cross_covariance(&fitted, &fitted, 1)?;
    let e = cross_covariance(&residuals, &residuals, 1)?;
    if (0..p).map(|j| e[j][j]).sum::<f64>() <= 1e-24 * n as f64 {
        return Err("MANOVA has no identifiable residual variation".into());
    }
    let whitening = symmetric_power(&e, -0.5)?;
    let roots = eigen(&multiply(&multiply(&whitening, &h)?, &whitening)?)?;
    let largest = roots.values[0].abs().max(f64::MIN_POSITIVE);
    if roots.values.iter().any(|v| *v < -1e-10 * largest) {
        return Err("MANOVA hypothesis roots are numerically invalid".into());
    }
    let eigenvalues = roots.values.iter().map(|v| v.max(0.)).collect::<Vec<_>>();
    let log_wilks = -eigenvalues.iter().map(|v| v.ln_1p()).sum::<f64>();
    let wilks = log_wilks.exp();
    let pillai = eigenvalues.iter().map(|v| v / (1. + v)).sum::<f64>();
    let hotelling = finite(eigenvalues.iter().sum(), "Hotelling-Lawley trace")?;
    let roy = eigenvalues[0];
    let response = p as f64;
    let effect = xcols as f64;
    let residual_df = (n - xcols - 1) as f64;
    let s = response.min(effect);
    let m = ((response - effect).abs() - 1.) / 2.;
    let nn = (residual_df - response - 1.) / 2.;
    let mut f_values = vec![None; 4];
    let mut numerator = vec![None; 4];
    let mut denominator = vec![None; 4];
    let mut p_values = vec![None; 4];
    let exact = s == 1.;
    let candidates = if exact {
        let df1 = response * effect;
        let df2 = residual_df - response + 1.;
        vec![(hotelling * df2 / df1, df1, df2); 4]
    } else {
        let t = if response * response + effect * effect > 5. {
            ((response * response * effect * effect - 4.)
                / (response * response + effect * effect - 5.))
                .sqrt()
        } else {
            1.
        };
        let r = residual_df - (response - effect + 1.) / 2.;
        let u = (response * effect - 2.) / 4.;
        let df2 = r * t - 2. * u;
        let pq = response * effect;
        let maximum = response.max(effect);
        let roy_df2 = residual_df - maximum + effect;
        vec![
            ((-log_wilks / t).exp_m1() * df2 / pq, pq, df2),
            (
                pillai / (s - pillai) * (2. * nn + s + 1.) / (2. * m + s + 1.),
                s * (2. * m + s + 1.),
                s * (2. * nn + s + 1.),
            ),
            (
                2. * (s * nn + 1.) * hotelling / (s * s * (2. * m + s + 1.)),
                s * (2. * m + s + 1.),
                2. * (s * nn + 1.),
            ),
            (roy * roy_df2 / maximum, maximum, roy_df2),
        ]
    };
    for (index, (statistic, df1, df2)) in candidates.into_iter().enumerate() {
        if df1 > 0. && df2 > 0. && statistic.is_finite() && statistic >= 0. {
            f_values[index] = Some(statistic);
            numerator[index] = Some(df1);
            denominator[index] = Some(df2);
            p_values[index] = Some(f_tail(statistic, df1, df2)?);
        }
    }
    let original = |a: &[Vec<f64>]| -> Result<Vec<Vec<f64>>, String> {
        (0..p)
            .map(|i| {
                (0..p)
                    .map(|j| finite(a[i][j] * yf.scales[i] * yf.scales[j], "original SSCP"))
                    .collect()
            })
            .collect()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"whole_model_numeric_design_manova\",\"n_samples\":{n},\"x_columns\":{xcols},\"y_columns\":{p},\"hypothesis_df\":{xcols},\"error_df\":{residual_df},\"hypothesis_sscp\":{},\"error_sscp\":{},\"eigenvalues\":{},\"wilks_lambda\":{wilks},\"log_wilks_lambda\":{log_wilks},\"pillai_trace\":{pillai},\"hotelling_lawley_trace\":{hotelling},\"roy_largest_root\":{roy},\"test_order\":[\"wilks\",\"pillai\",\"hotelling_lawley\",\"roy\"],\"f_statistics\":{},\"numerator_df\":{},\"denominator_df\":{},\"p_values\":{},\"f_is_exact\":{exact},\"roy_p_value_is_lower_bound\":{},\"f_reference\":\"{}\",\"null_inference\":\"undefined_degrees_or_unrepresentable_F\"}}",matrix(&original(&h)?),matrix(&original(&e)?),numbers(&eigenvalues),optional_numbers(&f_values),optional_numbers(&numerator),optional_numbers(&denominator),optional_numbers(&p_values),!exact,if exact{"rank_one_exact_F"}else{"JMP_19_2_approximate_F_including_Roy_bound"}))
}
fn canonical_inference(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || args.len() != 3 {
        return Err("canonical_correlation_inference requires sorted unregularized correlations as one-column rows and [n_samples,p,q]".into());
    }
    let n = integer(args[0], 3, 1_000_000, "sample count")?;
    let p = integer(args[1], 1, 64, "X variables")?;
    let q = integer(args[2], 1, 64, "Y variables")?;
    if n <= p + q || data.len() != p.min(q) {
        return Err("CCA inference needs n>p+q and all min(p,q) canonical correlations".into());
    }
    let correlations = data.iter().map(|r| r[0]).collect::<Vec<_>>();
    if correlations.iter().any(|v| *v < 0. || *v > 1.)
        || correlations.windows(2).any(|v| v[0] < v[1])
    {
        return Err("canonical correlations must descend within [0,1]".into());
    }
    let factor = n as f64 - 1. - (p + q + 1) as f64 / 2.;
    let mut lambdas = Vec::new();
    let mut statistics = Vec::new();
    let mut dfs = Vec::new();
    let mut probabilities = Vec::new();
    for j in 0..correlations.len() {
        let log_lambda = correlations[j..]
            .iter()
            .map(|v| (-v * v).ln_1p())
            .sum::<f64>();
        let df = ((p - j) * (q - j)) as f64;
        lambdas.push(log_lambda.exp());
        dfs.push(df);
        if log_lambda.is_infinite() {
            statistics.push(None);
            probabilities.push(0.);
        } else {
            let statistic = finite(-factor * log_lambda, "Bartlett statistic")?;
            statistics.push(Some(statistic));
            probabilities.push(chi_square_tail(statistic, df)?);
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"canonical_correlation_inference\",\"method\":\"bartlett_asymptotic_chi_square\",\"exact\":false,\"n_samples\":{n},\"x_variables\":{p},\"y_variables\":{q},\"correlations\":{},\"wilks_lambdas\":{},\"bartlett_factor\":{factor},\"bartlett_chi_square\":{},\"degrees_of_freedom\":{},\"p_values\":{},\"null_statistic\":\"positive_infinity_for_a_perfect_correlation\",\"assumption\":\"complete_unregularized_sample_canonical_correlations_under_multivariate_normality\"}}",numbers(&correlations),numbers(&lambdas),optional_numbers(&statistics),numbers(&dfs),numbers(&probabilities)))
}
fn mca(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let n = data.len();
    let variables = data[0].len();
    if n < 2
        || !(2..=32).contains(&variables)
        || args.len() > 1
        || data
            .iter()
            .flatten()
            .any(|v| v.fract() != 0. || v.abs() > 9_007_199_254_740_991.)
    {
        return Err("MCA requires >=2rows,2..32integer-coded categorical variables and [dimensions=min(2,n-1,total_categories-variables)]".into());
    }
    let mut levels = Vec::with_capacity(variables);
    let mut category_variables = Vec::new();
    let mut category_levels = Vec::new();
    let mut indices = vec![vec![0; variables]; n];
    for j in 0..variables {
        let mut values = data.iter().map(|r| r[j]).collect::<Vec<_>>();
        values.sort_by(f64::total_cmp);
        values.dedup_by(|a, b| *a == *b);
        let offset = category_levels.len();
        if offset + values.len() > 64 {
            return Err("MCA supports at most64totalobservedcategories".into());
        }
        for (i, row) in data.iter().enumerate() {
            indices[i][j] = offset + values.iter().position(|v| *v == row[j]).unwrap();
        }
        category_variables.extend(vec![j as f64; values.len()]);
        category_levels.extend_from_slice(&values);
        levels.push(values);
    }
    let categories = category_levels.len();
    let max_axes = (categories - variables).min(n - 1);
    if max_axes == 0 {
        return Err("MCA has no nonconstant categorical dimensions".into());
    }
    let dimensions = integer(
        args.first().copied().unwrap_or(max_axes.min(2) as f64),
        1,
        max_axes,
        "MCA dimensions",
    )?;
    work(n, categories, categories)?;
    let mut counts = vec![0.; categories];
    let mut burt = vec![vec![0.; categories]; categories];
    for row in &indices {
        for a in row {
            counts[*a] += 1.;
            for b in row {
                burt[*a][*b] += 1.;
            }
        }
    }
    let masses = counts
        .iter()
        .map(|v| v / (n * variables) as f64)
        .collect::<Vec<_>>();
    let gram = (0..categories)
        .map(|i| {
            (0..categories)
                .map(|j| {
                    (burt[i][j] - counts[i] * counts[j] / n as f64)
                        / (variables as f64 * counts[i].sqrt() * counts[j].sqrt())
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let e = eigen(&gram)?;
    if e.values.iter().any(|v| *v < -1e-10) {
        return Err("MCA indicator Gram is numerically indefinite".into());
    }
    let inertias = e.values[..max_axes]
        .iter()
        .map(|v| if *v <= 1e-12 { 0. } else { *v })
        .collect::<Vec<_>>();
    let category_coordinates = (0..categories)
        .map(|i| {
            (0..dimensions)
                .map(|j| e.vectors[i][j] * inertias[j].sqrt() / masses[i].sqrt())
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let row_coordinates = indices
        .iter()
        .map(|row| {
            (0..dimensions)
                .map(|j| {
                    row.iter()
                        .map(|category| e.vectors[*category][j] / masses[*category].sqrt())
                        .sum::<f64>()
                        / variables as f64
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    // Numeric null axes have no identified individual scores; set them to zero.
    let row_coordinates = row_coordinates
        .into_iter()
        .map(|row| {
            row.iter()
                .enumerate()
                .map(|(j, v)| if inertias[j] == 0. { 0. } else { *v })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let total = inertias.iter().sum::<f64>();
    let adjusted = inertias
        .iter()
        .map(|v| {
            if *v > 1. / variables as f64 {
                ((variables as f64 / (variables - 1) as f64) * (v - 1. / variables as f64)).powi(2)
            } else {
                0.
            }
        })
        .collect::<Vec<_>>();
    let adjusted_total = adjusted.iter().sum::<f64>();
    let greenacre = (variables as f64 / (variables - 1) as f64)
        * (inertias.iter().map(|v| v * v).sum::<f64>()
            - (categories - variables) as f64 / (variables * variables) as f64);
    let proportions = |denominator: f64| {
        adjusted
            .iter()
            .map(|v| {
                if denominator <= 1e-12 {
                    None
                } else {
                    Some(v / denominator)
                }
            })
            .collect::<Vec<_>>()
    };
    let singular = inertias.iter().map(|v| v.sqrt()).collect::<Vec<_>>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"complete_indicator_burt_multiple_correspondence\",\"n_samples\":{n},\"variables\":{variables},\"total_categories\":{categories},\"dimensions\":{dimensions},\"category_variables\":{},\"category_levels\":{},\"category_counts\":{},\"category_masses\":{},\"burt_table\":{},\"inertias\":{},\"singular_values\":{},\"total_inertia\":{total},\"row_coordinates\":{},\"category_coordinates\":{},\"benzecri_inertias\":{},\"benzecri_proportions\":{},\"greenacre_total_inertia\":{},\"greenacre_proportions\":{},\"inertia_convention\":\"indicator_principal_inertias_equal_standardized_Burt_eigenvalues\",\"scaling\":\"principal\",\"null_adjusted_proportion\":\"zero_adjusted_total\"}}",numbers(&category_variables),numbers(&category_levels),numbers(&counts),numbers(&masses),matrix(&burt),numbers(&inertias),numbers(&singular),matrix(&row_coordinates),matrix(&category_coordinates),numbers(&adjusted),optional_numbers(&proportions(adjusted_total)),greenacre.max(0.),optional_numbers(&proportions(greenacre))))
}
fn ml_factor(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 7 || !(2..=24).contains(&data[0].len()) {
        return Err("factor_analysis_ml needs 2..24variables and [factors,rotation0none/1varimax=1,maxiter1..1000=200,tol=1e-8,uniqueness_floor=1e-6,input0rows/1correlation=0,n_samples_for_correlation=0]".into());
    }
    let p = data[0].len();
    let k = integer(args[0], 1, p - 1, "factor count")?;
    let df = ((p - k) * (p - k)) as i64 - (p + k) as i64;
    if df < 0 {
        return Err("ML factor model has negative degrees of freedom and is unidentified".into());
    }
    let df = df / 2;
    let rotate = unit(args.get(1).copied().unwrap_or(1.), "varimax")?;
    let maxiter = integer(
        args.get(2).copied().unwrap_or(200.),
        1,
        1000,
        "ML iterations",
    )?;
    let tol = tolerance(args.get(3).copied().unwrap_or(1e-8))?;
    let lower = args.get(4).copied().unwrap_or(1e-6);
    if !(1e-10..=1.).contains(&lower) {
        return Err("uniqueness floor must be in [1e-10,1]".into());
    }
    let correlation_input = unit(args.get(5).copied().unwrap_or(0.), "correlation input")?;
    let n_provided = integer(
        args.get(6).copied().unwrap_or(0.),
        0,
        1_000_000,
        "sample count",
    )?;
    if !correlation_input && n_provided != 0 {
        return Err("raw-row ML factor sample count is determined from rows; leave n_samples_for_correlation=0".into());
    }
    let n = if correlation_input {
        n_provided
    } else {
        data.len()
    };
    if n != 0 && n <= p {
        return Err("ML factor inference/row fitting needs n_samples>variables".into());
    }
    work(p, p, p * maxiter)?;
    let frame = if correlation_input {
        Frame {
            x: vec![],
            origins: vec![0.; p],
            offsets: vec![0.; p],
            scales: vec![1.; p],
        }
    } else {
        Frame::new(data, true)?
    };
    let r = if correlation_input {
        if data.len() != p {
            return Err("correlation input must be square".into());
        }
        for (i, row) in data.iter().enumerate() {
            if (row[i] - 1.).abs() > 1e-10 {
                return Err("correlation diagonal must be one".into());
            }
            for (j, value) in row.iter().enumerate() {
                if (*value - data[j][i]).abs() > 1e-10 || value.abs() > 1. + 1e-10 {
                    return Err("correlation input must be symmetric within [-1,1]".into());
                }
            }
        }
        data.to_vec()
    } else {
        cross_covariance(&frame.x, &frame.x, data.len())?
    };
    let r_inverse = symmetric_power(&r, -1.)?;
    let mut uniqueness = (0..p)
        .map(|j| (1. / r_inverse[j][j]).clamp(lower, 1.))
        .collect::<Vec<_>>();
    let mut current = factor_profile(&r, &uniqueness, k)?;
    let mut inverse_hessian = identity(p);
    let mut converged = false;
    let mut iterations = 0;
    let mut line_search_failed = false;
    let mut history = vec![current.criterion.max(0.)];
    for iteration in 0..maxiter {
        let projected = projected_gradient(&uniqueness, &current.gradient, lower);
        let maximum = projected.iter().map(|v| v.abs()).fold(0_f64, f64::max);
        if maximum <= tol {
            converged = true;
            break;
        }
        let mut direction = inverse_hessian
            .iter()
            .map(|row| -dot(row, &projected))
            .collect::<Vec<_>>();
        // Restrict quasi-Newton coupling to the current free subspace. An active
        // bound has zero projected gradient but can otherwise acquire a costly
        // inward direction from off-diagonal inverse-Hessian entries.
        let active = (0..p)
            .map(|j| active_bound(uniqueness[j], current.gradient[j], lower))
            .collect::<Vec<_>>();
        for j in 0..p {
            if active[j] {
                direction[j] = 0.;
            }
        }
        if dot(&direction, &projected) >= 0. {
            direction = projected.iter().map(|v| -v).collect();
            inverse_hessian = identity(p);
        }
        let mut alpha = 1.;
        let mut accepted = None;
        for _ in 0..60 {
            let candidate = uniqueness
                .iter()
                .zip(&direction)
                .map(|(x, d)| (x + alpha * d).clamp(lower, 1.))
                .collect::<Vec<_>>();
            let step = candidate
                .iter()
                .zip(&uniqueness)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
            let slope = dot(&current.gradient, &step);
            if slope < 0. {
                let evaluation = factor_profile(&r, &candidate, k)?;
                if evaluation.criterion
                    <= current.criterion
                        + 1e-4 * slope
                        + 64. * f64::EPSILON * (1. + current.criterion.abs())
                {
                    accepted = Some((candidate, evaluation, step));
                    break;
                }
            }
            alpha *= 0.5;
        }
        let Some((candidate, evaluation, step)) = accepted else {
            line_search_failed = true;
            break;
        };
        let mut change = evaluation
            .gradient
            .iter()
            .zip(&current.gradient)
            .map(|(a, b)| a - b)
            .collect::<Vec<_>>();
        let active_changed =
            (0..p).any(|j| active[j] != active_bound(candidate[j], evaluation.gradient[j], lower));
        for j in 0..p {
            if active[j] {
                change[j] = 0.;
            }
        }
        let curvature = dot(&step, &change);
        if !active_changed
            && curvature > 1e-12 * norm(step.iter().copied()) * norm(change.iter().copied())
            && curvature > 0.
        {
            let hy = inverse_hessian
                .iter()
                .map(|row| dot(row, &change))
                .collect::<Vec<_>>();
            let coefficient = (1. + dot(&change, &hy) / curvature) / curvature;
            for i in 0..p {
                for j in 0..p {
                    inverse_hessian[i][j] = finite(
                        inverse_hessian[i][j] + coefficient * step[i] * step[j]
                            - (step[i] * hy[j] + hy[i] * step[j]) / curvature,
                        "ML inverse Hessian",
                    )?;
                }
            }
        } else {
            inverse_hessian = identity(p);
        }
        uniqueness = candidate;
        current = evaluation;
        iterations = iteration + 1;
        history.push(current.criterion.max(0.));
    }
    let projected = projected_gradient(&uniqueness, &current.gradient, lower);
    let gradient_max = projected.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    if gradient_max <= tol {
        converged = true;
    }
    let boundary = uniqueness
        .iter()
        .any(|v| *v <= lower + 1e-8 * lower.max(1e-6) || *v >= 1. - 1e-8);
    let (loadings, rotation, rotation_converged) = if rotate && k > 1 {
        varimax(&current.loadings, maxiter, tol)?
    } else {
        (current.loadings.clone(), identity(k), true)
    };
    let common = multiply(&loadings, &transpose(&loadings))?;
    let communalities = (0..p).map(|j| common[j][j]).collect::<Vec<_>>();
    let mut reproduced = common;
    for j in 0..p {
        reproduced[j][j] += uniqueness[j];
    }
    let residual = norm(
        (0..p)
            .flat_map(|i| (0..p).map(move |j| (i, j)))
            .map(|(i, j)| r[i][j] - reproduced[i][j]),
    ) / ((p * p) as f64).sqrt();
    let score_weights = multiply(&symmetric_power(&reproduced, -1.)?, &loadings)?;
    let scores = if correlation_input {
        vec![]
    } else {
        multiply(&frame.x, &score_weights)?
    };
    let models = frame.score_models(&score_weights)?;
    let criterion = current.criterion.max(0.);
    let correction = n as f64 - 1. - (2 * p + 5) as f64 / 6. - 2. * k as f64 / 3.;
    let model_statistic =
        if n > 0 && df > 0 && correction > 0. && converged && rotation_converged && !boundary {
            Some(criterion * correction)
        } else {
            None
        };
    let model_p = match model_statistic {
        Some(statistic) => Some(chi_square_tail(statistic, df as f64)?),
        None => None,
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"bounded_uniqueness_gaussian_ML_factor_analysis\",\"factors\":{k},\"n_samples\":{n},\"model_df\":{df},\"uniqueness_floor\":{lower},\"converged\":{},\"optimizer_converged\":{converged},\"rotation_converged\":{rotation_converged},\"iterations\":{iterations},\"line_search_failed\":{line_search_failed},\"projected_gradient_max\":{gradient_max},\"boundary\":{boundary},\"criterion\":{criterion},\"criterion_history\":{},\"criterion_definition\":\"logdet_Sigma_plus_trace_R_Sigma_inverse_minus_logdet_R_minus_p\",\"loadings\":{},\"unrotated_loadings\":{},\"rotation_matrix\":{},\"uniquenesses\":{},\"communalities\":{},\"reproduced_correlation\":{},\"residual_rms\":{residual},\"model_chi_square\":{},\"model_p_value\":{},\"model_test\":\"Bartlett_corrected_asymptotic_chi_square_null_on_boundary_nonconvergence_or_no_df_n\",\"scores\":{},\"score_predict_parameters\":{},\"score_input\":\"{}\",\"score_method\":\"fitted_model_covariance_inverse_loadings\",\"optimizer\":\"projected_BFGS_SMC_start_local_optimum_not_guaranteed_global\"}}",converged&&rotation_converged,numbers(&history),matrix(&loadings),matrix(&current.loadings),matrix(&rotation),numbers(&uniqueness),numbers(&communalities),matrix(&reproduced),nullable(model_statistic),nullable(model_p),matrix(&scores),matrix(&models),if correlation_input{"standardized_variables"}else{"original_variables"}))
}
struct FactorProfile {
    criterion: f64,
    gradient: Vec<f64>,
    loadings: Vec<Vec<f64>>,
}
fn factor_profile(r: &[Vec<f64>], uniqueness: &[f64], k: usize) -> Result<FactorProfile, String> {
    let p = r.len();
    let scales = uniqueness.iter().map(|v| v.sqrt()).collect::<Vec<_>>();
    let scaled = (0..p)
        .map(|i| (0..p).map(|j| r[i][j] / scales[i] / scales[j]).collect())
        .collect::<Vec<Vec<f64>>>();
    let e = eigen(&scaled)?;
    if e.values.iter().any(|v| *v <= 0.) {
        return Err("ML factor profile has nonpositive correlation roots".into());
    }
    let loadings = (0..p)
        .map(|i| {
            (0..k)
                .map(|j| scales[i] * e.vectors[i][j] * (e.values[j] - 1.).max(0.).sqrt())
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let mut criterion = 0.;
    for (j, value) in e.values.iter().enumerate() {
        if j >= k || *value < 1. {
            let delta = value - 1.;
            criterion += if delta.abs() < 1e-5 {
                delta * delta * (0.5 - delta / 3. + delta * delta / 4. - delta.powi(3) / 5.)
            } else {
                delta - delta.ln_1p()
            };
        }
    }
    let inverse = (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    (0..p)
                        .map(|a| {
                            e.vectors[i][a] * e.vectors[j][a]
                                / if a < k { e.values[a].max(1.) } else { 1. }
                        })
                        .sum::<f64>()
                        / scales[i]
                        / scales[j]
                })
                .collect()
        })
        .collect::<Vec<Vec<f64>>>();
    let product = multiply(&inverse, r)?;
    let gradient = (0..p)
        .map(|i| finite(inverse[i][i] - dot(&product[i], &inverse[i]), "ML gradient"))
        .collect::<Result<Vec<_>, _>>()?;
    Ok(FactorProfile {
        criterion: finite(criterion, "ML likelihood criterion")?,
        gradient,
        loadings,
    })
}
fn projected_gradient(uniqueness: &[f64], gradient: &[f64], lower: f64) -> Vec<f64> {
    uniqueness
        .iter()
        .zip(gradient)
        .map(|(x, g)| if active_bound(*x, *g, lower) { 0. } else { *g })
        .collect()
}
fn active_bound(x: f64, gradient: f64, lower: f64) -> bool {
    (x <= lower + 1e-14 && gradient > 0.) || (x >= 1. - 1e-14 && gradient < 0.)
}
fn nullable(value: Option<f64>) -> String {
    value.map_or("null".into(), |v| v.to_string())
}
fn optional_numbers(values: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        values
            .iter()
            .map(|v| nullable(*v))
            .collect::<Vec<_>>()
            .join(",")
    )
}

fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numeric range"))
    }
}
fn integer(x: f64, lo: usize, hi: usize, name: &str) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("{name} must be an integer in [{lo},{hi}]"))
    } else {
        Ok(x as usize)
    }
}
fn unit(x: f64, name: &str) -> Result<bool, String> {
    match x {
        0. => Ok(false),
        1. => Ok(true),
        _ => Err(format!("{name} must be 0 or 1")),
    }
}
fn tolerance(x: f64) -> Result<f64, String> {
    if (1e-12..=1e-2).contains(&x) {
        Ok(x)
    } else {
        Err("tolerance must be in [1e-12,1e-2]".into())
    }
}
fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.saturating_mul(b).saturating_mul(c) > MAX_WORK {
        Err("latent computation exceeds numerical work limit".into())
    } else {
        Ok(())
    }
}
fn numbers(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(v: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}
fn norm(v: impl Iterator<Item = f64>) -> f64 {
    v.fold(0_f64, f64::hypot)
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn mean(values: impl Iterator<Item = f64>) -> Result<f64, String> {
    let mut m = 0.;
    for (i, v) in values.enumerate() {
        m = finite(
            m + finite(v - m, "mean difference")? / (i + 1) as f64,
            "mean",
        )?;
    }
    Ok(m)
}
fn identity(p: usize) -> Vec<Vec<f64>> {
    (0..p)
        .map(|i| (0..p).map(|j| f64::from(i == j)).collect())
        .collect()
}
fn transpose(a: &[Vec<f64>]) -> Vec<Vec<f64>> {
    (0..a[0].len())
        .map(|j| a.iter().map(|r| r[j]).collect())
        .collect()
}
fn multiply(a: &[Vec<f64>], b: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let bt = transpose(b);
    a.iter()
        .map(|r| {
            bt.iter()
                .map(|c| finite(dot(r, c), "matrix product"))
                .collect()
        })
        .collect()
}

struct Frame {
    x: Vec<Vec<f64>>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
}
impl Frame {
    fn new(data: &[Vec<f64>], standardize: bool) -> Result<Self, String> {
        let p = data[0].len();
        if p > MAX_FEATURES || data.len() < 2 {
            return Err("fit needs >=2 rows and <=64 features".into());
        }
        let origins = data[0].clone();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|r| r[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = if standardize {
                    finite(
                        norm(data.iter().map(|r| (r[j] - origins[j]) - offsets[j]))
                            / (data.len() as f64).sqrt(),
                        "feature scale",
                    )?
                } else {
                    1.
                };
                Ok(if s == 0. { 1. } else { s })
            })
            .collect::<Result<Vec<_>, String>>()?;
        let x = data
            .iter()
            .map(|r| {
                (0..p)
                    .map(|j| {
                        finite(
                            ((r[j] - origins[j]) - offsets[j]) / scales[j],
                            "standardized feature",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<_>, _>>()?;
        Ok(Self {
            x,
            origins,
            offsets,
            scales,
        })
    }
    fn score_models(&self, weights: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
        transpose(weights)
            .iter()
            .map(|column| {
                let mut model = vec![0., 0., 0.];
                model.extend_from_slice(&self.origins);
                model.extend_from_slice(&self.offsets);
                for (v, s) in column.iter().zip(&self.scales) {
                    model.push(finite(v / s, "score coefficient")?);
                }
                Ok(model)
            })
            .collect()
    }
}
fn cross_covariance(
    a: &[Vec<f64>],
    b: &[Vec<f64>],
    denominator: usize,
) -> Result<Vec<Vec<f64>>, String> {
    let p = a[0].len();
    let q = b[0].len();
    work(a.len(), p, q)?;
    (0..p)
        .map(|i| {
            (0..q)
                .map(|j| {
                    finite(
                        a.iter().zip(b).map(|(x, y)| x[i] * y[j]).sum::<f64>() / denominator as f64,
                        "covariance",
                    )
                })
                .collect()
        })
        .collect()
}

struct Eigen {
    values: Vec<f64>,
    vectors: Vec<Vec<f64>>,
}
#[allow(clippy::needless_range_loop)]
fn eigen(matrix: &[Vec<f64>]) -> Result<Eigen, String> {
    let p = matrix.len();
    work(p, p, 64 * p)?;
    let scale = matrix
        .iter()
        .flatten()
        .map(|v| v.abs())
        .fold(0_f64, f64::max);
    if scale == 0. {
        return Ok(Eigen {
            values: vec![0.; p],
            vectors: identity(p),
        });
    }
    let mut a = matrix
        .iter()
        .map(|r| r.iter().map(|v| v / scale).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    let mut v = identity(p);
    let mut converged = p == 1;
    for _ in 0..64 {
        let mut maximum = 0_f64;
        for i in 0..p {
            for j in i + 1..p {
                let off = a[i][j];
                maximum = maximum.max(off.abs());
                if off.abs() <= 2e-15 {
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
                for k in 0..p {
                    if k != i && k != j {
                        let x = a[k][i];
                        let y = a[k][j];
                        a[k][i] = c * x - s * y;
                        a[i][k] = a[k][i];
                        a[k][j] = s * x + c * y;
                        a[j][k] = a[k][j];
                    }
                    let x = v[k][i];
                    let y = v[k][j];
                    v[k][i] = c * x - s * y;
                    v[k][j] = s * x + c * y;
                }
            }
        }
        if maximum <= 1e-13 {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err("symmetric eigensolver did not converge".into());
    }
    let mut order = (0..p).collect::<Vec<_>>();
    order.sort_by(|i, j| a[*j][*j].total_cmp(&a[*i][*i]).then(i.cmp(j)));
    let values = order
        .iter()
        .map(|j| finite(a[*j][*j] * scale, "eigenvalue"))
        .collect::<Result<Vec<_>, _>>()?;
    let mut vectors = (0..p)
        .map(|i| order.iter().map(|j| v[i][*j]).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for j in 0..p {
        let anchor = (0..p)
            .max_by(|i, k| {
                vectors[*i][j]
                    .abs()
                    .total_cmp(&vectors[*k][j].abs())
                    .then(k.cmp(i))
            })
            .unwrap();
        if vectors[anchor][j] < 0. {
            for row in &mut vectors {
                row[j] = -row[j];
            }
        }
    }
    Ok(Eigen { values, vectors })
}
fn symmetric_power(a: &[Vec<f64>], power: f64) -> Result<Vec<Vec<f64>>, String> {
    let e = eigen(a)?;
    let p = a.len();
    let threshold = e.values[0].abs() * 1e-12;
    if e.values.iter().any(|v| *v <= threshold) {
        return Err(
            "matrix is singular or numerically unidentified; increase regularization".into(),
        );
    }
    (0..p)
        .map(|i| {
            (0..p)
                .map(|j| {
                    finite(
                        (0..p)
                            .map(|k| e.vectors[i][k] * e.values[k].powf(power) * e.vectors[j][k])
                            .sum(),
                        "symmetric matrix power",
                    )
                })
                .collect()
        })
        .collect()
}

type RotatedLoadings = (Vec<Vec<f64>>, Vec<Vec<f64>>, bool);
fn varimax(loadings: &[Vec<f64>], maxiter: usize, tol: f64) -> Result<RotatedLoadings, String> {
    let p = loadings.len();
    let k = loadings[0].len();
    let scales = loadings
        .iter()
        .map(|r| norm(r.iter().copied()))
        .collect::<Vec<_>>();
    let mut l = loadings
        .iter()
        .zip(&scales)
        .map(|(r, s)| {
            r.iter()
                .map(|v| if *s == 0. { 0. } else { v / s })
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    let mut rotation = identity(k);
    let mut converged = false;
    for _ in 0..maxiter {
        let mut max_angle = 0_f64;
        for a in 0..k {
            for b in a + 1..k {
                let u = l
                    .iter()
                    .map(|r| r[a] * r[a] - r[b] * r[b])
                    .collect::<Vec<_>>();
                let v = l.iter().map(|r| 2. * r[a] * r[b]).collect::<Vec<_>>();
                let su = u.iter().sum::<f64>();
                let sv = v.iter().sum::<f64>();
                let aa = dot(&u, &u) - dot(&v, &v) - (su * su - sv * sv) / p as f64;
                let bb = 2. * dot(&u, &v) - 2. * su * sv / p as f64;
                let angle = 0.25 * bb.atan2(aa);
                max_angle = max_angle.max(angle.abs());
                let c = angle.cos();
                let s = angle.sin();
                for r in l.iter_mut().chain(&mut rotation) {
                    let x = r[a];
                    let y = r[b];
                    r[a] = c * x + s * y;
                    r[b] = -s * x + c * y;
                }
            }
        }
        if max_angle <= tol {
            converged = true;
            break;
        }
    }
    for (i, r) in l.iter_mut().enumerate() {
        for v in r {
            *v = finite(*v * scales[i], "rotated loading")?;
        }
    }
    Ok((l, rotation, converged))
}

fn log_gamma(z: f64) -> f64 {
    let c = [
        676.5203681218851,
        -1259.1392167224028,
        771.3234287776531,
        -176.6150291621406,
        12.507343278686905,
        -0.13857109526572012,
        9.984369578019572e-6,
        1.5056327351493116e-7,
    ];
    if z < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * z).sin().ln()
            - log_gamma(1. - z);
    }
    let z = z - 1.;
    let mut x = 0.9999999999998099;
    for (i, c) in c.iter().enumerate() {
        x += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * (2. * std::f64::consts::PI).ln() + (z + 0.5) * t.ln() - t + x.ln()
}

fn beta_fraction(a: f64, b: f64, x: f64) -> Result<f64, String> {
    let minimum = 1e-300;
    let mut c = 1.;
    let mut d = 1. - (a + b) * x / (a + 1.);
    if d.abs() < minimum {
        d = minimum;
    }
    d = 1. / d;
    let mut h = d;
    for m in 1..=10_000 {
        let m = m as f64;
        let m2 = 2. * m;
        let aa = m * (b - m) * x / ((a - 1. + m2) * (a + m2));
        d = 1. + aa * d;
        if d.abs() < minimum {
            d = minimum;
        }
        c = 1. + aa / c;
        if c.abs() < minimum {
            c = minimum;
        }
        d = 1. / d;
        h *= d * c;
        let aa = -(a + m) * (a + b + m) * x / ((a + m2) * (a + 1. + m2));
        d = 1. + aa * d;
        if d.abs() < minimum {
            d = minimum;
        }
        c = 1. + aa / c;
        if c.abs() < minimum {
            c = minimum;
        }
        d = 1. / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.).abs() < 2e-14 {
            return finite(h, "beta fraction");
        }
    }
    Err("incomplete beta failed to converge".into())
}

fn regularized_beta(x: f64, a: f64, b: f64) -> Result<f64, String> {
    if x <= 0. {
        return Ok(0.);
    }
    if x >= 1. {
        return Ok(1.);
    }
    let factor =
        (log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * x.ln() + b * (-x).ln_1p()).exp();
    let result = if x < (a + 1.) / (a + b + 2.) {
        factor * beta_fraction(a, b, x)? / a
    } else {
        1. - factor * beta_fraction(b, a, 1. - x)? / b
    };
    Ok(finite(result, "beta probability")?.clamp(0., 1.))
}

fn f_tail(f: f64, df1: f64, df2: f64) -> Result<f64, String> {
    regularized_beta(1. / (1. + df1 / df2 * f), df2 / 2., df1 / 2.)
}
fn chi_square_tail(statistic: f64, df: f64) -> Result<f64, String> {
    if statistic == 0. {
        return Ok(1.);
    }
    let a = df / 2.;
    let x = statistic / 2.;
    let prefactor = (a * x.ln() - x - log_gamma(a)).exp();
    if x < a + 1. {
        let mut term = 1. / a;
        let mut sum = term;
        for i in 1..=10000 {
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() <= sum.abs() * 2e-14 {
                return Ok((1. - prefactor * sum).clamp(0., 1.));
            }
        }
    } else {
        let minimum = 1e-300;
        let mut b = x + 1. - a;
        let mut c = 1. / minimum;
        let mut d = 1. / b;
        let mut h = d;
        for i in 1..=10000 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.;
            d = an * d + b;
            if d.abs() < minimum {
                d = minimum;
            }
            c = b + an / c;
            if c.abs() < minimum {
                c = minimum;
            }
            d = 1. / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() < 2e-14 {
                return Ok(finite(prefactor * h, "chi-square tail")?.clamp(0., 1.));
            }
        }
    }
    Err("chi-square tail did not converge".into())
}
