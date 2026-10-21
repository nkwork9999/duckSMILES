//! Explicit CPU functional-data preprocessing; no GUI, implicit model fitting or
//! claim of numerical identity to undocumented JMP preprocessing defaults.
const MAX_WORK: usize = 200_000_000;
#[derive(Clone)]
struct Point {
    x: f64,
    y: f64,
    row: usize,
}
#[derive(Clone)]
struct Curve {
    id: f64,
    role: usize,
    points: Vec<Point>,
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if !(3..=4).contains(&width)
        || data.len() > 100000
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
    {
        return Err("functional preprocessing requires finite rectangular [x,y,id,optional_role] data and bounded arguments".into());
    }
    let mut curves: Vec<Curve> = Vec::new();
    for (i, row) in data.iter().enumerate() {
        valid_id(row[2])?;
        let role = if width == 4 {
            integer(row[3], 0, 1)?
        } else {
            0
        };
        let g = if let Some(g) = curves.iter().position(|c| c.id == row[2]) {
            if curves[g].role != role {
                return Err("one curve cannot mix training and validation roles".into());
            }
            g
        } else {
            if curves.len() == 128 {
                return Err("at most128 functions are supported".into());
            }
            curves.push(Curve {
                id: row[2],
                role,
                points: Vec::new(),
            });
            curves.len() - 1
        };
        curves[g].points.push(Point {
            x: row[0],
            y: row[1],
            row: i,
        });
        if curves[g].points.len() > 10000 {
            return Err("at most10000 observations per function".into());
        }
    }
    for c in &mut curves {
        c.points
            .sort_by(|a, b| a.x.total_cmp(&b.x).then(a.row.cmp(&b.row)));
    }
    match operation {
        "functional_grid" => grid(&curves, args),
        "functional_bin" => bin(&curves, args),
        "functional_thin" => thin(&curves, args),
        "functional_snip" => snip(&curves, args),
        "functional_arwls" => arwls(&curves, args),
        "functional_baseline" => supplied(&curves, args),
        "functional_dtw_matrix" => dtw_matrix(&curves, args),
        "functional_golden_match" => golden(&curves, args),
        _ => Err("unknown functional preprocessing operation".into()),
    }
}
fn finite(x: f64, what: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{what} exceeds finite numeric range"))
    }
}
fn integer(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.fract() == 0. && x >= lo as f64 && x <= hi as f64 {
        Ok(x as usize)
    } else {
        Err(format!("integer parameter must be in[{lo},{hi}]"))
    }
}
fn valid_id(x: f64) -> Result<(), String> {
    if x.fract() == 0. && x.abs() <= 9_007_199_254_740_991. {
        Ok(())
    } else {
        Err("function IDs must be exact integers within absolute2^53-1".into())
    }
}
fn num(x: f64) -> String {
    if x.is_finite() {
        x.to_string()
    } else {
        "null".into()
    }
}
fn numbers(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(|x| num(*x)).collect::<Vec<_>>().join(",")
    )
}
fn matrix(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| numbers(r)).collect::<Vec<_>>().join(",")
    )
}
fn bools(x: &[bool]) -> String {
    format!(
        "[{}]",
        x.iter().map(bool::to_string).collect::<Vec<_>>().join(",")
    )
}
fn ids(curves: &[Curve]) -> Vec<f64> {
    curves.iter().map(|c| c.id).collect()
}
fn roles(curves: &[Curve]) -> Vec<f64> {
    curves.iter().map(|c| c.role as f64).collect()
}
fn strict(c: &Curve, minimum: usize) -> Result<(), String> {
    if c.points.len() < minimum || c.points.windows(2).any(|p| p[0].x >= p[1].x) {
        Err("function requires sufficient strictly increasing distinct x values".into())
    } else {
        Ok(())
    }
}
fn uniform(c: &Curve, minimum: usize) -> Result<(), String> {
    strict(c, minimum)?;
    let first = c.points[0].x;
    let span = finite(c.points.last().unwrap().x - first, "grid span")?;
    let step = span / (c.points.len() - 1) as f64;
    let tol = 1e-8 * step.abs()
        + 8. * f64::EPSILON * c.points.iter().map(|p| p.x.abs()).fold(0_f64, f64::max);
    if c.points
        .iter()
        .enumerate()
        .any(|(i, p)| ((p.x - first) - i as f64 * step).abs() > tol)
    {
        Err(
            "baseline algorithm requires an evenly spaced x grid; interpolate explicitly first"
                .into(),
        )
    } else {
        Ok(())
    }
}
fn interpolate(c: &Curve, x: f64, policy: usize) -> Result<Option<f64>, String> {
    let p = &c.points;
    let position = p.partition_point(|p| p.x < x);
    if position < p.len() && p[position].x == x {
        return Ok(Some(p[position].y));
    }
    if position == 0 || position == p.len() {
        return match policy {
            0 => Ok(None),
            1 => Ok(Some(if position == 0 {
                p[0].y
            } else {
                p.last().unwrap().y
            })),
            _ => Err("query lies outside the function interpolation support".into()),
        };
    }
    let lo = &p[position - 1];
    let hi = &p[position];
    let fraction = finite(
        (x - lo.x) / finite(hi.x - lo.x, "interpolation interval")?,
        "interpolation fraction",
    )?;
    Ok(Some(finite(
        (1. - fraction) * lo.y + fraction * hi.y,
        "linear interpolation",
    )?))
}
fn grid(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if !matches!(args.len(), 0 | 1 | 4) {
        return Err("grid args [],[size],or[size,lower,upper,outside_policy]".into());
    }
    for c in curves {
        strict(c, 2)?;
    }
    let train = curves.iter().filter(|c| c.role == 0).collect::<Vec<_>>();
    if train.is_empty() && args.len() != 4 {
        return Err("automatic grid needs at least one training function".into());
    }
    let longest = train.iter().map(|c| c.points.len()).max().unwrap_or(0);
    let n = if args.len() == 4 {
        integer(args[0], 2, 4096)?
    } else if args.is_empty() {
        longest.min(4096)
    } else {
        integer(args[0], 2, longest.min(4096))?
    };
    let grid_estimation = if args.len() == 4 {
        "supplied_bounds_and_size"
    } else {
        "training_functions_only"
    };
    let lo = if args.len() == 4 {
        args[1]
    } else {
        train
            .iter()
            .map(|c| c.points[0].x)
            .fold(f64::NEG_INFINITY, f64::max)
    };
    let hi = if args.len() == 4 {
        args[2]
    } else {
        train
            .iter()
            .map(|c| c.points.last().unwrap().x)
            .fold(f64::INFINITY, f64::min)
    };
    if lo >= hi {
        return Err("common grid needs a positive shared domain or explicit bounds".into());
    }
    let span = finite(hi - lo, "common grid span")?;
    let policy = if args.len() == 4 {
        integer(args[3], 0, 2)?
    } else {
        0
    };
    let grid = (0..n)
        .map(|i| {
            if i == 0 {
                lo
            } else if i == n - 1 {
                hi
            } else {
                lo + span * (i as f64 / (n - 1) as f64)
            }
        })
        .collect::<Vec<_>>();
    if grid.windows(2).any(|v| v[0] >= v[1]) {
        return Err("grid is finer than representable x spacing".into());
    }
    let values = curves
        .iter()
        .map(|c| {
            grid.iter()
                .map(|x| interpolate(c, *x, policy).map(|v| v.unwrap_or(f64::NAN)))
                .collect()
        })
        .collect::<Result<Vec<Vec<f64>>, String>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"training_domain_common_grid_linear_interpolation\",\"ids\":{},\"roles\":{},\"grid\":{},\"values\":{},\"outside_policy\":{policy},\"grid_estimation\":\"{grid_estimation}\"}}",numbers(&ids(curves)),numbers(&roles(curves)),numbers(&grid),matrix(&values)))
}
fn bin(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if !(3..=4).contains(&args.len()) {
        return Err("bin args[number,lower,upper,outside0reject/1drop=1]".into());
    }
    let n = integer(args[0], 1, 4096)?;
    let (lo, hi) = (args[1], args[2]);
    if lo >= hi {
        return Err("bin bounds must increase".into());
    }
    let span = finite(hi - lo, "bin span")?;
    let policy = integer(args.get(3).copied().unwrap_or(1.), 0, 1)?;
    let mut values = Vec::new();
    let mut counts = Vec::new();
    let mut dropped = Vec::new();
    for c in curves {
        let mut mean = vec![0.; n];
        let mut count = vec![0.; n];
        let mut removed = 0.;
        for p in &c.points {
            if p.x < lo || p.x > hi {
                if policy == 0 {
                    return Err("observation is outside explicit bin range".into());
                }
                removed += 1.;
                continue;
            }
            let index = if p.x == hi {
                n - 1
            } else {
                (((p.x - lo) / span) * n as f64).floor() as usize
            };
            let index = index.min(n - 1);
            count[index] += 1.;
            mean[index] = finite(
                mean[index] * ((count[index] - 1.) / count[index]) + p.y / count[index],
                "bin average",
            )?;
        }
        for i in 0..n {
            if count[i] == 0. {
                mean[i] = f64::NAN;
            }
        }
        values.push(mean);
        counts.push(count);
        dropped.push(removed);
    }
    let centers = (0..n)
        .map(|i| lo + span * ((i as f64 + 0.5) / n as f64))
        .collect::<Vec<_>>();
    let edges = (0..=n)
        .map(|i| {
            if i == n {
                hi
            } else {
                lo + span * (i as f64 / n as f64)
            }
        })
        .collect::<Vec<_>>();
    if edges.windows(2).any(|p| p[0] >= p[1]) || centers.windows(2).any(|p| p[0] >= p[1]) {
        return Err("bin edges/centers are finer than representable x spacing".into());
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"explicit_equal_width_function_bins\",\"ids\":{},\"roles\":{},\"bin_edges\":{},\"grid\":{},\"values\":{},\"counts\":{},\"dropped\":{},\"edge_contract\":\"left_closed_right_open_except_final_right_closed\"}}",numbers(&ids(curves)),numbers(&roles(curves)),numbers(&edges),numbers(&centers),matrix(&values),matrix(&counts),numbers(&dropped)))
}
fn thin(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if args.len() > 3 {
        return Err("thinning args[stride=2,offset=0,keep_last=0]".into());
    }
    let stride = integer(args.first().copied().unwrap_or(2.), 1, 10000)?;
    let offset = integer(args.get(1).copied().unwrap_or(0.), 0, stride - 1)?;
    let last = integer(args.get(2).copied().unwrap_or(0.), 0, 1)?;
    let mut x = Vec::new();
    let mut y = Vec::new();
    let mut rows = Vec::new();
    for c in curves {
        let selected = c
            .points
            .iter()
            .enumerate()
            .filter(|(i, _)| i % stride == offset || last == 1 && *i == c.points.len() - 1)
            .map(|(_, p)| p)
            .collect::<Vec<_>>();
        x.push(selected.iter().map(|p| p.x).collect());
        y.push(selected.iter().map(|p| p.y).collect());
        rows.push(selected.iter().map(|p| p.row as f64).collect());
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"sorted_per_function_stride_selection\",\"ids\":{},\"roles\":{},\"x\":{},\"values\":{},\"original_row_indices\":{},\"stride\":{stride},\"offset\":{offset}}}",numbers(&ids(curves)),numbers(&roles(curves)),matrix(&x),matrix(&y),matrix(&rows)))
}
fn baseline_json(
    curves: &[Curve],
    baselines: &[Vec<f64>],
    corrected: &[Vec<f64>],
    applied: &[bool],
    extra: &str,
) -> String {
    format!("{{\"version\":1,\"ids\":{},\"roles\":{},\"x\":{},\"original\":{},\"baselines\":{},\"corrected\":{},\"baseline_applied\":{},\"validation_contract\":\"unchanged_observations_zero_placeholder_baseline_not_estimated\",{extra}}}",numbers(&ids(curves)),numbers(&roles(curves)),matrix(&curves.iter().map(|c|c.points.iter().map(|p|p.x).collect()).collect::<Vec<_>>()),matrix(&curves.iter().map(|c|c.points.iter().map(|p|p.y).collect()).collect::<Vec<_>>()),matrix(baselines),matrix(corrected),bools(applied))
}
fn snip(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 3 {
        return Err(
            "SNIP args[max_halfwindow,transform0raw/1loglog=1,direction0increasing/1decreasing=0]"
                .into(),
        );
    }
    let span = integer(args[0], 1, 4999)?;
    let transform = integer(args.get(1).copied().unwrap_or(1.), 0, 1)?;
    let direction = integer(args.get(2).copied().unwrap_or(0.), 0, 1)?;
    if curves
        .iter()
        .filter(|c| c.role == 0)
        .map(|c| c.points.len())
        .sum::<usize>()
        .saturating_mul(span)
        > MAX_WORK
    {
        return Err("SNIP clipping work limit exceeded".into());
    }
    let mut baselines = Vec::new();
    let mut corrected = Vec::new();
    let mut applied = Vec::new();
    for c in curves {
        if c.role == 1 {
            baselines.push(vec![0.; c.points.len()]);
            corrected.push(c.points.iter().map(|p| p.y).collect());
            applied.push(false);
            continue;
        }
        uniform(c, 2 * span + 1)?;
        if transform == 1 && c.points.iter().any(|p| p.y < 0.) {
            return Err("Ryan log-log SNIP requires nonnegative observations".into());
        }
        let original = c
            .points
            .iter()
            .map(|p| {
                if transform == 1 {
                    p.y.ln_1p().ln_1p()
                } else {
                    p.y
                }
            })
            .collect::<Vec<_>>();
        let mut z = original.clone();
        for iteration in 0..span {
            let width = if direction == 0 {
                iteration + 1
            } else {
                span - iteration
            };
            let mut next = z.clone();
            for i in width..z.len() - width {
                next[i] = z[i].min(z[i - width] / 2. + z[i + width] / 2.);
            }
            z = next;
        }
        let baseline = z
            .iter()
            .enumerate()
            .map(|(i, v)| {
                let value = if *v == original[i] {
                    c.points[i].y
                } else if transform == 1 {
                    v.exp_m1().exp_m1().min(c.points[i].y).max(0.)
                } else {
                    *v
                };
                finite(value, "SNIP baseline")
            })
            .collect::<Result<Vec<_>, _>>()?;
        let correction = c
            .points
            .iter()
            .zip(&baseline)
            .map(|(p, b)| finite(p.y - b, "SNIP corrected signal"))
            .collect::<Result<Vec<_>, _>>()?;
        baselines.push(baseline);
        corrected.push(correction);
        applied.push(true);
    }
    Ok(baseline_json(curves,&baselines,&corrected,&applied,&format!("\"algorithm\":\"synchronous_SNIP_peak_clipping\",\"max_halfwindow\":{span},\"transform\":{transform},\"direction\":{direction},\"scope\":\"explicit_span_variant_not_original_FWHM_low_statistics_filter_or_JMP_defaults\"")))
}
fn supplied(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if !(1..=2).contains(&args.len()) {
        return Err("supplied baseline args[baseline_id,outside0null/1constant/2reject=2]".into());
    }
    valid_id(args[0])?;
    let policy = integer(args.get(1).copied().unwrap_or(2.), 0, 2)?;
    let baseline = curves
        .iter()
        .find(|c| c.id == args[0])
        .ok_or("unknown baseline function ID")?;
    if baseline.role == 1 {
        return Err("validation function cannot become an estimated training baseline".into());
    }
    strict(baseline, 2)?;
    let output = curves
        .iter()
        .filter(|c| c.id != args[0])
        .cloned()
        .collect::<Vec<_>>();
    if output.is_empty() {
        return Err("supplied baseline needs at least one other function".into());
    }
    let mut baselines = Vec::new();
    let mut corrected = Vec::new();
    let mut applied = Vec::new();
    for c in &output {
        if c.role == 1 {
            baselines.push(vec![0.; c.points.len()]);
            corrected.push(c.points.iter().map(|p| p.y).collect());
            applied.push(false);
            continue;
        }
        let mut b = Vec::new();
        let mut y = Vec::new();
        for p in &c.points {
            let base = interpolate(baseline, p.x, policy)?;
            b.push(base.unwrap_or(f64::NAN));
            y.push(if let Some(v) = base {
                finite(p.y - v, "baseline-subtracted signal")?
            } else {
                f64::NAN
            });
        }
        baselines.push(b);
        corrected.push(y);
        applied.push(true);
    }
    Ok(baseline_json(&output,&baselines,&corrected,&applied,&format!("\"algorithm\":\"supplied_baseline_linear_interpolation_subtraction\",\"baseline_id\":{},\"baseline_excluded\":true",args[0])))
}
// Lower bandwidth two Cholesky factor of W+lambda*D2' D2. No moving average,
// hidden diagonal regularization or dense cubic solve is substituted.
fn whittaker(y: &[f64], weights: &[f64], lambda: f64) -> Result<(Vec<f64>, f64), String> {
    let n = y.len();
    let mut diagonal = weights.to_vec();
    let mut a1 = vec![0.; n];
    let mut a2 = vec![0.; n];
    for k in 0..n - 2 {
        diagonal[k] += lambda;
        diagonal[k + 1] += 4. * lambda;
        diagonal[k + 2] += lambda;
        a1[k + 1] -= 2. * lambda;
        a1[k + 2] -= 2. * lambda;
        a2[k + 2] += lambda;
    }
    let mut diag = vec![0.; n];
    let mut l1 = vec![0.; n];
    let mut l2 = vec![0.; n];
    for i in 0..n {
        if i >= 2 {
            l2[i] = a2[i] / diag[i - 2];
        }
        if i >= 1 {
            l1[i] = (a1[i] - if i >= 2 { l2[i] * l1[i - 1] } else { 0. }) / diag[i - 1];
        }
        let pivot = diagonal[i] - l1[i] * l1[i] - l2[i] * l2[i];
        if pivot <= 0. || !pivot.is_finite() {
            return Err(
                "arWLS weighted penalty is numerically singular; no hidden ridge applied".into(),
            );
        }
        diag[i] = pivot.sqrt();
    }
    let rhs = y
        .iter()
        .zip(weights)
        .map(|(y, w)| y * w)
        .collect::<Vec<_>>();
    let mut z = rhs.clone();
    for i in 0..n {
        z[i] = (z[i]
            - if i >= 1 { l1[i] * z[i - 1] } else { 0. }
            - if i >= 2 { l2[i] * z[i - 2] } else { 0. })
            / diag[i];
    }
    for i in (0..n).rev() {
        z[i] = finite(
            (z[i]
                - if i + 1 < n { l1[i + 1] * z[i + 1] } else { 0. }
                - if i + 2 < n { l2[i + 2] * z[i + 2] } else { 0. })
                / diag[i],
            "arWLS penalized solution",
        )?;
    }
    let mut residual = 0_f64;
    let mut anorm = 0_f64;
    for i in 0..n {
        let mut az = diagonal[i] * z[i];
        let mut rownorm = diagonal[i].abs();
        for (j, a) in [
            (i.wrapping_sub(1), a1[i]),
            (i.wrapping_sub(2), a2[i]),
            (i + 1, if i + 1 < n { a1[i + 1] } else { 0. }),
            (i + 2, if i + 2 < n { a2[i + 2] } else { 0. }),
        ] {
            if j < n {
                az += a * z[j];
                rownorm += a.abs();
            }
        }
        residual = residual.max((az - rhs[i]).abs());
        anorm = anorm.max(rownorm);
    }
    let denominator = anorm * z.iter().map(|v| v.abs()).fold(0_f64, f64::max)
        + rhs.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    Ok((
        z,
        if denominator == 0. {
            0.
        } else {
            residual / denominator
        },
    ))
}
fn logistic_negative(x: f64) -> f64 {
    if x >= 0. {
        let e = (-x).exp();
        e / (1. + e)
    } else {
        1. / (1. + x.exp())
    }
}
fn arwls(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if args.len() > 4 {
        return Err(
            "arWLS args[lambda=1e5,relative_weight_tolerance=1e-6,maxiter=100,negative_std_ddof=0]"
                .into(),
        );
    }
    let lambda = args.first().copied().unwrap_or(1e5);
    let tol = args.get(1).copied().unwrap_or(1e-6);
    let maxiter = integer(args.get(2).copied().unwrap_or(100.), 1, 1000)?;
    let ddof = integer(args.get(3).copied().unwrap_or(0.), 0, 1)?;
    if !(1e-8..=1e10).contains(&lambda) || !(1e-12..=1e-2).contains(&tol) {
        return Err("arWLS lambda or weight tolerance out of supported range".into());
    }
    let ntrain = curves
        .iter()
        .filter(|c| c.role == 0)
        .map(|c| c.points.len())
        .sum::<usize>();
    if ntrain.saturating_mul(maxiter).saturating_mul(16) > MAX_WORK {
        return Err("arWLS iteration work limit exceeded".into());
    }
    let mut baselines = Vec::new();
    let mut corrected = Vec::new();
    let mut applied = Vec::new();
    let mut all_weights = Vec::new();
    let mut convergence = Vec::new();
    let mut iterations = Vec::new();
    let mut ratios = Vec::new();
    let mut backward = Vec::new();
    let mut reasons = Vec::new();
    for c in curves {
        let n = c.points.len();
        if c.role == 1 {
            baselines.push(vec![0.; n]);
            corrected.push(c.points.iter().map(|p| p.y).collect());
            all_weights.push(vec![0.; n]);
            applied.push(false);
            convergence.push("null".to_string());
            iterations.push(0.);
            ratios.push(f64::NAN);
            backward.push(f64::NAN);
            reasons.push("\"validation_not_fitted\"".to_string());
            continue;
        }
        uniform(c, 3)?;
        let origin = c.points[0].y;
        let differences = c
            .points
            .iter()
            .map(|p| finite(p.y - origin, "origin-relative intensity"))
            .collect::<Result<Vec<_>, _>>()?;
        let scale = differences.iter().map(|v| v.abs()).fold(0_f64, f64::max);
        let scale = if scale == 0. { 1. } else { scale };
        let y = differences.iter().map(|v| v / scale).collect::<Vec<_>>();
        let affine = y.iter().enumerate().all(|(i, v)| {
            (v - y[n - 1] * (i as f64 / (n - 1) as f64)).abs() <= 128. * f64::EPSILON
        });
        let mut weights = vec![1.; n];
        let mut z = y.clone();
        let mut converged = affine;
        let mut reason = if affine {
            "affine_within_numeric_tolerance"
        } else {
            "iteration_limit"
        };
        let mut count = 0;
        let mut ratio = if affine { 0. } else { f64::INFINITY };
        let mut error = 0.;
        if !affine {
            for iteration in 1..=maxiter {
                let (solution, err) = whittaker(&y, &weights, lambda)?;
                z = solution;
                error = err;
                count = iteration;
                let residual = y.iter().zip(&z).map(|(y, z)| y - z).collect::<Vec<_>>();
                let negative = residual
                    .iter()
                    .filter(|v| **v < 0.)
                    .copied()
                    .collect::<Vec<_>>();
                if negative.len() <= ddof.max(1) {
                    reason = "negative_noise_unidentified";
                    break;
                }
                let mean = negative.iter().sum::<f64>() / negative.len() as f64;
                let sd = (negative.iter().map(|v| (v - mean).powi(2)).sum::<f64>()
                    / (negative.len() - ddof) as f64)
                    .sqrt();
                if sd == 0. || !sd.is_finite() {
                    reason = "negative_noise_scale_zero";
                    break;
                }
                let proposed = residual
                    .iter()
                    .map(|v| logistic_negative(2. * ((*v - (-mean + 2. * sd)) / sd)))
                    .collect::<Vec<_>>();
                let norm = weights.iter().fold(0_f64, |a, b| a.hypot(*b));
                ratio = finite(
                    weights
                        .iter()
                        .zip(&proposed)
                        .fold(0_f64, |a, (u, v)| a.hypot(u - v))
                        / norm,
                    "relative arWLS weight change",
                )?;
                if ratio <= tol && error <= 1e-10 {
                    converged = true;
                    reason = "relative_weight_convergence";
                    break;
                }
                if error > 1e-10 {
                    reason = "linear_solver_backward_error";
                    break;
                }
                if iteration < maxiter {
                    weights = proposed;
                }
            }
        }
        let baseline = if affine {
            c.points.iter().map(|p| p.y).collect()
        } else {
            z.iter()
                .map(|z| finite(origin + scale * z, "arWLS baseline"))
                .collect::<Result<Vec<_>, _>>()?
        };
        let correction = c
            .points
            .iter()
            .zip(&baseline)
            .map(|(p, b)| finite(p.y - b, "arWLS corrected signal"))
            .collect::<Result<Vec<_>, _>>()?;
        baselines.push(baseline);
        corrected.push(correction);
        all_weights.push(weights);
        applied.push(true);
        convergence.push(converged.to_string());
        iterations.push(count as f64);
        ratios.push(ratio);
        backward.push(error);
        reasons.push(format!("\"{reason}\""));
    }
    Ok(baseline_json(curves,&baselines,&corrected,&applied,&format!("\"algorithm\":\"Baek2015_asymmetrically_reweighted_second_difference_penalized_LS\",\"weights\":{},\"converged\":[{}],\"iterations\":{},\"weight_relative_changes\":{},\"linear_backward_errors\":{},\"stop_reasons\":[{}],\"lambda\":{lambda},\"negative_std_ddof\":{ddof},\"weight_contract\":\"returned_weights_solve_returned_baseline;ratio_compares_proposed_next_weights\"",matrix(&all_weights),convergence.join(","),numbers(&iterations),numbers(&ratios),numbers(&backward),reasons.join(","))))
}
fn rms(values: &[f64]) -> f64 {
    let peak = values.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    if peak == 0. {
        0.
    } else {
        peak * values
            .iter()
            .map(|v| (v / peak).powi(2) / values.len() as f64)
            .sum::<f64>()
            .sqrt()
    }
}
fn normalized(c: &Curve, mode: usize) -> Result<Vec<f64>, String> {
    let y = c.points.iter().map(|p| p.y).collect::<Vec<_>>();
    if mode == 0 {
        return Ok(y);
    }
    let origin = y[0];
    let differences = y.iter().map(|v| v - origin).collect::<Vec<_>>();
    let offset = differences.iter().map(|v| v / y.len() as f64).sum::<f64>();
    let centered = differences.iter().map(|v| v - offset).collect::<Vec<_>>();
    let values = if centered.iter().all(|v| v.is_finite()) {
        centered
    } else {
        let peak = y.iter().map(|v| v.abs()).fold(0_f64, f64::max);
        let scaled = y.iter().map(|v| v / peak).collect::<Vec<_>>();
        let mean = scaled.iter().map(|v| v / y.len() as f64).sum::<f64>();
        scaled.iter().map(|v| v - mean).collect()
    };
    let peak = values.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    if peak == 0. {
        return Err("cannot standardize a constant function".into());
    }
    let scaled = values.iter().map(|v| v / peak).collect::<Vec<_>>();
    let scale = rms(&scaled);
    if scale == 0. || !scale.is_finite() {
        return Err("function scale is not numerically identified".into());
    }
    Ok(scaled.iter().map(|v| v / scale).collect())
}
fn dtw_distance(a: &[f64], b: &[f64], window: usize, metric: usize) -> Result<f64, String> {
    let (n, m) = (a.len(), b.len());
    if n.saturating_mul(m) > 4_000_000 || window < n.abs_diff(m) {
        return Err("DTW window cannot reach endpoints or pair exceeds4000000cells".into());
    }
    let origin = a[0];
    let shifted = a.iter().chain(b).map(|v| v - origin).collect::<Vec<_>>();
    let scale = if shifted.iter().all(|v| v.is_finite()) {
        shifted.iter().map(|v| v.abs()).fold(0_f64, f64::max)
    } else {
        a.iter().chain(b).map(|v| v.abs()).fold(0_f64, f64::max)
    };
    let scale = if scale == 0. { 1. } else { scale };
    let mut previous = vec![f64::INFINITY; m + 1];
    previous[0] = 0.;
    for i in 1..=n {
        let mut current = vec![f64::INFINITY; m + 1];
        for j in 1.max(i.saturating_sub(window))..=m.min(i.saturating_add(window)) {
            let delta = a[i - 1] - b[j - 1];
            let difference = if delta.is_finite() {
                delta / scale
            } else {
                a[i - 1] / scale - b[j - 1] / scale
            };
            let cost = if metric == 0 {
                difference * difference
            } else {
                difference.abs()
            };
            current[j] = cost + previous[j - 1].min(previous[j]).min(current[j - 1]);
        }
        previous = current;
    }
    finite(
        if metric == 0 {
            previous[m].sqrt() * scale
        } else {
            previous[m] * scale
        },
        "DTW distance",
    )
}
fn window(args: &[f64], position: usize, maxlen: usize) -> Result<usize, String> {
    let v = args.get(position).copied().unwrap_or(-1.);
    if v == -1. {
        Ok(maxlen)
    } else {
        integer(v, 0, maxlen)
    }
}
fn dtw_matrix(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if args.len() > 3 || curves.len() < 2 {
        return Err("DTW matrix needs >=2functions; args[window=-1,normalize0none/1population=0,cost0squared/1absolute=0]".into());
    }
    let maxlen = curves.iter().map(|c| c.points.len()).max().unwrap();
    let window = window(args, 0, maxlen)?;
    let normalize = integer(args.get(1).copied().unwrap_or(0.), 0, 1)?;
    let metric = integer(args.get(2).copied().unwrap_or(0.), 0, 1)?;
    let mut work = 0_usize;
    for i in 0..curves.len() {
        for j in 0..i {
            work = work.saturating_add(
                curves[i]
                    .points
                    .len()
                    .saturating_mul(curves[j].points.len()),
            );
        }
    }
    if work > MAX_WORK {
        return Err("all-pairs DTW work limit exceeded".into());
    }
    let values = curves
        .iter()
        .map(|c| normalized(c, normalize))
        .collect::<Result<Vec<_>, _>>()?;
    let mut distances = vec![vec![0.; curves.len()]; curves.len()];
    for i in 0..curves.len() {
        for j in 0..i {
            let d = dtw_distance(&values[i], &values[j], window, metric)?;
            distances[i][j] = d;
            distances[j][i] = d;
        }
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"all_pairs_Sakoe_Chiba_DTW\",\"ids\":{},\"roles\":{},\"distances\":{},\"window_radius\":{window},\"normalization\":{normalize},\"cost_metric\":{metric},\"nominal_cell_pairs\":{work},\"x_contract\":\"x_sort_order_only_no_physical_time_cost\"}}",numbers(&ids(curves)),numbers(&roles(curves)),matrix(&distances)))
}
fn golden(curves: &[Curve], args: &[f64]) -> Result<String, String> {
    if args.is_empty() || args.len() > 4 || curves.len() < 2 {
        return Err("golden matching needs target plus candidates; args[target_id,metric0RMSE/1MAE/2DTW=0,window=-1,normalize=0]".into());
    }
    valid_id(args[0])?;
    let target = curves
        .iter()
        .find(|c| c.id == args[0])
        .ok_or("unknown golden target function")?;
    let metric = integer(args.get(1).copied().unwrap_or(0.), 0, 2)?;
    let maxlen = curves.iter().map(|c| c.points.len()).max().unwrap();
    let window = window(args, 2, maxlen)?;
    let normalize = integer(args.get(3).copied().unwrap_or(0.), 0, 1)?;
    let t = normalized(target, normalize)?;
    let candidates = curves
        .iter()
        .filter(|c| c.id != target.id)
        .collect::<Vec<_>>();
    let work = candidates
        .iter()
        .map(|c| c.points.len().saturating_mul(t.len()))
        .fold(0_usize, usize::saturating_add);
    if metric == 2 && work > MAX_WORK {
        return Err("golden DTW work limit exceeded".into());
    }
    let mut scores = Vec::new();
    for c in &candidates {
        let y = normalized(c, normalize)?;
        let score = if metric == 2 {
            dtw_distance(&t, &y, window, 0)?
        } else {
            if c.points.len() != target.points.len()
                || c.points.iter().zip(&target.points).any(|(a, b)| a.x != b.x)
            {
                return Err(
                    "pointwise golden matching requires an identical explicit x grid".into(),
                );
            }
            let differences = y
                .iter()
                .zip(&t)
                .map(|(y, t)| finite(y - t, "golden difference"))
                .collect::<Result<Vec<_>, _>>()?;
            if metric == 0 {
                rms(&differences)
            } else {
                differences.iter().map(|v| v.abs() / y.len() as f64).sum()
            }
        };
        scores.push(finite(score, "golden curve score")?);
    }
    let mut ranking = (0..candidates.len()).collect::<Vec<_>>();
    ranking.sort_by(|i, j| scores[*i].total_cmp(&scores[*j]).then(i.cmp(j)));
    Ok(format!("{{\"version\":1,\"algorithm\":\"supplied_golden_curve_distance_ranking\",\"target_id\":{},\"candidate_ids\":{},\"candidate_roles\":{},\"scores\":{},\"ranking_ids\":{},\"metric\":{metric},\"normalization\":{normalize},\"model_fitted\":false,\"target_excluded\":true,\"scope\":\"distance_matching_not_FPC_profiler_or_golden_curve_optimization\"}}",target.id,numbers(&candidates.iter().map(|c|c.id).collect::<Vec<_>>()),numbers(&candidates.iter().map(|c|c.role as f64).collect::<Vec<_>>()),numbers(&scores),numbers(&ranking.iter().map(|i|candidates[*i].id).collect::<Vec<_>>())))
}
