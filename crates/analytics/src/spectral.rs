//! Bounded spectral preprocessing and common-grid functional analysis.
//!
//! Fixed Savitzky–Golay settings, global polynomial detrending, orthonormal Haar,
//! direct DFT and trapezoid-grid functional PCA are explicit numerical methods.
//! They do not reproduce JMP's automatic bandwidth/basis selection, sparse-grid
//! models, spline models, complete wavelet library, or fitted peak models.
//! Curves normally occupy rows and equally spaced samples occupy columns.

const MAX_ROWS: usize = 10_000;
const MAX_CELLS: usize = 1_000_000;
const MAX_WORK: usize = 50_000_000;
type Matrix = Vec<Vec<f64>>;
type PolynomialQr = (Matrix, Matrix);

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data.first().map_or(0, Vec::len);
    if data.is_empty()
        || data.len() > MAX_ROWS
        || p == 0
        || p > 1024
        || data.len() > MAX_CELLS / p
        || args.len() > 16384
        || data
            .iter()
            .any(|r| r.len() != p || r.iter().any(|v| !v.is_finite()))
        || args.iter().any(|v| !v.is_finite())
    {
        return Err("spectral data must be finite rectangular rows within 10000 curves, 1024 columns and 1000000 cells; arguments finite and <=16384".into());
    }
    match operation {
        "spectral_snv" => snv(data, args),
        "spectral_msc" => msc(data, args),
        "savitzky_golay" => savgol(data, args),
        "curve_derivative" => derivative(data, args),
        "spectral_detrend" => detrend(data, args),
        "haar_transform" => haar(data, args),
        "fourier_transform" => fourier(data, args),
        "dynamic_time_warping" => dtw(data, args),
        "functional_pca" => fpca(data, args),
        "peak_find" => peaks(data, args),
        _ => Err(format!("unknown spectral operation: {operation}")),
    }
}
fn arguments(a: &[f64], min: usize, max: usize) -> Result<(), String> {
    if (min..=max).contains(&a.len()) {
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
fn finite(x: f64, name: &str) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.checked_mul(b)
        .and_then(|v| v.checked_mul(c))
        .is_none_or(|v| v > MAX_WORK)
    {
        Err("spectral computation exceeds 50000000 work units".into())
    } else {
        Ok(())
    }
}
fn vector(v: &[f64]) -> Result<String, String> {
    if v.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite spectral result".into());
    }
    Ok(format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    ))
}
fn matrix(v: &[Vec<f64>]) -> Result<String, String> {
    Ok(format!(
        "[{}]",
        v.iter()
            .map(|r| vector(r))
            .collect::<Result<Vec<_>, _>>()?
            .join(",")
    ))
}
fn output(name: &str, rows: &[Vec<f64>], extra: &str) -> Result<String, String> {
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"{name}\",\"rows\":{}{extra}}}",
        matrix(rows)?
    ))
}
fn center(v: &[f64]) -> Result<(f64, f64, Vec<f64>), String> {
    let origin = v[0];
    let mut offset = 0.;
    for (i, &v) in v.iter().enumerate() {
        let d = finite(v - origin, "origin-relative spectrum")?;
        offset = finite(offset + (d - offset) / (i + 1) as f64, "spectrum mean")?;
    }
    let centered = v
        .iter()
        .map(|v| finite((v - origin) - offset, "centered spectrum"))
        .collect::<Result<Vec<_>, _>>()?;
    Ok((origin, offset, centered))
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64, String> {
    finite(
        a.iter().zip(b).map(|(a, b)| a * b).sum(),
        "spectral inner product",
    )
}
fn snv(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 1)?;
    let ddof = integer(args.first().copied().unwrap_or(1.), 0, 1, "variance ddof")?;
    if data[0].len() < 2 {
        return Err("SNV needs at least two samples per spectrum".into());
    }
    let mut rows = Vec::new();
    let mut means = Vec::new();
    let mut sd = Vec::new();
    for row in data {
        let (o, m, c) = center(row)?;
        let max = c.iter().map(|v| v.abs()).fold(0., f64::max);
        if max == 0. {
            return Err("SNV is undefined for a constant spectrum".into());
        }
        let scaled =
            (c.iter().map(|v| (v / max).powi(2)).sum::<f64>() / (row.len() - ddof) as f64).sqrt();
        let scale = finite(max * scaled, "spectrum standard deviation")?;
        rows.push(c.iter().map(|v| (v / max) / scaled).collect());
        means.push(finite(o + m, "spectrum mean")?);
        sd.push(scale);
    }
    output(
        "per_spectrum_standard_normal_variate",
        &rows,
        &format!(
            ",\"ddof\":{ddof},\"means\":{},\"standard_deviations\":{}",
            vector(&means)?,
            vector(&sd)?
        ),
    )
}
fn msc(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let p = data[0].len();
    if p < 2 || !args.is_empty() && args.len() != p {
        return Err("MSC needs at least two samples and either no arguments or one explicit reference value per column".into());
    }
    let reference = if args.is_empty() {
        (0..p)
            .map(|j| {
                let values = data.iter().map(|r| r[j]).collect::<Vec<_>>();
                let (o, m, _) = center(&values)?;
                finite(o + m, "reference spectrum")
            })
            .collect::<Result<Vec<_>, _>>()?
    } else {
        args.to_vec()
    };
    let (ro, rm, rc) = center(&reference)?;
    let scale = rc.iter().map(|v| v.abs()).fold(0., f64::max);
    if scale == 0. {
        return Err("MSC reference spectrum must have positive variance".into());
    }
    let r = rc.iter().map(|v| v / scale).collect::<Vec<_>>();
    let denominator = dot(&r, &r)?;
    let mut rows = Vec::new();
    let mut slopes = Vec::new();
    let mut intercepts = Vec::new();
    for row in data {
        let (o, m, c) = center(row)?;
        let ymax = c.iter().map(|v| v.abs()).fold(0., f64::max);
        if ymax == 0. {
            return Err("MSC slope is zero for a constant spectrum".into());
        }
        let ys = c.iter().map(|v| v / ymax).collect::<Vec<_>>();
        let slope = finite((dot(&r, &ys)? / denominator) * (ymax / scale), "MSC slope")?;
        if slope == 0. {
            return Err("MSC reference regression has zero slope".into());
        }
        let corrected = c
            .iter()
            .map(|v| finite(ro + (rm + v / slope), "MSC corrected spectrum"))
            .collect::<Result<Vec<_>, _>>()?;
        rows.push(corrected);
        slopes.push(slope);
        intercepts.push(finite(
            (o - slope * ro) + (m - slope * rm),
            "MSC intercept",
        )?);
    }
    output("multiplicative_scatter_correction_reference_regression",&rows,&format!(",\"reference_spectrum\":{},\"reference_parameters\":{},\"reference_estimator\":\"{}\",\"slopes\":{},\"intercepts\":{}",vector(&reference)?,vector(&reference)?,if args.is_empty(){"mean_of_supplied_curves"}else{"supplied_saved_reference"},vector(&slopes)?,vector(&intercepts)?))
}
// QR polynomial basis in [-1,1], columns ascending powers, with two passes.
fn polynomial_qr(n: usize, degree: usize) -> Result<PolynomialQr, String> {
    let x = (0..n)
        .map(|i| 2. * i as f64 / (n - 1) as f64 - 1.)
        .collect::<Vec<_>>();
    let m = degree + 1;
    let mut q = Vec::<Vec<f64>>::new();
    let mut r = vec![vec![0.; m]; m];
    let columns = (0..m)
        .map(|j| x.iter().map(|x| x.powi(j as i32)).collect::<Vec<_>>())
        .collect::<Vec<_>>();
    for (j, mut v) in columns.into_iter().enumerate() {
        let original = dot(&v, &v)?.sqrt();
        for _ in 0..2 {
            for i in 0..j {
                let amount = dot(&q[i], &v)?;
                r[i][j] += amount;
                for (v, q) in v.iter_mut().zip(&q[i]) {
                    *v -= amount * q;
                }
            }
        }
        let norm = dot(&v, &v)?.sqrt();
        if norm <= original * 1e-12 {
            return Err("polynomial grid is numerically rank deficient".into());
        }
        r[j][j] = norm;
        for v in &mut v {
            *v /= norm;
        }
        q.push(v);
    }
    Ok((q, r))
}
fn polynomial_fit(y: &[f64], q: &[Vec<f64>], r: &[Vec<f64>]) -> Result<Vec<f64>, String> {
    let m = q.len();
    let mut b = q.iter().map(|q| dot(q, y)).collect::<Result<Vec<_>, _>>()?;
    for i in (0..m).rev() {
        b[i] = finite(
            (b[i] - (i + 1..m).map(|j| r[i][j] * b[j]).sum::<f64>()) / r[i][i],
            "polynomial coefficient",
        )?;
    }
    Ok(b)
}
fn savgol(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 1, 4)?;
    let n = data[0].len();
    let window = integer(args[0], 3, n.min(101), "SG window")?;
    if window & 1 == 0 {
        return Err("SG window must be odd".into());
    }
    let degree = integer(
        args.get(1).copied().unwrap_or(2.),
        0,
        8.min(window - 1),
        "SG degree",
    )?;
    let derivative = integer(
        args.get(2).copied().unwrap_or(0.),
        0,
        degree,
        "SG derivative",
    )?;
    let delta = args.get(3).copied().unwrap_or(1.);
    if delta <= 0. {
        return Err("sample spacing must be positive".into());
    }
    work(data.len(), n, window)?;
    let (q, r) = polynomial_qr(window, degree)?;
    let half = (window / 2) as f64;
    let denominator = finite(
        (half * delta).powi(derivative as i32),
        "SG derivative spacing power",
    )?;
    if denominator == 0. {
        return Err("SG derivative spacing underflow".into());
    }
    let mut weights = Vec::new();
    for position in 0..window {
        let t = (position as f64 - half) / half;
        let mut c = vec![0.; degree + 1];
        for (j, v) in c.iter_mut().enumerate().skip(derivative) {
            let factorial = (j - derivative + 1..=j).fold(1., |a, k| a * k as f64);
            *v = factorial * t.powi((j - derivative) as i32) / denominator;
        }
        let mut z = c;
        for i in 0..z.len() {
            z[i] = finite(
                (z[i] - (0..i).map(|j| r[j][i] * z[j]).sum::<f64>()) / r[i][i],
                "SG derivative weights",
            )?;
        }
        weights.push(
            (0..window)
                .map(|i| (0..z.len()).map(|j| q[j][i] * z[j]).sum::<f64>())
                .collect::<Vec<_>>(),
        );
    }
    let mut rows = Vec::new();
    for row in data {
        let mut filtered = Vec::with_capacity(n);
        for i in 0..n {
            let start = i.saturating_sub(window / 2).min(n - window);
            let origin = row[start];
            let relative = row[start..start + window]
                .iter()
                .map(|v| finite(v - origin, "SG window value"))
                .collect::<Result<Vec<_>, _>>()?;
            let value = dot(&weights[i - start], &relative)?;
            filtered.push(finite(
                if derivative == 0 {
                    origin + value
                } else {
                    value
                },
                "SG output",
            )?);
        }
        rows.push(filtered);
    }
    output("fixed_savitzky_golay_interp_edges",&rows,&format!(",\"window\":{window},\"polynomial_degree\":{degree},\"derivative\":{derivative},\"sample_spacing\":{delta},\"boundary_mode\":\"interp\""))
}
fn derivative(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 2)?;
    let order = integer(
        args.first().copied().unwrap_or(1.),
        1,
        2,
        "derivative order",
    )?;
    let delta = args.get(1).copied().unwrap_or(1.);
    let n = data[0].len();
    if n < 3 || delta <= 0. {
        return Err("second-order gradient stencils need >=3 samples and positive spacing".into());
    }
    let mut rows = data.to_vec();
    for _ in 0..order {
        rows = rows
            .iter()
            .map(|y| {
                let mut out = vec![0.; n];
                out[0] = finite(
                    (2. * (y[1] - y[0]) - 0.5 * (y[2] - y[0])) / delta,
                    "left derivative",
                )?;
                for i in 1..n - 1 {
                    out[i] = finite(
                        0.5 * ((y[i + 1] - y[i]) + (y[i] - y[i - 1])) / delta,
                        "central derivative",
                    )?;
                }
                out[n - 1] = finite(
                    (2. * (y[n - 1] - y[n - 2]) - 0.5 * (y[n - 1] - y[n - 3])) / delta,
                    "right derivative",
                )?;
                Ok(out)
            })
            .collect::<Result<Vec<_>, String>>()?;
    }
    output(
        "repeated_second_order_central_gradient",
        &rows,
        &format!(",\"order\":{order},\"sample_spacing\":{delta},\"boundary_accuracy_order\":2"),
    )
}
fn detrend(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 1)?;
    let n = data[0].len();
    if n < 2 {
        return Err("detrending needs >=2 samples".into());
    }
    let degree = integer(
        args.first().copied().unwrap_or(1.),
        0,
        8.min(n - 1),
        "baseline degree",
    )?;
    work(data.len(), n, degree + 1)?;
    let (q, r) = polynomial_qr(n, degree)?;
    let mut rows = Vec::new();
    let mut baselines = Vec::new();
    let mut coefficients = Vec::new();
    for y in data {
        let origin = y[0];
        let relative = y
            .iter()
            .map(|v| finite(v - origin, "detrend centered value"))
            .collect::<Result<Vec<_>, _>>()?;
        let mut b = polynomial_fit(&relative, &q, &r)?;
        let mut baseline = Vec::new();
        let mut residual = Vec::new();
        for (i, &observed) in relative.iter().enumerate() {
            let t = 2. * i as f64 / (n - 1) as f64 - 1.;
            let fitted = b.iter().rev().fold(0., |a, c| a * t + c);
            baseline.push(finite(origin + fitted, "polynomial baseline")?);
            residual.push(finite(observed - fitted, "detrended value")?);
        }
        b[0] = finite(origin + b[0], "baseline intercept")?;
        rows.push(residual);
        baselines.push(baseline);
        coefficients.push(b);
    }
    output("global_polynomial_baseline_detrend",&rows,&format!(",\"degree\":{degree},\"polynomial_grid\":\"equally_spaced_minus1_to_plus1\",\"baselines\":{},\"coefficients\":{},\"coefficient_order\":\"ascending_power\"",matrix(&baselines)?,matrix(&coefficients)?))
}
fn haar(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 2)?;
    let inverse = integer(args.first().copied().unwrap_or(0.), 0, 1, "inverse flag")?;
    let n = data[0].len();
    if n < 2 || !n.is_power_of_two() {
        return Err("Haar requires power-of-two spectrum length >=2".into());
    }
    let maximum = n.ilog2() as usize;
    let levels = integer(
        args.get(1).copied().unwrap_or(maximum as f64),
        1,
        maximum,
        "Haar levels",
    )?;
    let q = std::f64::consts::FRAC_1_SQRT_2;
    let mut rows = Vec::new();
    for row in data {
        let mut v = row.clone();
        if inverse == 0 {
            let mut length = n;
            for _ in 0..levels {
                let mut next = vec![0.; length];
                for i in 0..length / 2 {
                    next[i] = finite(v[2 * i] * q + v[2 * i + 1] * q, "Haar average")?;
                    next[length / 2 + i] = finite(v[2 * i] * q - v[2 * i + 1] * q, "Haar detail")?;
                }
                v[..length].copy_from_slice(&next);
                length /= 2;
            }
        } else {
            let mut length = n >> (levels - 1);
            for _ in 0..levels {
                let mut next = vec![0.; length];
                for i in 0..length / 2 {
                    next[2 * i] = finite(v[i] * q + v[length / 2 + i] * q, "inverse Haar value")?;
                    next[2 * i + 1] =
                        finite(v[i] * q - v[length / 2 + i] * q, "inverse Haar value")?;
                }
                v[..length].copy_from_slice(&next);
                length *= 2;
            }
        }
        rows.push(v);
    }
    output("orthonormal_haar_mallat_transform",&rows,&format!(",\"inverse\":{inverse},\"levels\":{levels},\"coefficient_layout\":\"approximation_then_coarse_to_fine_details\""))
}
fn fourier(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 1)?;
    let inverse = integer(args.first().copied().unwrap_or(0.), 0, 1, "inverse flag")?;
    let width = data[0].len();
    if inverse == 1 && width & 1 != 0 {
        return Err(
            "inverse Fourier input packs real coefficients followed by imaginary coefficients"
                .into(),
        );
    }
    let n = if inverse == 0 { width } else { width / 2 };
    if n > 512 {
        return Err("direct Fourier transform supports at most512 samples, keeping packed inverse models within1024 columns".into());
    }
    work(data.len(), n, n)?;
    let mut real = Vec::new();
    let mut imaginary = Vec::new();
    for row in data {
        let mut re = vec![0.; n];
        let mut im = vec![0.; n];
        for k in 0..n {
            for t in 0..n {
                let angle = 2. * std::f64::consts::PI * k as f64 * t as f64 / n as f64;
                let c = angle.cos();
                let s = angle.sin();
                if inverse == 0 {
                    re[k] += row[t] * c;
                    im[k] -= row[t] * s;
                } else {
                    re[k] += (row[t] * c - row[n + t] * s) / n as f64;
                    im[k] += (row[t] * s + row[n + t] * c) / n as f64;
                }
            }
            finite(re[k], "Fourier real value")?;
            finite(im[k], "Fourier imaginary value")?;
        }
        if inverse == 0 {
            re.extend(&im);
        }
        real.push(re);
        imaginary.push(im);
    }
    output("direct_discrete_fourier_transform",&real,&format!(",\"inverse\":{inverse},\"normalization\":\"forward_unscaled_inverse_divide_n\",\"complex_layout\":\"all_real_then_all_imaginary\",\"imaginary_rows\":{}",matrix(&imaginary)?))
}
fn dtw(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 1)?;
    if data[0].len() != 2 {
        return Err("DTW rows must be [value,series_id0or1] in within-series order".into());
    }
    let mut series = [Vec::new(), Vec::new()];
    for row in data {
        if row[1] != 0. && row[1] != 1. {
            return Err("DTW series labels must be0or1".into());
        }
        series[row[1] as usize].push(row[0]);
    }
    let n = series[0].len();
    let m = series[1].len();
    if n == 0 || m == 0 || n.checked_mul(m).is_none_or(|v| v > 4_000_000) {
        return Err("DTW needs both series and at most4000000 cell pairs".into());
    }
    let window = integer(
        args.first().copied().unwrap_or(n.max(m) as f64),
        n.abs_diff(m),
        n.max(m),
        "DTW window radius",
    )?;
    let origin = series[0][0];
    let scale = series
        .iter()
        .flatten()
        .map(|v| finite(v - origin, "DTW origin-relative value"))
        .collect::<Result<Vec<_>, _>>()?
        .iter()
        .map(|v| v.abs())
        .fold(0., f64::max);
    let scale = if scale == 0. { 1. } else { scale };
    let a = series[0]
        .iter()
        .map(|v| (v - origin) / scale)
        .collect::<Vec<_>>();
    let b = series[1]
        .iter()
        .map(|v| (v - origin) / scale)
        .collect::<Vec<_>>();
    let mut costs = vec![f64::INFINITY; (n + 1) * (m + 1)];
    let mut previous = vec![0u8; costs.len()];
    costs[0] = 0.;
    for i in 1..=n {
        for j in 1.max(i.saturating_sub(window))..=m.min(i + window) {
            let idx = i * (m + 1) + j;
            let mut best = costs[(i - 1) * (m + 1) + j - 1];
            let mut direction = 0;
            if costs[(i - 1) * (m + 1) + j] < best {
                best = costs[(i - 1) * (m + 1) + j];
                direction = 1;
            }
            if costs[i * (m + 1) + j - 1] < best {
                best = costs[i * (m + 1) + j - 1];
                direction = 2;
            }
            costs[idx] = best + (a[i - 1] - b[j - 1]).powi(2);
            previous[idx] = direction;
        }
    }
    let final_cost = finite(costs[n * (m + 1) + m], "DTW endpoint cost")?;
    let distance = finite(final_cost.sqrt() * scale, "DTW distance")?;
    let mut path = Vec::new();
    let (mut i, mut j) = (n, m);
    while i > 0 && j > 0 {
        path.push(vec![(i - 1) as f64, (j - 1) as f64]);
        match previous[i * (m + 1) + j] {
            0 => {
                i -= 1;
                j -= 1;
            }
            1 => i -= 1,
            _ => j -= 1,
        }
    }
    if i != 0 || j != 0 {
        return Err("DTW path did not reach both origins".into());
    }
    path.reverse();
    Ok(format!("{{\"version\":1,\"algorithm\":\"sakoe_chiba_dtw_squared_cost_sqrt_distance\",\"distance\":{distance},\"cost_scale\":{scale},\"normalized_squared_cost\":{final_cost},\"window_radius\":{window},\"path\":{}}}",matrix(&path)?))
}
fn eigen(mut a: Vec<Vec<f64>>) -> Result<(Vec<f64>, Vec<Vec<f64>>), String> {
    let n = a.len();
    let mut v = vec![vec![0.; n]; n];
    for (i, r) in v.iter_mut().enumerate() {
        r[i] = 1.;
    }
    let max = (0..n).map(|j| a[j][j].abs()).fold(0., f64::max);
    if max == 0. {
        return Ok((vec![0.; n], v));
    }
    let threshold = max * 2e-14;
    let mut done = false;
    for _ in 0..80 {
        let mut largest = 0_f64;
        for i in 0..n {
            for j in i + 1..n {
                let off = a[i][j];
                largest = largest.max(off.abs());
                if off.abs() <= threshold {
                    continue;
                }
                let tau = (a[j][j] - a[i][i]) / (2. * off);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / (-tau + tau.hypot(1.))
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
                    .map(|(k, r)| (k, c * r[i] - s * r[j], s * r[i] + c * r[j]))
                    .collect::<Vec<_>>();
                for (k, u, w) in updates {
                    a[k][i] = u;
                    a[i][k] = u;
                    a[k][j] = w;
                    a[j][k] = w;
                }
                for row in &mut v {
                    let u = row[i];
                    let w = row[j];
                    row[i] = c * u - s * w;
                    row[j] = s * u + c * w;
                }
            }
        }
        if largest <= threshold {
            done = true;
            break;
        }
    }
    if !done {
        return Err("functional PCA eigensystem failed to converge".into());
    }
    let mut order = (0..n).collect::<Vec<_>>();
    order.sort_by(|&i, &j| a[j][j].total_cmp(&a[i][i]));
    Ok((
        order.iter().map(|&i| a[i][i]).collect(),
        (0..n)
            .map(|j| order.iter().map(|&i| v[j][i]).collect())
            .collect(),
    ))
}
fn fpca(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 1, 2)?;
    let n = data.len();
    let p = data[0].len();
    if n < 2 || p < 2 || n.min(p) > 128 {
        return Err(
            "grid functional PCA requires>=2curves,>=2samples and min(curves,samples)<=128".into(),
        );
    }
    let count = integer(args[0], 1, (n - 1).min(p), "functional components")?;
    let delta = args.get(1).copied().unwrap_or(1.);
    if delta <= 0. {
        return Err("functional grid spacing must be positive".into());
    }
    work(n, p, n.min(p))?;
    let weights = (0..p)
        .map(|j| {
            if j == 0 || j == p - 1 {
                delta / 2.
            } else {
                delta
            }
        })
        .collect::<Vec<_>>();
    if weights.contains(&0.) {
        return Err("quadrature weights underflow".into());
    }
    let roots = weights.iter().map(|v| v.sqrt()).collect::<Vec<_>>();
    let mut origins = Vec::new();
    let mut offsets = Vec::new();
    let mut means = Vec::new();
    let mut a = vec![vec![0.; p]; n];
    for j in 0..p {
        let values = data.iter().map(|r| r[j]).collect::<Vec<_>>();
        let (o, m, c) = center(&values)?;
        origins.push(o);
        offsets.push(m);
        means.push(finite(o + m, "functional mean")?);
        for i in 0..n {
            a[i][j] = finite(c[i] * roots[j], "weighted centered curve")?;
        }
    }
    let magnitude = a.iter().flatten().map(|v| v.abs()).fold(0., f64::max);
    if magnitude == 0. {
        return Err("functional PCA is undefined without between-curve variation".into());
    }
    for row in &mut a {
        for v in row {
            *v /= magnitude;
        }
    }
    let dual = n <= p;
    let size = n.min(p);
    let mut gram = vec![vec![0.; size]; size];
    for i in 0..size {
        for j in i..size {
            let value = if dual {
                dot(&a[i], &a[j])?
            } else {
                (0..n).map(|k| a[k][i] * a[k][j]).sum::<f64>()
            } / (n - 1) as f64;
            gram[i][j] = value;
            gram[j][i] = value;
        }
    }
    let (eigenvalues, vectors) = eigen(gram)?;
    let threshold = eigenvalues[0] * 1e-12;
    let rank = eigenvalues.iter().filter(|&&v| v > threshold).count();
    if count > rank {
        return Err("requested functional components exceed numerical covariance rank".into());
    }
    let total = eigenvalues.iter().filter(|&&v| v > 0.).sum::<f64>();
    let mut values = Vec::new();
    let mut proportions = Vec::new();
    let mut functions = Vec::new();
    let mut loading = Vec::new();
    for k in 0..count {
        let mut axis = if dual {
            let denominator = ((n - 1) as f64 * eigenvalues[k]).sqrt();
            (0..p)
                .map(|j| (0..n).map(|i| a[i][j] * vectors[i][k]).sum::<f64>() / denominator)
                .collect::<Vec<_>>()
        } else {
            (0..p).map(|j| vectors[j][k]).collect()
        };
        let pivot = (0..p)
            .max_by(|&i, &j| axis[i].abs().total_cmp(&axis[j].abs()))
            .unwrap();
        if axis[pivot] < 0. {
            for v in &mut axis {
                *v = -*v;
            }
        }
        values.push(finite(
            (magnitude * eigenvalues[k].sqrt()).powi(2),
            "functional eigenvalue",
        )?);
        proportions.push(eigenvalues[k] / total);
        functions.push(
            axis.iter()
                .zip(&roots)
                .map(|(v, w)| finite(v / w, "functional eigenfunction"))
                .collect::<Result<Vec<_>, _>>()?,
        );
        loading.push(axis);
    }
    let scores = a
        .iter()
        .map(|r| {
            loading
                .iter()
                .map(|v| finite(dot(r, v)? * magnitude, "functional score"))
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    let reconstructed = (0..n)
        .map(|i| {
            (0..p)
                .map(|j| {
                    finite(
                        origins[j]
                            + (offsets[j]
                                + (0..count)
                                    .map(|k| scores[i][k] * functions[k][j])
                                    .sum::<f64>()),
                        "functional reconstruction",
                    )
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!("{{\"version\":1,\"algorithm\":\"common_equally_spaced_grid_trapezoid_l2_fpca\",\"sample_spacing\":{delta},\"components\":{count},\"numerical_rank\":{rank},\"quadrature_weights\":{},\"mean_curve\":{},\"mean_origins\":{},\"mean_offsets\":{},\"eigenvalues\":{},\"explained_variance_proportions\":{},\"eigenfunctions\":{},\"scores\":{},\"reconstructed_rows\":{},\"scope\":\"raw_common_grid_no_basis_smoothing_or_sparse_curve_estimation\"}}",vector(&weights)?,vector(&means)?,vector(&origins)?,vector(&offsets)?,vector(&values)?,vector(&proportions)?,matrix(&functions)?,matrix(&scores)?,matrix(&reconstructed)?))
}
fn peaks(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    arguments(args, 0, 2)?;
    let n = data[0].len();
    if n < 3 {
        return Err("peak finding needs at least3samples".into());
    }
    let prominence = args.first().copied().unwrap_or(0.);
    if prominence < 0. {
        return Err("minimum prominence must be nonnegative".into());
    }
    let distance = integer(
        args.get(1).copied().unwrap_or(1.),
        1,
        n,
        "minimum peak distance",
    )?;
    let mut indices = Vec::new();
    let mut heights = Vec::new();
    let mut prominences = Vec::new();
    let mut lefts = Vec::new();
    let mut rights = Vec::new();
    let mut used = 0usize;
    for y in data {
        let mut candidates = Vec::new();
        let mut i = 1;
        while i + 1 < n {
            if y[i] > y[i - 1] {
                let start = i;
                while i + 1 < n && y[i + 1] == y[start] {
                    i += 1;
                }
                if i + 1 < n && y[i] > y[i + 1] {
                    candidates.push((start + i) / 2);
                }
            }
            i += 1;
        }
        candidates.sort_by(|&i, &j| y[j].total_cmp(&y[i]).then(i.cmp(&j)));
        let mut selected = Vec::<usize>::new();
        let mut blocked = vec![false; n];
        for p in candidates {
            if !blocked[p] {
                selected.push(p);
                let start = p.saturating_sub(distance - 1);
                let end = (p + distance).min(n);
                blocked[start..end].fill(true);
            }
        }
        selected.sort_unstable();
        let mut ix = Vec::new();
        let mut hs = Vec::new();
        let mut ps = Vec::new();
        let mut ls = Vec::new();
        let mut rs = Vec::new();
        for p in selected {
            let (mut left, mut right) = (p, p);
            let (mut lb, mut rb) = (p, p);
            let (mut lv, mut rv) = (y[p], y[p]);
            while left > 0 && y[left - 1] <= y[p] {
                left -= 1;
                used += 1;
                if y[left] < lv {
                    lv = y[left];
                    lb = left;
                }
            }
            while right + 1 < n && y[right + 1] <= y[p] {
                right += 1;
                used += 1;
                if y[right] < rv {
                    rv = y[right];
                    rb = right;
                }
            }
            if used > MAX_WORK {
                return Err("peak prominence search exceeds50000000steps".into());
            }
            let pr = finite(y[p] - lv.max(rv), "peak prominence")?;
            if pr >= prominence {
                ix.push(p as f64);
                hs.push(y[p]);
                ps.push(pr);
                ls.push(lb as f64);
                rs.push(rb as f64);
            }
        }
        indices.push(ix);
        heights.push(hs);
        prominences.push(ps);
        lefts.push(ls);
        rights.push(rs);
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"local_maxima_global_topographic_prominence\",\"minimum_prominence\":{prominence},\"minimum_distance\":{distance},\"plateau_location\":\"lower_middle_index\",\"distance_ties\":\"lowest_index_first\",\"peak_indices\":{},\"peak_heights\":{},\"prominences\":{},\"left_bases\":{},\"right_bases\":{},\"scope\":\"sampled_local_peaks_not_fitted_FDE_peak_models\"}}",matrix(&indices)?,matrix(&heights)?,matrix(&prominences)?,matrix(&lefts)?,matrix(&rights)?))
}
