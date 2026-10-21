//! Convex binary C-SVC, check-loss regression, and fixed-scale Huber M-regression.
//! All preprocessing is fitted on training rows only. Saved models retain origins
//! and mean offsets separately to avoid losing small variations at large origins.
//! SVM uses feasible pair-coordinate SMO with a reported dual KKT gap. Quantile
//! ADMM reports both residual criteria; iteration limits never imply convergence.

const MAX_WORK: usize = 200_000_000;
pub fn run(op: &str, d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let w = d.first().map_or(0, Vec::len);
    if w == 0
        || w > 33
        || d.len() > 10000
        || d.len().saturating_mul(w) > 1_000_000
        || d.iter()
            .any(|r| r.len() != w || r.iter().any(|v| !v.is_finite()))
        || a.len() > 16384
        || a.iter().any(|v| !v.is_finite())
    {
        return Err(
            "robust learning requires finite rectangular data, <=10000 rows and 33 columns".into(),
        );
    }
    let out = match op {
        "svm_fit" => svm(d, a),
        "svm_predict" => svm_predict(d, a),
        "quantile_regression" => quantile(d, a),
        "huber_regression" => huber(d, a),
        "robust_predict" => predict(d, a),
        _ => Err("unknown robust learning operation".into()),
    }?;
    if out
        .split([':', ',', '[', ']', '}', ' '])
        .any(|v| matches!(v, "NaN" | "inf" | "-inf"))
    {
        Err("robust learning result exceeds finite numerical range".into())
    } else {
        Ok(out)
    }
}
fn nums(x: &[f64]) -> String {
    format!(
        "[{}]",
        x.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn mat(x: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        x.iter().map(|r| nums(r)).collect::<Vec<_>>().join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    a.iter().fold(0_f64, |v, x| v.hypot(*x))
}
fn int(x: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if x.fract() != 0. || x < lo as f64 || x > hi as f64 {
        Err(format!("expected integer in [{lo},{hi}]"))
    } else {
        Ok(x as usize)
    }
}
fn tolerance(x: f64) -> Result<f64, String> {
    if (1e-12..=1e-2).contains(&x) {
        Ok(x)
    } else {
        Err("tolerance must be in [1e-12,1e-2]".into())
    }
}
fn inv(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let scales = (0..n).map(|i| a[i][i].sqrt()).collect::<Vec<_>>();
    if scales.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err("singular regression information".into());
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            let v =
                a[i][j] / scales[i] / scales[j] - (0..j).map(|k| l[i][k] * l[j][k]).sum::<f64>();
            if i == j {
                if v <= 1e-12 || !v.is_finite() {
                    return Err("rank deficient regression information".into());
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut out = vec![vec![0.; n]; n];
    for col in 0..n {
        let mut b = vec![0.; n];
        for i in 0..n {
            b[i] = ((if i == col { 1. } else { 0. }) - (0..i).map(|j| l[i][j] * b[j]).sum::<f64>())
                / l[i][i];
        }
        let mut x = vec![0.; n];
        for i in (0..n).rev() {
            x[i] = (b[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            out[i][col] = x[i] / scales[i] / scales[col];
        }
    }
    Ok(out)
}
fn product(a: &[Vec<f64>], b: &[Vec<f64>]) -> Vec<Vec<f64>> {
    (0..a.len())
        .map(|i| {
            (0..b[0].len())
                .map(|j| (0..b.len()).map(|k| a[i][k] * b[k][j]).sum())
                .collect()
        })
        .collect()
}
fn info(x: &[Vec<f64>], weights: &[f64]) -> Vec<Vec<f64>> {
    let q = x[0].len();
    (0..q)
        .map(|i| {
            (0..q)
                .map(|j| x.iter().zip(weights).map(|(r, w)| r[i] * r[j] * w).sum())
                .collect()
        })
        .collect()
}
fn wls(x: &[Vec<f64>], y: &[f64], w: &[f64]) -> Result<Vec<f64>, String> {
    let inverse = inv(&info(x, w))?;
    let rhs = (0..x[0].len())
        .map(|i| {
            x.iter()
                .zip(y)
                .zip(w)
                .map(|((r, y), w)| r[i] * y * w)
                .sum::<f64>()
        })
        .collect::<Vec<_>>();
    Ok(inverse.iter().map(|r| dot(r, &rhs)).collect())
}
#[derive(Clone)]
struct Standard {
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
}
impl Standard {
    fn fit(d: &[Vec<f64>], p: usize) -> Result<Self, String> {
        let origins = d[0][..p].to_vec();
        let mut offsets = vec![];
        let mut scales = vec![];
        for c in 0..p {
            let off = d
                .iter()
                .map(|r| (r[c] - origins[c]) / d.len() as f64)
                .sum::<f64>();
            let sd = d.iter().fold(0_f64, |v, r| {
                v.hypot(((r[c] - origins[c]) - off) / (d.len() as f64).sqrt())
            });
            if !off.is_finite() || !sd.is_finite() {
                return Err("standardization exceeds numeric range".into());
            }
            offsets.push(off);
            scales.push(if sd > 0. { sd } else { 1. });
        }
        Ok(Self {
            origins,
            offsets,
            scales,
        })
    }
    fn row(&self, r: &[f64]) -> Result<Vec<f64>, String> {
        let v = (0..self.origins.len())
            .map(|i| ((r[i] - self.origins[i]) - self.offsets[i]) / self.scales[i])
            .collect::<Vec<_>>();
        if v.iter().any(|v| !v.is_finite()) {
            Err("standardized features exceed numeric range".into())
        } else {
            Ok(v)
        }
    }
}
struct Regression {
    standard: Standard,
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
    origin: f64,
    offset: f64,
    scale: f64,
}
impl Regression {
    fn new(d: &[Vec<f64>]) -> Result<Self, String> {
        let p = d[0].len() - 1;
        if d.len() <= p + 1 {
            return Err("regression requires positive residual degrees of freedom".into());
        }
        let standard = Standard::fit(d, p)?;
        let x = d
            .iter()
            .map(|r| {
                let mut v = vec![1.];
                v.extend(standard.row(r)?);
                Ok(v)
            })
            .collect::<Result<Vec<_>, String>>()?;
        let origin = d[0][p];
        let offset = d
            .iter()
            .map(|r| (r[p] - origin) / d.len() as f64)
            .sum::<f64>();
        let scale = d.iter().fold(0_f64, |v, r| {
            v.hypot(((r[p] - origin) - offset) / (d.len() as f64).sqrt())
        });
        if !offset.is_finite() || !scale.is_finite() {
            return Err("response normalization exceeds range".into());
        }
        let scale = if scale > 0. { scale } else { 1. };
        let y = d
            .iter()
            .map(|r| ((r[p] - origin) - offset) / scale)
            .collect();
        Ok(Self {
            standard,
            x,
            y,
            origin,
            offset,
            scale,
        })
    }
    fn saved(&self, b: &[f64]) -> Vec<f64> {
        let mut a = vec![
            self.standard.origins.len() as f64,
            self.origin,
            self.offset,
            self.scale,
            b[0],
        ];
        a.extend(&self.standard.origins);
        a.extend(&self.standard.offsets);
        a.extend(&self.standard.scales);
        a.extend(&b[1..]);
        a
    }
    fn result(&self, b: &[f64]) -> (f64, Vec<f64>, Vec<f64>, Vec<f64>) {
        let coef = b[1..]
            .iter()
            .zip(&self.standard.scales)
            .map(|(b, s)| b * self.scale / s)
            .collect::<Vec<_>>();
        let intercept = self.origin + (self.offset + self.scale * b[0])
            - coef
                .iter()
                .enumerate()
                .map(|(i, b)| b * (self.standard.origins[i] + self.standard.offsets[i]))
                .sum::<f64>();
        let fitted = self
            .x
            .iter()
            .map(|r| self.origin + (self.offset + self.scale * dot(r, b)))
            .collect();
        let residual = self
            .x
            .iter()
            .zip(&self.y)
            .map(|(r, y)| (y - dot(r, b)) * self.scale)
            .collect();
        (intercept, coef, fitted, residual)
    }
}
fn predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() < 5 {
        return Err("portable regression parameters missing".into());
    }
    let p = int(a[0], 0, 32)?;
    if (p != d[0].len() && !(p == 0 && d[0].len() == 1))
        || a.len() != 5 + 4 * p
        || a[3] <= 0.
        || a[5 + 2 * p..5 + 3 * p].iter().any(|x| *x <= 0.)
    {
        return Err("invalid portable regression dimensions or scale".into());
    }
    let s = Standard {
        origins: a[5..5 + p].to_vec(),
        offsets: a[5 + p..5 + 2 * p].to_vec(),
        scales: a[5 + 2 * p..5 + 3 * p].to_vec(),
    };
    let fit = d
        .iter()
        .map(|r| {
            s.row(r)
                .map(|v| a[1] + (a[2] + a[3] * (a[4] + dot(&v, &a[5 + 3 * p..]))))
        })
        .collect::<Result<Vec<_>, _>>()?;
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"robust_linear_predict\",\"predictions\":{}}}",
        nums(&fit)
    ))
}
fn quantile(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() > 4 {
        return Err("quantile arguments are quantile,rho,max_iterations,tolerance".into());
    }
    let tau = *a.first().unwrap_or(&0.5);
    let rho = *a.get(1).unwrap_or(&1.);
    if tau <= 0. || tau >= 1. || !(1e-4..=1e4).contains(&rho) {
        return Err("quantile in (0,1),rho in [1e-4,1e4] required".into());
    }
    let maxiter = int(*a.get(2).unwrap_or(&10000.), 1, 50000)?;
    let tol = tolerance(*a.get(3).unwrap_or(&1e-7))?;
    let r = Regression::new(d)?;
    let n = d.len();
    let q = r.x[0].len();
    if n.saturating_mul(q).saturating_mul(maxiter) > MAX_WORK {
        return Err("quantile work limit exceeded".into());
    }
    let inverse = inv(&info(&r.x, &vec![1.; n]))?;
    let mut beta = wls(&r.x, &r.y, &vec![1.; n])?;
    let mut residual =
        r.x.iter()
            .zip(&r.y)
            .map(|(x, y)| y - dot(x, &beta))
            .collect::<Vec<_>>();
    let mut u = vec![0.; n];
    let mut primal = 0.;
    let mut dual = 0.;
    let mut pt = 0.;
    let mut dt = 0.;
    let mut converged = false;
    let mut iterations = 0;
    for it in 0..maxiter {
        iterations = it + 1;
        let rhs = (0..q)
            .map(|c| {
                r.x.iter()
                    .enumerate()
                    .map(|(i, x)| x[c] * (r.y[i] - residual[i] - u[i]))
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        beta = inverse.iter().map(|row| dot(row, &rhs)).collect();
        let xb = r.x.iter().map(|x| dot(x, &beta)).collect::<Vec<_>>();
        let old = residual.clone();
        for i in 0..n {
            let v = r.y[i] - xb[i] - u[i];
            residual[i] = if v > tau / rho {
                v - tau / rho
            } else if v < -(1. - tau) / rho {
                v + (1. - tau) / rho
            } else {
                0.
            };
            u[i] += xb[i] + residual[i] - r.y[i];
        }
        primal = norm(
            &(0..n)
                .map(|i| xb[i] + residual[i] - r.y[i])
                .collect::<Vec<_>>(),
        );
        let dv = (0..q)
            .map(|c| {
                rho * r
                    .x
                    .iter()
                    .enumerate()
                    .map(|(i, x)| x[c] * (residual[i] - old[i]))
                    .sum::<f64>()
            })
            .collect::<Vec<_>>();
        dual = norm(&dv);
        pt = tol * ((n as f64).sqrt() + norm(&xb).max(norm(&residual)).max(norm(&r.y)));
        let xu = (0..q)
            .map(|c| rho * r.x.iter().zip(&u).map(|(x, u)| x[c] * u).sum::<f64>())
            .collect::<Vec<_>>();
        dt = tol * ((q as f64).sqrt() + norm(&xu));
        if primal <= pt && dual <= dt {
            converged = true;
            break;
        }
    }
    let (intercept, coef, fit, res) = r.result(&beta);
    let objective = res
        .iter()
        .map(|e| if *e >= 0. { tau * e } else { (tau - 1.) * e })
        .sum::<f64>();
    Ok(format!("{{\"version\":1,\"algorithm\":\"quantile_regression_admm\",\"quantile\":{},\"rho\":{},\"intercept\":{},\"coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"objective\":{},\"primal_residual\":{},\"dual_residual\":{},\"primal_tolerance\":{},\"dual_tolerance\":{},\"convergence_units\":\"standardized_response_and_design\",\"converged\":{},\"iterations\":{},\"predict_parameters\":{}}}",tau,rho,intercept,nums(&coef),nums(&fit),nums(&res),objective,primal,dual,pt,dt,converged,iterations,nums(&r.saved(&beta))))
}
fn median(x: &[f64]) -> f64 {
    let mut a = x.to_vec();
    a.sort_by(f64::total_cmp);
    if a.len() & 1 == 0 {
        a[a.len() / 2 - 1] / 2. + a[a.len() / 2] / 2.
    } else {
        a[a.len() / 2]
    }
}
fn huber(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() > 4 {
        return Err("Huber arguments are delta,scale,max_iterations,tolerance".into());
    }
    let delta = *a.first().unwrap_or(&1.345);
    let requested_scale = *a.get(1).unwrap_or(&0.);
    if delta <= 0. || delta > 100. || requested_scale < 0. {
        return Err("Huber delta in (0,100], scale>=0 required".into());
    }
    let maxiter = int(*a.get(2).unwrap_or(&200.), 1, 5000)?;
    let tol = tolerance(*a.get(3).unwrap_or(&1e-8))?;
    let r = Regression::new(d)?;
    let n = d.len();
    let q = r.x[0].len();
    if n.saturating_mul(q * q).saturating_mul(maxiter) > MAX_WORK {
        return Err("Huber work limit exceeded".into());
    }
    let mut beta = wls(&r.x, &r.y, &vec![1.; n])?;
    let initial =
        r.x.iter()
            .zip(&r.y)
            .map(|(x, y)| y - dot(x, &beta))
            .collect::<Vec<_>>();
    let mid = median(&initial);
    let estimated =
        median(&initial.iter().map(|e| (e - mid).abs()).collect::<Vec<_>>()) / 0.6744897501960817;
    let mut scale = if requested_scale > 0. {
        requested_scale / r.scale
    } else {
        estimated
    };
    let exact = norm(&initial) <= 64. * f64::EPSILON * norm(&r.y).max(1.);
    if scale == 0. {
        if exact {
            scale = 1.;
        } else {
            return Err("OLS residual MAD is zero; specify a positive Huber scale".into());
        }
    }
    if !scale.is_finite() || scale <= 0. {
        return Err("Huber scale is outside finite range".into());
    }
    let mut weights = vec![1.; n];
    let mut converged = exact;
    let mut iterations = 0;
    if !exact {
        for it in 0..maxiter {
            iterations = it + 1;
            for (i, weight) in weights.iter_mut().enumerate() {
                let e = (r.y[i] - dot(&r.x[i], &beta)) / scale;
                *weight = if e.abs() <= delta {
                    1.
                } else {
                    delta / e.abs()
                };
            }
            let next = wls(&r.x, &r.y, &weights)?;
            let change = norm(
                &next
                    .iter()
                    .zip(&beta)
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>(),
            );
            beta = next;
            if change <= tol * (1. + norm(&beta)) {
                converged = true;
                break;
            }
        }
    }
    let errors =
        r.x.iter()
            .zip(&r.y)
            .map(|(x, y)| (y - dot(x, &beta)) / scale)
            .collect::<Vec<_>>();
    for (w, e) in weights.iter_mut().zip(&errors) {
        *w = if e.abs() <= delta {
            1.
        } else {
            delta / e.abs()
        };
    }
    let psi = errors
        .iter()
        .map(|e| e.clamp(-delta, delta))
        .collect::<Vec<_>>();
    let objective = errors
        .iter()
        .map(|e| {
            if e.abs() <= delta {
                0.5 * e * e
            } else {
                delta * (e.abs() - delta / 2.)
            }
        })
        .sum::<f64>();
    let score = (0..q)
        .map(|c| r.x.iter().zip(&psi).map(|(x, p)| x[c] * p).sum::<f64>())
        .collect::<Vec<_>>();
    let bread = info(
        &r.x,
        &errors
            .iter()
            .map(|e| if e.abs() <= delta { 1. } else { 0. })
            .collect::<Vec<_>>(),
    );
    // Restore score units before squaring. Squaring psi first can underflow
    // when a harmlessly large user scale leaves the estimator equal to OLS.
    let clipped_residual =
        r.x.iter()
            .zip(&r.y)
            .map(|(x, y)| {
                let e = y - dot(x, &beta);
                if (e / scale).abs() <= delta {
                    e
                } else {
                    e.signum() * (delta * scale)
                }
            })
            .collect::<Vec<_>>();
    let meat = info(
        &r.x,
        &clipped_residual.iter().map(|e| e * e).collect::<Vec<_>>(),
    );
    let covariance = inv(&bread).ok().map(|v| product(&product(&v, &meat), &v));
    let (cov, se) = if let Some(c) = covariance {
        let se = (1..q)
            .map(|i| c[i][i].max(0.).sqrt() * r.scale / r.standard.scales[i - 1])
            .collect::<Vec<_>>();
        (mat(&c), nums(&se))
    } else {
        ("null".into(), "null".into())
    };
    let (intercept, coef, fit, res) = r.result(&beta);
    Ok(format!("{{\"version\":1,\"algorithm\":\"huber_fixed_scale_irls\",\"delta\":{},\"robust_scale\":{},\"scale_method\":\"{}\",\"intercept\":{},\"coefficients\":{},\"fitted_values\":{},\"residuals\":{},\"weights\":{},\"objective\":{},\"score_norm\":{},\"standardized_covariance\":{},\"coefficient_standard_errors\":{},\"covariance_method\":\"fixed_scale_HC0_Huber_score_sandwich\",\"converged\":{},\"iterations\":{},\"predict_parameters\":{}}}",delta,if exact && requested_scale==0.{0.}else{scale*r.scale},if requested_scale>0.{"user_fixed"}else if exact{"exact_fit_zero_scale"}else{"OLS_residual_MAD_fixed"},intercept,nums(&coef),nums(&fit),nums(&res),nums(&weights),objective,norm(&score),cov,se,converged,iterations,nums(&r.saved(&beta))))
}
fn kernel(a: &[f64], b: &[f64], kind: usize, gamma: f64) -> f64 {
    if kind == 0 {
        dot(a, b)
    } else {
        (-gamma * a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum::<f64>()).exp()
    }
}
fn working(labels: &[f64], alpha: &[f64], f: &[f64], cost: f64) -> (usize, usize, f64, f64) {
    let eps = 32. * f64::EPSILON * cost;
    let mut up = (0, f64::NEG_INFINITY);
    let mut low = (0, f64::INFINITY);
    for i in 0..labels.len() {
        let t = labels[i] - f[i];
        if ((labels[i] > 0. && alpha[i] < cost - eps) || (labels[i] < 0. && alpha[i] > eps))
            && t > up.1
        {
            up = (i, t);
        }
        if ((labels[i] > 0. && alpha[i] > eps) || (labels[i] < 0. && alpha[i] < cost - eps))
            && t < low.1
        {
            low = (i, t);
        }
    }
    (up.0, low.0, up.1, low.1)
}
fn svm(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.is_empty() || a.len() > 6 || d[0].len() < 2 {
        return Err(
            "SVM requires feature/label rows and n_train,kernel,C,gamma,tolerance,max_iterations"
                .into(),
        );
    }
    let n = int(a[0], 2, d.len().min(512))?;
    let p = d[0].len() - 1;
    if p > 16 {
        return Err("SVM supports <=16 features".into());
    }
    let kind = int(*a.get(1).unwrap_or(&1.), 0, 1)?;
    let cost = *a.get(2).unwrap_or(&1.);
    let gamma = *a.get(3).unwrap_or(&(1. / p as f64));
    let tol = tolerance(*a.get(4).unwrap_or(&1e-6))?;
    let maxiter = int(*a.get(5).unwrap_or(&100000.), 1, 200000)?;
    if !(1e-8..=1e6).contains(&cost)
        || gamma <= 0.
        || gamma > 1e6
        || n.saturating_mul(maxiter) > MAX_WORK
        || n.saturating_mul(n).saturating_mul(p) > MAX_WORK
    {
        return Err("invalid SVM cost/gamma or work limit".into());
    }
    let labels = d[..n]
        .iter()
        .map(|r| match r[p] {
            0. => Ok(-1.),
            1. => Ok(1.),
            _ => Err("binary SVM training labels must be 0/1".to_string()),
        })
        .collect::<Result<Vec<_>, _>>()?;
    if !labels.contains(&-1.) || !labels.contains(&1.) {
        return Err("SVM training requires both classes".into());
    }
    let standard = Standard::fit(&d[..n], p)?;
    let x = d
        .iter()
        .map(|r| standard.row(r))
        .collect::<Result<Vec<_>, _>>()?;
    let k = (0..n)
        .map(|i| {
            (0..n)
                .map(|j| kernel(&x[i], &x[j], kind, gamma))
                .collect::<Vec<_>>()
        })
        .collect::<Vec<_>>();
    if k.iter().flatten().any(|v| !v.is_finite()) {
        return Err("SVM kernel exceeds numeric range".into());
    }
    let mut alpha = vec![0.; n];
    let mut fitted = vec![0.; n];
    let mut iterations = 0;
    for it in 0..maxiter {
        let (i, j, up, low) = working(&labels, &alpha, &fitted, cost);
        let gap = up - low;
        if gap <= tol {
            break;
        }
        let bound_i = if labels[i] > 0. {
            cost - alpha[i]
        } else {
            alpha[i]
        };
        let bound_j = if labels[j] > 0. {
            alpha[j]
        } else {
            cost - alpha[j]
        };
        let bound = bound_i.min(bound_j);
        let curvature = (k[i][i] + k[j][j] - 2. * k[i][j]).max(0.);
        let step = if curvature > 0. {
            (gap / curvature).min(bound)
        } else {
            bound
        };
        if step <= 0. || !step.is_finite() {
            break;
        }
        alpha[i] = (alpha[i] + labels[i] * step).clamp(0., cost);
        alpha[j] = (alpha[j] - labels[j] * step).clamp(0., cost);
        for t in 0..n {
            fitted[t] += step * (k[t][i] - k[t][j]);
        }
        iterations = it + 1;
    }
    // Recompute from dual variables, so serialized predictions and diagnostics use
    // exactly the same model even after many incremental gradient updates.
    fitted = (0..n)
        .map(|i| (0..n).map(|j| alpha[j] * labels[j] * k[i][j]).sum::<f64>())
        .collect();
    let (_, _, up, low) = working(&labels, &alpha, &fitted, cost);
    let gap = (up - low).max(0.);
    let converged = gap <= tol;
    let eps = 32. * f64::EPSILON * cost;
    let free = (0..n)
        .filter(|i| alpha[*i] > eps && alpha[*i] < cost - eps)
        .collect::<Vec<_>>();
    let bias = if free.is_empty() {
        up / 2. + low / 2.
    } else {
        free.iter()
            .map(|i| (labels[*i] - fitted[*i]) / free.len() as f64)
            .sum::<f64>()
    };
    let sv = (0..n).filter(|i| alpha[*i] > 0.).collect::<Vec<_>>();
    if 5 + 3 * p + sv.len() * (p + 1) > 16384 {
        return Err("SVM portable model exceeds argument limit".into());
    }
    let decision = x
        .iter()
        .enumerate()
        .map(|(i, row)| {
            if i < n {
                fitted[i] + bias
            } else {
                bias + sv
                    .iter()
                    .map(|j| alpha[*j] * labels[*j] * kernel(row, &x[*j], kind, gamma))
                    .sum::<f64>()
            }
        })
        .collect::<Vec<_>>();
    let pred = decision
        .iter()
        .map(|f| if *f >= 0. { 1. } else { 0. })
        .collect::<Vec<_>>();
    let norm_squared = alpha
        .iter()
        .zip(&labels)
        .zip(&fitted)
        .map(|((a, y), f)| a * y * f)
        .sum::<f64>()
        .max(0.);
    let dual = alpha.iter().sum::<f64>() - 0.5 * norm_squared;
    let primal = 0.5 * norm_squared
        + cost
            * (0..n)
                .map(|i| (1. - labels[i] * decision[i]).max(0.))
                .sum::<f64>();
    let equality = dot(&alpha, &labels).abs();
    let box_violation = alpha
        .iter()
        .map(|v| (-v).max(v - cost).max(0.))
        .fold(0_f64, f64::max);
    let mut saved = vec![p as f64, kind as f64, gamma, bias, sv.len() as f64];
    saved.extend(&standard.origins);
    saved.extend(&standard.offsets);
    saved.extend(&standard.scales);
    saved.extend(sv.iter().map(|i| alpha[*i] * labels[*i]));
    for i in &sv {
        saved.extend(&x[*i]);
    }
    Ok(format!("{{\"version\":1,\"algorithm\":\"binary_C_SVC_pair_coordinate_SMO\",\"kernel\":\"{}\",\"cost\":{},\"gamma\":{},\"n_train\":{},\"bias\":{},\"support_indices\":{},\"dual_alphas\":{},\"decision_values\":{},\"predictions\":{},\"dual_objective\":{},\"primal_objective\":{},\"duality_gap\":{},\"kkt_gap\":{},\"equality_residual\":{},\"box_violation\":{},\"converged\":{},\"iterations\":{},\"predict_parameters\":{}}}",if kind==0{"linear"}else{"rbf"},cost,gamma,n,bias,nums(&sv.iter().map(|i|*i as f64).collect::<Vec<_>>()),nums(&alpha),nums(&decision),nums(&pred),dual,primal,(primal-dual).max(0.),gap,equality,box_violation,converged,iterations,nums(&saved)))
}
fn svm_predict(d: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    if a.len() < 5 {
        return Err("saved SVM model missing".into());
    }
    let p = int(a[0], 1, 16)?;
    let kind = int(a[1], 0, 1)?;
    let n = int(a[4], 1, 512)?;
    if d[0].len() != p
        || a.len() != 5 + 3 * p + n * (p + 1)
        || a[2] <= 0.
        || a[5 + 2 * p..5 + 3 * p].iter().any(|v| *v <= 0.)
        || d.len().saturating_mul(n).saturating_mul(p) > MAX_WORK
    {
        return Err("invalid saved SVM model dimensions, scales or work".into());
    }
    let standard = Standard {
        origins: a[5..5 + p].to_vec(),
        offsets: a[5 + p..5 + 2 * p].to_vec(),
        scales: a[5 + 2 * p..5 + 3 * p].to_vec(),
    };
    let coef = &a[5 + 3 * p..5 + 3 * p + n];
    let support = &a[5 + 3 * p + n..];
    let f = d
        .iter()
        .map(|r| {
            standard.row(r).map(|x| {
                a[3] + (0..n)
                    .map(|i| coef[i] * kernel(&x, &support[i * p..(i + 1) * p], kind, a[2]))
                    .sum::<f64>()
            })
        })
        .collect::<Result<Vec<_>, _>>()?;
    let predictions = f
        .iter()
        .map(|v| if *v >= 0. { 1. } else { 0. })
        .collect::<Vec<_>>();
    Ok(format!(
        "{{\"version\":1,\"algorithm\":\"svm_predict\",\"decision_values\":{},\"predictions\":{}}}",
        nums(&f),
        nums(&predictions)
    ))
}
