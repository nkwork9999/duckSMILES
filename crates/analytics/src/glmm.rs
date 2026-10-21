//! Bounded Gaussian random-intercept Bernoulli/logit and Poisson/log GLMMs.
//! The objective is integrated likelihood evaluated by adaptive Gauss-Hermite
//! quadrature. No random slopes, correlated random effects, repeated covariance,
//! grouped binomial trials, offsets or non-Gaussian random effects are supplied.

const MAX_WORK: usize = 200_000_000;

pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 10000
        || width > 128
        || data.len() > 1_000_000 / width
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite()))
    {
        return Err("GLMM requires finite rectangular bounded data/parameters".into());
    }
    match operation {
        "glmm_binomial" => fit(data, args, 0),
        "glmm_poisson" => fit(data, args, 1),
        "glmm_predict" => predict(data, args),
        _ => Err("unknown GLMM operation".into()),
    }
}
struct Frame {
    x: Vec<Vec<f64>>,
    y: Vec<f64>,
    constants: Vec<f64>,
    groups: Vec<Vec<usize>>,
    ids: Vec<f64>,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    family: usize,
}
impl Frame {
    fn new(data: &[Vec<f64>], family: usize) -> Result<Self, String> {
        let p = data[0]
            .len()
            .checked_sub(2)
            .ok_or("GLMM rows need [features...,response,group]")?;
        let n = data.len();
        if p > 16 || n > 2000 || n <= p + 2 {
            return Err("GLMM needs n>features+2, n<=2000 and <=16 fixed features".into());
        }
        let origins = data[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| mean(data.iter().map(|r| r[j] - origins[j])))
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = finite(
                    norm(data.iter().map(|r| (r[j] - origins[j]) - offsets[j])) / (n as f64).sqrt(),
                    "feature scale",
                )?;
                if s == 0. {
                    Err("fixed design contains a constant feature".into())
                } else {
                    Ok(s)
                }
            })
            .collect::<Result<Vec<_>, String>>()?;
        let x = data
            .iter()
            .map(|r| {
                let mut row = vec![1.];
                for j in 0..p {
                    row.push(finite(
                        ((r[j] - origins[j]) - offsets[j]) / scales[j],
                        "standardized feature",
                    )?);
                }
                Ok(row)
            })
            .collect::<Result<Vec<_>, String>>()?;
        // Twice-reorthogonalized QR detects redundant fixed terms before fitting.
        let mut orthogonal: Vec<Vec<f64>> = Vec::new();
        for j in 0..=p {
            let mut column = x.iter().map(|r| r[j]).collect::<Vec<_>>();
            for _ in 0..2 {
                for q in &orthogonal {
                    let d = dot(&column, q);
                    for (v, u) in column.iter_mut().zip(q) {
                        *v -= d * u;
                    }
                }
            }
            let length = norm(column.iter().copied());
            if length <= 1e-10 * (n as f64).sqrt() {
                return Err("fixed design is rank deficient".into());
            }
            orthogonal.push(column.iter().map(|v| v / length).collect());
        }
        let mut ids = Vec::new();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        let mut y = Vec::new();
        let mut constants = Vec::new();
        for (i, row) in data.iter().enumerate() {
            let response = row[p];
            if (family == 0 && response != 0. && response != 1.)
                || (family == 1 && (response < 0. || response.fract() != 0. || response > 1e9))
            {
                return Err("GLMM response must be binary, or a Poisson integer in [0,1e9]".into());
            }
            let id = row[p + 1];
            valid_id(id)?;
            let g = if let Some(g) = ids.iter().position(|v| *v == id) {
                g
            } else {
                ids.push(id);
                groups.push(Vec::new());
                ids.len() - 1
            };
            groups[g].push(i);
            y.push(response);
            constants.push(if family == 1 {
                log_gamma(response + 1.)
            } else {
                0.
            });
        }
        if !(2..=200).contains(&groups.len()) || groups.iter().filter(|g| g.len() > 1).count() < 2 {
            return Err("GLMM needs 2..200 groups including at least two repeated groups".into());
        }
        let average = mean(y.iter().copied())?;
        if average == 0. || (family == 0 && average == 1.) {
            return Err("response has no finite intercept maximum likelihood estimate".into());
        }
        Ok(Self {
            x,
            y,
            constants,
            groups,
            ids,
            origins,
            offsets,
            scales,
            family,
        })
    }
    fn eta(&self, beta: &[f64]) -> Result<Vec<f64>, String> {
        self.x
            .iter()
            .map(|x| finite(dot(x, beta), "linear predictor"))
            .collect()
    }
}
// A fixed-effect recession direction also increases the group-integrated
// likelihood pointwise. Bernoulli requires signed Xd>=0 with one strict row;
// Poisson requires Xd=0 on positive counts and Xd<=0 on zeros, with one strict
// zero. The bounded LP maximizes their average improvement over |d_j|<=1.
// Splitting each direction into positive/negative parts provides a feasible
// all-zero slack basis, so no infeasible Phase I or fitted-tolerance heuristic
// is required. Bland indices resolve degenerate pivots deterministically.
#[allow(clippy::needless_range_loop)]
fn check_finite_mle(frame: &Frame) -> Result<(), String> {
    let d = frame.x[0].len();
    let variables = 2 * d;
    let mut constraints: Vec<Vec<f64>> = Vec::new();
    let mut objective = vec![0.; d];
    for (x, &y) in frame.x.iter().zip(&frame.y) {
        let scale = norm(x.iter().copied());
        let direction = if frame.family == 0 {
            x.iter()
                .map(|v| (2. * y - 1.) * v / scale)
                .collect::<Vec<_>>()
        } else {
            x.iter().map(|v| -v / scale).collect::<Vec<_>>()
        };
        constraints.push(direction.iter().map(|v| -v).collect());
        if frame.family == 1 && y > 0. {
            constraints.push(direction.clone());
        } else {
            for j in 0..d {
                objective[j] += direction[j];
            }
        }
    }
    let m = constraints.len() + variables;
    let mut table = vec![vec![0.; variables + 1]; m + 1];
    for (i, row) in constraints.iter().enumerate() {
        for j in 0..d {
            table[i][j] = row[j];
            table[i][j + d] = -row[j];
        }
    }
    for j in 0..variables {
        table[constraints.len() + j][j] = 1.;
        table[constraints.len() + j][variables] = 1.;
    }
    for j in 0..d {
        table[m][j] = -objective[j];
        table[m][j + d] = objective[j];
    }
    let mut basic = (variables..variables + m).collect::<Vec<_>>();
    let mut nonbasic = (0..variables).collect::<Vec<_>>();
    let mut resolved = false;
    for iteration in 0_usize..4096 {
        if (iteration + 1)
            .saturating_mul(m + 1)
            .saturating_mul(variables + 1)
            > MAX_WORK
        {
            return Err("finite-MLE separation LP exceeds numerical work limit".into());
        }
        let entering = (0..variables)
            .filter(|&j| table[m][j] < -1e-12)
            .min_by_key(|&j| nonbasic[j]);
        let Some(s) = entering else {
            resolved = true;
            break;
        };
        let leaving = (0..m).filter(|&i| table[i][s] > 1e-12).min_by(|&i, &k| {
            let a = table[i][variables] / table[i][s];
            let b = table[k][variables] / table[k][s];
            if (a - b).abs() <= 1e-12 * (1. + a.abs() + b.abs()) {
                basic[i].cmp(&basic[k])
            } else {
                a.total_cmp(&b)
            }
        });
        let r = leaving.ok_or("finite-MLE separation LP became numerically unbounded")?;
        let inverse = 1. / table[r][s];
        for i in 0..=m {
            if i != r {
                for j in 0..=variables {
                    if j != s {
                        table[i][j] -= table[r][j] * table[i][s] * inverse;
                    }
                }
            }
        }
        for j in 0..=variables {
            if j != s {
                table[r][j] *= inverse;
            }
        }
        for i in 0..=m {
            if i != r {
                table[i][s] *= -inverse;
            }
        }
        table[r][s] = inverse;
        std::mem::swap(&mut basic[r], &mut nonbasic[s]);
        if table.iter().flatten().any(|v| !v.is_finite()) {
            return Err("finite-MLE separation LP overflows".into());
        }
    }
    if !resolved {
        return Err("finite-MLE separation LP iteration limit; existence unconfirmed".into());
    }
    let mut solution = vec![0.; variables];
    for i in 0..m {
        if basic[i] < variables {
            solution[basic[i]] = table[i][variables];
        }
    }
    let direction = (0..d)
        .map(|j| solution[j] - solution[j + d])
        .collect::<Vec<_>>();
    let minimum = constraints
        .iter()
        .map(|r| -dot(r, &direction))
        .fold(f64::INFINITY, f64::min);
    let improvement = dot(&objective, &direction);
    if minimum < -1e-8 || table[m][variables] < -1e-8 {
        return Err("finite-MLE separation LP lost feasibility".into());
    }
    if improvement > 1e-10 || table[m][variables] > 1e-10 {
        return Err("fixed design has a complete/quasi or zero-count recession direction; no identified finite GLMM MLE".into());
    }
    Ok(())
}

fn valid_id(id: f64) -> Result<(), String> {
    if id.fract() != 0. || id.abs() > 9_007_199_254_740_991. {
        Err("group identifiers must be exact integers of absolute value <=2^53-1".into())
    } else {
        Ok(())
    }
}
struct Rule {
    nodes: Vec<f64>,
    weights: Vec<f64>,
}
fn normal_rule(order: usize) -> Result<Rule, String> {
    let mut jacobi = vec![vec![0.; order]; order];
    for j in 1..order {
        jacobi[j - 1][j] = (j as f64).sqrt();
        jacobi[j][j - 1] = jacobi[j - 1][j];
    }
    let e = eigen(&jacobi)?;
    let mut nodes = e.values;
    let mut weights = e.vectors[0].iter().map(|v| v * v).collect::<Vec<_>>();
    nodes.reverse();
    weights.reverse();
    let sum: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= sum;
        if *w <= 0. {
            return Err("Gaussian quadrature weight underflow".into());
        }
    }
    Ok(Rule { nodes, weights })
}
struct Posterior {
    log_integral: f64,
    mode: f64,
    mean: f64,
    variance: f64,
    nodes: Vec<f64>,
    weights: Vec<f64>,
}
fn conditional(frame: &Frame, group: &[usize], eta: &[f64], b: f64) -> f64 {
    group
        .iter()
        .map(|&i| {
            let t = eta[i] + b;
            if frame.family == 0 {
                if frame.y[i] == 1. {
                    -softplus(-t)
                } else {
                    -softplus(t)
                }
            } else if t > 709. {
                -f64::INFINITY
            } else {
                frame.y[i] * t - t.exp() - frame.constants[i]
            }
        })
        .sum()
}
fn mode_score(frame: &Frame, group: &[usize], eta: &[f64], b: f64, v: f64) -> (f64, f64) {
    let mut score = -b / v;
    let mut curvature = 1. / v;
    for &i in group {
        let t = eta[i] + b;
        let mu = if frame.family == 0 {
            sigmoid(t)
        } else {
            t.exp()
        };
        score += frame.y[i] - mu;
        curvature += if frame.family == 0 {
            mu * (1. - mu)
        } else {
            mu
        };
    }
    (score, curvature)
}
fn posterior(
    frame: &Frame,
    group: &[usize],
    eta: &[f64],
    v: f64,
    rule: &Rule,
) -> Result<Posterior, String> {
    if v == 0. {
        return Ok(Posterior {
            log_integral: finite(
                conditional(frame, group, eta, 0.),
                "boundary log likelihood",
            )?,
            mode: 0.,
            mean: 0.,
            variance: 0.,
            nodes: vec![0.],
            weights: vec![1.],
        });
    }
    let mut lo = -v.sqrt().max(1.);
    let mut hi = -lo;
    for _ in 0..60 {
        if mode_score(frame, group, eta, lo, v).0 > 0. {
            break;
        }
        lo *= 2.;
    }
    for _ in 0..60 {
        if mode_score(frame, group, eta, hi, v).0 < 0. {
            break;
        }
        hi *= 2.;
    }
    if mode_score(frame, group, eta, lo, v).0 <= 0. || mode_score(frame, group, eta, hi, v).0 >= 0.
    {
        return Err("posterior mode could not be bracketed".into());
    }
    let mut b = 0_f64.clamp(lo, hi);
    let mut root_ok = false;
    for _ in 0..100 {
        let (score, curvature) = mode_score(frame, group, eta, b, v);
        if score.abs()
            <= 1e-11 * (1. + group.iter().map(|&i| frame.y[i]).sum::<f64>() + b.abs() / v)
        {
            root_ok = true;
            break;
        }
        if score > 0. {
            lo = b;
        } else {
            hi = b;
        }
        let candidate = b + score / curvature;
        let next = if candidate.is_finite() && candidate > lo && candidate < hi {
            candidate
        } else {
            lo + (hi - lo) * 0.5
        };
        if (next - b).abs() <= 1e-12 * (1. + b.abs()) {
            b = next;
            root_ok = true;
            break;
        }
        b = next;
    }
    if !root_ok {
        return Err("posterior mode root did not converge".into());
    }
    let curvature = mode_score(frame, group, eta, b, v).1;
    if !curvature.is_finite() || curvature <= 0. {
        return Err("posterior curvature is not finite positive".into());
    }
    let nodes = rule
        .nodes
        .iter()
        .map(|z| b + z / curvature.sqrt())
        .collect::<Vec<_>>();
    let logs = nodes
        .iter()
        .zip(&rule.nodes)
        .zip(&rule.weights)
        .map(|((&point, &z), &w)| {
            w.ln() + conditional(frame, group, eta, point)
                - point * point / (2. * v)
                - 0.5 * (2. * std::f64::consts::PI * v).ln()
                + 0.5 * z * z
        })
        .collect::<Vec<_>>();
    let maximum = logs.iter().copied().fold(-f64::INFINITY, f64::max);
    if !maximum.is_finite() {
        return Err("adaptive quadrature likelihood is not finite".into());
    }
    let mut weights = logs.iter().map(|x| (x - maximum).exp()).collect::<Vec<_>>();
    let total: f64 = weights.iter().sum();
    for w in &mut weights {
        *w /= total;
    }
    let log_integral = finite(
        maximum + total.ln() + 0.5 * (2. * std::f64::consts::PI / curvature).ln(),
        "integrated log likelihood",
    )?;
    let posterior_mean = dot(&nodes, &weights);
    let variance = finite(
        nodes
            .iter()
            .zip(&weights)
            .map(|(x, w)| w * (x - posterior_mean).powi(2))
            .sum(),
        "posterior variance",
    )?;
    Ok(Posterior {
        log_integral,
        mode: b,
        mean: posterior_mean,
        variance,
        nodes,
        weights,
    })
}
fn likelihood(
    frame: &Frame,
    theta: &[f64],
    rule: &Rule,
    free_variance: bool,
) -> Result<f64, String> {
    let d = frame.x[0].len();
    let eta = frame.eta(&theta[..d])?;
    let variance = if free_variance { theta[d].exp_m1() } else { 0. };
    finite(
        -frame
            .groups
            .iter()
            .map(|g| posterior(frame, g, &eta, variance, rule).map(|p| p.log_integral))
            .collect::<Result<Vec<_>, _>>()?
            .iter()
            .sum::<f64>()
            / frame.y.len() as f64,
        "average negative log likelihood",
    )
}
struct Optimum {
    theta: Vec<f64>,
    value: f64,
    converged: bool,
    iterations: usize,
    gradient_max: f64,
    line_search_failed: bool,
}
fn active(x: f64, g: f64, lo: f64, hi: f64) -> bool {
    let slack = 1e-12 * (hi - lo).abs().min(1.);
    (x <= lo + slack && g > 0.) || (x >= hi - slack && g < 0.)
}
fn objective_gradient(
    frame: &Frame,
    theta: &[f64],
    rule: &Rule,
    free_variance: bool,
    lower: &[f64],
    upper: &[f64],
    value: f64,
) -> Result<Vec<f64>, String> {
    let d = frame.x[0].len();
    if !free_variance || theta[d] == 0. {
        let eta = frame.eta(&theta[..d])?;
        let mu = eta
            .iter()
            .map(|x| {
                if frame.family == 0 {
                    sigmoid(*x)
                } else {
                    x.exp()
                }
            })
            .collect::<Vec<_>>();
        let mut gradient = (0..d)
            .map(|j| {
                frame
                    .x
                    .iter()
                    .enumerate()
                    .map(|(i, r)| r[j] * (mu[i] - frame.y[i]))
                    .sum::<f64>()
                    / frame.y.len() as f64
            })
            .collect::<Vec<_>>();
        if free_variance {
            let score = frame
                .groups
                .iter()
                .map(|group| {
                    let s: f64 = group.iter().map(|&i| frame.y[i] - mu[i]).sum();
                    let information: f64 = group
                        .iter()
                        .map(|&i| {
                            if frame.family == 0 {
                                mu[i] * (1. - mu[i])
                            } else {
                                mu[i]
                            }
                        })
                        .sum();
                    -0.5 * (s * s - information)
                })
                .sum::<f64>()
                / frame.y.len() as f64;
            gradient.push(score);
        }
        if gradient.iter().any(|v| !v.is_finite()) {
            return Err("GLMM gradient overflows".into());
        }
        return Ok(gradient);
    }
    (0..theta.len())
        .map(|j| {
            let step = 2e-5 * (1. + theta[j].abs());
            let mut plus = theta.to_vec();
            let mut minus = theta.to_vec();
            plus[j] = (theta[j] + step).min(upper[j]);
            minus[j] = (theta[j] - step).max(lower[j]);
            let fp = if plus[j] == theta[j] {
                value
            } else {
                likelihood(frame, &plus, rule, free_variance)?
            };
            let fm = if minus[j] == theta[j] {
                value
            } else {
                likelihood(frame, &minus, rule, free_variance)?
            };
            finite((fp - fm) / (plus[j] - minus[j]), "likelihood gradient")
        })
        .collect()
}
fn optimize(
    frame: &Frame,
    initial: Vec<f64>,
    rule: &Rule,
    free_variance: bool,
    upper_variance: f64,
    maxiter: usize,
    tol: f64,
) -> Result<Optimum, String> {
    let d = initial.len();
    let mut lower = vec![-30.; d];
    let mut upper = vec![30.; d];
    if free_variance {
        lower[d - 1] = 0.;
        upper[d - 1] = upper_variance.ln_1p();
    }
    let mut theta = initial;
    let mut current = likelihood(frame, &theta, rule, free_variance)?;
    let mut gradient =
        objective_gradient(frame, &theta, rule, free_variance, &lower, &upper, current)?;
    let mut inverse = identity(d);
    let mut iterations = 0;
    let mut line_search_failed = false;
    for iteration in 0..maxiter {
        let active_set = (0..d)
            .map(|j| active(theta[j], gradient[j], lower[j], upper[j]))
            .collect::<Vec<_>>();
        let projected = (0..d)
            .map(|j| if active_set[j] { 0. } else { gradient[j] })
            .collect::<Vec<_>>();
        if projected.iter().map(|g| g.abs()).fold(0_f64, f64::max) <= tol {
            break;
        }
        let mut direction = inverse
            .iter()
            .map(|r| -dot(r, &projected))
            .collect::<Vec<_>>();
        for j in 0..d {
            if active_set[j] {
                direction[j] = 0.;
            }
        }
        if dot(&direction, &projected) >= 0. {
            direction = projected.iter().map(|g| -g).collect();
            inverse = identity(d);
        }
        let mut alpha = 1.;
        let mut accepted = None;
        for _ in 0..50 {
            let candidate = (0..d)
                .map(|j| (theta[j] + alpha * direction[j]).clamp(lower[j], upper[j]))
                .collect::<Vec<_>>();
            let step = (0..d).map(|j| candidate[j] - theta[j]).collect::<Vec<_>>();
            let slope = dot(&gradient, &step);
            if slope < 0. {
                if let Ok(value) = likelihood(frame, &candidate, rule, free_variance) {
                    if value <= current + 1e-4 * slope + 16. * f64::EPSILON * (1. + current.abs()) {
                        accepted = Some((candidate, step, value));
                        break;
                    }
                }
            }
            alpha *= 0.5;
        }
        let Some((candidate, step, value)) = accepted else {
            line_search_failed = true;
            break;
        };
        let next_gradient = objective_gradient(
            frame,
            &candidate,
            rule,
            free_variance,
            &lower,
            &upper,
            value,
        )?;
        let changed = (0..d)
            .any(|j| active_set[j] != active(candidate[j], next_gradient[j], lower[j], upper[j]));
        let change = (0..d)
            .map(|j| {
                if active_set[j] {
                    0.
                } else {
                    next_gradient[j] - gradient[j]
                }
            })
            .collect::<Vec<_>>();
        let curvature = dot(&step, &change);
        if !changed
            && curvature > 1e-12 * norm(step.iter().copied()) * norm(change.iter().copied())
            && curvature > 0.
        {
            let hy = inverse.iter().map(|r| dot(r, &change)).collect::<Vec<_>>();
            let coefficient = (1. + dot(&change, &hy) / curvature) / curvature;
            for i in 0..d {
                for j in 0..d {
                    inverse[i][j] = finite(
                        inverse[i][j] + coefficient * step[i] * step[j]
                            - (step[i] * hy[j] + hy[i] * step[j]) / curvature,
                        "inverse Hessian",
                    )?;
                }
            }
        } else {
            inverse = identity(d);
        }
        theta = candidate;
        gradient = next_gradient;
        current = value;
        iterations = iteration + 1;
    }
    let gradient_max = (0..d)
        .map(|j| {
            if active(theta[j], gradient[j], lower[j], upper[j]) {
                0.
            } else {
                gradient[j].abs()
            }
        })
        .fold(0_f64, f64::max);
    Ok(Optimum {
        theta,
        value: current,
        converged: gradient_max <= tol,
        iterations,
        gradient_max,
        line_search_failed,
    })
}
#[allow(clippy::needless_range_loop)]
fn fixed_covariance(
    frame: &Frame,
    theta: &[f64],
    rule: &Rule,
    boundary: bool,
) -> Result<Vec<Vec<f64>>, String> {
    let d = frame.x[0].len();
    let information = if boundary {
        let eta = frame.eta(&theta[..d])?;
        (0..d)
            .map(|j| {
                (0..d)
                    .map(|k| {
                        finite(
                            frame
                                .x
                                .iter()
                                .zip(&eta)
                                .map(|(x, t)| {
                                    let mu = if frame.family == 0 {
                                        sigmoid(*t)
                                    } else {
                                        t.exp()
                                    };
                                    x[j] * x[k]
                                        * if frame.family == 0 {
                                            mu * (1. - mu)
                                        } else {
                                            mu
                                        }
                                })
                                .sum(),
                            "fixed effect information",
                        )
                    })
                    .collect()
            })
            .collect::<Result<Vec<Vec<f64>>, String>>()?
    } else {
        let m = theta.len();
        let base = likelihood(frame, theta, rule, true)?;
        let mut hessian = vec![vec![0.; m]; m];
        for j in 0..m {
            let h = 2e-4 * (1. + theta[j].abs());
            if j == d && theta[j] <= h {
                return Err("variance too near boundary for joint information".into());
            }
            for k in 0..=j {
                let t = 2e-4 * (1. + theta[k].abs());
                let value = if j == k {
                    let mut a = theta.to_vec();
                    let mut b = theta.to_vec();
                    a[j] += h;
                    b[j] -= h;
                    (likelihood(frame, &a, rule, true)? - 2. * base
                        + likelihood(frame, &b, rule, true)?)
                        / (h * h)
                } else {
                    let mut pp = theta.to_vec();
                    let mut pm = theta.to_vec();
                    let mut mp = theta.to_vec();
                    let mut mm = theta.to_vec();
                    pp[j] += h;
                    pp[k] += t;
                    pm[j] += h;
                    pm[k] -= t;
                    mp[j] -= h;
                    mp[k] += t;
                    mm[j] -= h;
                    mm[k] -= t;
                    (likelihood(frame, &pp, rule, true)?
                        - likelihood(frame, &pm, rule, true)?
                        - likelihood(frame, &mp, rule, true)?
                        + likelihood(frame, &mm, rule, true)?)
                        / (4. * h * t)
                };
                hessian[j][k] =
                    finite(value * frame.y.len() as f64, "observed likelihood Hessian")?;
                hessian[k][j] = hessian[j][k];
            }
        }
        hessian
    };
    let inverse = symmetric_power(&information, -1.)?;
    let mut transform = identity(d);
    for j in 1..d {
        transform[j][j] = finite(1. / frame.scales[j - 1], "coefficient unit transform")?;
        transform[0][j] = finite(
            -frame.origins[j - 1] / frame.scales[j - 1]
                - frame.offsets[j - 1] / frame.scales[j - 1],
            "intercept unit transform",
        )?;
    }
    let fixed = inverse[..d]
        .iter()
        .map(|r| r[..d].to_vec())
        .collect::<Vec<_>>();
    multiply(&multiply(&transform, &fixed)?, &transpose(&transform))
}
fn fit(data: &[Vec<f64>], args: &[f64], family: usize) -> Result<String, String> {
    if args.len() > 4 {
        return Err(
            "GLMM args are [maxiter=200,tolerance=1e-6,odd_order=15,variance_upper=100]".into(),
        );
    }
    let maxiter = integer(
        args.first().copied().unwrap_or(200.),
        1,
        1000,
        "iteration limit",
    )?;
    let tol = args.get(1).copied().unwrap_or(1e-6);
    if !(1e-10..=1e-3).contains(&tol) {
        return Err("GLMM tolerance must be in [1e-10,1e-3]".into());
    }
    let order = integer(
        args.get(2).copied().unwrap_or(15.),
        3,
        31,
        "quadrature order",
    )?;
    if order % 2 != 1 {
        return Err("quadrature order must be odd".into());
    }
    let variance_upper = args.get(3).copied().unwrap_or(100.);
    if variance_upper <= 0. || variance_upper > 1000. {
        return Err("variance upper bound must be in (0,1000]".into());
    }
    let frame = Frame::new(data, family)?;
    let d = frame.x[0].len();
    work(data.len(), (d + 1) * order, 4 * maxiter)?;
    check_finite_mle(&frame)?;
    let rule = normal_rule(order)?;
    let reference_rule = normal_rule(order + 4)?;
    let average = mean(frame.y.iter().copied())?;
    let initial_intercept = if family == 0 {
        (average / (1. - average)).ln()
    } else {
        average.ln()
    };
    let mut initial = vec![0.; d];
    initial[0] = initial_intercept.clamp(-30., 30.);
    let fixed = optimize(&frame, initial, &rule, false, variance_upper, maxiter, tol)?;
    if family == 0 {
        let margins = frame
            .eta(&fixed.theta)?
            .iter()
            .zip(&frame.y)
            .map(|(t, y)| (2. * y - 1.) * t)
            .collect::<Vec<_>>();
        // Any fixed-effect direction with nonnegative signed margins and at
        // least one positive margin proves complete or quasi separation.
        if margins.iter().all(|m| *m >= 0.) && margins.iter().any(|m| *m > 1e-9) {
            return Err(
                "fixed effects exhibit a separating direction; finite Bernoulli MLE does not exist"
                    .into(),
            );
        }
    }
    let mut best = Optimum {
        theta: {
            let mut t = fixed.theta.clone();
            t.push(0.);
            t
        },
        value: fixed.value,
        converged: fixed.converged,
        iterations: fixed.iterations,
        gradient_max: fixed.gradient_max,
        line_search_failed: fixed.line_search_failed,
    };
    // A separately optimized boundary is always compared with interior starts.
    // Its one-sided variance score must also satisfy the constrained criterion.
    let mut lo = vec![-30.; d + 1];
    let mut hi = vec![30.; d + 1];
    lo[d] = 0.;
    hi[d] = variance_upper.ln_1p();
    let g = objective_gradient(&frame, &best.theta, &rule, true, &lo, &hi, best.value)?;
    best.gradient_max = best.gradient_max.max((-g[d]).max(0.));
    best.converged &= best.gradient_max <= tol;
    let mut starts_converged = usize::from(best.converged);
    for variance in [0.25_f64, 1., 4.] {
        let mut start = fixed.theta.clone();
        start.push(variance.min(variance_upper).ln_1p());
        if let Ok(candidate) = optimize(&frame, start, &rule, true, variance_upper, maxiter, tol) {
            starts_converged += usize::from(candidate.converged);
            if candidate.value < best.value - 1e-11 * (1. + best.value.abs())
                || ((candidate.value - best.value).abs() <= 1e-11 * (1. + best.value.abs())
                    && candidate.converged
                    && !best.converged)
            {
                best = candidate;
            }
        }
    }
    let variance = best.theta[d].exp_m1();
    let boundary = variance == 0.;
    let eta = frame.eta(&best.theta[..d])?;
    let low = frame
        .groups
        .iter()
        .map(|g| posterior(&frame, g, &eta, variance, &rule))
        .collect::<Result<Vec<_>, _>>()?;
    let post = frame
        .groups
        .iter()
        .map(|g| posterior(&frame, g, &eta, variance, &reference_rule))
        .collect::<Result<Vec<_>, _>>()?;
    let group_delta = low
        .iter()
        .zip(&post)
        .map(|(a, b)| (a.log_integral - b.log_integral).abs())
        .collect::<Vec<_>>();
    let log_likelihood = finite(post.iter().map(|p| p.log_integral).sum(), "log likelihood")?;
    let log_likelihood_order_delta = finite(
        (low.iter().map(|p| p.log_integral).sum::<f64>() - log_likelihood).abs(),
        "quadrature order difference",
    )?;
    let quadrature_tolerance = tol.max(1e-8) * (1. + log_likelihood.abs());
    let quadrature_converged = log_likelihood_order_delta <= quadrature_tolerance
        && group_delta
            .iter()
            .all(|v| *v <= quadrature_tolerance / frame.groups.len() as f64);
    let coefficient_bound = best.theta[..d].iter().any(|v| v.abs() >= 30. - 1e-6);
    let variance_upper_bound = best.theta[d] >= variance_upper.ln_1p() * (1. - 1e-8);
    let converged =
        best.converged && quadrature_converged && !coefficient_bound && !variance_upper_bound;
    let slopes = (0..d - 1)
        .map(|j| {
            finite(
                best.theta[j + 1] / frame.scales[j],
                "original unit coefficient",
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let intercept = finite(
        best.theta[0]
            - slopes
                .iter()
                .enumerate()
                .map(|(j, b)| b * frame.origins[j] + b * frame.offsets[j])
                .sum::<f64>(),
        "original unit intercept",
    )?;
    let means = post.iter().map(|p| p.mean).collect::<Vec<_>>();
    let modes = post.iter().map(|p| p.mode).collect::<Vec<_>>();
    let variances = post.iter().map(|p| p.variance).collect::<Vec<_>>();
    let saved_order = if boundary { 1 } else { order + 4 };
    let mut portable = vec![
        1.,
        family as f64,
        (d - 1) as f64,
        frame.groups.len() as f64,
        saved_order as f64,
        variance,
        best.theta[0],
    ];
    portable.extend(&frame.origins);
    portable.extend(&frame.offsets);
    portable.extend(&slopes);
    portable.extend(&frame.ids);
    portable.extend(&means);
    portable.extend(&modes);
    portable.extend(&variances);
    for p in &post {
        portable.extend(&p.nodes);
    }
    for p in &post {
        portable.extend(&p.weights);
    }
    if portable.len() > 16384 || portable.iter().any(|v| !v.is_finite()) {
        return Err("saved GLMM model exceeds finite portable bounds".into());
    }
    let covariance = fixed_covariance(&frame, &best.theta, &reference_rule, boundary).ok();
    let se = covariance.as_ref().map(|c| {
        c.iter()
            .enumerate()
            .map(|(j, r)| if r[j] >= 0. { Some(r[j].sqrt()) } else { None })
            .collect::<Vec<_>>()
    });
    let covariance_json = covariance.as_ref().map_or("null".into(), |c| matrix(c));
    let coefficient_se = se.as_ref().map_or("null".into(), |v| {
        format!(
            "[{}]",
            v[1..]
                .iter()
                .map(|x| optional(*x))
                .collect::<Vec<_>>()
                .join(",")
        )
    });
    let intercept_se = se.as_ref().and_then(|s| s[0]);
    let parameter_count = d + 1;
    Ok(format!("{{\"version\":1,\"algorithm\":\"adaptive_gauss_hermite_gaussian_random_intercept_ml\",\"family\":\"{}\",\"n_samples\":{},\"n_groups\":{},\"intercept\":{intercept},\"coefficients\":{},\"centered_intercept\":{},\"random_intercept_variance\":{variance},\"log_likelihood\":{log_likelihood},\"aic\":{},\"parameter_count\":{parameter_count},\"boundary\":{boundary},\"variance_upper_bound\":{variance_upper_bound},\"coefficient_guard_bound\":{coefficient_bound},\"optimizer_converged\":{},\"converged\":{converged},\"iterations\":{},\"projected_gradient_max\":{},\"line_search_failed\":{},\"converged_starts\":{starts_converged},\"quadrature_order\":{order},\"validation_order\":{},\"quadrature_tolerance\":{quadrature_tolerance},\"quadrature_converged\":{quadrature_converged},\"log_likelihood_order_delta\":{log_likelihood_order_delta},\"group_log_likelihood_order_deltas\":{},\"group_ids\":{},\"posterior_means\":{},\"posterior_modes\":{},\"posterior_variances\":{},\"fixed_effect_covariance\":{covariance_json},\"intercept_standard_error\":{},\"coefficient_standard_errors\":{coefficient_se},\"standard_error_method\":\"{}\",\"predict_parameters\":{},\"prediction_contract\":\"conditional=posterior_mean_random_effect;posterior=saved_posterior_integral;population=Gaussian_prior_integral\",\"quadrature_error_convention\":\"q_versus_q_plus_4_consistency_not_rigorous_bound\"}}",if family==0{"bernoulli_logit"}else{"poisson_log"},data.len(),frame.groups.len(),numbers(&slopes),best.theta[0],-2.*log_likelihood+2.*parameter_count as f64,best.converged,best.iterations,best.gradient_max,best.line_search_failed,order+4,numbers(&group_delta),numbers(&frame.ids),numbers(&means),numbers(&modes),numbers(&variances),optional(intercept_se),if covariance.is_none(){"unidentified_information_null"}else if boundary{"fixed_glm_information_at_variance_zero"}else{"inverse_joint_observed_integrated_likelihood_hessian"},numbers(&portable)))
}
fn predict(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 7 || args[0] != 1. {
        return Err("GLMM prediction needs unchanged version1 predict_parameters".into());
    }
    let family = integer(args[1], 0, 1, "family")?;
    let p = integer(args[2], 0, 16, "feature count")?;
    let groups = integer(args[3], 1, 200, "saved group count")?;
    let q = integer(args[4], 1, 35, "saved quadrature count")?;
    let v = args[5];
    let center = args[6];
    if !(0. ..=1000.).contains(&v)
        || args.len() != 7 + 3 * p + 4 * groups + 2 * groups * q
        || data[0].len() != p + 1
    {
        return Err("GLMM saved model shape or variance is invalid".into());
    }
    work(data.len(), q + 35, groups)?;
    let origin = &args[7..7 + p];
    let offset = &args[7 + p..7 + 2 * p];
    let slopes = &args[7 + 2 * p..7 + 3 * p];
    let start = 7 + 3 * p;
    let ids = &args[start..start + groups];
    let means = &args[start + groups..start + 2 * groups];
    let modes = &args[start + 2 * groups..start + 3 * groups];
    let vars = &args[start + 3 * groups..start + 4 * groups];
    for (i, &id) in ids.iter().enumerate() {
        valid_id(id)?;
        if ids[..i].contains(&id) {
            return Err("saved group IDs must be unique".into());
        }
    }
    let node_start = start + 4 * groups;
    let weight_start = node_start + groups * q;
    for g in 0..groups {
        let nodes = &args[node_start + g * q..node_start + (g + 1) * q];
        let weights = &args[weight_start + g * q..weight_start + (g + 1) * q];
        if vars[g] < 0.
            || weights.iter().any(|w| *w < 0.)
            || (weights.iter().sum::<f64>() - 1.).abs() > 1e-8
        {
            return Err("saved posterior variance or weights are invalid".into());
        }
        let normalized = weights
            .iter()
            .map(|w| w / weights.iter().sum::<f64>())
            .collect::<Vec<_>>();
        let m = finite(dot(nodes, &normalized), "saved posterior mean")?;
        let variance = nodes
            .iter()
            .zip(&normalized)
            .map(|(x, w)| w * (x - m).powi(2))
            .sum::<f64>();
        finite(variance, "saved posterior second moment")?;
        if (m - means[g]).abs() > 1e-7 * (1. + m.abs())
            || (variance - vars[g]).abs() > 1e-7 * (1. + variance)
            || (v == 0.
                && (means[g] != 0.
                    || modes[g] != 0.
                    || vars[g] != 0.
                    || nodes.iter().any(|x| *x != 0.)))
        {
            return Err("saved posterior moments do not match quadrature".into());
        }
    }
    let prior = normal_rule(31)?;
    let prior_check = normal_rule(35)?;
    let mut conditional_predictions = Vec::new();
    let mut posterior_predictions = Vec::new();
    let mut population_predictions = Vec::new();
    let mut known = Vec::new();
    let mut maximum_delta = 0_f64;
    for row in data {
        valid_id(row[p])?;
        let eta = finite(
            center
                + (0..p)
                    .map(|j| slopes[j] * ((row[j] - origin[j]) - offset[j]))
                    .sum::<f64>(),
            "prediction linear predictor",
        )?;
        let group = ids.iter().position(|id| *id == row[p]);
        let link = |x: f64| -> Result<f64, String> {
            finite(
                if family == 0 { sigmoid(x) } else { x.exp() },
                "GLMM predicted mean",
            )
        };
        let population = if family == 1 {
            link(eta + v / 2.)?
        } else if v == 0. {
            sigmoid(eta)
        } else {
            let a = prior
                .nodes
                .iter()
                .zip(&prior.weights)
                .map(|(z, w)| w * sigmoid(eta + v.sqrt() * z))
                .sum::<f64>();
            let b = prior_check
                .nodes
                .iter()
                .zip(&prior_check.weights)
                .map(|(z, w)| w * sigmoid(eta + v.sqrt() * z))
                .sum::<f64>();
            maximum_delta = maximum_delta.max((a - b).abs());
            b
        };
        let conditional = link(eta + group.map_or(0., |g| means[g]))?;
        let posterior = if let Some(g) = group {
            let nodes = &args[node_start + g * q..node_start + (g + 1) * q];
            let weights = &args[weight_start + g * q..weight_start + (g + 1) * q];
            let total: f64 = weights.iter().sum();
            if family == 1 {
                // Combine in log space before exp: a tail node may overflow
                // individually while its weighted predictive contribution is finite.
                let logarithms = nodes
                    .iter()
                    .zip(weights)
                    .filter(|(_, w)| **w > 0.)
                    .map(|(b, w)| b + (w / total).ln())
                    .collect::<Vec<_>>();
                let maximum = logarithms.iter().copied().fold(-f64::INFINITY, f64::max);
                let sum: f64 = logarithms.iter().map(|x| (x - maximum).exp()).sum();
                link(eta + maximum + sum.ln())?
            } else {
                let prediction: f64 = nodes
                    .iter()
                    .zip(weights)
                    .map(|(b, w)| (w / total) * sigmoid(eta + b))
                    .sum();
                finite(prediction, "posterior predicted probability").map(|x| x.clamp(0., 1.))?
            }
        } else {
            population
        };
        conditional_predictions.push(conditional);
        posterior_predictions.push(posterior);
        population_predictions.push(population);
        known.push(group.is_some());
    }
    Ok(format!("{{\"version\":1,\"conditional_predictions\":{},\"posterior_predictions\":{},\"population_predictions\":{},\"known_groups\":[{}],\"prior_quadrature_orders\":[31,35],\"prior_quadrature_order_delta\":{maximum_delta},\"prior_quadrature_converged\":{},\"posterior_quadrature_convention\":\"saved_discrete_posterior_from_fit\"}}",numbers(&conditional_predictions),numbers(&posterior_predictions),numbers(&population_predictions),known.iter().map(bool::to_string).collect::<Vec<_>>().join(","),maximum_delta<=1e-8))
}

fn optional(x: Option<f64>) -> String {
    x.map_or("null".into(), |v| v.to_string())
}
fn sigmoid(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
fn softplus(x: f64) -> f64 {
    if x > 0. {
        x + (-x).exp().ln_1p()
    } else {
        x.exp().ln_1p()
    }
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

fn work(a: usize, b: usize, c: usize) -> Result<(), String> {
    if a.saturating_mul(b).saturating_mul(c) > MAX_WORK {
        Err("GLMM computation exceeds numerical work limit".into())
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
