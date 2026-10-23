//! Shared-inflation natural-support regression with explicit boundary candidates.
//! Gamma uses JMP's shape=mu/sigma,scale=sigma convention. Model-based and HC0
//! covariance are named native contracts; exact JMP inference is not implied.
use crate::quality::log_gamma as lanczos_log_gamma;
use std::ops::{Add, Div, Mul, Neg, Sub};

type Matrix = Vec<Vec<f64>>;
// Integral of DLMF5.7.4, used without subtracting nearly equal logGamma values.
fn log_gamma1p(x: f64) -> f64 {
    if x < 0.01 {
        let c = [
            -0.5772156649015329,
            0.8224670334241132,
            -0.40068563438653143,
            0.27058080842778454,
            -0.20738555102867398,
            0.1695571769974082,
            -0.1440498967688461,
            0.12550966952474304,
            -0.11133426586956469,
            0.1000994575127818,
        ];
        x * c.iter().rev().fold(0., |v, c| v * x + c)
    } else {
        lanczos_log_gamma(x + 1.)
    }
}
fn log_gamma(x: f64) -> f64 {
    if x < 0.5 {
        log_gamma1p(x) - x.ln()
    } else {
        lanczos_log_gamma(x)
    }
}
// log(x^a*exp(-x)/Gamma(a+1)). The normalized lower-tail series starts at1.
// DLMF5.11.1 supplies the Stirling correction; a local log1p deviance series
// avoids subtracting terms of order a*log(a) near the large-shape mode.
fn gamma_log_normalizer(a: f64, x: f64, log_x: f64) -> f64 {
    if x == f64::INFINITY {
        return f64::NEG_INFINITY;
    }
    if a < 16. {
        return -x + a * log_x - log_gamma1p(a);
    }
    let delta = (x - a) / a;
    let deviance = if delta.abs() < 0.5 {
        let mut term = -delta * delta / 2.;
        let mut sum = term;
        for k in 3..=80 {
            term *= -delta * (k - 1) as f64 / k as f64;
            sum += term;
            if term.abs() <= 2e-16 * sum.abs() {
                break;
            }
        }
        sum
    } else {
        log_x - a.ln() - x / a + 1.
    };
    let inv = 1. / a;
    let inv2 = inv * inv;
    let correction = inv
        * (1. / 12.
            - inv2
                * (1. / 360.
                    - inv2
                        * (1. / 1260.
                            - inv2 * (1. / 1680. - inv2 * (1. / 1188. - inv2 * 691. / 360360.)))));
    a * deviance - 0.5 * (2. * std::f64::consts::PI * a).ln() - correction
}
const MAX_WORK: usize = 200_000_000;
struct Budget {
    used: usize,
}
impl Budget {
    fn charge(&mut self, n: usize) -> Result<(), String> {
        self.used = self.used.saturating_add(n);
        if self.used > MAX_WORK {
            Err("zero-inflation numerical work exceeds 200000000 units".into())
        } else {
            Ok(())
        }
    }
}
fn finite(x: f64) -> Result<f64, String> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err("zero-inflation result is not finite".into())
    }
}
fn integer(x: f64, max: usize) -> Result<usize, String> {
    if x >= 0. && x <= max as f64 && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err("expected bounded nonnegative integer".into())
    }
}
fn dispersion_family(f: usize) -> bool {
    matches!(f, 1 | 3 | 4)
}
fn trial_family(f: usize) -> bool {
    f <= 1
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
    x.max(0.) + (-x.abs()).exp().ln_1p()
}
fn logadd(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY {
        return b;
    }
    if b == f64::NEG_INFINITY {
        return a;
    }
    let m = a.max(b);
    m + ((a - m).exp() + (b - m).exp()).ln()
}
fn vector(a: &[f64]) -> String {
    format!(
        "[{}]",
        a.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(a: &Matrix) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| vector(r)).collect::<Vec<_>>().join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}

fn log_choose(n: usize, y: usize) -> f64 {
    log_gamma(n as f64 + 1.) - log_gamma(y as f64 + 1.) - log_gamma((n - y) as f64 + 1.)
}
fn base_log(
    f: usize,
    eta: f64,
    disp: f64,
    n: usize,
    y: f64,
    budget: &mut Budget,
) -> Result<f64, String> {
    budget.charge(32)?;
    if y < 0. || (f != 4 && (y.fract() != 0. || y > 1e6)) || (trial_family(f) && y > n as f64) {
        return Ok(f64::NEG_INFINITY);
    }
    let mu = if trial_family(f) {
        sigmoid(eta)
    } else {
        eta.exp()
    };
    if mu <= 0. || !mu.is_finite() {
        return Err("unrepresentable base mean".into());
    }
    if f != 4 && !trial_family(f) && mu > 1e6 {
        return Err("count mean exceeds 1000000".into());
    }
    let k = y as usize;
    let out = match f {
        0 => log_choose(n, k) - y * softplus(-eta) - (n as f64 - y) * softplus(eta),
        1 if disp == 0. => return base_log(0, eta, 0., n, y, budget),
        1 => {
            budget.charge(n.saturating_mul(8))?;
            let mut l = log_choose(n, k);
            for i in 0..k {
                l += (mu * (1. - disp) + i as f64 * disp).ln();
            }
            for i in 0..n - k {
                l += (sigmoid(-eta) * (1. - disp) + i as f64 * disp).ln();
            }
            for i in 0..n {
                l -= (1. - disp + i as f64 * disp).ln();
            }
            l
        }
        2 => -mu + y * eta - log_gamma(y + 1.),
        3 if disp == 0. => return base_log(2, eta, 0., n, y, budget),
        3 => {
            budget.charge(k.saturating_mul(8))?;
            let z = mu * disp;
            let ratio = if z.abs() < 1e-4 {
                1. - z / 2. + z * z / 3. - z * z * z / 4. + z.powi(4) / 5.
            } else {
                z.ln_1p() / z
            };
            let mut l = y * eta - log_gamma(y + 1.) - y * z.ln_1p() - mu * ratio;
            for i in 0..k {
                l += (i as f64 * disp).ln_1p();
            }
            l
        }
        4 => {
            if y == 0. {
                return Ok(f64::NEG_INFINITY);
            }
            let shape = mu / disp;
            if !(1e-8..=1e7).contains(&shape) {
                return Err("Gamma shape outside numerical range".into());
            }
            let value =
                gamma_log_normalizer(shape, y / disp, y.ln() - disp.ln()) + shape.ln() - y.ln();
            if value == f64::NEG_INFINITY {
                return Ok(value);
            }
            value
        }
        _ => return Err("unknown zero-inflation family".into()),
    };
    finite(out)
}

fn gamma_tails(a: f64, x: f64, log_x: f64, budget: &mut Budget) -> Result<(f64, f64), String> {
    if !(1e-8..=1e7).contains(&a) {
        return Err("Gamma tail shape outside numerical range".into());
    }
    if log_x == f64::NEG_INFINITY {
        return Ok((0., 1.));
    }
    if x == f64::INFINITY {
        return Ok((1., 0.));
    }
    // Preserve a*log(x) even when the positive quotient y/scale underflows.
    budget.charge(80)?;
    if a < 0.5 && x < 1. {
        // DLMF8.7.1's alternating expansion. Sum the correction separately
        // from1 and use log1p, preserving a tiny upper tail when a is small.
        let mut factorial_term = 1.;
        let mut correction = 0.;
        for n in 1..=10000 {
            budget.charge(24)?;
            factorial_term *= -x / n as f64;
            let term = factorial_term * a / (a + n as f64);
            correction += term;
            if term.abs() <= 2e-16 * correction.abs() {
                let lp = a * log_x - log_gamma1p(a) + correction.ln_1p();
                return Ok((lp.exp(), -lp.exp_m1()));
            }
        }
        return Err("small-shape Gamma tail iteration limit".into());
    }
    let lead = gamma_log_normalizer(a, x, log_x);
    if x < a + 1. {
        let mut term = 1.;
        let mut sum = term;
        for i in 1..=100000 {
            budget.charge(16)?;
            term *= x / (a + i as f64);
            sum += term;
            if term.abs() <= sum.abs() * 2e-15 {
                let lp = lead + sum.ln();
                return Ok((lp.exp().min(1.), -lp.exp_m1()));
            }
        }
    } else {
        let tiny = 1e-300;
        let mut b = x + 1. - a;
        let mut c = 1. / tiny;
        let mut d = 1. / b;
        let mut h = d;
        for i in 1..=100000 {
            budget.charge(24)?;
            let an = -(i as f64) * (i as f64 - a);
            b += 2.;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny.copysign(d);
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny.copysign(c);
            }
            d = 1. / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.).abs() <= 2e-15 {
                let lq = lead + a.ln() + h.ln();
                return Ok((-lq.exp_m1(), lq.exp().min(1.)));
            }
        }
    }
    Err("Gamma tail iteration limit".into())
}
fn recurrence(f: usize, eta: f64, disp: f64, n: usize, k: usize) -> f64 {
    let mu = if trial_family(f) {
        sigmoid(eta)
    } else {
        eta.exp()
    };
    let k = k as f64;
    match f {
        0 => (n as f64 - k).ln() - (k + 1.).ln() + eta,
        1 if disp == 0. => recurrence(0, eta, 0., n, k as usize),
        1 => {
            (n as f64 - k).ln() - (k + 1.).ln() + (mu * (1. - disp) + k * disp).ln()
                - (sigmoid(-eta) * (1. - disp) + (n as f64 - k - 1.) * disp).ln()
        }
        2 => mu.ln() - (k + 1.).ln(),
        3 if disp == 0. => recurrence(2, eta, 0., n, k as usize),
        3 => (1. + k * disp).ln() + mu.ln() - (k + 1.).ln() - (mu * disp).ln_1p(),
        _ => unreachable!(),
    }
}
fn count_tails(
    f: usize,
    eta: f64,
    disp: f64,
    n: usize,
    k: usize,
    budget: &mut Budget,
) -> Result<(f64, f64), String> {
    if trial_family(f) && k >= n {
        return Ok((1., 0.));
    }
    let mu = if trial_family(f) {
        sigmoid(eta)
    } else {
        eta.exp()
    };
    let mean = if trial_family(f) { n as f64 * mu } else { mu };
    if (k as f64) < mean {
        let mut term = base_log(f, eta, disp, n, 0., budget)?;
        let mut sum = term;
        for i in 0..k {
            budget.charge(16)?;
            term += recurrence(f, eta, disp, n, i);
            sum = logadd(sum, term);
        }
        let l = sum.min(0.);
        Ok((l.exp(), -l.exp_m1()))
    } else if trial_family(f) {
        let mut term = base_log(f, eta, disp, n, n as f64, budget)?;
        let mut sum = term;
        for i in (k + 1..n).rev() {
            budget.charge(16)?;
            term -= recurrence(f, eta, disp, n, i);
            sum = logadd(sum, term);
        }
        let l = sum.min(0.);
        Ok((-l.exp_m1(), l.exp()))
    } else {
        let mut term = base_log(f, eta, disp, n, (k + 1) as f64, budget)?;
        let mut sum = term;
        for i in k + 1..1_000_000 {
            budget.charge(16)?;
            term += recurrence(f, eta, disp, n, i);
            sum = logadd(sum, term);
            if term - sum < -36. && recurrence(f, eta, disp, n, i) < 0. {
                let l = sum.min(0.);
                return Ok((-l.exp_m1(), l.exp()));
            }
        }
        Err("count tail iteration limit".into())
    }
}
fn tails(
    f: usize,
    eta: f64,
    disp: f64,
    pi: f64,
    n: usize,
    y: f64,
    budget: &mut Budget,
) -> Result<(f64, f64), String> {
    if y < 0. {
        return Ok((0., 1.));
    }
    let (lo, hi) = if f == 4 {
        gamma_tails(eta.exp() / disp, y / disp, y.ln() - disp.ln(), budget)?
    } else {
        if y > 1e6 && !(trial_family(f) && y >= n as f64) {
            return Err("count CDF query exceeds 1000000".into());
        }
        count_tails(f, eta, disp, n, y.floor() as usize, budget)?
    };
    Ok((pi + (1. - pi) * lo, (1. - pi) * hi))
}
fn quantile(
    f: usize,
    eta: f64,
    disp: f64,
    pi: f64,
    n: usize,
    u: f64,
    budget: &mut Budget,
) -> Result<f64, String> {
    if u <= pi {
        return Ok(0.);
    }
    // Compare the base distribution's two conditional tails separately. Adding
    // pi first loses half an ulp of probability immediately above the zero atom.
    let lower = (u - pi) / (1. - pi);
    let upper = (1. - u) / (1. - pi);
    let reached = |value, budget: &mut Budget| -> Result<bool, String> {
        let (cdf, sf) = tails(f, eta, disp, 0., n, value, budget)?;
        Ok(if lower <= upper {
            cdf >= lower
        } else {
            sf <= upper
        })
    };
    if reached(0., budget)? {
        return Ok(0.);
    }
    if f == 4 {
        let smallest = f64::from_bits(1);
        let smallest_cdf = tails(f, eta, disp, 0., n, smallest, budget)?.0;
        if smallest_cdf > lower {
            return Err("positive Gamma quantile is below the smallest representable f64".into());
        }
        if smallest_cdf == lower {
            return Ok(smallest);
        }
        let mut hi = eta.exp().max(disp);
        for _ in 0..1024 {
            if reached(hi, budget)? {
                break;
            }
            hi = finite(hi * 2.)?;
        }
        if !reached(hi, budget)? {
            return Err("Gamma quantile not bracketed".into());
        }
        let mut lo = smallest.ln();
        hi = hi.ln();
        for _ in 0..180 {
            let mid = lo + (hi - lo) / 2.;
            if reached(mid.exp(), budget)? {
                hi = mid;
            } else {
                lo = mid;
            }
            if hi - lo <= 1e-12 {
                return finite(hi.exp());
            }
        }
        Err("Gamma quantile failed its log-domain convergence criterion".into())
    } else {
        let mut lo = 0usize;
        let mut hi = if trial_family(f) {
            n
        } else {
            (eta.exp() + 16.).ceil().min(1e6) as usize
        };
        while !reached(hi as f64, budget)? {
            if hi >= 1_000_000 {
                return Err("count quantile exceeds 1000000".into());
            }
            hi = (hi * 2).min(1_000_000);
        }
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            if reached(mid as f64, budget)? {
                hi = mid;
            } else {
                lo = mid + 1;
            }
        }
        Ok(lo as f64)
    }
}

// Triangular solves and symmetry averaging need cross-row matrix indexing.
#[allow(clippy::needless_range_loop)]
fn inverse_spd(a: &Matrix) -> Option<Matrix> {
    let n = a.len();
    let s = a
        .iter()
        .enumerate()
        .map(|(i, r)| r[i].sqrt())
        .collect::<Vec<_>>();
    if s.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return None;
    }
    let mut l = vec![vec![0.; n]; n];
    for i in 0..n {
        for j in 0..=i {
            if !a[i][j].is_finite() || (a[i][j] - a[j][i]).abs() > 1e-10 * s[i] * s[j] {
                return None;
            }
            let v = a[i][j] / s[i] / s[j] - dot(&l[i][..j], &l[j][..j]);
            if i == j {
                if v <= 1e-12 {
                    return None;
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut out = vec![vec![0.; n]; n];
    for k in 0..n {
        let mut x = vec![0.; n];
        for i in 0..n {
            x[i] = ((if i == k { 1. } else { 0. }) - dot(&l[i][..i], &x[..i])) / l[i][i];
        }
        for i in (0..n).rev() {
            x[i] = (x[i] - (i + 1..n).map(|j| l[j][i] * x[j]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            out[i][k] = x[i] / s[i] / s[k];
            if !out[i][k].is_finite() {
                return None;
            }
        }
    }
    for i in 0..n {
        for j in 0..i {
            let x = (out[i][j] + out[j][i]) / 2.;
            out[i][j] = x;
            out[j][i] = x;
        }
    }
    Some(out)
}
fn normal_quantile(p: f64) -> f64 {
    let a = [
        -39.69683028665376,
        220.9460984245205,
        -275.9285104469687,
        138.357751867269,
        -30.66479806614716,
        2.506628277459239,
    ];
    let b = [
        -54.47609879822406,
        161.5858368580409,
        -155.6989798598866,
        66.80131188771972,
        -13.28068155288572,
    ];
    let c = [
        -0.007784894002430293,
        -0.3223964580411365,
        -2.400758277161838,
        -2.549732539343734,
        4.374664141464968,
        2.938163982698783,
    ];
    let d = [
        0.007784695709041462,
        0.3224671290700398,
        2.445134137142996,
        3.754408661907416,
    ];
    if !(0.02425..=0.97575).contains(&p) {
        let q = (-2. * p.min(1. - p).ln()).sqrt();
        let x = (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5])
            / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.);
        if p > 0.5 {
            -x
        } else {
            x
        }
    } else {
        let q = p - 0.5;
        let r = q * q;
        (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q
            / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.)
    }
}
fn predict(data: &[Vec<f64>], a: &[f64], budget: &mut Budget) -> Result<String, String> {
    if a.len() < 10 || a[0] != 1. {
        return Err("invalid zero-inflation carrier version/length".into());
    }
    let f = integer(a[1], 4)?;
    let p = integer(a[2], 16)?;
    let dim = p + 2 + usize::from(dispersion_family(f));
    if a.len() != 10 + 4 * p + dim * dim || data[0].len() != p + 3 || !(0. < a[3] && a[3] < 1.) {
        return Err("invalid zero-inflation carrier dimensions/confidence".into());
    }
    let kind = integer(a[4], 1)?;
    let has = integer(a[5], 1)?;
    let boundary = integer(a[6], 1)?;
    let beta = &a[7..8 + p];
    let disp = a[8 + p];
    let pi = a[9 + p];
    if !(0. ..1.).contains(&pi)
        || (!dispersion_family(f) && disp != 0.)
        || (f == 1 && !(0. ..1.).contains(&disp))
        || (f == 3 && disp < 0.)
        || (f == 4 && disp <= 0.)
    {
        return Err("invalid dispersion/inflation carrier parameters".into());
    }
    if (pi == 0. || (matches!(f, 1 | 3) && disp == 0.)) && boundary != 1 {
        return Err("boundary flag inconsistent with model".into());
    }
    if has == 1 && boundary == 1 {
        return Err("boundary model cannot claim joint Wald covariance".into());
    }
    let origins = &a[10 + p..10 + 2 * p];
    let offsets = &a[10 + 2 * p..10 + 3 * p];
    let scales = &a[10 + 3 * p..10 + 4 * p];
    if scales.iter().any(|x| *x <= 0.) {
        return Err("saved feature scales must be positive".into());
    }
    let covariance = a[10 + 4 * p..]
        .chunks_exact(dim)
        .map(|r| r.to_vec())
        .collect::<Matrix>();
    budget.charge(dim.pow(3))?;
    if has == 1 {
        if inverse_spd(&covariance).is_none() {
            return Err("saved covariance must be symmetric positive definite".into());
        }
    } else if covariance.iter().flatten().any(|x| *x != 0.) {
        return Err("disabled covariance must be all zero".into());
    }
    let names = [
        "linear_predictor",
        "conditional_base_mean",
        "mixture_mean",
        "mixture_variance",
        "zero_probability",
        "query_probability_mass",
        "query_positive_density",
        "cdf",
        "survival",
        "quantile",
    ];
    let mut out = vec![Vec::new(); names.len()];
    let mut uncertainty = vec![];
    let mut intervals = vec![];
    let z = normal_quantile((1. + a[3]) / 2.);
    for row in data {
        let trials = integer(row[p], 1_000_000)?;
        if (trial_family(f) && (trials == 0 || (f == 1 && trials <= 1)))
            || (!trial_family(f) && trials != 1)
        {
            return Err("invalid prediction trials".into());
        }
        let y = row[p + 1];
        let u = row[p + 2];
        if !(0. < u && u < 1.) {
            return Err("quantile probability must be strictly between zero and one".into());
        }
        let mut x = vec![1.];
        for j in 0..p {
            x.push(finite(((row[j] - origins[j]) - offsets[j]) / scales[j])?);
        }
        let eta = finite(dot(beta, &x))?;
        let linkmean = if trial_family(f) {
            sigmoid(eta)
        } else {
            eta.exp()
        };
        if !linkmean.is_finite() || linkmean <= 0. || (!trial_family(f) && f != 4 && linkmean > 1e6)
        {
            return Err("prediction mean outside numerical range".into());
        }
        let mean = if trial_family(f) {
            trials as f64 * linkmean
        } else {
            linkmean
        };
        let var = match f {
            0 => mean * sigmoid(-eta),
            1 => mean * sigmoid(-eta) * (1. + (trials as f64 - 1.) * disp),
            2 => mean,
            3 => mean * (1. + mean * disp),
            4 => mean * disp,
            _ => unreachable!(),
        };
        let mixmean = finite((1. - pi) * mean)?;
        let mixvar = finite((1. - pi) * var + pi * (1. - pi) * mean * mean)?;
        let g0 = if f == 4 {
            0.
        } else {
            base_log(f, eta, disp, trials, 0., budget)?.exp()
        };
        let zero = pi + (1. - pi) * g0;
        let mut pmf = 0.;
        let mut density = 0.;
        if f == 4 {
            if y == 0. {
                pmf = pi;
            } else if y > 0. {
                density = finite((1. - pi) * base_log(f, eta, disp, trials, y, budget)?.exp())?;
            }
        } else if y == 0. {
            pmf = zero;
        } else {
            pmf = (1. - pi) * base_log(f, eta, disp, trials, y, budget)?.exp();
        }
        let (cdf, sf) = tails(f, eta, disp, pi, trials, y, budget)?;
        let q = quantile(f, eta, disp, pi, trials, u, budget)?;
        for (v, x) in out
            .iter_mut()
            .zip([eta, mean, mixmean, mixvar, zero, pmf, density, cdf, sf, q])
        {
            v.push(finite(x)?);
        }
        if has == 1 {
            let deriv = if trial_family(f) {
                mean * sigmoid(-eta)
            } else {
                mean
            };
            let mut g = x.iter().map(|x| (1. - pi) * deriv * x).collect::<Vec<_>>();
            if dispersion_family(f) {
                g.push(0.);
            }
            g.push(-mean);
            let v = finite(
                (0..dim)
                    .map(|i| {
                        (0..dim)
                            .map(|j| g[i] * covariance[i][j] * g[j])
                            .sum::<f64>()
                    })
                    .sum(),
            )?;
            if v < 0. {
                return Err("negative mean delta variance".into());
            }
            let se = v.sqrt();
            uncertainty.push(se.to_string());
            intervals.push(format!(
                "[{},{}]",
                finite(mixmean - z * se)?,
                finite(mixmean + z * se)?
            ));
        } else {
            uncertainty.push("null".into());
            intervals.push("null".into());
        }
    }
    let fields = names
        .iter()
        .zip(out.iter())
        .map(|(k, v)| format!("\"{k}\":{}", vector(v)))
        .collect::<Vec<_>>()
        .join(",");
    Ok(format!("{{\"version\":1,\"algorithm\":\"zi_predict\",{fields},\"mean_standard_error\":[{}],\"mean_confidence_interval\":[{}],\"covariance_method\":\"{}\",\"at_boundary\":{},\"jmp_equivalence_verified\":false,\"work_units\":{}}}",uncertainty.join(","),intervals.join(","),if kind==0{"observed_information"}else{"sandwich_HC0"},boundary==1,budget.used))
}

// Three local coordinates: linear predictor,transformed dispersion,and inflation.
// The full design Hessian is assembled from these row derivatives.
#[derive(Clone, Copy)]
struct Jet {
    v: f64,
    g: [f64; 3],
    h: [[f64; 3]; 3],
}
impl Jet {
    fn c(v: f64) -> Self {
        Self {
            v,
            g: [0.; 3],
            h: [[0.; 3]; 3],
        }
    }
    fn x(v: f64, i: usize) -> Self {
        let mut j = Self::c(v);
        j.g[i] = 1.;
        j
    }
    fn chain(self, v: f64, d: f64, dd: f64) -> Self {
        let mut o = Self::c(v);
        for i in 0..3 {
            o.g[i] = d * self.g[i];
            for j in 0..3 {
                o.h[i][j] = d * self.h[i][j] + dd * self.g[i] * self.g[j];
            }
        }
        o
    }
    fn exp(self) -> Self {
        let e = self.v.exp();
        self.chain(e, e, e)
    }
    fn ln(self) -> Self {
        self.chain(self.v.ln(), 1. / self.v, -1. / (self.v * self.v))
    }
    fn log1p(self) -> Self {
        self.chain(
            self.v.ln_1p(),
            1. / (1. + self.v),
            -1. / (1. + self.v).powi(2),
        )
    }
    fn logistic(self) -> Self {
        let p = sigmoid(self.v);
        let q = sigmoid(-self.v);
        self.chain(p, p * q, p * q * (q - p))
    }
    fn softplus(self) -> Self {
        let p = sigmoid(self.v);
        self.chain(softplus(self.v), p, p * sigmoid(-self.v))
    }
    fn lgamma(self) -> Self {
        let mut x = self.v;
        let mut d = 0.;
        let mut dd = 0.;
        while x < 12. {
            d -= 1. / x;
            dd += 1. / (x * x);
            x += 1.;
        }
        let z = 1. / x;
        let z2 = z * z;
        d += x.ln()
            - z / 2.
            - z2 * (1. / 12. - z2 * (1. / 120. - z2 * (1. / 252. - z2 * (1. / 240. - z2 / 132.))));
        dd += z
            + z2 / 2.
            + z * z2
                * (1. / 6. - z2 * (1. / 30. - z2 * (1. / 42. - z2 * (1. / 30. - z2 * 5. / 66.))));
        self.chain(log_gamma(self.v), d, dd)
    }
}
impl Add for Jet {
    type Output = Self;
    fn add(self, b: Self) -> Self {
        let mut o = Self::c(self.v + b.v);
        for i in 0..3 {
            o.g[i] = self.g[i] + b.g[i];
            for j in 0..3 {
                o.h[i][j] = self.h[i][j] + b.h[i][j];
            }
        }
        o
    }
}
impl Neg for Jet {
    type Output = Self;
    fn neg(self) -> Self {
        self * -1.
    }
}
impl Sub for Jet {
    type Output = Self;
    fn sub(self, b: Self) -> Self {
        self + -b
    }
}
impl Mul for Jet {
    type Output = Self;
    fn mul(self, b: Self) -> Self {
        let mut o = Self::c(self.v * b.v);
        for i in 0..3 {
            o.g[i] = self.g[i] * b.v + self.v * b.g[i];
            for j in 0..3 {
                o.h[i][j] = self.h[i][j] * b.v
                    + self.v * b.h[i][j]
                    + self.g[i] * b.g[j]
                    + b.g[i] * self.g[j];
            }
        }
        o
    }
}
impl Div for Jet {
    type Output = Self;
    fn div(self, b: Self) -> Self {
        self * b.chain(1. / b.v, -1. / b.v.powi(2), 2. / b.v.powi(3))
    }
}
impl Add<f64> for Jet {
    type Output = Self;
    fn add(self, b: f64) -> Self {
        self + Self::c(b)
    }
}
impl Sub<f64> for Jet {
    type Output = Self;
    fn sub(self, b: f64) -> Self {
        self - Self::c(b)
    }
}
impl Mul<f64> for Jet {
    type Output = Self;
    fn mul(self, b: f64) -> Self {
        let mut o = self;
        o.v *= b;
        for i in 0..3 {
            o.g[i] *= b;
            for j in 0..3 {
                o.h[i][j] *= b;
            }
        }
        o
    }
}
impl Div<f64> for Jet {
    type Output = Self;
    fn div(self, b: f64) -> Self {
        self * (1. / b)
    }
}
fn jet_logadd(a: Jet, b: Jet) -> Jet {
    let m = a.v.max(b.v);
    ((a - m).exp() + (b - m).exp()).ln() + m
}

struct Frame {
    family: usize,
    p: usize,
    x: Matrix,
    y: Vec<f64>,
    trials: Vec<usize>,
    weights: Vec<f64>,
    total: f64,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
}
impl Frame {
    fn new(family: usize, data: &[Vec<f64>]) -> Result<Self, String> {
        let p = data[0]
            .len()
            .checked_sub(3)
            .ok_or("fit rows need [features...,response,trials,frequency]")?;
        if p > 16 {
            return Err("zero-inflation fit supports at most16features".into());
        }
        let mut rows = vec![];
        let mut y = vec![];
        let mut trials = vec![];
        let mut weights = vec![];
        let mut total = 0.;
        for row in data {
            let w = integer(row[p + 2], 1_000_000_000)? as f64;
            if w == 0. {
                continue;
            }
            let n = integer(row[p + 1], 1_000_000)?;
            let v = row[p];
            if (trial_family(family) && (n == 0 || (family == 1 && n < 2)))
                || (!trial_family(family) && n != 1)
            {
                return Err("invalid fit trials".into());
            }
            if v < 0.
                || (family != 4 && (v.fract() != 0. || v > 1e6))
                || (trial_family(family) && v > n as f64)
            {
                return Err("response outside natural family support".into());
            }
            total += w;
            if total > 1e9 {
                return Err("total frequency exceeds1000000000".into());
            }
            rows.push(row);
            y.push(v);
            trials.push(n);
            weights.push(w);
        }
        let d = p + 2 + usize::from(dispersion_family(family));
        if rows.is_empty() || total <= d as f64 || y.iter().all(|x| *x == 0.) {
            return Err(
                "insufficient effective observations or all-zero unidentified model".into(),
            );
        }
        if trial_family(family) && y.iter().zip(&trials).all(|(y, n)| *y == *n as f64) {
            return Err("all-success binomial has no finite base MLE".into());
        }
        if family == 0 && p == 0 && trials.iter().all(|n| *n == 1) {
            return Err(
                "intercept-only Bernoulli inflation and base probability are confounded".into(),
            );
        }
        if family == 4 {
            let positive = y.iter().copied().filter(|x| *x > 0.).collect::<Vec<_>>();
            if positive.iter().all(|x| *x == positive[0]) {
                return Err(
                    "constant positive Gamma responses have an unbounded scale-zero likelihood"
                        .into(),
                );
            }
        }
        let origins = rows[0][..p].to_vec();
        let offsets = (0..p)
            .map(|j| {
                finite(
                    rows.iter()
                        .zip(&weights)
                        .map(|(r, w)| (r[j] - origins[j]) * (w / total))
                        .sum(),
                )
            })
            .collect::<Result<Vec<_>, _>>()?;
        let scales = (0..p)
            .map(|j| {
                let s = rows.iter().zip(&weights).fold(0_f64, |s, (r, w)| {
                    s.hypot(((r[j] - origins[j]) - offsets[j]) * (w / total).sqrt())
                });
                if s > 0. && s.is_finite() {
                    Ok(s)
                } else {
                    Err("constant or unrepresentable feature scale".to_string())
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let x = rows
            .iter()
            .map(|r| {
                let mut x = vec![1.];
                for j in 0..p {
                    x.push(((r[j] - origins[j]) - offsets[j]) / scales[j]);
                }
                x
            })
            .collect::<Matrix>();
        let mut q: Matrix = vec![];
        for j in 0..=p {
            let mut c = x
                .iter()
                .zip(&weights)
                .map(|(r, w)| r[j] * (w / total).sqrt())
                .collect::<Vec<_>>();
            for _ in 0..2 {
                for a in &q {
                    let v = dot(a, &c);
                    for (c, a) in c.iter_mut().zip(a) {
                        *c -= v * a;
                    }
                }
            }
            let n = c.iter().fold(0_f64, |n, x| n.hypot(*x));
            if n < 1e-10 {
                return Err("rank-deficient fixed feature design".into());
            }
            q.push(c.iter().map(|x| x / n).collect());
        }
        Ok(Self {
            family,
            p,
            x,
            y,
            trials,
            weights,
            total,
            origins,
            offsets,
            scales,
        })
    }
}
#[derive(Clone, Copy)]
struct Mode {
    disp: bool,
    pi: bool,
}
struct Eval {
    value: f64,
    g: Vec<f64>,
    h: Matrix,
    meat: Matrix,
}
fn row_log(
    f: usize,
    eta: Jet,
    disp: Jet,
    pi: Jet,
    n: usize,
    y: f64,
    budget: &mut Budget,
) -> Result<Jet, String> {
    budget.charge(96)?;
    if f == 4 && y == 0. {
        if pi.v == 0. {
            return Err("zero Gamma observation needs positive inflation".into());
        }
        return Ok(pi.ln());
    }
    let mu = if trial_family(f) {
        eta.logistic()
    } else {
        eta.exp()
    };
    if !mu.v.is_finite() || mu.v <= 0. || (!trial_family(f) && f != 4 && mu.v > 1e6) {
        return Err("fit mean outside numerical range".into());
    }
    let k = y as usize;
    let base = match f {
        0 => Jet::c(log_choose(n, k)) - (-eta).softplus() * y - eta.softplus() * (n as f64 - y),
        1 if disp.v == 0. => return row_log(0, eta, disp, pi, n, y, budget),
        1 => {
            budget.charge(n.saturating_mul(160))?;
            let mut l = Jet::c(log_choose(n, k));
            for i in 0..k {
                l = l + (mu * (Jet::c(1.) - disp) + disp * i as f64).ln();
            }
            for i in 0..n - k {
                l = l + ((-eta).logistic() * (Jet::c(1.) - disp) + disp * i as f64).ln();
            }
            for i in 0..n {
                l = l - (Jet::c(1.) - disp + disp * i as f64).ln();
            }
            l
        }
        2 => -mu + eta * y - Jet::c(log_gamma(y + 1.)),
        3 if disp.v == 0. => return row_log(2, eta, disp, pi, n, y, budget),
        3 => {
            budget.charge(k.saturating_mul(80))?;
            let z = mu * disp;
            let ratio = if z.v.abs() < 1e-4 {
                Jet::c(1.) - z / 2. + z * z / 3. - z * z * z / 4. + z * z * z * z / 5.
            } else {
                z.log1p() / z
            };
            let mut l = eta * y - Jet::c(log_gamma(y + 1.)) - z.log1p() * y - mu * ratio;
            for i in 0..k {
                l = l + (disp * i as f64).log1p();
            }
            l
        }
        4 => {
            let shape = mu / disp;
            if !(1e-8..=1e7).contains(&shape.v) {
                return Err("Gamma fit shape outside numerical range".into());
            }
            let mut result =
                (shape - 1.) * y.ln() - Jet::c(y) / disp - shape.lgamma() - shape * disp.ln();
            result.v = gamma_log_normalizer(shape.v, y / disp.v, y.ln() - disp.v.ln())
                + shape.v.ln()
                - y.ln();
            result
        }
        _ => unreachable!(),
    };
    if pi.v == 0. {
        return Ok(base);
    }
    let ordinary = (Jet::c(1.) - pi).ln() + base;
    Ok(if y == 0. {
        jet_logadd(pi.ln(), ordinary)
    } else {
        ordinary
    })
}
fn parameters(frame: &Frame, t: &[f64], mode: Mode) -> (Jet, Jet) {
    let disp = if mode.disp {
        let d = Jet::x(t[frame.p + 1], 1);
        if frame.family == 1 {
            d.logistic()
        } else {
            d.exp()
        }
    } else {
        Jet::c(0.)
    };
    let pi = if mode.pi {
        Jet::x(t[frame.p + 1 + usize::from(mode.disp)], 2).logistic()
    } else {
        Jet::c(0.)
    };
    (disp, pi)
}
fn evaluate(
    frame: &Frame,
    t: &[f64],
    mode: Mode,
    budget: &mut Budget,
    meat: bool,
) -> Result<Eval, String> {
    if t.iter().any(|v| !v.is_finite() || v.abs() > 60.) {
        return Err("optimizer parameter guard".into());
    }
    let d = t.len();
    let (disp, pi) = parameters(frame, t, mode);
    let mut e = Eval {
        value: 0.,
        g: vec![0.; d],
        h: vec![vec![0.; d]; d],
        meat: vec![vec![0.; d]; d],
    };
    for i in 0..frame.y.len() {
        budget.charge(d * d * 4)?;
        let eta = Jet::x(dot(&t[..frame.p + 1], &frame.x[i]), 0);
        let j = row_log(
            frame.family,
            eta,
            disp,
            pi,
            frame.trials[i],
            frame.y[i],
            budget,
        )?;
        if !j.v.is_finite()
            || j.g.iter().any(|x| !x.is_finite())
            || j.h.iter().flatten().any(|x| !x.is_finite())
        {
            return Err("nonfinite likelihood derivatives".into());
        }
        let w = frame.weights[i] / frame.total;
        e.value -= w * j.v;
        for a in 0..d {
            let (ca, fa) = if a <= frame.p {
                (0, frame.x[i][a])
            } else if mode.disp && a == frame.p + 1 {
                (1, 1.)
            } else {
                (2, 1.)
            };
            e.g[a] -= w * j.g[ca] * fa;
            for b in 0..d {
                let (cb, fb) = if b <= frame.p {
                    (0, frame.x[i][b])
                } else if mode.disp && b == frame.p + 1 {
                    (1, 1.)
                } else {
                    (2, 1.)
                };
                e.h[a][b] -= w * j.h[ca][cb] * fa * fb;
                if meat {
                    e.meat[a][b] += frame.weights[i] * j.g[ca] * j.g[cb] * fa * fb;
                }
            }
        }
    }
    finite(e.value)?;
    Ok(e)
}
struct Opt {
    t: Vec<f64>,
    mode: Mode,
    e: Eval,
    iterations: usize,
    converged: bool,
    line_failed: bool,
}
fn optimize(
    frame: &Frame,
    mut t: Vec<f64>,
    mode: Mode,
    maxiter: usize,
    tol: f64,
    budget: &mut Budget,
) -> Result<Opt, String> {
    let d = t.len();
    let mut h = vec![vec![0.; d]; d];
    for (i, r) in h.iter_mut().enumerate() {
        r[i] = 1.;
    }
    let mut e = evaluate(frame, &t, mode, budget, false)?;
    let mut iterations = 0;
    let mut line_failed = false;
    for _ in 0..maxiter {
        if e.g.iter().fold(0_f64, |m, x| m.max(x.abs())) <= tol {
            break;
        }
        iterations += 1;
        budget.charge(d.pow(3))?;
        let inv = inverse_spd(&e.h).unwrap_or_else(|| h.clone());
        let mut step = inv.iter().map(|r| -dot(r, &e.g)).collect::<Vec<_>>();
        if dot(&step, &e.g) >= 0. {
            step = e.g.iter().map(|x| -x).collect();
        }
        let norm = step.iter().fold(0_f64, |m, x| m.max(x.abs()));
        if norm > 10. {
            for v in &mut step {
                *v *= 10. / norm;
            }
        }
        let descent = dot(&step, &e.g);
        let mut alpha = 1.;
        let mut accepted = None;
        for _ in 0..40 {
            let candidate = t
                .iter()
                .zip(&step)
                .map(|(x, s)| x + alpha * s)
                .collect::<Vec<_>>();
            match evaluate(frame, &candidate, mode, budget, false) {
                Ok(ne) if ne.value <= e.value + 1e-4 * alpha * descent => {
                    accepted = Some((candidate, ne));
                    break;
                }
                Err(err) if budget.used > MAX_WORK => return Err(err),
                _ => alpha *= 0.5,
            }
        }
        let Some((next, ne)) = accepted else {
            line_failed = true;
            break;
        };
        let s = next.iter().zip(&t).map(|(a, b)| a - b).collect::<Vec<_>>();
        let y =
            ne.g.iter()
                .zip(&e.g)
                .map(|(a, b)| a - b)
                .collect::<Vec<_>>();
        let sy = dot(&s, &y);
        if sy > 1e-14 {
            let hy = h.iter().map(|r| dot(r, &y)).collect::<Vec<_>>();
            let yhy = dot(&y, &hy);
            for i in 0..d {
                for j in 0..d {
                    h[i][j] +=
                        (1. + yhy / sy) * s[i] * s[j] / sy - (s[i] * hy[j] + hy[i] * s[j]) / sy;
                }
            }
        }
        t = next;
        e = ne;
    }
    let converged = e.g.iter().fold(0_f64, |m, x| m.max(x.abs())) <= tol
        && !line_failed
        && t.iter().all(|x| x.abs() < 40.);
    Ok(Opt {
        t,
        mode,
        e,
        iterations,
        converged,
        line_failed,
    })
}
fn transform_cov(c: &Matrix, j: &Matrix) -> Matrix {
    let n = j.len();
    (0..n)
        .map(|a| {
            (0..n)
                .map(|b| {
                    (0..n)
                        .map(|i| (0..n).map(|k| j[a][i] * c[i][k] * j[b][k]).sum::<f64>())
                        .sum()
                })
                .collect()
        })
        .collect()
}
fn fit(
    family: usize,
    data: &[Vec<f64>],
    args: &[f64],
    budget: &mut Budget,
) -> Result<String, String> {
    if args.len() > 4 {
        return Err("fit args[maxiter,tolerance,confidence,covariance_kind]".into());
    }
    let maxiter = integer(args.first().copied().unwrap_or(600.), 5000)?;
    let tol = args.get(1).copied().unwrap_or(1e-9);
    let confidence = args.get(2).copied().unwrap_or(0.95);
    let kind = integer(args.get(3).copied().unwrap_or(0.), 1)?;
    if maxiter == 0 || !(1e-12..=1e-3).contains(&tol) || !(0. < confidence && confidence < 1.) {
        return Err("invalid fit iteration/tolerance/confidence".into());
    }
    let frame = Frame::new(family, data)?;
    let p = frame.p;
    let dim = p + 2 + usize::from(dispersion_family(family));
    let zeros = frame
        .y
        .iter()
        .zip(&frame.weights)
        .filter(|(y, _)| **y == 0.)
        .map(|(_, w)| w)
        .sum::<f64>()
        / frame.total;
    let mean = frame
        .y
        .iter()
        .zip(&frame.weights)
        .map(|(y, w)| y * w / frame.total)
        .sum::<f64>();
    let trials = frame
        .trials
        .iter()
        .zip(&frame.weights)
        .map(|(n, w)| *n as f64 * w / frame.total)
        .sum::<f64>();
    let mut best: Option<Opt> = None;
    let mut starts = 0;
    let disp_modes = if matches!(family, 1 | 3) {
        vec![false, true]
    } else {
        vec![dispersion_family(family)]
    };
    for free_disp in disp_modes {
        for free_pi in [false, true] {
            if (family == 4 && zeros > 0. && !free_pi) || (zeros == 0. && free_pi) {
                continue;
            }
            let ds: Vec<f64> = if free_disp {
                if family == 1 {
                    vec![0.05, 0.3, 0.7]
                } else {
                    vec![0.1, 0.5, 2.]
                }
            } else {
                vec![0.]
            };
            let ps: Vec<f64> = if free_pi {
                vec![
                    (zeros * 0.5).clamp(0.02, 0.8),
                    (zeros * 0.9).clamp(0.1, 0.9),
                ]
            } else {
                vec![0.]
            };
            for disp in &ds {
                for pi in &ps {
                    let mu = if trial_family(family) {
                        (mean / trials / (1. - pi)).clamp(0.01, 0.99)
                    } else {
                        mean / (1. - pi)
                    };
                    let mut t = vec![0.; p + 1];
                    t[0] = if trial_family(family) {
                        (mu / (1. - mu)).ln()
                    } else {
                        mu.ln()
                    };
                    if free_disp {
                        t.push(if family == 1 {
                            (disp / (1. - disp)).ln()
                        } else {
                            disp.ln()
                        });
                    }
                    if free_pi {
                        t.push((pi / (1. - pi)).ln());
                    }
                    starts += 1;
                    let mode = Mode {
                        disp: free_disp,
                        pi: free_pi,
                    };
                    let o = match optimize(&frame, t, mode, maxiter, tol, budget) {
                        Ok(o) => o,
                        Err(e) if budget.used > MAX_WORK => return Err(e),
                        Err(_) => continue,
                    };
                    let replace = best.as_ref().is_none_or(|b| {
                        o.e.value < b.e.value - 1e-11
                            || ((o.e.value - b.e.value).abs() <= 1e-11
                                && ((o.converged && !b.converged)
                                    || (o.converged == b.converged && o.t.len() < b.t.len())))
                    });
                    if replace {
                        best = Some(o);
                    }
                }
            }
        }
    }
    let best = best.ok_or("no finite zero-inflation likelihood candidate")?;
    let (dj, pj) = parameters(&frame, &best.t, best.mode);
    let disp = dj.v;
    let pi = pj.v;
    let guard = best.t.iter().any(|v| v.abs() >= 25.);
    let boundary = !best.mode.pi || (dispersion_family(family) && !best.mode.disp) || guard;
    let e = evaluate(&frame, &best.t, best.mode, budget, true)?;
    let mut normalized: Option<Matrix> = None;
    let mut model_based = None;
    let mut sandwich = None;
    let mut raw_j = vec![vec![0.; dim]; dim];
    for (i, r) in raw_j.iter_mut().enumerate() {
        r[i] = 1.;
    }
    for j in 0..p {
        raw_j[0][j + 1] = -(frame.origins[j] + frame.offsets[j]) / frame.scales[j];
        raw_j[j + 1][j + 1] = 1. / frame.scales[j];
    }
    if best.converged && !boundary {
        budget.charge(dim.pow(4) * 4)?;
        if let Some(inverse) = inverse_spd(&e.h) {
            let ct = inverse
                .iter()
                .map(|r| r.iter().map(|v| v / frame.total).collect())
                .collect::<Matrix>();
            let hc = transform_cov(&e.meat, &ct);
            let mut j = vec![vec![0.; dim]; dim];
            for (i, r) in j.iter_mut().enumerate() {
                r[i] = 1.;
            }
            if dispersion_family(family) {
                j[p + 1][p + 1] = if family == 1 {
                    disp * (1. - disp)
                } else {
                    disp
                };
            }
            j[dim - 1][dim - 1] = pi * (1. - pi);
            let cn = transform_cov(&ct, &j);
            let hn = transform_cov(&hc, &j);
            let selected = if kind == 0 { &cn } else { &hn };
            if selected.iter().flatten().all(|v| v.is_finite()) && inverse_spd(selected).is_some() {
                normalized = Some(selected.clone());
                model_based = Some(transform_cov(&cn, &raw_j));
                sandwich = Some(transform_cov(&hn, &raw_j));
            }
        }
    }
    let covariance = if kind == 0 {
        model_based.as_ref()
    } else {
        sandwich.as_ref()
    };
    let mut coefficients = vec![best.t[0]];
    for j in 0..p {
        let b = finite(best.t[j + 1] / frame.scales[j])?;
        coefficients[0] = finite(coefficients[0] - b * (frame.origins[j] + frame.offsets[j]))?;
        coefficients.push(b);
    }
    let mut saved = vec![
        1.,
        family as f64,
        p as f64,
        confidence,
        kind as f64,
        f64::from(normalized.is_some()),
        f64::from(boundary),
    ];
    saved.extend(&best.t[..p + 1]);
    saved.extend([disp, pi]);
    saved.extend(&frame.origins);
    saved.extend(&frame.offsets);
    saved.extend(&frame.scales);
    if let Some(c) = &normalized {
        saved.extend(c.iter().flatten());
    } else {
        saved.resize(saved.len() + dim * dim, 0.);
    }
    if saved.iter().any(|v| !v.is_finite()) {
        return Err("nonfinite saved zero-inflation model".into());
    }
    let mut physical = coefficients.clone();
    if dispersion_family(family) {
        physical.push(disp);
    }
    physical.push(pi);
    let (se, ci) = if let Some(c) = covariance {
        let se = (0..dim)
            .map(|i| finite(c[i][i].sqrt()))
            .collect::<Result<Vec<_>, _>>()?;
        let z = normal_quantile((1. + confidence) / 2.);
        let ci = physical
            .iter()
            .zip(&se)
            .map(|(v, se)| Ok(vec![finite(v - z * se)?, finite(v + z * se)?]))
            .collect::<Result<Matrix, String>>()?;
        (vector(&se), matrix(&ci))
    } else {
        ("null".into(), "null".into())
    };
    let cov = |c: Option<&Matrix>| c.map_or_else(|| "null".into(), matrix);
    let ll = finite(-e.value * frame.total)?;
    let aic = if best.converged && !boundary {
        finite(2. * dim as f64 - 2. * ll)?.to_string()
    } else {
        "null".into()
    };
    let bic = if best.converged && !boundary {
        finite(frame.total.ln() * dim as f64 - 2. * ll)?.to_string()
    } else {
        "null".into()
    };
    Ok(format!("{{\"version\":1,\"algorithm\":\"zi_fit\",\"family\":{family},\"n_samples\":{},\"effective_sample_size\":{},\"n_features\":{p},\"coefficients\":{},\"dispersion\":{disp},\"pi\":{pi},\"log_likelihood\":{ll},\"aic\":{aic},\"bic\":{bic},\"covariance\":{},\"model_based_covariance\":{},\"sandwich_HC0_covariance\":{},\"standard_errors\":{se},\"confidence_intervals\":{ci},\"covariance_method\":\"{}\",\"confidence\":{confidence},\"converged\":{},\"at_boundary\":{boundary},\"inference_valid\":{},\"finite_mle_existence_verified\":false,\"global_optimum_verified\":false,\"jmp_equivalence_verified\":false,\"gradient_max\":{},\"iterations\":{},\"starts\":{starts},\"line_search_failed\":{},\"predict_parameters\":{},\"work_units\":{}}}",frame.y.len(),frame.total,vector(&coefficients),cov(covariance),cov(model_based.as_ref()),cov(sandwich.as_ref()),if kind==0{"observed_information"}else{"sandwich_HC0"},best.converged,normalized.is_some(),e.g.iter().fold(0_f64,|m,x|m.max(x.abs())),best.iterations,best.line_failed,vector(&saved),budget.used))
}

pub fn run(op: &str, data: &[Vec<f64>], a: &[f64]) -> Result<String, String> {
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 100000
        || width > 1024
        || data.len() > 1_000_000 / width
        || a.len() > 16384
        || a.iter().any(|v| !v.is_finite())
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|x| !x.is_finite() || x.abs() > 1e150))
    {
        return Err("zero-inflation needs finite bounded rectangular data/parameters".into());
    }
    let mut budget = Budget { used: 0 };
    if op == "zi_predict" {
        return predict(data, a, &mut budget);
    }
    let family = match op {
        "zi_binomial" => 0,
        "zi_beta_binomial" => 1,
        "zi_poisson" => 2,
        "zi_negative_binomial" => 3,
        "zi_gamma" => 4,
        _ => return Err("unknown zero-inflation operation".into()),
    };
    fit(family, data, a, &mut budget)
}
