//! Joint unpenalized censored-response likelihoods in JMP 19.2 parameterizations.
//! Fit rows: [features..., lower, upper, kind, integer frequency]. Kinds are
//! 0 exact, 1 left, 2 right, 3 interval; unused endpoint placeholders are finite.
//! All likelihood constants are retained. Numerical optimization is local BFGS;
//! convergence and positive observed information do not prove a global optimum.

// Triangular matrix solves and rank updates require symmetric cross-row indexing.
#![allow(clippy::needless_range_loop)]
use std::cell::Cell;
const LN_2PI: f64 = 1.837_877_066_409_345_3;
const MAX_WORK: usize = 200_000_000;
type Matrix = Vec<Vec<f64>>;
thread_local! {static TOTAL_WORK:Cell<usize>=const {Cell::new(0)};}
fn charge(cost: usize) -> Result<(), String> {
    TOTAL_WORK.with(|counter| {
        let n = counter
            .get()
            .checked_add(cost)
            .ok_or("censored work overflow")?;
        counter.set(n);
        if n > MAX_WORK {
            Err("censored fit/prediction exceeds 200000000 cumulative numerical work units".into())
        } else {
            Ok(())
        }
    })
}
fn used_work() -> usize {
    TOTAL_WORK.with(Cell::get)
}
fn finite(v: f64, name: &str) -> Result<f64, String> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(format!("{name} exceeds finite numerical range"))
    }
}
fn int(v: f64, lo: usize, hi: usize, name: &str) -> Result<usize, String> {
    if v.is_finite() && v.fract() == 0. && v >= lo as f64 && v <= hi as f64 {
        Ok(v as usize)
    } else {
        Err(format!("{name} must be an integer in [{lo},{hi}]"))
    }
}
fn vector(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn matrix(v: &[Vec<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|r| vector(r)).collect::<Vec<_>>().join(",")
    )
}
fn optional(v: Option<f64>) -> String {
    v.map_or_else(|| "null".into(), |x| x.to_string())
}
fn optionals(v: &[Option<f64>]) -> String {
    format!(
        "[{}]",
        v.iter().map(|&x| optional(x)).collect::<Vec<_>>().join(",")
    )
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn weighted_scale(values: impl Iterator<Item = (f64, f64)>, total: f64) -> Result<f64, String> {
    let pairs = values.collect::<Vec<_>>();
    let largest = pairs.iter().map(|(x, _)| x.abs()).fold(0_f64, f64::max);
    if largest == 0. || !largest.is_finite() {
        return Err("constant/unrepresentable variation".into());
    }
    let variance = pairs
        .iter()
        .map(|(x, w)| (x / largest).powi(2) * (w / total))
        .sum::<f64>();
    let value = largest * variance.sqrt();
    if value <= 0. {
        return Err("scale underflows finite range".into());
    }
    finite(value, "weighted scale")
}
fn delta_error(g: &[f64], c: &[Vec<f64>]) -> Option<f64> {
    let gm = g.iter().map(|x| x.abs()).fold(0_f64, f64::max);
    if gm == 0. {
        return Some(0.);
    }
    let cm = c.iter().flatten().map(|x| x.abs()).fold(0_f64, f64::max);
    if cm <= 0. || !cm.is_finite() {
        return None;
    }
    let ng = g.iter().map(|x| x / gm).collect::<Vec<_>>();
    let product = c
        .iter()
        .map(|r| r.iter().zip(&ng).map(|(a, b)| (a / cm) * b).sum::<f64>())
        .collect::<Vec<_>>();
    let value = dot(&ng, &product);
    if value <= 0. || !value.is_finite() {
        return None;
    }
    let se = (gm.ln() + 0.5 * cm.ln() + 0.5 * value.ln()).exp();
    (se > 0. && se.is_finite()).then_some(se)
}
fn log1mexp(a: f64) -> Result<f64, String> {
    if a > 0. {
        return Err("invalid log probability".into());
    }
    Ok(if a < -std::f64::consts::LN_2 {
        (-a.exp()).ln_1p()
    } else {
        (-a.exp_m1()).ln()
    })
}
fn logdiff(a: f64, b: f64) -> Result<f64, String> {
    if b == f64::NEG_INFINITY {
        return Ok(a);
    }
    if a <= b {
        return Err("censor interval probability is numerically unresolved".into());
    }
    Ok(a + log1mexp(b - a)?)
}
fn loggamma(z: f64) -> f64 {
    let c = [
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        return std::f64::consts::PI.ln()
            - (std::f64::consts::PI * z).sin().ln()
            - loggamma(1. - z);
    }
    let z = z - 1.;
    let mut x = 0.999_999_999_999_809_9;
    for (i, c) in c.iter().enumerate() {
        x += c / (z + i as f64 + 1.);
    }
    let t = z + 7.5;
    0.5 * LN_2PI + (z + 0.5) * t.ln() - t + x.ln()
}
// Cancellation-free log Gamma(1+a), 0<a<1/4. The omitted zeta
// series tail is bounded by 2*a^33/(33*(1-a)); no rounded 1+a is formed.
fn loggamma1p_small(a: f64) -> f64 {
    const C: [f64; 31] = [
        0.8224670334241133,
        -0.4006856343865313,
        0.27058080842778454,
        -0.2073855510286741,
        0.16955717699740813,
        -0.14404989676884614,
        0.12550966952474307,
        -0.11133426586956473,
        0.10009945751278182,
        -0.09095401714582904,
        0.08335384054610899,
        -0.0769325164113522,
        0.07143294629536133,
        -0.06666870588242048,
        0.06250095514121302,
        -0.058823978658684564,
        0.055555767627403614,
        -0.052631679379616644,
        0.05000004769810169,
        -0.04761907033014222,
        0.04545455629320467,
        -0.04347826605304026,
        0.04166666915034121,
        -0.04000000119214014,
        0.03846153903467519,
        -0.037037037312989324,
        0.03571428584733336,
        -0.034482758684919304,
        0.03333333336437758,
        -0.03225806453115041,
        0.03125000000727597,
    ];
    let mut sum = 0.;
    for c in C.iter().rev() {
        sum = *c + a * sum;
    }
    a * (-0.577_215_664_901_532_9 + a * sum)
}
fn gamma_logq_fraction(a: f64, x: f64, lead: f64) -> Result<f64, String> {
    let mut b = x + 1. - a;
    let mut c = 1e300;
    let mut d = 1. / b;
    let mut h = d;
    for k in 1..=10000 {
        charge(20)?;
        let an = -(k as f64) * (k as f64 - a);
        b += 2.;
        d = an * d + b;
        if d.abs() < 1e-300 {
            d = 1e-300;
        }
        c = b + an / c;
        if c.abs() < 1e-300 {
            c = 1e-300;
        }
        d = 1. / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.).abs() <= 2e-15 {
            if h <= 0. {
                return Err("gamma continued fraction lost positivity".into());
            }
            return Ok((lead + h.ln()).min(0.));
        }
    }
    Err("incomplete gamma fraction did not converge".into())
}
// Returns both log tails; the smaller probability is calculated directly.
fn gamma_logs(a: f64, x: f64) -> Result<(f64, f64), String> {
    gamma_logs_logx(a, x, x.ln())
}
fn gamma_logs_logx(a: f64, x: f64, log_x: f64) -> Result<(f64, f64), String> {
    if a <= 0. || a > 1e6 || !a.is_finite() {
        return Err("gamma shape outside positive numerical range <=1e6".into());
    }
    if log_x == f64::NEG_INFINITY {
        return Ok((f64::NEG_INFINITY, 0.));
    }
    if x == f64::INFINITY {
        return Ok((0., f64::NEG_INFINITY));
    }
    if a < 0.25 && log_x <= 0. {
        // Integral expansion P=exp(a logx-logGamma(1+a))*(1+a*S),
        // S=sum[n>=1](-x)^n/(n!*(a+n)). Alternating decreasing terms
        // for x<=1 give a next-term remainder bound. Form Q with expm1,
        // avoiding subtraction from a rounded near-one lower probability.
        let prefactor = a * log_x - loggamma1p_small(a);
        let mut term = 1.;
        let mut sum = 0.;
        for n in 1..=10000 {
            charge(12)?;
            term *= -x / n as f64;
            let add = term / (a + n as f64);
            sum += add;
            if add.abs() <= sum.abs() * 1e-16 || term == 0. {
                let lp = prefactor + (a * sum).ln_1p();
                let q = -prefactor.exp_m1() - prefactor.exp() * a * sum;
                if q > 0. && q < 0.5 {
                    let lq = q.ln();
                    return Ok((log1mexp(lq)?, lq));
                }
                return Ok((lp.min(0.), log1mexp(lp.min(0.))?));
            }
        }
        return Err("incomplete gamma alternating series did not converge".into());
    }
    if x == 0. {
        // x underflowed after scaling, while logx and logP remain finite.
        let lp = (a * log_x - loggamma(a + 1.)).min(0.);
        return Ok((lp, log1mexp(lp)?));
    }
    let lead = -x + a * log_x - loggamma(a);
    if x < a + 1. {
        let mut term = 1. / a;
        let mut sum = term;
        for k in 1..=10000 {
            charge(12)?;
            term *= x / (a + k as f64);
            sum += term;
            if term.abs() <= sum.abs() * 2e-15 {
                let lp = (lead + sum.ln()).min(0.);
                if lp > -std::f64::consts::LN_2 {
                    let lq = gamma_logq_fraction(a, x, lead)?;
                    return Ok((log1mexp(lq)?, lq));
                }
                return Ok((lp, log1mexp(lp)?));
            }
        }
        return Err("incomplete gamma series did not converge".into());
    }
    let lq = gamma_logq_fraction(a, x, lead)?;
    Ok((log1mexp(lq)?, lq))
}
fn beta_fraction(a: f64, b: f64, x: f64) -> Result<f64, String> {
    let qab = a + b;
    let qap = a + 1.;
    let qam = a - 1.;
    let mut c = 1.;
    let mut d = 1. - qab * x / qap;
    if d.abs() < 1e-300 {
        d = 1e-300;
    }
    d = 1. / d;
    let mut h = d;
    for m in 1..=10000 {
        charge(40)?;
        let mf = m as f64;
        let m2 = 2. * mf;
        let aa = mf * (b - mf) * x / ((qam + m2) * (a + m2));
        d = 1. + aa * d;
        if d.abs() < 1e-300 {
            d = 1e-300;
        }
        c = 1. + aa / c;
        if c.abs() < 1e-300 {
            c = 1e-300;
        }
        d = 1. / d;
        h *= d * c;
        let aa = -(a + mf) * (qab + mf) * x / ((a + m2) * (qap + m2));
        d = 1. + aa * d;
        if d.abs() < 1e-300 {
            d = 1e-300;
        }
        c = 1. + aa / c;
        if c.abs() < 1e-300 {
            c = 1e-300;
        }
        d = 1. / d;
        let delta = d * c;
        h *= delta;
        if (delta - 1.).abs() < 2e-15 {
            return finite(h, "incomplete beta fraction");
        }
    }
    Err("incomplete beta did not converge".into())
}
fn beta_logs(a: f64, b: f64, x: f64) -> Result<(f64, f64), String> {
    if a <= 0. || b <= 0. || a > 1e6 || b > 1e6 {
        return Err("beta shapes outside positive numerical range <=1e6".into());
    }
    if x <= 0. {
        return Ok((f64::NEG_INFINITY, 0.));
    }
    if x >= 1. {
        return Ok((0., f64::NEG_INFINITY));
    }
    // Exact one-shape forms retain an arbitrarily small complementary shape.
    // Computing 1-SF after a rounded near-one incomplete beta loses such CDFs.
    if a == 1. {
        let lq = b * (-x).ln_1p();
        return Ok((log1mexp(lq)?, lq));
    }
    if b == 1. {
        let lp = a * x.ln();
        return Ok((lp, log1mexp(lp)?));
    }
    let lead = loggamma(a + b) - loggamma(a) - loggamma(b) + a * x.ln() + b * (-x).ln_1p();
    if x < (a + 1.) / (a + b + 2.) {
        let lp = (lead + beta_fraction(a, b, x)?.ln() - a.ln()).min(0.);
        if lp > -std::f64::consts::LN_2 {
            let lq = (lead + beta_fraction(b, a, 1. - x)?.ln() - b.ln()).min(0.);
            return Ok((log1mexp(lq)?, lq));
        }
        Ok((lp, log1mexp(lp)?))
    } else {
        let lq = (lead + beta_fraction(b, a, 1. - x)?.ln() - b.ln()).min(0.);
        if lq > -std::f64::consts::LN_2 {
            let lp = (lead + beta_fraction(a, b, x)?.ln() - a.ln()).min(0.);
            return Ok((lp, log1mexp(lp)?));
        }
        Ok((log1mexp(lq)?, lq))
    }
}
fn normal_logs(z: f64) -> Result<(f64, f64), String> {
    if z == 0. {
        return Ok((-std::f64::consts::LN_2, -std::f64::consts::LN_2));
    }
    if z.abs() > 1e150 {
        return Ok(if z < 0. {
            (f64::NEG_INFINITY, 0.)
        } else {
            (0., f64::NEG_INFINITY)
        });
    }
    let tail = gamma_logs(0.5, 0.5 * z * z)?.1 - std::f64::consts::LN_2;
    if z < 0. {
        Ok((tail, log1mexp(tail)?))
    } else {
        Ok((log1mexp(tail)?, tail))
    }
}
fn logistic(x: f64) -> f64 {
    if x >= 0. {
        1. / (1. + (-x).exp())
    } else {
        let e = x.exp();
        e / (1. + e)
    }
}
fn family_name(f: usize) -> &'static str {
    match f {
        0 => "normal",
        1 => "exponential",
        2 => "gamma",
        3 => "weibull",
        4 => "lognormal",
        _ => "beta",
    }
}
fn log_exponential_cdf(log_x: f64) -> f64 {
    if log_x < -35. {
        log_x
    } else {
        (-(-log_x.exp()).exp_m1()).ln()
    }
}
fn distribution(f: usize, eta: f64, ls: f64, y: f64) -> Result<(f64, f64, f64), String> {
    let sigma = ls.exp();
    if !sigma.is_finite() || sigma == 0. {
        return Err("scale outside numerical range".into());
    }
    if f == 0 {
        let z = (y - eta) / sigma;
        let (lc, lsf) = normal_logs(z)?;
        return Ok((-ls - 0.5 * LN_2PI - 0.5 * z * z, lc, lsf));
    }
    if f == 5 {
        let mu = logistic(eta);
        let a = mu / sigma;
        let b = logistic(-eta) / sigma;
        let (lc, lsf) = beta_logs(a, b, y)?;
        let lp = if y > 0. && y < 1. {
            (a - 1.) * y.ln() + (b - 1.) * (-y).ln_1p() + loggamma(a + b)
                - loggamma(a)
                - loggamma(b)
        } else {
            f64::NEG_INFINITY
        };
        return Ok((lp, lc, lsf));
    }
    if y <= 0. {
        return Ok((f64::NEG_INFINITY, f64::NEG_INFINITY, 0.));
    }
    if f == 1 {
        let z = finite((y.ln() - eta).exp(), "exponential scaled response")?;
        return Ok((-eta - z, log_exponential_cdf(y.ln() - eta), -z));
    }
    if f == 2 {
        let mu = finite(eta.exp(), "gamma mean")?;
        let a = mu / sigma;
        let log_x = y.ln() - ls;
        let (lc, lsf) = gamma_logs_logx(a, log_x.exp(), log_x)?;
        return Ok((
            (a - 1.) * y.ln() - y / sigma - loggamma(a) - a * ls,
            lc,
            lsf,
        ));
    }
    let z = (y.ln() - eta) / sigma;
    if f == 3 {
        let e = finite(z.exp(), "Weibull exponent")?;
        Ok((-y.ln() - ls + z - e, log_exponential_cdf(z), -e))
    } else {
        let (lc, lsf) = normal_logs(z)?;
        Ok((-y.ln() - ls - 0.5 * LN_2PI - 0.5 * z * z, lc, lsf))
    }
}
fn log_average<F: Fn(f64) -> Result<f64, String>>(fun: &F, order: usize) -> Result<f64, String> {
    let mut terms = Vec::with_capacity(order);
    for i in 0..order {
        let mut z = (std::f64::consts::PI * (i as f64 + 0.75) / (order as f64 + 0.5)).cos();
        let mut derivative = 0.;
        for _ in 0..32 {
            charge(order * 6)?;
            let mut p0 = 1.;
            let mut p1 = z;
            for k in 2..=order {
                let pk = ((2 * k - 1) as f64 * z * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = pk;
            }
            derivative = order as f64 * (z * p1 - p0) / (z * z - 1.);
            let next = z - p1 / derivative;
            if (next - z).abs() < 2e-16 {
                z = next;
                break;
            }
            z = next;
        }
        let weight = 1. / ((1. - z * z) * derivative * derivative);
        terms.push(finite(
            weight.ln() + fun(0.5 * (z + 1.))?,
            "interval quadrature term",
        )?);
    }
    let largest = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    finite(
        largest + terms.iter().map(|v| (v - largest).exp()).sum::<f64>().ln(),
        "interval quadrature",
    )
}
fn log_integral<F: Fn(f64) -> Result<f64, String>>(
    log_width: f64,
    log_pdf: f64,
    ratio: F,
) -> Result<f64, String> {
    let mut previous = log_average(&ratio, 16)?;
    for order in [32, 64, 128] {
        let current = log_average(&ratio, order)?;
        if (current - previous).abs() <= 2e-13 * (1. + current.abs()) {
            return finite(
                log_width + log_pdf + current,
                "interval integrated log probability",
            );
        }
        previous = current;
    }
    Err("interval Gauss-Legendre quadrature orders did not agree".into())
}
fn normal_interval(a: f64, log_delta: f64) -> Result<f64, String> {
    let delta = log_delta.exp();
    let base = -0.5 * LN_2PI - 0.5 * a * a;
    log_integral(log_delta, base, |u| {
        let v = delta * u;
        finite(-a * v - 0.5 * v * v, "normal interval density ratio")
    })
}
fn interval(f: usize, eta: f64, ls: f64, lo: f64, hi: f64, log_width: f64) -> Result<f64, String> {
    let sigma = ls.exp();
    if f == 1 {
        if hi <= 0. {
            return Err("interval has no positive exponential mass".into());
        }
        let lower = lo.max(0.);
        let width = if lo < 0. { hi } else { hi - lo };
        let hazard = if lower > 0. {
            (lower.ln() - eta).exp()
        } else {
            0.
        };
        return finite(
            -hazard + log_exponential_cdf(width.ln() - eta),
            "exponential interval log probability",
        );
    }
    if f == 3 {
        if lo <= 0. {
            return Ok(distribution(f, eta, ls, hi)?.1);
        }
        let log_h = (lo.ln() - eta) / sigma;
        let delta = ((hi - lo) / lo).ln_1p() / sigma;
        let log_change = if delta > 50. {
            delta + (-(-delta).exp()).ln_1p()
        } else {
            delta.exp_m1().ln()
        };
        return finite(
            -log_h.exp() + log_exponential_cdf(log_h + log_change),
            "Weibull interval log probability",
        );
    }
    let (_, a, b) = distribution(f, eta, ls, lo)?;
    let (_, c, d) = distribution(f, eta, ls, hi)?;
    let raw = if a < -std::f64::consts::LN_2 {
        logdiff(c, a)
    } else {
        logdiff(b, d)
    };
    if f == 0 {
        let aa = (lo - eta) / sigma;
        let ld = log_width - ls;
        if raw.is_ok() && ld.exp() * (1. + aa.abs()) > 1e-3 {
            return raw;
        }
        return normal_interval(aa, ld);
    }
    if f == 4 {
        if lo <= 0. {
            return Ok(c);
        }
        let aa = (lo.ln() - eta) / sigma;
        let delta = ((hi - lo) / lo).ln_1p() / sigma;
        if raw.is_ok() && delta * (1. + aa.abs()) > 1e-3 {
            return raw;
        }
        return normal_interval(aa, delta.ln());
    }
    if f == 2 {
        if lo <= 0. {
            return Ok(c);
        }
        let width = hi - lo;
        if raw.is_ok() && width / lo > 1e-4 {
            return raw;
        }
        let shape = eta.exp() / sigma;
        let base = (shape - 1.) * lo.ln() - lo / sigma - loggamma(shape) - shape * ls;
        return log_integral(width.ln(), base, |u| {
            finite(
                (shape - 1.) * ((width / lo) * u).ln_1p() - (width / sigma) * u,
                "gamma interval density ratio",
            )
        });
    }
    if lo <= 0. {
        return Ok(c);
    }
    if hi >= 1. {
        return Ok(b);
    }
    let width = hi - lo;
    if raw.is_ok() && width / lo.min(1. - hi) > 1e-4 {
        return raw;
    }
    let shape_a = logistic(eta) / sigma;
    let shape_b = logistic(-eta) / sigma;
    let base =
        (shape_a - 1.) * lo.ln() + (shape_b - 1.) * (-lo).ln_1p() + loggamma(shape_a + shape_b)
            - loggamma(shape_a)
            - loggamma(shape_b);
    log_integral(width.ln(), base, |u| {
        finite(
            (shape_a - 1.) * ((width / lo) * u).ln_1p()
                + (shape_b - 1.) * (-(width / (1. - lo)) * u).ln_1p(),
            "beta interval density ratio",
        )
    })
}
fn likelihood(
    f: usize,
    eta: f64,
    ls: f64,
    lo: f64,
    hi: f64,
    kind: usize,
    log_width: f64,
) -> Result<f64, String> {
    let v = match kind {
        0 => distribution(f, eta, ls, lo)?.0,
        1 => distribution(f, eta, ls, hi)?.1,
        2 => distribution(f, eta, ls, lo)?.2,
        _ => interval(f, eta, ls, lo, hi, log_width)?,
    };
    finite(v, "censored log likelihood")
}

fn inverse_spd(a: &[Vec<f64>]) -> Option<Matrix> {
    let n = a.len();
    let mut l = vec![vec![0.; n]; n];
    let max = a
        .iter()
        .enumerate()
        .map(|(i, r)| r[i].abs())
        .fold(0_f64, f64::max);
    if max == 0. || !max.is_finite() {
        return None;
    }
    for i in 0..n {
        for j in 0..=i {
            let mut v = a[i][j];
            for k in 0..j {
                v -= l[i][k] * l[j][k];
            }
            if i == j {
                if !v.is_finite() || v <= max * 1e-12 {
                    return None;
                }
                l[i][j] = v.sqrt();
            } else {
                l[i][j] = v / l[j][j];
            }
        }
    }
    let mut out = vec![vec![0.; n]; n];
    for col in 0..n {
        let mut y = vec![0.; n];
        for i in 0..n {
            let s = (0..i).map(|k| l[i][k] * y[k]).sum::<f64>();
            y[i] = (f64::from(i == col) - s) / l[i][i];
        }
        let mut x = vec![0.; n];
        for i in (0..n).rev() {
            x[i] = (y[i] - (i + 1..n).map(|k| l[k][i] * x[k]).sum::<f64>()) / l[i][i];
        }
        for i in 0..n {
            if !x[i].is_finite() {
                return None;
            }
            out[i][col] = x[i];
        }
    }
    Some(out)
}
struct Problem {
    family: usize,
    p: usize,
    d: usize,
    x: Matrix,
    lo: Vec<f64>,
    hi: Vec<f64>,
    kind: Vec<usize>,
    w: Vec<f64>,
    total: f64,
    origins: Vec<f64>,
    offsets: Vec<f64>,
    scales: Vec<f64>,
    yorigin: f64,
    yscale: f64,
    counts: [f64; 4],
    log_width: Vec<f64>,
}
impl Problem {
    fn new(family: usize, data: &[Vec<f64>]) -> Result<Self, String> {
        let width = data[0].len();
        if !(4..=20).contains(&width) {
            return Err(
                "censored fit requires 0..16 features then lower,upper,kind,frequency".into(),
            );
        }
        let p = width - 4;
        let d = p + 1 + usize::from(family != 1);
        let mut active = Vec::new();
        let mut total = 0.;
        let mut counts = [0.; 4];
        for r in data {
            let kind = int(r[p + 2], 0, 3, "censor kind")?;
            int(r[p + 3], 0, 1_000_000, "frequency")?;
            if r[p + 3] == 0. {
                continue;
            }
            let lo = r[p];
            let hi = r[p + 1];
            if kind == 0 && lo != hi || kind == 3 && lo >= hi {
                return Err("exact endpoints must match; interval endpoints must increase".into());
            }
            if family > 0
                && family < 5
                && (kind == 0 && lo <= 0. || kind == 1 && hi <= 0. || kind == 3 && hi <= 0.)
            {
                return Err(
                    "positive distribution requires observations/finite upper limits above zero"
                        .into(),
                );
            }
            if family == 5
                && (kind == 0 && !(lo > 0. && lo < 1.)
                    || kind == 1 && hi <= 0.
                    || kind == 2 && lo >= 1.
                    || kind == 3 && (hi <= 0. || lo >= 1.))
            {
                return Err(
                    "beta response/censoring interval has no positive probability on (0,1)".into(),
                );
            }
            active.push(r);
            total += r[p + 3];
            counts[kind] += r[p + 3];
        }
        if active.len() < p + 1 || total <= d as f64 || total > 1e9 {
            return Err("censored fit needs full-rank data, at least p+1 positive-frequency rows and frequency sum in (parameters,1e9]".into());
        }
        if counts[0] + counts[1] + counts[3] == 0. || counts[0] + counts[2] + counts[3] == 0. {
            return Err(
                "all observations censored in one direction have no finite unpenalized maximum"
                    .into(),
            );
        }
        // With no exact failures and an overlapping censor interval, a free
        // dispersion can collapse to zero while every observation has probability
        // tending to one. Positive-definite local information cannot prove a finite
        // maximum exists, so reject this directly identifiable recession case.
        if family != 1 && counts[0] == 0. {
            let mut lower = if family == 0 { f64::NEG_INFINITY } else { 0. };
            let mut upper = if family == 5 { 1. } else { f64::INFINITY };
            for r in &active {
                let kind = r[p + 2] as usize;
                if kind == 2 || kind == 3 {
                    lower = lower.max(r[p]);
                }
                if kind == 1 || kind == 3 {
                    upper = upper.min(r[p + 1]);
                }
            }
            if lower < upper {
                return Err("common censor interval permits zero-dispersion supremum without finite identified maximum".into());
            }
        }
        let n = active.len();
        let mut x = vec![vec![1.; p + 1]; n];
        let mut origins = vec![];
        let mut offsets = vec![];
        let mut scales = vec![];
        for j in 0..p {
            let origin = active[0][j];
            let offset = active
                .iter()
                .map(|r| (r[j] - origin) * r[p + 3] / total)
                .sum::<f64>();
            let scale = weighted_scale(
                active.iter().map(|r| ((r[j] - origin) - offset, r[p + 3])),
                total,
            )?;
            if !scale.is_finite() || scale <= 0. {
                return Err(
                    "constant/unrepresentable covariate is not identifiable with intercept".into(),
                );
            }
            for i in 0..n {
                x[i][j + 1] = ((active[i][j] - origin) - offset) / scale;
            }
            origins.push(origin);
            offsets.push(offset);
            scales.push(scale);
        }
        // Weighted QR rank test, without forming the squared-condition Gram matrix.
        let mut orthogonal: Matrix = vec![];
        for j in 0..=p {
            let mut col = (0..n)
                .map(|i| x[i][j] * (active[i][p + 3] / total).sqrt())
                .collect::<Vec<_>>();
            for _ in 0..2 {
                for q in &orthogonal {
                    let v = dot(&col, q);
                    for i in 0..n {
                        col[i] -= v * q[i];
                    }
                }
            }
            let norm = dot(&col, &col).sqrt();
            if norm <= 1e-10 {
                return Err("censored fixed design is rank deficient".into());
            }
            for v in &mut col {
                *v /= norm;
            }
            orthogonal.push(col);
        }
        let proxy = active
            .iter()
            .map(|r| match r[p + 2] as usize {
                0 | 2 => r[p],
                1 => r[p + 1],
                _ => r[p] + (r[p + 1] - r[p]) * 0.5,
            })
            .collect::<Vec<_>>();
        let (yorigin, yscale) = if family == 0 {
            let origin = proxy[0];
            let mean = proxy
                .iter()
                .zip(&active)
                .map(|(v, r)| (v - origin) * r[p + 3] / total)
                .sum::<f64>();
            let scale = weighted_scale(
                proxy
                    .iter()
                    .zip(&active)
                    .map(|(v, r)| ((v - origin) - mean, r[p + 3])),
                total,
            )?;
            if !scale.is_finite() || scale <= 0. {
                return Err("normal response variation is zero or unrepresentable".into());
            }
            (origin + mean, scale)
        } else {
            (0., 1.)
        };
        Ok(Self {
            family,
            p,
            d,
            x,
            lo: active.iter().map(|r| (r[p] - yorigin) / yscale).collect(),
            hi: active
                .iter()
                .map(|r| (r[p + 1] - yorigin) / yscale)
                .collect(),
            kind: active.iter().map(|r| r[p + 2] as usize).collect(),
            w: active.iter().map(|r| r[p + 3]).collect(),
            total,
            origins,
            offsets,
            scales,
            yorigin,
            yscale,
            counts,
            log_width: active
                .iter()
                .map(|r| {
                    if r[p + 2] == 3. {
                        (r[p + 1] - r[p]).ln() - yscale.ln()
                    } else {
                        0.
                    }
                })
                .collect(),
        })
    }
    fn value(&self, theta: &[f64]) -> Result<f64, String> {
        charge(self.x.len() * self.d)?;
        if theta.iter().any(|x| !x.is_finite() || x.abs() > 1e6)
            || self.family != 1 && theta[self.d - 1].abs() > 30.
        {
            return Err("censored numerical coefficient/scale guard reached".into());
        }
        let ls = if self.family == 1 {
            0.
        } else {
            theta[self.d - 1]
        };
        let mut sum = 0.;
        for i in 0..self.x.len() {
            let eta = dot(&self.x[i], &theta[..=self.p]);
            let l = likelihood(
                self.family,
                eta,
                ls,
                self.lo[i],
                self.hi[i],
                self.kind[i],
                self.log_width[i],
            )?;
            let l = if self.kind[i] == 0 && self.family == 0 {
                l - self.yscale.ln()
            } else {
                l
            };
            sum -= self.w[i] / self.total * l;
        }
        finite(sum, "average negative likelihood")
    }
    fn initial(&self) -> Result<Vec<f64>, String> {
        let mut y = vec![];
        for i in 0..self.lo.len() {
            let v = match self.kind[i] {
                0 | 2 => self.lo[i],
                1 => self.hi[i],
                _ => self.lo[i] + 0.5 * (self.hi[i] - self.lo[i]),
            };
            y.push(if self.family == 3 || self.family == 4 {
                v.max(1e-12).ln()
            } else {
                v
            });
        }
        let mean = y
            .iter()
            .zip(&self.w)
            .map(|(y, w)| y * w / self.total)
            .sum::<f64>();
        let variance = y
            .iter()
            .zip(&self.w)
            .map(|(y, w)| (y - mean).powi(2) * w / self.total)
            .sum::<f64>();
        let mut t = vec![0.; self.d];
        t[0] = match self.family {
            1 | 2 => mean.max(1e-8).ln(),
            5 => {
                let m = mean.clamp(1e-6, 1. - 1e-6);
                (m / (1. - m)).ln()
            }
            _ => mean,
        };
        if self.family != 1 {
            let sigma = match self.family {
                0 | 4 => variance.sqrt(),
                2 => variance / mean.max(1e-8),
                3 => variance.sqrt() * 6_f64.sqrt() / std::f64::consts::PI,
                5 => variance / (mean * (1. - mean) - variance).max(variance * 0.05),
                _ => 1.,
            };
            t[self.d - 1] = sigma.clamp(1e-6, 1e6).ln();
        }
        self.value(&t)?;
        Ok(t)
    }
}
fn gradient<F: Fn(&[f64]) -> Result<f64, String>>(f: &F, x: &[f64]) -> Result<Vec<f64>, String> {
    let mut out = vec![0.; x.len()];
    for j in 0..x.len() {
        let h = 6e-6 * (1. + x[j].abs());
        let mut a = x.to_vec();
        let mut b = x.to_vec();
        a[j] += h;
        b[j] -= h;
        out[j] = finite((f(&a)? - f(&b)?) / (2. * h), "likelihood gradient")?;
    }
    Ok(out)
}
struct Opt {
    theta: Vec<f64>,
    value: f64,
    gradient: f64,
    iterations: usize,
    converged: bool,
    line_failed: bool,
}
fn optimize(p: &Problem, initial: Vec<f64>, maxiter: usize, tol: f64) -> Result<Opt, String> {
    let n = initial.len();
    let mut x = initial;
    let mut v = p.value(&x)?;
    let mut g = gradient(&|x| p.value(x), &x)?;
    let mut h = vec![vec![0.; n]; n];
    for i in 0..n {
        h[i][i] = 1.;
    }
    let mut iterations = 0;
    let mut failed = false;
    for k in 0..maxiter {
        charge(n * n * 8)?;
        iterations = k;
        if g.iter().map(|v| v.abs()).fold(0_f64, f64::max) <= tol {
            break;
        }
        let mut direction = h.iter().map(|r| -dot(r, &g)).collect::<Vec<_>>();
        if dot(&direction, &g) >= -1e-15 * dot(&g, &g) {
            direction = g.iter().map(|v| -v).collect();
            for i in 0..n {
                for j in 0..n {
                    h[i][j] = f64::from(i == j);
                }
            }
        }
        let largest = direction.iter().map(|v| v.abs()).fold(0_f64, f64::max);
        if largest > 10. {
            for d in &mut direction {
                *d *= 10. / largest;
            }
        }
        let slope = dot(&g, &direction);
        let mut alpha = 1.;
        let mut accepted = None;
        for _ in 0..50 {
            let candidate = x
                .iter()
                .zip(&direction)
                .map(|(x, d)| x + alpha * d)
                .collect::<Vec<_>>();
            if let Ok(nv) = p.value(&candidate) {
                if nv <= v + 1e-4 * alpha * slope {
                    accepted = Some((candidate, nv));
                    break;
                }
            }
            alpha *= 0.5;
        }
        let Some((nx, nv)) = accepted else {
            failed = true;
            break;
        };
        let ng = gradient(&|x| p.value(x), &nx)?;
        let s = nx.iter().zip(&x).map(|(a, b)| a - b).collect::<Vec<_>>();
        let y = ng.iter().zip(&g).map(|(a, b)| a - b).collect::<Vec<_>>();
        let sy = dot(&s, &y);
        if sy > 1e-12 * dot(&s, &s).sqrt() * dot(&y, &y).sqrt() {
            let hy = h.iter().map(|r| dot(r, &y)).collect::<Vec<_>>();
            let yhy = dot(&y, &hy);
            for i in 0..n {
                for j in 0..n {
                    h[i][j] +=
                        (1. + yhy / sy) * s[i] * s[j] / sy - (s[i] * hy[j] + hy[i] * s[j]) / sy;
                }
            }
        }
        x = nx;
        v = nv;
        g = ng;
        iterations = k + 1;
    }
    let gn = g.iter().map(|v| v.abs()).fold(0_f64, f64::max);
    let guard = x.iter().take(p.p + 1).any(|v| v.abs() >= 999999.)
        || p.family != 1 && x[p.d - 1].abs() >= 29.999;
    Ok(Opt {
        theta: x,
        value: v,
        gradient: gn,
        iterations,
        converged: gn <= tol && !guard,
        line_failed: failed,
    })
}
fn hessian(p: &Problem, x: &[f64]) -> Result<Matrix, String> {
    let d = x.len();
    let f0 = p.value(x)?;
    let mut h = vec![vec![0.; d]; d];
    for i in 0..d {
        for j in 0..=i {
            let mut estimates = [0.; 2];
            for (level, e) in estimates.iter_mut().enumerate() {
                let hi = 2e-4 * (1. + x[i].abs()) / (level + 1) as f64;
                let hj = 2e-4 * (1. + x[j].abs()) / (level + 1) as f64;
                *e = if i == j {
                    let mut a = x.to_vec();
                    let mut b = x.to_vec();
                    a[i] += hi;
                    b[i] -= hi;
                    (p.value(&a)? - 2. * f0 + p.value(&b)?) / (hi * hi)
                } else {
                    let mut a = x.to_vec();
                    let mut b = x.to_vec();
                    let mut c = x.to_vec();
                    let mut z = x.to_vec();
                    a[i] += hi;
                    a[j] += hj;
                    b[i] += hi;
                    b[j] -= hj;
                    c[i] -= hi;
                    c[j] += hj;
                    z[i] -= hi;
                    z[j] -= hj;
                    (p.value(&a)? - p.value(&b)? - p.value(&c)? + p.value(&z)?) / (4. * hi * hj)
                };
            }
            let v = finite(
                (4. * estimates[1] - estimates[0]) / 3. * p.total,
                "observed information",
            )?;
            h[i][j] = v;
            h[j][i] = v;
        }
    }
    Ok(h)
}
fn normal_quantile(p: f64) -> Result<f64, String> {
    let mut lo = -40.;
    let mut hi = 40.;
    for _ in 0..100 {
        let mid = 0.5 * (lo + hi);
        if normal_logs(mid)?.0 < p.ln() {
            lo = mid
        } else {
            hi = mid
        }
    }
    Ok(0.5 * (lo + hi))
}
fn fitted(data: &[Vec<f64>], args: &[f64], family: usize) -> Result<String, String> {
    let p = Problem::new(family, data)?;
    if args.len() > 3 && args.len() != 3 + p.d {
        return Err("fit accepts [maxiter,tolerance,confidence,optional original coefficients and log_sigma]".into());
    }
    let maxiter = int(
        args.first().copied().unwrap_or(400.),
        0,
        2000,
        "max iterations",
    )?;
    let tol = args.get(1).copied().unwrap_or(1e-7);
    let confidence = args.get(2).copied().unwrap_or(0.95);
    if !(1e-10..=1e-3).contains(&tol) || !(0.5..0.999999).contains(&confidence) {
        return Err("tolerance must be 1e-10..1e-3 and confidence [0.5,0.999999)".into());
    }
    let mut initial = p.initial()?;
    if args.len() > 3 {
        let raw = &args[3..];
        initial[0] = (raw[0] - p.yorigin) / p.yscale;
        for j in 0..p.p {
            initial[j + 1] = raw[j + 1] * p.scales[j] / p.yscale;
            initial[0] += raw[j + 1] * (p.origins[j] + p.offsets[j]) / p.yscale;
        }
        if family != 1 {
            initial[p.d - 1] = raw[p.d - 1] - p.yscale.ln();
        }
    }
    let mut best = optimize(&p, initial.clone(), maxiter, tol)?;
    if family != 1 && args.len() <= 3 && maxiter > 0 {
        for offset in [-0.7, 0.7] {
            let mut start = initial.clone();
            start[p.d - 1] += offset;
            if let Ok(candidate) = optimize(&p, start, maxiter, tol) {
                if candidate.value < best.value {
                    best = candidate;
                }
            }
        }
    }
    let h = hessian(&p, &best.theta)?;
    charge(p.d * p.d * p.d * 4)?;
    let cov = if best.converged {
        inverse_spd(&h)
    } else {
        None
    };
    let mut transform = vec![vec![0.; p.d]; p.d];
    transform[0][0] = p.yscale;
    for j in 0..p.p {
        transform[0][j + 1] = -(p.origins[j] + p.offsets[j]) * p.yscale / p.scales[j];
        transform[j + 1][j + 1] = p.yscale / p.scales[j];
    }
    if family != 1 {
        transform[p.d - 1][p.d - 1] = 1.;
    }
    charge(p.d * p.d * p.d * p.d * 4)?;
    let raw_cov = cov.as_ref().map(|c| {
        let mut out = vec![vec![0.; p.d]; p.d];
        for i in 0..p.d {
            for j in 0..p.d {
                for k in 0..p.d {
                    for l in 0..p.d {
                        out[i][j] += transform[i][k] * c[k][l] * transform[j][l];
                    }
                }
            }
        }
        out
    });
    let raw_cov = raw_cov.filter(|c| {
        c.iter().flatten().all(|v| v.is_finite()) && c.iter().enumerate().all(|(i, r)| r[i] > 0.)
    });
    let raw = transform[..=p.p]
        .iter()
        .map(|r| dot(r, &best.theta))
        .collect::<Vec<_>>();
    let mut coefficients = raw;
    coefficients[0] += p.yorigin;
    for v in &coefficients {
        finite(*v, "original coefficient")?;
    }
    let se = raw_cov
        .as_ref()
        .map(|c| (0..p.d).map(|i| c[i][i].max(0.).sqrt()).collect::<Vec<_>>());
    let z = normal_quantile(0.5 + 0.5 * confidence)?;
    let mut estimates = coefficients.clone();
    if family != 1 {
        estimates.push(best.theta[p.d - 1] + p.yscale.ln())
    }
    let intervals = se.as_ref().map(|s| {
        estimates
            .iter()
            .zip(s)
            .map(|(v, s)| vec![v - z * s, v + z * s])
            .collect::<Vec<_>>()
    });
    let mut saved = vec![
        1.,
        family as f64,
        p.p as f64,
        confidence,
        p.yorigin,
        p.yscale,
        f64::from(cov.is_some()),
    ];
    saved.extend_from_slice(&best.theta[..=p.p]);
    saved.push(if family == 1 { 0. } else { best.theta[p.d - 1] });
    saved.extend_from_slice(&p.origins);
    saved.extend_from_slice(&p.offsets);
    saved.extend_from_slice(&p.scales);
    saved.extend(
        cov.clone()
            .unwrap_or_else(|| vec![vec![0.; p.d]; p.d])
            .iter()
            .flatten(),
    );
    if saved.len() > 16384 || saved.iter().any(|x| !x.is_finite()) {
        return Err("saved censored model is not representable at public ABI".into());
    }
    let ll = -best.value * p.total;
    let sigma = if family == 1 {
        None
    } else {
        Some(best.theta[p.d - 1].exp() * p.yscale)
    };
    let covariance = raw_cov
        .as_ref()
        .map_or_else(|| "null".into(), |c| matrix(c));
    let standard_errors = se.as_ref().map_or_else(|| "null".into(), |v| vector(v));
    let intervals = intervals
        .as_ref()
        .map_or_else(|| "null".into(), |v| matrix(v));
    Ok(format!("{{\"version\":1,\"algorithm\":\"censored_{}\",\"n_samples\":{},\"n_features\":{},\"effective_sample_size\":{},\"coefficients\":{},\"sigma\":{},\"log_sigma\":{},\"log_likelihood\":{},\"aic\":{},\"bic\":{},\"covariance\":{},\"standard_errors\":{},\"confidence_intervals\":{},\"confidence\":{},\"censor_counts\":{},\"converged\":{},\"information_identified\":{},\"inference_valid\":{},\"global_optimum_verified\":false,\"iterations\":{},\"gradient_max\":{},\"line_search_failed\":{},\"work_units\":{},\"predict_parameters\":{}}}",family_name(family),p.x.len(),p.p,p.total,vector(&coefficients),optional(sigma),optional(sigma.map(f64::ln)),ll,2.*p.d as f64-2.*ll,(p.d as f64)*p.total.ln()-2.*ll,covariance,standard_errors,intervals,confidence,vector(&p.counts),best.converged,cov.is_some(),best.converged&&raw_cov.is_some(),best.iterations,best.gradient,best.line_failed,used_work(),vector(&saved)))
}
fn mean(f: usize, eta: f64, ls: f64) -> Result<f64, String> {
    finite(
        match f {
            0 => eta,
            1 | 2 => eta.exp(),
            3 => (eta + loggamma(1. + ls.exp())).exp(),
            4 => (eta + 0.5 * (2. * ls).exp()).exp(),
            _ => logistic(eta),
        },
        "conditional mean",
    )
}
fn positive_quantile(v: f64, name: &str) -> Result<f64, String> {
    if v <= 0. {
        return Err(format!("{name} underflows positive distribution support"));
    }
    finite(v, name)
}
fn quantile(f: usize, eta: f64, ls: f64, p: f64) -> Result<f64, String> {
    if f == 0 {
        return finite(eta + ls.exp() * normal_quantile(p)?, "normal quantile");
    }
    if f == 1 {
        return positive_quantile((eta + (-(-p).ln_1p()).ln()).exp(), "exponential quantile");
    }
    if f == 3 {
        return positive_quantile(
            (eta + ls.exp() * (-(-p).ln_1p()).ln()).exp(),
            "Weibull quantile",
        );
    }
    if f == 4 {
        return positive_quantile(
            (eta + ls.exp() * normal_quantile(p)?).exp(),
            "lognormal quantile",
        );
    }
    let mut lower = -745.;
    let mut upper = if f == 5 { 0. } else { eta.max(ls) + 1. };
    if f == 5 {
        for _ in 0..8 {
            if distribution(f, eta, ls, logistic(upper))?.1 >= p.ln() {
                break;
            }
            upper = if upper == 0. { 1. } else { upper * 2. };
            if upper > 36.7 {
                return Err("beta quantile rounds outside representable open support".into());
            }
        }
    }
    if f == 2 {
        for _ in 0..1024 {
            if distribution(f, eta, ls, upper.exp())?.1 >= p.ln() {
                break;
            }
            upper += std::f64::consts::LN_2;
            if upper > 709. {
                return Err("gamma quantile cannot be bracketed in finite response range".into());
            }
        }
    }
    for _ in 0..140 {
        let mid = lower + 0.5 * (upper - lower);
        let y = if f == 5 { logistic(mid) } else { mid.exp() };
        if distribution(f, eta, ls, y)?.1 < p.ln() {
            lower = mid
        } else {
            upper = mid
        }
    }
    let q = if f == 5 {
        logistic(0.5 * (lower + upper))
    } else {
        (0.5 * (lower + upper)).exp()
    };
    if q <= 0. || f == 5 && q >= 1. {
        return Err("quantile is not representable inside distribution support".into());
    }
    finite(q, "distribution quantile")
}
fn prediction(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if args.len() < 10 || args[0] != 1. {
        return Err("prediction needs complete version-1 censored predict_parameters".into());
    }
    let f = int(args[1], 0, 5, "family")?;
    let p = int(args[2], 0, 16, "feature count")?;
    let d = p + 1 + usize::from(f != 1);
    if args.len() != 9 + 4 * p + d * d || data[0].len() != p + 2 {
        return Err("saved censored model length or prediction row width is invalid".into());
    }
    let confidence = args[3];
    let yo = args[4];
    let ys = args[5];
    let has_cov = int(args[6], 0, 1, "covariance flag")? == 1;
    if !(0.5..0.999999).contains(&confidence) || ys <= 0. || f != 0 && (yo != 0. || ys != 1.) {
        return Err("invalid saved response transform/confidence".into());
    }
    let coef = &args[7..8 + p];
    let ls = args[8 + p];
    if f == 1 && ls != 0. || ls.abs() > 30. {
        return Err("invalid saved scale".into());
    }
    let origins = &args[9 + p..9 + 2 * p];
    let offsets = &args[9 + 2 * p..9 + 3 * p];
    let scales = &args[9 + 3 * p..9 + 4 * p];
    if scales.iter().any(|v| *v <= 0.) {
        return Err("invalid saved feature scales".into());
    }
    let c = args[9 + 4 * p..]
        .chunks_exact(d)
        .map(<[f64]>::to_vec)
        .collect::<Vec<_>>();
    if has_cov {
        for i in 0..d {
            for j in 0..d {
                let scale = 1_f64.max(c[i][j].abs()).max(c[j][i].abs());
                if (c[i][j] / scale - c[j][i] / scale).abs()
                    > 1e-8 * (1. / scale + (c[i][j] / scale).abs() + (c[j][i] / scale).abs())
                {
                    return Err("saved covariance is not symmetric".into());
                }
            }
        }
        if inverse_spd(&c).is_none() {
            return Err("saved covariance is not positive definite".into());
        }
    } else if c.iter().flatten().any(|v| *v != 0.) {
        return Err("absent covariance must be encoded as zeros".into());
    }
    let z = normal_quantile(0.5 + 0.5 * confidence)?;
    let mut means = vec![];
    let mut cdf = vec![];
    let mut survival = vec![];
    let mut density = vec![];
    let mut quantiles = vec![];
    let mut mean_se = vec![];
    let mut mean_intervals = vec![];
    let mut quantile_se = vec![];
    for r in data {
        charge(d * d * 8 + 32)?;
        let probability = r[p + 1];
        if probability <= 0. || probability >= 1. {
            return Err("prediction quantile probabilities must lie strictly in (0,1)".into());
        }
        let mut x = vec![1.];
        for j in 0..p {
            x.push(finite(
                ((r[j] - origins[j]) - offsets[j]) / scales[j],
                "prediction feature",
            )?)
        }
        let eta = finite(dot(&x, coef), "prediction location")?;
        let y = (r[p] - yo) / ys;
        let (lp, lc, lsf) = distribution(f, eta, ls, y)?;
        let m = finite(yo + ys * mean(f, eta, ls)?, "predicted mean")?;
        let q = finite(
            yo + ys * quantile(f, eta, ls, probability)?,
            "predicted quantile",
        )?;
        means.push(m);
        cdf.push(lc.exp());
        survival.push(lsf.exp());
        density.push(finite(lp.exp() / ys, "prediction density")?);
        quantiles.push(q);
        let mut theta = coef.to_vec();
        if f != 1 {
            theta.push(ls)
        }
        let get = |t: &[f64], is_quantile: bool| -> Result<f64, String> {
            let e = dot(&x, &t[..=p]);
            let l = if f == 1 { 0. } else { t[d - 1] };
            if is_quantile {
                quantile(f, e, l, probability).map(|v| v * ys)
            } else {
                mean(f, e, l).map(|v| v * ys)
            }
        };
        let sem = if has_cov {
            let g = if f == 5 {
                let derivative = logistic(eta) * logistic(-eta) * ys;
                let mut g = x.iter().map(|v| v * derivative).collect::<Vec<_>>();
                g.push(0.);
                g
            } else {
                gradient(&|t| get(t, false), &theta)?
            };
            delta_error(&g, &c)
        } else {
            None
        };
        let seq = if has_cov {
            let g = gradient(&|t| get(t, true), &theta)?;
            delta_error(&g, &c)
        } else {
            None
        };
        mean_se.push(sem);
        quantile_se.push(seq);
        mean_intervals.push(sem.and_then(|s| {
            let lower = m - z * s;
            let upper = m + z * s;
            (lower.is_finite() && upper.is_finite()).then_some(vec![lower, upper])
        }));
    }
    let ci = format!(
        "[{}]",
        mean_intervals
            .iter()
            .map(|v| v.as_ref().map_or_else(|| "null".into(), |v| vector(v)))
            .collect::<Vec<_>>()
            .join(",")
    );
    Ok(format!("{{\"version\":1,\"algorithm\":\"censored_predict\",\"family\":\"{}\",\"means\":{},\"cdf\":{},\"survival\":{},\"density\":{},\"quantiles\":{},\"mean_standard_errors\":{},\"mean_confidence_intervals\":{},\"quantile_standard_errors\":{},\"confidence\":{},\"uncertainty\":\"observed_information_delta_normal_approximation\"}}",family_name(f),vector(&means),vector(&cdf),vector(&survival),vector(&density),vector(&quantiles),optionals(&mean_se),ci,optionals(&quantile_se),confidence))
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    TOTAL_WORK.with(|counter| counter.set(0));
    run_inner(operation, data, args)
}
fn run_inner(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    charge(
        data.len()
            .saturating_mul(data.first().map_or(0, Vec::len))
            .saturating_mul(8),
    )?;
    let width = data.first().map_or(0, Vec::len);
    if width == 0
        || data.len() > 100000
        || width > 1024
        || data.len() > 1_000_000 / width
        || data
            .iter()
            .any(|r| r.len() != width || r.iter().any(|v| !v.is_finite() || v.abs() > 1e150))
        || args.len() > 16384
        || args.iter().any(|v| !v.is_finite())
    {
        return Err(
            "censored regression requires bounded finite rectangular data/parameters".into(),
        );
    }
    match operation {
        "censored_normal" => fitted(data, args, 0),
        "censored_exponential" => fitted(data, args, 1),
        "censored_gamma" => fitted(data, args, 2),
        "censored_weibull" => fitted(data, args, 3),
        "censored_lognormal" => fitted(data, args, 4),
        "censored_beta" => fitted(data, args, 5),
        "censored_predict" => prediction(data, args),
        _ => Err("unknown censored regression operation".into()),
    }
}
