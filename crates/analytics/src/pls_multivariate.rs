//! Centered multivariate PLS (two-block regression deflation or SIMPLS),
//! indicator-response PLS-DA, VIP and portable conditional-score prediction.
//! Contract and independent references: tasks/jmp-pls-api.json.
use crate::{explore::symmetric_eigen, glm::student_quantile};
type Matrix = Vec<Vec<f64>>;
const MAX_WORK: usize = 200_000_000;
struct Budget(usize);
impl Budget {
    fn take(&mut self, work: usize) -> Result<(), String> {
        self.0 = self.0.checked_add(work).ok_or("PLS work overflow")?;
        if self.0 > MAX_WORK {
            return Err("PLS cumulative work budget exhausted".into());
        }
        Ok(())
    }
}
fn integer(v: f64, lo: usize, hi: usize) -> Result<usize, String> {
    if !v.is_finite() || v.fract() != 0. || v < lo as f64 || v > hi as f64 {
        return Err(format!("PLS integer must be in [{lo},{hi}]"));
    }
    Ok(v as usize)
}
fn flag(v: f64) -> Result<bool, String> {
    Ok(integer(v, 0, 1)? == 1)
}
fn finite(v: f64) -> Result<f64, String> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err("nonfinite PLS calculation".into())
    }
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(x, y)| x * y).sum()
}
fn norm(a: &[f64]) -> f64 {
    a.iter().fold(0_f64, |n, v| n.hypot(*v))
}
fn numbers(a: &[f64]) -> Result<String, String> {
    a.iter()
        .map(|v| Ok(finite(*v)?.to_string()))
        .collect::<Result<Vec<_>, String>>()
        .map(|v| format!("[{}]", v.join(",")))
}
fn matrix(a: &[Vec<f64>]) -> Result<String, String> {
    a.iter()
        .map(|r| numbers(r))
        .collect::<Result<Vec<_>, _>>()
        .map(|v| format!("[{}]", v.join(",")))
}
fn transpose(columns: &[Vec<f64>], rows: usize) -> Matrix {
    (0..rows)
        .map(|i| columns.iter().map(|c| c[i]).collect())
        .collect()
}
fn validate(data: &[Vec<f64>]) -> Result<usize, String> {
    let p = data.first().map_or(0, Vec::len);
    if p == 0
        || data.len() > 100_000
        || data.len() > 1_000_000 / p
        || data
            .iter()
            .any(|r| r.len() != p || r.iter().any(|v| !v.is_finite() || v.abs() > 1e150))
    {
        return Err(
            "PLS needs nonempty rectangular finite data; at most100000rows/1Mcells and abs1e150"
                .into(),
        );
    }
    Ok(p)
}
#[derive(Clone)]
struct Transform {
    origin: Vec<f64>,
    offset: Vec<f64>,
    scale: Vec<f64>,
}
impl Transform {
    fn indicators(data: &Matrix, scale: bool) -> Result<(Self, Matrix), String> {
        let n = data.len();
        let q = data[0].len();
        // Indicators have a known zero origin. Exact integer class counts avoid
        // order-dependent online-mean roundoff and 1 + (-2/3) cancellation,
        // which otherwise breaks genuine equal-prior ties in baseline PLS-DA.
        let counts = (0..q)
            .map(|j| data.iter().filter(|r| r[j] == 1.).count() as f64)
            .collect::<Vec<_>>();
        let offset = counts.iter().map(|v| v / n as f64).collect::<Vec<_>>();
        let sd = counts
            .iter()
            .map(|v| {
                if scale {
                    (v * (n as f64 - v) / (n as f64 * (n - 1) as f64)).sqrt()
                } else {
                    1.
                }
            })
            .collect::<Vec<_>>();
        if sd.iter().any(|v| *v <= 0.) {
            return Err("indicator training columns must have positive variance".into());
        }
        let t = Self {
            origin: vec![0.; q],
            offset,
            scale: sd,
        };
        let z = data.iter().map(|r| t.apply(r)).collect::<Result<_, _>>()?;
        Ok((t, z))
    }
    fn fit(data: &Matrix, scale: bool) -> Result<(Self, Matrix), String> {
        let n = data.len();
        let p = data[0].len();
        let origin = data[0].clone();
        let mut offset = vec![0.; p];
        for (i, row) in data.iter().enumerate() {
            for j in 0..p {
                offset[j] += ((row[j] - origin[j]) - offset[j]) / (i + 1) as f64;
            }
        }
        let mut sd = vec![1.; p];
        if scale {
            for j in 0..p {
                let v = data
                    .iter()
                    .fold(0_f64, |s, r| s.hypot((r[j] - origin[j]) - offset[j]))
                    / (n as f64 - 1.).sqrt();
                sd[j] = if v == 0. { 1. } else { finite(v)? };
            }
        }
        let t = Self {
            origin,
            offset,
            scale: sd,
        };
        let z = data.iter().map(|r| t.apply(r)).collect::<Result<_, _>>()?;
        Ok((t, z))
    }
    fn apply(&self, row: &[f64]) -> Result<Vec<f64>, String> {
        row.iter()
            .enumerate()
            .map(|(j, v)| finite(((v - self.origin[j]) - self.offset[j]) / self.scale[j]))
            .collect()
    }
    fn encode(&self, v: &mut Vec<f64>) {
        v.extend(&self.origin);
        v.extend(&self.offset);
        v.extend(&self.scale);
    }
    fn decode(v: &[f64], cursor: &mut usize, p: usize) -> Result<Self, String> {
        let origin = v[*cursor..*cursor + p].to_vec();
        *cursor += p;
        let offset = v[*cursor..*cursor + p].to_vec();
        *cursor += p;
        let scale = v[*cursor..*cursor + p].to_vec();
        *cursor += p;
        if scale.iter().any(|s| *s <= 0.) {
            return Err("PLS saved scales must be positive".into());
        }
        Ok(Self {
            origin,
            offset,
            scale,
        })
    }
}
struct Model {
    da: bool,
    p: usize,
    q: usize,
    a: usize,
    n: usize,
    confidence: f64,
    has_inference: bool,
    x: Transform,
    y: Transform,
    b: Matrix,
    r: Matrix,
    tt: Vec<f64>,
    residual_sd: Vec<f64>,
    labels: Vec<f64>,
}
impl Model {
    fn encode(&self) -> Vec<f64> {
        let mut v = vec![
            2.,
            f64::from(self.da),
            self.p as f64,
            self.q as f64,
            self.a as f64,
            self.n as f64,
            self.confidence,
            f64::from(self.has_inference),
        ];
        self.x.encode(&mut v);
        self.y.encode(&mut v);
        v.extend(self.b.iter().flatten());
        v.extend(self.r.iter().flatten());
        v.extend(&self.tt);
        v.extend(&self.residual_sd);
        v.extend(&self.labels);
        v
    }
    fn decode(v: &[f64], da: bool) -> Result<Self, String> {
        if v.len() < 8
            || v.iter().any(|x| !x.is_finite())
            || (v[0] != 1. && v[0] != 2.)
            || flag(v[1])? != da
        {
            return Err("PLS saved model version or kind mismatch".into());
        }
        let p = integer(v[2], 1, 128)?;
        let q = integer(v[3], 1, 32)?;
        let n = integer(v[5], 2, 100_000)?;
        let a = integer(v[4], 0, p.min(n - 1))?;
        let confidence = v[6];
        let has_inference = flag(v[7])?;
        if !(0. < confidence && confidence < 1.)
            || v.len() != 8 + 3 * p + 4 * q + p * q + p * a + a + if da { q } else { 0 }
            || has_inference != (!da && n > a + 1)
        {
            return Err("PLS saved model length or inference flags invalid".into());
        }
        let mut c = 8;
        let x = Transform::decode(v, &mut c, p)?;
        let y = Transform::decode(v, &mut c, q)?;
        let b = (0..p)
            .map(|_| {
                let z = v[c..c + q].to_vec();
                c += q;
                z
            })
            .collect();
        let r = (0..p)
            .map(|_| {
                let z = v[c..c + a].to_vec();
                c += a;
                z
            })
            .collect();
        let tt = v[c..c + a].to_vec();
        c += a;
        let mut residual_sd = v[c..c + q].to_vec();
        c += q;
        let labels = v[c..].to_vec();
        if tt.iter().any(|t| *t <= 0.)
            || residual_sd
                .iter()
                .any(|s| *s < 0. || (!has_inference && *s != 0.))
            || (da && (q < 2 || labels.windows(2).any(|r| r[0] >= r[1])))
        {
            return Err("invalid PLS score norm, residual uncertainty or class labels".into());
        }
        // Version1 stored variances; version2 stores residual standard deviations
        // so representable standard errors survive squared-variance underflow.
        if v[0] == 1. {
            for sd in &mut residual_sd {
                *sd = sd.sqrt();
            }
        }
        Ok(Self {
            da,
            p,
            q,
            a,
            n,
            confidence,
            has_inference,
            x,
            y,
            b,
            r,
            tt,
            residual_sd,
            labels,
        })
    }
    fn predict(&self, data: &[Vec<f64>], budget: &mut Budget) -> Result<String, String> {
        budget.take(data.len() * (self.p * (self.q + self.a) + self.q * 10))?;
        let z = data
            .iter()
            .map(|r| self.x.apply(r))
            .collect::<Result<Vec<_>, _>>()?;
        let scores = z
            .iter()
            .map(|r| {
                (0..self.a)
                    .map(|a| finite((0..self.p).map(|j| r[j] * self.r[j][a]).sum()))
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Matrix, _>>()?;
        let predictions = z
            .iter()
            .map(|r| {
                (0..self.q)
                    .map(|k| {
                        finite(
                            self.y.origin[k]
                                + (self.y.offset[k]
                                    + (0..self.p).map(|j| r[j] * self.b[j][k]).sum::<f64>()
                                        * self.y.scale[k]),
                        )
                    })
                    .collect::<Result<Vec<_>, _>>()
            })
            .collect::<Result<Matrix, _>>()?;
        let mut out = format!(
            "{{\"version\":1,\"algorithm\":\"{}\",\"predictions\":{},\"x_scores\":{}",
            if self.da {
                "plsda_predict"
            } else {
                "pls2_predict"
            },
            matrix(&predictions)?,
            matrix(&scores)?
        );
        if self.da {
            let predicted = predictions
                .iter()
                .map(|r| self.labels[argmax(r)])
                .collect::<Vec<_>>();
            out.push_str(&format!(",\"class_labels\":{},\"predicted_labels\":{},\"score_semantics\":\"raw_indicator_predictions_not_calibrated_probabilities\"",numbers(&self.labels)?,numbers(&predicted)?));
        } else if self.has_inference {
            let critical = student_quantile(self.confidence, (self.n - self.a - 1) as f64)?;
            let mut meanse = Vec::new();
            let mut indivse = Vec::new();
            let mut lower = Vec::new();
            let mut upper = Vec::new();
            let mut plower = Vec::new();
            let mut pupper = Vec::new();
            let mut leverage = Vec::new();
            for (s, y) in scores.iter().zip(&predictions) {
                let h = finite(
                    1. / self.n as f64
                        + s.iter()
                            .zip(&self.tt)
                            .map(|(v, tt)| (v / tt.sqrt()).powi(2))
                            .sum::<f64>(),
                )?;
                leverage.push(h);
                let m = self
                    .residual_sd
                    .iter()
                    .map(|v| finite(v * h.sqrt()))
                    .collect::<Result<Vec<_>, _>>()?;
                let p = self
                    .residual_sd
                    .iter()
                    .map(|v| finite(v * (1. + h).sqrt()))
                    .collect::<Result<Vec<_>, _>>()?;
                lower.push(y.iter().zip(&m).map(|(v, e)| v - critical * e).collect());
                upper.push(y.iter().zip(&m).map(|(v, e)| v + critical * e).collect());
                plower.push(y.iter().zip(&p).map(|(v, e)| v - critical * e).collect());
                pupper.push(y.iter().zip(&p).map(|(v, e)| v + critical * e).collect());
                meanse.push(m);
                indivse.push(p);
            }
            out.push_str(&format!(",\"leverage\":{},\"mean_standard_error\":{},\"individual_standard_error\":{},\"mean_ci_lower\":{},\"mean_ci_upper\":{},\"individual_ci_lower\":{},\"individual_ci_upper\":{},\"uncertainty\":\"normality_based_conditional_score_approximation\"",numbers(&leverage)?,matrix(&meanse)?,matrix(&indivse)?,matrix(&lower)?,matrix(&upper)?,matrix(&plower)?,matrix(&pupper)?));
        } else {
            out.push_str(",\"mean_standard_error\":null,\"individual_standard_error\":null,\"mean_ci_lower\":null,\"mean_ci_upper\":null,\"individual_ci_lower\":null,\"individual_ci_upper\":null");
        }
        out.push('}');
        Ok(out)
    }
}
fn argmax(v: &[f64]) -> usize {
    let mut best = 0;
    for j in 1..v.len() {
        if v[j] > v[best] {
            best = j;
        }
    }
    best
}
fn cross(x: &Matrix, y: &Matrix, budget: &mut Budget) -> Result<Matrix, String> {
    let p = x[0].len();
    let q = y[0].len();
    budget.take(x.len() * p * q * 2)?;
    (0..p)
        .map(|j| {
            (0..q)
                .map(|k| finite(x.iter().zip(y).map(|(a, b)| a[j] * b[k]).sum()))
                .collect()
        })
        .collect()
}
fn leading(s: &Matrix, budget: &mut Budget) -> Result<Option<Vec<f64>>, String> {
    let p = s.len();
    let q = s[0].len();
    let largest = s.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
    if largest == 0. {
        return Ok(None);
    }
    let d = p.min(q);
    budget.take(p * q * d + 100 * d * d * d)?;
    let a = s
        .iter()
        .map(|r| r.iter().map(|v| v / largest).collect::<Vec<_>>())
        .collect::<Matrix>();
    let gram = (0..d)
        .map(|i| {
            (0..d)
                .map(|j| {
                    if q <= p {
                        (0..p).map(|k| a[k][i] * a[k][j]).sum()
                    } else {
                        dot(&a[i], &a[j])
                    }
                })
                .collect()
        })
        .collect();
    let (e, v) = symmetric_eigen(gram)?;
    let k = argmax(&e);
    let mut w = if q <= p {
        (0..p)
            .map(|j| (0..q).map(|i| a[j][i] * v[i][k]).sum())
            .collect::<Vec<_>>()
    } else {
        (0..p).map(|j| v[j][k]).collect()
    };
    let length = norm(&w);
    if length == 0. {
        return Ok(None);
    }
    for v in &mut w {
        *v /= length;
    }
    let pivot = argmax(&w.iter().map(|v| v.abs()).collect::<Vec<_>>());
    if w[pivot] < 0. {
        for v in &mut w {
            *v = -*v;
        }
    }
    Ok(Some(w))
}
fn train(
    data: &[Vec<f64>],
    args: &[f64],
    da: bool,
    budget: &mut Budget,
    report: bool,
) -> Result<(Model, String), String> {
    let n = data.len();
    let width = data[0].len();
    if n < 2 || args.len() < if da { 1 } else { 2 } || args.len() > if da { 4 } else { 6 } {
        return Err("PLS fit argument count or n invalid".into());
    }
    let q = if da { 0 } else { integer(args[0], 1, 32)? };
    let p = width
        .checked_sub(if da { 1 } else { q })
        .ok_or("PLS response width exceeds data")?;
    if p == 0 || p > 128 {
        return Err("PLS needs1..128 numeric predictors".into());
    }
    let k = usize::from(!da);
    let requested = integer(args[k], usize::from(report), p.min(n - 1))?;
    let method = integer(args.get(k + 1).copied().unwrap_or(0.), 0, 1)?;
    let sx = flag(args.get(k + 2).copied().unwrap_or(1.))?;
    let sy = flag(args.get(k + 3).copied().unwrap_or(1.))?;
    let confidence = if da {
        0.95
    } else {
        args.get(5).copied().unwrap_or(0.95)
    };
    if !(0. < confidence && confidence < 1.) {
        return Err("PLS confidence must be strictly(0,1)".into());
    }
    let mut labels = if da {
        data.iter().map(|r| r[p]).collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    labels.sort_by(f64::total_cmp);
    labels.dedup_by(|a, b| *a == *b);
    let q = if da { labels.len() } else { q };
    if q > 32 || (da && q < 2) {
        return Err("PLS-DA needs2..32 classes".into());
    }
    let rawx = data.iter().map(|r| r[..p].to_vec()).collect::<Matrix>();
    let rawy = data
        .iter()
        .map(|r| {
            if da {
                labels.iter().map(|v| f64::from(*v == r[p])).collect()
            } else {
                r[p..].to_vec()
            }
        })
        .collect::<Matrix>();
    budget.take(n * (p + q) * 6)?;
    let (xt, x) = Transform::fit(&rawx, sx)?;
    let (yt, y) = if da {
        Transform::indicators(&rawy, sy)?
    } else {
        Transform::fit(&rawy, sy)?
    };
    let mut xr = x.clone();
    let mut yr = y.clone();
    let mut s = cross(&x, &y, budget)?;
    let initial = s.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max);
    let xss = x.iter().flatten().map(|v| v * v).sum::<f64>();
    let yss = y.iter().flatten().map(|v| v * v).sum::<f64>();
    finite(xss)?;
    finite(yss)?;
    let mut rotations: Matrix = Vec::new();
    let mut loadings: Matrix = Vec::new();
    let mut yloadings: Matrix = Vec::new();
    let mut weights: Matrix = Vec::new();
    let mut scores: Matrix = Vec::new();
    let mut orthogonal: Matrix = Vec::new();
    let mut tt = Vec::new();
    let mut explainedx = Vec::new();
    let mut explainedy = Vec::new();
    let mut b = vec![vec![0.; q]; p];
    for _ in 0..requested {
        if s.iter().flatten().map(|v| v.abs()).fold(0_f64, f64::max) <= initial * 1e-12 {
            break;
        }
        let Some(w) = leading(&s, budget)? else {
            break;
        };
        budget.take(n * (p + q) * 10 + p * p * 4)?;
        let mut r = w.clone();
        if method == 0 {
            for (oldr, oldp) in rotations.iter().zip(&loadings) {
                let z = dot(oldp, &w);
                for j in 0..p {
                    r[j] -= oldr[j] * z;
                }
            }
        }
        let mut t = x.iter().map(|row| dot(row, &r)).collect::<Vec<_>>();
        // Reorthogonalize scores and projections together to keep t=X*r.
        for _ in 0..2 {
            for ((oldr, oldt), oldtt) in rotations.iter().zip(&scores).zip(&tt) {
                let z = dot(oldt, &t) / oldtt;
                for j in 0..p {
                    r[j] -= oldr[j] * z;
                }
                for i in 0..n {
                    t[i] -= oldt[i] * z;
                }
            }
        }
        let mut tnorm = dot(&t, &t);
        if tnorm <= xss * 1e-24 {
            break;
        }
        if method == 1 {
            let norm = tnorm.sqrt();
            for v in &mut t {
                *v /= norm;
            }
            for v in &mut r {
                *v /= norm;
            }
            tnorm = dot(&t, &t);
        }
        let xp = if method == 0 { &xr } else { &x };
        let yp = if method == 0 { &yr } else { &y };
        let loading = (0..p)
            .map(|j| xp.iter().zip(&t).map(|(row, v)| row[j] * v).sum::<f64>() / tnorm)
            .collect::<Vec<_>>();
        let c = (0..q)
            .map(|j| yp.iter().zip(&t).map(|(row, v)| row[j] * v).sum::<f64>() / tnorm)
            .collect::<Vec<_>>();
        if method == 1 {
            let mut v = loading.clone();
            for _ in 0..2 {
                for old in &orthogonal {
                    let z = dot(old, &v);
                    for j in 0..p {
                        v[j] -= old[j] * z;
                    }
                }
            }
            let length = norm(&v);
            if length <= norm(&loading) * 1e-12 {
                break;
            }
            for v in &mut v {
                *v /= length;
            }
            let projection = (0..q)
                .map(|k| (0..p).map(|j| v[j] * s[j][k]).sum::<f64>())
                .collect::<Vec<_>>();
            for (row, weight) in s.iter_mut().zip(&v) {
                for (cell, z) in row.iter_mut().zip(&projection) {
                    *cell -= weight * z;
                }
            }
            orthogonal.push(v);
        }
        for j in 0..p {
            for k in 0..q {
                b[j][k] = finite(b[j][k] + r[j] * c[k])?;
            }
        }
        explainedx.push(tnorm * dot(&loading, &loading));
        explainedy.push(tnorm * dot(&c, &c));
        let mut vipweight = if method == 0 { w } else { r.clone() };
        let wn = norm(&vipweight);
        for v in &mut vipweight {
            *v /= wn;
        }
        rotations.push(r);
        loadings.push(loading.clone());
        yloadings.push(c.clone());
        weights.push(vipweight);
        scores.push(t.clone());
        tt.push(tnorm);
        if method == 0 {
            for i in 0..n {
                for j in 0..p {
                    xr[i][j] -= t[i] * loading[j];
                }
                for k in 0..q {
                    yr[i][k] -= t[i] * c[k];
                }
            }
            s = cross(&xr, &yr, budget)?;
        }
    }
    let a = tt.len();
    let fitted = (0..n)
        .map(|i| {
            (0..q)
                .map(|k| {
                    finite(
                        yt.origin[k]
                            + (yt.offset[k]
                                + dot(&x[i], &b.iter().map(|r| r[k]).collect::<Vec<_>>())
                                    * yt.scale[k]),
                    )
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .collect::<Result<Matrix, _>>()?;
    let mut sse = vec![0.; q];
    let mut residual_norms = vec![0.; q];
    for k in 0..q {
        let length = norm(
            &(0..n)
                .map(|i| {
                    ((rawy[i][k] - yt.origin[k]) - yt.offset[k])
                        - dot(&x[i], &b.iter().map(|r| r[k]).collect::<Vec<_>>()) * yt.scale[k]
                })
                .collect::<Vec<_>>(),
        );
        residual_norms[k] = finite(length)?;
        sse[k] = finite(length * length)?;
    }
    let has_inference = !da && n > a + 1;
    let residual_sd = if has_inference {
        residual_norms
            .iter()
            .map(|v| v / ((n - a - 1) as f64).sqrt())
            .collect()
    } else {
        vec![0.; q]
    };
    let model = Model {
        da,
        p,
        q,
        a,
        n,
        confidence,
        has_inference,
        x: xt,
        y: yt,
        b,
        r: transpose(&rotations, p),
        tt,
        residual_sd,
        labels,
    };
    if !report {
        return Ok((model, String::new()));
    }
    let total = explainedy.iter().sum::<f64>();
    let vip = if total > 0. {
        numbers(
            &(0..p)
                .map(|j| {
                    (p as f64
                        * weights
                            .iter()
                            .zip(&explainedy)
                            .map(|(w, ss)| ss * w[j] * w[j])
                            .sum::<f64>()
                        / total)
                        .max(0.)
                        .sqrt()
                })
                .collect::<Vec<_>>(),
        )?
    } else {
        "null".into()
    };
    let coefficients = (0..p)
        .map(|j| {
            (0..q)
                .map(|k| finite(model.b[j][k] * model.y.scale[k] / model.x.scale[j]))
                .collect()
        })
        .collect::<Result<Matrix, _>>()?;
    let intercept = (0..q)
        .map(|k| {
            finite(
                model.y.origin[k] + model.y.offset[k]
                    - (0..p)
                        .map(|j| (model.x.origin[j] + model.x.offset[j]) * coefficients[j][k])
                        .sum::<f64>(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let xratio = if xss > 0. {
        numbers(&explainedx.iter().map(|v| v / xss).collect::<Vec<_>>())?
    } else {
        "null".into()
    };
    let yratio = if yss > 0. {
        numbers(&explainedy.iter().map(|v| v / yss).collect::<Vec<_>>())?
    } else {
        "null".into()
    };
    let mut out=format!("{{\"version\":1,\"algorithm\":\"{}\",\"method\":\"{}\",\"n_samples\":{n},\"n_features\":{p},\"n_responses\":{q},\"components_requested\":{requested},\"components_used\":{a},\"coefficients\":{},\"intercepts\":{},\"coefficients_standardized\":{},\"x_weights\":{},\"x_rotations\":{},\"x_loadings\":{},\"y_loadings\":{},\"x_scores\":{},\"score_squared_norms\":{},\"x_explained_fraction\":{xratio},\"y_explained_fraction\":{yratio},\"vip\":{vip},\"fitted\":{},\"sse\":{},\"residual_df\":{},\"predict_parameters\":{},\"centering\":\"training_offset_safe\",\"scale_x\":{sx},\"scale_y\":{sy},\"scale_convention\":\"sample_ddof1_when_enabled_constant_scale1\",\"confidence\":{confidence},\"full_JMP_parity\":false",if da{"plsda_fit"}else{"pls2_fit"},if method==0{"NIPALS_SVD"}else{"SIMPLS_SVD"},matrix(&coefficients)?,numbers(&intercept)?,matrix(&model.b)?,matrix(&transpose(&weights,p))?,matrix(&model.r)?,matrix(&transpose(&loadings,p))?,matrix(&transpose(&yloadings,q))?,matrix(&transpose(&scores,n))?,numbers(&model.tt)?,matrix(&fitted)?,numbers(&sse)?,n-a-1,numbers(&model.encode())?);
    out.push_str(&format!(
        ",\"residual_standard_deviations\":{}",
        if has_inference {
            numbers(&model.residual_sd)?
        } else {
            "null".into()
        }
    ));
    if da {
        let predicted = fitted
            .iter()
            .map(|r| model.labels[argmax(r)])
            .collect::<Vec<_>>();
        let accuracy = predicted
            .iter()
            .zip(data)
            .filter(|(pred, row)| **pred == row[p])
            .count() as f64
            / n as f64;
        out.push_str(&format!(",\"class_labels\":{},\"predicted_labels\":{},\"training_accuracy\":{accuracy},\"score_semantics\":\"raw_indicator_predictions_not_calibrated_probabilities\"",numbers(&model.labels)?,numbers(&predicted)?));
    }
    out.push('}');
    Ok((model, out))
}
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let width = validate(data)?;
    if args.iter().any(|v| !v.is_finite()) {
        return Err("PLS parameters must be finite".into());
    }
    let mut budget = Budget(0);
    match operation {
        "pls2_fit" | "plsda_fit" => {
            train(data, args, operation == "plsda_fit", &mut budget, true).map(|(_, out)| out)
        }
        "pls2_validate" | "plsda_validate" => {
            validate_components(data, args, operation == "plsda_validate", &mut budget)
        }
        "pls2_predict" | "plsda_predict" => {
            let m = Model::decode(args, operation == "plsda_predict")?;
            if width != m.p {
                return Err("PLS prediction feature width mismatch".into());
            }
            m.predict(data, &mut budget)
        }
        _ => Err("unknown multivariate PLS operation".into()),
    }
}

fn evaluate_row(model: &Model, row: &[f64]) -> Result<(Vec<f64>, Vec<f64>), String> {
    let z = model.x.apply(&row[..model.p])?;
    let mut predictions = Vec::with_capacity(model.q);
    let mut errors = Vec::with_capacity(model.q);
    for k in 0..model.q {
        let delta =
            finite((0..model.p).map(|j| z[j] * model.b[j][k]).sum::<f64>() * model.y.scale[k])?;
        predictions.push(finite(model.y.origin[k] + (model.y.offset[k] + delta))?);
        let actual = if model.da {
            f64::from(row[model.p] == model.labels[k])
        } else {
            row[model.p + k]
        };
        errors.push(finite(
            ((actual - model.y.origin[k]) - model.y.offset[k]) - delta,
        )?);
    }
    Ok((predictions, errors))
}
fn nullable_matrix(a: &[Vec<Option<f64>>]) -> Result<String, String> {
    a.iter()
        .map(|row| {
            row.iter()
                .map(|v| match v {
                    Some(x) => Ok(finite(*x)?.to_string()),
                    None => Ok("null".into()),
                })
                .collect::<Result<Vec<_>, String>>()
                .map(|v| format!("[{}]", v.join(",")))
        })
        .collect::<Result<Vec<_>, _>>()
        .map(|v| format!("[{}]", v.join(",")))
}
fn validate_components(
    data: &[Vec<f64>],
    args: &[f64],
    da: bool,
    budget: &mut Budget,
) -> Result<String, String> {
    use std::collections::BTreeMap;
    let k = usize::from(!da);
    let width = data[0].len();
    let n = data.len();
    if args.len() < k + 1 || args.len() > k + 5 {
        return Err("PLS validation parameter count invalid".into());
    }
    let mut q = if da { 1 } else { integer(args[0], 1, 32)? };
    let p = width
        .checked_sub(q + 1)
        .ok_or("PLS validation response/role width invalid")?;
    if !(1..=128).contains(&p) {
        return Err("PLS validation needs1..128 predictors".into());
    }
    let requested = integer(args[k], 1, p.min(32))?;
    let method = integer(args.get(k + 1).copied().unwrap_or(0.), 0, 1)?;
    let sx = flag(args.get(k + 2).copied().unwrap_or(1.))?;
    let sy = flag(args.get(k + 3).copied().unwrap_or(1.))?;
    let mode = integer(args.get(k + 4).copied().unwrap_or(0.), 0, 1)?;
    let ids = data
        .iter()
        .map(|r| {
            let v = r[width - 1];
            if v.fract() != 0.
                || !(0. ..=9_007_199_254_740_991.).contains(&v)
                || (mode == 1 && v > 2.)
            {
                Err("validation fold/role must be an exact integer in its domain".into())
            } else {
                Ok(v as u64)
            }
        })
        .collect::<Result<Vec<_>, String>>()?;
    let reference_indices = (0..n)
        .filter(|i| mode == 0 || ids[*i] == 0)
        .collect::<Vec<_>>();
    let validation_indices = (0..n)
        .filter(|i| mode == 0 || ids[*i] == 1)
        .collect::<Vec<_>>();
    let test_indices = (0..n)
        .filter(|i| mode == 1 && ids[*i] == 2)
        .collect::<Vec<_>>();
    if reference_indices.len() < 2 || validation_indices.is_empty() {
        return Err("validation needs>=2training and>=1validation rows".into());
    }
    let mut labels = if da {
        reference_indices
            .iter()
            .map(|i| data[*i][p])
            .collect::<Vec<_>>()
    } else {
        Vec::new()
    };
    labels.sort_by(f64::total_cmp);
    labels.dedup_by(|a, b| *a == *b);
    if da {
        q = labels.len();
        if !(2..=32).contains(&q) || data.iter().any(|r| !labels.contains(&r[p])) {
            return Err(
                "PLS validation classes must be2..32 training levels and cover evaluation labels"
                    .into(),
            );
        }
    }
    if (requested + 1) * validation_indices.len() * q > 2_000_000 {
        return Err("PLS validation path exceeds2Mprediction cells".into());
    }
    let mut groups = BTreeMap::<u64, Vec<usize>>::new();
    for i in &validation_indices {
        groups
            .entry(if mode == 0 { ids[*i] } else { 1 })
            .or_default()
            .push(*i);
    }
    if mode == 0 && !(2..=1000).contains(&groups.len()) {
        return Err("K-fold needs2..1000nonempty groups".into());
    }
    for (id, valid) in &groups {
        let ntrain = if mode == 0 {
            n - valid.len()
        } else {
            reference_indices.len()
        };
        if ntrain < 2 || requested > ntrain - 1 {
            return Err(format!(
                "PLS validation group{id} has insufficient training rows for requested factors"
            ));
        }
    }
    budget.take((reference_indices.len() + validation_indices.len()) * q * 6)?;
    let actual = |i: usize| -> Vec<f64> {
        if da {
            labels.iter().map(|v| f64::from(*v == data[i][p])).collect()
        } else {
            data[i][p..width - 1].to_vec()
        }
    };
    let reference_y = reference_indices
        .iter()
        .map(|i| actual(*i))
        .collect::<Matrix>();
    let (_, centered) = Transform::fit(&reference_y, false)?;
    let reference_sd = (0..q)
        .map(|j| {
            finite(
                centered.iter().fold(0_f64, |s, r| s.hypot(r[j]))
                    / ((reference_y.len() - 1) as f64).sqrt(),
            )
        })
        .collect::<Result<Vec<_>, _>>()?;
    let variances = reference_sd
        .iter()
        .map(|v| finite(v * v))
        .collect::<Result<Vec<_>, _>>()?;
    if variances.iter().any(|v| *v <= 0.) {
        return Err("PLS rootmeanPRESS reference response variance must be positive".into());
    }
    let validation_y = validation_indices
        .iter()
        .map(|i| actual(*i))
        .collect::<Matrix>();
    let sst = if validation_y.len() > 1 {
        let (_, y) = Transform::fit(&validation_y, false)?;
        (0..q)
            .map(|j| y.iter().fold(0_f64, |s, r| s.hypot(r[j] / reference_sd[j])))
            .collect::<Vec<f64>>()
    } else {
        vec![0.; q]
    };
    let mut positions = vec![usize::MAX; n];
    for (j, i) in validation_indices.iter().enumerate() {
        positions[*i] = j;
    }
    let mut path = Vec::new();
    let mut press = Vec::new();
    let mut normalized = Vec::new();
    let mut root = Vec::new();
    let mut q2 = Vec::new();
    let mut accuracy = Vec::new();
    let mut provenance = Vec::new();
    let mut used_path = Vec::new();
    for component in 0..=requested {
        let mut predictions = vec![vec![0.; q]; validation_indices.len()];
        let mut raw = vec![0_f64; q];
        let mut standardized = vec![0_f64; q];
        let mut normalized_raw = vec![0_f64; q];
        let mut used = Vec::new();
        for (id, valid) in &groups {
            budget.take(n * width)?;
            let train_rows = (0..n)
                .filter(|i| {
                    if mode == 0 {
                        ids[*i] != *id
                    } else {
                        ids[*i] == 0
                    }
                })
                .map(|i| data[i][..width - 1].to_vec())
                .collect::<Matrix>();
            let config = if da {
                vec![
                    component as f64,
                    method as f64,
                    f64::from(sx),
                    f64::from(sy),
                ]
            } else {
                vec![
                    q as f64,
                    component as f64,
                    method as f64,
                    f64::from(sx),
                    f64::from(sy),
                    0.95,
                ]
            };
            let (model, _) = train(&train_rows, &config, da, budget, false)?;
            if da && model.labels != labels {
                return Err(
                    "every PLS-DA training split must contain every evaluation class".into(),
                );
            }
            used.push(model.a as f64);
            budget.take(valid.len() * p * q * 2)?;
            let mut group_press = vec![0_f64; q];
            let mut group_normalized = vec![0_f64; q];
            for i in valid {
                let (pred, errors) = evaluate_row(&model, &data[*i])?;
                predictions[positions[*i]] = pred;
                for j in 0..q {
                    group_press[j] = finite(group_press[j].hypot(errors[j]))?;
                    group_normalized[j] =
                        finite(group_normalized[j].hypot(errors[j] / reference_sd[j]))?;
                }
            }
            for j in 0..q {
                raw[j] = finite(raw[j].hypot(group_press[j]))?;
                standardized[j] = finite(
                    standardized[j].hypot(group_normalized[j] / (valid.len() as f64).sqrt()),
                )?;
                normalized_raw[j] = finite(normalized_raw[j].hypot(group_normalized[j]))?;
            }
            if component == 0 {
                provenance.push(format!("{{\"group\":{id},\"training_count\":{},\"validation_count\":{},\"feature_origins\":{},\"feature_mean_offsets\":{},\"feature_scales\":{},\"response_origins\":{},\"response_mean_offsets\":{},\"response_scales\":{}}}",model.n,valid.len(),numbers(&model.x.origin)?,numbers(&model.x.offset)?,numbers(&model.x.scale)?,numbers(&model.y.origin)?,numbers(&model.y.offset)?,numbers(&model.y.scale)?));
            }
        }
        let divisor = if groups.len() > 1 {
            (groups.len() - 1) as f64
        } else {
            1.
        };
        // Keep norms until after normalization: squaring a tiny raw response
        // residual first can erase a representable standardized PRESS/root.
        root.push(finite(norm(&standardized) / (divisor * q as f64).sqrt())?);
        q2.push(
            normalized_raw
                .iter()
                .zip(&sst)
                .map(|(v, ss)| {
                    if *ss > 0. {
                        Some(1. - (v / ss).powi(2))
                    } else {
                        None
                    }
                })
                .collect::<Vec<_>>(),
        );
        for v in &mut raw {
            *v = finite(*v * *v)?;
        }
        for v in &mut standardized {
            *v = finite((*v / divisor.sqrt()).powi(2))?;
        }
        if da {
            accuracy.push(
                validation_indices
                    .iter()
                    .zip(&predictions)
                    .filter(|(i, pred)| labels[argmax(pred)] == data[**i][p])
                    .count() as f64
                    / validation_indices.len() as f64,
            );
        }
        path.push(predictions);
        press.push(raw);
        normalized.push(standardized);
        used_path.push(used);
    }
    let selected = argmax(&root.iter().map(|v| -*v).collect::<Vec<_>>());
    let train_rows = reference_indices
        .iter()
        .map(|i| data[*i][..width - 1].to_vec())
        .collect::<Matrix>();
    budget.take(train_rows.len() * width)?;
    let config = if da {
        vec![selected as f64, method as f64, f64::from(sx), f64::from(sy)]
    } else {
        vec![
            q as f64,
            selected as f64,
            method as f64,
            f64::from(sx),
            f64::from(sy),
            0.95,
        ]
    };
    let (model, _) = train(&train_rows, &config, da, budget, false)?;
    let mut tests = Vec::new();
    let mut testpress = vec![0.; q];
    budget.take(test_indices.len() * p * q * 2)?;
    for i in &test_indices {
        let (pred, errors) = evaluate_row(&model, &data[*i])?;
        tests.push(pred);
        for j in 0..q {
            testpress[j] = finite(testpress[j] + errors[j] * errors[j])?;
        }
    }
    let cube = path
        .iter()
        .map(|a| matrix(a))
        .collect::<Result<Vec<_>, _>>()?
        .join(",");
    let mut out=format!("{{\"version\":1,\"algorithm\":\"{}\",\"mode\":\"{}\",\"components\":{},\"components_used_by_group\":{},\"selected_components\":{selected},\"selected_model_components_used\":{},\"validation_indices\":{},\"test_indices\":{},\"predictions_by_components\":[{cube}],\"raw_press_per_response\":{},\"normalized_press_per_response\":{},\"root_mean_press\":{},\"pooled_predictive_r2_per_response\":{},\"reference_response_variances\":{},\"split_training_transforms\":[{}],\"test_predictions\":{},\"test_press_per_response\":{},\"predict_parameters\":{},\"selection_rule\":\"minimum_documented_root_mean_PRESS_smallest_exact_tie\",\"variance_reference\":\"{}\",\"training_transform_scope\":\"each_training_split_only\",\"pooled_R2_is_not_asserted_JMP_Q2\":true,\"full_JMP_parity\":false",if da{"plsda_validate"}else{"pls2_validate"},if mode==0{"supplied_folds"}else{"supplied_roles"},numbers(&(0..=requested).map(|v|v as f64).collect::<Vec<_>>())?,matrix(&used_path)?,model.a,numbers(&validation_indices.iter().map(|v|*v as f64).collect::<Vec<_>>())?,numbers(&test_indices.iter().map(|v|*v as f64).collect::<Vec<_>>())?,matrix(&press)?,matrix(&normalized)?,numbers(&root)?,nullable_matrix(&q2)?,numbers(&variances)?,provenance.join(","),matrix(&tests)?,if test_indices.is_empty(){"null".into()}else{numbers(&testpress)?},numbers(&model.encode())?,if mode==0{"whole_response_for_evaluation_only"}else{"training_role0_only"});
    if da {
        let predicted = path[selected]
            .iter()
            .map(|r| labels[argmax(r)])
            .collect::<Vec<_>>();
        out.push_str(&format!(",\"class_labels\":{},\"validation_accuracy_by_components\":{},\"selected_validation_labels\":{},\"score_semantics\":\"raw_indicator_predictions_not_calibrated_probabilities\"",numbers(&labels)?,numbers(&accuracy)?,numbers(&predicted)?));
    }
    out.push('}');
    Ok(out)
}
