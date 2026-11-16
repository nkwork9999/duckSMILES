//! Explicit-input materials descriptors and liquid-viscosity correlations.
//! No phase, ionic radius, pressure, or non-Newtonian behavior is inferred.
fn nums(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn positive(x: f64) -> Result<f64, String> {
    if x.is_finite() && x > 0. {
        Ok(x)
    } else {
        Err("expected a finite positive value".into())
    }
}
fn fractions(d: &[Vec<f64>]) -> Result<(), String> {
    if d.iter().any(|r| !(0. ..=1.).contains(&r[0]))
        || (d.iter().map(|r| r[0]).sum::<f64>() - 1.).abs() > 1e-10
    {
        return Err(
            "mole fractions must be in [0,1] and sum to one; no implicit normalization".into(),
        );
    }
    Ok(())
}
pub fn run(op: &str, d: &[Vec<f64>], p: &[f64]) -> Result<String, String> {
    let w = d.first().map_or(0, Vec::len);
    if d.is_empty()
        || d.len() > 10000
        || w == 0
        || w > 128
        || d.iter()
            .any(|r| r.len() != w || r.iter().any(|x| !x.is_finite()))
        || p.iter().any(|x| !x.is_finite())
    {
        return Err(
            "materials methods require finite rectangular data, <=10000 rows, <=128 columns".into(),
        );
    }
    let out = match op {
        "perovskite_descriptors" => {
            if w != 6 || !p.is_empty() {
                return Err("expected [rA,rB,rX,zA,zB,zX] rows and [] parameters".into());
            }
            let mut rows = Vec::new();
            for r in d {
                for radius in &r[..3] {
                    positive(*radius)?;
                }
                if r[3] <= 0. || r[4] <= 0. || r[5] >= 0. {
                    return Err("ABX3 requires positive cation and negative anion charges".into());
                }
                let t = (r[0] + r[2]) / (std::f64::consts::SQRT_2 * (r[1] + r[2]));
                let mu = r[1] / r[2];
                let charge = r[3] + r[4] + 3. * r[5];
                rows.push(format!("{{\"radii\":{},\"charges\":{},\"goldschmidt_t\":{t},\"octahedral_mu\":{mu},\"charge_sum\":{charge},\"charge_balanced\":{}}}", nums(&r[..3]), nums(&r[3..]), charge.abs() <= 1e-10));
            }
            format!("{{\"schema\":1,\"method\":\"ABX3-explicit-site-v1\",\"radius_unit\":\"consistent_user_unit\",\"stability_prediction\":false,\"rows\":[{}]}}", rows.join(","))
        }
        "mixture_descriptors" => {
            if w < 2 || !p.is_empty() {
                return Err("expected [mole_fraction,features...] rows and [] parameters".into());
            }
            fractions(d)?;
            let means: Vec<f64> = (1..w)
                .map(|j| d.iter().map(|r| r[0] * r[j]).sum())
                .collect();
            let sd: Vec<f64> = (1..w)
                .map(|j| {
                    d.iter()
                        .map(|r| r[0] * (r[j] - means[j - 1]).powi(2))
                        .sum::<f64>()
                        .sqrt()
                })
                .collect();
            let entropy: f64 = d
                .iter()
                .filter(|r| r[0] > 0.)
                .map(|r| -r[0] * r[0].ln())
                .sum();
            format!("{{\"schema\":1,\"method\":\"mole-fraction-moments-v1\",\"means\":{},\"population_std\":{},\"mixing_entropy_over_R\":{entropy}}}", nums(&means), nums(&sd))
        }
        "viscosity_arrhenius" => {
            if w != 1 || p.len() != 2 {
                return Err("expected [temperature_K] rows; [A,B_K], ln(eta/(Pa s))=A+B/T".into());
            }
            let mut eta = Vec::new();
            for r in d {
                eta.push(positive((p[0] + p[1] / positive(r[0])?).exp())?);
            }
            format!("{{\"schema\":1,\"method\":\"Andrade-Arrhenius-v1\",\"temperature_K\":{},\"parameters\":{},\"viscosity_Pa_s\":{}}}", nums(&d.iter().map(|r| r[0]).collect::<Vec<_>>()), nums(p), nums(&eta))
        }
        "viscosity_fit" => {
            if w != 2 || d.len() < 3 || !p.is_empty() {
                return Err(
                    "fit requires >=3 [temperature_K,viscosity_Pa_s] rows and [] parameters".into(),
                );
            }
            let n = d.len() as f64;
            let mut x = Vec::new();
            let mut y = Vec::new();
            for r in d {
                x.push(1. / positive(r[0])?);
                y.push(positive(r[1])?.ln());
            }
            let xm = x.iter().sum::<f64>() / n;
            let ym = y.iter().sum::<f64>() / n;
            let ss: f64 = x.iter().map(|v| (v - xm).powi(2)).sum();
            if ss <= f64::EPSILON * x.iter().map(|v| v * v).sum::<f64>() {
                return Err("temperature range is rank deficient".into());
            }
            let b: f64 = x
                .iter()
                .zip(&y)
                .map(|(x, y)| (x - xm) * (y - ym))
                .sum::<f64>()
                / ss;
            let a = ym - b * xm;
            let residuals: Vec<f64> = x.iter().zip(&y).map(|(x, y)| y - (a + b * x)).collect();
            let sse: f64 = residuals.iter().map(|v| v * v).sum();
            let rmse = (sse / n).sqrt();
            let residual_sd = (sse / (n - 2.)).sqrt();
            let tmin = d.iter().map(|r| r[0]).fold(f64::INFINITY, f64::min);
            let tmax = d.iter().map(|r| r[0]).fold(f64::NEG_INFINITY, f64::max);
            let ea = b * 8.31446261815324;
            format!("{{\"schema\":1,\"method\":\"OLS-ln-Pa-s-vs-inverse-K-v1\",\"n\":{},\"parameters\":[{a},{b}],\"activation_energy_J_mol\":{ea},\"temperature_range_K\":[{tmin},{tmax}],\"rmse_log\":{rmse},\"residual_std_log\":{residual_sd},\"residuals_log\":{}}}",d.len(),nums(&residuals))
        }
        "viscosity_mix" => {
            if w != 2 || !(p.is_empty() || (p.len() == 1 && d.len() == 2)) {
                return Err("expected [mole_fraction,viscosity_Pa_s] rows; [] ideal log mix or [G12] for a binary Grunberg-Nissan mixture".into());
            }
            fractions(d)?;
            let mut ln_eta = 0.;
            for r in d {
                ln_eta += r[0] * positive(r[1])?.ln();
            }
            if !p.is_empty() {
                ln_eta += d[0][0] * d[1][0] * p[0];
            }
            let eta = positive(ln_eta.exp())?;
            format!("{{\"schema\":1,\"method\":\"log-mole-fraction-Grunberg-Nissan-v1\",\"G12\":{},\"viscosity_Pa_s\":{eta},\"same_temperature_and_pressure_required\":true}}",if p.is_empty() {"null".into()} else {p[0].to_string()})
        }
        _ => return Err("unknown materials operation".into()),
    };
    if out
        .split([':', ',', '[', ']', '}', ' '])
        .any(|v| matches!(v, "NaN" | "inf" | "-inf"))
    {
        return Err("materials result exceeds finite range".into());
    }
    Ok(out)
}
