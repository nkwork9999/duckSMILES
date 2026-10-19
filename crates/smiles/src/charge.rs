//! Pinned RDKit2026.03.6 Gasteiger-Marsili 12-iteration charge model.
//! BSD-3-Clause upstream notice applies; see RDKIT-DATA-LICENSE.txt.
use crate::electronic::Context;
use crate::parser::Molecule;
fn charges(mol: &Molecule, context: &Context) -> (Vec<f64>, Vec<f64>, bool) {
    use crate::electronic::Hybridization as H;
    let n = mol.atoms.len();
    let adj = mol.adjacency();
    let mut q = vec![0f64; n];
    let mut hq = vec![0f64; n];
    // Same-element sites separated by two conjugated bonds share the initial
    // formal charge; preserve repeated markers and input traversal order.
    for (i, a) in mol.atoms.iter().enumerate() {
        let mut formal = a.charge as f64;
        if formal.abs() > f64::EPSILON && q[i].abs() < f64::EPSILON {
            let mut marker = vec![i];
            for &(middle, b1) in &adj[i] {
                if !context.conjugated[b1] {
                    continue;
                }
                for &(j, b2) in &adj[middle] {
                    if b1 != b2 && context.conjugated[b2] && context.z[i] == context.z[j] {
                        formal += mol.atoms[j].charge as f64;
                        marker.push(j);
                    }
                }
            }
            let value = formal / marker.len() as f64;
            for j in marker {
                q[j] = value;
            }
        }
    }
    let parameters: Vec<[f64; 3]> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            let mode = match context.hybridization[i] {
                H::Sp3 => "sp3",
                H::Sp2 => "sp2",
                H::Sp => "sp",
                _ => {
                    if context.z[i] == 1 {
                        "*"
                    } else if context.z[i] == 16 {
                        match adj[i].iter().filter(|&&(j, _)| context.z[j] == 8).count() {
                            2 => "so2",
                            1 => "so",
                            _ => "sp3",
                        }
                    } else {
                        ""
                    }
                }
            };
            crate::charge_data::parameters(&a.symbol, mode)
        })
        .collect();
    let ion: Vec<f64> = parameters
        .iter()
        .enumerate()
        .map(|(i, p)| {
            if context.z[i] == 1 {
                20.02
            } else {
                p[0] + p[1] + p[2]
            }
        })
        .collect();
    let hp = crate::charge_data::parameters("H", "*");
    let mut energy = vec![0.; n];
    let mut damp = 0.5;
    for _ in 0..12 {
        for i in 0..n {
            energy[i] = parameters[i][0] + q[i] * (parameters[i][1] + parameters[i][2] * q[i]);
        }
        for i in 0..n {
            let mut dq = 0.;
            for &(j, _) in &adj[i] {
                let dx = energy[j] - energy[i];
                let sign = if dx < 0. { 0. } else { 1. };
                // Keep this expression rather than selecting the denominator:
                // the pinned model's 0*NaN / unknown-parameter behavior matters.
                dq += dx / (sign * (ion[i] - ion[j]) + ion[j]);
            }
            let nh = mol.atoms[i].hydrogen as f64;
            if nh > 0. {
                let h = hq[i] / nh;
                let e = hp[0] + h * (hp[1] + hp[2] * h);
                let dx = e - energy[i];
                let sign = if dx < 0. { 0. } else { 1. };
                let dqh = dx / (sign * (ion[i] - 20.02) + 20.02);
                dq += nh * dqh;
                hq[i] -= nh * dqh * damp;
            }
            q[i] += damp * dq;
        }
        damp *= 0.5;
    }
    let parameterized = parameters.iter().all(|p| p.iter().any(|&v| v != 0.));
    (q, hq, parameterized)
}

fn py_min(a: f64, b: f64) -> f64 {
    if b < a {
        b
    } else {
        a
    }
}
fn py_max(a: f64, b: f64) -> f64 {
    if b > a {
        b
    } else {
        a
    }
}
pub(crate) struct Values {
    pub(crate) atoms: Vec<f64>,
    pub(crate) parameterized: bool,
    radicals: f64,
    min: f64,
    max: f64,
    vsa: Option<[f64; 14]>,
}
impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        let value = match id {
            "NumRadicalElectrons" => self.radicals,
            "MinPartialCharge" => self.min,
            "MaxPartialCharge" => self.max,
            "MinAbsPartialCharge" => py_min(self.min.abs(), self.max.abs()),
            "MaxAbsPartialCharge" => py_max(self.min.abs(), self.max.abs()),
            _ => {
                let i = id
                    .strip_prefix("PEOE_VSA")?
                    .parse::<usize>()
                    .ok()?
                    .checked_sub(1)?;
                *self.vsa.as_ref()?.get(i)?
            }
        };
        value.is_finite().then_some(value)
    }
}
pub(crate) fn calculate(
    mol: &Molecule,
    context: &Context,
    areas: Option<&crate::surface::Areas>,
) -> Values {
    let (q, _hq, parameterized) = charges(mol, context);
    let mut min = 500.;
    let mut max = -500.;
    for &charge in &q {
        min = py_min(charge, min);
        max = py_max(charge, max);
    }
    let vsa = areas.filter(|a| a.atoms.len() == q.len()).map(|areas| {
        let mut v = [0.; 14];
        for (&charge, &area) in q.iter().zip(&areas.atoms) {
            v[crate::surface::bin_of(&crate::charge_data::BINS, charge)] += area;
        }
        v
    });
    let radicals = context.states.iter().map(|s| s.radicals as f64).sum();
    Values {
        atoms: q,
        parameterized,
        radicals,
        min,
        max,
        vsa,
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::BondOrder;
    use std::collections::BTreeMap;
    fn reference_contexts() -> BTreeMap<&'static str, Molecule> {
        use crate::parser::{Atom, Bond};
        include_str!("../tests/fixtures/charge_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<_> = l.split('\t').collect();
                let atoms = f[1]
                    .split(';')
                    .map(|a| {
                        let a: Vec<_> = a.split(',').collect();
                        Atom {
                            symbol: a[0].into(),
                            hydrogen: a[1].parse().unwrap(),
                            charge: a[2].parse().unwrap(),
                            aromatic: a[3] == "1",
                            in_bracket: a[4] == "1",
                            ..Default::default()
                        }
                    })
                    .collect();
                let bonds: Vec<_> = f[2]
                    .split(';')
                    .filter(|s| !s.is_empty())
                    .map(|b| {
                        let b: Vec<_> = b.split(',').collect();
                        let order = match b[2] {
                            "1.0" => BondOrder::Single,
                            "2.0" => BondOrder::Double,
                            "3.0" => BondOrder::Triple,
                            "1.5" => BondOrder::Aromatic,
                            _ => panic!("unexpected bond"),
                        };
                        Bond::new(b[0].parse().unwrap(), b[1].parse().unwrap(), order)
                    })
                    .collect();
                (
                    f[0],
                    Molecule {
                        atoms,
                        bond_count: bonds.len() as i32,
                        bonds,
                    },
                )
            })
            .collect()
    }

    #[test]
    fn all_charge_scalars_match_independent_sanitized_contexts() {
        let contexts = reference_contexts();
        let mut previous = "";
        let mut actual = None;
        let mut errors = Vec::new();
        for l in include_str!("../tests/fixtures/descriptor_charge.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = l.split('\t').collect();
            let expected: f64 = f[2].parse().unwrap();
            if previous != f[0] {
                let mol = &contexts[f[0]];
                let c = Context::new(mol).unwrap();
                let areas = crate::surface::areas(mol);
                actual = Some(calculate(mol, &c, areas.as_ref()));
                previous = f[0];
            }
            let value = actual.as_ref().unwrap().value(f[1]);
            let tolerance = if f[1] == "NumRadicalElectrons" {
                0.
            } else {
                1e-9
            };
            let ok = if expected.is_finite() {
                value.is_some_and(|v| (v - expected).abs() <= tolerance)
            } else {
                value.is_none()
            };
            if !ok {
                errors.push(format!("{} {}: {value:?} vs {expected}", f[0], f[1]));
            }
        }
        assert!(
            errors.is_empty(),
            "{} scalar errors:\n{}",
            errors.len(),
            errors
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
    #[test]
    fn peoe_bins_keep_right_inclusion_and_nan_last() {
        let b = &crate::charge_data::BINS;
        for (i, &edge) in b.iter().enumerate() {
            assert_eq!(crate::surface::bin_of(b, edge - 1e-12), i);
            assert_eq!(crate::surface::bin_of(b, edge), i + 1);
            assert_eq!(crate::surface::bin_of(b, edge + 1e-12), i + 1);
        }
        assert_eq!(crate::surface::bin_of(b, f64::NAN), 13);
        assert_eq!(crate::surface::bin_of(b, f64::NEG_INFINITY), 0);
        assert_eq!(crate::surface::bin_of(b, f64::INFINITY), 13);
    }
    fn agrees(a: f64, b: f64) -> bool {
        if b.is_finite() {
            a.is_finite() && (a - b).abs() <= 1e-12
        } else {
            a.is_nan() && b.is_nan() || a == b
        }
    }
    #[test]
    fn atomic_gasteiger_and_attached_h_charges_match_independent_oracle() {
        let contexts = reference_contexts();
        let mut previous = "";
        let mut actual = (Vec::new(), Vec::new(), false);
        let mut errors = Vec::new();
        for line in include_str!("../tests/fixtures/charge_atomic.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let i: usize = f[1].parse().unwrap();
            if previous != f[0] {
                let m = &contexts[f[0]];
                actual = charges(m, &Context::new(m).unwrap());
                previous = f[0];
            }
            let q: f64 = f[2].parse().unwrap();
            let h: f64 = f[3].parse().unwrap();
            if !agrees(actual.0[i], q) || !agrees(actual.1[i], h) {
                errors.push(format!(
                    "{} atom{i}: q {} vs {q}; H {} vs {h}",
                    f[0], actual.0[i], actual.1[i]
                ));
            }
        }
        assert!(
            errors.is_empty(),
            "{} atomic errors:\n{}",
            errors.len(),
            errors
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
    #[test]
    fn full_hybridization_radical_and_conjugation_context_matches_oracle() {
        use crate::electronic::Hybridization as H;
        let contexts = reference_contexts();
        let mut errors = Vec::new();
        for line in include_str!("../tests/fixtures/charge_atomic.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let i: usize = f[1].parse().unwrap();
            let c = Context::new(&contexts[f[0]]).unwrap();
            let hyb = match c.hybridization[i] {
                H::S => "S",
                H::Sp => "SP",
                H::Sp2 => "SP2",
                H::Sp3 => "SP3",
                H::Sp2d => "SP2D",
                H::Sp3d => "SP3D",
                H::Sp3d2 => "SP3D2",
                H::Unspecified => "UNSPECIFIED",
            };
            if hyb != f[4] || c.states[i].radicals != f[5].parse::<i32>().unwrap() {
                errors.push(format!(
                    "{} atom{i}: {hyb}/{} vs {}/{}",
                    f[0], c.states[i].radicals, f[4], f[5]
                ));
            }
        }
        for line in include_str!("../tests/fixtures/charge_conjugation.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let i: usize = f[1].parse().unwrap();
            let c = Context::new(&contexts[f[0]]).unwrap();
            if c.conjugated[i] != (f[2] == "1") {
                errors.push(format!(
                    "{} bond{i}: conjugation {} vs {}",
                    f[0], c.conjugated[i], f[2]
                ));
            }
        }
        assert!(
            errors.is_empty(),
            "{} context errors:\n{}",
            errors.len(),
            errors
                .iter()
                .take(20)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
}
