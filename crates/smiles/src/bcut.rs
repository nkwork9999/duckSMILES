//! Pinned BCUT2D: heavy-atom Burden matrices with four atomic properties.
use crate::{electronic::Context, parser::Molecule};
fn hydrogen_suppressed(mol: &Molecule, context: &Context) -> Molecule {
    use crate::parser::BondOrder;
    let adj = mol.adjacency();
    let mut remap = vec![usize::MAX; mol.atoms.len()];
    let mut atoms = Vec::new();
    for (i, atom) in mol.atoms.iter().enumerate() {
        if atom.symbol == "H" {
            continue;
        }
        remap[i] = atoms.len();
        let mut a = atom.clone();
        let removed = adj[i]
            .iter()
            .filter(|&&(j, _)| mol.atoms[j].symbol == "H")
            .count() as i32;
        if removed > 0 {
            let vals = crate::graph_data::VALENCES[context.z[i] - 1];
            let nondefault = vals.iter().skip(1).any(|&v| v == context.states[i].valence);
            if a.in_bracket || a.aromatic || a.chirality.is_some() || nondefault {
                a.hydrogen += removed;
            } else {
                let accum: f64 = adj[i]
                    .iter()
                    .filter(|&&(j, _)| mol.atoms[j].symbol != "H")
                    .map(|&(_, b)| mol.bonds[b].valence_contribution(i))
                    .sum();
                let aromatic = adj[i].iter().any(|&(j, b)| {
                    mol.atoms[j].symbol != "H" && mol.bonds[b].order == BondOrder::Aromatic
                });
                let valence =
                    crate::electronic::integral_valence(context.z[i], a.charge, aromatic, accum);
                let effective = (context.z[i] as i32 - a.charge).clamp(1, 118) as usize;
                let allowed = crate::graph_data::VALENCES[effective - 1];
                a.hydrogen = if aromatic {
                    (allowed[0] - valence).max(0)
                } else {
                    allowed
                        .iter()
                        .copied()
                        .find(|&v| v >= valence)
                        .map_or(0, |v| v - valence)
                };
            }
        }
        a.nbr_order.clear();
        atoms.push(a);
    }
    let bonds: Vec<_> = mol
        .bonds
        .iter()
        .filter(|b| remap[b.a] != usize::MAX && remap[b.b] != usize::MAX)
        .map(|b| {
            let mut b = *b;
            b.a = remap[b.a];
            b.b = remap[b.b];
            b
        })
        .collect();
    Molecule {
        atoms,
        bond_count: bonds.len() as i32,
        bonds,
    }
}

pub(crate) struct Values([Option<f64>; 8]);
impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        let i = match id {
            "BCUT2D_MWHI" => 0,
            "BCUT2D_MWLOW" => 1,
            "BCUT2D_CHGHI" => 2,
            "BCUT2D_CHGLO" => 3,
            "BCUT2D_LOGPHI" => 4,
            "BCUT2D_LOGPLOW" => 5,
            "BCUT2D_MRHI" => 6,
            "BCUT2D_MRLOW" => 7,
            _ => return None,
        };
        self.0[i]
    }
}
fn from_properties(
    mol: &Molecule,
    charge: &crate::charge::Values,
    crippen: &[(f64, f64)],
) -> Values {
    use crate::parser::BondOrder;
    let n = mol.atoms.len();
    if n == 0 {
        return Values([Some(0.); 8]);
    }
    // The registered wrapper catches the throwing BCUT vector function, so an
    // unknown Gasteiger mode makes all eight fields undefined, including MW.
    if !charge.parameterized || charge.atoms.len() != n || crippen.len() != n {
        return Values([None; 8]);
    }
    let masses: Option<Vec<f64>> = mol
        .atoms
        .iter()
        .map(|a| match a.isotope.filter(|&i| i != 0) {
            Some(i) => crate::weights::isotope_mass(&a.symbol, i),
            None => crate::weights::atomic_weight(&a.symbol),
        })
        .collect();
    let Some(masses) = masses else {
        return Values([None; 8]);
    };
    let mut matrix = vec![0.001; n * n];
    for b in &mol.bonds {
        let v = match b.order {
            BondOrder::Dative => return Values([None; 8]),
            BondOrder::Single => 1.,
            BondOrder::Double => 0.7071067811865475,
            BondOrder::Triple => 0.5773502691896258,
            BondOrder::Aromatic => 0.8164965809277261,
        };
        matrix[b.a * n + b.b] = v;
        matrix[b.b * n + b.a] = v;
    }
    let mut out = [None; 8];
    for k in 0..4 {
        for i in 0..n {
            matrix[i * n + i] = match k {
                0 => masses[i],
                1 => charge.atoms[i],
                2 => crippen[i].0,
                _ => crippen[i].1,
            };
        }
        if let Some((lo, hi)) = crate::spectral::extrema(&matrix, n) {
            out[k * 2] = Some(hi);
            out[k * 2 + 1] = Some(lo);
        }
    }
    Values(out)
}
pub(crate) fn calculate(
    mol: &Molecule,
    context: &Context,
    charge: &crate::charge::Values,
    crippen: &[(f64, f64)],
) -> Values {
    if !mol.atoms.iter().any(|a| a.symbol == "H") {
        return from_properties(mol, charge, crippen);
    }
    let heavy = hydrogen_suppressed(mol, context);
    if heavy.atoms.is_empty() {
        return Values([Some(0.); 8]);
    }
    let Some(c) = Context::new(&heavy) else {
        return Values([None; 8]);
    };
    let q = crate::charge::calculate(&heavy, &c, None);
    let crippen = crate::logp_crippen::assign_contributions(&heavy);
    from_properties(&heavy, &q, &crippen)
}
#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::{Atom, Bond, BondOrder};
    use std::collections::BTreeMap;
    fn contexts(data: &'static str) -> BTreeMap<&'static str, Molecule> {
        data.lines()
            .filter(|l| !l.starts_with('#'))
            .map(|l| {
                let f: Vec<_> = l.split('\t').collect();
                let atoms = f[1]
                    .split(';')
                    .filter(|a| !a.is_empty())
                    .map(|a| {
                        let v: Vec<_> = a.split(',').collect();
                        Atom {
                            symbol: v[0].into(),
                            hydrogen: v[1].parse().unwrap(),
                            charge: v[2].parse().unwrap(),
                            aromatic: v[3] == "1",
                            in_bracket: v[4] == "1",
                            isotope: if v[5] == "0" {
                                None
                            } else {
                                Some(v[5].parse().unwrap())
                            },
                            ..Default::default()
                        }
                    })
                    .collect();
                let bonds: Vec<_> = f[2]
                    .split(';')
                    .filter(|b| !b.is_empty())
                    .map(|b| {
                        let v: Vec<_> = b.split(',').collect();
                        let o = match v[2] {
                            "1.0" => BondOrder::Single,
                            "2.0" => BondOrder::Double,
                            "3.0" => BondOrder::Triple,
                            "1.5" => BondOrder::Aromatic,
                            _ => panic!("unexpected bond"),
                        };
                        Bond::new(v[0].parse().unwrap(), v[1].parse().unwrap(), o)
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
    fn registered_bcut_values_match_independent_sanitized_contexts() {
        let original = contexts(include_str!("../tests/fixtures/bcut_original_contexts.tsv"));
        let mut previous = "";
        let mut result = None;
        let mut errors = Vec::new();
        for l in include_str!("../tests/fixtures/descriptor_bcut.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = l.split('\t').collect();
            let expected: f64 = f[2].parse().unwrap();
            if previous != f[0] {
                let m = &original[f[0]];
                let c = Context::new(m).unwrap();
                let q = crate::charge::calculate(m, &c, None);
                let cr = crate::logp_crippen::assign_contributions(m);
                result = Some(calculate(m, &c, &q, &cr));
                previous = f[0];
            }
            let v = result.as_ref().unwrap().value(f[1]);
            if !(if expected.is_finite() {
                v.is_some_and(|v| (v - expected).abs() <= 1e-9)
            } else {
                v.is_none()
            }) {
                errors.push(format!("{} {}: {v:?} vs {expected}", f[0], f[1]));
            }
        }
        assert!(
            errors.is_empty(),
            "{} errors:\n{}",
            errors.len(),
            errors
                .iter()
                .take(25)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }
    #[test]
    fn remove_all_hydrogens_matches_independent_graph_and_attached_counts() {
        let original = contexts(include_str!("../tests/fixtures/bcut_original_contexts.tsv"));
        let reference = contexts(include_str!("../tests/fixtures/bcut_contexts.tsv"));
        let mut errors = Vec::new();
        for (s, m) in original {
            let actual = hydrogen_suppressed(&m, &Context::new(&m).unwrap());
            let expected = &reference[s];
            if actual.atoms.len() != expected.atoms.len()
                || actual.bonds.len() != expected.bonds.len()
            {
                errors.push(format!("{s}: graph shape differs"));
                continue;
            }
            for (i, (a, b)) in actual.atoms.iter().zip(&expected.atoms).enumerate() {
                if (
                    a.symbol.as_str(),
                    a.hydrogen,
                    a.charge,
                    a.aromatic,
                    a.in_bracket,
                    a.isotope,
                ) != (
                    b.symbol.as_str(),
                    b.hydrogen,
                    b.charge,
                    b.aromatic,
                    b.in_bracket,
                    b.isotope,
                ) {
                    errors.push(format!(
                        "{s} atom{i}: H{}/{} flag{}/{}",
                        a.hydrogen, b.hydrogen, a.in_bracket, b.in_bracket
                    ));
                }
            }
            for (a, b) in actual.bonds.iter().zip(&expected.bonds) {
                if (a.a, a.b, a.order) != (b.a, b.b, b.order) {
                    errors.push(format!("{s}: bond differs"));
                }
            }
        }
        assert!(
            errors.is_empty(),
            "{} errors:\n{}",
            errors.len(),
            errors.join("\n")
        );
    }
}
