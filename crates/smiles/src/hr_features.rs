//! Bounded MS-FINDER HR grammar adapted to structural fingerprints (DM005).
//! Source: Tsugawa et al. 2016, doi:10.1021/acs.analchem.6b00770, Fig.1.
//! Single-bond CNOPS boundaries only; no empirical frequency/BDE weights.
//! This is not a complete fragmentation simulator or an observed ion assignment.
use crate::ms_features::{self, Fragment, ELEMENTS};
use crate::parser::{parse, BondOrder, Molecule};
use std::collections::{BTreeMap, BTreeSet};

pub const SPEC: &str = "duck-hr10-v1";
pub const NAMES: [&str; 10] = [
    "hr_positive_ion_mass",
    "hr_positive_ion_formula",
    "hr_positive_loss_formula",
    "hr_positive_cut_environment",
    "hr_positive_ion_path",
    "hr_negative_ion_mass",
    "hr_negative_ion_formula",
    "hr_negative_loss_formula",
    "hr_negative_cut_environment",
    "hr_negative_ion_path",
];
// NIST 2022 CODATA electron mass in u; neutral-H minus proton is not identical.
pub const ELECTRON_MASS: f64 = 0.0005485799090441;
type SparseCounts = BTreeMap<usize, u32>;

#[derive(Clone, Debug)]
pub struct Ion {
    pub fragment: usize,
    pub polarity: i8,
    pub hydrogen_shift: i8,
    pub formula: [u16; 13],
    pub mz: f64,
}

#[derive(Clone, Debug)]
pub struct Features {
    pub channels: Vec<Vec<(usize, u32)>>,
    pub ions: Vec<Ion>,
    pub fragments: Vec<Fragment>,
    pub eligible_fragments: usize,
    pub unsupported_boundary_fragments: usize,
}

/// Changes from the uncapped fragment retaining the parent's original H atoms.
pub fn first_shifts(element: &str, polarity: i8) -> Result<&'static [i8], String> {
    match (polarity, element) {
        (1, "C") => Ok(&[0]),
        (1, "N" | "O") => Ok(&[2]),
        (1, "P" | "S") => Ok(&[0, 2]),
        (-1, "C" | "P") => Ok(&[-2, 0]),
        (-1, "N" | "O") => Ok(&[0]),
        (-1, "S") => Ok(&[-1, 0]),
        (1 | -1, _) => Err(format!("HR has no CNOPS rule for {element}")),
        _ => Err("polarity must be +1 or -1".into()),
    }
}

fn boundary(m: &Molecule, f: &Fragment) -> Vec<(usize, usize)> {
    m.bonds
        .iter()
        .enumerate()
        .filter_map(|(i, b)| {
            let a = f.atoms.binary_search(&b.a).is_ok();
            let z = f.atoms.binary_search(&b.b).is_ok();
            (a != z).then_some((i, if a { b.a } else { b.b }))
        })
        .collect()
}

fn shifts(elements: &[&str], polarity: i8) -> Result<BTreeSet<i8>, String> {
    if elements.is_empty() || elements.len() > 2 {
        return Err("HR version1 requires one or two boundary bonds".into());
    }
    let mut result = BTreeSet::new();
    for first in 0..elements.len() {
        for &delta in first_shifts(elements[first], polarity)? {
            if elements.len() == 1 {
                result.insert(delta);
            } else {
                // P3/P4 or N4/N5: every subsequent single-bond boundary ±1H.
                first_shifts(elements[1 - first], polarity)?;
                result.insert(delta - 1);
                result.insert(delta + 1);
            }
        }
    }
    Ok(result)
}

fn ion_formula(f: &Fragment, delta: i8, parent: &[u16; 13], polarity: i8) -> Option<[u16; 13]> {
    let hydrogens = i32::from(f.formula[0]) + i32::from(delta);
    let available = i32::from(parent[0]) + i32::from(polarity);
    if hydrogens < 0 || hydrogens > available {
        return None;
    }
    let mut formula = f.formula;
    formula[0] = u16::try_from(hydrogens).ok()?;
    Some(formula)
}

fn ion_mass(formula: &[u16; 13], polarity: i8) -> f64 {
    formula
        .iter()
        .zip(ELEMENTS)
        .map(|(&count, element)| {
            f64::from(count)
                * crate::weights::monoisotopic_mass(element).expect("validated element")
        })
        .sum::<f64>()
        - f64::from(polarity) * ELECTRON_MASS
}

fn emit(counts: &mut SparseCounts, width: usize, words: &[u64]) {
    let bin = (ms_features::hash(words) % width as u64) as usize;
    *counts.entry(bin).or_default() += 1;
}

pub fn calculate(smiles: &str, width: usize) -> Result<Features, String> {
    let m = parse(smiles).ok_or("parse failed")?;
    if m.atoms.iter().any(|atom| atom.charge != 0) {
        return Err("HR version1 requires zero formal charge on every parent atom".into());
    }
    let original = ms_features::calculate(smiles, width)?;
    let mut channels = vec![SparseCounts::new(); 10];
    let mut ions = Vec::new();
    let boundaries: Vec<_> = original.fragments.iter().map(|f| boundary(&m, f)).collect();
    let eligible: Vec<_> = boundaries
        .iter()
        .map(|b| {
            !b.is_empty()
                && b.len() <= 2
                && b.iter().all(|&(bond, atom)| {
                    matches!(m.bonds[bond].order, BondOrder::Single)
                        && first_shifts(&m.atoms[atom].symbol, 1).is_ok()
                })
        })
        .collect();
    for polarity in [1_i8, -1_i8] {
        let offset = if polarity == 1 { 0 } else { 5 };
        let mut by_fragment = vec![Vec::new(); original.fragments.len()];
        for (index, f) in original.fragments.iter().enumerate() {
            if !eligible[index] {
                continue;
            }
            let elements: Vec<_> = boundaries[index]
                .iter()
                .map(|&(_, atom)| m.atoms[atom].symbol.as_str())
                .collect();
            for delta in shifts(&elements, polarity)? {
                let Some(formula) = ion_formula(f, delta, &original.formula, polarity) else {
                    continue;
                };
                let mz = ion_mass(&formula, polarity);
                if !mz.is_finite() || mz <= 0.0 {
                    continue;
                }
                let ion = Ion {
                    fragment: index,
                    polarity,
                    hydrogen_shift: delta,
                    formula,
                    mz,
                };
                emit(&mut channels[offset], width, &[1, mz.round() as u64]);
                let mut token = vec![2];
                token.extend(formula.iter().map(|&n| u64::from(n)));
                emit(&mut channels[offset + 1], width, &token);
                let mut loss = vec![3];
                for (element, &amount) in formula.iter().enumerate() {
                    let parent = i32::from(original.formula[element])
                        + if element == 0 { i32::from(polarity) } else { 0 };
                    let difference = parent - i32::from(amount);
                    if difference < 0 {
                        return Err("invalid parent/ion elemental conservation".into());
                    }
                    loss.push(difference as u64);
                }
                emit(&mut channels[offset + 2], width, &loss);
                emit(
                    &mut channels[offset + 3],
                    width,
                    &[
                        4,
                        mz.round() as u64,
                        f.environment,
                        (delta + 8) as u64,
                        u64::from(f.cuts),
                    ],
                );
                by_fragment[index].push(ions.len());
                ions.push(ion);
            }
        }
        for (parent_index, parent) in original.fragments.iter().enumerate() {
            if parent.cuts != 1 || !eligible[parent_index] {
                continue;
            }
            let (first_bond, charged_endpoint) = boundaries[parent_index][0];
            for (child_index, child) in original.fragments.iter().enumerate() {
                if child.cuts != 2
                    || !eligible[child_index]
                    || child.atoms.binary_search(&charged_endpoint).is_err()
                    || !child
                        .atoms
                        .iter()
                        .all(|atom| parent.atoms.binary_search(atom).is_ok())
                    || !boundaries[child_index]
                        .iter()
                        .any(|&(b, _)| b == first_bond)
                {
                    continue;
                }
                let mut unique_paths = BTreeSet::new();
                for &pi in &by_fragment[parent_index] {
                    for &ci in &by_fragment[child_index] {
                        let a = &ions[pi];
                        let b = &ions[ci];
                        if (b.hydrogen_shift - a.hydrogen_shift).abs() != 1
                            || b.formula
                                .iter()
                                .zip(a.formula)
                                .any(|(&child, parent)| child > parent)
                            || a.mz <= b.mz
                        {
                            continue;
                        }
                        unique_paths.insert((a.hydrogen_shift, b.hydrogen_shift));
                    }
                }
                for (pdelta, cdelta) in unique_paths {
                    let pformula = ion_formula(parent, pdelta, &original.formula, polarity)
                        .ok_or("invalid path parent")?;
                    let cformula = ion_formula(child, cdelta, &original.formula, polarity)
                        .ok_or("invalid path child")?;
                    let pmass = ion_mass(&pformula, polarity);
                    let cmass = ion_mass(&cformula, polarity);
                    emit(
                        &mut channels[offset + 4],
                        width,
                        &[
                            5,
                            pmass.round() as u64,
                            cmass.round() as u64,
                            (pmass - cmass).round() as u64,
                        ],
                    );
                }
            }
        }
    }
    let eligible_fragments = eligible.iter().filter(|&&x| x).count();
    Ok(Features {
        channels: channels
            .into_iter()
            .map(|c| c.into_iter().collect())
            .collect(),
        ions,
        eligible_fragments,
        unsupported_boundary_fragments: original.fragments.len() - eligible_fragments,
        fragments: original.fragments,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn published_initial_rules_and_invalid_domain() {
        assert_eq!(first_shifts("C", 1).unwrap(), &[0]);
        assert_eq!(first_shifts("N", 1).unwrap(), &[2]);
        assert_eq!(first_shifts("O", 1).unwrap(), &[2]);
        assert_eq!(first_shifts("S", -1).unwrap(), &[-1, 0]);
        assert_eq!(first_shifts("P", -1).unwrap(), &[-2, 0]);
        assert!(first_shifts("Cl", 1).is_err());
        assert!(first_shifts("C", 0).is_err());
        assert!(calculate("[NH3+]CC(=O)[O-]", 2048).is_err());
        assert!(calculate("[13CH3]O", 2048).is_err());
        assert!(calculate("C.C", 2048).is_err());
        assert!(calculate("CO", 1).is_err());
    }

    #[test]
    fn two_cut_order_is_a_union_not_double_counted() {
        assert_eq!(shifts(&["C", "N"], 1).unwrap(), BTreeSet::from([-1, 1, 3]));
        assert_eq!(
            shifts(&["C", "N"], 1).unwrap(),
            shifts(&["N", "C"], 1).unwrap()
        );
        assert_eq!(shifts(&["C", "C"], 1).unwrap(), BTreeSet::from([-1, 1]));
        assert_eq!(
            shifts(&["O", "C"], -1).unwrap(),
            BTreeSet::from([-3, -1, 1])
        );
    }

    #[test]
    fn nitrogen_positive_and_sulfur_negative_ion_mass() {
        let methylamine = calculate("CN", 2048).unwrap();
        let ammonium = methylamine
            .ions
            .iter()
            .find(|ion| {
                ion.polarity == 1
                    && ion.formula[2] == 1
                    && ion.formula[1] == 0
                    && ion.formula[0] == 4
            })
            .unwrap();
        let expected = crate::weights::monoisotopic_mass("N").unwrap()
            + 4.0 * crate::weights::monoisotopic_mass("H").unwrap()
            - ELECTRON_MASS;
        assert!((ammonium.mz - expected).abs() < 1e-12);
        let methanethiol = calculate("CS", 2048).unwrap();
        assert!(methanethiol.ions.iter().any(|ion| ion.polarity == -1
            && ion.formula[4] == 1
            && ion.formula[1] == 0
            && ion.hydrogen_shift == -1));
    }

    #[test]
    fn elemental_conservation_and_unsupported_boundaries() {
        for smiles in [
            "CO",
            "CCN",
            "CC(=O)O",
            "C1CCCCC1",
            "COP(=O)(O)O",
            "CS(=O)(=O)O",
        ] {
            let x = calculate(smiles, 2048).unwrap();
            let parent = ms_features::calculate(smiles, 2048).unwrap().formula;
            for ion in &x.ions {
                for (element, &amount) in ion.formula.iter().enumerate() {
                    let cap = i32::from(parent[element])
                        + if element == 0 {
                            i32::from(ion.polarity)
                        } else {
                            0
                        };
                    assert!(i32::from(amount) <= cap);
                }
                assert!(ion.mz > 0.0);
            }
        }
        let aromatic = calculate("c1ccccc1", 2048).unwrap();
        assert_eq!(aromatic.eligible_fragments, 0);
        assert!(aromatic.unsupported_boundary_fragments > 0);
        assert!(aromatic.ions.is_empty());
        let carbonyl = calculate("CC=O", 2048).unwrap();
        assert!(carbonyl.unsupported_boundary_fragments > 0);
    }

    #[test]
    fn representation_invariance_and_nonempty_paths() {
        for (a, b) in [
            ("CCN", "NCC"),
            ("CC(=O)O", "OC(C)=O"),
            ("CC1CCCCC1", "C1CCC(C)CC1"),
        ] {
            assert_eq!(
                calculate(a, 2048).unwrap().channels,
                calculate(b, 2048).unwrap().channels
            );
        }
        let x = calculate("CCNCCO", 2048).unwrap();
        assert!(!x.channels[4].is_empty());
        assert!(!x.channels[9].is_empty());
    }
}
