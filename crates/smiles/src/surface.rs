//! Labute approximate surface and VSA aggregation; no coordinates or 3D model.
//! Pinned RDKit MolSurf.cpp / EState_VSA.py; see RDKIT-DATA-LICENSE.txt.
use crate::parser::{BondOrder, Molecule, RingInfo};
use crate::surface_data::{AREA_BINS, ESTATE_BINS, LOGP_BINS, MR_BINS, RADII};

pub(crate) struct Areas {
    pub atoms: Vec<f64>,
    pub hydrogen: f64,
}

#[cfg(test)]
pub(crate) fn areas(mol: &Molecule) -> Option<Areas> {
    areas_with_ring_info(mol, &mol.ring_info())
}

pub(crate) fn areas_with_ring_info(mol: &Molecule, rings: &RingInfo) -> Option<Areas> {
    let aromatic = mol.aromatic_bond_flags(rings);
    let radii: Vec<f64> = mol
        .atoms
        .iter()
        .map(|atom| {
            crate::periodic_data::ELEMENTS
                .iter()
                .position(|e| e.0 == atom.symbol)
                .map(|z| RADII[z])
        })
        .collect::<Option<_>>()?;
    let mut atoms = vec![0.0; radii.len()];
    for (index, bond) in mol.bonds.iter().enumerate() {
        let (ri, rj) = (radii[bond.a], radii[bond.b]);
        let correction = if aromatic[index] {
            0.1
        } else {
            match bond.order {
                BondOrder::Aromatic => 0.0,
                BondOrder::Single | BondOrder::Dative => 0.0,
                BondOrder::Double => 0.2,
                BondOrder::Triple => 0.3,
            }
        };
        let bij = ri + rj - correction;
        let dij = (ri - rj).abs().max(bij).min(ri + rj);
        atoms[bond.a] += rj * rj - (ri - dij) * (ri - dij) / dij;
        atoms[bond.b] += ri * ri - (rj - dij) * (rj - dij) / dij;
    }
    // The registered default has one aggregate H correction at every graph
    // atom, independent of its attached-H count. This is the pinned Labute
    // approximation, not explicit-H expansion or a 3D solvent model.
    let rh = RADII[0];
    let mut hydrogen = 0.0;
    for (i, &ri) in radii.iter().enumerate() {
        let dij = (ri - rh).abs().max(ri + rh).min(ri + rh);
        atoms[i] += rh * rh - (ri - dij) * (ri - dij) / dij;
        hydrogen += ri * ri - (rh - dij) * (rh - dij) / dij;
    }
    for (value, &ri) in atoms.iter_mut().zip(&radii) {
        *value = std::f64::consts::PI * ri * (4.0 * ri - *value);
    }
    if hydrogen.abs() > 1e-4 {
        hydrogen = std::f64::consts::PI * rh * (4.0 * rh - hydrogen);
    }
    Some(Areas { atoms, hydrogen })
}

pub(crate) fn bin_of(bounds: &[f64], value: f64) -> usize {
    // Match bisect_right / upper_bound, including the pinned NaN behavior.
    bounds.partition_point(|bound| !(value < *bound))
}

pub(crate) struct Values {
    labute: f64,
    mr: [f64; 10],
    logp: [f64; 12],
    estate_vsa: Option<[f64; 11]>,
    vsa_estate: Option<[f64; 10]>,
}

pub(crate) fn calculate(
    mol: &Molecule,
    crippen: &[(f64, f64)],
    estate: Option<&[f64]>,
    area: &Areas,
) -> Option<Values> {
    if crippen.len() != mol.atoms.len() {
        return None;
    }
    let mut mr = [0.0; 10];
    let mut logp = [0.0; 12];
    for (&(lp, refractivity), &weight) in crippen.iter().zip(&area.atoms) {
        mr[bin_of(&MR_BINS, refractivity)] += weight;
        logp[bin_of(&LOGP_BINS, lp)] += weight;
    }
    let electronic = estate.filter(|v| v.len() == mol.atoms.len()).map(|estate| {
        let mut estate_vsa = [0.0; 11];
        let mut vsa_estate = [0.0; 10];
        for (&index, &weight) in estate.iter().zip(&area.atoms) {
            estate_vsa[bin_of(&ESTATE_BINS, index)] += weight;
            vsa_estate[bin_of(&AREA_BINS, weight)] += index;
        }
        (estate_vsa, vsa_estate)
    });
    Some(Values {
        labute: area.atoms.iter().sum::<f64>() + area.hydrogen,
        mr,
        logp,
        estate_vsa: electronic.map(|v| v.0),
        vsa_estate: electronic.map(|v| v.1),
    })
}

impl Values {
    pub(crate) fn value(&self, id: &str) -> Option<f64> {
        fn entry(id: &str, prefix: &str, values: &[f64]) -> Option<f64> {
            let index = id
                .strip_prefix(prefix)?
                .parse::<usize>()
                .ok()?
                .checked_sub(1)?;
            values.get(index).copied()
        }
        let result = if id == "LabuteASA" {
            Some(self.labute)
        } else if id.starts_with("SMR_VSA") {
            entry(id, "SMR_VSA", &self.mr)
        } else if id.starts_with("SlogP_VSA") {
            entry(id, "SlogP_VSA", &self.logp)
        } else if id.starts_with("EState_VSA") {
            self.estate_vsa
                .as_ref()
                .and_then(|v| entry(id, "EState_VSA", v))
        } else if id.starts_with("VSA_EState") {
            self.vsa_estate
                .as_ref()
                .and_then(|v| entry(id, "VSA_EState", v))
        } else {
            None
        };
        result.filter(|v| v.is_finite())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn atomic_surface_and_crippen_contributions_match_independent_reference() {
        let mut previous = "";
        let mut area = None;
        let mut typed = Vec::new();
        let mut differences = Vec::new();
        for line in include_str!("../tests/fixtures/surface_atomic.tsv")
            .lines()
            .filter(|s| !s.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let (smiles, kind) = (fields[0], fields[1]);
            let index: usize = fields[2].parse().unwrap();
            let expected: f64 = fields[3].parse().unwrap();
            if previous != smiles {
                let mol = crate::parser::parse(smiles).expect(smiles);
                area = areas(&mol);
                typed = crate::logp_crippen::calc_with_atom_contributions(&mol).atoms;
                previous = smiles;
            }
            let actual = match kind {
                "hydrogen" => area.as_ref().map(|a| a.hydrogen),
                "area" => area.as_ref().and_then(|a| a.atoms.get(index).copied()),
                "logp" => typed.get(index).map(|v| v.0),
                "mr" => typed.get(index).map(|v| v.1),
                _ => panic!("unknown reference field"),
            };
            let agrees = actual.is_some_and(|v| {
                if expected.is_finite() {
                    (v - expected).abs() <= 1e-9
                } else {
                    v.is_nan()
                }
            });
            if !agrees {
                differences.push(format!(
                    "{smiles} {kind}[{index}]: {actual:?} vs {expected}"
                ));
            }
        }
        assert!(
            differences.is_empty(),
            "{} atomic differences; first examples:\n{}",
            differences.len(),
            differences
                .iter()
                .take(25)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    #[test]
    fn every_surface_bin_keeps_its_lower_boundary() {
        use crate::surface_data::{AREA_BINS, ESTATE_BINS, LOGP_BINS, MR_BINS};
        for bounds in [
            &AREA_BINS[..],
            &ESTATE_BINS[..],
            &LOGP_BINS[..],
            &MR_BINS[..],
        ] {
            for (i, &edge) in bounds.iter().enumerate() {
                assert_eq!(bin_of(bounds, edge), i + 1, "exact lower boundary {edge}");
                assert_eq!(bin_of(bounds, edge - 1e-9), i, "below {edge}");
                assert_eq!(bin_of(bounds, edge + 1e-9), i + 1, "above {edge}");
            }
            assert_eq!(bin_of(bounds, f64::NEG_INFINITY), 0);
            assert_eq!(bin_of(bounds, f64::INFINITY), bounds.len());
            // Python bisect_right / C++ upper_bound puts NaN at the end.
            assert_eq!(bin_of(bounds, f64::NAN), bounds.len());
        }
    }
}
