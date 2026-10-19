//! Parse-once, catalog-ordered descriptor values. Null means undefined, never a placeholder zero.
pub use crate::descriptor_catalog::{DescriptorSpec, CATALOG_JSON, SPECS};
use crate::{
    parser::{BondOrder, Molecule, RingInfo},
    workbench::BasicDescriptors,
};

pub fn calculate(mol: &Molecule, base: &BasicDescriptors, mr: f64) -> Vec<Option<f64>> {
    let info = mol.ring_info();
    let crippen = crate::logp_crippen::calc_with_atom_contributions(mol);
    calculate_with_ring_info(mol, base, mr, &crippen.atoms, &info)
}

pub(crate) fn calculate_with_ring_info(
    mol: &Molecule,
    base: &BasicDescriptors,
    mr: f64,
    crippen: &[(f64, f64)],
    info: &RingInfo,
) -> Vec<Option<f64>> {
    let classes = crate::ring_classes_with_info(mol, info);
    let structural = crate::structural_descriptors::calculate(mol, info);
    let extras = crate::extra_descriptors::calculate(mol);
    let information = crate::information::calculate(mol);
    let morgan = crate::morgan_density::calculate(mol, info);
    let qed = crate::qed::calculate(mol, info, base);
    let fragment_counts = crate::fragment_descriptors::calculate(mol, info);
    let estate = crate::estate::indices(mol);
    let areas = crate::surface::areas_with_ring_info(mol, info);
    let surface = areas
        .as_ref()
        .and_then(|a| crate::surface::calculate(mol, crippen, estate.as_deref(), a));
    let electronic = crate::electronic::Context::with_ring_info(mol, info);
    let stereo = electronic
        .as_ref()
        .map(|c| crate::stereo::legacy(mol, info, c, &structural.bridgehead_atoms));
    let hydrogen_transfers = stereo
        .as_ref()
        .map_or(0, |s| s.implicit_hydrogen_transfers(mol, info));
    let sps = electronic.as_ref().zip(stereo.as_ref()).and_then(|(c, s)| {
        crate::stereo::spacial_score(mol, info, c, &structural.bridgehead_atoms, s)
    });
    let graph = electronic
        .as_ref()
        .and_then(|c| crate::graph_descriptors::calculate(mol, c));
    let charge = electronic
        .as_ref()
        .map(|c| crate::charge::calculate(mol, c, areas.as_ref()));
    let bcut = electronic
        .as_ref()
        .zip(charge.as_ref())
        .map(|(c, q)| crate::bcut::calculate(mol, c, q, crippen));
    let count = |symbol: &str| mol.atoms.iter().filter(|a| a.symbol == symbol).count() as f64;
    let hetero = mol
        .atoms
        .iter()
        .filter(|a| a.symbol != "H" && a.symbol != "C")
        .count() as f64;
    let heavy_mass = mol
        .atoms
        .iter()
        .filter(|a| a.symbol != "H")
        .try_fold(0.0, |sum, a| {
            let mass = match a.isotope.filter(|i| *i != 0) {
                Some(i) => crate::weights::isotope_mass(&a.symbol, i),
                None => crate::weights::atomic_weight(&a.symbol),
            }?;
            Some(sum + mass)
        });
    let valence_electrons = mol.atoms.iter().try_fold(0i64, |sum, a| {
        let z = crate::periodic_data::ELEMENTS
            .iter()
            .position(|e| e.0 == a.symbol)?;
        Some(
            sum + crate::descriptor_catalog::OUTER_ELECTRONS[z] as i64 - a.charge as i64
                + a.hydrogen as i64,
        )
    });
    // Lipinski NHOHCount counts attached hydrogen atoms, unlike NumHDonors' sites.
    let nhoh: i32 = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|(_, a)| matches!(a.symbol.as_str(), "N" | "O"))
        .map(|(i, a)| {
            a.hydrogen
                + mol
                    .neighbors(i)
                    .iter()
                    .filter(|(j, _)| mol.atoms[*j].symbol == "H")
                    .count() as i32
        })
        .sum();
    SPECS
        .iter()
        .map(|s| {
            let value = match s.id {
                "qed" => qed,
                "FpDensityMorgan1" => morgan?[0],
                "FpDensityMorgan2" => morgan?[1],
                "FpDensityMorgan3" => morgan?[2],
                "MolWt" => base.mw,
                "ExactMolWt" => base.exact_mass,
                "MolLogP" => base.logp,
                "TPSA" => base.tpsa,
                "NumHDonors" => base.hbd as f64,
                "NumHAcceptors" => base.hba as f64,
                "HeavyAtomCount" => base.heavy_atoms as f64,
                "ds_BondCount" => base.bonds as f64,
                "ds_FormalCharge" => base.formal_charge as f64,
                "ds_FragmentCount" => base.fragments as f64,
                "NumRotatableBonds" => base.rotatable_bonds as f64,
                "FractionCSP3" => base.fraction_csp3,
                "RingCount" => classes.len() as f64,
                "NumAromaticRings" => classes.iter().filter(|r| r.aromatic).count() as f64,
                "NumAliphaticRings" => classes.iter().filter(|r| !r.aromatic).count() as f64,
                "NumSaturatedRings" => classes.iter().filter(|r| r.saturated).count() as f64,
                "NumAromaticHeterocycles" => classes
                    .iter()
                    .filter(|r| r.aromatic && r.has_hetero)
                    .count() as f64,
                "NumAromaticCarbocycles" => classes
                    .iter()
                    .filter(|r| r.aromatic && !r.has_hetero)
                    .count() as f64,
                "NumAliphaticHeterocycles" => classes
                    .iter()
                    .filter(|r| !r.aromatic && r.has_hetero)
                    .count() as f64,
                "NumAliphaticCarbocycles" => classes
                    .iter()
                    .filter(|r| !r.aromatic && !r.has_hetero)
                    .count() as f64,
                "NumSaturatedHeterocycles" => classes
                    .iter()
                    .filter(|r| r.saturated && r.has_hetero)
                    .count() as f64,
                "NumSaturatedCarbocycles" => classes
                    .iter()
                    .filter(|r| r.saturated && !r.has_hetero)
                    .count() as f64,
                "NumHeterocycles" => classes.iter().filter(|r| r.has_hetero).count() as f64,
                "HeavyAtomMolWt" => return heavy_mass.filter(|v| v.is_finite()),
                "NumHeteroatoms" => hetero,
                "NOCount" => count("N") + count("O"),
                "NHOHCount" => nhoh as f64,
                "NumValenceElectrons" => return valence_electrons.map(|n| n as f64),
                "ds_ExplicitHydrogens" => {
                    (crate::num_explicit_hydrogens(mol) - hydrogen_transfers) as f64
                }
                "ds_ImplicitHydrogens" => {
                    (crate::num_implicit_hydrogens(mol) + hydrogen_transfers) as f64
                }
                "ds_TotalHydrogens" => crate::num_total_hydrogens(mol) as f64,
                "ds_CountC" => count("C"),
                "ds_CountN" => count("N"),
                "ds_CountO" => count("O"),
                "ds_CountS" => count("S"),
                "ds_CountP" => count("P"),
                "ds_CountF" => count("F"),
                "ds_CountCl" => count("Cl"),
                "ds_CountBr" => count("Br"),
                "ds_CountI" => count("I"),
                "ds_HalogenCount" => count("F") + count("Cl") + count("Br") + count("I"),
                "ds_SingleBonds" => crate::num_bonds_with_order(mol, BondOrder::Single) as f64,
                "ds_DoubleBonds" => crate::num_bonds_with_order(mol, BondOrder::Double) as f64,
                "ds_TripleBonds" => crate::num_bonds_with_order(mol, BondOrder::Triple) as f64,
                "ds_AromaticBonds" => crate::num_bonds_with_order(mol, BondOrder::Aromatic) as f64,
                "ds_RingAtoms" => info.atom_in_ring.iter().filter(|v| **v).count() as f64,
                "ds_RingBonds" => info.bond_in_ring.iter().filter(|v| **v).count() as f64,
                "ds_LargestRingSize" => info.rings.iter().map(Vec::len).max().unwrap_or(0) as f64,
                "ds_AromaticAtoms" => mol.atoms.iter().filter(|a| a.aromatic).count() as f64,
                "MolMR" => mr,
                "SPS" => return sps,
                "Ipc" => return information.as_ref().map(|value| value.weighted),
                "AvgIpc" => return information.as_ref().map(|value| value.average),
                "MaxEStateIndex" => {
                    return estate
                        .as_ref()
                        .and_then(|v| v.iter().copied().reduce(f64::max));
                }
                "MinEStateIndex" => {
                    return estate
                        .as_ref()
                        .and_then(|v| v.iter().copied().reduce(f64::min));
                }
                "MaxAbsEStateIndex" => {
                    return estate
                        .as_ref()
                        .and_then(|v| v.iter().map(|v| v.abs()).reduce(f64::max));
                }
                "MinAbsEStateIndex" => {
                    return estate
                        .as_ref()
                        .and_then(|v| v.iter().map(|v| v.abs()).reduce(f64::min));
                }
                _ => {
                    return crate::fragment_descriptors::index(s.id)
                        .and_then(|i| fragment_counts[i])
                        .or_else(|| surface.as_ref().and_then(|v| v.value(s.id)))
                        .or_else(|| graph.as_ref().and_then(|v| v.value(s.id)))
                        .or_else(|| charge.as_ref().and_then(|v| v.value(s.id)))
                        .or_else(|| bcut.as_ref().and_then(|v| v.value(s.id)))
                        .or_else(|| structural.value(s.id))
                        .or_else(|| stereo.as_ref().and_then(|v| v.value(s.id)))
                        .or_else(|| extras.value(s.id));
                }
            };
            value.is_finite().then_some(value)
        })
        .collect()
}

/// # Safety
/// A non-null output must address `capacity` writable, disjoint bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_workbench_catalog(output: *mut u8, capacity: usize) -> i32 {
    if output.is_null() && capacity != 0 {
        return -1;
    }
    if !output.is_null() && capacity >= CATALOG_JSON.len() {
        // SAFETY: caller owns writable capacity bytes; static JSON cannot overlap them.
        unsafe {
            std::ptr::copy_nonoverlapping(CATALOG_JSON.as_ptr(), output, CATALOG_JSON.len());
        }
    }
    CATALOG_JSON.len() as i32
}
