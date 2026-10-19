//! Pinned chemical preparation for workbench calculations only.
//! Coordination includes a bounded canonical subset, not the full sanitizer.
//! Source text, stereo annotations and the raw graph are retained.
use crate::parser::{BondOrder, Molecule};

/// Prepared-only strict explicit-valence check, matching pinned Atom.cpp.
/// Read-only: reject invalid inputs without changing a valid graph's values.
/// This is not Atom.hasValenceViolation, which applies additional charge rules.
pub(crate) fn check_explicit_valence(mol: &Molecule) -> Option<()> {
    let numbers = crate::electronic::atomic_numbers(mol)?;
    let mut sums: Vec<_> = mol
        .atoms
        .iter()
        .map(|a| {
            if a.in_bracket && !a.prepared_h_is_implicit {
                f64::from(a.hydrogen)
            } else {
                0.
            }
        })
        .collect();
    let mut aromatic: Vec<_> = mol.atoms.iter().map(|a| a.aromatic).collect();
    for bond in &mol.bonds {
        for i in [bond.a, bond.b] {
            sums[i] += bond.valence_contribution(i);
            aromatic[i] |= bond.order == BondOrder::Aromatic;
        }
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        let z = numbers[i];
        let original = crate::graph_data::VALENCES[z - 1];
        // The strict property-cache path clamps the effective atomic number;
        // checkIt=true in hasValenceViolation has a different contract.
        let effective = if original == [-1] {
            z
        } else {
            (z as i64 - i64::from(a.charge)).clamp(0, 118) as usize
        };
        let allowed = if effective == 0 {
            &[-1][..]
        } else {
            crate::graph_data::VALENCES[effective - 1]
        };
        let default = allowed[0];
        let mut accum = sums[i];
        if aromatic[i] && default >= 0 && accum > f64::from(default) {
            let mut lower = default;
            for &v in allowed {
                if v < 0 || f64::from(v) > accum {
                    break;
                }
                lower = v;
            }
            if accum - f64::from(lower) <= 1.5 {
                accum = f64::from(lower);
            }
        }
        let explicit = (accum + 0.1).round() as i64;
        let mut maximum = *allowed.last()?;
        let mut offset = 0i64;
        if (effective > 16 && matches!(z, 15 | 16)) || (effective > 34 && matches!(z, 33 | 34)) {
            maximum = *original.last()?;
            offset = -i64::from(a.charge);
        }
        if z == 1 && a.charge == -1 {
            maximum = 2;
        }
        // -1 is the reference's unbounded high-valence sentinel.
        if maximum >= 0 && *original.last()? >= 0 && explicit + offset > i64::from(maximum) {
            return None;
        }
    }
    Some(())
}

/// Fill implicit H at atoms incident to aromatic bonds using numeric1.5 bond
/// valence and Atom.cpp rounding. Bracket H and raw-parser graphs are retained.
/// An aromatic bond need not imply an aromatic atom or an aromatic ring.
pub(crate) fn aromatic_hydrogens(mol: &mut Molecule) {
    let mut sums = vec![0.; mol.atoms.len()];
    let mut touched = vec![false; mol.atoms.len()];
    for bond in &mol.bonds {
        for index in [bond.a, bond.b] {
            sums[index] += bond.valence_contribution(index);
            touched[index] |= bond.order == BondOrder::Aromatic;
        }
    }
    for (index, atom) in mol.atoms.iter_mut().enumerate() {
        if atom.in_bracket || !touched[index] {
            continue;
        }
        let Some(z) = crate::periodic_data::ELEMENTS
            .iter()
            .position(|e| e.0 == atom.symbol)
            .map(|z| z + 1)
        else {
            continue;
        };
        let original = crate::graph_data::VALENCES[z - 1];
        let effective = if original == [-1] {
            z
        } else {
            (z as i32 - atom.charge).clamp(1, 118) as usize
        };
        let default = crate::graph_data::VALENCES[effective - 1][0];
        let explicit = crate::electronic::integral_valence(z, atom.charge, true, sums[index]);
        atom.hydrogen = if default < 0 {
            0
        } else {
            (default - explicit).max(0)
        };
    }
}

/// Aromaticity.cpp default-valence and radical eligibility, scoped to prepared
/// workbench graphs. Raw parser perception retains its existing contract.
pub(crate) fn aromatic_donations(
    mol: &Molecule,
    rings: &crate::parser::RingInfo,
) -> Option<Vec<Option<u32>>> {
    let z = crate::electronic::atomic_numbers(mol)?;
    let states = crate::electronic::atom_states(mol, &z);
    let adj = mol.adjacency();
    let mut donors: Vec<Option<u32>> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            // Default aromaticity admits higher-row exceptions only Se/Te.
            if z[i] > 18 && !matches!(z[i], 34 | 52) {
                return None;
            }
            let default = crate::graph_data::VALENCES[z[i] - 1][0];
            let effective = z[i] as i32 - a.charge;
            if default > 0
                && (!(1..=118).contains(&effective)
                    || states[i].valence > crate::graph_data::VALENCES[effective as usize - 1][0])
            {
                return None;
            }
            if states[i].radicals > 0 && (z[i] != 6 || a.charge != 0) {
                return None;
            }
            // Aromaticity.cpp countAtomElec: sigma coordination, default
            // valence, lone-pair electrons and cached radicals. In particular
            // a saturated boron with zero spare electrons is NOT an eligible
            // empty-orbital donor merely because its element is B.
            let degree = adj[i]
                .iter()
                .filter(|&&(_, b)| mol.bonds[b].valence_contribution(i) > 0.)
                .count() as i32
                + a.hydrogen;
            if default <= 1 || degree > 3 {
                return None;
            }
            let outer = crate::descriptor_catalog::OUTER_ELECTRONS[z[i] - 1];
            let mut electrons =
                default - degree + (outer - default - a.charge).max(0) - states[i].radicals;
            let explicit_h = if a.in_bracket && !a.prepared_h_is_implicit {
                a.hydrogen
            } else {
                0
            };
            let explicit_valence = states[i].valence - a.hydrogen + explicit_h;
            if electrons > 1 && explicit_valence - adj[i].len() as i32 > 1 {
                electrons = 1;
            }
            if electrons < 0 {
                return None;
            }
            let multiple = adj[i]
                .iter()
                .filter(|&&(_, b)| {
                    matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Triple)
                })
                .count();
            if multiple > 1 {
                return None;
            }
            let exocyclic = adj[i]
                .iter()
                .find(|&&(_, b)| {
                    !rings.bond_in_ring[b] && mol.bonds[b].valence_contribution(i) >= 2.
                })
                .map(|&(j, _)| j);
            let cyclic_multiple = adj[i]
                .iter()
                .any(|&(_, b)| rings.bond_in_ring[b] && mol.bonds[b].valence_contribution(i) >= 2.);
            let steals = exocyclic.is_some_and(|j| {
                let other = crate::descriptor_catalog::OUTER_ELECTRONS[z[j] - 1];
                other > outer || (other == outer && z[j] < z[i])
            });
            match electrons {
                0 => {
                    if exocyclic.is_some() {
                        Some(0)
                    } else if cyclic_multiple {
                        Some(1)
                    } else {
                        None
                    }
                }
                1 => {
                    if exocyclic.is_some() {
                        Some(u32::from(!steals))
                    } else if states[i].valence - a.hydrogen != degree - a.hydrogen {
                        Some(1)
                    } else if a.charge == 1 {
                        Some(0)
                    } else {
                        None
                    }
                }
                _ => {
                    if steals {
                        electrons -= 1;
                    }
                    Some(if electrons % 2 == 1 { 1 } else { 2 })
                }
            }
        })
        .collect();
    // Pinned Aromaticity.cpp caches each donor at its first ring encounter.
    // Neutral two-coordinate O/S without multiple bonds are ether bridges in
    // macrocycles (>=9), rather than two-electron aromatic donors. Small-ring
    // O/S such as furan/thiophene retain their donor behavior.
    let mut seen = vec![false; mol.atoms.len()];
    for ring in &rings.rings {
        for &bi in ring {
            for i in [mol.bonds[bi].a, mol.bonds[bi].b] {
                if seen[i] {
                    continue;
                }
                seen[i] = true;
                let a = &mol.atoms[i];
                if ring.len() >= 9
                    && donors[i] == Some(2)
                    && matches!(a.symbol.as_str(), "O" | "S")
                    && a.charge == 0
                    && adj[i].len() == 2
                    && adj[i].iter().all(|&(_, b)| {
                        !matches!(mol.bonds[b].order, BondOrder::Double | BondOrder::Triple)
                    })
                {
                    donors[i] = None;
                }
            }
        }
    }
    Some(donors)
}

fn explicit_valence_twice(mol: &Molecule, atom: usize, adj: &[Vec<(usize, usize)>]) -> i64 {
    let a = &mol.atoms[atom];
    // Bracket H counts are explicit. Implicit H must not turn neutral N with
    // four explicit bonds into a five-valence cleanup candidate.
    let h = if a.in_bracket {
        i64::from(a.hydrogen.max(0)) * 2
    } else {
        0
    };
    h + adj[atom]
        .iter()
        .map(|&(_, b)| match mol.bonds[b].order {
            BondOrder::Single => 2,
            BondOrder::Double => 4,
            BondOrder::Triple => 6,
            BondOrder::Aromatic => 3,
            BondOrder::Dative => {
                if mol.bonds[b].b == atom {
                    2
                } else {
                    0
                }
            }
        })
        .sum::<i64>()
}

// QueryOps MDL M-query membership, not a periodic-table "transition metal" set.
fn is_metal(z: usize) -> bool {
    (1..=118).contains(&z)
        && !matches!(
            z,
            1 | 2
                | 5
                | 6
                | 7
                | 8
                | 9
                | 10
                | 14
                | 15
                | 16
                | 17
                | 18
                | 33
                | 34
                | 35
                | 36
                | 52
                | 53
                | 54
                | 85
                | 86
        )
}

/// Infer donor→metal bonds following pinned MolOps.cpp. Multiple metal choices
/// require validated canonical ordering; unsupported rank contexts are refused.
pub(crate) fn coordination(mol: &mut Molecule) -> Option<bool> {
    let z = crate::electronic::atomic_numbers(mol)?;
    let adj = mol.adjacency();
    let mut conversions = Vec::new();
    for i in 0..mol.atoms.len() {
        if is_metal(z[i]) || matches!(z[i], 1 | 2 | 9 | 10) {
            continue;
        }
        let effective = z[i] as i32 - mol.atoms[i].charge;
        if !(1..=118).contains(&effective) {
            continue;
        }
        let max = *crate::graph_data::VALENCES[effective as usize - 1].last()?;
        if max <= 0 {
            continue;
        }
        // Atom.cpp applies aromatic valence rounding to incident AROMATIC
        // bonds even when the source atom itself is uppercase/aliphatic.
        // Otherwise neutral N with two ':' edges plus one metal SINGLE is
        // falsely considered overvalent and converted to a dative bond.
        let aromatic_valence = mol.atoms[i].aromatic
            || adj[i]
                .iter()
                .any(|&(_, b)| mol.bonds[b].order == BondOrder::Aromatic);
        let ev = crate::electronic::integral_valence(
            z[i],
            mol.atoms[i].charge,
            aromatic_valence,
            explicit_valence_twice(mol, i, &adj) as f64 / 2.,
        );
        if ev < max
            || (ev == max
                && !(mol.atoms[i].aromatic && adj[i].len() as i32 + mol.atoms[i].hydrogen == 4))
        {
            continue;
        }
        let candidates: Vec<_> = adj[i]
            .iter()
            .copied()
            .filter(|&(j, b)| is_metal(z[j]) && mol.bonds[b].order == BondOrder::Single)
            .collect();
        if !candidates.is_empty() {
            conversions.push((i, candidates));
        }
    }
    let mut chosen = Vec::new();
    if conversions.iter().any(|(_, c)| c.len() > 1) {
        if let Some(ranks) = crate::workbench_canonical::ranks(mol) {
            conversions.sort_by_key(|(donor, _)| ranks[*donor]);
            let mut dative_counts = vec![0usize; mol.atoms.len()];
            for bond in &mol.bonds {
                if bond.order == BondOrder::Dative {
                    dative_counts[bond.a] += 1;
                    dative_counts[bond.b] += 1;
                }
            }
            for (donor, candidates) in &conversions {
                let &(metal, bi) = candidates.iter().min_by_key(|&&(metal, _)| {
                    (dative_counts[metal], std::cmp::Reverse(ranks[metal]))
                })?;
                chosen.push((*donor, metal, bi));
                dative_counts[*donor] += 1;
                dative_counts[metal] += 1;
            }
        } else if conversions.len() == 1 {
            // Stereo tags remain intact. If there is only one conversion and
            // one candidate's INITIAL canonical prefix is strictly greatest,
            // no neighbor or stereo refinement can change that comparison.
            // The bracket/single-bond scope also guarantees zero initial dative
            // counts. Ties still require the complete canonical algorithm.
            let keys = crate::workbench_canonical::initial_keys(mol)?;
            let (donor, candidates) = &conversions[0];
            let &(metal, bi) = candidates.iter().max_by_key(|&&(m, _)| keys[m])?;
            if candidates
                .iter()
                .filter(|&&(m, _)| keys[m] == keys[metal])
                .count()
                != 1
            {
                return None;
            }
            chosen.push((*donor, metal, bi));
        } else {
            return None;
        }
    } else {
        chosen.extend(
            conversions
                .iter()
                .map(|(donor, c)| (*donor, c[0].0, c[0].1)),
        );
    }
    for &(donor, metal, bi) in &chosen {
        let bond = &mut mol.bonds[bi];
        bond.a = donor;
        bond.b = metal;
        bond.order = BondOrder::Dative;
    }
    Some(!chosen.is_empty())
}

/// Pinned cleanupChirality clears an out-of-range geometry tag without rejecting
/// the molecule. Explicit permutation grammar is validated during preparation.
pub(crate) fn non_tetrahedral_chirality(mol: &mut Molecule) {
    if !mol
        .atoms
        .iter()
        .any(|a| a.non_tetrahedral_geometry().is_some())
    {
        return;
    }
    let adj = mol.adjacency();
    for (i, atom) in mol.atoms.iter_mut().enumerate() {
        if let Some(geometry) = atom.non_tetrahedral_geometry() {
            let total_degree = adj[i].len() as i32 + atom.hydrogen;
            if !(2..=geometry.maximum_degree()).contains(&total_degree) {
                atom.chirality = None;
            }
        }
    }
}

pub(crate) fn apply(mol: &mut Molecule) -> bool {
    let adj = mol.adjacency();
    let mut changed = false;
    let mut changed_atoms = vec![false; mol.atoms.len()];
    let nitrogens: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .filter(|&(i, a)| {
            a.symbol == "N" && a.charge == 0 && explicit_valence_twice(mol, i, &adj) == 10
        })
        .map(|(i, _)| i)
        .collect();

    for &i in &nitrogens {
        if let Some(&(j, b)) = adj[i].iter().find(|&&(j, b)| {
            mol.atoms[j].symbol == "O"
                && mol.atoms[j].charge == 0
                && mol.bonds[b].order == BondOrder::Double
        }) {
            mol.bonds[b].order = BondOrder::Single;
            mol.atoms[i].charge = 1;
            mol.atoms[j].charge = -1;
            changed_atoms[i] = true;
            changed_atoms[j] = true;
            changed = true;
        }
    }
    // Use the original candidate list, even if the first pass already changed
    // a center. This preserves the upstream O=N#N two-pass charge behavior.
    for &i in &nitrogens {
        if let Some(&(j, b)) = adj[i].iter().find(|&&(j, b)| {
            mol.atoms[j].symbol == "N"
                && mol.atoms[j].charge == 0
                && mol.bonds[b].order == BondOrder::Triple
        }) {
            mol.bonds[b].order = BondOrder::Double;
            mol.atoms[i].charge = 1;
            mol.atoms[j].charge = -1;
            changed_atoms[i] = true;
            changed_atoms[j] = true;
            changed = true;
        }
    }
    for i in 0..mol.atoms.len() {
        if mol.atoms[i].charge != 0 {
            continue;
        }
        let valence = explicit_valence_twice(mol, i, &adj);
        if mol.atoms[i].symbol == "P" && valence == 10 && adj[i].len() == 3 {
            let oxygen = adj[i]
                .iter()
                .filter(|&&(j, b)| {
                    mol.atoms[j].symbol == "O"
                        && mol.atoms[j].charge == 0
                        && mol.bonds[b].order == BondOrder::Double
                })
                .last()
                .copied();
            let carbon_or_nitrogen = adj[i].iter().any(|&(j, b)| {
                matches!(mol.atoms[j].symbol.as_str(), "C" | "N")
                    && adj[j].len() >= 2
                    && mol.bonds[b].order == BondOrder::Double
            });
            if let Some((j, b)) = oxygen.filter(|_| carbon_or_nitrogen) {
                mol.bonds[b].order = BondOrder::Single;
                mol.atoms[i].charge = 1;
                mol.atoms[j].charge = -1;
                changed_atoms[i] = true;
                changed_atoms[j] = true;
                changed = true;
            }
        } else if matches!(mol.atoms[i].symbol.as_str(), "Cl" | "Br" | "I")
            && matches!(valence, 6 | 10 | 14)
            && adj[i].iter().all(|&(j, _)| mol.atoms[j].symbol == "O")
        {
            let mut charge = 0;
            for &(j, b) in &adj[i] {
                if mol.bonds[b].order == BondOrder::Double {
                    mol.bonds[b].order = BondOrder::Single;
                    mol.atoms[j].charge = -1;
                    changed_atoms[i] = true;
                    changed_atoms[j] = true;
                    charge += 1;
                    changed = true;
                }
            }
            mol.atoms[i].charge = charge;
        }
    }
    for (i, touched) in changed_atoms.into_iter().enumerate() {
        if !touched || mol.atoms[i].in_bracket {
            continue;
        }
        // Property-cache recomputation follows cleanup upstream. A center
        // changed in both N passes can gain implicit H; retain explicit H.
        let Some(z) = crate::periodic_data::ELEMENTS
            .iter()
            .position(|e| e.0 == mol.atoms[i].symbol)
            .map(|z| z + 1)
        else {
            continue;
        };
        let effective = z as i32 - mol.atoms[i].charge;
        if !(1..=118).contains(&effective) {
            continue;
        }
        let valence = (explicit_valence_twice(mol, i, &adj) + 1) / 2;
        mol.atoms[i].hydrogen = crate::graph_data::VALENCES[effective as usize - 1]
            .iter()
            .copied()
            .find(|&v| i64::from(v) >= valence)
            .map_or(0, |v| (i64::from(v) - valence) as i32);
    }
    changed
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn prepared_ring_membership_and_first_sizes_match_independent_reference() {
        let mut membership_errors = 0;
        let mut first_size_errors = 0;
        let mut tied_order_diagnostics = 0;
        let mut details = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/prepared_rings.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let (source, expected) = line.split_once('\t').unwrap();
            let (mol, rings) = crate::workbench::prepare(source).expect(source);
            let actual: Vec<String> = rings
                .rings
                .iter()
                .map(|ring| {
                    let mut edges: Vec<_> = ring
                        .iter()
                        .map(|&bi| {
                            let b = &mol.bonds[bi];
                            format!("{}-{}", b.a.min(b.b), b.a.max(b.b))
                        })
                        .collect();
                    edges.sort();
                    edges.join(",")
                })
                .collect();
            let mut actual_sets = actual.clone();
            actual_sets.sort();
            let mut expected_sets: Vec<_> = expected
                .split(';')
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect();
            expected_sets.sort();
            if actual_sets != expected_sets {
                membership_errors += 1;
                if details.len() < 12 {
                    details.push(format!(
                        "{}: ring sets {} vs {}",
                        source.chars().take(100).collect::<String>(),
                        actual.len(),
                        expected_sets.len()
                    ));
                }
            }
            let reference_order: Vec<_> = expected
                .split(';')
                .filter(|s| !s.is_empty())
                .map(str::to_owned)
                .collect();
            let first_sizes = |keys: &[String]| {
                let mut sizes = vec![0; mol.atoms.len()];
                for key in keys {
                    let edges: Vec<_> = key.split(',').collect();
                    for edge in &edges {
                        let (a, b) = edge.split_once('-').unwrap();
                        for i in [a.parse::<usize>().unwrap(), b.parse::<usize>().unwrap()] {
                            if sizes[i] == 0 {
                                sizes[i] = edges.len();
                            }
                        }
                    }
                }
                sizes
            };
            if first_sizes(&actual) != first_sizes(&reference_order) {
                first_size_errors += 1;
                if details.len() < 12 {
                    details.push(format!("source {count}: first ring sizes differ"));
                }
            }
            // std::sort does not specify equal-size ties. These are diagnostic:
            // donor first-size vectors are checked here; all source stereo
            // masks/factors are independently checked by the SPS fixture.
            tied_order_diagnostics += usize::from(actual.join(";") != expected);
            count += 1;
        }
        assert_eq!(count, 658);
        eprintln!(
            "{count} ring references; {tied_order_diagnostics} internal tied-order differences (diagnostic only)"
        );
        assert_eq!(
            (membership_errors, first_size_errors),
            (0, 0),
            "first membership differences:\n{}",
            details.join("\n")
        );
    }

    #[test]
    fn workbench_aromatic_preparation_preserves_raw_graph_consumer_semantics() {
        let source = "C1:C:C:C:C:C:1";
        let (written, _) = crate::native_graph::graph_as_written_with_ring_info(source).unwrap();
        assert!(
            written.atoms.iter().all(|a| !a.aromatic),
            "uppercase atoms remain aliphatic as written"
        );
        let legacy = crate::native_graph::graph(source).unwrap();
        assert!(
            legacy.atoms.iter().all(|a| a.aromatic),
            "legacy perceived graph keeps its contract"
        );
        assert_eq!(crate::workbench::analyze(source).unwrap().formula, "C6H6");
        for (source, aromatic_bonds) in [
            ("[n+]1ccccc1", 6),
            ("CP1(C)=CC=CN=C1C", 6),
            ("c1cc[asH]c1", 5),
        ] {
            let raw = crate::native_graph::graph(source).unwrap();
            assert_eq!(
                raw.bonds
                    .iter()
                    .filter(|b| b.order == BondOrder::Aromatic)
                    .count(),
                aromatic_bonds
            );
            let row = crate::workbench::analyze(source).unwrap();
            let field = crate::descriptor_profile::SPECS
                .iter()
                .position(|s| s.id == "ds_AromaticBonds")
                .unwrap();
            assert_eq!(row.descriptor_values[field], Some(0.));
            let raw_after = crate::native_graph::graph(source).unwrap();
            assert_eq!(
                raw_after
                    .bonds
                    .iter()
                    .filter(|b| b.order == BondOrder::Aromatic)
                    .count(),
                aromatic_bonds
            );
        }
    }

    #[test]
    fn directed_preparation_matches_independent_atomic_bond_and_ring_contexts() {
        use crate::electronic::Hybridization as H;
        let mut errors = Vec::new();
        for line in include_str!("../tests/fixtures/coordination_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let (mut mol, mut rings) =
                crate::native_graph::graph_as_written_with_ring_info(f[0]).unwrap();
            apply(&mut mol);
            if coordination(&mut mol).expect(f[0]) {
                rings = mol.ring_info();
            }
            mol.perceive_aromaticity_with_ring_info(&rings);
            let context = crate::electronic::Context::new(&mol).unwrap();
            let atoms: Vec<_> = f[1].split(';').collect();
            assert_eq!(atoms.len(), mol.atoms.len());
            for (i, (a, expected)) in mol.atoms.iter().zip(atoms).enumerate() {
                let e: Vec<_> = expected.split(',').collect();
                let hybrid = match context.hybridization[i] {
                    H::Unspecified => 0,
                    H::S => 1,
                    H::Sp => 2,
                    H::Sp2 => 3,
                    H::Sp3 => 4,
                    H::Sp2d => 5,
                    H::Sp3d => 6,
                    H::Sp3d2 => 7,
                };
                let actual = format!(
                    "{},{},{},{},{},{},{},{}",
                    a.symbol,
                    a.hydrogen,
                    a.charge,
                    usize::from(a.aromatic),
                    a.isotope.unwrap_or(0),
                    context.states[i].valence,
                    hybrid,
                    context.states[i].radicals
                );
                // noImplicit and bracket syntax are representation choices;
                // compare actual H, total valence and electronic state.
                let wanted = format!(
                    "{},{},{},{},{},{},{},{}",
                    e[0], e[1], e[2], e[3], e[5], e[6], e[7], e[8]
                );
                if actual != wanted {
                    errors.push(format!("{} atom{i}: {actual} vs {wanted}", f[0]));
                }
            }
            // Parsers append ring-closure edges at different positions.
            // Compare undirected endpoint pairs, preserving donor→metal order
            // only for dative edges, and map rings to endpoint pairs too.
            let bonds: Vec<Vec<usize>> = f[2]
                .split(';')
                .filter(|s| !s.is_empty())
                .map(|s| s.split(',').map(|v| v.parse().unwrap()).collect())
                .collect();
            assert_eq!(bonds.len(), mol.bonds.len());
            let pair = |a: usize, b: usize| (a.min(b), a.max(b));
            let mut expected_bonds = Vec::new();
            for b in &bonds {
                let (a, z) = if b[2] == 17 {
                    (b[0], b[1])
                } else {
                    pair(b[0], b[1])
                };
                expected_bonds.push((a, z, b[2], b[3]));
            }
            let mut actual_bonds = Vec::new();
            for (i, b) in mol.bonds.iter().enumerate() {
                let kind = match b.order {
                    BondOrder::Single => 1,
                    BondOrder::Double => 2,
                    BondOrder::Triple => 3,
                    BondOrder::Aromatic => 12,
                    BondOrder::Dative => 17,
                };
                let (a, z) = if kind == 17 {
                    (b.a, b.b)
                } else {
                    pair(b.a, b.b)
                };
                actual_bonds.push((a, z, kind, usize::from(context.conjugated[i])));
            }
            expected_bonds.sort_unstable();
            actual_bonds.sort_unstable();
            if actual_bonds != expected_bonds {
                errors.push(format!(
                    "{} bonds: {actual_bonds:?} vs {expected_bonds:?}",
                    f[0]
                ));
            }
            let mut expected: Vec<Vec<(usize, usize)>> = f[3]
                .split(';')
                .filter(|s| !s.is_empty())
                .map(|s| {
                    s.split(',')
                        .map(|v| {
                            let b = &bonds[v.parse::<usize>().unwrap()];
                            pair(b[0], b[1])
                        })
                        .collect()
                })
                .collect();
            let mut actual: Vec<Vec<_>> = rings
                .rings
                .iter()
                .map(|r| {
                    r.iter()
                        .map(|&bi| pair(mol.bonds[bi].a, mol.bonds[bi].b))
                        .collect()
                })
                .collect();
            for r in &mut expected {
                r.sort_unstable();
            }
            for r in &mut actual {
                r.sort_unstable();
            }
            expected.sort_unstable();
            actual.sort_unstable();
            if actual != expected {
                errors.push(format!("{} rings: {actual:?} vs {expected:?}", f[0]));
            }
        }
        assert!(
            errors.is_empty(),
            "{} context errors:\n{}",
            errors.len(),
            errors.into_iter().take(30).collect::<Vec<_>>().join("\n")
        );
    }

    #[test]
    fn cleanup_is_idempotent_and_does_not_change_the_raw_graph_contract() {
        let source = "CN(=O)=O";
        let raw = crate::native_graph::graph(source).unwrap();
        assert_eq!(
            raw.atoms.iter().map(|a| a.charge).collect::<Vec<_>>(),
            [0, 0, 0, 0]
        );
        assert_eq!(
            raw.bonds
                .iter()
                .filter(|b| b.order == BondOrder::Double)
                .count(),
            2
        );
        let mut prepared = raw;
        assert!(apply(&mut prepared));
        // Independent literal graph C[N+]([O-])=O, preserving source indices.
        assert_eq!(
            prepared.atoms.iter().map(|a| a.charge).collect::<Vec<_>>(),
            [0, 1, -1, 0]
        );
        assert_eq!(
            prepared.bonds.iter().map(|b| b.order).collect::<Vec<_>>(),
            [BondOrder::Single, BondOrder::Single, BondOrder::Double]
        );
        assert!(
            !apply(&mut prepared),
            "a second cleanup must not change the graph"
        );
        let row = crate::workbench::analyze(source).unwrap();
        assert!((row.logp - -0.10710000000000008).abs() <= 1e-9);
        // Moving preparation into the global parser would violate this API.
        let raw_after = crate::native_graph::graph(source).unwrap();
        assert_eq!(
            raw_after.atoms.iter().map(|a| a.charge).collect::<Vec<_>>(),
            [0, 0, 0, 0]
        );
        assert_eq!(
            raw_after
                .bonds
                .iter()
                .filter(|b| b.order == BondOrder::Double)
                .count(),
            2
        );
    }
}
