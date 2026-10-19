//! Bounded analysis entry point for native workbench hosts.
// For simple undirected adjacency, sum |characteristic coefficients| is at
// most (1 + sqrt(2m/n))^n. These budgets bound weighted information below
// 3^640 * log2(641), within finite binary64. Byte/row ABI guards still apply.
pub const MAX_ATOMS: usize = 640;
pub const MAX_BONDS: usize = 1280;

pub struct BasicDescriptors {
    pub mw: f64,
    pub exact_mass: f64,
    pub formula: String,
    pub heavy_atoms: usize,
    pub bonds: usize,
    pub hbd: usize,
    pub hba: usize,
    pub tpsa: f64,
    pub fragments: usize,
    pub formal_charge: i32,
    pub logp: f64,
    pub rotatable_bonds: usize,
    pub fraction_csp3: f64,
    pub descriptor_values: Vec<Option<f64>>,
}

pub(crate) fn prepare(input: &str) -> Option<(crate::parser::Molecule, crate::parser::RingInfo)> {
    let mut mol = crate::workbench_cx::graph(input)?;
    crate::workbench_cleanup::apply(&mut mol);
    crate::workbench_cleanup::coordination(&mut mol)?;
    let rings = crate::workbench_rings::calculate(&mol)?;
    crate::workbench_cleanup::aromatic_hydrogens(&mut mol);
    crate::workbench_cleanup::check_explicit_valence(&mol)?;
    crate::workbench_kekulize::apply(&mut mol, &rings)?;
    let donations = crate::workbench_cleanup::aromatic_donations(&mol, &rings)?;
    mol.perceive_aromaticity_with_donations(&rings, Some(&donations));
    crate::workbench_cleanup::non_tetrahedral_chirality(&mut mol);
    Some((mol, rings))
}

pub fn analyze(input: &str) -> Option<BasicDescriptors> {
    let (mol, rings) = prepare(input)?;
    let mw = mol.molecular_weight();
    let exact_mass = mol.exact_mass();
    if !mw.is_finite() || !exact_mass.is_finite() {
        return None;
    }
    let crippen = crate::logp_crippen::calc_with_atom_contributions(&mol);
    // Reuse the parser's ring topology for acceptors, strict rotors, ring
    // classes and the shared fragment matcher.
    let mut result = BasicDescriptors {
        mw,
        exact_mass,
        formula: mol.formula(),
        heavy_atoms: mol.atoms.iter().filter(|a| a.symbol != "H").count(),
        bonds: mol.bonds.len(),
        hbd: crate::num_h_donors(&mol),
        hba: crate::num_h_acceptors_with_ring_bonds(&mol, &rings.bond_in_ring),
        tpsa: crate::tpsa::calc_tpsa_with_ring_info(&mol, &rings),
        fragments: mol.components().len(),
        formal_charge: mol.atoms.iter().map(|a| a.charge).sum(),
        logp: crippen.logp,
        rotatable_bonds: crate::num_rotatable_bonds_with_ring_bonds(&mol, &rings.bond_in_ring),
        fraction_csp3: crate::fraction_csp3(&mol),
        descriptor_values: Vec::new(),
    };
    result.descriptor_values = crate::descriptor_profile::calculate_with_ring_info(
        &mol,
        &result,
        crippen.mr,
        &crippen.atoms,
        &rings,
    );
    Some(result)
}

/// UTF-8, newline-delimited input; one result per source row. Required length is
/// returned without writing when capacity is insufficient. Negative = bad ABI input.
///
/// # Safety
/// Non-null pointers must address readable/writable regions of the given lengths.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_workbench_analyze(
    input: *const u8,
    len: usize,
    out: *mut u8,
    capacity: usize,
) -> i32 {
    std::panic::catch_unwind(|| {
        if input.is_null() || len > 1_048_576 || (out.is_null() && capacity != 0) { return -1; }
        // SAFETY: the C caller guarantees readable input and a disjoint writable output.
        let bytes = unsafe { std::slice::from_raw_parts(input, len) };
        let Ok(text) = std::str::from_utf8(bytes) else { return -1; };
        let row_count = if text.is_empty() { 0 } else { text.split('\n').count() };
        if row_count > 10_000 { return -2; }
        use std::fmt::Write;
        let mut json = String::from("{\"schema\":1,\"rows\":[");
        for (index, line) in text.split('\n').take(row_count).enumerate() {
            if index != 0 { json.push(','); }
            match analyze(line.trim()) {
                None => { let _ = write!(json, "{{\"index\":{index},\"status\":\"invalid_or_unsupported\"}}"); }
                Some(d) => {
                    let _ = write!(json, "{{\"index\":{index},\"status\":\"ok\",\"formula\":\"{}\",\"mw\":{},\"exact_mass\":{},\"heavy_atoms\":{},\"bonds\":{},\"formal_charge\":{},\"fragments\":{},\"hbd\":{},\"hba\":{},\"tpsa\":{},\"logp\":{},\"rotatable_bonds\":{},\"fraction_csp3\":{}}}",
                        d.formula, d.mw, d.exact_mass, d.heavy_atoms, d.bonds, d.formal_charge,
                        d.fragments, d.hbd, d.hba, d.tpsa, d.logp, d.rotatable_bonds, d.fraction_csp3);
                    json.pop(); // append the catalog-ordered vector to the existing object
                    json.push_str(",\"descriptor_values\":[");
                    for (i, value) in d.descriptor_values.iter().enumerate() {
                        if i > 0 { json.push(','); }
                        match value { Some(value) => { let _ = write!(json, "{value}"); }, None => json.push_str("null") }
                    }
                    json.push_str("]}");
                }
            }
        }
        json.push_str("]}");
        let Ok(required) = i32::try_from(json.len()) else { return -2; };
        if capacity >= json.len() && !out.is_null() {
            // SAFETY: caller owns capacity bytes; the string owns a disjoint source.
            unsafe { std::ptr::copy_nonoverlapping(json.as_ptr(), out, json.len()); }
        }
        required
    }).unwrap_or(-3)
}

#[cfg(test)]
mod tests {
    #[test]
    fn large_input_domain_matches_declared_budgets_without_changing_raw_graph() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/large_domain.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let actual = super::prepare(f[1]);
            let expected = f[2] == "1";
            if actual.is_some() != expected {
                differences.push(format!(
                    "{}: accepted={} expected={expected}",
                    f[0],
                    actual.is_some()
                ));
            }
            if let Some((mol, _)) = actual {
                assert_eq!(
                    mol.atoms.len(),
                    f[3].parse::<usize>().unwrap(),
                    "{} atoms",
                    f[0]
                );
                assert_eq!(
                    mol.bonds.len(),
                    f[4].parse::<usize>().unwrap(),
                    "{} bonds",
                    f[0]
                );
                assert_eq!(mol.formula(), f[5], "{} formula", f[0]);
            }
            assert!(
                crate::native_graph::graph(f[1]).is_none(),
                "raw graph must keep its smaller domain"
            );
            count += 1;
        }
        assert_eq!(count, 6);
        assert!(differences.is_empty(), "{}", differences.join("\n"));
    }

    #[test]
    fn explicit_valence_acceptance_matches_pinned_sanitization() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/sanitization_valence.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let actual = super::prepare(fields[0]);
            let expected = fields[1] == "1";
            if actual.is_some() != expected {
                differences.push(format!(
                    "{}: accepted={} expected={expected}",
                    fields[0],
                    actual.is_some()
                ));
            } else if let Some((mol, _)) = actual {
                if mol.formula() != fields[2] {
                    differences.push(format!(
                        "{}: formula={} expected={}",
                        fields[0],
                        mol.formula(),
                        fields[2]
                    ));
                }
            }
            count += 1;
        }
        assert_eq!(count, 14603);
        assert!(
            differences.is_empty(),
            "{} acceptance/formula differences, first40:\n{}",
            differences.len(),
            differences
                .iter()
                .take(40)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    #[test]
    #[ignore = "external unselected100k source is supplied explicitly"]
    fn supplied_population_preparation_acceptance() {
        let source = std::env::var("CHEMOB_ACCEPTANCE_SMILES").expect("external source path");
        let text = std::fs::read_to_string(source).unwrap();
        let inputs: Vec<_> = text.lines().collect();
        assert_eq!(inputs.len(), 100_000);
        let unique: std::collections::HashSet<_> = inputs.iter().copied().collect();
        assert_eq!(
            unique.len(),
            inputs.len(),
            "original distinct structures required"
        );
        let unsupported: Vec<_> = inputs
            .iter()
            .enumerate()
            .filter_map(|(i, s)| super::prepare(s.trim()).is_none().then_some(i))
            .collect();
        println!("POPULATION_ACCEPTANCE_JSON {{\"rows\":{},\"accepted\":{},\"unsupported_indices\":{:?}}}",
            inputs.len(), inputs.len()-unsupported.len(), unsupported);
    }

    #[test]
    fn prepared_kekule_graphs_match_independent_sanitization() {
        let mut differences = Vec::new();
        let mut count = 0;
        for line in include_str!("../tests/fixtures/kekulize_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            let Some((mol, rings)) = super::prepare(f[0]) else {
                differences.push(format!("{}: unexpectedly unsupported", f[0]));
                continue;
            };
            let z = crate::electronic::atomic_numbers(&mol).unwrap();
            let states = crate::electronic::atom_states(&mol, &z);
            let mut fields = vec![mol.formula()];
            for column in 0..6 {
                fields.push(
                    mol.atoms
                        .iter()
                        .enumerate()
                        .map(|(i, a)| {
                            let explicit = crate::profile_h_is_explicit(a);
                            match column {
                                0 => a.hydrogen,
                                1 => i32::from(a.aromatic),
                                2 => {
                                    if explicit {
                                        a.hydrogen
                                    } else {
                                        0
                                    }
                                }
                                3 => {
                                    if explicit {
                                        0
                                    } else {
                                        a.hydrogen
                                    }
                                }
                                4 => states[i].valence,
                                _ => states[i].radicals,
                            }
                            .to_string()
                        })
                        .collect::<Vec<_>>()
                        .join(","),
                );
            }
            let flags = mol.aromatic_bond_flags(&rings);
            let mut edges: Vec<_> = mol
                .bonds
                .iter()
                .enumerate()
                .map(|(i, b)| {
                    use crate::parser::BondOrder::*;
                    let order = match b.order {
                        Single => 1,
                        Double => 2,
                        Triple => 3,
                        Aromatic => 12,
                        Dative => 17,
                    };
                    (b.a.min(b.b), b.a.max(b.b), order, i32::from(flags[i]))
                })
                .collect();
            edges.sort_unstable();
            fields.push(
                edges
                    .iter()
                    .map(|(a, b, o, f)| format!("{a},{b},{o},{f}"))
                    .collect::<Vec<_>>()
                    .join(";"),
            );
            let labels = [
                "formula",
                "attached H",
                "aromatic atoms",
                "explicit H",
                "implicit H",
                "valence",
                "radicals",
                "typed bonds",
            ];
            for ((actual, expected), label) in fields.iter().zip(&f[1..]).zip(labels) {
                if actual != expected {
                    differences.push(format!("{} {label}: {actual} vs {expected}", f[0]));
                }
            }
            count += 1;
        }
        assert_eq!(count, 176);
        assert!(
            differences.is_empty(),
            "{} differences, first30:\n{}",
            differences.len(),
            differences
                .iter()
                .take(30)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    #[test]
    fn impossible_aromatic_kekule_inputs_are_issue_rows() {
        for source in include_str!("../tests/fixtures/kekulize_invalid.txt").lines() {
            assert!(
                super::prepare(source).is_none(),
                "invalid aromatic input {source}"
            );
        }
    }

    #[test]
    fn dative_inputs_preserve_endpoint_direction_hydrogens_and_raw_grammar() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/dative_input.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let f: Vec<_> = line.split('\t').collect();
            assert!(
                crate::native_graph::graph(f[0]).is_none(),
                "legacy raw API must retain its grammar"
            );
            let (m, rings) = super::prepare(f[0]).expect(f[0]);
            assert_eq!(m.formula(), f[1], "{} formula", f[0]);
            let h = m
                .atoms
                .iter()
                .map(|a| a.hydrogen.to_string())
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(h, f[2], "{} attached H", f[0]);
            let mut edges: Vec<_> = m
                .bonds
                .iter()
                .filter(|b| b.order == crate::parser::BondOrder::Dative)
                .map(|b| (b.a, b.b))
                .collect();
            edges.sort_unstable();
            assert_eq!(
                edges
                    .iter()
                    .map(|(a, b)| format!("{a},{b}"))
                    .collect::<Vec<_>>()
                    .join(";"),
                f[3],
                "{} donor→acceptor",
                f[0]
            );
            assert_eq!(
                rings.rings.len().to_string(),
                f[4],
                "{} covalent ring count",
                f[0]
            );
            count += 1;
        }
        assert_eq!(count, 13);
        for source in ["N->", "N<-", "->[Cu]", "N-><-[Cu]", "N<->[Cu]", "N->(C)"] {
            assert!(super::prepare(source).is_none(), "malformed arrow {source}");
        }
    }

    #[test]
    fn prepared_hydrogen_vectors_match_reference_and_raw_contract_is_retained() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/valence_h_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let (mol, _) = super::prepare(fields[0]).expect(fields[0]);
            let hydrogens = mol
                .atoms
                .iter()
                .map(|a| a.hydrogen.to_string())
                .collect::<Vec<_>>()
                .join(",");
            assert_eq!(hydrogens, fields[1], "{} H distribution", fields[0]);
            count += 1;
        }
        assert_eq!(count, 41);
        let before = crate::native_graph::graph("C:C").unwrap().formula();
        assert_eq!(before, "C2H6");
        assert_eq!(super::analyze("C:C").unwrap().formula, "C2H4");
        assert_eq!(crate::native_graph::graph("C:C").unwrap().formula(), before);
    }

    #[test]
    fn blank_final_structure_is_retained_as_an_issue_row() {
        let input = b"CCO\n";
        let mut output = vec![0u8; 4096];
        let n = unsafe {
            super::ds_workbench_analyze(
                input.as_ptr(),
                input.len(),
                output.as_mut_ptr(),
                output.len(),
            )
        };
        assert!(n > 0);
        let json = std::str::from_utf8(&output[..n as usize]).unwrap();
        assert!(json.contains("\"index\":1,\"status\":\"invalid_or_unsupported\""));
    }
}

/// Cache namespace; bump whenever descriptor or graph semantics change.
pub const REVISION: &str = "ducksmiles-workbench-v42-cluster291-rdkit2026036";

/// # Safety
/// Non-null output must address capacity writable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_workbench_revision(out: *mut u8, capacity: usize) -> i32 {
    if out.is_null() && capacity != 0 {
        return -1;
    }
    if capacity >= REVISION.len() && !out.is_null() {
        // SAFETY: caller provides the documented writable region; source is static.
        unsafe {
            std::ptr::copy_nonoverlapping(REVISION.as_ptr(), out, REVISION.len());
        }
    }
    REVISION.len() as i32
}
