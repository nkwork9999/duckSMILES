//! Topological polar surface area (TPSA) descriptor.
//!
//! This follows RDKit's `MolSurf.cpp` TPSA atom-contribution implementation:
//! count each hetero atom's local bond environment, assign an Ertl-style
//! contribution, then sum the contributions. The public DuckDB function uses
//! RDKit's default scope: nitrogen and oxygen atoms only.

use crate::parser::{BondOrder, Molecule};

#[derive(Clone, Copy, Debug, Default)]
struct AtomEnv {
    n_nbrs: i32,
    n_sing: i32,
    n_doub: i32,
    n_trip: i32,
    n_arom: i32,
    n_hs: i32,
    in_three_ring: bool,
}

/// Compute RDKit-default TPSA: N/O atom contributions, excluding S/P.
pub fn calc_tpsa(mol: &Molecule) -> f64 {
    calc_tpsa_with_options(mol, false)
}

fn calc_tpsa_with_options(mol: &Molecule, include_s_and_p: bool) -> f64 {
    if mol.atoms.is_empty() {
        return f64::NAN;
    }

    let envs = atom_environments(mol);
    mol.atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            atom_contribution(atom.symbol.as_str(), atom.charge, &envs[i], include_s_and_p)
        })
        .sum()
}

fn atom_environments(mol: &Molecule) -> Vec<AtomEnv> {
    let mut envs = vec![AtomEnv::default(); mol.atoms.len()];

    for bond in &mol.bonds {
        let a_sym = mol.atoms[bond.a].symbol.as_str();
        let b_sym = mol.atoms[bond.b].symbol.as_str();

        if a_sym == "H" {
            envs[bond.b].n_hs += 1;
            continue;
        }
        if b_sym == "H" {
            envs[bond.a].n_hs += 1;
            continue;
        }

        envs[bond.a].n_nbrs += 1;
        envs[bond.b].n_nbrs += 1;
        match bond.order {
            BondOrder::Single => {
                envs[bond.a].n_sing += 1;
                envs[bond.b].n_sing += 1;
            }
            BondOrder::Double => {
                envs[bond.a].n_doub += 1;
                envs[bond.b].n_doub += 1;
            }
            BondOrder::Triple => {
                envs[bond.a].n_trip += 1;
                envs[bond.b].n_trip += 1;
            }
            BondOrder::Aromatic => {
                envs[bond.a].n_arom += 1;
                envs[bond.b].n_arom += 1;
            }
        }
    }

    for (i, atom) in mol.atoms.iter().enumerate() {
        envs[i].n_hs += atom.hydrogen.max(0);
        envs[i].in_three_ring = is_atom_in_three_membered_ring(mol, i);
    }

    envs
}

fn is_atom_in_three_membered_ring(mol: &Molecule, atom_idx: usize) -> bool {
    let neighbors = mol
        .neighbors(atom_idx)
        .into_iter()
        .map(|(idx, _)| idx)
        .filter(|idx| mol.atoms[*idx].symbol != "H")
        .collect::<Vec<_>>();

    for i in 0..neighbors.len() {
        for j in (i + 1)..neighbors.len() {
            if are_bonded(mol, neighbors[i], neighbors[j]) {
                return true;
            }
        }
    }
    false
}

fn are_bonded(mol: &Molecule, a: usize, b: usize) -> bool {
    mol.bonds
        .iter()
        .any(|bond| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
}

fn atom_contribution(symbol: &str, charge: i32, env: &AtomEnv, include_s_and_p: bool) -> f64 {
    match symbol {
        "N" => nitrogen_contribution(charge, env),
        "O" => oxygen_contribution(charge, env),
        "P" if include_s_and_p => phosphorus_contribution(charge, env),
        "S" if include_s_and_p => sulfur_contribution(charge, env),
        _ => 0.0,
    }
}

fn nitrogen_contribution(chg: i32, e: &AtomEnv) -> f64 {
    let mut tmp = -1.0;
    match e.n_nbrs {
        1 => {
            if e.n_hs == 0 && chg == 0 && e.n_trip == 1 {
                tmp = 23.79;
            } else if e.n_hs == 1 && chg == 0 && e.n_doub == 1 {
                tmp = 23.85;
            } else if e.n_hs == 2 && chg == 0 && e.n_sing == 1 {
                tmp = 26.02;
            } else if e.n_hs == 2 && chg == 1 && e.n_doub == 1 {
                tmp = 25.59;
            } else if e.n_hs == 3 && chg == 1 && e.n_sing == 1 {
                tmp = 27.64;
            }
        }
        2 => {
            if e.n_hs == 0 && chg == 0 && e.n_sing == 1 && e.n_doub == 1 {
                tmp = 12.36;
            } else if e.n_hs == 0 && chg == 0 && e.n_trip == 1 && e.n_doub == 1 {
                tmp = 13.60;
            } else if e.n_hs == 1 && chg == 0 && e.n_sing == 2 && e.in_three_ring {
                tmp = 21.94;
            } else if e.n_hs == 1 && chg == 0 && e.n_sing == 2 && !e.in_three_ring {
                tmp = 12.03;
            } else if e.n_hs == 0 && chg == 1 && e.n_trip == 1 && e.n_sing == 1 {
                tmp = 4.36;
            } else if e.n_hs == 1 && chg == 1 && e.n_doub == 1 && e.n_sing == 1 {
                tmp = 13.97;
            } else if e.n_hs == 2 && chg == 1 && e.n_sing == 2 {
                tmp = 16.61;
            } else if e.n_hs == 0 && chg == 0 && e.n_arom == 2 {
                tmp = 12.89;
            } else if e.n_hs == 1 && chg == 0 && e.n_arom == 2 {
                tmp = 15.79;
            } else if e.n_hs == 1 && chg == 1 && e.n_arom == 2 {
                tmp = 14.14;
            }
        }
        3 => {
            if e.n_hs == 0 && chg == 0 && e.n_sing == 3 && e.in_three_ring {
                tmp = 3.01;
            } else if e.n_hs == 0 && chg == 0 && e.n_sing == 3 && !e.in_three_ring {
                tmp = 3.24;
            } else if e.n_hs == 0 && chg == 0 && e.n_sing == 1 && e.n_doub == 2 {
                tmp = 11.68;
            } else if e.n_hs == 0 && chg == 1 && e.n_sing == 2 && e.n_doub == 1 {
                tmp = 3.01;
            } else if e.n_hs == 1 && chg == 1 && e.n_sing == 3 {
                tmp = 4.44;
            } else if e.n_hs == 0 && chg == 0 && e.n_arom == 3 {
                tmp = 4.41;
            } else if e.n_hs == 0 && chg == 0 && e.n_sing == 1 && e.n_arom == 2 {
                tmp = 4.93;
            } else if e.n_hs == 0 && chg == 0 && e.n_doub == 1 && e.n_arom == 2 {
                tmp = 8.39;
            } else if e.n_hs == 0 && chg == 1 && e.n_arom == 3 {
                tmp = 4.10;
            } else if e.n_hs == 0 && chg == 1 && e.n_sing == 1 && e.n_arom == 2 {
                tmp = 3.88;
            }
        }
        4 => {
            if e.n_hs == 0 && e.n_sing == 4 && chg == 1 {
                tmp = 0.0;
            }
        }
        _ => {}
    }

    if tmp < 0.0 {
        tmp = 30.5 - (e.n_nbrs as f64) * 8.2 + (e.n_hs as f64) * 1.5;
        if tmp < 0.0 {
            tmp = 0.0;
        }
    }
    tmp
}

fn oxygen_contribution(chg: i32, e: &AtomEnv) -> f64 {
    let mut tmp = -1.0;
    match e.n_nbrs {
        1 => {
            if e.n_hs == 0 && chg == 0 && e.n_doub == 1 {
                tmp = 17.07;
            } else if e.n_hs == 1 && chg == 0 && e.n_sing == 1 {
                tmp = 20.23;
            } else if e.n_hs == 0 && chg == -1 && e.n_sing == 1 {
                tmp = 23.06;
            }
        }
        2 => {
            if e.n_hs == 0 && chg == 0 && e.n_sing == 2 && e.in_three_ring {
                tmp = 12.53;
            } else if e.n_hs == 0 && chg == 0 && e.n_sing == 2 && !e.in_three_ring {
                tmp = 9.23;
            } else if e.n_hs == 0 && chg == 0 && e.n_arom == 2 {
                tmp = 13.14;
            }
        }
        _ => {}
    }

    if tmp < 0.0 {
        tmp = 28.5 - (e.n_nbrs as f64) * 8.6 + (e.n_hs as f64) * 1.5;
        if tmp < 0.0 {
            tmp = 0.0;
        }
    }
    tmp
}

fn phosphorus_contribution(chg: i32, e: &AtomEnv) -> f64 {
    match e.n_nbrs {
        2 if e.n_hs == 0 && chg == 0 && e.n_sing == 1 && e.n_doub == 1 => 34.14,
        3 if e.n_hs == 0 && chg == 0 && e.n_sing == 3 => 13.59,
        3 if e.n_hs == 1 && chg == 0 && e.n_sing == 2 && e.n_doub == 1 => 23.47,
        4 if e.n_hs == 0 && chg == 0 && e.n_sing == 3 && e.n_doub == 1 => 9.81,
        _ => 0.0,
    }
}

fn sulfur_contribution(chg: i32, e: &AtomEnv) -> f64 {
    match e.n_nbrs {
        1 if e.n_hs == 0 && chg == 0 && e.n_doub == 1 => 32.09,
        1 if e.n_hs == 1 && chg == 0 && e.n_sing == 1 => 38.80,
        2 if e.n_hs == 0 && chg == 0 && e.n_sing == 2 => 25.30,
        2 if e.n_hs == 0 && chg == 0 && e.n_arom == 2 => 28.24,
        3 if e.n_hs == 0 && chg == 0 && e.n_arom == 2 && e.n_doub == 1 => 21.70,
        3 if e.n_hs == 0 && chg == 0 && e.n_sing == 2 && e.n_doub == 1 => 19.21,
        4 if e.n_hs == 0 && chg == 0 && e.n_sing == 2 && e.n_doub == 2 => 8.38,
        _ => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;
    use crate::test_fixtures;

    fn tpsa_of(smiles: &str) -> f64 {
        calc_tpsa(&parse(smiles).expect("parse smiles"))
    }

    fn assert_close(smiles: &str, expected: f64) {
        let got = tpsa_of(smiles);
        assert!(
            (got - expected).abs() < 0.01,
            "TPSA({}) = {:.4}, expected {:.4}",
            smiles,
            got,
            expected
        );
    }

    #[test]
    fn tpsa_named_fixture_regression() {
        let cases = [
            ("methane", test_fixtures::METHANE, 0.00),
            ("benzene", test_fixtures::BENZENE, 0.00),
            ("water", test_fixtures::WATER, 31.50),
            ("ethanol", test_fixtures::ETHANOL, 20.23),
            ("phenol", test_fixtures::PHENOL, 20.23),
            ("aspirin", test_fixtures::ASPIRIN, 63.60),
            ("salicylic acid", test_fixtures::SALICYLIC_ACID, 57.53),
            ("paracetamol", test_fixtures::PARACETAMOL, 49.33),
            ("ibuprofen", test_fixtures::IBUPROFEN, 37.30),
            ("cholesterol", test_fixtures::CHOLESTEROL, 20.23),
        ];

        for (name, smiles, expected) in cases {
            let got = tpsa_of(smiles);
            assert!(
                (got - expected).abs() < 0.01,
                "TPSA fixture {name} ({smiles}) = {got:.4}, expected {expected:.4}"
            );
        }
    }

    #[test]
    fn tpsa_simple_non_fixture_oxygen_cases() {
        assert_close("CC(=O)C", 17.07);
        assert_close("CC(=O)O", 37.30);
    }

    #[test]
    fn tpsa_nitrogen_cases() {
        assert_close("C#N", 23.79);
        assert_close("N", 35.00);
        assert_close("Nc1ccccc1", 26.02);
        assert_close("c1ccncc1", 12.89);
    }

    #[test]
    fn tpsa_three_membered_ring_cases() {
        assert_close("COC", 9.23);
        assert_close("C1CO1", 12.53);
    }

    #[test]
    fn tpsa_default_excludes_s_and_p() {
        let mol = parse("CS").unwrap();
        assert_eq!(calc_tpsa(&mol), 0.0);
        assert_eq!(calc_tpsa_with_options(&mol, true), 38.80);
    }
}
