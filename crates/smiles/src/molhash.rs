use crate::parser::{Atom, Bond, BondOrder, Molecule};
use crate::scaffold::{generic_scaffold_smiles, murcko_scaffold_smiles};

const METHODS: &[&str] = &[
    "canonical_smiles",
    "formula",
    "net_charge",
    "degree_vector",
    "atom_bond_counts",
    "element_graph",
    "bond_order_graph",
    "anonymous_graph",
    "murcko_scaffold",
    "generic_scaffold",
];

pub fn methods_json() -> String {
    let mut out = String::from("[");
    for (i, method) in METHODS.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push('"');
        out.push_str(method);
        out.push('"');
    }
    out.push(']');
    out
}

pub fn mol_hash(mol: &Molecule, method: &str) -> Option<String> {
    let method = normalize_method(method);
    match method.as_str() {
        "canonical_smiles" | "canonical" | "smiles" => Some(mol.canonical_smiles()),
        "formula" => Some(mol.formula()),
        "net_charge" | "charge" => Some(mol.total_charge().to_string()),
        "degree_vector" | "degreevector" => Some(degree_vector(mol)),
        "atom_bond_counts" | "counts" => Some(format!(
            "{}:{}",
            mol.heavy_atom_count(),
            heavy_bond_count(mol)
        )),
        "element_graph" | "elementgraph" => Some(element_graph(mol).canonical_smiles()),
        "bond_order_graph" | "bondordergraph" => Some(bond_order_graph(mol).canonical_smiles()),
        "anonymous_graph" | "anonymousgraph" => Some(anonymous_graph(mol).canonical_smiles()),
        "murcko_scaffold" | "murcko" | "scaffold" => Some(murcko_scaffold_smiles(mol)),
        "generic_scaffold" | "generic" | "framework" => Some(generic_scaffold_smiles(mol)),
        _ => None,
    }
}

fn normalize_method(method: &str) -> String {
    method
        .trim()
        .chars()
        .filter_map(|ch| match ch {
            '-' | ' ' | '.' => Some('_'),
            '_' => Some('_'),
            ch if ch.is_ascii_alphanumeric() => Some(ch.to_ascii_lowercase()),
            _ => None,
        })
        .collect()
}

fn clean_atom(symbol: &str) -> Atom {
    Atom {
        symbol: symbol.to_string(),
        hydrogen: 0,
        charge: 0,
        aromatic: false,
        in_bracket: false,
        isotope: None,
        atom_map: None,
        chirality: None,
    }
}

fn heavy_index_map(mol: &Molecule) -> Vec<usize> {
    let mut map = vec![usize::MAX; mol.atoms.len()];
    let mut next = 0;
    for (idx, atom) in mol.atoms.iter().enumerate() {
        if atom.symbol != "H" {
            map[idx] = next;
            next += 1;
        }
    }
    map
}

fn graph_from<F>(mol: &Molecule, atom_symbol: F, keep_order: bool) -> Molecule
where
    F: Fn(&Atom) -> String,
{
    let index_map = heavy_index_map(mol);
    let mut atoms = Vec::new();
    for atom in &mol.atoms {
        if atom.symbol != "H" {
            atoms.push(clean_atom(&atom_symbol(atom)));
        }
    }

    let mut bonds = Vec::new();
    for bond in &mol.bonds {
        let a = index_map[bond.a];
        let b = index_map[bond.b];
        if a != usize::MAX && b != usize::MAX {
            bonds.push(Bond {
                a,
                b,
                order: if keep_order {
                    bond.order
                } else {
                    BondOrder::Single
                },
            });
        }
    }

    Molecule {
        bond_count: bonds.len() as i32,
        atoms,
        bonds,
    }
}

fn element_graph(mol: &Molecule) -> Molecule {
    graph_from(mol, |atom| atom.symbol.clone(), false)
}

fn bond_order_graph(mol: &Molecule) -> Molecule {
    let index_map = heavy_index_map(mol);
    let mut atoms = Vec::new();
    for atom in &mol.atoms {
        if atom.symbol == "H" {
            continue;
        }
        atoms.push(Atom {
            symbol: atom.symbol.clone(),
            hydrogen: 0,
            charge: 0,
            aromatic: atom.aromatic,
            in_bracket: false,
            isotope: None,
            atom_map: None,
            chirality: None,
        });
    }

    let mut bonds = Vec::new();
    for bond in &mol.bonds {
        let a = index_map[bond.a];
        let b = index_map[bond.b];
        if a != usize::MAX && b != usize::MAX {
            bonds.push(Bond {
                a,
                b,
                order: bond.order,
            });
        }
    }

    Molecule {
        bond_count: bonds.len() as i32,
        atoms,
        bonds,
    }
}

fn anonymous_graph(mol: &Molecule) -> Molecule {
    graph_from(mol, |_| "C".to_string(), false)
}

fn heavy_bond_count(mol: &Molecule) -> usize {
    mol.bonds
        .iter()
        .filter(|bond| mol.atoms[bond.a].symbol != "H" && mol.atoms[bond.b].symbol != "H")
        .count()
}

fn degree_vector(mol: &Molecule) -> String {
    let mut degrees = Vec::new();
    for idx in 0..mol.atoms.len() {
        if mol.atoms[idx].symbol == "H" {
            continue;
        }
        let degree = mol
            .neighbors(idx)
            .iter()
            .filter(|(nbr, _)| mol.atoms[*nbr].symbol != "H")
            .count();
        degrees.push(degree);
    }
    degrees.sort_unstable();
    degrees
        .iter()
        .map(|degree| degree.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn molhash_core_methods_are_stable() {
        let mol = parse("CC(=O)[O-].[Na+]").unwrap();
        assert_eq!(mol_hash(&mol, "formula").unwrap(), "C2H3NaO2");
        assert_eq!(mol_hash(&mol, "net_charge").unwrap(), "0");
        assert_eq!(mol_hash(&mol, "degree_vector").unwrap(), "0,1,1,1,3");
        assert_eq!(mol_hash(&mol, "atom_bond_counts").unwrap(), "5:3");
        assert_eq!(mol_hash(&mol, "element_graph").unwrap(), "C(C)(O)O.[Na]");
        assert_eq!(mol_hash(&mol, "anonymous_graph").unwrap(), "C.C(C)(C)C");
    }

    #[test]
    fn molhash_scaffold_methods_reuse_existing_scaffolds() {
        let mol = parse("Cc1ccccc1").unwrap();
        assert_eq!(mol_hash(&mol, "murcko_scaffold").unwrap(), "c1ccccc1");
        assert_eq!(mol_hash(&mol, "generic_scaffold").unwrap(), "C1CCCCC1");
        assert!(methods_json().contains("element_graph"));
    }
}
