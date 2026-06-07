use crate::parser::{Atom, BondOrder, Molecule};

pub fn largest_fragment_smiles(mol: &Molecule) -> String {
    selected_components(mol, LargestMode::Largest)
        .into_iter()
        .next()
        .and_then(|component| mol.subgraph_from_atoms(&component))
        .map(|frag| frag.canonical_smiles())
        .unwrap_or_default()
}

pub fn strip_salts_smiles(mol: &Molecule) -> String {
    let components = selected_components(mol, LargestMode::KeepOrganic);
    canonical_join(mol, &components)
}

pub fn neutralize_charges_smiles(mol: &Molecule) -> String {
    let mut normalized = mol.clone();
    neutralize_charges_in_place(&mut normalized);
    normalized.canonical_smiles()
}

pub fn normalize_smiles(mol: &Molecule) -> String {
    let mut normalized = mol.clone();
    neutralize_charges_in_place(&mut normalized);
    normalized.perceive_aromaticity();
    normalized.canonical_smiles()
}

pub fn fragment_parent_smiles(mol: &Molecule) -> String {
    let stripped = selected_components(mol, LargestMode::KeepOrganic);
    let parent_component = if stripped.is_empty() {
        selected_components(mol, LargestMode::Largest)
    } else {
        choose_largest(stripped, mol)
            .into_iter()
            .collect::<Vec<_>>()
    };
    let Some(component) = parent_component.into_iter().next() else {
        return String::new();
    };
    let Some(mut parent) = mol.subgraph_from_atoms(&component) else {
        return String::new();
    };
    neutralize_charges_in_place(&mut parent);
    parent.perceive_aromaticity();
    parent.canonical_smiles()
}

#[derive(Clone, Copy)]
enum LargestMode {
    Largest,
    KeepOrganic,
}

fn selected_components(mol: &Molecule, mode: LargestMode) -> Vec<Vec<usize>> {
    let components = mol.components();
    if components.is_empty() {
        return Vec::new();
    }
    match mode {
        LargestMode::Largest => choose_largest(components, mol).into_iter().collect(),
        LargestMode::KeepOrganic => {
            let organic = components
                .iter()
                .filter(|component| component_has_carbon(mol, component))
                .cloned()
                .collect::<Vec<_>>();
            if organic.is_empty() {
                choose_largest(components, mol).into_iter().collect()
            } else {
                organic
            }
        }
    }
}

fn choose_largest(components: Vec<Vec<usize>>, mol: &Molecule) -> Option<Vec<usize>> {
    components.into_iter().max_by(|a, b| {
        component_score(mol, a)
            .cmp(&component_score(mol, b))
            .then_with(|| b.cmp(a))
    })
}

fn component_score(mol: &Molecule, component: &[usize]) -> (usize, usize, i64, usize) {
    let heavy = component
        .iter()
        .filter(|&&idx| mol.atoms[idx].symbol != "H")
        .count();
    let carbon = component
        .iter()
        .filter(|&&idx| mol.atoms[idx].symbol == "C")
        .count();
    let mass_milli = mol
        .subgraph_from_atoms(component)
        .map(|frag| (frag.molecular_weight() * 1000.0).round() as i64)
        .unwrap_or(0);
    (carbon, heavy, mass_milli, component.len())
}

fn component_has_carbon(mol: &Molecule, component: &[usize]) -> bool {
    component.iter().any(|&idx| mol.atoms[idx].symbol == "C")
}

fn canonical_join(mol: &Molecule, components: &[Vec<usize>]) -> String {
    let mut parts = components
        .iter()
        .filter_map(|component| mol.subgraph_from_atoms(component))
        .map(|fragment| fragment.canonical_smiles())
        .collect::<Vec<_>>();
    parts.sort();
    parts.join(".")
}

fn neutralize_charges_in_place(mol: &mut Molecule) {
    for idx in 0..mol.atoms.len() {
        if should_skip_ion_neutralization(&mol.atoms[idx]) {
            continue;
        }

        if mol.atoms[idx].charge < 0 {
            let add_h = (-mol.atoms[idx].charge) as i32;
            mol.atoms[idx].hydrogen += add_h;
            mol.atoms[idx].charge = 0;
            mol.atoms[idx].in_bracket = false;
            continue;
        }

        if mol.atoms[idx].charge > 0 {
            let remove_h = mol.atoms[idx].charge.min(mol.atoms[idx].hydrogen.max(0));
            mol.atoms[idx].hydrogen -= remove_h;
            mol.atoms[idx].charge -= remove_h;
            if mol.atoms[idx].charge == 0 {
                mol.atoms[idx].in_bracket = false;
            }
        }
    }

    normalize_nitro(mol);
}

fn should_skip_ion_neutralization(atom: &Atom) -> bool {
    matches!(
        atom.symbol.as_str(),
        "Li" | "Na"
            | "K"
            | "Rb"
            | "Cs"
            | "Mg"
            | "Ca"
            | "Sr"
            | "Ba"
            | "Al"
            | "Zn"
            | "Fe"
            | "Cu"
            | "F"
            | "Cl"
            | "Br"
            | "I"
    )
}

fn normalize_nitro(mol: &mut Molecule) {
    for n_idx in 0..mol.atoms.len() {
        if mol.atoms[n_idx].symbol != "N" {
            continue;
        }
        let oxygen_neighbors = mol
            .neighbors(n_idx)
            .into_iter()
            .filter(|(nbr, _)| mol.atoms[*nbr].symbol == "O")
            .collect::<Vec<_>>();
        let has_double_o = oxygen_neighbors
            .iter()
            .any(|(_, order)| *order == BondOrder::Double);
        let has_single_o = oxygen_neighbors
            .iter()
            .any(|(_, order)| *order == BondOrder::Single);
        if has_double_o && has_single_o && mol.atoms[n_idx].charge == 0 {
            mol.atoms[n_idx].charge = 1;
            if let Some((o_idx, _)) = oxygen_neighbors
                .iter()
                .find(|(_, order)| *order == BondOrder::Single)
            {
                mol.atoms[*o_idx].charge = -1;
                mol.atoms[*o_idx].hydrogen = 0;
                mol.atoms[*o_idx].in_bracket = true;
            }
            mol.atoms[n_idx].in_bracket = true;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn largest_fragment_prefers_organic_heavy_component() {
        let mol = parse("CC(=O)[O-].[Na+]").unwrap();
        assert_eq!(largest_fragment_smiles(&mol), "C(C)(=O)[O-]");
        assert_eq!(strip_salts_smiles(&mol), "C(C)(=O)[O-]");
    }

    #[test]
    fn neutralize_and_parent_are_canonical() {
        let mol = parse("CC(=O)[O-].[Na+]").unwrap();
        assert_eq!(neutralize_charges_smiles(&mol), "C(C)(=O)O.[Na+]");
        assert_eq!(fragment_parent_smiles(&mol), "C(C)(=O)O");

        let chloride = parse("[Cl-]").unwrap();
        assert_eq!(neutralize_charges_smiles(&chloride), "[Cl-]");
    }

    #[test]
    fn strip_salts_keeps_multiple_organic_fragments() {
        let mol = parse("CCO.CN.[Cl-]").unwrap();
        assert_eq!(strip_salts_smiles(&mol), "C(C)O.CN");
    }
}
