use crate::parser::{BondOrder, Molecule};
use crate::scaffold::murcko_scaffold_mol;

#[derive(Debug, Clone, Default)]
struct Mapping {
    pairs: Vec<(usize, usize)>,
}

pub fn mcs_smarts(a: &Molecule, b: &Molecule) -> String {
    let mapping = maximum_common_mapping(a, b);
    mapping_to_smarts(a, &mapping)
}

pub fn mcs_json(a: &Molecule, b: &Molecule) -> String {
    let mapping = maximum_common_mapping(a, b);
    let smarts = mapping_to_smarts(a, &mapping);
    let mut atoms_a = mapping.pairs.iter().map(|(x, _)| *x).collect::<Vec<_>>();
    atoms_a.sort_unstable();
    let mut atoms_b = mapping.pairs.iter().map(|(_, y)| *y).collect::<Vec<_>>();
    atoms_b.sort_unstable();
    let bond_count = mapped_bond_count(a, &mapping);

    let mut out = String::from("{\"num_atoms\":");
    out.push_str(&mapping.pairs.len().to_string());
    out.push_str(",\"num_bonds\":");
    out.push_str(&bond_count.to_string());
    out.push_str(",\"smarts\":\"");
    out.push_str(&json_escape(&smarts));
    out.push_str("\",\"atoms_a\":");
    push_usize_array_1based(&mut out, &atoms_a);
    out.push_str(",\"atoms_b\":");
    push_usize_array_1based(&mut out, &atoms_b);
    out.push_str(",\"atom_pairs\":[");
    for (i, (a_idx, b_idx)) in mapping.pairs.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"a\":");
        out.push_str(&(a_idx + 1).to_string());
        out.push_str(",\"b\":");
        out.push_str(&(b_idx + 1).to_string());
        out.push('}');
    }
    out.push_str("]}");
    out
}

pub fn scaffold_network_json(mol: &Molecule) -> String {
    let original = mol.canonical_smiles();
    let murcko = murcko_scaffold_mol(mol);
    let murcko_smiles = murcko
        .as_ref()
        .map(|scaffold| scaffold.canonical_smiles())
        .unwrap_or_default();
    let generic = murcko.as_ref().map(genericize).unwrap_or_default();
    let ring_fragments = murcko
        .as_ref()
        .map(ring_fragment_smiles)
        .unwrap_or_default();

    let mut nodes = Vec::new();
    push_node(&mut nodes, "molecule", &original);
    if !murcko_smiles.is_empty() {
        push_node(&mut nodes, "murcko_scaffold", &murcko_smiles);
    }
    if !generic.is_empty() {
        push_node(&mut nodes, "generic_scaffold", &generic);
    }
    for ring in &ring_fragments {
        push_node(&mut nodes, "ring_system", ring);
    }
    nodes.sort();
    nodes.dedup();

    let mut edges = Vec::new();
    if !murcko_smiles.is_empty() {
        push_edge(&mut edges, &original, &murcko_smiles, "murcko");
    }
    if !generic.is_empty() && !murcko_smiles.is_empty() {
        push_edge(&mut edges, &murcko_smiles, &generic, "generic");
    }
    for ring in &ring_fragments {
        if !murcko_smiles.is_empty() && ring != &murcko_smiles {
            push_edge(&mut edges, &murcko_smiles, ring, "ring_system");
        }
    }
    edges.sort();
    edges.dedup();

    let mut out = String::from("{\"nodes\":[");
    for (i, node) in nodes.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(node);
    }
    out.push_str("],\"edges\":[");
    for (i, edge) in edges.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(edge);
    }
    out.push_str("]}");
    out
}

fn maximum_common_mapping(a: &Molecule, b: &Molecule) -> Mapping {
    let (small, large, swapped) = if a.heavy_atom_count() <= b.heavy_atom_count() {
        (a, b, false)
    } else {
        (b, a, true)
    };

    let mut best = Mapping::default();
    for s_idx in 0..small.atoms.len() {
        if small.atoms[s_idx].symbol == "H" {
            continue;
        }
        for l_idx in 0..large.atoms.len() {
            if large.atoms[l_idx].symbol == "H" || !atoms_compatible(small, large, s_idx, l_idx) {
                continue;
            }
            let mut mapping = vec![usize::MAX; small.atoms.len()];
            let mut used_large = vec![false; large.atoms.len()];
            mapping[s_idx] = l_idx;
            used_large[l_idx] = true;
            search_connected(small, large, &mut mapping, &mut used_large, &mut best);
        }
    }

    if swapped {
        Mapping {
            pairs: best.pairs.into_iter().map(|(x, y)| (y, x)).collect(),
        }
    } else {
        best
    }
}

fn search_connected(
    small: &Molecule,
    large: &Molecule,
    mapping: &mut [usize],
    used_large: &mut [bool],
    best: &mut Mapping,
) {
    update_best(mapping, best);

    let frontier = frontier_atoms(small, mapping);
    for s_idx in frontier {
        for l_idx in 0..large.atoms.len() {
            if used_large[l_idx]
                || large.atoms[l_idx].symbol == "H"
                || !atoms_compatible(small, large, s_idx, l_idx)
                || !compatible_with_mapped_neighbors(small, large, mapping, s_idx, l_idx)
            {
                continue;
            }
            mapping[s_idx] = l_idx;
            used_large[l_idx] = true;
            search_connected(small, large, mapping, used_large, best);
            used_large[l_idx] = false;
            mapping[s_idx] = usize::MAX;
        }
    }
}

fn update_best(mapping: &[usize], best: &mut Mapping) {
    let mut pairs = mapping
        .iter()
        .enumerate()
        .filter_map(|(idx, &mapped)| {
            if mapped == usize::MAX {
                None
            } else {
                Some((idx, mapped))
            }
        })
        .collect::<Vec<_>>();
    pairs.sort_unstable();
    let candidate = Mapping { pairs };
    if candidate.pairs.len() > best.pairs.len()
        || (candidate.pairs.len() == best.pairs.len() && candidate.pairs < best.pairs)
    {
        *best = candidate;
    }
}

fn frontier_atoms(mol: &Molecule, mapping: &[usize]) -> Vec<usize> {
    let mut frontier = Vec::new();
    for (idx, mapped) in mapping.iter().enumerate() {
        if *mapped != usize::MAX || mol.atoms[idx].symbol == "H" {
            continue;
        }
        if mol
            .neighbors(idx)
            .iter()
            .any(|(nbr, _)| mapping[*nbr] != usize::MAX)
        {
            frontier.push(idx);
        }
    }
    frontier.sort_unstable();
    frontier
}

fn atoms_compatible(a: &Molecule, b: &Molecule, a_idx: usize, b_idx: usize) -> bool {
    let aa = &a.atoms[a_idx];
    let bb = &b.atoms[b_idx];
    aa.symbol == bb.symbol && aa.aromatic == bb.aromatic && aa.charge == bb.charge
}

fn compatible_with_mapped_neighbors(
    small: &Molecule,
    large: &Molecule,
    mapping: &[usize],
    s_idx: usize,
    l_idx: usize,
) -> bool {
    for (s_nbr, s_order) in small.neighbors(s_idx) {
        let l_nbr = mapping[s_nbr];
        if l_nbr == usize::MAX {
            continue;
        }
        let Some(l_order) = bond_between(large, l_idx, l_nbr) else {
            return false;
        };
        if !bond_compatible(s_order, l_order) {
            return false;
        }
    }
    true
}

fn bond_between(mol: &Molecule, a: usize, b: usize) -> Option<BondOrder> {
    mol.bonds.iter().find_map(|bond| {
        if (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a) {
            Some(bond.order)
        } else {
            None
        }
    })
}

fn bond_compatible(a: BondOrder, b: BondOrder) -> bool {
    a == b
        || (a == BondOrder::Aromatic && b == BondOrder::Single)
        || (a == BondOrder::Single && b == BondOrder::Aromatic)
}

fn mapped_bond_count(mol: &Molecule, mapping: &Mapping) -> usize {
    let atoms = mapping
        .pairs
        .iter()
        .map(|(idx, _)| *idx)
        .collect::<Vec<_>>();
    mol.bonds
        .iter()
        .filter(|bond| atoms.contains(&bond.a) && atoms.contains(&bond.b))
        .count()
}

fn mapping_to_smarts(mol: &Molecule, mapping: &Mapping) -> String {
    if mapping.pairs.is_empty() {
        return String::new();
    }
    let atoms = mapping
        .pairs
        .iter()
        .map(|(idx, _)| *idx)
        .collect::<Vec<_>>();
    mol.subgraph_from_atoms(&atoms)
        .map(|subgraph| subgraph.canonical_smiles())
        .unwrap_or_default()
}

fn genericize(mol: &Molecule) -> String {
    let mut generic = mol.clone();
    for atom in &mut generic.atoms {
        atom.symbol = "C".to_string();
        atom.hydrogen = 0;
        atom.charge = 0;
        atom.aromatic = false;
        atom.in_bracket = false;
        atom.isotope = None;
        atom.atom_map = None;
        atom.chirality = None;
    }
    for bond in &mut generic.bonds {
        bond.order = BondOrder::Single;
    }
    generic.canonical_smiles()
}

fn ring_fragment_smiles(mol: &Molecule) -> Vec<String> {
    let ring = mol.ring_info();
    let mut fragments = Vec::new();
    for ring_bonds in ring.rings {
        let mut atoms = Vec::new();
        for bond_idx in ring_bonds {
            if let Some(bond) = mol.bonds.get(bond_idx) {
                atoms.push(bond.a);
                atoms.push(bond.b);
            }
        }
        atoms.sort_unstable();
        atoms.dedup();
        if let Some(subgraph) = mol.subgraph_from_atoms(&atoms) {
            fragments.push(subgraph.canonical_smiles());
        }
    }
    fragments.sort();
    fragments.dedup();
    fragments
}

fn push_node(nodes: &mut Vec<String>, kind: &str, smiles: &str) {
    let mut node = String::from("{\"kind\":\"");
    node.push_str(kind);
    node.push_str("\",\"smiles\":\"");
    node.push_str(&json_escape(smiles));
    node.push_str("\"}");
    nodes.push(node);
}

fn push_edge(edges: &mut Vec<String>, from: &str, to: &str, kind: &str) {
    let mut edge = String::from("{\"from\":\"");
    edge.push_str(&json_escape(from));
    edge.push_str("\",\"to\":\"");
    edge.push_str(&json_escape(to));
    edge.push_str("\",\"kind\":\"");
    edge.push_str(kind);
    edge.push_str("\"}");
    edges.push(edge);
}

fn push_usize_array_1based(out: &mut String, values: &[usize]) {
    out.push('[');
    for (i, value) in values.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str(&(value + 1).to_string());
    }
    out.push(']');
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            ch if ch <= '\u{1f}' => out.push_str(&format!("\\u{:04x}", ch as u32)),
            _ => out.push(ch),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    #[test]
    fn mcs_finds_connected_common_subgraph() {
        let a = parse("CC(=O)O").unwrap();
        let b = parse("CC(=O)N").unwrap();
        assert_eq!(mcs_smarts(&a, &b), "C(C)=O");
        let json = mcs_json(&a, &b);
        assert!(json.contains("\"num_atoms\":3"));
        assert!(json.contains("\"smarts\":\"C(C)=O\""));
    }

    #[test]
    fn scaffold_network_reports_scaffold_and_rings() {
        let mol = parse("Cc1ccccc1").unwrap();
        let json = scaffold_network_json(&mol);
        assert!(json.contains("\"murcko_scaffold\""));
        assert!(json.contains("c1ccccc1"));
        assert!(json.contains("\"generic_scaffold\""));
    }
}
