use std::collections::{BTreeMap, HashMap};
use crate::weights;

const ORGANIC_SUBSET: &[&str] = &["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];

fn is_organic(sym: &str) -> bool {
    ORGANIC_SUBSET.contains(&sym)
}

fn is_aromatic_char(c: char) -> bool {
    matches!(c, 'b' | 'c' | 'n' | 'o' | 'p' | 's')
}

fn default_valence(sym: &str, aromatic: bool) -> i32 {
    match sym {
        "B" => if aromatic { 2 } else { 3 },
        "C" => if aromatic { 3 } else { 4 },
        "N" => if aromatic { 2 } else { 3 },
        "O" => 2,
        "P" => 3,
        "S" => 2,
        "F" | "Cl" | "Br" | "I" => 1,
        _ => 0,
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Aromatic,
}

#[derive(Clone, Copy, Debug)]
pub struct Bond {
    pub a: usize,
    pub b: usize,
    pub order: BondOrder,
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub symbol: String,
    pub hydrogen: i32,
    pub charge: i32,
    pub aromatic: bool,
    pub in_bracket: bool,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            symbol: String::new(),
            hydrogen: 0,
            charge: 0,
            aromatic: false,
            in_bracket: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub bond_count: i32,
}

impl Molecule {
    /// Return (neighbor_index, bond_order) pairs for the given atom.
    pub fn neighbors(&self, idx: usize) -> Vec<(usize, BondOrder)> {
        let mut result = Vec::new();
        for b in &self.bonds {
            if b.a == idx {
                result.push((b.b, b.order));
            } else if b.b == idx {
                result.push((b.a, b.order));
            }
        }
        result
    }

    /// Return a new Molecule with explicit H atom vertices added.
    /// Each heavy atom's implicit H count is materialized as separate H atoms
    /// connected by single bonds. Used by SMARTS matchers that need to match
    /// [#1] patterns (e.g., Wildman-Crippen H types).
    pub fn with_explicit_hydrogens(&self) -> Molecule {
        let mut atoms = self.atoms.clone();
        let mut bonds = self.bonds.clone();
        let original_count = self.atoms.len();

        for i in 0..original_count {
            let h_count = self.atoms[i].hydrogen.max(0) as usize;
            for _ in 0..h_count {
                let h_idx = atoms.len();
                atoms.push(Atom {
                    symbol: "H".to_string(),
                    hydrogen: 0,
                    charge: 0,
                    aromatic: false,
                    in_bracket: true,
                });
                bonds.push(Bond {
                    a: i,
                    b: h_idx,
                    order: BondOrder::Single,
                });
            }
            // Zero out the implicit H count on the heavy atom since H's are now explicit
            atoms[i].hydrogen = 0;
        }

        let bond_count = bonds.len() as i32;
        Molecule {
            atoms,
            bonds,
            bond_count,
        }
    }

    pub fn element_counts(&self) -> BTreeMap<String, i32> {
        let mut counts = BTreeMap::new();
        for atom in &self.atoms {
            *counts.entry(atom.symbol.clone()).or_insert(0) += 1;
            if atom.hydrogen > 0 {
                *counts.entry("H".to_string()).or_insert(0) += atom.hydrogen;
            }
        }
        counts
    }

    /// Hill system: C first, H second, rest alphabetical
    pub fn formula(&self) -> String {
        let mut counts = self.element_counts();
        if counts.is_empty() {
            return String::new();
        }

        let mut result = String::new();
        let mut append = |elem: &str, count: i32| {
            result.push_str(elem);
            if count > 1 {
                result.push_str(&count.to_string());
            }
        };

        if let Some(c_count) = counts.remove("C") {
            append("C", c_count);
            if let Some(h_count) = counts.remove("H") {
                append("H", h_count);
            }
        }

        for (elem, count) in &counts {
            append(elem, *count);
        }
        result
    }

    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.symbol != "H").count()
    }

    pub fn molecular_weight(&self) -> f64 {
        let counts = self.element_counts();
        let mut weight = 0.0;
        for (elem, count) in &counts {
            if let Some(w) = weights::atomic_weight(elem) {
                weight += w * (*count as f64);
            }
        }
        weight
    }

    pub fn exact_mass(&self) -> f64 {
        let counts = self.element_counts();
        let mut mass = 0.0;
        for (elem, count) in &counts {
            if let Some(m) = weights::monoisotopic_mass(elem) {
                mass += m * (*count as f64);
            }
        }
        mass
    }

    /// Serialize molecule to SMILES with every atom written in bracket form
    /// (verbose form: each atom is `[ELEM]`, `[ELEM+chg]`, or just `[H]`).
    /// Bond orders are shown explicitly (`=`, `#`, `:`); single bonds are implicit.
    /// Disconnected components are joined with `.`.
    /// Output is valid SMILES that round-trips through `parse`.
    pub fn to_smiles_verbose(&self) -> String {
        if self.atoms.is_empty() {
            return String::new();
        }
        let n = self.atoms.len();
        let mut visited = vec![false; n];
        let mut ring_id_for_pair: std::collections::HashMap<(usize, usize), usize> =
            std::collections::HashMap::new();
        let mut next_ring_id: usize = 1;
        let mut output = String::new();
        let mut first = true;
        for start in 0..n {
            if visited[start] {
                continue;
            }
            if !first {
                output.push('.');
            }
            first = false;
            self.dfs_smiles(
                start,
                None,
                &mut visited,
                &mut ring_id_for_pair,
                &mut next_ring_id,
                &mut output,
            );
        }
        output
    }

    fn dfs_smiles(
        &self,
        atom_idx: usize,
        from: Option<usize>,
        visited: &mut [bool],
        ring_id_for_pair: &mut std::collections::HashMap<(usize, usize), usize>,
        next_ring_id: &mut usize,
        output: &mut String,
    ) {
        visited[atom_idx] = true;
        let atom = &self.atoms[atom_idx];
        // Bracket atom: [<sym>(<charge>)]
        output.push('[');
        if atom.aromatic {
            output.push_str(&atom.symbol.to_lowercase());
        } else {
            output.push_str(&atom.symbol);
        }
        if atom.charge > 0 {
            output.push('+');
            if atom.charge > 1 {
                output.push_str(&atom.charge.to_string());
            }
        } else if atom.charge < 0 {
            output.push('-');
            if atom.charge < -1 {
                output.push_str(&(-atom.charge).to_string());
            }
        }
        output.push(']');

        // Collect neighbors except the atom we came from. Snapshot now; visited[] may
        // mutate during recursion below.
        let pending: Vec<(usize, BondOrder)> = self
            .neighbors(atom_idx)
            .into_iter()
            .filter(|(n_idx, _)| Some(*n_idx) != from)
            .collect();

        // Emit ring-closure digits for already-visited neighbors at function entry.
        for (other_idx, order) in pending.iter().filter(|(idx, _)| visited[*idx]) {
            emit_ring_closure(atom_idx, *other_idx, *order, ring_id_for_pair, next_ring_id, output);
        }

        // Iterate the originally-unvisited neighbors in order. A neighbor that becomes
        // visited mid-iteration (via a sibling's DFS) is emitted as a ring closure;
        // a still-unvisited neighbor is recursed into. The LAST still-unvisited
        // neighbor is emitted without parens.
        let unvisited_at_entry: Vec<(usize, BondOrder)> = pending
            .iter()
            .filter(|(idx, _)| !visited[*idx])
            .copied()
            .collect();
        for (i, (next_idx, order)) in unvisited_at_entry.iter().enumerate() {
            if visited[*next_idx] {
                // A sibling DFS already visited this atom → emit a ring closure here.
                emit_ring_closure(
                    atom_idx,
                    *next_idx,
                    *order,
                    ring_id_for_pair,
                    next_ring_id,
                    output,
                );
                continue;
            }
            // Is there any still-unvisited atom AFTER this one in our list?
            let any_after = unvisited_at_entry
                .iter()
                .skip(i + 1)
                .any(|(idx, _)| !visited[*idx]);
            let is_last = !any_after;
            if !is_last {
                output.push('(');
            }
            match order {
                BondOrder::Double => output.push('='),
                BondOrder::Triple => output.push('#'),
                BondOrder::Aromatic => output.push(':'),
                BondOrder::Single => {}
            }
            self.dfs_smiles(
                *next_idx,
                Some(atom_idx),
                visited,
                ring_id_for_pair,
                next_ring_id,
                output,
            );
            if !is_last {
                output.push(')');
            }
        }
    }
}

fn emit_ring_closure(
    a: usize,
    b: usize,
    order: BondOrder,
    ring_id_for_pair: &mut std::collections::HashMap<(usize, usize), usize>,
    next_ring_id: &mut usize,
    output: &mut String,
) {
    let key = if a < b { (a, b) } else { (b, a) };
    let id = *ring_id_for_pair.entry(key).or_insert_with(|| {
        let id = *next_ring_id;
        *next_ring_id += 1;
        id
    });
    // Bond symbol prefix (only non-aromatic non-single needs explicit symbol on closure)
    match order {
        BondOrder::Double => output.push('='),
        BondOrder::Triple => output.push('#'),
        BondOrder::Aromatic => {}
        BondOrder::Single => {}
    }
    if id < 10 {
        output.push((b'0' + id as u8) as char);
    } else {
        output.push('%');
        if id < 100 {
            output.push((b'0' + (id / 10) as u8) as char);
            output.push((b'0' + (id % 10) as u8) as char);
        }
    }
}

/// Parse bracket atom: [NH3+], [Fe+2], [13C@@H], etc.
fn parse_bracket_atom(chars: &[char], pos: &mut usize) -> Option<Atom> {
    let mut atom = Atom { in_bracket: true, ..Default::default() };

    // Skip isotope
    while *pos < chars.len() && chars[*pos].is_ascii_digit() {
        *pos += 1;
    }
    if *pos >= chars.len() {
        return None;
    }

    // Element symbol
    if chars[*pos] == '*' {
        atom.symbol = "*".to_string();
        *pos += 1;
    } else if chars[*pos].is_ascii_lowercase() {
        atom.aromatic = true;
        let mut sym = chars[*pos].to_uppercase().to_string();
        *pos += 1;
        if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
            sym.push(chars[*pos]);
            *pos += 1;
        }
        atom.symbol = sym;
    } else if chars[*pos].is_ascii_uppercase() {
        let mut sym = String::new();
        sym.push(chars[*pos]);
        *pos += 1;
        if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
            sym.push(chars[*pos]);
            *pos += 1;
        }
        atom.symbol = sym;
    } else {
        return None;
    }

    // Skip chirality (@, @@)
    while *pos < chars.len() && chars[*pos] == '@' {
        *pos += 1;
    }

    // Explicit H
    if *pos < chars.len() && chars[*pos] == 'H' {
        *pos += 1;
        if *pos < chars.len() && chars[*pos].is_ascii_digit() {
            atom.hydrogen = (chars[*pos] as i32) - ('0' as i32);
            *pos += 1;
        } else {
            atom.hydrogen = 1;
        }
    }

    // Charge
    if *pos < chars.len() && (chars[*pos] == '+' || chars[*pos] == '-') {
        let sign: i32 = if chars[*pos] == '+' { 1 } else { -1 };
        *pos += 1;
        if *pos < chars.len() && chars[*pos].is_ascii_digit() {
            atom.charge = sign * ((chars[*pos] as i32) - ('0' as i32));
            *pos += 1;
        } else {
            atom.charge = sign;
            let ch = if sign > 0 { '+' } else { '-' };
            while *pos < chars.len() && chars[*pos] == ch {
                atom.charge += sign;
                *pos += 1;
            }
        }
    }

    if *pos >= chars.len() || chars[*pos] != ']' {
        return None;
    }
    *pos += 1;
    Some(atom)
}

fn resolve_bond_order(
    explicit: Option<BondOrder>,
    next_order: i32,
    a_aromatic: bool,
    b_aromatic: bool,
) -> BondOrder {
    if let Some(o) = explicit {
        return o;
    }
    if a_aromatic && b_aromatic {
        return BondOrder::Aromatic;
    }
    match next_order {
        2 => BondOrder::Double,
        3 => BondOrder::Triple,
        _ => BondOrder::Single,
    }
}

/// Parse a SMILES string into a Molecule. Returns None on invalid input.
pub fn parse(smi: &str) -> Option<Molecule> {
    if smi.is_empty() {
        return None;
    }

    let chars: Vec<char> = smi.chars().collect();
    let mut atoms: Vec<Atom> = Vec::new();
    let mut bonds: Vec<Bond> = Vec::new();
    let mut bond_count: i32 = 0;
    let mut degree: Vec<i32> = Vec::new();
    let mut branch_stack: Vec<i32> = Vec::new();
    let mut ring_openings: HashMap<i32, i32> = HashMap::new();
    let mut prev_atom: i32 = -1;
    let mut next_bond_order: i32 = 1;
    let mut explicit_bond: Option<BondOrder> = None;
    let mut pos = 0;

    while pos < chars.len() {
        let c = chars[pos];

        // Branch
        if c == '(' {
            if prev_atom < 0 { return None; }
            branch_stack.push(prev_atom);
            pos += 1;
            continue;
        }
        if c == ')' {
            prev_atom = branch_stack.pop()?;
            next_bond_order = 1;
            pos += 1;
            continue;
        }

        // Bond symbols
        if matches!(c, '=' | '#' | '-' | ':' | '/' | '\\') {
            next_bond_order = match c {
                '=' => 2,
                '#' => 3,
                _ => 1,
            };
            explicit_bond = Some(match c {
                '=' => BondOrder::Double,
                '#' => BondOrder::Triple,
                ':' => BondOrder::Aromatic,
                _ => BondOrder::Single,
            });
            pos += 1;
            continue;
        }

        // Dot (disconnected fragments)
        if c == '.' {
            prev_atom = -1;
            pos += 1;
            continue;
        }

        // Ring closure
        if c.is_ascii_digit() || c == '%' {
            let ring_num;
            if c == '%' {
                pos += 1;
                if pos + 1 >= chars.len() || !chars[pos].is_ascii_digit() || !chars[pos + 1].is_ascii_digit() {
                    return None;
                }
                ring_num = ((chars[pos] as i32) - ('0' as i32)) * 10
                         + ((chars[pos + 1] as i32) - ('0' as i32));
                pos += 2;
            } else {
                ring_num = (c as i32) - ('0' as i32);
                pos += 1;
            }

            if let Some(other) = ring_openings.remove(&ring_num) {
                let a = other as usize;
                let b = prev_atom as usize;
                let order = resolve_bond_order(
                    explicit_bond,
                    next_bond_order,
                    atoms[a].aromatic,
                    atoms[b].aromatic,
                );
                bonds.push(Bond { a, b, order });
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                next_bond_order = 1;
                explicit_bond = None;
            } else {
                ring_openings.insert(ring_num, prev_atom);
            }
            continue;
        }

        // Bracket atom
        if c == '[' {
            pos += 1;
            let atom = parse_bracket_atom(&chars, &mut pos)?;
            let idx = atoms.len() as i32;
            atoms.push(atom);
            degree.push(0);
            if prev_atom >= 0 {
                let a = prev_atom as usize;
                let b = idx as usize;
                let order = resolve_bond_order(
                    explicit_bond,
                    next_bond_order,
                    atoms[a].aromatic,
                    atoms[b].aromatic,
                );
                bonds.push(Bond { a, b, order });
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                next_bond_order = 1;
                explicit_bond = None;
            }
            prev_atom = idx;
            continue;
        }

        // Aromatic atom
        if is_aromatic_char(c) {
            let atom = Atom {
                symbol: c.to_uppercase().to_string(),
                aromatic: true,
                ..Default::default()
            };
            let idx = atoms.len() as i32;
            atoms.push(atom);
            degree.push(0);
            if prev_atom >= 0 {
                let a = prev_atom as usize;
                let b = idx as usize;
                let order = resolve_bond_order(
                    explicit_bond,
                    next_bond_order,
                    atoms[a].aromatic,
                    atoms[b].aromatic,
                );
                bonds.push(Bond { a, b, order });
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                next_bond_order = 1;
                explicit_bond = None;
            }
            prev_atom = idx;
            pos += 1;
            continue;
        }

        // Organic subset atom
        if c.is_ascii_uppercase() {
            let mut sym = String::new();
            sym.push(c);
            if pos + 1 < chars.len() && chars[pos + 1].is_ascii_lowercase() {
                let two = format!("{}{}", c, chars[pos + 1]);
                if is_organic(&two) {
                    sym = two;
                    pos += 1;
                }
            }
            pos += 1;
            if !is_organic(&sym) {
                return None;
            }

            let atom = Atom {
                symbol: sym,
                ..Default::default()
            };
            let idx = atoms.len() as i32;
            atoms.push(atom);
            degree.push(0);
            if prev_atom >= 0 {
                let a = prev_atom as usize;
                let b = idx as usize;
                let order = resolve_bond_order(
                    explicit_bond,
                    next_bond_order,
                    atoms[a].aromatic,
                    atoms[b].aromatic,
                );
                bonds.push(Bond { a, b, order });
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                next_bond_order = 1;
                explicit_bond = None;
            }
            prev_atom = idx;
            continue;
        }

        // Unknown character
        return None;
    }

    if !branch_stack.is_empty() || !ring_openings.is_empty() || atoms.is_empty() {
        return None;
    }

    // Compute implicit H for non-bracket organic subset atoms
    for (i, atom) in atoms.iter_mut().enumerate() {
        if atom.in_bracket {
            continue;
        }
        let val = default_valence(&atom.symbol, atom.aromatic);
        let implicit_h = (val - degree[i]).max(0);
        atom.hydrogen = implicit_h;
    }

    Some(Molecule { atoms, bonds, bond_count })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_water() {
        let mol = parse("O").unwrap();
        assert_eq!(mol.atoms.len(), 1);
        assert_eq!(mol.atoms[0].symbol, "O");
        assert_eq!(mol.atoms[0].hydrogen, 2);
    }

    #[test]
    fn test_parse_ethanol() {
        let mol = parse("CCO").unwrap();
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bond_count, 2);
    }

    #[test]
    fn test_parse_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(mol.atoms.len(), 6);
        assert_eq!(mol.bond_count, 6);
        assert!(mol.atoms[0].aromatic);
    }

    // =================================================================
    // Bond/adjacency tests (added with feature/logP)
    // =================================================================

    #[test]
    fn test_bonds_ethanol() {
        let mol = parse("CCO").unwrap();
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.bonds[0].a, 0);
        assert_eq!(mol.bonds[0].b, 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Single);
        assert_eq!(mol.bonds[1].a, 1);
        assert_eq!(mol.bonds[1].b, 2);
    }

    #[test]
    fn test_bonds_double() {
        let mol = parse("C=O").unwrap();
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Double);
    }

    #[test]
    fn test_bonds_triple() {
        let mol = parse("C#N").unwrap();
        assert_eq!(mol.bonds.len(), 1);
        assert_eq!(mol.bonds[0].order, BondOrder::Triple);
    }

    #[test]
    fn test_bonds_benzene_aromatic() {
        let mol = parse("c1ccccc1").unwrap();
        assert_eq!(mol.bonds.len(), 6);
        for b in &mol.bonds {
            assert_eq!(b.order, BondOrder::Aromatic, "benzene bond should be aromatic");
        }
    }

    #[test]
    fn test_bonds_ring_closure() {
        let mol = parse("C1CCCCC1").unwrap();
        assert_eq!(mol.bonds.len(), 6);
        // Last bond is the ring closure between atom 5 and atom 0
        let last = mol.bonds.last().unwrap();
        assert!((last.a == 0 && last.b == 5) || (last.a == 5 && last.b == 0));
    }

    #[test]
    fn test_neighbors() {
        let mol = parse("CCO").unwrap();
        let n0 = mol.neighbors(0);
        assert_eq!(n0.len(), 1);
        assert_eq!(n0[0].0, 1);
        let n1 = mol.neighbors(1);
        assert_eq!(n1.len(), 2);
    }

    #[test]
    fn test_bonds_branch() {
        // Isobutane CC(C)C: atom 1 has 3 neighbors
        let mol = parse("CC(C)C").unwrap();
        assert_eq!(mol.bonds.len(), 3);
        assert_eq!(mol.neighbors(1).len(), 3);
    }

    #[test]
    fn test_bonds_disconnected() {
        let mol = parse("[Na+].[Cl-]").unwrap();
        assert_eq!(mol.bonds.len(), 0);
    }

    // =================================================================
    // Explicit hydrogen expansion tests (for SMARTS [#1] matching)
    // =================================================================

    #[test]
    fn test_explicit_h_methane() {
        let mol = parse("C").unwrap().with_explicit_hydrogens();
        // C + 4 H atoms
        assert_eq!(mol.atoms.len(), 5);
        assert_eq!(mol.bonds.len(), 4);
        // The original carbon now has hydrogen count 0
        assert_eq!(mol.atoms[0].hydrogen, 0);
        // All other atoms are H
        for i in 1..5 {
            assert_eq!(mol.atoms[i].symbol, "H");
        }
    }

    #[test]
    fn test_explicit_h_water() {
        let mol = parse("O").unwrap().with_explicit_hydrogens();
        assert_eq!(mol.atoms.len(), 3); // O + 2 H
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.atoms[0].symbol, "O");
        assert_eq!(mol.atoms[1].symbol, "H");
        assert_eq!(mol.atoms[2].symbol, "H");
    }

    #[test]
    fn test_explicit_h_ethanol() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        // 3 heavy + (3 + 2 + 1) H = 9 atoms
        assert_eq!(mol.atoms.len(), 9);
        // Heavy bonds (2) + H bonds (6) = 8
        assert_eq!(mol.bonds.len(), 8);
    }

    #[test]
    fn test_explicit_h_neighbors() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        // CH3 (atom 0): should have 1 heavy neighbor (C) + 3 H neighbors = 4
        assert_eq!(mol.neighbors(0).len(), 4);
        // CH2 (atom 1): 2 heavy + 2 H = 4
        assert_eq!(mol.neighbors(1).len(), 4);
        // OH (atom 2): 1 heavy + 1 H = 2
        assert_eq!(mol.neighbors(2).len(), 2);
    }

    #[test]
    fn test_to_smiles_verbose_methane() {
        let s = parse("C").unwrap().with_explicit_hydrogens().to_smiles_verbose();
        // [C]([H])([H])([H])[H]
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 5);
        assert_eq!(parsed.heavy_atom_count(), 1);
    }

    #[test]
    fn test_to_smiles_verbose_ethanol() {
        let s = parse("CCO").unwrap().with_explicit_hydrogens().to_smiles_verbose();
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 9);
        assert_eq!(parsed.heavy_atom_count(), 3);
    }

    #[test]
    fn test_to_smiles_verbose_water() {
        let s = parse("O").unwrap().with_explicit_hydrogens().to_smiles_verbose();
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 3);
        assert_eq!(parsed.heavy_atom_count(), 1);
    }

    #[test]
    fn test_to_smiles_verbose_benzene() {
        let s = parse("c1ccccc1").unwrap().with_explicit_hydrogens().to_smiles_verbose();
        let parsed = parse(&s).expect("round-trip parse");
        // 6 C + 6 H = 12 atoms; round-trip preserves heavy atom count
        assert_eq!(parsed.heavy_atom_count(), 6);
    }

    #[test]
    fn test_to_smiles_verbose_idempotent() {
        // Calling with_explicit_hydrogens twice produces the same result
        let mol1 = parse("CCO").unwrap().with_explicit_hydrogens();
        let mol2 = mol1.clone().with_explicit_hydrogens();
        assert_eq!(mol1.atoms.len(), mol2.atoms.len());
        assert_eq!(mol1.bonds.len(), mol2.bonds.len());
    }

    #[test]
    fn test_to_smiles_verbose_disconnected() {
        let s = parse("[Na+].[Cl-]").unwrap().to_smiles_verbose();
        // Should contain a dot separating the two fragments
        assert!(s.contains('.'));
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 2);
    }

    #[test]
    fn test_parse_charged() {
        let mol = parse("[NH4+]").unwrap();
        assert_eq!(mol.atoms[0].symbol, "N");
        assert_eq!(mol.atoms[0].charge, 1);
        assert_eq!(mol.atoms[0].hydrogen, 4);
    }

    #[test]
    fn test_parse_fragments() {
        let mol = parse("[Na+].[Cl-]").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.bond_count, 0);
    }

    #[test]
    fn test_parse_double_bond() {
        let mol = parse("C=C").unwrap();
        assert_eq!(mol.bond_count, 1);
        assert_eq!(mol.atoms[0].hydrogen, 2); // valence 4 - 2(double) = 2
    }

    #[test]
    fn test_parse_triple_bond() {
        let mol = parse("C#N").unwrap();
        assert_eq!(mol.bond_count, 1);
        assert_eq!(mol.atoms[0].hydrogen, 1); // valence 4 - 3(triple) = 1
    }

    #[test]
    fn test_invalid_empty() {
        assert!(parse("").is_none());
    }

    #[test]
    fn test_invalid_unclosed_branch() {
        assert!(parse("C(C").is_none());
    }

    #[test]
    fn test_invalid_unclosed_ring() {
        assert!(parse("C1CC").is_none());
    }

    #[test]
    fn test_naphthalene() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        assert_eq!(mol.heavy_atom_count(), 10);
    }

    #[test]
    fn test_cyclohexane() {
        let mol = parse("C1CCCCC1").unwrap();
        assert_eq!(mol.heavy_atom_count(), 6);
        assert_eq!(mol.bond_count, 6);
    }

    // =================================================================
    // Bulk molecule coverage tests
    // =================================================================

    /// Helper: parse must succeed, check formula and heavy atom count
    fn check(smiles: &str, expected_formula: &str, expected_heavy: usize) {
        let mol = parse(smiles).unwrap_or_else(|| panic!("Failed to parse: {}", smiles));
        assert_eq!(mol.formula(), expected_formula, "Formula mismatch for {}", smiles);
        assert_eq!(mol.heavy_atom_count(), expected_heavy, "Heavy atom count mismatch for {}", smiles);
    }

    /// Helper: parse must succeed (formula not checked, just parsability + heavy atoms)
    fn check_parse(smiles: &str, expected_heavy: usize) {
        let mol = parse(smiles).unwrap_or_else(|| panic!("Failed to parse: {}", smiles));
        assert_eq!(mol.heavy_atom_count(), expected_heavy, "Heavy atom count mismatch for {}", smiles);
    }

    // --- Simple organic molecules ---

    #[test]
    fn test_methane() { check("C", "CH4", 1); }

    #[test]
    fn test_ethane() { check("CC", "C2H6", 2); }

    #[test]
    fn test_propane() { check("CCC", "C3H8", 3); }

    #[test]
    fn test_butane() { check("CCCC", "C4H10", 4); }

    #[test]
    fn test_isobutane() { check("CC(C)C", "C4H10", 4); }

    #[test]
    fn test_neopentane() { check("CC(C)(C)C", "C5H12", 5); }

    #[test]
    fn test_ethylene() { check("C=C", "C2H4", 2); }

    #[test]
    fn test_acetylene() { check("C#C", "C2H2", 2); }

    #[test]
    fn test_propylene() { check("CC=C", "C3H6", 3); }

    #[test]
    fn test_1_3_butadiene() { check("C=CC=C", "C4H6", 4); }

    // --- Alcohols, aldehydes, ketones, acids ---

    #[test]
    fn test_methanol() { check("CO", "CH4O", 2); }

    #[test]
    fn test_formaldehyde() { check("C=O", "CH2O", 2); }

    #[test]
    fn test_acetone() { check("CC(=O)C", "C3H6O", 4); }

    #[test]
    fn test_acetic_acid() { check("CC(=O)O", "C2H4O2", 4); }

    #[test]
    fn test_formic_acid() { check("O=CO", "CH2O2", 3); }

    #[test]
    fn test_glycerol() { check("OCC(O)CO", "C3H8O3", 6); }

    // --- Amines, amides ---

    #[test]
    fn test_methylamine() { check("CN", "CH5N", 2); }

    #[test]
    fn test_dimethylamine() { check("CNC", "C2H7N", 3); }

    #[test]
    fn test_trimethylamine() { check("CN(C)C", "C3H9N", 4); }

    #[test]
    fn test_urea() { check("NC(=O)N", "CH4N2O", 4); }

    #[test]
    fn test_acetamide() { check("CC(=O)N", "C2H5NO", 4); }

    // --- Halogens ---

    #[test]
    fn test_chloromethane() { check("CCl", "CH3Cl", 2); }

    #[test]
    fn test_bromomethane() { check("CBr", "CH3Br", 2); }

    #[test]
    fn test_iodomethane() { check("CI", "CH3I", 2); }

    #[test]
    fn test_fluoromethane() { check("CF", "CH3F", 2); }

    #[test]
    fn test_dichloromethane() { check("ClCCl", "CH2Cl2", 3); }

    #[test]
    fn test_chloroform() { check("ClC(Cl)Cl", "CHCl3", 4); }

    #[test]
    fn test_carbon_tetrachloride() { check("ClC(Cl)(Cl)Cl", "CCl4", 5); }

    // --- Aromatic compounds ---

    #[test]
    fn test_toluene() { check("Cc1ccccc1", "C7H8", 7); }

    #[test]
    fn test_phenol() { check("Oc1ccccc1", "C6H6O", 7); }

    #[test]
    fn test_aniline() { check("Nc1ccccc1", "C6H7N", 7); }

    #[test]
    fn test_benzoic_acid() { check("OC(=O)c1ccccc1", "C7H6O2", 9); }

    #[test]
    fn test_nitrobenzene() { check("c1ccc([N+](=O)[O-])cc1", "C6H5NO2", 9); }

    #[test]
    fn test_styrene() { check("C=Cc1ccccc1", "C8H8", 8); }

    #[test]
    fn test_biphenyl() { check("c1ccc(-c2ccccc2)cc1", "C12H10", 12); }

    #[test]
    fn test_anthracene() { check("c1ccc2cc3ccccc3cc2c1", "C14H10", 14); }

    // --- Heterocycles ---

    #[test]
    fn test_pyridine() { check("c1ccncc1", "C5H5N", 6); }

    #[test]
    fn test_pyrrole() { check("c1cc[nH]c1", "C4H5N", 5); }

    #[test]
    fn test_furan() { check("c1ccoc1", "C4H4O", 5); }

    #[test]
    fn test_thiophene() { check("c1ccsc1", "C4H4S", 5); }

    #[test]
    fn test_imidazole() { check("c1c[nH]cn1", "C3H4N2", 5); }

    #[test]
    fn test_indole() { check("c1ccc2[nH]ccc2c1", "C8H7N", 9); }

    #[test]
    fn test_quinoline() { check("c1ccc2ncccc2c1", "C9H7N", 10); }

    // --- Rings ---

    #[test]
    fn test_cyclopropane() { check("C1CC1", "C3H6", 3); }

    #[test]
    fn test_cyclobutane() { check("C1CCC1", "C4H8", 4); }

    #[test]
    fn test_cyclopentane() { check("C1CCCC1", "C5H10", 5); }

    #[test]
    fn test_cycloheptane() { check("C1CCCCCC1", "C7H14", 7); }

    #[test]
    fn test_cyclooctane() { check("C1CCCCCCC1", "C8H16", 8); }

    // --- Drugs / bioactive molecules ---

    #[test]
    fn test_aspirin() {
        check("CC(=O)Oc1ccccc1C(=O)O", "C9H8O4", 13);
    }

    #[test]
    fn test_ibuprofen() {
        check("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "C13H18O2", 15);
    }

    #[test]
    fn test_paracetamol() {
        check("CC(=O)Nc1ccc(O)cc1", "C8H9NO2", 11);
    }

    #[test]
    fn test_caffeine() {
        check("Cn1c(=O)c2c(ncn2C)n(C)c1=O", "C8H10N4O2", 14);
    }

    #[test]
    fn test_nicotine() {
        check_parse("CN1CCC[C@@H]1c1cccnc1", 12);
    }

    #[test]
    fn test_dopamine() {
        check("NCCc1ccc(O)c(O)c1", "C8H11NO2", 11);
    }

    #[test]
    fn test_serotonin() {
        check("NCCc1c[nH]c2ccc(O)cc12", "C10H12N2O", 13);
    }

    #[test]
    fn test_adrenaline() {
        check("CNC[C@H](O)c1ccc(O)c(O)c1", "C9H13NO3", 13);
    }

    #[test]
    fn test_penicillin_g_core() {
        check_parse("CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C", 23);
    }

    #[test]
    fn test_cholesterol() {
        check_parse("C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C", 28);
    }

    #[test]
    fn test_glucose() {
        check_parse("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", 12); // 6C + 6O = 12 heavy
    }

    #[test]
    fn test_sucrose() {
        check_parse("OC[C@H]1OC(O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O", 23); // 12C + 11O = 23 heavy
    }

    // --- Amino acids ---

    #[test]
    fn test_glycine() { check("NCC(=O)O", "C2H5NO2", 5); }

    #[test]
    fn test_alanine() { check("C[C@@H](N)C(=O)O", "C3H7NO2", 6); }

    #[test]
    fn test_valine() { check("CC(C)[C@@H](N)C(=O)O", "C5H11NO2", 8); }

    #[test]
    fn test_leucine() { check("CC(C)C[C@@H](N)C(=O)O", "C6H13NO2", 9); }

    #[test]
    fn test_phenylalanine() { check_parse("[C@@H](Cc1ccccc1)(N)C(=O)O", 12); }

    #[test]
    fn test_tryptophan() { check_parse("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", 15); }

    #[test]
    fn test_cysteine() { check("N[C@@H](CS)C(=O)O", "C3H7NO2S", 7); }

    #[test]
    fn test_methionine() { check("CSCC[C@@H](N)C(=O)O", "C5H11NO2S", 9); }

    #[test]
    fn test_proline() { check_parse("OC(=O)[C@@H]1CCCN1", 8); }

    #[test]
    fn test_histidine() { check_parse("N[C@@H](Cc1c[nH]cn1)C(=O)O", 11); }

    // --- Sulfur compounds ---

    #[test]
    fn test_dimethyl_sulfoxide() { check("CS(=O)C", "C2H6OS", 4); }

    #[test]
    fn test_thioacetone() { check("CC(=S)C", "C3H6S", 4); }

    #[test]
    fn test_methanethiol() { check("CS", "CH4S", 2); }

    // --- Phosphorus ---

    #[test]
    fn test_trimethylphosphine() { check("CP(C)C", "C3H9P", 4); }

    // --- Boron ---

    #[test]
    fn test_borane() { check_parse("[BH3]", 1); }

    #[test]
    fn test_phenylboronic_acid() { check("OB(O)c1ccccc1", "C6H7BO2", 9); }

    // --- Charged species ---

    #[test]
    fn test_ammonium() { check_parse("[NH4+]", 1); }

    #[test]
    fn test_hydroxide() { check_parse("[OH-]", 1); }

    #[test]
    fn test_acetate() { check_parse("CC(=O)[O-]", 4); }

    #[test]
    fn test_sodium_chloride() { check_parse("[Na+].[Cl-]", 2); }

    #[test]
    fn test_calcium_chloride() { check_parse("[Ca+2].[Cl-].[Cl-]", 3); }

    #[test]
    fn test_sulfate() { check_parse("[O-]S(=O)(=O)[O-]", 5); }

    // --- Stereo SMILES (chirality/cis-trans) — just check parsability ---

    #[test]
    fn test_cis_2_butene() { check_parse(r"C/C=C\C", 4); }

    #[test]
    fn test_trans_2_butene() { check_parse("C/C=C/C", 4); }

    #[test]
    fn test_l_alanine() { check_parse("[C@@H](N)(C)C(=O)O", 6); }

    #[test]
    fn test_d_alanine() { check_parse("[C@H](N)(C)C(=O)O", 6); }

    // --- Multi-ring / fused systems ---

    #[test]
    fn test_adamantane() { check("C1C2CC3CC1CC(C2)C3", "C10H16", 10); }

    #[test]
    fn test_cubane() { check("C12C3C4C1C5C3C4C25", "C8H8", 8); }

    #[test]
    fn test_decalin() { check("C1CCC2CCCCC2C1", "C10H18", 10); }

    #[test]
    fn test_fluorene() { check("c1ccc2c(c1)Cc1ccccc1-2", "C13H10", 13); }

    // --- Two-digit ring closures ---

    #[test]
    fn test_two_digit_ring() {
        // %10 ring closure
        let mol = parse("C%10CCCCCCCCC%10").unwrap();
        assert_eq!(mol.heavy_atom_count(), 10);
    }

    // --- Edge cases ---

    #[test]
    fn test_single_bracket_atom() { check_parse("[Cu]", 1); }

    #[test]
    fn test_isotope_bracket() { check_parse("[13CH4]", 1); }

    #[test]
    fn test_wildcard_bracket() { check_parse("[*]", 1); }

    #[test]
    fn test_deep_branch() {
        check("C(C(C(C(C)C)C)C)C", "C9H20", 9);
    }

    #[test]
    fn test_long_chain() {
        // C20 alkane
        let smiles = "C".repeat(20);
        let mol = parse(&smiles).unwrap();
        assert_eq!(mol.heavy_atom_count(), 20);
    }

    // =================================================================
    // Nucleobases
    // =================================================================

    #[test]
    fn test_adenine() { check("c1nc(N)c2nc[nH]c2n1", "C5H5N5", 10); }

    #[test]
    fn test_guanine() { check("c1nc2c(n1)[nH]c(=O)n2N", "C4H4N5O", 10); }

    #[test]
    fn test_cytosine() { check("c1cnc(=O)[nH]c1N", "C4H5N3O", 8); }

    #[test]
    fn test_thymine() { check("Cc1c[nH]c(=O)[nH]c1=O", "C5H6N2O2", 9); }

    #[test]
    fn test_uracil() { check("c1c[nH]c(=O)[nH]c1=O", "C4H4N2O2", 8); }

    // =================================================================
    // Steroids (beyond cholesterol)
    // =================================================================

    #[test]
    fn test_estradiol() {
        check_parse("C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O", 20);
    }

    #[test]
    fn test_progesterone() {
        check_parse("CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@@]34C)C", 23);
    }

    #[test]
    fn test_cortisol() {
        check_parse("C[C@@]12C[C@@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO", 26);
    }

    #[test]
    fn test_testosterone() {
        check_parse("C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@@H]2O", 21);
    }

    // =================================================================
    // Anticancer drugs
    // =================================================================

    #[test]
    fn test_5_fluorouracil() {
        check("Fc1c[nH]c(=O)[nH]c1=O", "C4H3FN2O2", 9);
    }

    #[test]
    fn test_methotrexate() {
        check_parse("CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1", 33);
    }

    #[test]
    fn test_cisplatin() {
        // [Pt] in bracket
        check_parse("[NH3][Pt]([NH3])(Cl)Cl", 5);
    }

    #[test]
    fn test_doxorubicin_core() {
        check_parse("O=C1c2cccc(O)c2C(=O)c2c(O)cc(O)cc21", 19);
    }

    #[test]
    fn test_tamoxifen() {
        check_parse("CC/C(=C(\\c1ccccc1)c1ccccc1)c1ccc(OCCN(C)C)cc1", 28);
    }

    // =================================================================
    // Vitamins
    // =================================================================

    #[test]
    fn test_ascorbic_acid() {
        // Vitamin C
        check_parse("OC[C@H](O)[C@H]1OC(=O)C(O)=C1O", 12);
    }

    #[test]
    fn test_retinol() {
        // Vitamin A
        check_parse("CC1=C(/C=C/C(C)=C/C=C/C(C)=C/CO)C(C)(C)CCC1", 21);
    }

    #[test]
    fn test_pyridoxine() {
        // Vitamin B6
        check("Cc1ncc(CO)c(CO)c1O", "C8H11NO3", 12);
    }

    #[test]
    fn test_niacin() {
        // Vitamin B3 / nicotinic acid
        check("OC(=O)c1cccnc1", "C6H5NO2", 9);
    }

    #[test]
    fn test_riboflavin_core() {
        check_parse("Cc1cc2nc3c(=O)[nH]c(=O)nc3n(C)c2cc1C", 19);
    }

    // =================================================================
    // Sugars (monosaccharides)
    // =================================================================

    #[test]
    fn test_fructose() {
        check_parse("OC[C@H]1OC(O)(CO)[C@@H](O)[C@@H]1O", 12);
    }

    #[test]
    fn test_galactose() {
        check_parse("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@H]1O", 12);
    }

    #[test]
    fn test_ribose() {
        check_parse("OC[C@H]1OC(O)[C@H](O)[C@@H]1O", 10);
    }

    #[test]
    fn test_mannose() {
        check_parse("OC[C@H]1OC(O)[C@@H](O)[C@@H](O)[C@@H]1O", 12);
    }

    #[test]
    fn test_deoxyribose() {
        check_parse("OC[C@H]1OC(O)C[C@@H]1O", 9);
    }

    // =================================================================
    // Fatty acids
    // =================================================================

    #[test]
    fn test_palmitic_acid() {
        // C16:0
        check("CCCCCCCCCCCCCCCC(=O)O", "C16H32O2", 18);
    }

    #[test]
    fn test_stearic_acid() {
        // C18:0
        check("CCCCCCCCCCCCCCCCCC(=O)O", "C18H36O2", 20);
    }

    #[test]
    fn test_oleic_acid() {
        // C18:1 cis-9
        check_parse(r"CCCCCCCC/C=C\CCCCCCCC(=O)O", 20);
    }

    #[test]
    fn test_linoleic_acid() {
        // C18:2
        check_parse(r"CCCCCC=CCC=CCCCCCCCC(=O)O", 20);
    }

    #[test]
    fn test_arachidonic_acid() {
        // C20:4
        check_parse(r"CCCCCC=CCC=CCC=CCC=CCCCC(=O)O", 22);
    }

    #[test]
    fn test_dha() {
        // DHA C22:6
        check_parse("CC=CCC=CCC=CCC=CCC=CCC=CCCC(=O)O", 23);
    }

    // =================================================================
    // Organometallics / coordination compounds
    // =================================================================

    #[test]
    fn test_ferrocene_simplified() {
        // Simplified — two Cp fragments + Fe
        check_parse("[Fe+2].[c-]1cccc1.[c-]1cccc1", 11);
    }

    #[test]
    fn test_hemoglobin_fe() {
        // Just Fe bracket
        check_parse("[Fe+2]", 1);
    }

    #[test]
    fn test_zinc_acetate() {
        check_parse("[Zn+2].CC([O-])=O.CC([O-])=O", 9);
    }

    #[test]
    fn test_silver_ion() {
        check_parse("[Ag+]", 1);
    }

    #[test]
    fn test_mercury_chloride() {
        check_parse("Cl[Hg]Cl", 3);
    }

    #[test]
    fn test_copper_sulfate() {
        check_parse("[Cu+2].[O-]S([O-])(=O)=O", 6);
    }

    // =================================================================
    // Macrocycles / crown ethers / porphyrins
    // =================================================================

    #[test]
    fn test_12_crown_4() {
        check("C1COCCOCCOCCO1", "C8H16O4", 12);
    }

    #[test]
    fn test_18_crown_6() {
        check("C1COCCOCCOCCOCCOCCO1", "C12H24O6", 18);
    }

    #[test]
    fn test_cyclooctadecane() {
        // Large ring, 18-membered
        let smiles = "C1CCCCCCCCCCCCCCCCC1";
        let mol = parse(smiles).unwrap();
        assert_eq!(mol.heavy_atom_count(), 18);
    }

    #[test]
    fn test_porphyrin_core() {
        // Porphine — basic porphyrin macrocycle
        check_parse("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1[nH]2)[nH]5)n4)[nH]3", 24);
    }

    #[test]
    fn test_cyclodextrin_fragment() {
        // Maltose-like fragment (2 glucose units)
        check_parse("OC[C@H]1OC(O[C@@H]2[C@@H](O)[C@H](O)[C@@H](O)OC2CO)[C@H](O)[C@@H](O)[C@@H]1O", 23);
    }

    // =================================================================
    // Polymers / fragments
    // =================================================================

    #[test]
    fn test_ethylene_glycol() {
        check("OCCO", "C2H6O2", 4);
    }

    #[test]
    fn test_peg_trimer() {
        check("OCCOCCOCCO", "C6H14O4", 10);
    }

    #[test]
    fn test_lactic_acid() {
        check("C[C@@H](O)C(=O)O", "C3H6O3", 6);
    }

    #[test]
    fn test_caprolactam() {
        check("O=C1CCCCCN1", "C6H11NO", 8);
    }

    #[test]
    fn test_styrene_monomer() {
        check("C=Cc1ccccc1", "C8H8", 8);
    }

    #[test]
    fn test_vinyl_chloride() {
        check("C=CCl", "C2H3Cl", 3);
    }

    #[test]
    fn test_acrylonitrile() {
        check("C=CC#N", "C3H3N", 4);
    }

    // =================================================================
    // Neurotransmitters / signaling molecules
    // =================================================================

    #[test]
    fn test_gaba() {
        check("NCCCC(=O)O", "C4H9NO2", 7);
    }

    #[test]
    fn test_glutamate() {
        check("N[C@@H](CCC(=O)O)C(=O)O", "C5H9NO4", 10);
    }

    #[test]
    fn test_acetylcholine() {
        check_parse("CC(=O)OCC[N+](C)(C)C", 10);
    }

    #[test]
    fn test_histamine() {
        check("NCCc1c[nH]cn1", "C5H9N3", 8);
    }

    #[test]
    fn test_melatonin() {
        check("CC(=O)NCCc1c[nH]c2ccc(OC)cc12", "C13H16N2O2", 17);
    }

    // =================================================================
    // Antibiotics
    // =================================================================

    #[test]
    fn test_sulfanilamide() {
        check("Nc1ccc(cc1)S(=O)(=O)N", "C6H8N2O2S", 11);
    }

    #[test]
    fn test_chloramphenicol() {
        check_parse("O[C@@H](C(=O)NCc1ccc([N+](=O)[O-])cc1)[C@H](O)CO", 19);
    }

    #[test]
    fn test_trimethoprim() {
        check_parse("COc1cc(Cc2cnc(N)nc2N)cc(OC)c1OC", 21);
    }

    // =================================================================
    // Painkillers / NSAIDs
    // =================================================================

    #[test]
    fn test_naproxen() {
        check_parse("COc1ccc2cc(ccc2c1)[C@@H](C)C(=O)O", 17);
    }

    #[test]
    fn test_diclofenac() {
        check_parse("OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl", 19);
    }

    #[test]
    fn test_morphine() {
        check_parse("CN1CC[C@]23c4c5ccc(O)c4O[C@H]2C(=C[C@@H]1[C@@H]3O)C5", 20);
    }

    #[test]
    fn test_codeine() {
        check_parse("COc1ccc2C[C@H]3N(C)CC[C@@]45c2c1O[C@H]4[C@@H](O)C=C[C@@H]35", 22);
    }

    // =================================================================
    // Antidepressants / psychoactive
    // =================================================================

    #[test]
    fn test_fluoxetine() {
        // Prozac
        check_parse("CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1", 22);
    }

    #[test]
    fn test_diazepam() {
        check_parse("CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21", 20);
    }

    #[test]
    fn test_lsd_core() {
        // Ergoline core simplified
        check_parse("CCN(CC)C(=O)[C@H]1CN(C)[C@@H]2Cc3c[nH]c4cccc(C2=C1)c34", 24);
    }

    // =================================================================
    // Dyes / pigments
    // =================================================================

    #[test]
    fn test_indigo() {
        check_parse("O=C1/C(=C2\\Nc3ccccc3C2=O)Nc2ccccc21", 20);
    }

    #[test]
    fn test_methyl_orange_core() {
        // Azobenzene core
        check_parse("c1ccc(/N=N/c2ccccc2)cc1", 14);
    }

    #[test]
    fn test_fluorescein_core() {
        check_parse("OC(=O)c1ccccc1-c1c2ccc(=O)cc-2oc2cc(O)ccc12", 25);
    }

    // =================================================================
    // Explosives / energetic materials
    // =================================================================

    #[test]
    fn test_tnt() {
        check_parse("Cc1c(cc([N+](=O)[O-])cc1[N+](=O)[O-])[N+](=O)[O-]", 16);
    }

    #[test]
    fn test_nitroglycerin() {
        check("O[N+](=O)OCC(O[N+](=O)O)CO[N+](=O)O", "C3H8N3O9", 15);
    }

    // =================================================================
    // Pesticides / herbicides
    // =================================================================

    #[test]
    fn test_ddt() {
        check_parse("ClC(Cl)C(c1ccc(Cl)cc1)c1ccc(Cl)cc1", 18);
    }

    #[test]
    fn test_glyphosate() {
        check("OC(=O)CNCP(=O)(O)O", "C3H8NO5P", 10);
    }

    // =================================================================
    // Solvents / industrial chemicals
    // =================================================================

    #[test]
    fn test_dmf() {
        check("CN(C)C=O", "C3H7NO", 5);
    }

    #[test]
    fn test_dmso() {
        check("CS(=O)C", "C2H6OS", 4);
    }

    #[test]
    fn test_thf() {
        check("C1CCOC1", "C4H8O", 5);
    }

    #[test]
    fn test_dioxane() {
        check("C1COCCO1", "C4H8O2", 6);
    }

    #[test]
    fn test_acetonitrile() {
        check("CC#N", "C2H3N", 3);
    }

    #[test]
    fn test_dmac() {
        check("CN(C)C(=O)C", "C4H9NO", 6);
    }

    #[test]
    fn test_nmp() {
        check("CN1CCCC1=O", "C5H9NO", 7);
    }

    // =================================================================
    // Remaining amino acids (complete the 20)
    // =================================================================

    #[test]
    fn test_isoleucine() { check("CC[C@H](C)[C@@H](N)C(=O)O", "C6H13NO2", 9); }

    #[test]
    fn test_serine() { check("N[C@@H](CO)C(=O)O", "C3H7NO3", 7); }

    #[test]
    fn test_threonine() { check("C[C@@H](O)[C@@H](N)C(=O)O", "C4H9NO3", 8); }

    #[test]
    fn test_aspartate() { check("N[C@@H](CC(=O)O)C(=O)O", "C4H7NO4", 9); }

    #[test]
    fn test_asparagine() { check("N[C@@H](CC(=O)N)C(=O)O", "C4H8N2O3", 9); }

    #[test]
    fn test_glutamine() { check("N[C@@H](CCC(=O)N)C(=O)O", "C5H10N2O3", 10); }

    #[test]
    fn test_lysine() { check("NCCCC[C@@H](N)C(=O)O", "C6H14N2O2", 10); }

    #[test]
    fn test_arginine() { check_parse("N[C@@H](CCCNC(=N)N)C(=O)O", 12); }

    #[test]
    fn test_tyrosine() { check_parse("N[C@@H](Cc1ccc(O)cc1)C(=O)O", 13); }

    // =================================================================
    // Edge cases: very large / complex
    // =================================================================

    #[test]
    fn test_c60_fragment() {
        // Corannulene (bowl-shaped PAH, C60 fragment)
        check_parse("c1cc2ccc3ccc4ccc5ccc1c1c2c3c4c51", 20);
    }

    #[test]
    fn test_c30_chain() {
        let smiles = "C".repeat(30);
        let mol = parse(&smiles).unwrap();
        assert_eq!(mol.heavy_atom_count(), 30);
    }

    #[test]
    fn test_c50_chain() {
        let smiles = "C".repeat(50);
        let mol = parse(&smiles).unwrap();
        assert_eq!(mol.heavy_atom_count(), 50);
    }

    #[test]
    fn test_many_branches() {
        // Star-shaped: central C with 4 chains
        check("C(CCCC)(CCCC)(CCCC)CCCC", "C17H36", 17);
    }

    #[test]
    fn test_spiro_compound() {
        // Spiro[4.5]decane
        check("C1CCC2(CC1)CCCCC2", "C11H20", 11);
    }

    #[test]
    fn test_bridged_bicyclic() {
        // Norbornane (bicyclo[2.2.1]heptane)
        check("C1CC2CC1CC2", "C7H12", 7);
    }

    #[test]
    fn test_multiple_stereocenters() {
        // 4 stereocenters + 2 terminal C
        check_parse("[C@@H](O)([C@H](O)[C@@H](O)[C@H](O)C)C", 10);
    }
}
