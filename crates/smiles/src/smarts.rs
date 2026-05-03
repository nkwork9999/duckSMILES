//! Minimal SMARTS parser and matcher.
//!
//! Supports the subset needed by Wildman-Crippen LogP (RDKit's Crippen.txt):
//! - Element symbols (C, Cl, etc.) and aromatic (c, n, ...)
//! - `[#n]` atomic number, `[Hn]`, `[Xn]`, `[+n]`, `[-n]`
//! - Operators: `;` (AND_LO), `,` (OR), implicit AND_HI, `!` (NOT)
//! - `A` / `a` (aliphatic / aromatic wildcards), `*` (any)
//! - Atom lists via `,`: `[N,O,P]`
//! - Bonds: `-`, `=`, `#`, `:`, `~`, default = single-or-aromatic
//! - Branches `(...)` — tree-shaped patterns only (no ring closures).

use crate::parser::{BondOrder, Molecule};

// =============================================================================
// Types
// =============================================================================

#[derive(Clone, Debug)]
pub enum Primitive {
    Any,                          // *
    AnyAliphatic,                 // A
    AnyAromatic,                  // a
    AliphaticElement(String),     // uppercase or [C]
    AromaticElement(String),      // lowercase or [c]
    AtomicNum(u8),                // #<n>
    TotalH(u8),                   // H or H<n>
    Connections(u8),              // X<n>
    PositiveCharge(u8),           // + or +<n>
    NegativeCharge(u8),           // - or -<n>
}

#[derive(Clone, Debug)]
pub enum Predicate {
    Prim(Primitive),
    And(Vec<Predicate>),
    Or(Vec<Predicate>),
    Not(Box<Predicate>),
}

#[derive(Clone, Debug)]
pub struct AtomPattern {
    pub pred: Predicate,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondSpec {
    Default,   // single OR aromatic (no explicit symbol)
    Single,    // -
    Double,    // =
    Triple,    // #
    Aromatic,  // :
    Any,       // ~
}

#[derive(Clone, Debug)]
pub struct BondEdge {
    pub a: usize,
    pub b: usize,
    pub spec: BondSpec,
}

#[derive(Clone, Debug)]
pub struct Pattern {
    pub atoms: Vec<AtomPattern>,
    pub bonds: Vec<BondEdge>,
}

// =============================================================================
// Parser
// =============================================================================

pub fn parse_smarts(s: &str) -> Option<Pattern> {
    if s.is_empty() {
        return None;
    }
    let chars: Vec<char> = s.chars().collect();
    let mut atoms: Vec<AtomPattern> = Vec::new();
    let mut bonds: Vec<BondEdge> = Vec::new();
    let mut branch_stack: Vec<i32> = Vec::new();
    let mut prev_atom: i32 = -1;
    let mut pending_bond: Option<BondSpec> = None;
    let mut pos = 0;

    while pos < chars.len() {
        let c = chars[pos];

        if c == '(' {
            if prev_atom < 0 {
                return None;
            }
            branch_stack.push(prev_atom);
            pos += 1;
            continue;
        }
        if c == ')' {
            prev_atom = branch_stack.pop()?;
            pos += 1;
            continue;
        }

        // Bond symbols
        if matches!(c, '-' | '=' | '#' | ':' | '~' | '/' | '\\') {
            pending_bond = Some(match c {
                '-' => BondSpec::Single,
                '=' => BondSpec::Double,
                '#' => BondSpec::Triple,
                ':' => BondSpec::Aromatic,
                '~' => BondSpec::Any,
                _ => BondSpec::Single, // / \ — directional, treat as single
            });
            pos += 1;
            continue;
        }

        // Bracket atom
        if c == '[' {
            pos += 1;
            let pred = parse_bracket(&chars, &mut pos)?;
            let idx = add_atom(&mut atoms, &mut bonds, prev_atom, &mut pending_bond, pred);
            prev_atom = idx as i32;
            continue;
        }

        // Wildcards outside brackets
        if c == '*' {
            let idx = add_atom(
                &mut atoms,
                &mut bonds,
                prev_atom,
                &mut pending_bond,
                Predicate::Prim(Primitive::Any),
            );
            prev_atom = idx as i32;
            pos += 1;
            continue;
        }
        if c == 'A' {
            let idx = add_atom(
                &mut atoms,
                &mut bonds,
                prev_atom,
                &mut pending_bond,
                Predicate::Prim(Primitive::AnyAliphatic),
            );
            prev_atom = idx as i32;
            pos += 1;
            continue;
        }
        if c == 'a' {
            let idx = add_atom(
                &mut atoms,
                &mut bonds,
                prev_atom,
                &mut pending_bond,
                Predicate::Prim(Primitive::AnyAromatic),
            );
            prev_atom = idx as i32;
            pos += 1;
            continue;
        }

        // Element (uppercase or lowercase) outside brackets
        if c.is_ascii_alphabetic() {
            let aromatic = c.is_ascii_lowercase();
            let mut sym = c.to_ascii_uppercase().to_string();
            pos += 1;
            if !aromatic && pos < chars.len() && chars[pos].is_ascii_lowercase() {
                let combined = format!("{}{}", sym, chars[pos]);
                if is_known_element(&combined) {
                    sym = combined;
                    pos += 1;
                }
            }
            let prim = if aromatic {
                Primitive::AromaticElement(sym)
            } else {
                Primitive::AliphaticElement(sym)
            };
            let idx = add_atom(
                &mut atoms,
                &mut bonds,
                prev_atom,
                &mut pending_bond,
                Predicate::Prim(prim),
            );
            prev_atom = idx as i32;
            continue;
        }

        return None; // unknown character
    }

    if !branch_stack.is_empty() || atoms.is_empty() {
        return None;
    }

    Some(Pattern { atoms, bonds })
}

fn add_atom(
    atoms: &mut Vec<AtomPattern>,
    bonds: &mut Vec<BondEdge>,
    prev_atom: i32,
    pending_bond: &mut Option<BondSpec>,
    pred: Predicate,
) -> usize {
    let idx = atoms.len();
    atoms.push(AtomPattern { pred });
    if prev_atom >= 0 {
        bonds.push(BondEdge {
            a: prev_atom as usize,
            b: idx,
            spec: pending_bond.take().unwrap_or(BondSpec::Default),
        });
    }
    idx
}

fn is_known_element(s: &str) -> bool {
    // Two-letter organic-subset-like elements used in Crippen.txt
    matches!(
        s,
        "He" | "Li" | "Be" | "Ne" | "Na" | "Mg" | "Al" | "Si" | "Ar" | "Ca"
            | "Sc" | "Ti" | "Cr" | "Mn" | "Fe" | "Co" | "Ni" | "Cu" | "Zn"
            | "Ga" | "Ge" | "As" | "Se" | "Br" | "Kr" | "Rb" | "Sr" | "Zr"
            | "Nb" | "Mo" | "Tc" | "Ru" | "Rh" | "Pd" | "Ag" | "Cd" | "In"
            | "Sn" | "Sb" | "Te" | "Xe" | "Cs" | "Ba" | "Hf" | "Ta" | "Re"
            | "Os" | "Ir" | "Pt" | "Au" | "Hg" | "Tl" | "Pb" | "Bi" | "Po"
            | "At" | "Cl"
    )
}

// Bracket content grammar (precedence low → high):
//   and_lo := or_expr (';' or_expr)*
//   or_expr := and_hi (',' and_hi)*
//   and_hi := unary+                    (implicit concatenation)
//   unary := '!' unary | primitive
fn parse_bracket(chars: &[char], pos: &mut usize) -> Option<Predicate> {
    let pred = parse_and_lo(chars, pos)?;
    if *pos >= chars.len() || chars[*pos] != ']' {
        return None;
    }
    *pos += 1;
    Some(pred)
}

fn parse_and_lo(chars: &[char], pos: &mut usize) -> Option<Predicate> {
    let mut items = vec![parse_or(chars, pos)?];
    while *pos < chars.len() && chars[*pos] == ';' {
        *pos += 1;
        items.push(parse_or(chars, pos)?);
    }
    Some(if items.len() == 1 {
        items.pop().unwrap()
    } else {
        Predicate::And(items)
    })
}

fn parse_or(chars: &[char], pos: &mut usize) -> Option<Predicate> {
    let mut items = vec![parse_and_hi(chars, pos)?];
    while *pos < chars.len() && chars[*pos] == ',' {
        *pos += 1;
        items.push(parse_and_hi(chars, pos)?);
    }
    Some(if items.len() == 1 {
        items.pop().unwrap()
    } else {
        Predicate::Or(items)
    })
}

fn parse_and_hi(chars: &[char], pos: &mut usize) -> Option<Predicate> {
    let mut items = vec![parse_unary(chars, pos)?];
    loop {
        if *pos >= chars.len() {
            break;
        }
        let c = chars[*pos];
        if matches!(c, ';' | ',' | ']' | '&') {
            break;
        }
        items.push(parse_unary(chars, pos)?);
    }
    Some(if items.len() == 1 {
        items.pop().unwrap()
    } else {
        Predicate::And(items)
    })
}

fn parse_unary(chars: &[char], pos: &mut usize) -> Option<Predicate> {
    if *pos >= chars.len() {
        return None;
    }
    if chars[*pos] == '!' {
        *pos += 1;
        Some(Predicate::Not(Box::new(parse_unary(chars, pos)?)))
    } else {
        Some(Predicate::Prim(parse_primitive(chars, pos)?))
    }
}

fn parse_primitive(chars: &[char], pos: &mut usize) -> Option<Primitive> {
    if *pos >= chars.len() {
        return None;
    }
    let c = chars[*pos];

    // Skip leading isotope digits (ignored — Crippen doesn't care)
    if c.is_ascii_digit() {
        while *pos < chars.len() && chars[*pos].is_ascii_digit() {
            *pos += 1;
        }
        return parse_primitive(chars, pos);
    }

    match c {
        '*' => {
            *pos += 1;
            Some(Primitive::Any)
        }
        'A' => {
            *pos += 1;
            if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
                let sym = format!("A{}", chars[*pos]);
                if is_known_element(&sym) {
                    *pos += 1;
                    return Some(Primitive::AliphaticElement(sym));
                }
            }
            Some(Primitive::AnyAliphatic)
        }
        'a' => {
            *pos += 1;
            Some(Primitive::AnyAromatic)
        }
        '#' => {
            *pos += 1;
            let n = read_number(chars, pos)?;
            Some(Primitive::AtomicNum(n as u8))
        }
        'H' => {
            *pos += 1;
            let n = if *pos < chars.len() && chars[*pos].is_ascii_digit() {
                let d = chars[*pos] as u8 - b'0';
                *pos += 1;
                d
            } else {
                1
            };
            Some(Primitive::TotalH(n))
        }
        'X' => {
            *pos += 1;
            let d = chars.get(*pos)?.to_digit(10)? as u8;
            *pos += 1;
            Some(Primitive::Connections(d))
        }
        '+' => {
            *pos += 1;
            let n = if *pos < chars.len() && chars[*pos].is_ascii_digit() {
                let d = chars[*pos] as u8 - b'0';
                *pos += 1;
                d
            } else {
                1
            };
            Some(Primitive::PositiveCharge(n))
        }
        '-' => {
            *pos += 1;
            let n = if *pos < chars.len() && chars[*pos].is_ascii_digit() {
                let d = chars[*pos] as u8 - b'0';
                *pos += 1;
                d
            } else {
                1
            };
            Some(Primitive::NegativeCharge(n))
        }
        c if c.is_ascii_uppercase() => {
            let mut sym = String::from(c);
            *pos += 1;
            if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
                let combined = format!("{}{}", sym, chars[*pos]);
                if is_known_element(&combined) {
                    sym = combined;
                    *pos += 1;
                }
            }
            Some(Primitive::AliphaticElement(sym))
        }
        c if c.is_ascii_lowercase() => {
            *pos += 1;
            Some(Primitive::AromaticElement(c.to_ascii_uppercase().to_string()))
        }
        _ => None,
    }
}

fn read_number(chars: &[char], pos: &mut usize) -> Option<u32> {
    if *pos >= chars.len() || !chars[*pos].is_ascii_digit() {
        return None;
    }
    let mut n: u32 = 0;
    while *pos < chars.len() && chars[*pos].is_ascii_digit() {
        n = n * 10 + (chars[*pos] as u32 - '0' as u32);
        *pos += 1;
    }
    Some(n)
}

// =============================================================================
// Matcher
// =============================================================================

/// Try to match `pat` with `pat.atoms[0]` anchored at `start` in `mol`.
pub fn match_at(pat: &Pattern, mol: &Molecule, start: usize) -> bool {
    if pat.atoms.is_empty() || start >= mol.atoms.len() {
        return false;
    }
    if !atom_matches(&pat.atoms[0].pred, mol, start) {
        return false;
    }
    if pat.atoms.len() == 1 {
        return true;
    }

    let mut mapping = vec![usize::MAX; pat.atoms.len()];
    let mut used = vec![false; mol.atoms.len()];
    mapping[0] = start;
    used[start] = true;
    match_recurse(pat, mol, 1, &mut mapping, &mut used)
}

fn match_recurse(
    pat: &Pattern,
    mol: &Molecule,
    next: usize,
    mapping: &mut [usize],
    used: &mut [bool],
) -> bool {
    if next >= pat.atoms.len() {
        return true;
    }

    // Find a bond in the pattern that connects `next` to an already-mapped atom
    let found = pat.bonds.iter().find_map(|b| {
        if b.a == next && mapping[b.b] != usize::MAX {
            Some((b.b, b.spec))
        } else if b.b == next && mapping[b.a] != usize::MAX {
            Some((b.a, b.spec))
        } else {
            None
        }
    });
    let (pred_pat_idx, bond_spec) = match found {
        Some(v) => v,
        None => return false,
    };

    let anchor = mapping[pred_pat_idx];

    for (nbr_idx, bond_order) in mol.neighbors(anchor) {
        if used[nbr_idx] {
            continue;
        }
        if !bond_matches(bond_spec, bond_order) {
            continue;
        }
        if !atom_matches(&pat.atoms[next].pred, mol, nbr_idx) {
            continue;
        }
        mapping[next] = nbr_idx;
        used[nbr_idx] = true;
        if match_recurse(pat, mol, next + 1, mapping, used) {
            return true;
        }
        mapping[next] = usize::MAX;
        used[nbr_idx] = false;
    }
    false
}

fn bond_matches(spec: BondSpec, order: BondOrder) -> bool {
    match spec {
        BondSpec::Single => order == BondOrder::Single,
        BondSpec::Double => order == BondOrder::Double,
        BondSpec::Triple => order == BondOrder::Triple,
        BondSpec::Aromatic => order == BondOrder::Aromatic,
        BondSpec::Default => matches!(order, BondOrder::Single | BondOrder::Aromatic),
        BondSpec::Any => true,
    }
}

fn atom_matches(pred: &Predicate, mol: &Molecule, idx: usize) -> bool {
    match pred {
        Predicate::Prim(p) => prim_matches(p, mol, idx),
        Predicate::And(items) => items.iter().all(|p| atom_matches(p, mol, idx)),
        Predicate::Or(items) => items.iter().any(|p| atom_matches(p, mol, idx)),
        Predicate::Not(inner) => !atom_matches(inner, mol, idx),
    }
}

fn prim_matches(prim: &Primitive, mol: &Molecule, idx: usize) -> bool {
    let atom = &mol.atoms[idx];
    match prim {
        Primitive::Any => true,
        Primitive::AnyAliphatic => !atom.aromatic,
        Primitive::AnyAromatic => atom.aromatic,
        Primitive::AliphaticElement(sym) => !atom.aromatic && &atom.symbol == sym,
        Primitive::AromaticElement(sym) => atom.aromatic && &atom.symbol == sym,
        Primitive::AtomicNum(n) => atomic_num(&atom.symbol) == *n,
        Primitive::TotalH(n) => total_h(mol, idx) == *n as i32,
        Primitive::Connections(n) => total_connections(mol, idx) == *n as usize,
        Primitive::PositiveCharge(n) => atom.charge == *n as i32,
        Primitive::NegativeCharge(n) => atom.charge == -(*n as i32),
    }
}

/// Sum of implicit H count + explicit H neighbors (for SMARTS H<n> matching).
fn total_h(mol: &Molecule, idx: usize) -> i32 {
    let implicit = mol.atoms[idx].hydrogen.max(0);
    let explicit: i32 = mol
        .neighbors(idx)
        .iter()
        .filter(|(nbr, _)| mol.atoms[*nbr].symbol == "H")
        .count() as i32;
    implicit + explicit
}

/// Total connection count (heavy neighbors + implicit H, or all neighbors if H expanded).
fn total_connections(mol: &Molecule, idx: usize) -> usize {
    let direct_neighbors = mol.neighbors(idx).len();
    let implicit_h = mol.atoms[idx].hydrogen.max(0) as usize;
    direct_neighbors + implicit_h
}

fn atomic_num(sym: &str) -> u8 {
    match sym {
        "H" => 1, "He" => 2, "Li" => 3, "Be" => 4, "B" => 5, "C" => 6, "N" => 7,
        "O" => 8, "F" => 9, "Ne" => 10, "Na" => 11, "Mg" => 12, "Al" => 13, "Si" => 14,
        "P" => 15, "S" => 16, "Cl" => 17, "Ar" => 18, "K" => 19, "Ca" => 20,
        "Sc" => 21, "Ti" => 22, "V" => 23, "Cr" => 24, "Mn" => 25, "Fe" => 26,
        "Co" => 27, "Ni" => 28, "Cu" => 29, "Zn" => 30, "Ga" => 31, "Ge" => 32,
        "As" => 33, "Se" => 34, "Br" => 35, "Kr" => 36, "Rb" => 37, "Sr" => 38,
        "Y" => 39, "Zr" => 40, "Nb" => 41, "Mo" => 42, "Tc" => 43, "Ru" => 44,
        "Rh" => 45, "Pd" => 46, "Ag" => 47, "Cd" => 48, "In" => 49, "Sn" => 50,
        "Sb" => 51, "Te" => 52, "I" => 53, "Xe" => 54, "Cs" => 55, "Ba" => 56,
        "Hf" => 72, "Ta" => 73, "W" => 74, "Re" => 75, "Os" => 76, "Ir" => 77,
        "Pt" => 78, "Au" => 79, "Hg" => 80, "Tl" => 81, "Pb" => 82, "Bi" => 83,
        "Po" => 84, "At" => 85,
        _ => 0,
    }
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn matches_anywhere(smarts: &str, smiles: &str) -> bool {
        let pat = parse_smarts(smarts).expect("parse smarts");
        let mol = parse(smiles).expect("parse smiles");
        (0..mol.atoms.len()).any(|i| match_at(&pat, &mol, i))
    }

    fn count_matches(smarts: &str, smiles: &str) -> usize {
        let pat = parse_smarts(smarts).expect("parse smarts");
        let mol = parse(smiles).expect("parse smiles");
        (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count()
    }

    // --------- Parser basics ---------

    #[test]
    fn parse_simple() {
        let p = parse_smarts("C").unwrap();
        assert_eq!(p.atoms.len(), 1);
        assert_eq!(p.bonds.len(), 0);
    }

    #[test]
    fn parse_two_atoms_default_bond() {
        let p = parse_smarts("CC").unwrap();
        assert_eq!(p.atoms.len(), 2);
        assert_eq!(p.bonds.len(), 1);
        assert_eq!(p.bonds[0].spec, BondSpec::Default);
    }

    #[test]
    fn parse_explicit_bonds() {
        assert_eq!(parse_smarts("C=O").unwrap().bonds[0].spec, BondSpec::Double);
        assert_eq!(parse_smarts("C#N").unwrap().bonds[0].spec, BondSpec::Triple);
        assert_eq!(parse_smarts("c:c").unwrap().bonds[0].spec, BondSpec::Aromatic);
        assert_eq!(parse_smarts("C-C").unwrap().bonds[0].spec, BondSpec::Single);
    }

    #[test]
    fn parse_branch() {
        let p = parse_smarts("C(O)N").unwrap();
        assert_eq!(p.atoms.len(), 3);
        assert_eq!(p.bonds.len(), 2);
    }

    #[test]
    fn parse_bracket_elements() {
        assert!(parse_smarts("[C]").is_some());
        assert!(parse_smarts("[c]").is_some());
        assert!(parse_smarts("[#6]").is_some());
        assert!(parse_smarts("[Cl]").is_some());
        assert!(parse_smarts("[Fe]").is_some());
    }

    #[test]
    fn parse_bracket_combos() {
        assert!(parse_smarts("[CH3]").is_some());
        assert!(parse_smarts("[CH2X4]").is_some());
        assert!(parse_smarts("[NH+0]").is_some());
        assert!(parse_smarts("[O-]").is_some());
    }

    #[test]
    fn parse_crippen_patterns() {
        // Representative Crippen.txt patterns
        let samples = [
            "[CH4]",
            "[CH2](C)C",
            "[CH3][N,O,P,S,F,Cl,Br,I]",
            "[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]",
            "[c](:a)(:a):a",
            "[c](:a)(:a)-C",
            "[NH2+0][A;!#1]",
            "[NH+0]([A;!#1])[A;!#1]",
            "[NH3,NH2,NH;+,+2,+3]",
            "[O]=[#7,#8]",
            "[#9,#17,#35,#53;-]",
            "[S;-,-2,-3,-4,+1,+2,+3,+5,+6]",
            "[!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]",
            "[#1][#6,#1]",
        ];
        for s in samples {
            assert!(parse_smarts(s).is_some(), "failed to parse: {}", s);
        }
    }

    #[test]
    fn parse_bad_inputs() {
        assert!(parse_smarts("").is_none());
        assert!(parse_smarts("[").is_none());
        assert!(parse_smarts("[C").is_none());
        assert!(parse_smarts("C(").is_none());
    }

    // --------- Element matching ---------

    #[test]
    fn match_aliphatic_element() {
        assert!(matches_anywhere("C", "CCO"));
        assert!(matches_anywhere("O", "CCO"));
        assert!(!matches_anywhere("N", "CCO"));
    }

    #[test]
    fn match_aromatic_vs_aliphatic() {
        assert!(matches_anywhere("c", "c1ccccc1"));
        assert!(!matches_anywhere("c", "CCO"));
        assert!(!matches_anywhere("C", "c1ccccc1"));
    }

    #[test]
    fn match_atomic_num() {
        assert!(matches_anywhere("[#6]", "CCO"));
        assert!(matches_anywhere("[#8]", "CCO"));
        assert!(!matches_anywhere("[#7]", "CCO"));
        // #6 matches both aliphatic and aromatic C
        assert!(matches_anywhere("[#6]", "c1ccccc1"));
    }

    #[test]
    fn match_any_wildcards() {
        assert!(matches_anywhere("*", "CCO"));
        assert!(matches_anywhere("[A]", "CCO"));
        assert!(matches_anywhere("[a]", "c1ccccc1"));
        assert!(!matches_anywhere("[a]", "CCO"));
    }

    // --------- Hydrogens & connections ---------

    #[test]
    fn match_total_h() {
        // CCO: atoms = CH3, CH2, OH
        assert_eq!(count_matches("[CH3]", "CCO"), 1);
        assert_eq!(count_matches("[CH2]", "CCO"), 1);
        assert_eq!(count_matches("[H1]", "CCO"), 1); // the O
        assert_eq!(count_matches("[OH]", "CCO"), 1);
    }

    #[test]
    fn match_connections() {
        // CCO: CH3 (X4), CH2 (X4), OH (X2)
        assert_eq!(count_matches("[CX4]", "CCO"), 2);
        assert_eq!(count_matches("[OX2]", "CCO"), 1);
        // CC#N: atom 0 C = 1 heavy + 3 H = X4, atom 1 C = 2 heavy + 0 H = X2
        assert_eq!(count_matches("[CX2]", "CC#N"), 1);
    }

    // --------- Charges ---------

    #[test]
    fn match_charges() {
        assert!(matches_anywhere("[O-]", "CC(=O)[O-]"));
        assert!(matches_anywhere("[NH4+]", "[NH4+]"));
        assert!(matches_anywhere("[N+0]", "CCN")); // neutral N
        assert!(!matches_anywhere("[N+]", "CCN"));
    }

    // --------- Logical operators ---------

    #[test]
    fn match_not() {
        assert!(matches_anywhere("[!#1]", "CCO"));
        assert!(!matches_anywhere("[!C;!N;!O;!#1]", "CCO"));
    }

    #[test]
    fn match_or_list() {
        assert!(matches_anywhere("[N,O,P]", "CCO"));
        assert!(!matches_anywhere("[N,P,S]", "CCO"));
    }

    #[test]
    fn match_and_lo() {
        // A AND NOT hydrogen = non-H aliphatic
        assert!(matches_anywhere("[A;!#1]", "CCO"));
    }

    #[test]
    fn match_nh_charge_combo() {
        // [NH3,NH2,NH;+,+2,+3]: one of NH1/NH2/NH3 and positive charge.
        // [NH4+] has H=4 → should NOT match.
        assert!(!matches_anywhere("[NH3,NH2,NH;+,+2,+3]", "[NH4+]"));
        // [NH3+] has H=3 → matches.
        assert!(matches_anywhere("[NH3,NH2,NH;+,+2,+3]", "[NH3+]"));
    }

    // --------- Bond matching ---------

    #[test]
    fn match_double_bond() {
        assert!(matches_anywhere("C=O", "CC(=O)O"));
        assert!(!matches_anywhere("C=O", "CCO"));
    }

    #[test]
    fn match_triple_bond() {
        assert!(matches_anywhere("C#N", "CC#N"));
        assert!(!matches_anywhere("C#N", "CCN"));
    }

    #[test]
    fn match_aromatic_bond() {
        assert!(matches_anywhere("c:c", "c1ccccc1"));
    }

    // --------- Crippen-shaped patterns ---------

    #[test]
    fn crippen_c1_methane() {
        // C1: [CH4]
        let pat = parse_smarts("[CH4]").unwrap();
        let mol = parse("C").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_c1_ethane() {
        // C1: [CH3]C
        let pat = parse_smarts("[CH3]C").unwrap();
        let mol = parse("CC").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_c3_methanol() {
        // C3: [CH3][N,O,P,S,F,Cl,Br,I]
        let pat = parse_smarts("[CH3][N,O,P,S,F,Cl,Br,I]").unwrap();
        let mol = parse("CO").unwrap();
        assert!(match_at(&pat, &mol, 0)); // C has H3 and O neighbor
        assert!(!match_at(&pat, &mol, 1)); // O root would not match [CH3]
    }

    #[test]
    fn crippen_c18_benzene() {
        // C18: [cH] — aromatic C with 1 H
        assert_eq!(count_matches("[cH]", "c1ccccc1"), 6);
    }

    #[test]
    fn crippen_c19_fused_aromatic() {
        // C19: [c](:a)(:a):a — aromatic C bonded to 3 aromatic via aromatic bonds
        // Naphthalene has 2 ring-junction carbons matching this
        let count = count_matches("[c](:a)(:a):a", "c1ccc2ccccc2c1");
        assert_eq!(count, 2, "expected 2 ring-junction aromatic C");
    }

    #[test]
    fn crippen_o5_nitro_oxygen() {
        // O5: [O]=[#7,#8]
        let pat = parse_smarts("[O]=[#7,#8]").unwrap();
        let mol = parse("O=N").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_fallback_cs() {
        // CS: [#6] — any C, fallback
        let pat = parse_smarts("[#6]").unwrap();
        let mol = parse("CCO").unwrap();
        assert!(match_at(&pat, &mol, 0));
        assert!(match_at(&pat, &mol, 1));
        assert!(!match_at(&pat, &mol, 2)); // O
    }

    #[test]
    fn crippen_halogen_list() {
        // F / Cl / Br / I = [#9-0] / [#17-0] / [#35-0] / [#53-0]
        assert!(matches_anywhere("[#9-0]", "CF"));
        assert!(matches_anywhere("[#17-0]", "CCl"));
        assert!(matches_anywhere("[#35-0]", "CBr"));
    }

    #[test]
    fn crippen_o2_hydroxyl() {
        // O2: [OH,OH2]
        assert!(matches_anywhere("[OH,OH2]", "CCO")); // ethanol OH
        assert!(matches_anywhere("[OH,OH2]", "O"));   // water
    }

    // --------- Multi-atom traversal ---------

    #[test]
    fn match_three_atom_chain() {
        // [O]=C-N in urea-like fragment
        let pat = parse_smarts("O=CN").unwrap();
        let mol = parse("NC(=O)N").unwrap();
        // The O should match; traverse double bond to C, then single bond to N
        assert!(mol.atoms.iter().enumerate().any(|(i, _)| match_at(&pat, &mol, i)));
    }

    #[test]
    fn match_branch() {
        // C(=O)O pattern on acetic acid CC(=O)O
        let pat = parse_smarts("C(=O)O").unwrap();
        let mol = parse("CC(=O)O").unwrap();
        // The carboxyl C (atom 1) matches
        assert!(match_at(&pat, &mol, 1));
    }

    // --------- H-expanded patterns (Crippen H types) ---------

    #[test]
    fn match_hydrogen_with_explicit_h() {
        // [#1] should match explicit H atoms after expansion
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        // CCO has 3+2+1 = 6 H atoms
        assert_eq!(count, 6);
    }

    #[test]
    fn match_ch3_after_expansion() {
        // Critical: [CH3] must still match C atoms after H expansion
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[CH3]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 1); // only atom 0 (methyl C)
    }

    #[test]
    fn match_connections_after_expansion() {
        // [CX4] should still correctly identify sp3 carbons
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[CX4]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 2); // both C's
    }

    #[test]
    fn match_crippen_h1_pattern() {
        // H1: [#1][#6,#1] — H bonded to C or H
        // In CCO, all H's are bonded to C (atoms 0,1) or O (atom 2).
        // The 5 H's on C atoms match H1; the 1 H on O doesn't.
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1][#6,#1]").unwrap();
        let matches: Vec<usize> = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .collect();
        assert_eq!(matches.len(), 5);
    }

    #[test]
    fn match_crippen_h2_pattern() {
        // H2 (partial): [#1]O[CX4,c] — H bonded to O bonded to sp3 C or aromatic c
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1]O[CX4,c]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 1); // the single OH hydrogen
    }
}
