//! Minimal SMARTS parser and matcher.
//!
//! Originally written for the subset needed by Wildman-Crippen LogP
//! (RDKit's Crippen.txt). Extended for MACCS keys, which additionally need:
//! - Ring-closure bonds (`C1CCCCC1`, `*1~*~*~*~1`)
//! - Ring-membership atom queries (`[R]`, `[R0]`, `[#16R]`, `[!#6R]`)
//! - Ring / non-ring bond specs (`@`, `!@`) and not-aromatic bond (`!:`)
//! - Recursive SMARTS environments (`[$(...)]`)
//! - Uniquified match counting (RDKit `SubstructMatch(..., uniquify=true)`)
//!
//! Supported atom features:
//! - Element symbols (C, Cl, etc.) and aromatic (c, n, ...)
//! - `[#n]` atomic number, `[Hn]`, `[Xn]`, `[+n]`, `[-n]`, `[R]`/`[R0]`
//! - Operators: `;` (AND_LO), `,` (OR), implicit AND_HI, `!` (NOT)
//! - `A` / `a` (aliphatic / aromatic wildcards), `*` (any)
//! - Atom lists via `,`: `[N,O,P]`
//! - Recursive environments `$(...)`
//!
//! Supported bond features:
//! - `-`, `=`, `#`, `:`, `~`, `@`, `!@`, `!:`, default = single-or-aromatic
//! - Ring closures and branches `(...)`.

use crate::parser::{BondOrder, Molecule, RingInfo};

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
    InRing(bool),                 // R (true) or R0 (false)
    Recursive(Box<Pattern>),      // $(...)
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

/// A single bond constraint. Bond queries in MACCS are simple enough to model
/// as one base order constraint plus optional ring/aromaticity modifiers,
/// each possibly negated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BondSpec {
    Default,    // single OR aromatic (no explicit symbol)
    Single,     // -
    Double,     // =
    Triple,     // #
    Aromatic,   // :
    NotAromatic, // !:
    Any,        // ~
    InRing,     // @
    NotInRing,  // !@
}

/// A bond constraint is one or more `BondSpec`s that must ALL hold. This models
/// combinations like `=@` (double AND in-ring, MACCS bit 26). A lone `Default`
/// (empty meaning) is represented as `specs == [Default]`.
#[derive(Clone, Debug)]
pub struct BondEdge {
    pub a: usize,
    pub b: usize,
    pub specs: Vec<BondSpec>,
}

impl BondEdge {
    fn new(a: usize, b: usize, specs: Vec<BondSpec>) -> Self {
        let specs = if specs.is_empty() {
            vec![BondSpec::Default]
        } else {
            specs
        };
        BondEdge { a, b, specs }
    }
}

#[derive(Clone, Debug)]
pub struct Pattern {
    pub atoms: Vec<AtomPattern>,
    pub bonds: Vec<BondEdge>,
    /// True if any bond spec or atom query depends on ring membership, so the
    /// matcher must compute `RingInfo` for the molecule. Patterns without ring
    /// dependence (e.g. all of Crippen.txt) skip that work entirely.
    pub needs_ring_info: bool,
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
    // ring-closure bookkeeping: digit -> (atom_idx, pending bond specs at open)
    let mut ring_open: std::collections::HashMap<u32, (usize, Vec<BondSpec>)> =
        std::collections::HashMap::new();
    let mut prev_atom: i32 = -1;
    // Accumulated bond constraints for the next bond (ALL must hold), e.g. `=@`.
    let mut pending_bond: Vec<BondSpec> = Vec::new();
    let mut needs_ring_info = false;
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

        // Bond symbols (including `!@` / `!:`). A leading `!` before a bond
        // symbol is a negated bond constraint. Consecutive bond symbols
        // accumulate as an AND (e.g. `=@` = double AND in-ring).
        if c == '!' && pos + 1 < chars.len() && matches!(chars[pos + 1], '@' | ':') {
            pending_bond.push(match chars[pos + 1] {
                '@' => BondSpec::NotInRing,
                ':' => BondSpec::NotAromatic,
                _ => unreachable!(),
            });
            needs_ring_info = true;
            pos += 2;
            continue;
        }
        if matches!(c, '-' | '=' | '#' | ':' | '~' | '/' | '\\' | '@') {
            let spec = match c {
                '-' => BondSpec::Single,
                '=' => BondSpec::Double,
                '#' => BondSpec::Triple,
                ':' => BondSpec::Aromatic,
                '~' => BondSpec::Any,
                '@' => BondSpec::InRing,
                _ => BondSpec::Single, // / \ — directional, treat as single
            };
            if spec == BondSpec::InRing {
                needs_ring_info = true;
            }
            pending_bond.push(spec);
            pos += 1;
            continue;
        }

        // Ring-closure digit(s): `1`..`9` or `%nn`.
        if c.is_ascii_digit() || c == '%' {
            if prev_atom < 0 {
                return None;
            }
            let ring_num = if c == '%' {
                pos += 1;
                let d1 = chars.get(pos)?.to_digit(10)?;
                pos += 1;
                let d2 = chars.get(pos)?.to_digit(10)?;
                pos += 1;
                d1 * 10 + d2
            } else {
                pos += 1;
                c.to_digit(10).unwrap()
            };
            let specs = std::mem::take(&mut pending_bond);
            match ring_open.remove(&ring_num) {
                Some((other, open_specs)) => {
                    // Closing the ring: use whichever side carried an explicit
                    // bond spec (SMILES allows the spec on either side).
                    let final_specs = if !specs.is_empty() { specs } else { open_specs };
                    if final_specs
                        .iter()
                        .any(|s| matches!(s, BondSpec::InRing | BondSpec::NotInRing))
                    {
                        needs_ring_info = true;
                    }
                    bonds.push(BondEdge::new(other, prev_atom as usize, final_specs));
                }
                None => {
                    ring_open.insert(ring_num, (prev_atom as usize, specs));
                }
            }
            continue;
        }

        // Bracket atom
        if c == '[' {
            pos += 1;
            let (pred, ring_dep) = parse_bracket(&chars, &mut pos)?;
            if ring_dep {
                needs_ring_info = true;
            }
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

    if !branch_stack.is_empty() || !ring_open.is_empty() || atoms.is_empty() {
        return None;
    }

    Some(Pattern {
        atoms,
        bonds,
        needs_ring_info,
    })
}

fn add_atom(
    atoms: &mut Vec<AtomPattern>,
    bonds: &mut Vec<BondEdge>,
    prev_atom: i32,
    pending_bond: &mut Vec<BondSpec>,
    pred: Predicate,
) -> usize {
    let idx = atoms.len();
    atoms.push(AtomPattern { pred });
    if prev_atom >= 0 {
        let specs = std::mem::take(pending_bond);
        bonds.push(BondEdge::new(prev_atom as usize, idx, specs));
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
//
// Returns (predicate, ring_dependent) where ring_dependent is true if any
// primitive inside depends on ring membership.
fn parse_bracket(chars: &[char], pos: &mut usize) -> Option<(Predicate, bool)> {
    let (pred, ring_dep) = parse_and_lo(chars, pos)?;
    if *pos >= chars.len() || chars[*pos] != ']' {
        return None;
    }
    *pos += 1;
    Some((pred, ring_dep))
}

fn parse_and_lo(chars: &[char], pos: &mut usize) -> Option<(Predicate, bool)> {
    let (first, mut rd) = parse_or(chars, pos)?;
    let mut items = vec![first];
    while *pos < chars.len() && chars[*pos] == ';' {
        *pos += 1;
        let (p, r) = parse_or(chars, pos)?;
        rd |= r;
        items.push(p);
    }
    Some((
        if items.len() == 1 {
            items.pop().unwrap()
        } else {
            Predicate::And(items)
        },
        rd,
    ))
}

fn parse_or(chars: &[char], pos: &mut usize) -> Option<(Predicate, bool)> {
    let (first, mut rd) = parse_and_hi(chars, pos)?;
    let mut items = vec![first];
    while *pos < chars.len() && chars[*pos] == ',' {
        *pos += 1;
        let (p, r) = parse_and_hi(chars, pos)?;
        rd |= r;
        items.push(p);
    }
    Some((
        if items.len() == 1 {
            items.pop().unwrap()
        } else {
            Predicate::Or(items)
        },
        rd,
    ))
}

fn parse_and_hi(chars: &[char], pos: &mut usize) -> Option<(Predicate, bool)> {
    let (first, mut rd) = parse_unary(chars, pos)?;
    let mut items = vec![first];
    loop {
        if *pos >= chars.len() {
            break;
        }
        let c = chars[*pos];
        if matches!(c, ';' | ',' | ']') {
            break;
        }
        if c == '&' {
            // explicit high-precedence AND — same as implicit concatenation
            *pos += 1;
            continue;
        }
        let (p, r) = parse_unary(chars, pos)?;
        rd |= r;
        items.push(p);
    }
    Some((
        if items.len() == 1 {
            items.pop().unwrap()
        } else {
            Predicate::And(items)
        },
        rd,
    ))
}

fn parse_unary(chars: &[char], pos: &mut usize) -> Option<(Predicate, bool)> {
    if *pos >= chars.len() {
        return None;
    }
    if chars[*pos] == '!' {
        *pos += 1;
        let (inner, rd) = parse_unary(chars, pos)?;
        Some((Predicate::Not(Box::new(inner)), rd))
    } else {
        let (prim, rd) = parse_primitive(chars, pos)?;
        Some((Predicate::Prim(prim), rd))
    }
}

fn parse_primitive(chars: &[char], pos: &mut usize) -> Option<(Primitive, bool)> {
    if *pos >= chars.len() {
        return None;
    }
    let c = chars[*pos];

    // Skip leading isotope digits (ignored — Crippen/MACCS don't care)
    if c.is_ascii_digit() {
        while *pos < chars.len() && chars[*pos].is_ascii_digit() {
            *pos += 1;
        }
        return parse_primitive(chars, pos);
    }

    match c {
        '$' => {
            // Recursive SMARTS: $(...)
            *pos += 1;
            if *pos >= chars.len() || chars[*pos] != '(' {
                return None;
            }
            *pos += 1;
            // Capture balanced parentheses
            let start = *pos;
            let mut depth = 1;
            while *pos < chars.len() && depth > 0 {
                match chars[*pos] {
                    '(' => depth += 1,
                    ')' => depth -= 1,
                    _ => {}
                }
                if depth == 0 {
                    break;
                }
                *pos += 1;
            }
            if depth != 0 {
                return None;
            }
            let sub: String = chars[start..*pos].iter().collect();
            *pos += 1; // consume ')'
            let pat = parse_smarts(&sub)?;
            let rd = pat.needs_ring_info;
            Some((Primitive::Recursive(Box::new(pat)), rd))
        }
        '*' => {
            *pos += 1;
            Some((Primitive::Any, false))
        }
        'R' => {
            *pos += 1;
            // R0 means "not in ring"; R or R<n>0... we only need R / R0.
            if *pos < chars.len() && chars[*pos] == '0' {
                *pos += 1;
                Some((Primitive::InRing(false), true))
            } else {
                // optional ring-count digits (treat any R<n>, n>0, as "in ring")
                while *pos < chars.len() && chars[*pos].is_ascii_digit() {
                    *pos += 1;
                }
                Some((Primitive::InRing(true), true))
            }
        }
        'A' => {
            *pos += 1;
            if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
                let sym = format!("A{}", chars[*pos]);
                if is_known_element(&sym) {
                    *pos += 1;
                    return Some((Primitive::AliphaticElement(sym), false));
                }
            }
            Some((Primitive::AnyAliphatic, false))
        }
        'a' => {
            *pos += 1;
            Some((Primitive::AnyAromatic, false))
        }
        '#' => {
            *pos += 1;
            let n = read_number(chars, pos)?;
            Some((Primitive::AtomicNum(n as u8), false))
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
            Some((Primitive::TotalH(n), false))
        }
        'X' => {
            *pos += 1;
            let d = chars.get(*pos)?.to_digit(10)? as u8;
            *pos += 1;
            Some((Primitive::Connections(d), false))
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
            Some((Primitive::PositiveCharge(n), false))
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
            Some((Primitive::NegativeCharge(n), false))
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
            Some((Primitive::AliphaticElement(sym), false))
        }
        c if c.is_ascii_lowercase() => {
            *pos += 1;
            Some((
                Primitive::AromaticElement(c.to_ascii_uppercase().to_string()),
                false,
            ))
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

/// Per-molecule ring data, computed lazily and cached by raw pointer identity.
/// MACCS matching reuses one molecule across ~150 patterns, so we avoid
/// recomputing ring perception for every pattern.
struct MatchCtx<'a> {
    mol: &'a Molecule,
    ring: Option<RingInfo>,
    adj: Vec<Vec<(usize, usize)>>,
}

impl<'a> MatchCtx<'a> {
    fn new(mol: &'a Molecule, needs_ring: bool) -> Self {
        MatchCtx {
            mol,
            ring: if needs_ring { Some(mol.ring_info()) } else { None },
            adj: mol.adjacency(),
        }
    }

    /// Bond index between two adjacent atoms (None if not bonded).
    fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.adj[a].iter().find(|(n, _)| *n == b).map(|(_, bi)| *bi)
    }
}

/// Try to match `pat` with `pat.atoms[0]` anchored at `start` in `mol`.
/// Preserves the original semantics: returns true iff a subgraph isomorphism
/// of the pattern exists with pattern-atom 0 mapped to `start`.
pub fn match_at(pat: &Pattern, mol: &Molecule, start: usize) -> bool {
    let ctx = MatchCtx::new(mol, pat.needs_ring_info);
    match_at_ctx(pat, &ctx, start)
}

fn match_at_ctx(pat: &Pattern, ctx: &MatchCtx, start: usize) -> bool {
    if pat.atoms.is_empty() || start >= ctx.mol.atoms.len() {
        return false;
    }
    if !atom_matches(&pat.atoms[0].pred, ctx, start) {
        return false;
    }
    if pat.atoms.len() == 1 {
        return true;
    }

    let mut mapping = vec![usize::MAX; pat.atoms.len()];
    let mut used = vec![false; ctx.mol.atoms.len()];
    mapping[0] = start;
    used[start] = true;
    let mut found = false;
    // emit returns false to stop at the first complete match.
    match_recurse(pat, ctx, 1, &mut mapping, &mut used, &mut |_| {
        found = true;
        false
    });
    found
}

/// Enumerate all matches of `pat` anchored at `start`, invoking `emit` with the
/// atom mapping for each. `emit` returns false to stop early.
fn for_each_match_at(
    pat: &Pattern,
    ctx: &MatchCtx,
    start: usize,
    emit: &mut dyn FnMut(&[usize]) -> bool,
) {
    if pat.atoms.is_empty() || start >= ctx.mol.atoms.len() {
        return;
    }
    if !atom_matches(&pat.atoms[0].pred, ctx, start) {
        return;
    }
    if pat.atoms.len() == 1 {
        emit(&[start]);
        return;
    }
    let mut mapping = vec![usize::MAX; pat.atoms.len()];
    let mut used = vec![false; ctx.mol.atoms.len()];
    mapping[0] = start;
    used[start] = true;
    match_recurse(pat, ctx, 1, &mut mapping, &mut used, emit);
}

/// Recursive VF-style extension. Calls `emit(mapping)` on every complete match;
/// returns false from the whole recursion once `emit` asks to stop.
fn match_recurse(
    pat: &Pattern,
    ctx: &MatchCtx,
    next: usize,
    mapping: &mut [usize],
    used: &mut [bool],
    emit: &mut dyn FnMut(&[usize]) -> bool,
) -> bool {
    if next >= pat.atoms.len() {
        // Complete mapping — verify ALL pattern bonds (including ring closures)
        // are satisfied, then emit.
        for b in &pat.bonds {
            let ma = mapping[b.a];
            let mb = mapping[b.b];
            match ctx.bond_between(ma, mb) {
                Some(bi) => {
                    let order = ctx.mol.bonds[bi].order;
                    if !bond_specs_match(&b.specs, order, ctx, bi) {
                        return true; // continue search (this mapping invalid)
                    }
                }
                None => return true, // not bonded → invalid mapping, keep searching
            }
        }
        return emit(mapping);
    }

    // Find a pattern bond connecting `next` to an already-mapped atom.
    let found = pat.bonds.iter().find_map(|b| {
        if b.a == next && mapping[b.b] != usize::MAX {
            Some((b.b, &b.specs))
        } else if b.b == next && mapping[b.a] != usize::MAX {
            Some((b.a, &b.specs))
        } else {
            None
        }
    });
    let (pred_pat_idx, bond_specs) = match found {
        Some(v) => v,
        None => return true, // disconnected pattern atom — should not happen for our SMARTS
    };

    let anchor = mapping[pred_pat_idx];

    for &(nbr_idx, bi) in &ctx.adj[anchor] {
        if used[nbr_idx] {
            continue;
        }
        let order = ctx.mol.bonds[bi].order;
        if !bond_specs_match(bond_specs, order, ctx, bi) {
            continue;
        }
        if !atom_matches(&pat.atoms[next].pred, ctx, nbr_idx) {
            continue;
        }
        mapping[next] = nbr_idx;
        used[nbr_idx] = true;
        if !match_recurse(pat, ctx, next + 1, mapping, used, emit) {
            return false; // emit asked to stop
        }
        mapping[next] = usize::MAX;
        used[nbr_idx] = false;
    }
    true
}

/// All specs must hold (AND), e.g. `=@` = Double AND InRing.
fn bond_specs_match(specs: &[BondSpec], order: BondOrder, ctx: &MatchCtx, bond_idx: usize) -> bool {
    specs.iter().all(|&s| bond_matches(s, order, ctx, bond_idx))
}

fn bond_matches(spec: BondSpec, order: BondOrder, ctx: &MatchCtx, bond_idx: usize) -> bool {
    let in_ring = ctx
        .ring
        .as_ref()
        .map(|r| r.bond_in_ring[bond_idx])
        .unwrap_or(false);
    match spec {
        BondSpec::Single => order == BondOrder::Single,
        BondSpec::Double => order == BondOrder::Double,
        BondSpec::Triple => order == BondOrder::Triple,
        BondSpec::Aromatic => order == BondOrder::Aromatic,
        BondSpec::NotAromatic => order != BondOrder::Aromatic,
        BondSpec::Default => matches!(order, BondOrder::Single | BondOrder::Aromatic),
        BondSpec::Any => true,
        BondSpec::InRing => in_ring,
        BondSpec::NotInRing => !in_ring,
    }
}

fn atom_matches(pred: &Predicate, ctx: &MatchCtx, idx: usize) -> bool {
    match pred {
        Predicate::Prim(p) => prim_matches(p, ctx, idx),
        Predicate::And(items) => items.iter().all(|p| atom_matches(p, ctx, idx)),
        Predicate::Or(items) => items.iter().any(|p| atom_matches(p, ctx, idx)),
        Predicate::Not(inner) => !atom_matches(inner, ctx, idx),
    }
}

fn prim_matches(prim: &Primitive, ctx: &MatchCtx, idx: usize) -> bool {
    let mol = ctx.mol;
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
        Primitive::InRing(want) => {
            let in_ring = ctx
                .ring
                .as_ref()
                .map(|r| r.atom_in_ring[idx])
                .unwrap_or(false);
            in_ring == *want
        }
        Primitive::Recursive(sub) => {
            // The recursive environment must match with its atom 0 anchored at idx.
            match_at_ctx(sub, ctx, idx)
        }
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

// =============================================================================
// Public match helpers (presence + uniquified counting)
// =============================================================================

/// True if the pattern matches anywhere in the molecule.
pub fn matches_mol(pat: &Pattern, mol: &Molecule) -> bool {
    let ctx = MatchCtx::new(mol, pat.needs_ring_info);
    (0..mol.atoms.len()).any(|i| match_at_ctx(pat, &ctx, i))
}

/// Count of unique matches, mirroring RDKit `SubstructMatch(..., uniquify=true)`.
///
/// RDKit's `uniquify` deduplicates matches whose *set* of matched atoms is
/// identical (order-independent), which is what the MACCS count thresholds
/// were tuned against. We enumerate every anchored embedding and dedupe by the
/// sorted set of mapped molecule atoms.
pub fn count_unique(pat: &Pattern, mol: &Molecule) -> usize {
    unique_matches(pat, mol).len()
}

/// Return unique substructure matches as pattern-atom ordered molecule atom
/// mappings. Atom indices are 0-based internally; SQL JSON output converts
/// them to 1-based indices.
pub fn unique_matches(pat: &Pattern, mol: &Molecule) -> Vec<Vec<usize>> {
    let ctx = MatchCtx::new(mol, pat.needs_ring_info);
    let mut seen: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
    let mut matches = Vec::new();
    for start in 0..mol.atoms.len() {
        for_each_match_at(pat, &ctx, start, &mut |mapping| {
            let mut key: Vec<usize> = mapping.to_vec();
            key.sort_unstable();
            if seen.insert(key) {
                matches.push(mapping.to_vec());
            }
            true
        });
    }
    matches.sort();
    matches
}

pub fn unique_matches_json(pat: &Pattern, mol: &Molecule) -> String {
    let matches = unique_matches(pat, mol);
    let mut out = String::from("[");
    for (i, mapping) in matches.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"match\":");
        out.push_str(&(i + 1).to_string());
        out.push_str(",\"atom_indices\":");
        push_usize_array_1based(&mut out, mapping);
        out.push_str(",\"atoms\":[");
        for (query_idx, &target_idx) in mapping.iter().enumerate() {
            if query_idx > 0 {
                out.push(',');
            }
            out.push_str("{\"query_atom\":");
            out.push_str(&(query_idx + 1).to_string());
            out.push_str(",\"target_atom\":");
            out.push_str(&(target_idx + 1).to_string());
            out.push_str(",\"symbol\":\"");
            out.push_str(&json_escape(&mol.atoms[target_idx].symbol));
            out.push_str("\"}");
        }
        out.push_str("]}");
    }
    out.push(']');
    out
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
        matches_mol(&pat, &mol)
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
        assert_eq!(p.bonds[0].specs, vec![BondSpec::Default]);
    }

    #[test]
    fn parse_explicit_bonds() {
        assert_eq!(parse_smarts("C=O").unwrap().bonds[0].specs, vec![BondSpec::Double]);
        assert_eq!(parse_smarts("C#N").unwrap().bonds[0].specs, vec![BondSpec::Triple]);
        assert_eq!(parse_smarts("c:c").unwrap().bonds[0].specs, vec![BondSpec::Aromatic]);
        assert_eq!(parse_smarts("C-C").unwrap().bonds[0].specs, vec![BondSpec::Single]);
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
        assert_eq!(count_matches("[CH3]", "CCO"), 1);
        assert_eq!(count_matches("[CH2]", "CCO"), 1);
        assert_eq!(count_matches("[H1]", "CCO"), 1);
        assert_eq!(count_matches("[OH]", "CCO"), 1);
    }

    #[test]
    fn match_connections() {
        assert_eq!(count_matches("[CX4]", "CCO"), 2);
        assert_eq!(count_matches("[OX2]", "CCO"), 1);
        assert_eq!(count_matches("[CX2]", "CC#N"), 1);
    }

    // --------- Charges ---------

    #[test]
    fn match_charges() {
        assert!(matches_anywhere("[O-]", "CC(=O)[O-]"));
        assert!(matches_anywhere("[NH4+]", "[NH4+]"));
        assert!(matches_anywhere("[N+0]", "CCN"));
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
        assert!(matches_anywhere("[A;!#1]", "CCO"));
    }

    #[test]
    fn match_nh_charge_combo() {
        assert!(!matches_anywhere("[NH3,NH2,NH;+,+2,+3]", "[NH4+]"));
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
        let pat = parse_smarts("[CH4]").unwrap();
        let mol = parse("C").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_c1_ethane() {
        let pat = parse_smarts("[CH3]C").unwrap();
        let mol = parse("CC").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_c3_methanol() {
        let pat = parse_smarts("[CH3][N,O,P,S,F,Cl,Br,I]").unwrap();
        let mol = parse("CO").unwrap();
        assert!(match_at(&pat, &mol, 0));
        assert!(!match_at(&pat, &mol, 1));
    }

    #[test]
    fn crippen_c18_benzene() {
        assert_eq!(count_matches("[cH]", "c1ccccc1"), 6);
    }

    #[test]
    fn crippen_c19_fused_aromatic() {
        let count = count_matches("[c](:a)(:a):a", "c1ccc2ccccc2c1");
        assert_eq!(count, 2, "expected 2 ring-junction aromatic C");
    }

    #[test]
    fn crippen_o5_nitro_oxygen() {
        let pat = parse_smarts("[O]=[#7,#8]").unwrap();
        let mol = parse("O=N").unwrap();
        assert!(match_at(&pat, &mol, 0));
    }

    #[test]
    fn crippen_fallback_cs() {
        let pat = parse_smarts("[#6]").unwrap();
        let mol = parse("CCO").unwrap();
        assert!(match_at(&pat, &mol, 0));
        assert!(match_at(&pat, &mol, 1));
        assert!(!match_at(&pat, &mol, 2));
    }

    #[test]
    fn crippen_halogen_list() {
        assert!(matches_anywhere("[#9-0]", "CF"));
        assert!(matches_anywhere("[#17-0]", "CCl"));
        assert!(matches_anywhere("[#35-0]", "CBr"));
    }

    #[test]
    fn crippen_o2_hydroxyl() {
        assert!(matches_anywhere("[OH,OH2]", "CCO"));
        assert!(matches_anywhere("[OH,OH2]", "O"));
    }

    // --------- Multi-atom traversal ---------

    #[test]
    fn match_three_atom_chain() {
        let pat = parse_smarts("O=CN").unwrap();
        let mol = parse("NC(=O)N").unwrap();
        assert!(mol.atoms.iter().enumerate().any(|(i, _)| match_at(&pat, &mol, i)));
    }

    #[test]
    fn match_branch() {
        let pat = parse_smarts("C(=O)O").unwrap();
        let mol = parse("CC(=O)O").unwrap();
        assert!(match_at(&pat, &mol, 1));
    }

    // --------- H-expanded patterns (Crippen H types) ---------

    #[test]
    fn match_hydrogen_with_explicit_h() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 6);
    }

    #[test]
    fn match_ch3_after_expansion() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[CH3]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 1);
    }

    #[test]
    fn match_connections_after_expansion() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[CX4]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 2);
    }

    #[test]
    fn match_crippen_h1_pattern() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1][#6,#1]").unwrap();
        let matches: Vec<usize> = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .collect();
        assert_eq!(matches.len(), 5);
    }

    #[test]
    fn match_crippen_h2_pattern() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let pat = parse_smarts("[#1]O[CX4,c]").unwrap();
        let count = (0..mol.atoms.len())
            .filter(|&i| match_at(&pat, &mol, i))
            .count();
        assert_eq!(count, 1);
    }

    // --------- NEW: ring closures, ring queries, recursive, counting ---------

    #[test]
    fn ring_closure_benzene_pattern() {
        // 6-membered any-ring pattern matches benzene.
        assert!(matches_anywhere("*1~*~*~*~*~*~1", "c1ccccc1"));
        // 4-membered ring pattern should NOT match benzene.
        assert!(!matches_anywhere("*1~*~*~*~1", "c1ccccc1"));
    }

    #[test]
    fn ring_closure_small_ring() {
        assert!(matches_anywhere("*1~*~*~1", "C1CC1")); // cyclopropane
        assert!(!matches_anywhere("*1~*~*~1", "CCC")); // propane (no ring)
    }

    #[test]
    fn ring_atom_query() {
        // [R] matches ring atoms only
        assert_eq!(count_matches("[R]", "c1ccccc1"), 6);
        assert_eq!(count_matches("[R]", "CCO"), 0);
        // ring sulfur [#16R]
        assert!(matches_anywhere("[#16R]", "C1CCSCC1"));
        assert!(!matches_anywhere("[#16R]", "CCSCC"));
        // [!#6R] = ring atom that's not carbon
        assert!(matches_anywhere("[!#6R]", "c1ccncc1")); // pyridine N in ring
        assert!(!matches_anywhere("[!#6R]", "c1ccccc1")); // benzene: all ring C
    }

    #[test]
    fn ring_bond_specs() {
        // *@*!@*@* : ring bond, then non-ring bond, then ring bond
        // biphenyl-like: two rings joined by a single (non-ring) bond
        assert!(matches_anywhere("*@*!@*@*", "c1ccccc1-c1ccccc1"));
        // benzene alone has no non-ring bond between ring atoms → no match
        assert!(!matches_anywhere("*@*!@*@*", "c1ccccc1"));
    }

    #[test]
    fn not_in_ring_bond() {
        // *!@* matches any non-ring bond
        assert!(matches_anywhere("[#6]!@[#6]", "CC")); // ethane single non-ring bond
        assert!(!matches_anywhere("[#6]!@[#6]", "C1CC1")); // cyclopropane: all C-C in ring
    }

    #[test]
    fn recursive_smarts() {
        // [$([CH3])] matches a methyl carbon via recursion
        assert_eq!(count_matches("[$([CH3])]", "CCO"), 1);
        // [$(*~[CH2]~[CH2]~*)] anchored at atoms that start such a chain
        let pat = parse_smarts("[$(*~[CH2]~[CH2]~*)]").unwrap();
        let mol = parse("CCCC").unwrap();
        assert!((0..mol.atoms.len()).any(|i| match_at(&pat, &mol, i)));
    }

    #[test]
    fn count_unique_uniquify() {
        // [#8] on a molecule with 2 oxygens → count 2
        let pat = parse_smarts("[#8]").unwrap();
        let mol = parse("OCCO").unwrap();
        assert_eq!(count_unique(&pat, &mol), 2);
        // [CH3] on ethane: 2 methyls
        let pat = parse_smarts("[CH3]").unwrap();
        let mol = parse("CC").unwrap();
        assert_eq!(count_unique(&pat, &mol), 2);
    }

    #[test]
    fn count_unique_dedups_symmetric() {
        // [#6]~[#6] on ethane: one undirected bond → uniquify gives 1, not 2.
        let pat = parse_smarts("[#6]~[#6]").unwrap();
        let mol = parse("CC").unwrap();
        assert_eq!(count_unique(&pat, &mol), 1);
    }

    #[test]
    fn unique_matches_preserve_pattern_order_after_dedup() {
        let pat = parse_smarts("[#6]~[#6]").unwrap();
        let mol = parse("CCC").unwrap();
        assert_eq!(unique_matches(&pat, &mol), vec![vec![0, 1], vec![1, 2]]);

        let pat = parse_smarts("C=O").unwrap();
        let mol = parse("CC(=O)O").unwrap();
        assert_eq!(unique_matches(&pat, &mol), vec![vec![1, 2]]);
    }

    #[test]
    fn unique_matches_json_reports_atom_mapping() {
        let pat = parse_smarts("C=O").unwrap();
        let mol = parse("CC(=O)O").unwrap();
        assert_eq!(
            unique_matches_json(&pat, &mol),
            "[{\"match\":1,\"atom_indices\":[2,3],\"atoms\":[{\"query_atom\":1,\"target_atom\":2,\"symbol\":\"C\"},{\"query_atom\":2,\"target_atom\":3,\"symbol\":\"O\"}]}]"
        );
    }

    #[test]
    fn not_aromatic_bond() {
        // [#16]!:*:* — sulfur attached by non-aromatic bond to an aromatic chain
        // thiophene-substituted: Cc1ccccc1 has methyl by non-aromatic bond
        assert!(matches_anywhere("[#6]!:*:*", "Cc1ccccc1"));
    }
}
