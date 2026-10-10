use crate::weights;
use std::collections::{BTreeMap, HashMap};

const ORGANIC_SUBSET: &[&str] = &["B", "C", "N", "O", "P", "S", "F", "Cl", "Br", "I"];

fn is_organic(sym: &str) -> bool {
    ORGANIC_SUBSET.contains(&sym)
}

fn is_aromatic_char(c: char) -> bool {
    matches!(c, 'b' | 'c' | 'n' | 'o' | 'p' | 's')
}

fn default_valence(sym: &str, aromatic: bool) -> i32 {
    match sym {
        "B" => {
            if aromatic {
                2
            } else {
                3
            }
        }
        "C" => {
            if aromatic {
                3
            } else {
                4
            }
        }
        "N" => {
            if aromatic {
                2
            } else {
                3
            }
        }
        "O" => 2,
        "P" => 3,
        "S" => 2,
        "F" | "Cl" | "Br" | "I" => 1,
        _ => 0,
    }
}

fn bond_valence(order: BondOrder) -> i32 {
    match order {
        BondOrder::Single | BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
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
    /// Cis/trans direction of a single bond next to a double bond, as written
    /// from `a` to `b`: `1` for `/`, `-1` for `\`, `0` when unspecified.
    /// Reading the bond in the other direction inverts the sign.
    pub direction: i8,
}

impl Bond {
    pub fn new(a: usize, b: usize, order: BondOrder) -> Self {
        Self {
            a,
            b,
            order,
            direction: 0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct Atom {
    pub symbol: String,
    pub hydrogen: i32,
    pub charge: i32,
    pub aromatic: bool,
    pub in_bracket: bool,
    pub isotope: Option<u16>,
    pub atom_map: Option<u32>,
    pub chirality: Option<String>,
    /// Neighbours in the order they were written in the input, which is the
    /// reference frame `@`/`@@` is defined against. `-1` marks the slot of an
    /// implicit hydrogen inside brackets, `-2` an unresolved ring-closure slot.
    /// Empty when the atom carries no chirality.
    pub nbr_order: Vec<i32>,
}

impl Default for Atom {
    fn default() -> Self {
        Self {
            symbol: String::new(),
            hydrogen: 0,
            charge: 0,
            aromatic: false,
            in_bracket: false,
            isotope: None,
            atom_map: None,
            chirality: None,
            nbr_order: Vec::new(),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Molecule {
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub bond_count: i32,
}

/// Dense 0-based ranks of `keys`: equal keys share a rank, and a rank equals
/// the number of strictly smaller keys, so the result depends on the key values
/// alone and not on the order they were supplied in.
fn ranks_from_keys<K: Ord>(keys: &[K]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..keys.len()).collect();
    order.sort_by(|&a, &b| keys[a].cmp(&keys[b]));
    let mut rank = vec![0usize; keys.len()];
    let mut current = 0;
    for i in 0..order.len() {
        if i > 0 && keys[order[i]] != keys[order[i - 1]] {
            current = i;
        }
        rank[order[i]] = current;
    }
    rank
}

/// Index of the lowest set bit of a GF(2) row, or `None` for a zero row.
fn pivot_of(row: &[u64]) -> Option<usize> {
    row.iter()
        .enumerate()
        .find_map(|(w, v)| (*v != 0).then(|| w * 64 + v.trailing_zeros() as usize))
}

/// Reduce `vector` against a row-reduced span; true when it reduces to zero,
/// meaning the cycle is a sum of the cycles already in the span.
fn in_span(span: &[Vec<u64>], vector: &[u64]) -> bool {
    let mut vector = vector.to_vec();
    for row in span {
        let Some(p) = pivot_of(row) else { continue };
        if vector[p / 64] >> (p % 64) & 1 == 1 {
            for (w, value) in row.iter().enumerate() {
                vector[w] ^= value;
            }
        }
    }
    vector.iter().all(|w| *w == 0)
}

fn insert_into_span(span: &mut Vec<Vec<u64>>, vector: Vec<u64>) {
    let mut vector = vector;
    for row in span.iter() {
        let Some(p) = pivot_of(row) else { continue };
        if vector[p / 64] >> (p % 64) & 1 == 1 {
            for (w, value) in row.iter().enumerate() {
                vector[w] ^= value;
            }
        }
    }
    if vector.iter().any(|w| *w != 0) {
        span.push(vector);
    }
}

fn class_count(rank: &[usize]) -> usize {
    let mut seen: Vec<usize> = rank.to_vec();
    seen.sort_unstable();
    seen.dedup();
    seen.len()
}

fn bond_order_code(order: BondOrder) -> u8 {
    match order {
        BondOrder::Single => 0,
        BondOrder::Aromatic => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
    }
}

/// Parity of the permutation taking `reference` to `emitted`. Both list the
/// same neighbours; an odd permutation inverts a tetrahedral centre, which is
/// what `@` has to be flipped for.
fn permutation_is_odd(reference: &[i32], emitted: &[i32]) -> bool {
    if reference.len() != emitted.len() || reference.len() < 3 {
        return false;
    }
    let mut perm: Vec<usize> = Vec::with_capacity(emitted.len());
    let mut used = vec![false; reference.len()];
    for value in emitted {
        match reference
            .iter()
            .enumerate()
            .position(|(i, r)| r == value && !used[i])
        {
            Some(pos) => {
                used[pos] = true;
                perm.push(pos);
            }
            // Neighbour sets do not line up (an implicit hydrogen appeared or
            // vanished); leave the tag alone rather than invent a parity.
            None => return false,
        }
    }
    let mut swaps = 0;
    let mut perm = perm;
    for i in 0..perm.len() {
        while perm[i] != i {
            let j = perm[i];
            perm.swap(i, j);
            swaps += 1;
        }
    }
    swaps % 2 == 1
}

/// Total pi-electron count of a ring, or `None` when any atom disqualifies it.
fn sum_pi(contributions: &[Option<u32>]) -> Option<u32> {
    let mut total = 0;
    for c in contributions {
        total += (*c)?;
    }
    Some(total)
}

/// Result of ring perception over a molecule.
#[derive(Clone, Debug)]
pub struct RingInfo {
    /// Per-atom: true if the atom lies on at least one ring.
    pub atom_in_ring: Vec<bool>,
    /// Per-bond (indexed like `Molecule::bonds`): true if the bond is a ring bond.
    pub bond_in_ring: Vec<bool>,
    /// A smallest-set-of-smallest-rings, each ring as a list of bond indices.
    pub rings: Vec<Vec<usize>>,
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
                    ..Default::default()
                });
                bonds.push(Bond {
                    a: i,
                    b: h_idx,
                    order: BondOrder::Single,
                    direction: 0,
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

        // RDKit appends the net charge: `+`, `-`, `+2`, `-2`, …
        let charge = self.total_charge();
        match charge {
            0 => {}
            1 => result.push('+'),
            -1 => result.push('-'),
            c if c > 0 => result.push_str(&format!("+{}", c)),
            c => result.push_str(&format!("-{}", -c)),
        }
        result
    }

    /// Re-derive every atom's hydrogen count from its bonds, the way the parser
    /// does for atoms written without brackets. Used after a molecule has been
    /// rebuilt (scaffold extraction, element flattening) and the counts carried
    /// over from the parent no longer fit the new valences.
    pub fn recompute_implicit_hydrogens(&mut self) {
        let mut bonded = vec![0i32; self.atoms.len()];
        for bond in &self.bonds {
            let v = bond_valence(bond.order);
            bonded[bond.a] += v;
            bonded[bond.b] += v;
        }
        for (idx, atom) in self.atoms.iter_mut().enumerate() {
            atom.hydrogen = (default_valence(&atom.symbol, atom.aromatic) - bonded[idx]).max(0);
        }
    }

    pub fn heavy_atom_count(&self) -> usize {
        self.atoms.iter().filter(|a| a.symbol != "H").count()
    }

    pub fn total_charge(&self) -> i32 {
        self.atoms.iter().map(|atom| atom.charge).sum()
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

    /// Build an adjacency list of (neighbor, bond_index) for every atom.
    /// Bond index refers into `self.bonds`. Used by ring perception and by
    /// SMARTS ring-bond (`@` / `!@`) matching.
    pub fn adjacency(&self) -> Vec<Vec<(usize, usize)>> {
        let mut adj = vec![Vec::new(); self.atoms.len()];
        for (bi, b) in self.bonds.iter().enumerate() {
            if b.a < adj.len() && b.b < adj.len() {
                adj[b.a].push((b.b, bi));
                adj[b.b].push((b.a, bi));
            }
        }
        adj
    }

    /// Number of disconnected fragments (connected components) over the
    /// atom/bond graph. Matches RDKit `MolOps::getMolFrags(...).size()`,
    /// which is what MACCS bit 166 (`> 1` fragment) keys on.
    pub fn fragment_count(&self) -> usize {
        let n = self.atoms.len();
        if n == 0 {
            return 0;
        }
        let adj = self.adjacency();
        let mut seen = vec![false; n];
        let mut comps = 0;
        for start in 0..n {
            if seen[start] {
                continue;
            }
            comps += 1;
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(u) = stack.pop() {
                for &(v, _) in &adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        stack.push(v);
                    }
                }
            }
        }
        comps
    }

    /// Perceive aromaticity from Kekule input with a Hueckel 4n+2 pi-electron
    /// count over each smallest ring, then over each fused ring system.
    ///
    /// Every decision is taken against a snapshot of the *original* bond orders,
    /// so a ring is never judged differently depending on whether a fused
    /// neighbour happened to be rewritten first. Rings that arrive already
    /// aromatic (lowercase SMILES) are left untouched.
    pub fn perceive_aromaticity(&mut self) {
        let ring_info = self.ring_info();
        if ring_info.rings.is_empty() {
            return;
        }
        let orig: Vec<Bond> = self.bonds.clone();
        let mut system_atoms = vec![false; self.atoms.len()];
        let mut system_bonds = vec![false; self.bonds.len()];

        let ring_atoms: Vec<Vec<usize>> = ring_info
            .rings
            .iter()
            .map(|ring| Self::ring_atom_list(&orig, ring))
            .collect();

        // --- pass 1: each smallest ring on its own ---
        let mut aromatic = vec![false; ring_info.rings.len()];
        let mut pi_per_atom: Vec<Vec<Option<u32>>> = Vec::with_capacity(ring_info.rings.len());
        for (ri, ring) in ring_info.rings.iter().enumerate() {
            let atoms = &ring_atoms[ri];
            let contributions: Vec<Option<u32>> = atoms
                .iter()
                .map(|&idx| self.pi_contribution(idx, ring, &orig))
                .collect();
            let usable = atoms.len() == ring.len() && atoms.len() >= 3;
            let already_aromatic = ring
                .iter()
                .all(|&bi| matches!(orig[bi].order, BondOrder::Aromatic));
            if usable && !already_aromatic {
                if let Some(total) = sum_pi(&contributions) {
                    if total >= 2 && total % 4 == 2 {
                        aromatic[ri] = true;
                    }
                }
            }
            pi_per_atom.push(contributions);
        }

        // --- pass 2: fused ring systems (azulene-like) ---
        // Rings sharing a bond form a system; if every atom of the system is a
        // valid sp2 contributor and the system total is 4n+2, the whole system
        // is aromatic even when no individual ring is.
        for system in Self::fused_ring_systems(&ring_info.rings) {
            if system.len() < 2 || system.iter().all(|&ri| aromatic[ri]) {
                continue;
            }
            let mut atoms: Vec<usize> = Vec::new();
            for &ri in &system {
                for &idx in &ring_atoms[ri] {
                    if !atoms.contains(&idx) {
                        atoms.push(idx);
                    }
                }
            }
            let bonds: Vec<usize> = {
                let mut v: Vec<usize> = Vec::new();
                for &ri in &system {
                    for &bi in &ring_info.rings[ri] {
                        if !v.contains(&bi) {
                            v.push(bi);
                        }
                    }
                }
                v
            };
            let contributions: Vec<Option<u32>> = atoms
                .iter()
                .map(|&idx| self.pi_contribution(idx, &bonds, &orig))
                .collect();
            let already_aromatic = bonds
                .iter()
                .all(|&bi| matches!(orig[bi].order, BondOrder::Aromatic));
            if already_aromatic {
                continue;
            }
            if let Some(total) = sum_pi(&contributions) {
                if total >= 2 && total % 4 == 2 {
                    for &idx in &atoms {
                        system_atoms[idx] = true;
                    }
                    // Aromatic atoms need not share an aromatic fusion bond.
                    // Only the perimeter of a fused 4n+2 system is aromatic;
                    // internal bonds are aromatic only via an individual ring.
                    for &bi in &bonds {
                        if system
                            .iter()
                            .filter(|&&ri| ring_info.rings[ri].contains(&bi))
                            .count()
                            == 1
                        {
                            system_bonds[bi] = true;
                        }
                    }
                }
            }
        }

        for (ri, is_aromatic) in aromatic.iter().enumerate() {
            if !is_aromatic {
                continue;
            }
            for &idx in &ring_atoms[ri] {
                self.atoms[idx].aromatic = true;
            }
            for &bi in &ring_info.rings[ri] {
                self.bonds[bi].order = BondOrder::Aromatic;
            }
        }
        for (idx, aromatic) in system_atoms.into_iter().enumerate() {
            self.atoms[idx].aromatic |= aromatic;
        }
        for (bi, aromatic) in system_bonds.into_iter().enumerate() {
            if aromatic {
                self.bonds[bi].order = BondOrder::Aromatic;
            }
        }
    }

    /// Unique atom indices touched by a ring given as a list of bond indices.
    fn ring_atom_list(bonds: &[Bond], ring: &[usize]) -> Vec<usize> {
        let mut atoms = Vec::with_capacity(ring.len());
        for &bi in ring {
            let Some(bond) = bonds.get(bi) else {
                continue;
            };
            for idx in [bond.a, bond.b] {
                if !atoms.contains(&idx) {
                    atoms.push(idx);
                }
            }
        }
        atoms
    }

    /// Group rings into fused systems: two rings belong to the same system when
    /// they share at least one bond, transitively.
    fn fused_ring_systems(rings: &[Vec<usize>]) -> Vec<Vec<usize>> {
        let n = rings.len();
        let mut parent: Vec<usize> = (0..n).collect();
        fn find(parent: &mut Vec<usize>, mut x: usize) -> usize {
            while parent[x] != x {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            x
        }
        for i in 0..n {
            for j in (i + 1)..n {
                if rings[i].iter().any(|bi| rings[j].contains(bi)) {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }
        let mut systems: Vec<Vec<usize>> = vec![Vec::new(); n];
        for i in 0..n {
            let root = find(&mut parent, i);
            systems[root].push(i);
        }
        systems.into_iter().filter(|s| !s.is_empty()).collect()
    }

    /// Pi-electron contribution of one ring atom, or `None` when the atom
    /// disqualifies the ring (sp3 centre, triple bond, unsupported element).
    ///
    /// `ring` is the bond-index set the atom is being judged within, so a bond
    /// counts as "exocyclic" relative to that particular ring or ring system.
    fn pi_contribution(&self, idx: usize, ring: &[usize], orig: &[Bond]) -> Option<u32> {
        let atom = self.atoms.get(idx)?;
        // A neutral bracket nitrogen with an incomplete valence shell is a
        // radical, not a pyrrolic lone-pair donor (C1=C[N]C=C1).
        if atom.in_bracket && atom.symbol == "N" && atom.charge == 0 {
            let valence = atom.hydrogen
                + orig
                    .iter()
                    .filter(|b| b.a == idx || b.b == idx)
                    .map(|b| bond_valence(b.order))
                    .sum::<i32>();
            if !atom.aromatic && valence < 3 {
                return None;
            }
        }
        if !matches!(
            atom.symbol.as_str(),
            "B" | "C" | "N" | "O" | "P" | "S" | "Se" | "Te" | "As"
        ) {
            return None;
        }

        let mut in_ring_double = false;
        let mut in_ring_aromatic = false;
        let mut exocyclic_double_to: Option<&str> = None;
        for (bi, bond) in orig.iter().enumerate() {
            let other = if bond.a == idx {
                bond.b
            } else if bond.b == idx {
                bond.a
            } else {
                continue;
            };
            if matches!(bond.order, BondOrder::Triple) {
                return None;
            }
            if ring.contains(&bi) {
                match bond.order {
                    BondOrder::Double => in_ring_double = true,
                    BondOrder::Aromatic => in_ring_aromatic = true,
                    _ => {}
                }
            } else if matches!(bond.order, BondOrder::Double) {
                exocyclic_double_to = Some(self.atoms[other].symbol.as_str());
            }
        }

        // An exocyclic double bond to a more electronegative atom pulls the pi
        // pair away from the ring, leaving an empty p orbital that contributes
        // nothing — this is what keeps quinones non-aromatic while the caffeine
        // pyrimidinedione ring stays aromatic. A double bond to another carbon
        // still donates one electron into the ring system, which is what makes
        // the naphthalene half of acenaphthylene come out aromatic.
        if let Some(partner) = exocyclic_double_to {
            return match partner {
                "O" | "N" | "S" | "Se" | "Te" => Some(0),
                _ => Some(1),
            };
        }
        if in_ring_double {
            return Some(1);
        }
        if in_ring_aromatic {
            return Some(match atom.symbol.as_str() {
                "O" | "S" | "Se" | "Te" => 2,
                "N" | "P" | "As"
                    if atom.charge == 0
                        && (atom.hydrogen > 0
                            || orig.iter().filter(|b| b.a == idx || b.b == idx).count() == 3) =>
                {
                    2
                }
                _ => 1,
            });
        }

        // No double bond at this atom: it can only take part through a lone
        // pair (or an empty orbital when positively charged).
        match (atom.symbol.as_str(), atom.charge) {
            (_, c) if c > 0 => Some(0),
            ("C", c) if c < 0 => Some(2),
            ("C", _) => None,
            ("B", _) => Some(0),
            ("N" | "P" | "As", _) => Some(2),
            ("O" | "S" | "Se" | "Te", _) => Some(2),
            _ => None,
        }
    }

    /// Ring perception. Returns, for the whole molecule:
    ///   - `atom_in_ring`: per-atom flag (true if the atom lies on any cycle)
    ///   - `bond_in_ring`: per-bond flag (true if the bond lies on any cycle)
    ///   - `rings`: a Smallest-Set-of-Smallest-Rings as lists of bond indices
    ///
    /// A bond is a ring bond iff it is *not* a bridge: removing it leaves its
    /// endpoints still connected. Ring atoms are the endpoints of ring bonds.
    /// The SSSR is built greedily by collecting, for each ring bond, the
    /// shortest cycle through it (BFS shortest path between its endpoints in the
    /// graph with that bond removed), then deduplicating by bond-set. This is
    /// sufficient for MACCS, which only needs ring sizes and aromatic-ring
    /// counts, not a canonical SSSR.
    pub fn ring_info(&self) -> RingInfo {
        let n = self.atoms.len();
        let nb = self.bonds.len();
        let adj = self.adjacency();

        // --- ring bonds = non-bridges ---
        // A bond is a ring bond iff removing it leaves its endpoints connected.
        // Molecules are tiny, so a per-bond connectivity check (BFS with that
        // bond skipped) is both simple and unambiguously correct.
        let mut bond_in_ring = vec![false; nb];
        for (bi, b) in self.bonds.iter().enumerate() {
            if Self::connected_without(b.a, b.b, bi, &adj, n) {
                bond_in_ring[bi] = true;
            }
        }

        let mut atom_in_ring = vec![false; n];
        for (bi, b) in self.bonds.iter().enumerate() {
            if bond_in_ring[bi] {
                atom_in_ring[b.a] = true;
                atom_in_ring[b.b] = true;
            }
        }

        // --- symmetrised smallest set of smallest rings ---
        //
        // The rings kept are the *relevant* cycles: those that belong to at
        // least one minimum cycle basis. A cycle qualifies when it cannot be
        // written as a sum, over GF(2) in the bond space, of strictly shorter
        // cycles. This is what makes bicyclo[2.2.2]octane report three
        // six-membered rings rather than the cyclomatic number of two, and it
        // is the same set RDKit's symmetrised SSSR returns.
        let mut candidates: Vec<Vec<usize>> = Vec::new();
        let mut seen_keys: std::collections::HashSet<Vec<usize>> = std::collections::HashSet::new();
        let push_candidate =
            |cycle: Vec<usize>,
             candidates: &mut Vec<Vec<usize>>,
             seen: &mut std::collections::HashSet<Vec<usize>>| {
                let mut key = cycle.clone();
                key.sort_unstable();
                if seen.insert(key) {
                    candidates.push(cycle);
                }
            };

        for (bi, b) in self.bonds.iter().enumerate() {
            if !bond_in_ring[bi] {
                continue;
            }
            if let Some(cycle) = self.shortest_cycle_through(b.a, b.b, bi, &adj) {
                push_candidate(cycle, &mut candidates, &mut seen_keys);
            }
        }
        // Horton-style candidates: for every ring atom and every ring bond, the
        // cycle made of the two shortest paths from that atom to the bond's
        // ends. This is what turns up the bridging rings a per-bond scan misses.
        for root in 0..n {
            if !atom_in_ring[root] {
                continue;
            }
            for (bi, b) in self.bonds.iter().enumerate() {
                if !bond_in_ring[bi] {
                    continue;
                }
                if let Some(cycle) = self.cycle_through_via(root, b.a, b.b, bi, &adj) {
                    push_candidate(cycle, &mut candidates, &mut seen_keys);
                }
            }
        }

        candidates.sort_by(|a, b| a.len().cmp(&b.len()).then_with(|| a.cmp(b)));

        let words = nb.div_ceil(64).max(1);
        let to_vector = |cycle: &[usize]| {
            let mut vector = vec![0u64; words];
            for &bi in cycle {
                vector[bi / 64] ^= 1u64 << (bi % 64);
            }
            vector
        };
        // Row-reduced span of every ring accepted at a strictly smaller size.
        let mut smaller_span: Vec<Vec<u64>> = Vec::new();
        let mut rings: Vec<Vec<usize>> = Vec::new();
        let mut pending: Vec<Vec<u64>> = Vec::new();
        let mut current_size = 0usize;

        for cycle in candidates {
            if cycle.len() != current_size {
                for row in pending.drain(..) {
                    insert_into_span(&mut smaller_span, row);
                }
                current_size = cycle.len();
            }
            let vector = to_vector(&cycle);
            if !in_span(&smaller_span, &vector) {
                pending.push(vector);
                rings.push(cycle);
            }
        }

        RingInfo {
            atom_in_ring,
            bond_in_ring,
            rings,
        }
    }

    /// Shortest cycle that runs through bond `(a, b)` and also visits `root`:
    /// the shortest `root..a` and `root..b` paths joined by that bond, when the
    /// two paths share no vertex other than `root`.
    fn cycle_through_via(
        &self,
        root: usize,
        a: usize,
        b: usize,
        bond: usize,
        adj: &[Vec<(usize, usize)>],
    ) -> Option<Vec<usize>> {
        let (path_a, atoms_a) = self.shortest_path(root, a, bond, adj)?;
        let (path_b, atoms_b) = self.shortest_path(root, b, bond, adj)?;
        if atoms_a
            .iter()
            .filter(|x| **x != root)
            .any(|x| atoms_b.contains(x))
        {
            return None;
        }
        let mut cycle = path_a;
        cycle.extend(path_b);
        cycle.push(bond);
        cycle.sort_unstable();
        cycle.dedup();
        Some(cycle)
    }

    /// BFS shortest path, returning its bond indices and the atoms it visits.
    fn shortest_path(
        &self,
        from: usize,
        to: usize,
        skip: usize,
        adj: &[Vec<(usize, usize)>],
    ) -> Option<(Vec<usize>, Vec<usize>)> {
        let n = self.atoms.len();
        let mut prev: Vec<(usize, usize)> = vec![(usize::MAX, usize::MAX); n];
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(from);
        visited[from] = true;
        while let Some(u) = queue.pop_front() {
            if u == to {
                break;
            }
            for &(v, bi) in &adj[u] {
                if bi == skip || visited[v] {
                    continue;
                }
                visited[v] = true;
                prev[v] = (u, bi);
                queue.push_back(v);
            }
        }
        if !visited[to] {
            return None;
        }
        let mut bonds = Vec::new();
        let mut atoms = vec![to];
        let mut cur = to;
        while cur != from {
            let (p, bi) = prev[cur];
            if p == usize::MAX {
                return None;
            }
            bonds.push(bi);
            atoms.push(p);
            cur = p;
        }
        Some((bonds, atoms))
    }

    /// True if `a` can reach `b` in the graph with bond `skip` removed.
    fn connected_without(
        a: usize,
        b: usize,
        skip: usize,
        adj: &[Vec<(usize, usize)>],
        n: usize,
    ) -> bool {
        let mut visited = vec![false; n];
        let mut stack = vec![a];
        visited[a] = true;
        while let Some(u) = stack.pop() {
            if u == b {
                return true;
            }
            for &(v, bi) in &adj[u] {
                if bi == skip || visited[v] {
                    continue;
                }
                visited[v] = true;
                stack.push(v);
            }
        }
        false
    }

    /// BFS shortest path between `a` and `b` in the graph with bond `skip`
    /// removed, returning the cycle's bond indices (the path bonds plus `skip`).
    fn shortest_cycle_through(
        &self,
        a: usize,
        b: usize,
        skip: usize,
        adj: &[Vec<(usize, usize)>],
    ) -> Option<Vec<usize>> {
        let n = self.atoms.len();
        let mut prev: Vec<(usize, usize)> = vec![(usize::MAX, usize::MAX); n]; // (prev_atom, bond)
        let mut visited = vec![false; n];
        let mut queue = std::collections::VecDeque::new();
        queue.push_back(a);
        visited[a] = true;
        while let Some(u) = queue.pop_front() {
            if u == b {
                break;
            }
            for &(v, bi) in &adj[u] {
                if bi == skip || visited[v] {
                    continue;
                }
                visited[v] = true;
                prev[v] = (u, bi);
                queue.push_back(v);
            }
        }
        if !visited[b] {
            return None;
        }
        let mut bonds = vec![skip];
        let mut cur = b;
        while cur != a {
            let (p, bi) = prev[cur];
            if bi == usize::MAX {
                return None;
            }
            bonds.push(bi);
            cur = p;
        }
        Some(bonds)
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
            emit_ring_closure(
                atom_idx,
                *other_idx,
                *order,
                ring_id_for_pair,
                next_ring_id,
                output,
            );
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

    /// Canonical atom ranking: a permutation of `0..atoms.len()` derived only
    /// from the molecular graph, never from the order atoms happened to be
    /// written in. Equivalent atoms are separated by successive refinement of
    /// their neighbourhoods, and any classes that survive refinement (genuinely
    /// symmetric atoms) are split one at a time so the ranking is total.
    pub fn canonical_ranks(&self) -> Vec<usize> {
        let n = self.atoms.len();
        if n == 0 {
            return Vec::new();
        }
        let adj = self.adjacency();

        // Initial invariant: everything about an atom that a re-spelling of the
        // same molecule cannot change.
        let initial: Vec<(usize, String, i32, i32, bool, u32)> = (0..n)
            .map(|i| {
                let atom = &self.atoms[i];
                (
                    adj[i].len(),
                    atom.symbol.clone(),
                    atom.charge,
                    atom.hydrogen,
                    atom.aromatic,
                    atom.isotope.map(u32::from).unwrap_or(0),
                )
            })
            .collect();
        let mut rank = ranks_from_keys(&initial);

        loop {
            let refined = self.refine_ranks(&rank, &adj);
            if class_count(&refined) == class_count(&rank) {
                rank = refined;
                break;
            }
            rank = refined;
        }

        // Break ties until the ranking is a strict order. Members of a surviving
        // class are indistinguishable to the refinement, so picking any of them
        // yields the same output string.
        while class_count(&rank) < n {
            let target = (0..n)
                .map(|i| rank[i])
                .filter(|r| (0..n).filter(|&i| rank[i] == *r).count() > 1)
                .min()
                .expect("a duplicated rank exists");
            let pick = (0..n)
                .filter(|&i| rank[i] == target)
                .min()
                .expect("class is non-empty");
            for r in rank.iter_mut() {
                if *r > target {
                    *r += 1;
                }
            }
            for (i, r) in rank.iter_mut().enumerate() {
                if *r == target && i != pick {
                    *r += 1;
                }
            }
            loop {
                let refined = self.refine_ranks(&rank, &adj);
                if class_count(&refined) == class_count(&rank) {
                    rank = refined;
                    break;
                }
                rank = refined;
            }
        }
        rank
    }

    fn refine_ranks(&self, rank: &[usize], adj: &[Vec<(usize, usize)>]) -> Vec<usize> {
        let keys: Vec<(usize, Vec<(u8, usize)>)> = (0..self.atoms.len())
            .map(|i| {
                let mut nbrs: Vec<(u8, usize)> = adj[i]
                    .iter()
                    .map(|&(v, bi)| (bond_order_code(self.bonds[bi].order), rank[v]))
                    .collect();
                nbrs.sort_unstable();
                (rank[i], nbrs)
            })
            .collect();
        ranks_from_keys(&keys)
    }

    /// Canonical SMILES: identical for every spelling of the same molecule.
    ///
    /// Atoms are ordered by [`Molecule::canonical_ranks`], the traversal starts
    /// at the root that produces the lexicographically smallest string, and
    /// tetrahedral parity is recomputed against the neighbour order actually
    /// emitted, so re-spelling a molecule can never flip its stereocentres.
    pub fn canonical_smiles(&self) -> String {
        if self.atoms.is_empty() {
            return String::new();
        }
        let ranks = self.canonical_ranks();
        let mut molecule = self.clone();
        molecule.normalize_bond_directions(&ranks);
        let mut rendered = Vec::new();
        for component in molecule.components() {
            let mut best: Option<String> = None;
            for &start in &component {
                let candidate = molecule.component_smiles_from(start, &component, &ranks);
                if best.as_ref().map(|s| candidate < *s).unwrap_or(true) {
                    best = Some(candidate);
                }
            }
            if let Some(s) = best {
                rendered.push(s);
            }
        }
        rendered.sort();
        rendered.join(".")
    }

    /// Rewrite `/` and `\` marks into a canonical frame.
    ///
    /// Cis/trans is a relation between the two double-bond substituents, and
    /// flipping every mark around one double bond describes the same molecule.
    /// Anchoring the `/` on the lowest-ranked substituent of the lowest-ranked
    /// double-bond carbon fixes that freedom, so `F/C=C\F` and `C(=C/F)/F` stop
    /// canonicalising to mirror-image strings.
    fn normalize_bond_directions(&mut self, ranks: &[usize]) {
        let mut assigned = vec![0i8; self.bonds.len()];
        let directional_at = |mol: &Molecule, atom: usize, skip: usize| -> Vec<(usize, usize)> {
            mol.bonds
                .iter()
                .enumerate()
                .filter_map(|(bi, bond)| {
                    if bi == skip || bond.direction == 0 || bond.order != BondOrder::Single {
                        return None;
                    }
                    if bond.a == atom {
                        Some((bond.b, bi))
                    } else if bond.b == atom {
                        Some((bond.a, bi))
                    } else {
                        None
                    }
                })
                .collect()
        };
        let outward = |mol: &Molecule, bi: usize, atom: usize| -> i8 {
            if mol.bonds[bi].a == atom {
                mol.bonds[bi].direction
            } else {
                -mol.bonds[bi].direction
            }
        };

        for di in 0..self.bonds.len() {
            if self.bonds[di].order != BondOrder::Double {
                continue;
            }
            let (mut u, mut v) = (self.bonds[di].a, self.bonds[di].b);
            if ranks[v] < ranks[u] {
                std::mem::swap(&mut u, &mut v);
            }
            let at_u = directional_at(self, u, di);
            let at_v = directional_at(self, v, di);
            if at_u.is_empty() || at_v.is_empty() {
                continue;
            }
            let ref_u = *at_u.iter().min_by_key(|(other, _)| ranks[*other]).unwrap();
            let ref_v = *at_v.iter().min_by_key(|(other, _)| ranks[*other]).unwrap();
            let cis = outward(self, ref_u.1, u) == outward(self, ref_v.1, v);

            for (anchor, refs, sign) in [(u, &at_u, 1i8), (v, &at_v, if cis { 1 } else { -1 })] {
                let reference = if anchor == u { ref_u.1 } else { ref_v.1 };
                for &(_, bi) in refs {
                    let s = if bi == reference { sign } else { -sign };
                    assigned[bi] = if self.bonds[bi].a == anchor { s } else { -s };
                }
            }
        }

        for (bi, bond) in self.bonds.iter_mut().enumerate() {
            bond.direction = assigned[bi];
        }
    }

    pub fn components(&self) -> Vec<Vec<usize>> {
        let n = self.atoms.len();
        let adj = self.adjacency();
        let mut seen = vec![false; n];
        let mut components = Vec::new();
        for start in 0..n {
            if seen[start] {
                continue;
            }
            let mut comp = Vec::new();
            let mut stack = vec![start];
            seen[start] = true;
            while let Some(u) = stack.pop() {
                comp.push(u);
                for &(v, _) in &adj[u] {
                    if !seen[v] {
                        seen[v] = true;
                        stack.push(v);
                    }
                }
            }
            comp.sort_unstable();
            components.push(comp);
        }
        components
    }

    pub fn subgraph_from_atoms(&self, atoms_to_keep: &[usize]) -> Option<Molecule> {
        let mut keep = vec![false; self.atoms.len()];
        for &idx in atoms_to_keep {
            if idx >= keep.len() {
                return None;
            }
            keep[idx] = true;
        }

        let mut index_map = vec![usize::MAX; self.atoms.len()];
        let mut atoms = Vec::new();
        for (idx, atom) in self.atoms.iter().enumerate() {
            if !keep[idx] {
                continue;
            }
            let mut atom = atom.clone();
            // Bonds cut by the subgraph become hydrogens, so a fragment carved
            // out of a larger molecule still has sane valences. Whole
            // components lose nothing and are copied unchanged.
            let lost: i32 = self
                .bonds
                .iter()
                .filter(|bond| (bond.a == idx && !keep[bond.b]) || (bond.b == idx && !keep[bond.a]))
                .map(|bond| bond_valence(bond.order))
                .sum();
            if lost > 0 {
                atom.hydrogen += lost;
                atom.chirality = None;
                atom.nbr_order.clear();
            }
            index_map[idx] = atoms.len();
            atoms.push(atom);
        }
        if atoms.is_empty() {
            return None;
        }
        for atom in &mut atoms {
            atom.nbr_order = atom
                .nbr_order
                .iter()
                .map(|&slot| {
                    if slot < 0 || index_map[slot as usize] == usize::MAX {
                        slot
                    } else {
                        index_map[slot as usize] as i32
                    }
                })
                .collect();
        }

        let mut bonds = Vec::new();
        for bond in &self.bonds {
            let a = index_map[bond.a];
            let b = index_map[bond.b];
            if a != usize::MAX && b != usize::MAX {
                bonds.push(Bond {
                    a,
                    b,
                    order: bond.order,
                    direction: bond.direction,
                });
            }
        }

        Some(Molecule {
            bond_count: bonds.len() as i32,
            atoms,
            bonds,
        })
    }

    fn component_smiles_from(&self, start: usize, component: &[usize], ranks: &[usize]) -> String {
        let mut in_component = vec![false; self.atoms.len()];
        for &idx in component {
            in_component[idx] = true;
        }
        // A first pass classifies which bonds close rings, so the opening digit
        // can be written straight after its atom (`c1ccccc1`) instead of being
        // discovered only when the traversal comes back around.
        let mut seen = vec![false; self.atoms.len()];
        let mut back_edges = vec![false; self.bonds.len()];
        self.collect_back_edges(
            start,
            usize::MAX,
            ranks,
            &in_component,
            &mut seen,
            &mut back_edges,
        );

        let mut visited = vec![false; self.atoms.len()];
        let mut ring_id_for_pair = std::collections::HashMap::new();
        let mut next_ring_id = 1;
        let mut output = String::new();
        self.dfs_canonical_smiles(
            start,
            None,
            ranks,
            &in_component,
            &back_edges,
            &mut visited,
            &mut ring_id_for_pair,
            &mut next_ring_id,
            &mut output,
        );
        output
    }

    /// Neighbours of `atom_idx` in canonical order, as `(neighbour, bond index)`,
    /// skipping the bond the traversal arrived on.
    fn ordered_neighbors(
        &self,
        atom_idx: usize,
        from_bond: usize,
        ranks: &[usize],
        in_component: &[bool],
    ) -> Vec<(usize, usize)> {
        let mut out: Vec<(usize, usize)> = Vec::new();
        for (bi, bond) in self.bonds.iter().enumerate() {
            if bi == from_bond {
                continue;
            }
            let other = if bond.a == atom_idx {
                bond.b
            } else if bond.b == atom_idx {
                bond.a
            } else {
                continue;
            };
            if in_component[other] {
                out.push((other, bi));
            }
        }
        out.sort_by_key(|(idx, bi)| (ranks[*idx], bond_order_code(self.bonds[*bi].order), *bi));
        out
    }

    fn collect_back_edges(
        &self,
        atom_idx: usize,
        from_bond: usize,
        ranks: &[usize],
        in_component: &[bool],
        seen: &mut [bool],
        back_edges: &mut [bool],
    ) {
        seen[atom_idx] = true;
        for (nbr, bi) in self.ordered_neighbors(atom_idx, from_bond, ranks, in_component) {
            if back_edges[bi] {
                continue;
            }
            if seen[nbr] {
                back_edges[bi] = true;
            } else {
                self.collect_back_edges(nbr, bi, ranks, in_component, seen, back_edges);
            }
        }
    }

    /// Index of the bond joining two atoms, if any.
    fn bond_between(&self, a: usize, b: usize) -> Option<usize> {
        self.bonds
            .iter()
            .position(|bond| (bond.a == a && bond.b == b) || (bond.a == b && bond.b == a))
    }

    /// Cis/trans marker for the bond `from -> to`, inverted when the bond is
    /// traversed against the direction it was stored in.
    fn direction_symbol(&self, from: usize, to: usize) -> &'static str {
        let Some(bi) = self.bond_between(from, to) else {
            return "";
        };
        let bond = self.bonds[bi];
        if bond.order != BondOrder::Single || bond.direction == 0 {
            return "";
        }
        let dir = if bond.a == from {
            bond.direction
        } else {
            -bond.direction
        };
        if dir > 0 {
            "/"
        } else {
            "\\"
        }
    }

    #[allow(clippy::too_many_arguments)]
    fn dfs_canonical_smiles(
        &self,
        atom_idx: usize,
        from: Option<(usize, usize)>,
        ranks: &[usize],
        in_component: &[bool],
        back_edges: &[bool],
        visited: &mut [bool],
        ring_id_for_pair: &mut std::collections::HashMap<(usize, usize), usize>,
        next_ring_id: &mut usize,
        output: &mut String,
    ) {
        visited[atom_idx] = true;
        let from_bond = from.map(|(_, bi)| bi).unwrap_or(usize::MAX);
        let pending = self.ordered_neighbors(atom_idx, from_bond, ranks, in_component);

        // Ring-closure digits have to sit right after the atom, ahead of any
        // branch, so the emitted neighbour order is: the atom we came from, the
        // bracketed hydrogen, the ring closures, then the branches.
        let (closures, branches): (Vec<_>, Vec<_>) =
            pending.into_iter().partition(|(_, bi)| back_edges[*bi]);

        let mut emitted: Vec<i32> = Vec::new();
        if let Some((f, _)) = from {
            emitted.push(f as i32);
        }
        let atom = &self.atoms[atom_idx];
        if !self.can_write_bare(atom_idx) && atom.hydrogen > 0 {
            emitted.push(-1);
        }
        for (idx, _) in closures.iter().chain(branches.iter()) {
            emitted.push(*idx as i32);
        }

        output.push_str(&self.atom_token(atom_idx, &emitted));

        for (other_idx, bi) in &closures {
            emit_ring_closure(
                atom_idx,
                *other_idx,
                self.bonds[*bi].order,
                ring_id_for_pair,
                next_ring_id,
                output,
            );
        }

        for (i, (next_idx, bi)) in branches.iter().enumerate() {
            let is_last = i + 1 == branches.len();
            if !is_last {
                output.push('(');
            }
            let symbol = bond_smiles_symbol(self.bonds[*bi].order);
            if symbol.is_empty() {
                output.push_str(self.direction_symbol(atom_idx, *next_idx));
            } else {
                output.push_str(symbol);
            }
            self.dfs_canonical_smiles(
                *next_idx,
                Some((atom_idx, *bi)),
                ranks,
                in_component,
                back_edges,
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

    /// Whether the atom can be written without brackets: it must be in the
    /// organic subset, plain, and its hydrogen count must be exactly what the
    /// bare form implies — otherwise a pyrrole nitrogen would silently lose its
    /// hydrogen on the way out.
    fn can_write_bare(&self, idx: usize) -> bool {
        let atom = &self.atoms[idx];
        if atom.charge != 0
            || atom.symbol == "H"
            || !is_organic(&atom.symbol)
            || atom.isotope.is_some()
            || atom.atom_map.is_some()
            || atom.chirality.is_some()
        {
            return false;
        }
        let bonded: i32 = self
            .neighbors(idx)
            .into_iter()
            .map(|(_, order)| bond_valence(order))
            .sum();
        let implied = (default_valence(&atom.symbol, atom.aromatic) - bonded).max(0);
        implied == atom.hydrogen
    }

    /// Atom token for the canonical writer. `emitted` is the neighbour order
    /// this atom is about to be written with (`-1` marking its bracketed
    /// hydrogen); tetrahedral parity is corrected against the order the input
    /// used, so `@` and `@@` follow the molecule rather than the spelling.
    fn atom_token(&self, idx: usize, emitted: &[i32]) -> String {
        let atom = &self.atoms[idx];
        if self.can_write_bare(idx) {
            return if atom.aromatic {
                atom.symbol.to_ascii_lowercase()
            } else {
                atom.symbol.clone()
            };
        }

        let chirality = atom.chirality.as_deref().map(|tag| {
            if permutation_is_odd(&atom.nbr_order, emitted) {
                if tag == "@" {
                    "@@"
                } else {
                    "@"
                }
            } else if tag == "@" {
                "@"
            } else {
                "@@"
            }
        });

        let mut out = String::new();
        out.push('[');
        if let Some(isotope) = atom.isotope {
            out.push_str(&isotope.to_string());
        }
        if atom.aromatic {
            out.push_str(&atom.symbol.to_ascii_lowercase());
        } else {
            out.push_str(&atom.symbol);
        }
        if let Some(tag) = chirality {
            out.push_str(tag);
        }
        if atom.hydrogen > 0 {
            out.push('H');
            if atom.hydrogen > 1 {
                out.push_str(&atom.hydrogen.to_string());
            }
        }
        if atom.charge > 0 {
            out.push('+');
            if atom.charge > 1 {
                out.push_str(&atom.charge.to_string());
            }
        } else if atom.charge < 0 {
            out.push('-');
            if atom.charge < -1 {
                out.push_str(&(-atom.charge).to_string());
            }
        }
        if let Some(map) = atom.atom_map {
            out.push(':');
            out.push_str(&map.to_string());
        }
        out.push(']');
        out
    }

    /// Classify atom `idx` by hybridization: SP / SP2 / SP3.
    pub fn hybridization(&self, idx: usize) -> crate::conformer::params::Hybridization {
        use crate::conformer::params::Hybridization;
        let nbrs = self.neighbors(idx);
        if nbrs.iter().any(|(_, o)| *o == BondOrder::Triple) {
            return Hybridization::SP;
        }
        if self.atoms[idx].aromatic
            || nbrs
                .iter()
                .any(|(_, o)| matches!(o, BondOrder::Double | BondOrder::Aromatic))
        {
            return Hybridization::SP2;
        }
        Hybridization::SP3
    }
}

fn bond_smiles_symbol(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Single | BondOrder::Aromatic => "",
        BondOrder::Double => "=",
        BondOrder::Triple => "#",
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
    let mut atom = Atom {
        in_bracket: true,
        ..Default::default()
    };

    let isotope_start = *pos;
    while *pos < chars.len() && chars[*pos].is_ascii_digit() {
        *pos += 1;
    }
    if *pos > isotope_start {
        let isotope: String = chars[isotope_start..*pos].iter().collect();
        atom.isotope = isotope.parse::<u16>().ok();
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

    if *pos < chars.len() && chars[*pos] == '@' {
        let start = *pos;
        *pos += 1;
        if *pos < chars.len() && chars[*pos] == '@' {
            *pos += 1;
        }
        atom.chirality = Some(chars[start..*pos].iter().collect());
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

    if *pos < chars.len() && chars[*pos] == ':' {
        *pos += 1;
        let start = *pos;
        while *pos < chars.len() && chars[*pos].is_ascii_digit() {
            *pos += 1;
        }
        if start == *pos {
            return None;
        }
        let atom_map: String = chars[start..*pos].iter().collect();
        atom.atom_map = atom_map.parse::<u32>().ok();
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
    // Ring bond number -> (opening atom, its slot in `nbr_order`, direction).
    let mut ring_openings: HashMap<i32, (i32, usize, i8)> = HashMap::new();
    let mut prev_atom: i32 = -1;
    let mut next_bond_order: i32 = 1;
    let mut explicit_bond: Option<BondOrder> = None;
    let mut next_direction: i8 = 0;
    let mut pos = 0;

    while pos < chars.len() {
        let c = chars[pos];

        // Branch
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
            next_direction = match c {
                '/' => 1,
                '\\' => -1,
                _ => 0,
            };
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
                if pos + 1 >= chars.len()
                    || !chars[pos].is_ascii_digit()
                    || !chars[pos + 1].is_ascii_digit()
                {
                    return None;
                }
                ring_num = ((chars[pos] as i32) - ('0' as i32)) * 10
                    + ((chars[pos + 1] as i32) - ('0' as i32));
                pos += 2;
            } else {
                ring_num = (c as i32) - ('0' as i32);
                pos += 1;
            }

            if prev_atom < 0 {
                return None;
            }
            if let Some((other, slot, open_direction)) = ring_openings.remove(&ring_num) {
                let a = other as usize;
                let b = prev_atom as usize;
                let order = resolve_bond_order(
                    explicit_bond,
                    next_bond_order,
                    atoms[a].aromatic,
                    atoms[b].aromatic,
                );
                let mut bond = Bond::new(a, b, order);
                bond.direction = if open_direction != 0 {
                    open_direction
                } else {
                    // A direction written at the closing digit reads from the
                    // closing atom, so it flips when stored as a -> b.
                    -next_direction
                };
                bonds.push(bond);
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                atoms[a].nbr_order[slot] = b as i32;
                atoms[b].nbr_order.push(a as i32);
                next_bond_order = 1;
                explicit_bond = None;
                next_direction = 0;
            } else {
                let idx = prev_atom as usize;
                atoms[idx].nbr_order.push(-2);
                let slot = atoms[idx].nbr_order.len() - 1;
                ring_openings.insert(ring_num, (prev_atom, slot, next_direction));
                next_direction = 0;
            }
            continue;
        }

        // Bracket atom
        if c == '[' {
            pos += 1;
            let atom = parse_bracket_atom(&chars, &mut pos)?;
            let has_implicit_h = atom.hydrogen > 0;
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
                let mut bond = Bond::new(a, b, order);
                bond.direction = next_direction;
                bonds.push(bond);
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                atoms[a].nbr_order.push(b as i32);
                atoms[b].nbr_order.push(a as i32);
                next_bond_order = 1;
                explicit_bond = None;
                next_direction = 0;
            }
            // Inside brackets the hydrogen sits immediately after the preceding
            // atom in the neighbour order that `@`/`@@` refers to.
            if has_implicit_h {
                atoms[idx as usize].nbr_order.push(-1);
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
                let mut bond = Bond::new(a, b, order);
                bond.direction = next_direction;
                bonds.push(bond);
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                atoms[a].nbr_order.push(b as i32);
                atoms[b].nbr_order.push(a as i32);
                next_bond_order = 1;
                explicit_bond = None;
                next_direction = 0;
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
                let mut bond = Bond::new(a, b, order);
                bond.direction = next_direction;
                bonds.push(bond);
                bond_count += 1;
                degree[a] += next_bond_order;
                degree[b] += next_bond_order;
                atoms[a].nbr_order.push(b as i32);
                atoms[b].nbr_order.push(a as i32);
                next_bond_order = 1;
                explicit_bond = None;
                next_direction = 0;
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

    let mut mol = Molecule {
        atoms,
        bonds,
        bond_count,
    };
    // Implicit hydrogens are filled in from the bond orders *as written*. Doing
    // this before aromatic perception is what keeps a Kekule pyrrole nitrogen
    // (`C1=CNC=C1`, valence 3, two single bonds -> 1 H) equivalent to `[nH]`;
    // deriving it from the perceived aromatic valence would silently drop that H.
    let mut degree = vec![0; mol.atoms.len()];
    for bond in &mol.bonds {
        let valence = bond_valence(bond.order);
        degree[bond.a] += valence;
        degree[bond.b] += valence;
    }

    for (i, atom) in mol.atoms.iter_mut().enumerate() {
        if atom.in_bracket {
            continue;
        }
        let val = default_valence(&atom.symbol, atom.aromatic);
        let implicit_h = (val - degree[i]).max(0);
        atom.hydrogen = implicit_h;
    }

    mol.perceive_aromaticity();

    Some(mol)
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
            assert_eq!(
                b.order,
                BondOrder::Aromatic,
                "benzene bond should be aromatic"
            );
        }
    }

    #[test]
    fn test_kekule_benzene_perceives_aromaticity() {
        let mol = parse("C1=CC=CC=C1").unwrap();
        assert_eq!(mol.formula(), "C6H6");
        assert!(mol.atoms.iter().all(|atom| atom.aromatic));
        assert!(mol
            .bonds
            .iter()
            .all(|bond| bond.order == BondOrder::Aromatic));
        assert_eq!(mol.canonical_smiles(), "c1ccccc1");
    }

    #[test]
    fn test_kekule_pyridine_perceives_aromaticity() {
        let mol = parse("C1=CC=NC=C1").unwrap();
        assert_eq!(mol.formula(), "C5H5N");
        assert!(mol.atoms.iter().all(|atom| atom.aromatic));
        assert_eq!(mol.canonical_smiles(), "c1ccccn1");
    }

    #[test]
    fn test_nonaromatic_cyclohexane_stays_aliphatic() {
        let mol = parse("C1CCCCC1").unwrap();
        assert_eq!(mol.formula(), "C6H12");
        assert!(mol.atoms.iter().all(|atom| !atom.aromatic));
        assert!(mol.bonds.iter().all(|bond| bond.order == BondOrder::Single));
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
        let s = parse("C")
            .unwrap()
            .with_explicit_hydrogens()
            .to_smiles_verbose();
        // [C]([H])([H])([H])[H]
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 5);
        assert_eq!(parsed.heavy_atom_count(), 1);
    }

    #[test]
    fn test_to_smiles_verbose_ethanol() {
        let s = parse("CCO")
            .unwrap()
            .with_explicit_hydrogens()
            .to_smiles_verbose();
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 9);
        assert_eq!(parsed.heavy_atom_count(), 3);
    }

    #[test]
    fn test_to_smiles_verbose_water() {
        let s = parse("O")
            .unwrap()
            .with_explicit_hydrogens()
            .to_smiles_verbose();
        let parsed = parse(&s).expect("round-trip parse");
        assert_eq!(parsed.atoms.len(), 3);
        assert_eq!(parsed.heavy_atom_count(), 1);
    }

    #[test]
    fn test_to_smiles_verbose_benzene() {
        let s = parse("c1ccccc1")
            .unwrap()
            .with_explicit_hydrogens()
            .to_smiles_verbose();
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
        assert_eq!(
            mol.formula(),
            expected_formula,
            "Formula mismatch for {}",
            smiles
        );
        assert_eq!(
            mol.heavy_atom_count(),
            expected_heavy,
            "Heavy atom count mismatch for {}",
            smiles
        );
    }

    /// Helper: parse must succeed (formula not checked, just parsability + heavy atoms)
    fn check_parse(smiles: &str, expected_heavy: usize) {
        let mol = parse(smiles).unwrap_or_else(|| panic!("Failed to parse: {}", smiles));
        assert_eq!(
            mol.heavy_atom_count(),
            expected_heavy,
            "Heavy atom count mismatch for {}",
            smiles
        );
    }

    // --- Ring perception & fragments ---

    #[test]
    fn ring_benzene() {
        let mol = parse("c1ccccc1").unwrap();
        let ri = mol.ring_info();
        assert!(
            ri.atom_in_ring.iter().all(|&x| x),
            "all benzene atoms in ring"
        );
        assert_eq!(ri.rings.len(), 1, "benzene has 1 ring");
        assert_eq!(ri.rings[0].len(), 6, "benzene ring has 6 bonds");
    }

    #[test]
    fn ring_acyclic_has_no_rings() {
        let mol = parse("CCO").unwrap();
        let ri = mol.ring_info();
        assert!(ri.atom_in_ring.iter().all(|&x| !x));
        assert!(ri.bond_in_ring.iter().all(|&x| !x));
        assert_eq!(ri.rings.len(), 0);
    }

    #[test]
    fn ring_naphthalene_two_rings() {
        let mol = parse("c1ccc2ccccc2c1").unwrap();
        let ri = mol.ring_info();
        assert_eq!(ri.rings.len(), 2, "naphthalene has 2 rings");
        // 10 atoms all in ring
        assert_eq!(ri.atom_in_ring.iter().filter(|&&x| x).count(), 10);
    }

    #[test]
    fn ring_substituent_not_in_ring() {
        // toluene: the methyl C is not in the ring
        let mol = parse("Cc1ccccc1").unwrap();
        let ri = mol.ring_info();
        assert!(!ri.atom_in_ring[0], "methyl carbon not in ring");
        assert_eq!(ri.atom_in_ring.iter().filter(|&&x| x).count(), 6);
    }

    #[test]
    fn fragments_single_and_salt() {
        assert_eq!(parse("CCO").unwrap().fragment_count(), 1);
        assert_eq!(parse("[Na+].[Cl-]").unwrap().fragment_count(), 2);
        assert_eq!(parse("CC(=O)[O-].[Na+]").unwrap().fragment_count(), 2);
    }

    // --- Simple organic molecules ---

    #[test]
    fn test_methane() {
        check("C", "CH4", 1);
    }

    #[test]
    fn test_ethane() {
        check("CC", "C2H6", 2);
    }

    #[test]
    fn test_propane() {
        check("CCC", "C3H8", 3);
    }

    #[test]
    fn test_butane() {
        check("CCCC", "C4H10", 4);
    }

    #[test]
    fn test_isobutane() {
        check("CC(C)C", "C4H10", 4);
    }

    #[test]
    fn test_neopentane() {
        check("CC(C)(C)C", "C5H12", 5);
    }

    #[test]
    fn test_ethylene() {
        check("C=C", "C2H4", 2);
    }

    #[test]
    fn test_acetylene() {
        check("C#C", "C2H2", 2);
    }

    #[test]
    fn test_propylene() {
        check("CC=C", "C3H6", 3);
    }

    #[test]
    fn test_1_3_butadiene() {
        check("C=CC=C", "C4H6", 4);
    }

    // --- Alcohols, aldehydes, ketones, acids ---

    #[test]
    fn test_methanol() {
        check("CO", "CH4O", 2);
    }

    #[test]
    fn test_formaldehyde() {
        check("C=O", "CH2O", 2);
    }

    #[test]
    fn test_acetone() {
        check("CC(=O)C", "C3H6O", 4);
    }

    #[test]
    fn test_acetic_acid() {
        check("CC(=O)O", "C2H4O2", 4);
    }

    #[test]
    fn test_formic_acid() {
        check("O=CO", "CH2O2", 3);
    }

    #[test]
    fn test_glycerol() {
        check("OCC(O)CO", "C3H8O3", 6);
    }

    // --- Amines, amides ---

    #[test]
    fn test_methylamine() {
        check("CN", "CH5N", 2);
    }

    #[test]
    fn test_dimethylamine() {
        check("CNC", "C2H7N", 3);
    }

    #[test]
    fn test_trimethylamine() {
        check("CN(C)C", "C3H9N", 4);
    }

    #[test]
    fn test_urea() {
        check("NC(=O)N", "CH4N2O", 4);
    }

    #[test]
    fn test_acetamide() {
        check("CC(=O)N", "C2H5NO", 4);
    }

    // --- Halogens ---

    #[test]
    fn test_chloromethane() {
        check("CCl", "CH3Cl", 2);
    }

    #[test]
    fn test_bromomethane() {
        check("CBr", "CH3Br", 2);
    }

    #[test]
    fn test_iodomethane() {
        check("CI", "CH3I", 2);
    }

    #[test]
    fn test_fluoromethane() {
        check("CF", "CH3F", 2);
    }

    #[test]
    fn test_dichloromethane() {
        check("ClCCl", "CH2Cl2", 3);
    }

    #[test]
    fn test_chloroform() {
        check("ClC(Cl)Cl", "CHCl3", 4);
    }

    #[test]
    fn test_carbon_tetrachloride() {
        check("ClC(Cl)(Cl)Cl", "CCl4", 5);
    }

    // --- Aromatic compounds ---

    #[test]
    fn test_toluene() {
        check("Cc1ccccc1", "C7H8", 7);
    }

    #[test]
    fn test_phenol() {
        check("Oc1ccccc1", "C6H6O", 7);
    }

    #[test]
    fn test_aniline() {
        check("Nc1ccccc1", "C6H7N", 7);
    }

    #[test]
    fn test_benzoic_acid() {
        check("OC(=O)c1ccccc1", "C7H6O2", 9);
    }

    #[test]
    fn test_nitrobenzene() {
        check("c1ccc([N+](=O)[O-])cc1", "C6H5NO2", 9);
    }

    #[test]
    fn test_styrene() {
        check("C=Cc1ccccc1", "C8H8", 8);
    }

    #[test]
    fn test_biphenyl() {
        check("c1ccc(-c2ccccc2)cc1", "C12H10", 12);
    }

    #[test]
    fn test_anthracene() {
        check("c1ccc2cc3ccccc3cc2c1", "C14H10", 14);
    }

    // --- Heterocycles ---

    #[test]
    fn test_pyridine() {
        check("c1ccncc1", "C5H5N", 6);
    }

    #[test]
    fn test_pyrrole() {
        check("c1cc[nH]c1", "C4H5N", 5);
    }

    #[test]
    fn test_furan() {
        check("c1ccoc1", "C4H4O", 5);
    }

    #[test]
    fn test_thiophene() {
        check("c1ccsc1", "C4H4S", 5);
    }

    #[test]
    fn test_imidazole() {
        check("c1c[nH]cn1", "C3H4N2", 5);
    }

    #[test]
    fn test_indole() {
        check("c1ccc2[nH]ccc2c1", "C8H7N", 9);
    }

    #[test]
    fn test_quinoline() {
        check("c1ccc2ncccc2c1", "C9H7N", 10);
    }

    // --- Rings ---

    #[test]
    fn test_cyclopropane() {
        check("C1CC1", "C3H6", 3);
    }

    #[test]
    fn test_cyclobutane() {
        check("C1CCC1", "C4H8", 4);
    }

    #[test]
    fn test_cyclopentane() {
        check("C1CCCC1", "C5H10", 5);
    }

    #[test]
    fn test_cycloheptane() {
        check("C1CCCCCC1", "C7H14", 7);
    }

    #[test]
    fn test_cyclooctane() {
        check("C1CCCCCCC1", "C8H16", 8);
    }

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
        check_parse(
            "CC1([C@@H](N2[C@H](S1)[C@@H](C2=O)NC(=O)Cc3ccccc3)C(=O)O)C",
            23,
        );
    }

    #[test]
    fn test_cholesterol() {
        check_parse(
            "C[C@H](CCCC(C)C)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CC=C4[C@@]3(CC[C@@H](C4)O)C)C",
            28,
        );
    }

    #[test]
    fn test_glucose() {
        check_parse("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", 12); // 6C + 6O = 12 heavy
    }

    #[test]
    fn test_sucrose() {
        check_parse(
            "OC[C@H]1OC(O[C@@]2(CO)O[C@H](CO)[C@@H](O)[C@@H]2O)[C@H](O)[C@@H](O)[C@@H]1O",
            23,
        ); // 12C + 11O = 23 heavy
    }

    // --- Amino acids ---

    #[test]
    fn test_glycine() {
        check("NCC(=O)O", "C2H5NO2", 5);
    }

    #[test]
    fn test_alanine() {
        check("C[C@@H](N)C(=O)O", "C3H7NO2", 6);
    }

    #[test]
    fn test_valine() {
        check("CC(C)[C@@H](N)C(=O)O", "C5H11NO2", 8);
    }

    #[test]
    fn test_leucine() {
        check("CC(C)C[C@@H](N)C(=O)O", "C6H13NO2", 9);
    }

    #[test]
    fn test_phenylalanine() {
        check_parse("[C@@H](Cc1ccccc1)(N)C(=O)O", 12);
    }

    #[test]
    fn test_tryptophan() {
        check_parse("N[C@@H](Cc1c[nH]c2ccccc12)C(=O)O", 15);
    }

    #[test]
    fn test_cysteine() {
        check("N[C@@H](CS)C(=O)O", "C3H7NO2S", 7);
    }

    #[test]
    fn test_methionine() {
        check("CSCC[C@@H](N)C(=O)O", "C5H11NO2S", 9);
    }

    #[test]
    fn test_proline() {
        check_parse("OC(=O)[C@@H]1CCCN1", 8);
    }

    #[test]
    fn test_histidine() {
        check_parse("N[C@@H](Cc1c[nH]cn1)C(=O)O", 11);
    }

    // --- Sulfur compounds ---

    #[test]
    fn test_dimethyl_sulfoxide() {
        check("CS(=O)C", "C2H6OS", 4);
    }

    #[test]
    fn test_thioacetone() {
        check("CC(=S)C", "C3H6S", 4);
    }

    #[test]
    fn test_methanethiol() {
        check("CS", "CH4S", 2);
    }

    // --- Phosphorus ---

    #[test]
    fn test_trimethylphosphine() {
        check("CP(C)C", "C3H9P", 4);
    }

    // --- Boron ---

    #[test]
    fn test_borane() {
        check_parse("[BH3]", 1);
    }

    #[test]
    fn test_phenylboronic_acid() {
        check("OB(O)c1ccccc1", "C6H7BO2", 9);
    }

    // --- Charged species ---

    #[test]
    fn test_ammonium() {
        check_parse("[NH4+]", 1);
    }

    #[test]
    fn test_hydroxide() {
        check_parse("[OH-]", 1);
    }

    #[test]
    fn test_acetate() {
        check_parse("CC(=O)[O-]", 4);
    }

    #[test]
    fn test_sodium_chloride() {
        check_parse("[Na+].[Cl-]", 2);
    }

    #[test]
    fn test_calcium_chloride() {
        check_parse("[Ca+2].[Cl-].[Cl-]", 3);
    }

    #[test]
    fn test_sulfate() {
        check_parse("[O-]S(=O)(=O)[O-]", 5);
    }

    // --- Stereo SMILES (chirality/cis-trans) — just check parsability ---

    #[test]
    fn test_cis_2_butene() {
        check_parse(r"C/C=C\C", 4);
    }

    #[test]
    fn test_trans_2_butene() {
        check_parse("C/C=C/C", 4);
    }

    #[test]
    fn test_l_alanine() {
        check_parse("[C@@H](N)(C)C(=O)O", 6);
    }

    #[test]
    fn test_d_alanine() {
        check_parse("[C@H](N)(C)C(=O)O", 6);
    }

    // --- Multi-ring / fused systems ---

    #[test]
    fn test_adamantane() {
        check("C1C2CC3CC1CC(C2)C3", "C10H16", 10);
    }

    #[test]
    fn test_cubane() {
        check("C12C3C4C1C5C3C4C25", "C8H8", 8);
    }

    #[test]
    fn test_decalin() {
        check("C1CCC2CCCCC2C1", "C10H18", 10);
    }

    #[test]
    fn test_fluorene() {
        check("c1ccc2c(c1)Cc1ccccc1-2", "C13H10", 13);
    }

    // --- Two-digit ring closures ---

    #[test]
    fn test_two_digit_ring() {
        // %10 ring closure
        let mol = parse("C%10CCCCCCCCC%10").unwrap();
        assert_eq!(mol.heavy_atom_count(), 10);
    }

    // --- Edge cases ---

    #[test]
    fn test_single_bracket_atom() {
        check_parse("[Cu]", 1);
    }

    #[test]
    fn test_isotope_bracket() {
        check_parse("[13CH4]", 1);
    }

    #[test]
    fn test_wildcard_bracket() {
        check_parse("[*]", 1);
    }

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
    fn test_adenine() {
        check("c1nc(N)c2nc[nH]c2n1", "C5H5N5", 10);
    }

    #[test]
    fn test_guanine() {
        check("c1nc2c(n1)[nH]c(=O)n2N", "C4H4N5O", 10);
    }

    #[test]
    fn test_cytosine() {
        check("c1cnc(=O)[nH]c1N", "C4H5N3O", 8);
    }

    #[test]
    fn test_thymine() {
        check("Cc1c[nH]c(=O)[nH]c1=O", "C5H6N2O2", 9);
    }

    #[test]
    fn test_uracil() {
        check("c1c[nH]c(=O)[nH]c1=O", "C4H4N2O2", 8);
    }

    // =================================================================
    // Steroids (beyond cholesterol)
    // =================================================================

    #[test]
    fn test_estradiol() {
        check_parse("C[C@]12CC[C@H]3[C@@H](CCc4cc(O)ccc43)[C@@H]1CC[C@@H]2O", 20);
    }

    #[test]
    fn test_progesterone() {
        check_parse(
            "CC(=O)[C@H]1CC[C@@H]2[C@@]1(CC[C@H]3[C@H]2CCC4=CC(=O)CC[C@@]34C)C",
            23,
        );
    }

    #[test]
    fn test_cortisol() {
        check_parse(
            "C[C@@]12C[C@@H](O)[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@]2(O)C(=O)CO",
            26,
        );
    }

    #[test]
    fn test_testosterone() {
        check_parse(
            "C[C@]12CC[C@H]3[C@@H](CCC4=CC(=O)CC[C@@]43C)[C@@H]1CC[C@@H]2O",
            21,
        );
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
        check_parse(
            "CN(Cc1cnc2nc(N)nc(N)c2n1)c1ccc(C(=O)N[C@@H](CCC(=O)O)C(=O)O)cc1",
            33,
        );
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
        check_parse(
            "OC[C@H]1OC(O[C@@H]2[C@@H](O)[C@H](O)[C@@H](O)OC2CO)[C@H](O)[C@@H](O)[C@@H]1O",
            23,
        );
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
        check_parse(
            "COc1ccc2C[C@H]3N(C)CC[C@@]45c2c1O[C@H]4[C@@H](O)C=C[C@@H]35",
            22,
        );
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
        // The nitrate groups need their [O-]; written without it the molecule
        // really does carry a +3 charge, which the formula now reports.
        check(
            "[O-][N+](=O)OCC(O[N+]([O-])=O)CO[N+]([O-])=O",
            "C3H5N3O9",
            15,
        );
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
    fn test_isoleucine() {
        check("CC[C@H](C)[C@@H](N)C(=O)O", "C6H13NO2", 9);
    }

    #[test]
    fn test_serine() {
        check("N[C@@H](CO)C(=O)O", "C3H7NO3", 7);
    }

    #[test]
    fn test_threonine() {
        check("C[C@@H](O)[C@@H](N)C(=O)O", "C4H9NO3", 8);
    }

    #[test]
    fn test_aspartate() {
        check("N[C@@H](CC(=O)O)C(=O)O", "C4H7NO4", 9);
    }

    #[test]
    fn test_asparagine() {
        check("N[C@@H](CC(=O)N)C(=O)O", "C4H8N2O3", 9);
    }

    #[test]
    fn test_glutamine() {
        check("N[C@@H](CCC(=O)N)C(=O)O", "C5H10N2O3", 10);
    }

    #[test]
    fn test_lysine() {
        check("NCCCC[C@@H](N)C(=O)O", "C6H14N2O2", 10);
    }

    #[test]
    fn test_arginine() {
        check_parse("N[C@@H](CCCNC(=N)N)C(=O)O", 12);
    }

    #[test]
    fn test_tyrosine() {
        check_parse("N[C@@H](Cc1ccc(O)cc1)C(=O)O", 13);
    }

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

    // ── hybridization tests ───────────────────────────────────────────────────

    #[test]
    fn hybridization_sp3_ethane() {
        use crate::conformer::params::Hybridization;
        let mol = parse("CC").unwrap();
        assert_eq!(mol.hybridization(0), Hybridization::SP3);
        assert_eq!(mol.hybridization(1), Hybridization::SP3);
    }

    #[test]
    fn hybridization_sp2_ethylene() {
        use crate::conformer::params::Hybridization;
        let mol = parse("C=C").unwrap();
        assert_eq!(mol.hybridization(0), Hybridization::SP2);
        assert_eq!(mol.hybridization(1), Hybridization::SP2);
    }

    #[test]
    fn hybridization_sp_acetylene() {
        use crate::conformer::params::Hybridization;
        let mol = parse("C#C").unwrap();
        assert_eq!(mol.hybridization(0), Hybridization::SP);
        assert_eq!(mol.hybridization(1), Hybridization::SP);
    }

    #[test]
    fn hybridization_sp2_benzene() {
        use crate::conformer::params::Hybridization;
        let mol = parse("c1ccccc1").unwrap();
        for i in 0..6 {
            assert_eq!(mol.hybridization(i), Hybridization::SP2, "atom {i}");
        }
    }

    #[test]
    fn hybridization_mixed_aspirin() {
        // CC(=O)Oc1ccccc1C(=O)O
        use crate::conformer::params::Hybridization;
        let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
        // atom 0 = methyl C (SP3)
        assert_eq!(mol.hybridization(0), Hybridization::SP3);
        // atom 1 = carbonyl C (SP2, has C=O)
        assert_eq!(mol.hybridization(1), Hybridization::SP2);
    }
}
