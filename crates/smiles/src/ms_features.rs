//! DM003 research descriptors. Versioned graph hypotheses, not an ion simulator.
//! Atom-order independent; no stereochemical/CIP equivalence claim. Fragment
//! formulas retain parent H counts. No automatic valence completion or capping.
use crate::parser::{parse, BondOrder, Molecule};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
pub const SPEC: &str = "duck-ms18-v1";
pub const NAMES: [&str; 18] = [
    "single_mass",
    "fragment_formula",
    "bond_types",
    "edge_environment",
    "mass_environment",
    "neutral_losses",
    "charge_sites",
    "double_mass",
    "ring_mass",
    "pathways",
    "functional_groups",
    "group_distances",
    "hetero_placement",
    "ring_topology",
    "conjugation_branch",
    "composition",
    "properties",
    "path_fingerprint",
];
pub const ELEMENTS: [&str; 13] = [
    "H", "C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "B", "Si", "Se",
];
pub const H: f64 = 1.00782503223;
pub const PROTON: f64 = 1.007276466621;
pub(crate) type Adj = Vec<Vec<(usize, usize, u64)>>;
type Counts = BTreeMap<usize, u32>;
#[derive(Clone, Debug)]
pub struct Fragment {
    pub atoms: Vec<usize>,
    pub formula: [u16; 13],
    pub mass: f64,
    pub cuts: u8,
    pub ring: bool,
    pub basic: u16,
    pub acidic: u16,
    pub environment: u64,
}
#[derive(Clone, Debug)]
pub struct Features {
    pub channels: Vec<Vec<(usize, u32)>>,
    pub properties: Vec<f64>,
    pub fragments: Vec<Fragment>,
    pub formula: [u16; 13],
    pub mass: f64,
}
pub(crate) fn hash(v: &[u64]) -> u64 {
    let mut h = 14695981039346656037u64;
    for x in v {
        for b in x.to_le_bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(1099511628211);
        }
    }
    h
}
pub(crate) fn elem(s: &str) -> Result<usize, String> {
    ELEMENTS
        .iter()
        .position(|x| *x == s)
        .ok_or_else(|| format!("unsupported element {s}"))
}
pub(crate) fn order(o: BondOrder) -> u64 {
    match o {
        BondOrder::Single => 1,
        BondOrder::Double => 2,
        BondOrder::Triple => 3,
        BondOrder::Aromatic => 4,
    }
}
pub(crate) fn adjacency(m: &Molecule) -> Adj {
    let mut a = vec![vec![]; m.atoms.len()];
    for (i, b) in m.bonds.iter().enumerate() {
        a[b.a].push((b.b, i, order(b.order)));
        a[b.b].push((b.a, i, order(b.order)));
    }
    a
}
pub(crate) fn distances(a: &Adj, root: usize, cuts: &[usize]) -> Vec<usize> {
    let mut d = vec![usize::MAX; a.len()];
    d[root] = 0;
    let mut q = VecDeque::from([root]);
    while let Some(u) = q.pop_front() {
        for &(v, b, _) in &a[u] {
            if !cuts.contains(&b) && d[v] == usize::MAX {
                d[v] = d[u] + 1;
                q.push_back(v);
            }
        }
    }
    d
}
fn components(a: &Adj, cuts: &[usize]) -> Vec<Vec<usize>> {
    let mut seen = vec![false; a.len()];
    let mut result = vec![];
    for root in 0..a.len() {
        if seen[root] {
            continue;
        }
        let d = distances(a, root, cuts);
        let mut side = vec![];
        for (i, x) in d.iter().enumerate() {
            if *x != usize::MAX {
                seen[i] = true;
                side.push(i);
            }
        }
        result.push(side);
    }
    result
}
fn emit(c: &mut Counts, width: usize, v: &[u64]) {
    *c.entry((hash(v) % width as u64) as usize).or_default() += 1;
}
fn formula(m: &Molecule, ids: &[usize]) -> Result<[u16; 13], String> {
    let mut f = [0; 13];
    for &i in ids {
        let a = &m.atoms[i];
        f[elem(&a.symbol)?] += 1;
        f[0] += a.hydrogen as u16;
    }
    Ok(f)
}
fn mass(f: &[u16; 13]) -> f64 {
    f.iter()
        .zip(ELEMENTS)
        .map(|(&c, s)| f64::from(c) * crate::weights::monoisotopic_mass(s).unwrap())
        .sum()
}
// Named graph predicates, independent of proton affinity or pKa. Multiple tags
// may coexist at one atom. Tags are anchored to the atom, never SMARTS match order.
pub(crate) fn groups(m: &Molecule, a: &Adj) -> Vec<Vec<u64>> {
    let carbonyl: Vec<bool> = (0..a.len())
        .map(|i| {
            m.atoms[i].symbol == "C"
                && a[i]
                    .iter()
                    .any(|&(j, _, o)| o == 2 && m.atoms[j].symbol == "O")
        })
        .collect();
    m.atoms
        .iter()
        .enumerate()
        .map(|(i, x)| {
            let mut g = vec![];
            if x.symbol == "O" && x.hydrogen > 0 {
                g.push(1);
            }
            if x.symbol == "O" && x.hydrogen > 0 && a[i].iter().any(|&(j, _, _)| carbonyl[j]) {
                g.push(2);
            }
            if x.symbol == "O"
                && x.hydrogen == 0
                && a[i].len() == 2
                && a[i].iter().all(|x| x.2 == 1)
            {
                g.push(3);
            }
            if x.symbol == "N" {
                let amide = a[i].iter().any(|&(j, _, o)| o == 1 && carbonyl[j]);
                let nitrile = a[i].iter().any(|x| x.2 == 3);
                if amide {
                    g.push(5);
                } else if nitrile {
                    g.push(6);
                } else if x.aromatic {
                    g.push(7);
                } else {
                    g.push(4);
                }
            }
            if carbonyl[i] {
                g.push(8);
                if a[i]
                    .iter()
                    .any(|&(j, _, o)| o == 1 && m.atoms[j].symbol == "O")
                {
                    g.push(9);
                }
                if a[i]
                    .iter()
                    .any(|&(j, _, o)| o == 1 && m.atoms[j].symbol == "N")
                {
                    g.push(10);
                }
            }
            if x.symbol == "S"
                && a[i]
                    .iter()
                    .any(|&(j, _, o)| o == 2 && m.atoms[j].symbol == "O")
            {
                g.push(11);
            }
            if x.symbol == "P" && a[i].iter().any(|&(j, _, _)| m.atoms[j].symbol == "O") {
                g.push(12);
            }
            if x.symbol == "S" && x.hydrogen > 0 {
                g.push(13);
            }
            if x.charge > 0 {
                g.push(14);
            }
            if x.charge < 0 {
                g.push(15);
            }
            if matches!(x.symbol.as_str(), "F" | "Cl" | "Br" | "I") {
                g.push(16);
            }
            g
        })
        .collect()
}
pub(crate) fn basic(g: &[u64], m: &Molecule, i: usize) -> bool {
    (g.contains(&4) || (g.contains(&7) && m.atoms[i].hydrogen == 0)) && m.atoms[i].charge >= 0
}
pub(crate) fn acidic(g: &[u64]) -> bool {
    g.iter().any(|x| [2, 12, 13, 15].contains(x))
}
fn frag(
    m: &Molecule,
    ids: Vec<usize>,
    cuts: u8,
    ring: bool,
    env: u64,
    g: &[Vec<u64>],
) -> Result<Fragment, String> {
    let f = formula(m, &ids)?;
    Ok(Fragment {
        mass: mass(&f),
        formula: f,
        basic: ids.iter().filter(|&&i| basic(&g[i], m, i)).count() as u16,
        acidic: ids.iter().filter(|&&i| acidic(&g[i])).count() as u16,
        atoms: ids,
        cuts,
        ring,
        environment: env,
    })
}
/// All two-edge cuts are considered within the explicit size bound; none are
/// truncated by atom/bond order. Isotopes/charged/disconnected molecules fail
/// explicitly because this neutral-fragment hypothesis has no valid convention
/// for them. These are research-scope bounds, not parser capabilities.
pub fn calculate(smiles: &str, width: usize) -> Result<Features, String> {
    if !(64..=65536).contains(&width) {
        return Err("width outside 64..65536".into());
    }
    let m = parse(smiles).ok_or("parse failed")?;
    let n = m.atoms.len();
    let nb = m.bonds.len();
    if n == 0 || n > 128 || nb > 192 {
        return Err("graph exceeds research bounds".into());
    }
    for x in &m.atoms {
        elem(&x.symbol)?;
        if x.isotope.is_some() {
            return Err("isotopic mass not supported by ms18-v1".into());
        }
        if x.hydrogen < 0 {
            return Err("negative H count".into());
        }
    }
    if m.atoms.iter().map(|a| a.charge).sum::<i32>() != 0 {
        return Err("neutral parent required".into());
    }
    let adj = adjacency(&m);
    if components(&adj, &[]).len() != 1 {
        return Err("connected parent required".into());
    }
    let f = formula(&m, &(0..n).collect::<Vec<_>>())?;
    let totalmass = mass(&f);
    let g = groups(&m, &adj);
    let ds: Vec<_> = (0..n).map(|i| distances(&adj, i, &[])).collect();
    let mut c = vec![Counts::new(); 18];
    let mut cycle = vec![0; nb];
    let mut singles = vec![];
    let mut labels: Vec<Vec<u64>> = vec![(0..n)
        .map(|i| {
            hash(&[
                elem(&m.atoms[i].symbol).unwrap() as u64,
                m.atoms[i].hydrogen as u64,
                m.atoms[i].charge as u64,
                m.atoms[i].aromatic as u64,
                adj[i].len() as u64,
            ])
        })
        .collect()];
    for r in 1..=3 {
        let prev = &labels[r - 1];
        let next = (0..n)
            .map(|i| {
                let mut ns: Vec<_> = adj[i]
                    .iter()
                    .map(|&(j, _, o)| hash(&[o, prev[j]]))
                    .collect();
                ns.sort_unstable();
                let mut v = vec![r as u64, prev[i]];
                v.extend(ns);
                hash(&v)
            })
            .collect();
        labels.push(next);
    }
    for (bi, b) in m.bonds.iter().enumerate() {
        let d = distances(&adj, b.a, &[bi]);
        if d[b.b] != usize::MAX {
            cycle[bi] = d[b.b] + 1;
        }
        let mut e = [
            elem(&m.atoms[b.a].symbol)? as u64,
            elem(&m.atoms[b.b].symbol)? as u64,
        ];
        e.sort_unstable();
        emit(
            &mut c[2],
            width,
            &[3, e[0], e[1], order(b.order), (cycle[bi] > 0) as u64],
        );
        for (r, ls) in labels.iter().enumerate().skip(1) {
            let mut p = [ls[b.a], ls[b.b]];
            p.sort_unstable();
            emit(&mut c[3], width, &[4, r as u64, p[0], p[1], order(b.order)]);
        }
        if cycle[bi] > 0 {
            emit(
                &mut c[13],
                width,
                &[14, cycle[bi] as u64, e[0], e[1], order(b.order)],
            );
            continue;
        }
        for side in components(&adj, &[bi]) {
            let root = if side.binary_search(&b.a).is_ok() {
                b.a
            } else {
                b.b
            };
            let other = if root == b.a { b.b } else { b.a };
            let env = hash(&[labels[1][root], labels[1][other], order(b.order)]);
            let s = frag(&m, side, 1, false, env, &g)?;
            if (1.0..1000.0).contains(&s.mass) {
                emit(&mut c[0], width, &[604, s.mass.round() as u64]);
            }
            let mut fv = vec![2];
            fv.extend(s.formula.iter().map(|&v| u64::from(v)));
            emit(&mut c[1], width, &fv);
            for (el, &cnt) in s.formula.iter().enumerate() {
                if cnt > 0 {
                    emit(&mut c[1], width, &[20, el as u64, u64::from(cnt)]);
                }
            }
            emit(
                &mut c[4],
                width,
                &[
                    5,
                    s.mass.round() as u64,
                    elem(&m.atoms[root].symbol)? as u64,
                    elem(&m.atoms[other].symbol)? as u64,
                    order(b.order),
                ],
            );
            emit(&mut c[4], width, &[51, (s.mass / 5.0).round() as u64, env]);
            emit(
                &mut c[6],
                width,
                &[
                    7,
                    s.mass.round() as u64,
                    s.basic.min(4) as u64,
                    s.acidic.min(4) as u64,
                ],
            );
            singles.push(s);
        }
    }
    // Deduplicate physical connected atom sets; retain all shortest two-cut
    // pathways as features below. A fragment previously reachable with one cut
    // is not counted as a novel two-cut fragment.
    let mut seen: BTreeSet<Vec<usize>> = singles.iter().map(|s| s.atoms.clone()).collect();
    let mut doubles = vec![];
    for b1 in 0..nb {
        for b2 in b1 + 1..nb {
            let comps = components(&adj, &[b1, b2]);
            if comps.len() < 2 {
                continue;
            }
            for ids in comps {
                if seen.contains(&ids) {
                    continue;
                }
                // A fragment's actual boundary must use both removed bonds.
                let mut boundary = vec![];
                for &bi in &[b1, b2] {
                    let b = &m.bonds[bi];
                    if ids.binary_search(&b.a).is_ok() != ids.binary_search(&b.b).is_ok() {
                        boundary.push(bi);
                    }
                }
                if boundary.len() != 2 {
                    continue;
                }
                let mut ev = vec![];
                for &bi in &boundary {
                    let b = &m.bonds[bi];
                    let (u, v) = if ids.binary_search(&b.a).is_ok() {
                        (b.a, b.b)
                    } else {
                        (b.b, b.a)
                    };
                    ev.push(hash(&[labels[1][u], labels[1][v], order(b.order)]));
                }
                ev.sort_unstable();
                let ring = cycle[b1] > 0 && cycle[b2] > 0;
                let s = frag(&m, ids.clone(), 2, ring, hash(&ev), &g)?;
                emit(&mut c[7], width, &[8, s.mass.round() as u64]);
                if ring {
                    emit(&mut c[8], width, &[9, s.mass.round() as u64]);
                }
                // Parent graph -> one-cut fragment -> two-cut child. Only real
                // subset paths qualify; ring-only two-cut events have no such parent.
                for parent in &singles {
                    if s.atoms
                        .iter()
                        .all(|x| parent.atoms.binary_search(x).is_ok())
                    {
                        emit(
                            &mut c[9],
                            width,
                            &[
                                10,
                                parent.mass.round() as u64,
                                s.mass.round() as u64,
                                (parent.mass - s.mass).round() as u64,
                            ],
                        );
                    }
                }
                seen.insert(ids);
                doubles.push(s);
            }
        }
    }
    for (i, gs) in g.iter().enumerate() {
        for &tag in gs {
            emit(&mut c[10], width, &[11, tag]);
        }
        let ringdegree = adj[i].iter().filter(|&&(_, bi, _)| cycle[bi] > 0).count();
        emit(
            &mut c[13],
            width,
            &[140, ringdegree as u64, m.atoms[i].aromatic as u64],
        );
        let multiple = adj[i].iter().filter(|x| x.2 > 1).count();
        let adjacent_pi = adj[i]
            .iter()
            .filter(|&&(j, _, _)| adj[j].iter().any(|x| x.2 > 1))
            .count();
        emit(
            &mut c[14],
            width,
            &[
                15,
                elem(&m.atoms[i].symbol)? as u64,
                adj[i].len() as u64,
                multiple as u64,
                adjacent_pi as u64,
                m.atoms[i].aromatic as u64,
            ],
        );
        if m.atoms[i].symbol != "C" && m.atoms[i].symbol != "H" {
            emit(&mut c[12], width, &[13, labels[1][i]]);
        }
        for j in i + 1..n {
            for &u in gs {
                for &v in &g[j] {
                    emit(
                        &mut c[11],
                        width,
                        &[12, u.min(v), u.max(v), ds[i][j] as u64],
                    );
                }
            }
            if m.atoms[i].symbol != "C" && m.atoms[j].symbol != "C" {
                let mut p = [
                    elem(&m.atoms[i].symbol)? as u64,
                    elem(&m.atoms[j].symbol)? as u64,
                ];
                p.sort_unstable();
                emit(&mut c[12], width, &[131, p[0], p[1], ds[i][j] as u64]);
            }
            let mut p = [labels[0][i], labels[0][j]];
            p.sort_unstable();
            emit(&mut c[17], width, &[18, p[0], p[1], ds[i][j] as u64]);
        }
    }
    // Loss availability is gated by local functional group and parent formula.
    // Each feature is a hypothesis, not a claim that the loss occurs.
    let has = |tag| g.iter().any(|x| x.contains(&tag));
    let losses = [
        (1, 18.010564684, has(1) && f[0] >= 2 && f[3] >= 1),
        (2, 17.026549101, has(4) && f[0] >= 3 && f[2] >= 1),
        (3, 27.994914620, has(8) && f[1] >= 1 && f[3] >= 1),
        (4, 43.989829240, has(2) && f[1] >= 1 && f[3] >= 2),
    ];
    for (tag, loss, valid) in losses {
        if valid {
            emit(&mut c[5], width, &[6, tag]);
            emit(
                &mut c[5],
                width,
                &[61, tag, (totalmass - loss).round() as u64],
            );
        }
    }
    for (el, &cnt) in f.iter().enumerate() {
        emit(&mut c[15], width, &[16, el as u64, u64::from(cnt)]);
    }
    let dbe = 1.0
        + f[1] as f64
        + (f[2] as f64 - f[0] as f64 - f[6..10].iter().map(|&x| x as f64).sum::<f64>()) / 2.0;
    let dbe_valid = f[5] == 0 && f[10] == 0 && f[11] == 0;
    let props = vec![
        totalmass,
        crate::tpsa::calc_tpsa(&m),
        crate::logp_crippen::calc_logp(&m),
        crate::num_h_donors(&m) as f64,
        crate::num_h_acceptors(&m) as f64,
        crate::num_rotatable_bonds(&m) as f64,
        crate::fraction_csp3(&m),
        n as f64,
        nb as f64,
        if dbe_valid { dbe } else { 0.0 },
        f64::from(dbe_valid),
        m.atoms.iter().filter(|a| a.aromatic).count() as f64,
        cycle.iter().filter(|&&x| x > 0).count() as f64,
    ];
    for (i, &v) in props.iter().enumerate() {
        if !v.is_finite() {
            return Err("nonfinite native descriptor".into());
        }
        emit(
            &mut c[16],
            width,
            &[
                17,
                i as u64,
                ((v * 1e6).round() / 1e6 * 2.0).round() as i64 as u64,
            ],
        );
    }
    let mut fragments = singles;
    fragments.extend(doubles);
    Ok(Features {
        channels: c.into_iter().map(|x| x.into_iter().collect()).collect(),
        properties: props,
        fragments,
        formula: f,
        mass: totalmass,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cut_mass_is_conserved() {
        let a = calculate("CCO", 2048).unwrap();
        assert!((a.mass - 46.0418648).abs() < 1e-6);
        assert_eq!(a.fragments.iter().filter(|x| x.cuts == 1).count(), 4);
        for s in a.fragments.iter().filter(|x| x.cuts == 1) {
            assert!(a
                .fragments
                .iter()
                .any(|t| t.cuts == 1 && (s.mass + t.mass - a.mass).abs() < 1e-8));
        }
    }
    #[test]
    fn rings_need_two_cuts() {
        let a = calculate("C1CCCCC1", 2048).unwrap();
        assert!(a.channels[0].is_empty());
        assert!(!a.channels[8].is_empty());
        assert!(a.fragments.iter().all(|x| x.cuts == 2 && x.ring));
        assert!(a.channels[9].is_empty());
    }
    #[test]
    fn atom_order_invariant() {
        for (a, b) in [
            ("CCOC(=O)NCC", "CCNC(=O)OCC"),
            ("Oc1ccccc1C", "Cc1ccccc1O"),
            ("C1CCCCC1", "C1(CCCCC1)"),
        ] {
            let a = calculate(a, 2048).unwrap();
            let b = calculate(b, 2048).unwrap();
            assert_eq!(a.channels, b.channels);
            assert_eq!(a.formula, b.formula);
        }
    }
    #[test]
    fn descriptor_bin_boundary_is_order_invariant() {
        let a = calculate("COCC1(C(N)=O)CCCN1C(=O)c1ccc(OC(C)C)nc1", 2048).unwrap();
        let b = calculate("c1(ccc(nc1)OC(C)C)C(=O)N1CCCC1(COC)C(=O)N", 2048).unwrap();
        assert_eq!(a.channels[16], b.channels[16]);
    }
    #[test]
    fn invalid_scope_fails_explicitly() {
        for s in ["[13CH4]", "CC.O", "[NH4+]", "not smiles"] {
            assert!(calculate(s, 2048).is_err());
        }
    }
    #[test]
    fn actual_boundary_prevents_duplicate_cuts() {
        let a = calculate("CCCC", 2048).unwrap();
        let mut ids = BTreeSet::new();
        for f in &a.fragments {
            assert!(ids.insert(f.atoms.clone()));
        }
        assert_eq!(a.fragments.iter().filter(|f| f.cuts == 2).count(), 3);
    }
}
