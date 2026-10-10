//! Versioned research fingerprints, independent of the stable Morgan C ABI.
//! Counts retain repeated atom environments; FNV-1a u64 hashes are deterministic,
//! not collision-free. No RDKit bit compatibility or stereochemical invariance
//! beyond the parser's non-stereo graph is claimed. @/@@ are deliberately omitted.
use crate::parser::{parse, BondOrder, Molecule};
use std::collections::{BTreeMap, VecDeque};

pub const SPEC: &str = "duck-morgan-casmi-v1";
pub const SPEC_V2: &str = "duck-morgan-casmi-v2-local-plus-context";
pub const VARIANTS: [&str; 6] = ["count2", "count3", "cycle2", "shell2", "cut2", "hybrid2"];
type Adj = Vec<Vec<(usize, usize, u64)>>;

fn hash(values: &[u64]) -> u64 {
    let mut h = 14695981039346656037u64;
    for v in values {
        for b in v.to_le_bytes() {
            h = (h ^ u64::from(b)).wrapping_mul(1099511628211);
        }
    }
    h
}
fn symbol(s: &str) -> u64 {
    hash(&s.bytes().map(u64::from).collect::<Vec<_>>())
}
fn distances(adj: &Adj, start: usize, blocked: Option<usize>) -> Vec<usize> {
    let mut d = vec![usize::MAX; adj.len()];
    d[start] = 0;
    let mut q = VecDeque::from([start]);
    while let Some(u) = q.pop_front() {
        for &(v, bi, _) in &adj[u] {
            if Some(bi) != blocked && d[v] == usize::MAX {
                d[v] = d[u] + 1;
                q.push_back(v);
            }
        }
    }
    d
}
fn adjacency(m: &Molecule) -> Adj {
    let mut a = vec![vec![]; m.atoms.len()];
    for (i, b) in m.bonds.iter().enumerate() {
        let o = match b.order {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Aromatic => 4,
        };
        a[b.a].push((b.b, i, o));
        a[b.b].push((b.a, i, o));
    }
    a
}

/// V2 keeps the local radius-3 vector unchanged in [0,width), and emits
/// independently pooled context in [width,2*width). Unlike V1, a remote change
/// cannot replace every local atom label. Channels must be normalized separately.
/// Four contexts: edge-cycle sizes, atom-pair distance shells, coarse bridge-side
/// composition, and bridge-side neutral atomic masses (NOT predicted ion peaks).
pub fn fingerprints_v2(smiles: &str, width: usize) -> Result<Vec<Vec<(usize, u32)>>, String> {
    let base = fingerprints(smiles, width)?[1].clone();
    let mol = parse(smiles).ok_or("parse failed")?;
    let adj = adjacency(&mol);
    let mut ctx: Vec<BTreeMap<usize, u32>> = vec![BTreeMap::new(); 4];
    fn emit(m: &mut BTreeMap<usize, u32>, v: &[u64], width: usize) {
        *m.entry(width + (hash(v) % width as u64) as usize)
            .or_default() += 1;
    }
    for (i, a) in mol.atoms.iter().enumerate() {
        let d = distances(&adj, i, None);
        for (j, b) in mol.atoms.iter().enumerate().skip(i + 1) {
            if d[j] != usize::MAX {
                let mut pair = [symbol(&a.symbol), symbol(&b.symbol)];
                pair.sort_unstable();
                emit(
                    &mut ctx[1],
                    &[601, pair[0], pair[1], d[j].min(12) as u64],
                    width,
                );
            }
        }
    }
    for (bi, b) in mol.bonds.iter().enumerate() {
        let d = distances(&adj, b.a, Some(bi));
        if d[b.b] != usize::MAX {
            let mut pair = [
                symbol(&mol.atoms[b.a].symbol),
                symbol(&mol.atoms[b.b].symbol),
            ];
            pair.sort_unstable();
            emit(
                &mut ctx[0],
                &[602, pair[0], pair[1], (d[b.b] + 1) as u64],
                width,
            );
        } else {
            for root in [b.a, b.b] {
                let side = distances(&adj, root, Some(bi));
                let mut counts = [0u64; 6];
                let mut mass = Some(0.0);
                for (i, a) in mol.atoms.iter().enumerate() {
                    if side[i] == usize::MAX {
                        continue;
                    }
                    let k = match a.symbol.as_str() {
                        "C" => 0,
                        "N" => 1,
                        "O" => 2,
                        "S" => 3,
                        "P" => 4,
                        _ => 5,
                    };
                    counts[k] += 1;
                    mass = mass.and_then(|x| {
                        if a.isotope.is_some() {
                            None
                        } else {
                            crate::weights::monoisotopic_mass(&a.symbol)
                                .map(|m| x + m + f64::from(a.hydrogen) * 1.00782503)
                        }
                    });
                }
                for (element, count) in counts.iter().enumerate() {
                    if *count > 0 {
                        emit(
                            &mut ctx[2],
                            &[
                                603,
                                element as u64,
                                (*count).min(16),
                                symbol(&mol.atoms[root].symbol),
                            ],
                            width,
                        );
                    }
                }
                if let Some(m) = mass {
                    if (1.0..1000.0).contains(&m) {
                        emit(&mut ctx[3], &[604, m.round() as u64], width);
                    }
                }
            }
        }
    }
    Ok(ctx
        .into_iter()
        .map(|m| base.iter().copied().chain(m).collect())
        .collect())
}

/// Compute all variants together. Limit work explicitly; unsupported inputs fail
/// rather than silently receiving a zero fingerprint. Width 16..=65536.
pub fn fingerprints(smiles: &str, width: usize) -> Result<Vec<Vec<(usize, u32)>>, String> {
    if !(16..=65536).contains(&width) {
        return Err("width outside 16..65536".into());
    }
    let mol = parse(smiles).ok_or("parse failed")?;
    let n = mol.atoms.len();
    if n == 0 || n > 256 {
        return Err("atom count outside 1..256".into());
    }
    let adj = adjacency(&mol);
    let all_d: Vec<_> = (0..n).map(|i| distances(&adj, i, None)).collect();
    let mut cycles = vec![vec![]; n];
    let mut cuts = vec![vec![]; n];
    for (bi, b) in mol.bonds.iter().enumerate() {
        let d = distances(&adj, b.a, Some(bi));
        if d[b.b] != usize::MAX {
            // Shortest cycle through this edge; independent of SSSR choices.
            cycles[b.a].push((d[b.b] + 1) as u64);
            cycles[b.b].push((d[b.b] + 1) as u64);
        } else {
            // Each bridge-side composition, excluding disconnected salts.
            for (root, other) in [(b.a, b.b), (b.b, b.a)] {
                let side = distances(&adj, root, Some(bi));
                let mut composition: BTreeMap<(u64, u64), u64> = BTreeMap::new();
                let mut h = 0i64;
                let mut charge = 0i64;
                for (i, a) in mol.atoms.iter().enumerate() {
                    if side[i] != usize::MAX {
                        *composition
                            .entry((symbol(&a.symbol), u64::from(a.isotope.unwrap_or(0))))
                            .or_default() += 1;
                        h += i64::from(a.hydrogen);
                        charge += i64::from(a.charge);
                    }
                }
                let mut v = vec![
                    301,
                    h as u64,
                    charge as u64,
                    symbol(&mol.atoms[other].symbol),
                ];
                for ((s, iso), count) in composition {
                    v.extend([s, iso, count]);
                }
                cuts[root].push(hash(&v));
            }
        }
    }
    let base: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, a)| {
            hash(&[
                101,
                symbol(&a.symbol),
                a.hydrogen as u64,
                a.charge as u64,
                u64::from(a.isotope.unwrap_or(0)),
                u64::from(a.aromatic),
                adj[i].len() as u64,
                u64::from(!cycles[i].is_empty()),
            ])
        })
        .collect();
    let mut result = Vec::new();
    for variant in VARIANTS {
        let mut current = base.clone();
        for i in 0..n {
            let mut v = vec![201, base[i]];
            if ["cycle2", "shell2", "cut2", "hybrid2"].contains(&variant) {
                let mut sizes = cycles[i].clone();
                sizes.sort_unstable();
                v.extend([202, sizes.len() as u64]);
                v.extend(sizes);
            }
            if ["shell2", "hybrid2"].contains(&variant) {
                // Beyond the local Morgan radius: composition per graph-distance
                // shell, capped at 8; disconnected components never contribute.
                let mut shell = Vec::new();
                for (j, &label) in base.iter().enumerate() {
                    if i != j && all_d[i][j] != usize::MAX {
                        shell.push(hash(&[all_d[i][j].min(8) as u64, label]));
                    }
                }
                shell.sort_unstable();
                v.extend([203, shell.len() as u64]);
                v.extend(shell);
            }
            if ["cut2", "hybrid2"].contains(&variant) {
                let mut c = cuts[i].clone();
                c.sort_unstable();
                v.extend([204, c.len() as u64]);
                v.extend(c);
            }
            current[i] = hash(&v);
        }
        let mut features: BTreeMap<usize, u32> = BTreeMap::new();
        let radius = if variant == "count3" { 3 } else { 2 };
        for layer in 0..=radius {
            for &label in &current {
                *features
                    .entry((hash(&[401, layer, label]) % width as u64) as usize)
                    .or_default() += 1;
            }
            let next = (0..n)
                .map(|i| {
                    let mut nbr: Vec<_> = adj[i]
                        .iter()
                        .map(|&(j, _, o)| hash(&[o, current[j]]))
                        .collect();
                    nbr.sort_unstable();
                    let mut v = vec![501, layer, current[i], nbr.len() as u64];
                    v.extend(nbr);
                    hash(&v)
                })
                .collect();
            current = next;
        }
        result.push(features.into_iter().collect());
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn representation_invariant() {
        for (a, b) in [
            ("CCO", "OCC"),
            ("CC(C)O", "OC(C)C"),
            ("c1ccccc1O", "Oc1ccccc1"),
            ("C1CC2CCC1C2", "C1C2CCC(C1)C2"),
            ("CCO.[Na+]", "[Na+].OCC"),
        ] {
            assert_eq!(fingerprints(a, 2048), fingerprints(b, 2048), "{a} / {b}");
        }
    }
    #[test]
    fn counts_and_isotopes_are_retained() {
        assert_ne!(fingerprints("CCO", 2048), fingerprints("[13CH3]CO", 2048));
        let a = fingerprints("C", 2048).unwrap();
        let b = fingerprints("C.C", 2048).unwrap();
        for (x, y) in a.iter().zip(b.iter()) {
            assert_eq!(x.iter().map(|&(i, c)| (i, c * 2)).collect::<Vec<_>>(), *y);
        }
    }
    #[test]
    fn nonlocal_ring_context_resolves_local_ambiguity() {
        let a = fingerprints("C1CCCCCC1", 65536).unwrap();
        let b = fingerprints("C1CCCCCCC1", 65536).unwrap();
        // Local environments identical as sets; explicit cycle size separates.
        assert_eq!(
            a[0].iter().map(|x| x.0).collect::<Vec<_>>(),
            b[0].iter().map(|x| x.0).collect::<Vec<_>>()
        );
        assert_ne!(
            a[2].iter().map(|x| x.0).collect::<Vec<_>>(),
            b[2].iter().map(|x| x.0).collect::<Vec<_>>()
        );
    }
    #[test]
    fn explicit_limits() {
        assert!(fingerprints("CC", 0).is_err());
        assert!(fingerprints("C1CC", 2048).is_err());
        assert!(fingerprints(&"C".repeat(257), 2048).is_err());
    }

    #[test]
    fn v2_preserves_local_channel_and_atom_order() {
        let local = fingerprints("CCOC", 2048).unwrap()[1].clone();
        let a = fingerprints_v2("CCOC", 2048).unwrap();
        assert_eq!(a, fingerprints_v2("COCC", 2048).unwrap());
        for fp in a {
            assert_eq!(
                fp.into_iter().filter(|x| x.0 < 2048).collect::<Vec<_>>(),
                local
            );
        }
    }
}
