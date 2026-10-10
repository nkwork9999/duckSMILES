//! Fifty separately defined paper-inspired research descriptors (DM004).
//! These are adaptations/novel combinations, not fifty published algorithms.
//! See the experiment registry for sources, definitions and applicability.
use crate::ms_features::{self as ms, adjacency, distances, elem, groups, hash, order};
use crate::parser::parse;
use std::collections::{BTreeMap, BTreeSet};
pub const SPEC: &str = "duck-paper50-v2";
pub const NAMES: [&str; 50] = [
    "fragment_mass_defect",
    "fragment_kendrick_ch2",
    "loss_mass_defect",
    "mass_nitrogen_joint",
    "mass_oxygen_joint",
    "mass_sulfur_phosphorus_joint",
    "fragment_halogen_class",
    "fragment_dbe_profile",
    "fragment_van_krevelen",
    "element_retention",
    "cut_size_balance",
    "boundary_bond_orders",
    "boundary_element_pairs",
    "boundary_degrees",
    "boundary_aromaticity",
    "basic_retention",
    "acidic_retention",
    "charge_site_cut_distance",
    "donor_cut_distance",
    "acceptor_proxy_cut_distance",
    "fragment_cycle_rank",
    "fragment_aromatic_fraction",
    "hetero_distance_inverse",
    "formal_charge_separation",
    "fragment_branch_fraction",
    "boundary_branch_asymmetry",
    "fragment_environment_entropy",
    "equivalent_cleavage_multiplicity",
    "formula_degeneracy",
    "nominal_mass_formula_ambiguity",
    "path_loss_formula",
    "path_loss_mass_defect",
    "path_dbe_change",
    "path_hetero_retention",
    "two_cut_distance",
    "ring_cut_cycle_sizes",
    "boundary_acid_base_pair",
    "acetal_link_mass",
    "carbonyl_hetero_link_mass",
    "sulfonyl_phosphoryl_link_mass",
    "fragment_water_loss",
    "fragment_ammonia_loss",
    "fragment_co_loss",
    "fragment_co2_loss",
    "double_water_loss",
    "fragment_complement_formula",
    "nonlocal_environment_mean_distance",
    "aromatic_hetero_distance",
    "functional_triangles",
    "fragment_mass_entropy",
];
type Counts = BTreeMap<usize, u32>;
fn emit(c: &mut [Counts], i: usize, w: usize, v: &[u64]) {
    let mut x = vec![1000 + i as u64];
    x.extend_from_slice(v);
    *c[i].entry((hash(&x) % w as u64) as usize).or_default() += 1;
}
fn quant(v: f64, scale: f64) -> u64 {
    ((v * 1e7).round() / 1e7 * scale).round() as i64 as u64
}
fn dbe2(f: &[u16; 13]) -> Option<i64> {
    if f[5] + f[10] + f[11] > 0 {
        None
    } else {
        Some(
            2 + 2 * i64::from(f[1]) + i64::from(f[2])
                - i64::from(f[0])
                - f[6..10].iter().map(|&x| i64::from(x)).sum::<i64>(),
        )
    }
}
fn retained_acid(
    mol: &crate::parser::Molecule,
    a: &ms::Adj,
    g: &[Vec<u64>],
    inside: &[bool],
) -> bool {
    (0..inside.len()).any(|i| {
        inside[i]
            && g[i].contains(&2)
            && a[i].iter().any(|&(j, _, _)| {
                inside[j]
                    && g[j].contains(&8)
                    && a[j]
                        .iter()
                        .any(|&(k, _, o)| inside[k] && o == 2 && mol.atoms[k].symbol == "O")
            })
    })
}
fn entropy(counts: impl Iterator<Item = u32>) -> f64 {
    let x: Vec<_> = counts.collect();
    let total = x.iter().map(|&x| x as f64).sum::<f64>();
    if total == 0.0 {
        return 0.0;
    }
    x.iter()
        .map(|&v| {
            let p = v as f64 / total;
            -p * p.ln()
        })
        .sum()
}
/// Supported domain exactly matches ms18-v1. All graph cuts are exhaustively
/// bounded by that API. Width controls hashing only, not enumeration.
pub fn calculate(smiles: &str, width: usize) -> Result<Vec<Vec<(usize, u32)>>, String> {
    if !(64..=65536).contains(&width) {
        return Err("width outside64..65536".into());
    }
    let base = ms::calculate(smiles, 2048)?;
    let mol = parse(smiles).ok_or("parse failed")?;
    let a = adjacency(&mol);
    let n = a.len();
    let g = groups(&mol, &a);
    let ds: Vec<_> = (0..n).map(|i| distances(&a, i, &[])).collect();
    let mut c = vec![Counts::new(); 50];
    let labels: Vec<_> = (0..n)
        .map(|i| {
            let mut neighbors: Vec<_> = a[i]
                .iter()
                .map(|&(j, _, o)| {
                    hash(&[
                        elem(&mol.atoms[j].symbol).unwrap() as u64,
                        o,
                        mol.atoms[j].aromatic as u64,
                    ])
                })
                .collect();
            neighbors.sort_unstable();
            let mut v = vec![
                elem(&mol.atoms[i].symbol).unwrap() as u64,
                mol.atoms[i].hydrogen as u64,
                mol.atoms[i].charge as u64,
            ];
            v.extend(neighbors);
            hash(&v)
        })
        .collect();
    let cycle: Vec<_> = mol
        .bonds
        .iter()
        .enumerate()
        .map(|(i, b)| {
            let d = distances(&a, b.a, &[i])[b.b];
            if d == usize::MAX {
                0
            } else {
                d + 1
            }
        })
        .collect();
    let mut formulas: BTreeMap<[u16; 13], u32> = BTreeMap::new();
    let mut masses: BTreeMap<i64, u32> = BTreeMap::new();
    let mut mf: BTreeMap<i64, BTreeSet<[u16; 13]>> = BTreeMap::new();
    let mut envcount: BTreeMap<u64, u32> = BTreeMap::new();
    let parentbasic = (0..n).filter(|&i| ms::basic(&g[i], &mol, i)).count();
    let parentacid = (0..n).filter(|&i| ms::acidic(&g[i])).count();
    for f in &base.fragments {
        *formulas.entry(f.formula).or_default() += 1;
        *masses.entry(f.mass.round() as i64).or_default() += 1;
        mf.entry(f.mass.round() as i64)
            .or_default()
            .insert(f.formula);
        *envcount.entry(f.environment).or_default() += 1;
        let mass = f.mass.round() as u64;
        let defect = f.mass - f.mass.round();
        let km = f.mass * 14.0 / (12.0 + 2.0 * ms::H);
        let loss = base.mass - f.mass;
        emit(&mut c, 0, width, &[mass / 20, quant(defect, 1000.)]);
        emit(&mut c, 1, width, &[quant(km.round() - km, 1000.)]);
        emit(
            &mut c,
            2,
            width,
            &[quant(loss - loss.round(), 1000.), f.cuts as u64],
        );
        emit(&mut c, 3, width, &[mass, f.formula[2] as u64]);
        emit(&mut c, 4, width, &[mass, f.formula[3] as u64]);
        emit(
            &mut c,
            5,
            width,
            &[mass, f.formula[4] as u64, f.formula[5] as u64],
        );
        emit(
            &mut c,
            6,
            width,
            &[
                mass / 10,
                f.formula[6] as u64,
                f.formula[7] as u64,
                f.formula[8] as u64,
                f.formula[9] as u64,
            ],
        );
        if let Some(v) = dbe2(&f.formula) {
            emit(&mut c, 7, width, &[mass / 10, v as u64]);
        }
        if f.formula[1] > 0 {
            emit(
                &mut c,
                8,
                width,
                &[
                    quant(f.formula[0] as f64 / f.formula[1] as f64, 10.),
                    quant(f.formula[3] as f64 / f.formula[1] as f64, 10.),
                ],
            );
        }
        for el in 0..13 {
            if base.formula[el] > 0 {
                emit(
                    &mut c,
                    9,
                    width,
                    &[
                        el as u64,
                        quant(f.formula[el] as f64 / base.formula[el] as f64, 10.),
                        f.cuts as u64,
                    ],
                );
            }
        }
        emit(
            &mut c,
            10,
            width,
            &[quant(f.atoms.len() as f64 / n as f64, 10.), f.cuts as u64],
        );
        let mut inside = vec![false; n];
        for &i in &f.atoms {
            inside[i] = true;
        }
        let mut boundary = vec![];
        let mut internal = 0;
        let mut aromatic = 0;
        let mut branch = 0;
        let mut environments: BTreeMap<u64, u32> = BTreeMap::new();
        for &i in &f.atoms {
            aromatic += usize::from(mol.atoms[i].aromatic);
            branch += usize::from(a[i].iter().filter(|&&(j, _, _)| inside[j]).count() > 2);
            *environments.entry(labels[i]).or_default() += 1;
        }
        for (bi, b) in mol.bonds.iter().enumerate() {
            if inside[b.a] && inside[b.b] {
                internal += 1;
            } else if inside[b.a] != inside[b.b] {
                boundary.push((
                    bi,
                    if inside[b.a] { b.a } else { b.b },
                    if inside[b.a] { b.b } else { b.a },
                ));
            }
        }
        let mut orders: Vec<_> = boundary
            .iter()
            .map(|&(bi, _, _)| order(mol.bonds[bi].order))
            .collect();
        orders.sort_unstable();
        emit(&mut c, 11, width, &orders);
        let mut pairs = vec![];
        let mut degrees = vec![];
        let mut arom = vec![];
        for &(bi, u, v) in &boundary {
            pairs.push(hash(&[
                elem(&mol.atoms[u].symbol)? as u64,
                elem(&mol.atoms[v].symbol)? as u64,
                order(mol.bonds[bi].order),
            ]));
            degrees.push(hash(&[a[u].len() as u64, a[v].len() as u64]));
            arom.push(hash(&[
                mol.atoms[u].aromatic as u64,
                mol.atoms[v].aromatic as u64,
            ]));
            emit(
                &mut c,
                25,
                width,
                &[a[u].len() as u64, a[v].len() as u64, mass / 10],
            );
            emit(
                &mut c,
                36,
                width,
                &[
                    ms::basic(&g[u], &mol, u) as u64,
                    ms::acidic(&g[u]) as u64,
                    ms::basic(&g[v], &mol, v) as u64,
                    ms::acidic(&g[v]) as u64,
                    mass / 10,
                ],
            );
            let acetal = |i: usize| {
                mol.atoms[i].symbol == "C"
                    && a[i]
                        .iter()
                        .filter(|&&(j, _, o)| o == 1 && mol.atoms[j].symbol == "O")
                        .count()
                        >= 2
            };
            if (acetal(u) && mol.atoms[v].symbol == "O")
                || (acetal(v) && mol.atoms[u].symbol == "O")
            {
                emit(&mut c, 37, width, &[mass, f.cuts as u64]);
            }
            if (g[u].contains(&8) && matches!(mol.atoms[v].symbol.as_str(), "O" | "N"))
                || (g[v].contains(&8) && matches!(mol.atoms[u].symbol.as_str(), "O" | "N"))
            {
                emit(
                    &mut c,
                    38,
                    width,
                    &[
                        mass,
                        elem(&mol.atoms[u].symbol)? as u64,
                        elem(&mol.atoms[v].symbol)? as u64,
                    ],
                );
            }
            if g[u].iter().chain(g[v].iter()).any(|x| [11, 12].contains(x)) {
                emit(
                    &mut c,
                    39,
                    width,
                    &[
                        mass,
                        elem(&mol.atoms[u].symbol)? as u64,
                        elem(&mol.atoms[v].symbol)? as u64,
                    ],
                );
            }
        }
        pairs.sort_unstable();
        degrees.sort_unstable();
        arom.sort_unstable();
        emit(&mut c, 12, width, &pairs);
        emit(&mut c, 13, width, &degrees);
        emit(&mut c, 14, width, &arom);
        emit(
            &mut c,
            15,
            width,
            &[
                f.basic as u64,
                parentbasic as u64 - f.basic as u64,
                mass / 10,
            ],
        );
        emit(
            &mut c,
            16,
            width,
            &[
                f.acidic as u64,
                parentacid as u64 - f.acidic as u64,
                mass / 10,
            ],
        );
        for &(channel, kind) in &[(17, 0), (18, 1), (19, 2)] {
            let mut dd = vec![];
            for &i in &f.atoms {
                let valid = match kind {
                    0 => ms::basic(&g[i], &mol, i) || ms::acidic(&g[i]),
                    1 => {
                        matches!(mol.atoms[i].symbol.as_str(), "N" | "O" | "S")
                            && mol.atoms[i].hydrogen > 0
                    }
                    _ => crate::is_h_acceptor_atom(&mol, i),
                };
                if valid {
                    if let Some(d) = boundary.iter().map(|&(_, u, _)| ds[i][u]).min() {
                        dd.push(d);
                    }
                }
            }
            for d in dd {
                emit(&mut c, channel, width, &[d as u64, mass / 20]);
            }
        }
        emit(
            &mut c,
            20,
            width,
            &[(internal + 1 - f.atoms.len()) as u64, mass / 10],
        );
        emit(
            &mut c,
            21,
            width,
            &[
                quant(aromatic as f64 / f.atoms.len() as f64, 10.),
                mass / 10,
            ],
        );
        for (ii, &i) in f.atoms.iter().enumerate() {
            for &j in f.atoms.iter().skip(ii + 1) {
                if !matches!(mol.atoms[i].symbol.as_str(), "C" | "H")
                    && !matches!(mol.atoms[j].symbol.as_str(), "C" | "H")
                {
                    let mut p = [
                        elem(&mol.atoms[i].symbol)? as u64,
                        elem(&mol.atoms[j].symbol)? as u64,
                    ];
                    p.sort_unstable();
                    emit(
                        &mut c,
                        22,
                        width,
                        &[p[0], p[1], quant(1. / ds[i][j] as f64, 100.)],
                    );
                }
                if mol.atoms[i].charge != 0 && mol.atoms[j].charge != 0 {
                    emit(
                        &mut c,
                        23,
                        width,
                        &[
                            (mol.atoms[i].charge * mol.atoms[j].charge) as u64,
                            ds[i][j] as u64,
                        ],
                    );
                }
            }
        }
        emit(
            &mut c,
            24,
            width,
            &[quant(branch as f64 / f.atoms.len() as f64, 10.), mass / 10],
        );
        emit(
            &mut c,
            26,
            width,
            &[
                quant(entropy(environments.values().copied()), 10.),
                f.cuts as u64,
            ],
        );
        if boundary.len() == 2 {
            let (b1, u1, _) = boundary[0];
            let (b2, u2, _) = boundary[1];
            emit(&mut c, 34, width, &[ds[u1][u2] as u64, mass / 10]);
            if f.ring {
                let mut sizes = [cycle[b1] as u64, cycle[b2] as u64];
                sizes.sort_unstable();
                emit(&mut c, 35, width, &[sizes[0], sizes[1], mass / 10]);
            }
        }
        let hydroxyl = f
            .atoms
            .iter()
            .filter(|&&i| g[i].contains(&1) && a[i].iter().any(|&(j, _, _)| inside[j]))
            .count();
        let amine = f
            .atoms
            .iter()
            .any(|&i| g[i].contains(&4) && mol.atoms[i].hydrogen > 0);
        let carbonyl = f.atoms.iter().any(|&i| {
            g[i].contains(&8)
                && a[i]
                    .iter()
                    .any(|&(j, _, o)| inside[j] && o == 2 && mol.atoms[j].symbol == "O")
        });
        let acid = retained_acid(&mol, &a, &g, &inside);
        for (chan, loss, valid) in [
            (40, 18.010564684, hydroxyl > 0 && f.formula[0] >= 2),
            (41, 17.026549101, amine && f.formula[0] >= 3),
            (42, 27.994914620, carbonyl),
            (43, 43.989829240, acid && f.formula[3] >= 2),
            (
                44,
                36.021129368,
                hydroxyl >= 2 && f.formula[0] >= 4 && f.formula[3] >= 2,
            ),
        ] {
            if valid && f.mass > loss {
                emit(
                    &mut c,
                    chan,
                    width,
                    &[(f.mass - loss).round() as u64, f.cuts as u64],
                );
            }
        }
        let mut complementary: Vec<u64> = f.formula.iter().map(|&x| x as u64).collect();
        complementary.extend((0..13).map(|i| (base.formula[i] - f.formula[i]) as u64));
        emit(&mut c, 45, width, &complementary);
    }
    for (env, count) in envcount {
        emit(&mut c, 27, width, &[env, count as u64]);
    }
    for (formula, count) in formulas {
        let mut v = vec![count as u64];
        v.extend(formula.iter().map(|&x| x as u64));
        emit(&mut c, 28, width, &v);
    }
    for (mass, formulas) in mf {
        emit(&mut c, 29, width, &[mass as u64, formulas.len() as u64]);
    }
    let singles: Vec<_> = base.fragments.iter().filter(|f| f.cuts == 1).collect();
    for child in base.fragments.iter().filter(|f| f.cuts == 2) {
        for parent in &singles {
            if child
                .atoms
                .iter()
                .all(|i| parent.atoms.binary_search(i).is_ok())
            {
                let loss: Vec<_> = (0..13)
                    .map(|i| (parent.formula[i] - child.formula[i]) as u64)
                    .collect();
                emit(&mut c, 30, width, &loss);
                let lm = parent.mass - child.mass;
                emit(
                    &mut c,
                    31,
                    width,
                    &[quant(lm - lm.round(), 1000.), lm.round() as u64 / 10],
                );
                if let (Some(a), Some(b)) = (dbe2(&parent.formula), dbe2(&child.formula)) {
                    emit(
                        &mut c,
                        32,
                        width,
                        &[(a - b) as u64, child.mass.round() as u64 / 10],
                    );
                }
                for el in [2, 3, 4, 5] {
                    emit(
                        &mut c,
                        33,
                        width,
                        &[
                            el as u64,
                            parent.formula[el] as u64,
                            child.formula[el] as u64,
                        ],
                    );
                }
            }
        }
    }
    let neighbors: Vec<BTreeSet<usize>> = (0..n)
        .map(|i| std::iter::once(i).chain(a[i].iter().map(|x| x.0)).collect())
        .collect();
    for i in 0..n {
        for j in i + 1..n {
            if !neighbors[i].is_subset(&neighbors[j]) && !neighbors[j].is_subset(&neighbors[i]) {
                let sum: usize = neighbors[i]
                    .iter()
                    .map(|&u| neighbors[j].iter().map(|&v| ds[u][v]).sum::<usize>())
                    .sum();
                let count = neighbors[i].len() * neighbors[j].len();
                let mut ls = [labels[i], labels[j]];
                ls.sort_unstable();
                emit(
                    &mut c,
                    46,
                    width,
                    &[ls[0], ls[1], ((8 * sum + count) / (2 * count)) as u64],
                );
            }
            for (u, v) in [(i, j), (j, i)] {
                if mol.atoms[u].aromatic && !matches!(mol.atoms[v].symbol.as_str(), "C" | "H") {
                    emit(
                        &mut c,
                        47,
                        width,
                        &[
                            labels[u],
                            elem(&mol.atoms[v].symbol)? as u64,
                            ds[u][v] as u64,
                        ],
                    );
                }
            }
        }
    }
    let anchors: Vec<_> = (0..n).filter(|&i| !g[i].is_empty()).collect();
    for x in 0..anchors.len() {
        for y in x + 1..anchors.len() {
            for z in y + 1..anchors.len() {
                let v = [anchors[x], anchors[y], anchors[z]];
                if ds[v[0]][v[1]].max(ds[v[1]][v[2]]).max(ds[v[0]][v[2]]) > 8 {
                    continue;
                }
                let mut codes = vec![];
                for p in [
                    [0, 1, 2],
                    [0, 2, 1],
                    [1, 0, 2],
                    [1, 2, 0],
                    [2, 0, 1],
                    [2, 1, 0],
                ] {
                    let [i, j, k] = p.map(|p| v[p]);
                    codes.push([
                        g[i][0],
                        g[j][0],
                        g[k][0],
                        ds[i][j] as u64,
                        ds[j][k] as u64,
                        ds[i][k] as u64,
                    ]);
                }
                codes.sort_unstable();
                emit(&mut c, 48, width, &codes[0]);
            }
        }
    }
    emit(
        &mut c,
        49,
        width,
        &[
            quant(entropy(masses.values().copied()), 10.),
            base.fragments.len() as u64,
            masses.len() as u64,
        ],
    );
    Ok(c.into_iter().map(|x| x.into_iter().collect()).collect())
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn co2_loss_requires_own_carbonyl_oxygen() {
        let mol = parse("O=C(O)CCO").unwrap();
        let a = adjacency(&mol);
        let g = groups(&mol, &a);
        let mut inside = vec![true; mol.atoms.len()];
        assert!(retained_acid(&mol, &a, &g, &inside));
        inside[0] = false;
        assert!(!retained_acid(&mol, &a, &g, &inside));
    }
    #[test]
    fn fifty_defined_channels() {
        let a = calculate("CCOC(=O)NCC", 1024).unwrap();
        assert_eq!(a.len(), 50);
        assert!(a.iter().flatten().all(|&(i, n)| i < 1024 && n > 0));
    }
    #[test]
    fn atom_order() {
        for (a, b) in [
            ("CCOC(=O)NCC", "CCNC(=O)OCC"),
            ("Oc1ccccc1C", "Cc1ccccc1O"),
            ("CC1(CO)CCN1C(=O)C", "OCC1(C)CCN1C(C)=O"),
        ] {
            assert_eq!(calculate(a, 1024).unwrap(), calculate(b, 1024).unwrap());
        }
    }
    #[test]
    fn loss_requires_functionality() {
        let a = calculate("CCCC", 1024).unwrap();
        for c in &a[40..45] {
            assert!(c.is_empty());
        }
    }
    #[test]
    fn known_kendrick_homologues() {
        let a = 14.01565006446 * 14. / (12. + 2. * ms::H);
        assert!((a - 14.).abs() < 1e-9);
    }
}
