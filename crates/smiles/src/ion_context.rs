//! DM006: 24 HR-conditioned graph fingerprints per polarity and eight direct
//! spectrum descriptors. Custom MS-FINDER/FIORA/FraGNNet-inspired definitions,
//! not their trained models. See foundation experiments/DM006/PROTOCOL.md.
use crate::hr_features::{self, Features as HrFeatures};
use crate::ms_features::{self, Fragment, ELEMENTS};
use crate::parser::{parse, Molecule};
use std::collections::{BTreeMap, BTreeSet};

pub const SPEC: &str = "duck-ion-context48-v1";
pub const MATCH_SPEC: &str = "duck-hr-match8-v1";
pub const NAMES: [&str; 24] = [
    "ion_cut_depth",
    "ion_hydrogen_shift",
    "ion_ring_opening",
    "ion_formula_presence",
    "ion_formula_degeneracy",
    "ion_nominal_ambiguity",
    "ion_loss_mass_pair",
    "ion_hetero_composition",
    "ion_formula_dbe2",
    "ion_mass_defect",
    "loss_mass_by_shift",
    "loss_mass_defect",
    "loss_hetero_composition",
    "loss_small_formula",
    "ion_oriented_boundary",
    "ion_boundary_local_state",
    "ion_acid_base_retention",
    "ion_functional_retention",
    "ion_boundary_separation",
    "ion_size_hshift",
    "path_loss_formula",
    "path_shift_transition",
    "path_environment_transition",
    "path_degree",
];
pub const MATCH_NAMES: [&str; 8] = [
    "hr_peak_recall",
    "hr_intensity_recall",
    "hr_mass_precision",
    "hr_peak_f1",
    "hr_single_intensity",
    "hr_double_intensity",
    "hr_ring_intensity",
    "hr_path_intensity",
];

#[derive(Clone, Debug)]
pub struct Features {
    /// Positive24 channels followed by negative24 channels.
    pub channels: Vec<Vec<(usize, u32)>>,
    pub hr: HrFeatures,
    /// Directed parent/child indices in hr.ions, including both polarities.
    pub paths: Vec<(usize, usize)>,
}

fn boundary(m: &Molecule, f: &Fragment) -> Vec<(usize, usize, usize)> {
    m.bonds
        .iter()
        .enumerate()
        .filter_map(|(bi, b)| {
            let a = f.atoms.binary_search(&b.a).is_ok();
            let z = f.atoms.binary_search(&b.b).is_ok();
            (a != z).then_some(if a { (bi, b.a, b.b) } else { (bi, b.b, b.a) })
        })
        .collect()
}

fn paths(m: &Molecule, hr: &HrFeatures) -> Vec<(usize, usize)> {
    let boundaries: Vec<_> = hr.fragments.iter().map(|f| boundary(m, f)).collect();
    let mut by_fragment = vec![Vec::new(); hr.fragments.len()];
    for (i, ion) in hr.ions.iter().enumerate() {
        by_fragment[ion.fragment].push(i);
    }
    let mut edges = Vec::new();
    for (pi, p) in hr.fragments.iter().enumerate() {
        if p.cuts != 1 || by_fragment[pi].is_empty() {
            continue;
        }
        let (first, charged, _) = boundaries[pi][0];
        for (ci, c) in hr.fragments.iter().enumerate() {
            if c.cuts != 2
                || by_fragment[ci].is_empty()
                || c.atoms.len() >= p.atoms.len()
                || c.atoms.binary_search(&charged).is_err()
                || !c.atoms.iter().all(|v| p.atoms.binary_search(v).is_ok())
                || !boundaries[ci].iter().any(|&(bi, _, _)| bi == first)
            {
                continue;
            }
            for &a in &by_fragment[pi] {
                for &b in &by_fragment[ci] {
                    let u = &hr.ions[a];
                    let v = &hr.ions[b];
                    if u.polarity == v.polarity
                        && (u.hydrogen_shift - v.hydrogen_shift).abs() == 1
                        && u.mz > v.mz
                        && u.formula.iter().zip(v.formula).all(|(&x, y)| x >= y)
                    {
                        edges.push((a, b));
                    }
                }
            }
        }
    }
    edges
}

fn neutral_mass(f: &[u16; 13]) -> f64 {
    f.iter()
        .zip(ELEMENTS)
        .map(|(&n, e)| {
            f64::from(n) * crate::weights::monoisotopic_mass(e).expect("HR validated element")
        })
        .sum()
}
fn zigzag(value: i64) -> u64 {
    if value >= 0 {
        (value as u64) * 2
    } else {
        value.unsigned_abs() * 2 - 1
    }
}
fn coarse(mass: f64) -> u64 {
    (mass / 10.).floor() as u64
}
fn defect(mass: f64) -> u64 {
    zigzag(((mass - mass.round()) * 100.).round() as i64)
}
fn emit(c: &mut BTreeMap<usize, u32>, width: usize, words: &[u64]) {
    *c.entry((ms_features::hash(words) % width as u64) as usize)
        .or_default() += 1;
}
fn formula_token(prefix: u64, f: &[u16; 13]) -> Vec<u64> {
    let mut out = vec![prefix];
    out.extend(f.iter().map(|&x| u64::from(x)));
    out
}
fn small_loss(f: &[u16; 13]) -> Option<usize> {
    let motifs = [
        (2, 0, 0, 1, 0, 0),
        (3, 0, 1, 0, 0, 0),
        (0, 1, 0, 1, 0, 0),
        (0, 1, 0, 2, 0, 0),
        (2, 1, 0, 1, 0, 0),
        (2, 0, 0, 0, 1, 0),
        (0, 0, 0, 2, 1, 0),
        (3, 0, 0, 4, 0, 1),
    ];
    if f[6..].iter().any(|&x| x != 0) {
        return None;
    }
    motifs
        .iter()
        .position(|&(h, c, n, o, s, p)| f[..6] == [h, c, n, o, s, p])
}

pub fn calculate(smiles: &str, width: usize) -> Result<Features, String> {
    let hr = hr_features::calculate(smiles, width)?;
    let m = parse(smiles).ok_or("parse failed")?;
    let edges = paths(&m, &hr);
    let adj = ms_features::adjacency(&m);
    let groups = ms_features::groups(&m, &adj);
    let distances: Vec<_> = (0..m.atoms.len())
        .map(|i| ms_features::distances(&adj, i, &[]))
        .collect();
    let boundaries: Vec<_> = hr.fragments.iter().map(|f| boundary(&m, f)).collect();
    let mut parent = [0u16; 13];
    for atom in &m.atoms {
        parent[ms_features::elem(&atom.symbol)?] += 1;
        parent[0] += atom.hydrogen as u16;
    }
    let mut parent_groups = [0u64; 17];
    for tags in &groups {
        for &g in tags {
            parent_groups[g as usize] += 1;
        }
    }
    let parent_basic = (0..m.atoms.len())
        .filter(|&i| ms_features::basic(&groups[i], &m, i))
        .count() as u64;
    let parent_acidic = groups.iter().filter(|g| ms_features::acidic(g)).count() as u64;
    let mut c = vec![BTreeMap::new(); 48];
    let mut degrees = vec![[0u64; 2]; hr.ions.len()];
    for &(a, b) in &edges {
        degrees[a][1] += 1;
        degrees[b][0] += 1;
    }
    for mode in [1i8, -1i8] {
        let offset = if mode == 1 { 0 } else { 24 };
        let mut formulas: BTreeMap<[u16; 13], (f64, u64)> = BTreeMap::new();
        let mut nominal: BTreeMap<u64, BTreeSet<[u16; 13]>> = BTreeMap::new();
        for (ii, ion) in hr
            .ions
            .iter()
            .enumerate()
            .filter(|(_, i)| i.polarity == mode)
        {
            let f = &hr.fragments[ion.fragment];
            let mz = ion.mz;
            let nm = mz.round() as u64;
            let cm = coarse(mz);
            let shift = zigzag(i64::from(ion.hydrogen_shift));
            let bnd = &boundaries[ion.fragment];
            let mut loss = [0u16; 13];
            for k in 0..13 {
                let n = i32::from(parent[k]) + if k == 0 { i32::from(mode) } else { 0 }
                    - i32::from(ion.formula[k]);
                loss[k] = u16::try_from(n).map_err(|_| "negative conserved loss")?;
            }
            let lm = neutral_mass(&loss);
            emit(&mut c[offset], width, &[0, nm, u64::from(f.cuts)]);
            emit(&mut c[offset + 1], width, &[1, nm, shift]);
            emit(
                &mut c[offset + 2],
                width,
                &[2, nm, u64::from(f.ring), u64::from(f.cuts)],
            );
            formulas.entry(ion.formula).or_insert((mz, 0)).1 += 1;
            nominal.entry(nm).or_default().insert(ion.formula);
            emit(&mut c[offset + 6], width, &[6, cm, coarse(lm)]);
            emit(
                &mut c[offset + 7],
                width,
                &[
                    7,
                    cm,
                    ion.formula[2].into(),
                    ion.formula[3].into(),
                    ion.formula[4].into(),
                    ion.formula[5].into(),
                ],
            );
            if [5, 10, 11, 12].iter().all(|&i| ion.formula[i] == 0) {
                let z = &ion.formula;
                let dbe2 = 2 * i64::from(z[1]) + 2 + i64::from(z[2])
                    - i64::from(z[0])
                    - z[6..10].iter().map(|&n| i64::from(n)).sum::<i64>();
                emit(&mut c[offset + 8], width, &[8, cm, zigzag(dbe2)]);
            }
            emit(&mut c[offset + 9], width, &[9, cm, defect(mz)]);
            emit(&mut c[offset + 10], width, &[10, lm.round() as u64, shift]);
            emit(&mut c[offset + 11], width, &[11, coarse(lm), defect(lm)]);
            emit(
                &mut c[offset + 12],
                width,
                &[
                    12,
                    coarse(lm),
                    loss[2].into(),
                    loss[3].into(),
                    loss[4].into(),
                    loss[5].into(),
                ],
            );
            if let Some(tag) = small_loss(&loss) {
                emit(
                    &mut c[offset + 13],
                    width,
                    &[13, tag as u64, u64::from(f.cuts)],
                );
            }
            let mut ep = Vec::new();
            let mut local = Vec::new();
            for &(_, u, v) in bnd {
                let a = &m.atoms[u];
                let b = &m.atoms[v];
                ep.push(ms_features::hash(&[
                    ms_features::elem(&a.symbol)? as u64,
                    ms_features::elem(&b.symbol)? as u64,
                ]));
                local.push(ms_features::hash(&[
                    adj[u].len() as u64,
                    adj[v].len() as u64,
                    u64::from(a.aromatic),
                    u64::from(b.aromatic),
                    a.hydrogen as u64,
                    b.hydrogen as u64,
                ]));
            }
            ep.sort_unstable();
            local.sort_unstable();
            let mut token = vec![14, cm, shift];
            token.extend(ep);
            emit(&mut c[offset + 14], width, &token);
            let mut token = vec![15, cm];
            token.extend(local);
            emit(&mut c[offset + 15], width, &token);
            emit(
                &mut c[offset + 16],
                width,
                &[
                    16,
                    cm,
                    f.basic.into(),
                    parent_basic - u64::from(f.basic),
                    f.acidic.into(),
                    parent_acidic - u64::from(f.acidic),
                ],
            );
            let mut retained = [0u64; 17];
            for &i in &f.atoms {
                for &tag in &groups[i] {
                    retained[tag as usize] += 1;
                }
            }
            for tag in 1..17 {
                if parent_groups[tag] > 0 {
                    emit(
                        &mut c[offset + 17],
                        width,
                        &[
                            17,
                            cm,
                            tag as u64,
                            retained[tag],
                            parent_groups[tag] - retained[tag],
                        ],
                    );
                }
            }
            if bnd.len() == 2 {
                emit(
                    &mut c[offset + 18],
                    width,
                    &[18, cm, distances[bnd[0].1][bnd[1].1] as u64, shift],
                );
            }
            emit(
                &mut c[offset + 19],
                width,
                &[19, cm, (10 * f.atoms.len() / m.atoms.len()) as u64, shift],
            );
            emit(
                &mut c[offset + 23],
                width,
                &[23, cm, degrees[ii][0].min(32), degrees[ii][1].min(32)],
            );
        }
        for (formula, (mz, count)) in formulas {
            emit(&mut c[offset + 3], width, &formula_token(3, &formula));
            emit(
                &mut c[offset + 4],
                width,
                &[4, mz.round() as u64, count.min(32)],
            );
        }
        for (mass, fs) in nominal {
            emit(
                &mut c[offset + 5],
                width,
                &[5, mass, (fs.len() as u64).min(32)],
            );
        }
        for &(ai, bi) in &edges {
            let a = &hr.ions[ai];
            let b = &hr.ions[bi];
            if a.polarity != mode {
                continue;
            }
            let mut loss = [0u16; 13];
            for (i, n) in loss.iter_mut().enumerate() {
                *n = a.formula[i] - b.formula[i];
            }
            emit(&mut c[offset + 20], width, &formula_token(20, &loss));
            emit(
                &mut c[offset + 21],
                width,
                &[
                    21,
                    coarse(a.mz),
                    coarse(b.mz),
                    zigzag(a.hydrogen_shift.into()),
                    zigzag(b.hydrogen_shift.into()),
                ],
            );
            emit(
                &mut c[offset + 22],
                width,
                &[
                    22,
                    hr.fragments[a.fragment].environment,
                    hr.fragments[b.fragment].environment,
                    (a.mz - b.mz).round() as u64,
                ],
            );
        }
    }
    Ok(Features {
        channels: c.into_iter().map(|x| x.into_iter().collect()).collect(),
        hr,
        paths: edges,
    })
}

/// Eight bounded direct evidence descriptors, without empirical weights.
pub fn spectral_scores(
    hr: &HrFeatures,
    edges: &[(usize, usize)],
    mode: i8,
    precursor: f64,
    mz: &[f64],
    intensity: &[f64],
) -> Result<[f64; 8], String> {
    if !matches!(mode, 1 | -1)
        || !precursor.is_finite()
        || precursor <= 0.
        || mz.len() != intensity.len()
    {
        return Err("invalid mode, precursor or spectrum array lengths".into());
    }
    let mut peaks: Vec<_> = mz
        .iter()
        .copied()
        .zip(intensity.iter().copied())
        .filter(|&(m, i)| {
            m.is_finite() && i.is_finite() && m > 0. && i > 0. && m < precursor - 0.02
        })
        .collect();
    peaks.sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.total_cmp(&b.1)));
    let mut merged: Vec<(f64, f64)> = Vec::new();
    for (m, i) in peaks {
        if let Some(last) = merged.last_mut().filter(|v| v.0 == m) {
            last.1 += i;
        } else {
            merged.push((m, i));
        }
    }
    if merged.is_empty() {
        return Ok([0.; 8]);
    }
    if merged.iter().any(|(_, i)| !i.is_finite()) {
        return Err("merged intensity overflow".into());
    }
    let total: f64 = merged.iter().map(|p| p.1.sqrt()).sum();
    if !total.is_finite() {
        return Err("intensity normalization overflow".into());
    }
    let weights: Vec<_> = merged.iter().map(|p| p.1.sqrt() / total).collect();
    let nearest = |mass: f64| -> Option<usize> {
        let right = merged.partition_point(|p| p.0 < mass);
        let candidates = [
            right.checked_sub(1),
            (right < merged.len()).then_some(right),
        ];
        candidates
            .into_iter()
            .flatten()
            .min_by(|&a, &b| {
                (merged[a].0 - mass)
                    .abs()
                    .total_cmp(&(merged[b].0 - mass).abs())
                    .then(a.cmp(&b))
            })
            .filter(|&j| (merged[j].0 - mass).abs() <= 0.005_f64.max(mass * 20e-6))
    };
    let mut explained = vec![[false; 5]; merged.len()];
    let mut matches = vec![None; hr.ions.len()];
    let mut unique: BTreeMap<u64, f64> = BTreeMap::new();
    for (i, ion) in hr.ions.iter().enumerate() {
        if ion.polarity != mode || ion.mz >= precursor - 0.02 {
            continue;
        }
        unique
            .entry((ion.mz * 1e5).round() as u64)
            .and_modify(|m| *m = m.min(ion.mz))
            .or_insert(ion.mz);
        if let Some(j) = nearest(ion.mz) {
            matches[i] = Some(j);
            explained[j][0] = true;
            let f = &hr.fragments[ion.fragment];
            explained[j][usize::from(f.cuts)] = true;
            if f.ring {
                explained[j][3] = true;
            }
        }
    }
    for &(a, b) in edges {
        let (Some(u), Some(v)) = (matches[a], matches[b]) else {
            continue;
        };
        if u != v {
            explained[u][4] = true;
            explained[v][4] = true;
        }
    }
    let recall = explained.iter().filter(|x| x[0]).count() as f64 / merged.len() as f64;
    let precision = if unique.is_empty() {
        0.
    } else {
        unique.values().filter(|&&m| nearest(m).is_some()).count() as f64 / unique.len() as f64
    };
    let covered = |ch: usize| {
        explained
            .iter()
            .zip(&weights)
            .filter(|(x, _)| x[ch])
            .map(|(_, w)| w)
            .sum::<f64>()
    };
    Ok([
        recall,
        covered(0),
        precision,
        if recall + precision > 0. {
            2. * recall * precision / (recall + precision)
        } else {
            0.
        },
        covered(1),
        covered(2),
        covered(3),
        covered(4),
    ])
}

/// Convenience matcher avoiding construction of the48 additional fingerprints.
pub fn match_smiles(
    smiles: &str,
    mode: i8,
    precursor: f64,
    mz: &[f64],
    intensity: &[f64],
) -> Result<[f64; 8], String> {
    let hr = hr_features::calculate(smiles, 2048)?;
    let m = parse(smiles).ok_or("parse failed")?;
    let edges = paths(&m, &hr);
    spectral_scores(&hr, &edges, mode, precursor, mz, intensity)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn invariant_channels_and_path_degrees() {
        let a = calculate("CC(O)CN", 2048).unwrap();
        let b = calculate("NCC(O)C", 2048).unwrap();
        assert_eq!(a.channels, b.channels);
        assert_eq!(a.channels.len(), 48);
        assert!(!a.paths.is_empty());
        assert!(!a.channels[20].is_empty());
        assert!(!a.channels[44].is_empty());
        for &(u, v) in &a.paths {
            assert_eq!(a.hr.ions[u].polarity, a.hr.ions[v].polarity);
            assert!(a.hr.ions[u].mz > a.hr.ions[v].mz);
        }
    }
    #[test]
    fn preserves_hr_domain_and_zero_ring_paths() {
        assert!(calculate("[NH3+]CC(=O)[O-]", 2048).is_err());
        assert!(calculate("C.C", 2048).is_err());
        assert!(calculate("CO", 32).is_err());
        let x = calculate("C1CCCCC1", 2048).unwrap();
        assert!(x.paths.is_empty());
        assert!(x.channels[20].is_empty());
        assert!(!x.channels[2].is_empty());
    }
    #[test]
    fn motifs_are_exact_conserved_loss_not_containment() {
        let mut f = [0; 13];
        f[0] = 2;
        f[3] = 1;
        assert_eq!(small_loss(&f), Some(0));
        f[1] = 1;
        assert_eq!(small_loss(&f), Some(4));
        f[1] = 2;
        assert_eq!(small_loss(&f), None);
        assert_eq!(zigzag(-1), 1);
        assert_eq!(zigzag(1), 2);
    }
    #[test]
    fn matched_coverage_and_duplicate_invariance() {
        let x = calculate("CC(O)CN", 2048).unwrap();
        let i = x.hr.ions.iter().position(|i| i.polarity == 1).unwrap();
        let mass = x.hr.ions[i].mz;
        let a = spectral_scores(&x.hr, &x.paths, 1, 200., &[mass], &[4.]).unwrap();
        let b = spectral_scores(&x.hr, &x.paths, 1, 200., &[mass, mass], &[1., 3.]).unwrap();
        assert_eq!(a, b);
        assert_eq!(a[0], 1.);
        assert_eq!(a[1], 1.);
        assert!(a[2] > 0.);
        let mut duplicated = x.hr.clone();
        duplicated.ions.extend(x.hr.ions.clone());
        let d = spectral_scores(&duplicated, &x.paths, 1, 200., &[mass], &[4.]).unwrap();
        assert_eq!(a, d);
        let noise = spectral_scores(&x.hr, &x.paths, 1, 200., &[mass, 199.], &[1., 1.]).unwrap();
        assert_eq!(noise[0], 0.5);
        assert_eq!(noise[1], 0.5);
    }
    #[test]
    fn matches_handle_empty_invalid_and_precursor() {
        let x = calculate("CN", 2048).unwrap();
        assert_eq!(
            spectral_scores(&x.hr, &x.paths, 1, 30., &[], &[]).unwrap(),
            [0.; 8]
        );
        assert_eq!(
            spectral_scores(
                &x.hr,
                &x.paths,
                1,
                30.,
                &[30., f64::NAN, -1.],
                &[1., 1., 1.]
            )
            .unwrap(),
            [0.; 8]
        );
        assert!(spectral_scores(&x.hr, &x.paths, 1, 30., &[1.], &[]).is_err());
        assert!(spectral_scores(&x.hr, &x.paths, 0, 30., &[], &[]).is_err());
        assert!(spectral_scores(&x.hr, &x.paths, 1, f64::INFINITY, &[], &[]).is_err());
    }
}
