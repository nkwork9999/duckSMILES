//! Prepared-only, atom-index traversal of pinned RDKit Kekulize.cpp.
//! BSD-3-Clause; see RDKIT-DATA-LICENSE.txt. No dummy/query/wedge atoms are
//! admitted by this SMILES entry point. Canonical fallback is not implemented;
//! exhausting the reference's 100 backtracks returns an unsupported row.
use crate::parser::{BondOrder, Molecule, RingInfo};
use std::collections::VecDeque;

pub(crate) fn apply(mol: &mut Molecule, rings: &RingInfo) -> Option<()> {
    let aromatic: Vec<_> = mol
        .bonds
        .iter()
        .map(|b| b.order == BondOrder::Aromatic)
        .collect();
    if !aromatic.iter().any(|&a| a) && !mol.atoms.iter().any(|a| a.aromatic) {
        return Some(());
    }
    let adj = mol.adjacency();
    let z = crate::electronic::atomic_numbers(mol)?;
    let valences: Vec<_> = (0..mol.atoms.len())
        .map(|i| crate::electronic::total_valence(mol, i))
        .collect();
    let mut processed = vec![false; rings.rings.len()];
    for start in 0..rings.rings.len() {
        if processed[start] {
            continue;
        }
        // Fused systems share any bond, unlike the later aromaticity rule.
        let mut queue = vec![start];
        processed[start] = true;
        let mut atoms = vec![false; mol.atoms.len()];
        while let Some(ri) = queue.pop() {
            for &bi in &rings.rings[ri] {
                atoms[mol.bonds[bi].a] = true;
                atoms[mol.bonds[bi].b] = true;
            }
            for (rj, other) in rings.rings.iter().enumerate() {
                if !processed[rj] && other.iter().any(|b| rings.rings[ri].contains(b)) {
                    processed[rj] = true;
                    queue.push(rj);
                }
            }
        }
        let members: Vec<_> = (0..atoms.len()).filter(|&i| atoms[i]).collect();
        let is_aromatic = |i: usize, m: &Molecule| {
            m.atoms[i].aromatic || adj[i].iter().any(|&(_, b)| aromatic[b])
        };
        if !members.iter().any(|&i| is_aromatic(i, mol)) {
            continue;
        }
        let mut candidates = vec![false; atoms.len()];
        let mut done = Vec::new();
        let mut make_single = vec![false; mol.bonds.len()];
        for &i in &members {
            if !is_aromatic(i, mol) {
                done.push(i);
                continue;
            }
            let a = &mol.atoms[i];
            let mut sbo = a.hydrogen;
            let mut ignored = 0;
            for &(_, bi) in &adj[i] {
                let b = &mol.bonds[bi];
                if aromatic[bi] {
                    sbo += 1;
                    make_single[bi] = true;
                } else {
                    let contribution = b.valence_contribution(i).round() as i32;
                    sbo += contribution;
                    ignored += i32::from(contribution == 0);
                }
            }
            let vals = crate::graph_data::VALENCES[z[i] - 1];
            let mut charge = a.charge;
            if matches!(z[i], 3..=5 | 11..=13 | 19..=22 | 30..=32 | 37..=41 | 48..=51 | 55..=61 | 72..=73 | 80..=83 | 87..=93 | 104..=118)
            {
                charge = -charge;
            }
            if z[i] == 6 && charge > 0 {
                charge = -charge;
            }
            let mut dv = vals[0] + charge;
            for &v in vals.iter().skip(1) {
                if valences[i] <= dv || v <= 0 {
                    break;
                }
                dv = v + charge;
            }
            // SMILES radicals are assigned *after* Kekulize. Do not use our
            // derived final-state radical count while choosing candidates.
            let degree = adj[i].len() as i32 + if a.in_bracket { 0 } else { a.hydrogen } - ignored;
            if valences[i] == 5
                && sbo == 4
                && dv == 3
                && degree == 3
                && charge == 0
                && a.hydrogen == 0
                && matches!(z[i], 7 | 15 | 33)
            {
                dv = 5;
            }
            if degree < dv && (dv == sbo + 1 || (a.in_bracket && dv == sbo + 2)) {
                candidates[i] = true;
            }
        }
        for (bi, single) in make_single.into_iter().enumerate() {
            if single {
                mol.bonds[bi].order = BondOrder::Single;
            }
        }
        assign(mol, &adj, &aromatic, &atoms, &members, candidates, done)?;
    }
    for (i, a) in mol.atoms.iter_mut().enumerate() {
        if !a.aromatic {
            continue;
        }
        if !rings.atom_in_ring[i] {
            return None;
        }
        a.aromatic = false;
        if a.charge == 0
            && a.hydrogen == 1
            && a.in_bracket
            && matches!(a.symbol.as_str(), "N" | "P")
        {
            // Retain bracket/source provenance while distinguishing the
            // prepared implicit H from a literal [NH] in a nonaromatic input.
            a.prepared_h_is_implicit = true;
        }
    }
    if (0..mol.atoms.len()).any(|i| crate::electronic::total_valence(mol, i) != valences[i]) {
        return None;
    }
    Some(())
}

fn assign(
    mol: &mut Molecule,
    adj: &[Vec<(usize, usize)>],
    aromatic: &[bool],
    in_play: &[bool],
    members: &[usize],
    mut candidates: Vec<bool>,
    mut done: Vec<usize>,
) -> Option<()> {
    let mut queue = VecDeque::new();
    let mut options: Vec<Option<VecDeque<usize>>> = vec![None; mol.atoms.len()];
    let mut moves = Vec::new();
    let mut added = vec![false; mol.bonds.len()];
    let mut backtracks = 0;
    while done.len() < members.len() || !queue.is_empty() {
        let current = queue
            .pop_front()
            .or_else(|| members.iter().copied().find(|i| !done.contains(i)))?;
        done.push(current);
        let mut choices = if let Some(saved) = &options[current] {
            saved.clone()
        } else {
            let mut neighbors: Vec<_> = adj[current]
                .iter()
                .filter(|&&(j, _)| in_play[j] && !done.contains(&j))
                .copied()
                .collect();
            neighbors.sort_unstable_by_key(|&(j, _)| j);
            let mut choices = VecDeque::new();
            for (j, bi) in neighbors {
                if !queue.contains(&j) {
                    queue.push_back(j);
                }
                if candidates[current] && candidates[j] && aromatic[bi] {
                    choices.push_back(j);
                }
            }
            choices
        };
        if !candidates[current] {
            continue;
        }
        if let Some(next) = choices.pop_front() {
            let bi = adj[current].iter().find(|&&(j, _)| j == next)?.1;
            mol.bonds[bi].order = BondOrder::Double;
            mol.bonds[bi].direction = 0;
            candidates[current] = false;
            candidates[next] = false;
            added[bi] = true;
            if options[current].is_some() {
                if choices.is_empty() {
                    options[current] = None;
                    moves.pop();
                } else {
                    options[current] = Some(choices);
                }
            } else if !choices.is_empty() {
                moves.push(current);
                options[current] = Some(choices);
            }
        } else {
            let &last = moves.last()?;
            if backtracks >= 100 {
                return None;
            }
            let offset = done.iter().position(|&i| i == last)?;
            for &i in done[offset..].iter().rev() {
                queue.push_front(i);
            }
            done.truncate(offset);
            for (bi, b) in mol.bonds.iter_mut().enumerate() {
                if added[bi] && !done.contains(&b.a) && !done.contains(&b.b) {
                    added[bi] = false;
                    b.order = BondOrder::Single;
                    candidates[b.a] = true;
                    candidates[b.b] = true;
                }
            }
            backtracks += 1;
        }
    }
    Some(())
}
