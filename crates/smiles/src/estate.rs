//! Electrotopological state, Hall/Mohney/Kier, JCICS 31, 76–81 (1991).
//! Pinned to RDKit 2026.03.6 EState.py's definition, including explicit H nodes.
//! Upstream: Copyright (C) 2002–2006 Greg Landrum and Rational Discovery LLC.
//! See RDKIT-DATA-LICENSE.txt (BSD-3-Clause).
use crate::parser::Molecule;

pub(crate) fn indices(mol: &Molecule) -> Option<Vec<f64>> {
    let n = mol.atoms.len();
    if n == 0 {
        return None;
    }
    let adj = mol.adjacency();
    let mut intrinsic = Vec::with_capacity(n);
    for (i, atom) in mol.atoms.iter().enumerate() {
        let z = crate::periodic_data::ELEMENTS
            .iter()
            .position(|e| e.0 == atom.symbol)?
            + 1;
        let period = match z {
            1..=2 => 1,
            3..=10 => 2,
            11..=18 => 3,
            19..=36 => 4,
            37..=54 => 5,
            55..=86 => 6,
            _ => 7,
        };
        // RDKit GetTotalNumHs() excludes separate hydrogen neighbours here.
        // The formal charge is not subtracted from this definition's dv.
        let dv = crate::descriptor_catalog::OUTER_ELECTRONS[z - 1] - atom.hydrogen;
        intrinsic.push(if adj[i].is_empty() {
            0.0
        } else {
            (4.0 / f64::from(period * period) * f64::from(dv) + 1.0) / adj[i].len() as f64
        });
    }
    let mut accum = vec![0.0; n];
    let mut distances = vec![usize::MAX; n];
    let mut queue = Vec::with_capacity(n);
    for i in 0..n {
        distances.fill(usize::MAX);
        distances[i] = 0;
        queue.clear();
        queue.push(i);
        let mut front = 0;
        while front < queue.len() {
            let a = queue[front];
            front += 1;
            for &(b, _) in &adj[a] {
                if distances[b] == usize::MAX {
                    distances[b] = distances[a] + 1;
                    queue.push(b);
                }
            }
        }
        for j in i + 1..n {
            if distances[j] != usize::MAX {
                let separation = (distances[j] + 1) as f64;
                let delta = (intrinsic[i] - intrinsic[j]) / (separation * separation);
                accum[i] += delta;
                accum[j] -= delta;
            }
        }
    }
    Some(intrinsic.iter().zip(accum).map(|(a, b)| a + b).collect())
}
