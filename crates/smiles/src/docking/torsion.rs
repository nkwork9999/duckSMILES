use crate::parser::{BondOrder, Molecule};

// ── Torsion tree (flexible-docking fragment tree) ────────────────────────────
//
// A ligand is decomposed into rigid fragments connected by rotatable bonds.
// Each rotatable bond defines one torsional degree of freedom. Rotating a
// torsion rotates the "leaf-side" set of atoms about the bond axis while the
// "root-side" stays fixed.
//
// Root selection: atom 0 of the (explicit-H) molecule. The root fragment is
// held rigid by the global pose (translation + quaternion); torsions move
// everything downstream of it.
//
// A bond is rotatable iff:
//   - it is a single bond,
//   - it is NOT part of any ring (rotating a ring bond would break the ring),
//   - neither endpoint is terminal *in the heavy-atom sense*: each endpoint
//     must have at least one neighbour OTHER than the partner that is not a
//     hydrogen. This excludes X–H bonds and terminal methyl/amino rotations
//     (which are chemically degenerate for scoring purposes).
//
// Moving sets are computed by removing the candidate bond and asking which
// connected component the root lands in; the *other* component is the moving
// (leaf) side. Bonds are ordered root→leaf by the graph distance of their
// pivot atom from the root so that, when applied in order, a parent rotation
// carries its children's axes before the children rotate.

#[derive(Clone, Debug)]
pub struct RotatableBond {
    /// Root-side endpoint (the pivot; stays put under this torsion).
    pub pivot: usize,
    /// Leaf-side endpoint (defines the axis direction together with `pivot`).
    pub moving_end: usize,
    /// All atom indices that rotate when this torsion changes (includes
    /// `moving_end`, excludes `pivot`).
    pub moving: Vec<usize>,
}

#[derive(Clone, Debug, Default)]
pub struct TorsionTree {
    pub bonds: Vec<RotatableBond>,
}

impl TorsionTree {
    pub fn n_tors(&self) -> usize {
        self.bonds.len()
    }
}

// ── adjacency helpers ──────────────────────────────────────────────────────

fn adjacency(mol: &Molecule) -> Vec<Vec<usize>> {
    let n = mol.atoms.len();
    let mut adj = vec![Vec::new(); n];
    for b in &mol.bonds {
        adj[b.a].push(b.b);
        adj[b.b].push(b.a);
    }
    adj
}

/// BFS distances from `root` over the full molecule graph.
fn bfs_distances(adj: &[Vec<usize>], root: usize) -> Vec<i32> {
    let mut dist = vec![-1_i32; adj.len()];
    let mut queue = std::collections::VecDeque::new();
    dist[root] = 0;
    queue.push_back(root);
    while let Some(u) = queue.pop_front() {
        for &v in &adj[u] {
            if dist[v] < 0 {
                dist[v] = dist[u] + 1;
                queue.push_back(v);
            }
        }
    }
    dist
}

/// Atoms reachable from `start` without traversing the edge (avoid_a, avoid_b).
fn component_excluding_bond(
    adj: &[Vec<usize>],
    start: usize,
    avoid_a: usize,
    avoid_b: usize,
) -> Vec<usize> {
    let mut seen = vec![false; adj.len()];
    let mut stack = vec![start];
    seen[start] = true;
    let mut out = Vec::new();
    while let Some(u) = stack.pop() {
        out.push(u);
        for &v in &adj[u] {
            // skip the cut edge in both directions
            if (u == avoid_a && v == avoid_b) || (u == avoid_b && v == avoid_a) {
                continue;
            }
            if !seen[v] {
                seen[v] = true;
                stack.push(v);
            }
        }
    }
    out
}

/// True if `atom` has at least one neighbour other than `partner` that is not H.
fn has_heavy_neighbour_besides(mol: &Molecule, atom: usize, partner: usize, adj: &[Vec<usize>]) -> bool {
    adj[atom]
        .iter()
        .any(|&nbr| nbr != partner && mol.atoms[nbr].symbol != "H")
}

// ── builder ────────────────────────────────────────────────────────────────

/// Build the torsion tree for an explicit-H molecule. `root` is the atom held
/// fixed by the global rigid pose (callers pass 0).
pub fn build_torsion_tree(mol: &Molecule, root: usize) -> TorsionTree {
    let adj = adjacency(mol);
    let ring = mol.ring_info();
    let dist = bfs_distances(&adj, root);

    let mut bonds: Vec<RotatableBond> = Vec::new();

    for (bi, b) in mol.bonds.iter().enumerate() {
        if b.order != BondOrder::Single {
            continue;
        }
        if ring.bond_in_ring.get(bi).copied().unwrap_or(false) {
            continue;
        }
        // both endpoints must be non-terminal in the heavy sense
        if !has_heavy_neighbour_besides(mol, b.a, b.b, &adj)
            || !has_heavy_neighbour_besides(mol, b.b, b.a, &adj)
        {
            continue;
        }
        // Skip bonds the root can't reach (disconnected fragments / salts).
        if dist[b.a] < 0 || dist[b.b] < 0 {
            continue;
        }

        // Pivot = the endpoint closer to the root; moving_end = the farther one.
        let (pivot, moving_end) = if dist[b.a] <= dist[b.b] {
            (b.a, b.b)
        } else {
            (b.b, b.a)
        };

        // Moving set = component containing moving_end after cutting the bond.
        let moving = component_excluding_bond(&adj, moving_end, pivot, moving_end);

        // Guard: if the moving set somehow contains the root, this bond is a
        // ring-closure the ring perception missed — skip it.
        if moving.contains(&root) {
            continue;
        }

        bonds.push(RotatableBond { pivot, moving_end, moving });
    }

    // Order root→leaf by pivot distance so parent rotations precede children.
    bonds.sort_by_key(|rb| dist[rb.pivot]);

    TorsionTree { bonds }
}

// ── Rodrigues rotation of a point about an axis through a base point ─────────

/// Rotate `p` about the line through `base` with unit direction `axis` by
/// `angle` radians. Returns the rotated point.
pub fn rotate_about_axis(p: [f64; 3], base: [f64; 3], axis: [f64; 3], angle: f64) -> [f64; 3] {
    let (s, c) = angle.sin_cos();
    let one_c = 1.0 - c;
    // translate so base is origin
    let v = [p[0] - base[0], p[1] - base[1], p[2] - base[2]];
    let (kx, ky, kz) = (axis[0], axis[1], axis[2]);
    // k × v
    let cross = [
        ky * v[2] - kz * v[1],
        kz * v[0] - kx * v[2],
        kx * v[1] - ky * v[0],
    ];
    // k · v
    let dot = kx * v[0] + ky * v[1] + kz * v[2];
    // Rodrigues: v*cos + (k×v)*sin + k*(k·v)*(1-cos)
    let r = [
        v[0] * c + cross[0] * s + kx * dot * one_c,
        v[1] * c + cross[1] * s + ky * dot * one_c,
        v[2] * c + cross[2] * s + kz * dot * one_c,
    ];
    [r[0] + base[0], r[1] + base[1], r[2] + base[2]]
}

/// Apply a set of torsion angles to a reference conformation, returning new
/// internal coordinates. `coords` is consumed as the starting conformation;
/// each torsion `i` rotates `tree.bonds[i].moving` about its (current) axis by
/// `angles[i]` radians. Applied root→leaf (tree is pre-ordered).
pub fn apply_torsions(
    ref_coords: &[[f64; 3]],
    tree: &TorsionTree,
    angles: &[f64],
) -> Vec<[f64; 3]> {
    let mut coords = ref_coords.to_vec();
    for (i, rb) in tree.bonds.iter().enumerate() {
        let angle = angles.get(i).copied().unwrap_or(0.0);
        if angle == 0.0 {
            continue;
        }
        let base = coords[rb.pivot];
        let end = coords[rb.moving_end];
        let mut axis = [end[0] - base[0], end[1] - base[1], end[2] - base[2]];
        let norm = (axis[0] * axis[0] + axis[1] * axis[1] + axis[2] * axis[2]).sqrt();
        if norm < 1e-9 {
            continue;
        }
        axis = [axis[0] / norm, axis[1] / norm, axis[2] / norm];
        for &atom in &rb.moving {
            coords[atom] = rotate_about_axis(coords[atom], base, axis, angle);
        }
    }
    coords
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn mol(smi: &str) -> Molecule {
        parse(smi).unwrap().with_explicit_hydrogens()
    }

    #[test]
    fn ethane_has_no_rotatable_bond() {
        // C-C with only H's on each side → both ends terminal in heavy sense.
        let m = mol("CC");
        let tree = build_torsion_tree(&m, 0);
        assert_eq!(tree.n_tors(), 0, "ethane C-C should not be rotatable");
    }

    #[test]
    fn butane_has_one_rotatable_bond() {
        // CCCC: only the central C-C is rotatable (terminal C-C's are methyls).
        let m = mol("CCCC");
        let tree = build_torsion_tree(&m, 0);
        assert_eq!(tree.n_tors(), 1, "butane should have exactly 1 rotatable bond");
    }

    #[test]
    fn benzene_ring_bonds_not_rotatable() {
        let m = mol("c1ccccc1");
        let tree = build_torsion_tree(&m, 0);
        assert_eq!(tree.n_tors(), 0, "aromatic ring bonds must not be rotatable");
    }

    #[test]
    fn biphenyl_central_bond_rotatable() {
        // Two phenyls joined by a single bond → 1 rotatable bond.
        let m = mol("c1ccccc1-c1ccccc1");
        let tree = build_torsion_tree(&m, 0);
        assert_eq!(tree.n_tors(), 1, "biphenyl central bond should be rotatable");
    }

    #[test]
    fn moving_set_excludes_root() {
        let m = mol("CCCC");
        let tree = build_torsion_tree(&m, 0);
        for rb in &tree.bonds {
            assert!(!rb.moving.contains(&0), "moving set must not contain root");
            assert!(rb.moving.contains(&rb.moving_end));
        }
    }

    #[test]
    fn rotate_about_z_axis_90deg() {
        // (1,0,0) rotated 90° about z through origin → (0,1,0)
        let p = rotate_about_axis([1.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 1.0],
            std::f64::consts::FRAC_PI_2);
        assert!((p[0] - 0.0).abs() < 1e-9, "x={}", p[0]);
        assert!((p[1] - 1.0).abs() < 1e-9, "y={}", p[1]);
        assert!((p[2] - 0.0).abs() < 1e-9, "z={}", p[2]);
    }

    #[test]
    fn apply_torsions_preserves_bond_lengths() {
        // Rotating a torsion must not change any covalent bond length.
        let m = mol("CCCC");
        let tree = build_torsion_tree(&m, 0);
        // Make a simple stretched-out conformation.
        let n = m.atoms.len();
        let coords: Vec<[f64; 3]> = (0..n)
            .map(|i| [i as f64 * 1.5, (i % 2) as f64 * 0.5, 0.0])
            .collect();
        let angles = vec![1.0_f64; tree.n_tors()];
        let out = apply_torsions(&coords, &tree, &angles);

        for b in &m.bonds {
            let before = {
                let d = [
                    coords[b.a][0] - coords[b.b][0],
                    coords[b.a][1] - coords[b.b][1],
                    coords[b.a][2] - coords[b.b][2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            };
            let after = {
                let d = [
                    out[b.a][0] - out[b.b][0],
                    out[b.a][1] - out[b.b][1],
                    out[b.a][2] - out[b.b][2],
                ];
                (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt()
            };
            assert!((before - after).abs() < 1e-9,
                "bond {}-{} length changed: {before} → {after}", b.a, b.b);
        }
    }

    #[test]
    fn zero_torsion_is_identity() {
        let m = mol("CCCC");
        let tree = build_torsion_tree(&m, 0);
        let n = m.atoms.len();
        let coords: Vec<[f64; 3]> = (0..n).map(|i| [i as f64, 0.0, 0.0]).collect();
        let angles = vec![0.0_f64; tree.n_tors()];
        let out = apply_torsions(&coords, &tree, &angles);
        for i in 0..n {
            for k in 0..3 {
                assert!((coords[i][k] - out[i][k]).abs() < 1e-12);
            }
        }
    }
}
