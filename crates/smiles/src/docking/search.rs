use crate::docking::affinity::AffinityMap;
use crate::docking::atomtype::VinaType;
use crate::docking::pose::{apply_pose_flex, centre_coords, Pose, Quat};
use crate::docking::score::{repulsion, vdw, W_REPULSION};
use crate::docking::torsion::{apply_torsions, build_torsion_tree, TorsionTree};

// ── Docking Search ────────────────────────────────────────────────────────────
// Monte Carlo search with local gradient descent over the full pose:
//   translation (3) + orientation quaternion (4) + torsions (n_rot).
//
// Algorithm:
//  1. Generate `n_runs` random initial poses within the grid box.
//  2. For each pose, run local gradient descent (numerical gradient) over all
//     degrees of freedom (translation, rotation, AND torsions).
//  3. Score = inter-molecular (affinity map) + intra-molecular ligand strain.
//  4. Return top-k poses sorted by score.
//
// Phase 5 adds torsional flexibility: each rotatable bond is an extra DoF and
// `apply_pose_flex` bends the conformation before placing it.
// Phase 6 adds intra-molecular scoring: a repulsion-only strain term over
// 1-5-and-beyond atom pairs penalises self-clashing (folded) conformations so
// the optimiser cannot cheat by collapsing the ligand into the protein.

// ── LCG PRNG ─────────────────────────────────────────────────────────────────

struct Lcg { state: u64 }
impl Lcg {
    fn new(seed: u64) -> Self { Lcg { state: seed.wrapping_add(1) } }
    fn next_u64(&mut self) -> u64 {
        self.state = self.state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        self.state
    }
    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    /// Sample from [-1, 1]
    fn next_pm1(&mut self) -> f64 { self.next_f64() * 2.0 - 1.0 }
}

// ── Ligand bundle (reference conformation + flexibility metadata) ────────────

/// Everything the search needs about a ligand: its centred reference
/// conformation, per-atom Vina types, the torsion tree, the intra-molecular
/// non-bonded pair list, and the rotatable-bond count for the flexibility
/// penalty.
pub struct Ligand {
    pub ref_coords: Vec<[f64; 3]>,
    pub types: Vec<VinaType>,
    pub tree: TorsionTree,
    /// Heavy-atom pairs separated by >= 4 bonds (1-5 and beyond), scored for
    /// intramolecular strain (repulsion only).
    pub intra_pairs: Vec<(usize, usize)>,
    pub n_rot: u32,
}

impl Ligand {
    /// Build a rigid ligand (no torsions, no intra pairs) — used by the simple
    /// `dock()` entry point and its tests.
    fn rigid(ref_coords: Vec<[f64; 3]>, types: Vec<VinaType>, n_rot: u32) -> Self {
        Ligand {
            ref_coords,
            types,
            tree: TorsionTree::default(),
            intra_pairs: Vec::new(),
            n_rot,
        }
    }
}

// ── Intra-molecular strain (repulsion only) ─────────────────────────────────

/// Repulsion-only internal energy over the precomputed 1-5+ pair list. Uses the
/// internal (torsioned) coordinates; translation/rotation invariant.
fn intra_energy(coords: &[[f64; 3]], types: &[VinaType], pairs: &[(usize, usize)]) -> f64 {
    let mut e = 0.0;
    for &(i, j) in pairs {
        let dx = coords[i][0] - coords[j][0];
        let dy = coords[i][1] - coords[j][1];
        let dz = coords[i][2] - coords[j][2];
        let r = (dx * dx + dy * dy + dz * dz).sqrt();
        let d = r - vdw(types[i]) - vdw(types[j]);
        e += W_REPULSION * repulsion(d);
    }
    e
}

// ── pose <-> vector ─────────────────────────────────────────────────────────

fn random_pose(rng: &mut Lcg, center: [f64; 3], box_size: [f64; 3], n_tors: usize) -> Pose {
    let tx = center[0] + rng.next_pm1() * box_size[0];
    let ty = center[1] + rng.next_pm1() * box_size[1];
    let tz = center[2] + rng.next_pm1() * box_size[2];

    let (a, b, c, d) = (rng.next_pm1(), rng.next_pm1(), rng.next_pm1(), rng.next_pm1());
    let q = Quat { w: a, x: b, y: c, z: d }.normalise();

    let torsions = (0..n_tors)
        .map(|_| rng.next_pm1() * std::f64::consts::PI)
        .collect();

    Pose { translation: [tx, ty, tz], orientation: q, torsions }
}

fn pose_to_vec(pose: &Pose) -> Vec<f64> {
    let mut v = vec![
        pose.translation[0], pose.translation[1], pose.translation[2],
        pose.orientation.w, pose.orientation.x, pose.orientation.y, pose.orientation.z,
    ];
    v.extend_from_slice(&pose.torsions);
    v
}

fn vec_to_pose(v: &[f64]) -> Pose {
    let q = Quat { w: v[3], x: v[4], y: v[5], z: v[6] }.normalise();
    Pose {
        translation: [v[0], v[1], v[2]],
        orientation: q,
        torsions: v[7..].to_vec(),
    }
}

// ── pose scoring (inter via map + intra strain) ─────────────────────────────

const GD_ITERS: usize = 200;
const GD_H: f64 = 0.01;
const GD_TOL: f64 = 1e-4;
const OUT_OF_GRID: f64 = 1e9;

fn score_pose(map: &AffinityMap, lig: &Ligand, pose: &Pose) -> Option<f64> {
    // Inter-molecular: place the (bent) ligand and look up the maps.
    let coords = apply_pose_flex(&lig.ref_coords, pose, &lig.tree);
    let inter = map.score_ligand(&coords, &lig.types, lig.n_rot)?;

    // Intra-molecular: compute strain from the internal (bent) conformation.
    let intra = if lig.intra_pairs.is_empty() {
        0.0
    } else {
        let internal = apply_torsions(&lig.ref_coords, &lig.tree, &pose.torsions);
        intra_energy(&internal, &lig.types, &lig.intra_pairs)
    };

    Some(inter + intra)
}

fn local_optimise(map: &AffinityMap, lig: &Ligand, pose: &Pose, margin: f64) -> (Pose, f64) {
    let mut x = pose_to_vec(pose);
    let n = x.len();
    let mut step = 0.1_f64;

    let t_lo = [
        map.origin[0] + margin,
        map.origin[1] + margin,
        map.origin[2] + margin,
    ];
    let t_hi = [
        map.origin[0] + (map.n[0].saturating_sub(2)) as f64 * map.spacing - margin,
        map.origin[1] + (map.n[1].saturating_sub(2)) as f64 * map.spacing - margin,
        map.origin[2] + (map.n[2].saturating_sub(2)) as f64 * map.spacing - margin,
    ];
    let clamp_t = |v: &mut Vec<f64>| {
        for k in 0..3 {
            v[k] = v[k].clamp(t_lo[k], t_hi[k]);
        }
    };
    clamp_t(&mut x);

    let eval = |x: &[f64]| -> f64 {
        let p = vec_to_pose(x);
        score_pose(map, lig, &p).unwrap_or(OUT_OF_GRID)
    };

    let mut e = eval(&x);

    for _ in 0..GD_ITERS {
        let mut grad = vec![0.0_f64; n];
        for i in 0..n {
            let orig = x[i];
            x[i] = orig + GD_H;
            let ep = eval(&x);
            x[i] = orig - GD_H;
            let em = eval(&x);
            x[i] = orig;
            grad[i] = (ep - em) / (2.0 * GD_H);
        }

        let gnorm = grad.iter().map(|g| g * g).sum::<f64>().sqrt();
        if gnorm < GD_TOL { break; }

        let slope: f64 = grad.iter().map(|g| g * g).sum::<f64>();
        let mut alpha = step;
        let mut accepted = false;
        for _ in 0..10 {
            let xt: Vec<f64> = x.iter().zip(&grad).map(|(xi, gi)| xi - alpha * gi).collect();
            let et = eval(&xt);
            if et < e - 1e-4 * alpha * slope {
                x = xt;
                e = et;
                accepted = true;
                break;
            }
            alpha *= 0.5;
        }
        if accepted {
            step = (step * 1.2).min(0.5);
            clamp_t(&mut x);
        }
    }

    (vec_to_pose(&x), e)
}

// ── DockResult ────────────────────────────────────────────────────────────────

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct DockResult {
    pub score: f64,            // kcal/mol
    pub pose: Pose,
    pub coords: Vec<[f64; 3]>, // final ligand atom coordinates
}

// ── core docking loop ───────────────────────────────────────────────────────

fn dock_ligand(lig: &Ligand, map: &AffinityMap, n_runs: usize, seed: u64) -> Vec<DockResult> {
    let mut rng = Lcg::new(seed);
    let center = map.center;

    // Ligand radius (max atom distance from its centroid). The conformation can
    // rotate freely, so the worst-case footprint is a sphere of this radius;
    // keep the centroid far enough from the grid edge that no atom escapes.
    let radius = lig
        .ref_coords
        .iter()
        .map(|c| (c[0] * c[0] + c[1] * c[1] + c[2] * c[2]).sqrt())
        .fold(0.0_f64, f64::max);
    let margin = radius + 1.0;

    // Initial-placement box: keep the whole ligand inside the grid from the
    // start so the first score is finite (non-zero gradient to optimise).
    let half = |k: usize| (map.n[k] as f64 * map.spacing) / 2.0;
    let box_size = [
        (half(0) - margin).max(0.0),
        (half(1) - margin).max(0.0),
        (half(2) - margin).max(0.0),
    ];
    let n_tors = lig.tree.n_tors();

    let mut results: Vec<DockResult> = Vec::new();
    for _ in 0..n_runs {
        let init_pose = random_pose(&mut rng, center, box_size, n_tors);
        let (opt_pose, score) = local_optimise(map, lig, &init_pose, margin);
        if score < OUT_OF_GRID * 0.1 {
            let coords = apply_pose_flex(&lig.ref_coords, &opt_pose, &lig.tree);
            results.push(DockResult { score, pose: opt_pose, coords });
        }
    }

    results.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
    results
}

// ── dock (rigid, direct-coords entry) ───────────────────────────────────────
// Backwards-compatible rigid docking: caller supplies coords + types directly.

pub fn dock(
    ref_coords: &[[f64; 3]],
    lig_types: &[VinaType],
    map: &AffinityMap,
    n_rot: u32,
    n_runs: usize,
    seed: u64,
) -> Vec<DockResult> {
    let lig = Ligand::rigid(ref_coords.to_vec(), lig_types.to_vec(), n_rot);
    dock_ligand(&lig, map, n_runs, seed)
}

// ── dock_from_smiles (flexible) ─────────────────────────────────────────────
// SMILES → 3D conformer → typed → torsion tree → flexible docking.

pub fn dock_from_smiles(
    smiles: &str,
    map: &AffinityMap,
    n_runs: usize,
    seed: u64,
) -> Option<Vec<DockResult>> {
    let lig = build_ligand(smiles, seed)?;
    Some(dock_ligand(&lig, map, n_runs, seed))
}

/// Construct a flexible `Ligand` from SMILES: conformer + types + torsion tree
/// + intramolecular pair list.
fn build_ligand(smiles: &str, seed: u64) -> Option<Ligand> {
    use crate::conformer::generate_conformer;
    use crate::docking::atomtype::assign_types;
    use crate::parser::parse;

    let mol_bare = parse(smiles)?;
    let mol = mol_bare.with_explicit_hydrogens();
    let coords_flat = generate_conformer(smiles, seed)?;

    let n = mol.atoms.len();
    let mut ref_coords: Vec<[f64; 3]> = (0..n)
        .map(|i| [coords_flat[3 * i], coords_flat[3 * i + 1], coords_flat[3 * i + 2]])
        .collect();
    centre_coords(&mut ref_coords);

    let types = assign_types(&mol);
    let tree = build_torsion_tree(&mol, 0);
    let n_rot = tree.n_tors() as u32;
    let intra_pairs = intra_pair_list(&mol);

    Some(Ligand { ref_coords, types, tree, intra_pairs, n_rot })
}

/// Heavy-atom pairs separated by >= 4 bonds (graph distance), used for the
/// intramolecular strain term. Hydrogens are excluded for stability.
fn intra_pair_list(mol: &crate::parser::Molecule) -> Vec<(usize, usize)> {
    let n = mol.atoms.len();
    // adjacency
    let mut adj = vec![Vec::new(); n];
    for b in &mol.bonds {
        adj[b.a].push(b.b);
        adj[b.b].push(b.a);
    }
    let is_heavy = |i: usize| mol.atoms[i].symbol != "H";

    let mut pairs = Vec::new();
    for src in 0..n {
        if !is_heavy(src) { continue; }
        // BFS distances from src
        let mut dist = vec![-1_i32; n];
        let mut q = std::collections::VecDeque::new();
        dist[src] = 0;
        q.push_back(src);
        while let Some(u) = q.pop_front() {
            for &v in &adj[u] {
                if dist[v] < 0 {
                    dist[v] = dist[u] + 1;
                    q.push_back(v);
                }
            }
        }
        for dst in (src + 1)..n {
            if !is_heavy(dst) { continue; }
            if dist[dst] >= 4 {
                pairs.push((src, dst));
            }
        }
    }
    pairs
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::docking::affinity::AffinityMap;

    fn make_map() -> AffinityMap {
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::OA];
        AffinityMap::build(&prot_coords, &prot_types, [0.0, 0.0, 0.0], [8.0, 8.0, 8.0], 0.5)
    }

    #[test]
    fn dock_returns_results() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64]];
        let lig_types = vec![VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 3, 42);
        assert!(!results.is_empty(), "dock should return at least one result");
    }

    #[test]
    fn dock_results_sorted_by_score() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64]];
        let lig_types = vec![VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 5, 1);
        for w in results.windows(2) {
            assert!(w[0].score <= w[1].score, "results not sorted: {} > {}", w[0].score, w[1].score);
        }
    }

    #[test]
    fn dock_coords_length_matches_ligand() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0]];
        let lig_types = vec![VinaType::C, VinaType::C];
        let results = dock(&ref_coords, &lig_types, &map, 0, 3, 0);
        for r in &results {
            assert_eq!(r.coords.len(), 2);
        }
    }

    #[test]
    fn dock_best_hd_near_oa() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64]];
        let lig_types = vec![VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 20, 99);
        assert!(!results.is_empty(), "should find at least one pose");
        if let Some(best) = results.first() {
            let [x, y, z] = best.coords[0];
            let r = (x * x + y * y + z * z).sqrt();
            assert!(r < 8.0, "best pose outside grid: r={r:.2}");
            assert!(best.score < 0.0, "HD-OA should give negative score: {}", best.score);
        }
    }

    #[test]
    fn dock_simple_two_atom_ligand() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0]];
        let lig_types = vec![VinaType::C, VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 10, 0);
        assert!(!results.is_empty(), "dock should find at least one pose");
        assert!(!results[0].score.is_nan());
    }

    #[test]
    fn dock_from_smiles_water() {
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::NS];
        let map = AffinityMap::build(&prot_coords, &prot_types, [0.0, 0.0, 0.0], [8.0, 8.0, 8.0], 0.5);
        let results = dock_from_smiles("O", &map, 10, 7);
        assert!(results.is_some(), "dock_from_smiles should succeed for water");
        assert!(!results.unwrap().is_empty(), "should find at least one valid pose");
    }

    // ── Phase 5/6 specific tests ────────────────────────────────────────────

    #[test]
    fn build_ligand_butane_has_torsion() {
        // butane: 1 rotatable bond → ligand carries 1 torsion DoF.
        let lig = build_ligand("CCCC", 1).expect("butane ligand builds");
        assert_eq!(lig.n_rot, 1, "butane should expose 1 rotatable bond");
        assert_eq!(lig.tree.n_tors(), 1);
    }

    #[test]
    fn build_ligand_intra_pairs_nonempty_for_pentane() {
        // pentane heavy atoms C0..C4: C0–C4 are 4 bonds apart → at least 1 pair.
        let lig = build_ligand("CCCCC", 2).expect("pentane ligand builds");
        assert!(!lig.intra_pairs.is_empty(), "pentane should have >=1 intra (1-5) pair");
        // All intra pairs must be heavy-heavy.
        for &(i, j) in &lig.intra_pairs {
            assert_ne!(lig.types[i], VinaType::H);
            assert_ne!(lig.types[j], VinaType::H);
        }
    }

    #[test]
    fn flexible_dock_runs_and_scores_finite() {
        // A small flexible molecule docked against an OA atom — must complete
        // and yield a finite best score with all torsions applied.
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::OA];
        let map = AffinityMap::build(&prot_coords, &prot_types, [0.0, 0.0, 0.0], [10.0, 10.0, 10.0], 0.5);
        let results = dock_from_smiles("CCO", &map, 10, 3).expect("ethanol docks");
        assert!(!results.is_empty(), "flexible dock should find a pose");
        assert!(results[0].score.is_finite());
    }

    #[test]
    fn flexible_dock_with_torsion_optimises_dofs() {
        // 1,2-dichloroethane: C-C is rotatable → docking optimises a torsion DoF
        // and the returned pose carries one torsion angle.
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::OA];
        let map = AffinityMap::build(&prot_coords, &prot_types, [0.0, 0.0, 0.0], [12.0, 12.0, 12.0], 0.5);
        let lig = build_ligand("ClCCCl", 5).expect("builds");
        assert_eq!(lig.n_rot, 1, "ClCCCl has one rotatable C-C bond");
        let results = dock_ligand(&lig, &map, 8, 5);
        assert!(!results.is_empty(), "flexible dock should produce a pose");
        assert_eq!(results[0].pose.torsions.len(), 1, "pose should carry one torsion DoF");
        assert!(results[0].score.is_finite());
    }

    #[test]
    fn intra_energy_penalises_overlap() {
        // Two C atoms on top of each other → large positive repulsion.
        let coords = vec![[0.0, 0.0, 0.0_f64], [0.1, 0.0, 0.0]];
        let types = vec![VinaType::C, VinaType::C];
        let pairs = vec![(0usize, 1usize)];
        let e = intra_energy(&coords, &types, &pairs);
        assert!(e > 0.0, "overlapping atoms must incur positive strain: {e}");
        // Well-separated → zero.
        let coords2 = vec![[0.0, 0.0, 0.0_f64], [6.0, 0.0, 0.0]];
        let e2 = intra_energy(&coords2, &types, &pairs);
        assert_eq!(e2, 0.0, "separated atoms incur no repulsion");
    }
}
