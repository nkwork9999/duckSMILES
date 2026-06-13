use crate::docking::affinity::AffinityMap;
use crate::docking::atomtype::VinaType;
use crate::docking::pose::{apply_pose, centre_coords, Pose, Quat};

// ── Docking Search ────────────────────────────────────────────────────────────
// Monte Carlo search with local gradient descent.
//
// Algorithm:
//  1. Generate `n_runs` random initial poses within the grid box.
//  2. For each pose, run local gradient descent (numerical gradient).
//  3. Collect best pose per run.
//  4. Return top-k poses sorted by score.
//
// Torsional flexibility: Phase 4 searches translation+rotation only (rigid
// docking). Torsion angles are set to 0 for all rotatable bonds.

// ── LCG PRNG (reuse from embed.rs pattern) ───────────────────────────────────

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

// ── Random pose generator ─────────────────────────────────────────────────────

fn random_pose(rng: &mut Lcg, center: [f64; 3], box_size: [f64; 3], n_tors: usize) -> Pose {
    // Random translation within box
    let tx = center[0] + rng.next_pm1() * box_size[0];
    let ty = center[1] + rng.next_pm1() * box_size[1];
    let tz = center[2] + rng.next_pm1() * box_size[2];

    // Random unit quaternion via Gaussian method (Shoemake 1992 approximation)
    let (a, b, c, d) = (rng.next_pm1(), rng.next_pm1(), rng.next_pm1(), rng.next_pm1());
    let q = Quat { w: a, x: b, y: c, z: d }.normalise();

    // Random torsions in [-π, π]
    let torsions = (0..n_tors)
        .map(|_| rng.next_pm1() * std::f64::consts::PI)
        .collect();

    Pose { translation: [tx, ty, tz], orientation: q, torsions }
}

// ── Numerical gradient descent ────────────────────────────────────────────────
// Perturbs each DoF with step h and evaluates the score.
// Uses a simple back-tracking line search for step size.

const GD_ITERS: usize = 200;
const GD_H: f64 = 0.01;     // finite difference step (Å or rad)
const GD_TOL: f64 = 1e-4;   // gradient norm tolerance

fn score_pose(
    map: &AffinityMap,
    ref_coords: &[[f64; 3]],
    lig_types: &[VinaType],
    pose: &Pose,
    n_rot: u32,
) -> Option<f64> {
    let coords = apply_pose(ref_coords, pose);
    map.score_ligand(&coords, lig_types, n_rot)
}

fn pose_to_vec(pose: &Pose) -> Vec<f64> {
    let mut v = vec![
        pose.translation[0],
        pose.translation[1],
        pose.translation[2],
        pose.orientation.w,
        pose.orientation.x,
        pose.orientation.y,
        pose.orientation.z,
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

fn local_optimise(
    map: &AffinityMap,
    ref_coords: &[[f64; 3]],
    lig_types: &[VinaType],
    pose: &Pose,
    n_rot: u32,
) -> (Pose, f64) {
    let mut x = pose_to_vec(pose);
    let n = x.len();
    let mut step = 0.1_f64;

    // Bounds for translation (clamp to keep pose center inside grid, 2 Å margin)
    let margin = 2.0;
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
        score_pose(map, ref_coords, lig_types, &p, n_rot).unwrap_or(1e9)
    };

    let mut e = eval(&x);

    for _ in 0..GD_ITERS {
        // Numerical gradient
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

        // Armijo line search
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
    pub score: f64,           // kcal/mol
    pub pose: Pose,
    pub coords: Vec<[f64; 3]>, // final ligand atom coordinates
}

// ── dock ──────────────────────────────────────────────────────────────────────
// Main docking entry point.
//
// ref_coords: Phase 1 conformer with CoM at origin
// lig_types:  Vina atom types for each ligand atom
// map:        pre-computed affinity map (built from protein in Phase 3)
// n_rot:      number of rotatable bonds in ligand
// n_runs:     number of independent MC runs (typically 10-50)
// seed:       random seed

pub fn dock(
    ref_coords: &[[f64; 3]],
    lig_types: &[VinaType],
    map: &AffinityMap,
    n_rot: u32,
    n_runs: usize,
    seed: u64,
) -> Vec<DockResult> {
    let mut rng = Lcg::new(seed);
    let center = map.center;
    // Use grid half-extents as box size for random placement
    let box_size = [
        map.n[0] as f64 * map.spacing / 2.0 * 0.8,
        map.n[1] as f64 * map.spacing / 2.0 * 0.8,
        map.n[2] as f64 * map.spacing / 2.0 * 0.8,
    ];
    let n_tors = 0; // Phase 4: rigid docking only

    let mut results: Vec<DockResult> = Vec::new();

    for _ in 0..n_runs {
        let init_pose = random_pose(&mut rng, center, box_size, n_tors);
        let (opt_pose, score) = local_optimise(map, ref_coords, lig_types, &init_pose, n_rot);
        if score < 1e8 {
            let coords = apply_pose(ref_coords, &opt_pose);
            results.push(DockResult { score, pose: opt_pose, coords });
        }
    }

    // Sort by score (ascending = most negative = best)
    results.sort_by(|a, b| a.score.partial_cmp(&b.score).unwrap());
    results
}

// ── dock_from_smiles ──────────────────────────────────────────────────────────
// Convenience wrapper: SMILES string → 3D conformer → typed → docked.

pub fn dock_from_smiles(
    smiles: &str,
    map: &AffinityMap,
    n_runs: usize,
    seed: u64,
) -> Option<Vec<DockResult>> {
    use crate::conformer::generate_conformer;
    use crate::docking::atomtype::assign_types;
    use crate::parser::parse;

    let mol_bare = parse(smiles)?;
    let mol = mol_bare.with_explicit_hydrogens();
    let mut coords_flat = generate_conformer(smiles, seed)?;

    // Reshape flat Vec<f64> into Vec<[f64;3]>
    let n = mol.atoms.len();
    let mut ref_coords: Vec<[f64; 3]> = (0..n)
        .map(|i| [coords_flat[3 * i], coords_flat[3 * i + 1], coords_flat[3 * i + 2]])
        .collect();
    centre_coords(&mut ref_coords);
    coords_flat.clear();

    let lig_types = assign_types(&mol);

    // Count rotatable bonds (single bonds not in ring, not C=X)
    let n_rot = mol.bonds.iter().filter(|b| {
        use crate::parser::BondOrder;
        if b.order != BondOrder::Single { return false; }
        if mol.atoms[b.a].symbol == "H" || mol.atoms[b.b].symbol == "H" { return false; }
        true
    }).count() as u32;

    Some(dock(&ref_coords, &lig_types, map, n_rot, n_runs, seed))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::docking::affinity::AffinityMap;

    fn make_map() -> AffinityMap {
        // Single OA protein atom, grid around it
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::OA];
        AffinityMap::build(
            &prot_coords,
            &prot_types,
            [0.0, 0.0, 0.0],
            [8.0, 8.0, 8.0],
            0.5,
        )
    }

    #[test]
    fn dock_returns_results() {
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64]]; // single atom ligand
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
        // HD atom should be attracted to OA protein atom.
        // Best pose should end up inside the grid and have a negative score.
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64]];
        let lig_types = vec![VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 20, 99);
        assert!(!results.is_empty(), "should find at least one pose");
        if let Some(best) = results.first() {
            let [x, y, z] = best.coords[0];
            let r = (x * x + y * y + z * z).sqrt();
            // Best pose should be within the grid (< 8 Å) and attracted
            assert!(r < 8.0, "best pose outside grid: r={r:.2}");
            assert!(best.score < 0.0, "HD-OA should give negative score: {}", best.score);
        }
    }

    #[test]
    fn dock_simple_two_atom_ligand() {
        // Directly dock a 2-atom ligand (C + HD) to OA protein — no conformer needed
        let map = make_map();
        let ref_coords = vec![[0.0, 0.0, 0.0_f64], [1.0, 0.0, 0.0]];
        let lig_types = vec![VinaType::C, VinaType::HD];
        let results = dock(&ref_coords, &lig_types, &map, 0, 10, 0);
        assert!(!results.is_empty(), "dock should find at least one pose");
        assert!(!results[0].score.is_nan());
    }

    #[test]
    fn dock_from_smiles_water() {
        // Water is tiny (3 atoms, max ~1 Å from CoM) — definitely fits in size=8 grid
        let prot_coords = vec![[0.0_f64, 0.0, 0.0]];
        let prot_types = vec![VinaType::NS];
        let map = AffinityMap::build(
            &prot_coords, &prot_types,
            [0.0, 0.0, 0.0], [8.0, 8.0, 8.0], 0.5,
        );
        let results = dock_from_smiles("O", &map, 10, 7);
        assert!(results.is_some(), "dock_from_smiles should succeed for water");
        assert!(!results.unwrap().is_empty(), "should find at least one valid pose");
    }
}
