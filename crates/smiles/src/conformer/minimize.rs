use crate::conformer::params::{bond_length, ideal_angle, vdw_radius};
use crate::parser::Molecule;

// ── UFF-lite force field ───────────────────────────────────────────────────────
//
// Three energy terms:
//   E_bond  = Σ_{1-2} k_b (r − r₀)²          k_b = 700 kcal/(mol·Å²)
//   E_angle = Σ_{1-3} k_a (θ − θ₀)²          k_a = 100 kcal/(mol·rad²)
//   E_vdw   = Σ_{i<j, 1-4+} ε [(r_min/r)¹² − 2(r_min/r)⁶]  (Lennard-Jones)
//
// Gradient descent with adaptive step-size (Armijo back-tracking line search).
// Terminates when |∇E|_∞ < tol or max_iter reached.

const K_BOND: f64 = 700.0; // kcal / (mol·Å²)
const K_ANGLE: f64 = 100.0; // kcal / (mol·rad²)
const K_OOP: f64 = 400.0; // kcal/(mol·Å⁶) — sp2 out-of-plane (improper) penalty
const VDW_EPS: f64 = 0.10; // kcal/mol  (uniform ε for simplicity)
const MAX_ITER: usize = 2000;
const GRAD_TOL: f64 = 0.2; // kcal/(mol·Å) — convergence tolerance

// ── sp2 planarity (improper) centres ─────────────────────────────────────────
// Each sp2 atom with exactly three neighbours must be coplanar with them. We
// enforce this with a penalty on the signed tetrahedron volume V = u·(v×w),
// where u,v,w are the vectors from the centre to its three neighbours; V = 0 iff
// the centre lies in the plane of the three. This is the ETKDG/MMFF "ET"
// (experimental-torsion / planarity) ingredient missing from raw DG + UFF —
// without it aromatic rings and amide/carbonyl groups pucker out of plane.

#[inline]
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Collect (centre, n0, n1, n2) for every sp2 atom that has exactly three
/// neighbours (aromatic carbons, carbonyl/imine C, amide N, …).
fn sp2_impropers(mol: &Molecule) -> Vec<(usize, usize, usize, usize)> {
    use crate::conformer::params::Hybridization;
    let mut out = Vec::new();
    for c in 0..mol.atoms.len() {
        if mol.hybridization(c) != Hybridization::SP2 {
            continue;
        }
        let nbrs = mol.neighbors(c);
        if nbrs.len() == 3 {
            out.push((c, nbrs[0].0, nbrs[1].0, nbrs[2].0));
        }
    }
    out
}

// ── helpers ───────────────────────────────────────────────────────────────────

#[inline]
fn coord(c: &[f64], i: usize) -> (f64, f64, f64) {
    (c[3 * i], c[3 * i + 1], c[3 * i + 2])
}

#[inline]
fn dist_vec(c: &[f64], i: usize, j: usize) -> (f64, f64, f64, f64) {
    let (xi, yi, zi) = coord(c, i);
    let (xj, yj, zj) = coord(c, j);
    let dx = xi - xj;
    let dy = yi - yj;
    let dz = zi - zj;
    let r = (dx * dx + dy * dy + dz * dz).sqrt().max(1e-8);
    (dx, dy, dz, r)
}

// ── energy + gradient ─────────────────────────────────────────────────────────

fn energy_grad(
    mol: &Molecule,
    c: &[f64],
    grad: &mut [f64],
    impropers: &[(usize, usize, usize, usize)],
) -> f64 {
    let n = mol.atoms.len();
    for g in grad.iter_mut() {
        *g = 0.0;
    }
    let mut e = 0.0;

    // ── Bond (1-2) terms ──────────────────────────────────────────────────
    for b in &mol.bonds {
        let (i, j) = (b.a, b.b);
        let r0 = bond_length(&mol.atoms[i].symbol, &mol.atoms[j].symbol, b.order);
        let (dx, dy, dz, r) = dist_vec(c, i, j);
        let delta = r - r0;
        e += K_BOND * delta * delta;
        let df = 2.0 * K_BOND * delta / r;
        grad[3 * i] += df * dx;
        grad[3 * i + 1] += df * dy;
        grad[3 * i + 2] += df * dz;
        grad[3 * j] -= df * dx;
        grad[3 * j + 1] -= df * dy;
        grad[3 * j + 2] -= df * dz;
    }

    // ── Angle (1-3) terms ─────────────────────────────────────────────────
    for center in 0..n {
        let nbrs = mol.neighbors(center);
        let hyb = mol.hybridization(center);
        let theta0 = ideal_angle(&mol.atoms[center].symbol, hyb).to_radians();

        for ai in 0..nbrs.len() {
            for bi in (ai + 1)..nbrs.len() {
                let (a, _) = nbrs[ai];
                let (b, _) = nbrs[bi];

                // vectors from center → a and center → b
                let (dxa, dya, dza, ra) = dist_vec(c, center, a);
                let (dxb, dyb, dzb, rb) = dist_vec(c, center, b);

                // negate: these go center→a, center→b; we need a→center, b→center
                // but for dot product it's the same sign
                let dot = -dxa * (-dxb) - dya * (-dyb) - dza * (-dzb);
                // Actually: vectors u = a - center, v = b - center
                let ux = -dxa;
                let uy = -dya;
                let uz = -dza; // atom a - center (but dist_vec returns center-a so negate)
                let vx = -dxb;
                let vy = -dyb;
                let vz = -dzb;
                let dot2 = ux * vx + uy * vy + uz * vz;
                let cos_theta = (dot2 / (ra * rb)).clamp(-1.0 + 1e-10, 1.0 - 1e-10);
                let theta = cos_theta.acos();
                let delta = theta - theta0;
                e += K_ANGLE * delta * delta;

                // Gradient via ∂θ/∂r (chain rule through acos)
                let sin_theta = theta.sin().max(1e-8);
                let df = 2.0 * K_ANGLE * delta / sin_theta;

                // ∂cos_θ/∂u_k = v_k/(ra·rb) - cos_θ·u_k/ra²
                // ∂θ/∂u_k    = −∂cos_θ/∂u_k / sin_θ
                let inv_ra_rb = 1.0 / (ra * rb);
                let inv_ra2 = 1.0 / (ra * ra);
                let inv_rb2 = 1.0 / (rb * rb);

                let ga_x = vx * inv_ra_rb - cos_theta * ux * inv_ra2;
                let ga_y = vy * inv_ra_rb - cos_theta * uy * inv_ra2;
                let ga_z = vz * inv_ra_rb - cos_theta * uz * inv_ra2;

                let gb_x = ux * inv_ra_rb - cos_theta * vx * inv_rb2;
                let gb_y = uy * inv_ra_rb - cos_theta * vy * inv_rb2;
                let gb_z = uz * inv_ra_rb - cos_theta * vz * inv_rb2;

                // ∂E/∂a_k = df * ∂θ/∂u_k = -df * ∂cos/∂u_k / sin (handled above)
                // Note: u = a - center → ∂u/∂a = I, ∂u/∂center = -I
                grad[3 * a] -= df * ga_x;
                grad[3 * a + 1] -= df * ga_y;
                grad[3 * a + 2] -= df * ga_z;
                grad[3 * b] -= df * gb_x;
                grad[3 * b + 1] -= df * gb_y;
                grad[3 * b + 2] -= df * gb_z;
                grad[3 * center] += df * (ga_x + gb_x);
                grad[3 * center + 1] += df * (ga_y + gb_y);
                grad[3 * center + 2] += df * (ga_z + gb_z);

                let _ = dot;
            }
        }
    }

    // ── VDW (1-4+ pairs) terms ────────────────────────────────────────────
    // Build 1-2 and 1-3 exclusion sets
    let mut exclude = vec![false; n * n];
    for b in &mol.bonds {
        let (i, j) = (b.a, b.b);
        exclude[i * n + j] = true;
        exclude[j * n + i] = true;
    }
    for center in 0..n {
        let nbrs = mol.neighbors(center);
        for ai in 0..nbrs.len() {
            for bi in (ai + 1)..nbrs.len() {
                let (a, _) = nbrs[ai];
                let (b, _) = nbrs[bi];
                exclude[a * n + b] = true;
                exclude[b * n + a] = true;
            }
        }
    }

    for i in 0..n {
        for j in (i + 1)..n {
            if exclude[i * n + j] {
                continue;
            }
            let r_min =
                (vdw_radius(&mol.atoms[i].symbol) + vdw_radius(&mol.atoms[j].symbol)) * 0.5;
            let (dx, dy, dz, r) = dist_vec(c, i, j);
            let ratio = r_min / r;
            let ratio6 = ratio * ratio * ratio * ratio * ratio * ratio;
            let ratio12 = ratio6 * ratio6;
            e += VDW_EPS * (ratio12 - 2.0 * ratio6);
            // dE/dr = ε (-12 r_min^12 / r^13 + 12 r_min^6 / r^7)
            //       = 12ε/r (r_min/r)^6 [1 - (r_min/r)^6] ... wait
            // dE/dr = 12ε (ratio6 - ratio12) / r
            let de_dr = 12.0 * VDW_EPS * (ratio6 - ratio12) / r;
            let f = de_dr / r;
            grad[3 * i] += f * dx;
            grad[3 * i + 1] += f * dy;
            grad[3 * i + 2] += f * dz;
            grad[3 * j] -= f * dx;
            grad[3 * j + 1] -= f * dy;
            grad[3 * j + 2] -= f * dz;
        }
    }

    // ── sp2 out-of-plane (improper) terms ─────────────────────────────────
    // V = u·(v×w);  E = K_OOP·V²;  ∂V/∂a = v×w, ∂V/∂b = w×u, ∂V/∂d = u×v,
    // ∂V/∂centre = −(those three).
    for &(ctr, a, b, d) in impropers {
        let cc = coord(c, ctr);
        let ca = coord(c, a);
        let cb = coord(c, b);
        let cd = coord(c, d);
        let u = [ca.0 - cc.0, ca.1 - cc.1, ca.2 - cc.2];
        let v = [cb.0 - cc.0, cb.1 - cc.1, cb.2 - cc.2];
        let w = [cd.0 - cc.0, cd.1 - cc.1, cd.2 - cc.2];
        let vxw = cross3(v, w);
        let vol = u[0] * vxw[0] + u[1] * vxw[1] + u[2] * vxw[2];
        e += K_OOP * vol * vol;
        let factor = 2.0 * K_OOP * vol;
        let ga = vxw; // ∂V/∂a
        let gb = cross3(w, u); // ∂V/∂b
        let gd = cross3(u, v); // ∂V/∂d
        for k in 0..3 {
            grad[3 * a + k] += factor * ga[k];
            grad[3 * b + k] += factor * gb[k];
            grad[3 * d + k] += factor * gd[k];
            grad[3 * ctr + k] -= factor * (ga[k] + gb[k] + gd[k]);
        }
    }

    e
}

// ── Gradient descent with Armijo line search ─────────────────────────────────

pub fn minimize(mol: &Molecule, coords: &mut Vec<f64>) -> f64 {
    let n3 = coords.len();
    let mut grad = vec![0.0_f64; n3];
    let mut step = 0.01_f64;
    let impropers = sp2_impropers(mol);
    let mut scratch = vec![0.0_f64; n3];

    for _ in 0..MAX_ITER {
        let e0 = energy_grad(mol, coords, &mut grad, &impropers);

        // Convergence check: max |gradient component|
        let g_max = grad.iter().map(|g| g.abs()).fold(0.0_f64, f64::max);
        if g_max < GRAD_TOL {
            break;
        }

        // Armijo back-tracking
        let mut alpha = step;
        let slope: f64 = grad.iter().map(|g| g * g).sum::<f64>();
        let trial: Vec<f64> = (0..n3)
            .map(|k| coords[k] - alpha * grad[k])
            .collect();
        let mut e_trial = energy_grad(mol, &trial, &mut scratch, &impropers);
        let mut iters = 0;
        while e_trial > e0 - 1e-4 * alpha * slope && iters < 20 {
            alpha *= 0.5;
            let t: Vec<f64> = (0..n3).map(|k| coords[k] - alpha * grad[k]).collect();
            e_trial = energy_grad(mol, &t, &mut scratch, &impropers);
            iters += 1;
        }

        // Accept step
        for k in 0..n3 {
            coords[k] -= alpha * grad[k];
        }
        // Grow step if we accepted on first try
        if iters == 0 {
            step = (step * 1.2).min(0.1);
        }
    }

    // Final energy at the converged geometry.
    energy_grad(mol, coords, &mut grad, &impropers)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::conformer::bounds::build_bounds;
    use crate::conformer::embed::embed;
    use crate::conformer::smoothing::smooth;
    use crate::parser::parse;

    fn gen_conformer(smi: &str, seed: u64) -> (Vec<f64>, Molecule) {
        let mol = parse(smi).unwrap().with_explicit_hydrogens();
        let mut bm = build_bounds(&mol);
        smooth(&mut bm);
        let mut coords = embed(&bm, seed).unwrap();
        minimize(&mol, &mut coords);
        (coords, mol)
    }

    #[test]
    fn minimize_no_nan() {
        let (coords, _) = gen_conformer("CCO", 0);
        for &v in &coords {
            assert!(v.is_finite(), "NaN/Inf after minimize: {v}");
        }
    }

    #[test]
    fn bond_lengths_after_minimize_ethanol() {
        let (coords, mol) = gen_conformer("CCO", 1);
        for b in &mol.bonds {
            let dx = coords[3 * b.a] - coords[3 * b.b];
            let dy = coords[3 * b.a + 1] - coords[3 * b.b + 1];
            let dz = coords[3 * b.a + 2] - coords[3 * b.b + 2];
            let r = (dx * dx + dy * dy + dz * dz).sqrt();
            let r0 = bond_length(
                &mol.atoms[b.a].symbol,
                &mol.atoms[b.b].symbol,
                b.order,
            );
            assert!(
                (r - r0).abs() < 0.3,
                "bond {}-{} ({}-{}) dist={r:.3} ideal={r0:.3}",
                b.a, b.b, mol.atoms[b.a].symbol, mol.atoms[b.b].symbol
            );
        }
    }

    #[test]
    fn minimize_aspirin_finishes() {
        // Just verify no panic and all coords finite.
        let (coords, _) = gen_conformer("CC(=O)Oc1ccccc1C(=O)O", 42);
        for &v in &coords {
            assert!(v.is_finite());
        }
    }

    /// Max perpendicular deviation of ring points (in ring order) from their
    /// best-fit plane, using Newell's robust polygon normal.
    fn planarity_deviation(pts: &[[f64; 3]]) -> f64 {
        let n = pts.len();
        let mut cen = [0.0; 3];
        for p in pts { for k in 0..3 { cen[k] += p[k]; } }
        for k in 0..3 { cen[k] /= n as f64; }
        // Newell's method: robust normal for a (roughly planar) polygon.
        let mut nrm = [0.0_f64; 3];
        for i in 0..n {
            let a = pts[i];
            let b = pts[(i + 1) % n];
            nrm[0] += (a[1] - b[1]) * (a[2] + b[2]);
            nrm[1] += (a[2] - b[2]) * (a[0] + b[0]);
            nrm[2] += (a[0] - b[0]) * (a[1] + b[1]);
        }
        let m = (nrm[0]*nrm[0]+nrm[1]*nrm[1]+nrm[2]*nrm[2]).sqrt().max(1e-12);
        for k in 0..3 { nrm[k] /= m; }
        pts.iter()
            .map(|p| ((p[0]-cen[0])*nrm[0] + (p[1]-cen[1])*nrm[1] + (p[2]-cen[2])*nrm[2]).abs())
            .fold(0.0_f64, f64::max)
    }

    #[test]
    fn benzene_ring_is_planar_after_minimize() {
        // The six aromatic carbons must come out essentially coplanar thanks to
        // the sp2 out-of-plane term. Use the real ensemble path (flattest of N).
        use crate::conformer::generate_conformer;
        let coords = generate_conformer("c1ccccc1", 7).expect("benzene embeds");
        let mol = crate::parser::parse("c1ccccc1").unwrap().with_explicit_hydrogens();
        let ring: Vec<[f64; 3]> = (0..mol.atoms.len())
            .filter(|&i| mol.atoms[i].symbol == "C")
            .map(|i| [coords[3*i], coords[3*i+1], coords[3*i+2]])
            .collect();
        assert_eq!(ring.len(), 6);
        let dev = planarity_deviation(&ring);
        // "lite" target: raw DG embeds benzene with ~0.45 Å pucker; the ring
        // planarity bounds + sp2 out-of-plane term bring it to ≈0.1 Å (a ~4×
        // improvement). Sub-0.01 Å flatness would need L-BFGS + full MMFF.
        assert!(dev < 0.15, "benzene ring not planar: max deviation {dev:.3} Å");
    }

    #[test]
    fn no_vdw_clash_after_minimize() {
        // After minimisation, non-bonded atom pairs should be ≥ 1.0 Å apart.
        let (coords, mol) = gen_conformer("CCCCC", 7);
        let n = mol.atoms.len();
        for i in 0..n {
            for j in (i + 2)..n {
                let dx = coords[3 * i] - coords[3 * j];
                let dy = coords[3 * i + 1] - coords[3 * j + 1];
                let dz = coords[3 * i + 2] - coords[3 * j + 2];
                let r = (dx * dx + dy * dy + dz * dz).sqrt();
                assert!(
                    r >= 1.0,
                    "VDW clash between atoms {i} and {j}: r={r:.3}"
                );
            }
        }
    }
}
