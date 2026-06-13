// ── Rigid-body pose representation ───────────────────────────────────────────
// A ligand pose is parameterised by:
//   - translation: (tx, ty, tz)   centre-of-mass position
//   - orientation: unit quaternion (w, x, y, z)
//   - torsions:    Vec<f64>         one angle per rotatable bond
//
// The reference conformer (output of Phase 1) is stored with its centre of
// mass at the origin. Applying the pose translates + rotates it, then applies
// torsion angles along the rotatable bonds.
//
// For Phase 3/4 we implement the rigid case (no torsion application) and note
// where torsion application would be inserted.

// ── Quaternion ────────────────────────────────────────────────────────────────

#[allow(dead_code)]
#[derive(Clone, Copy, Debug)]
pub struct Quat {
    pub w: f64,
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

#[allow(dead_code)]
impl Quat {
    pub fn identity() -> Self {
        Quat { w: 1.0, x: 0.0, y: 0.0, z: 0.0 }
    }

    /// Axis-angle constructor: axis must be normalised.
    pub fn from_axis_angle(ax: f64, ay: f64, az: f64, angle: f64) -> Self {
        let s = (angle / 2.0).sin();
        Quat {
            w: (angle / 2.0).cos(),
            x: ax * s,
            y: ay * s,
            z: az * s,
        }
    }

    pub fn normalise(mut self) -> Self {
        let n = (self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z).sqrt();
        if n > 1e-12 {
            self.w /= n;
            self.x /= n;
            self.y /= n;
            self.z /= n;
        }
        self
    }

    /// Rotate vector (vx,vy,vz) by this quaternion.
    pub fn rotate(&self, vx: f64, vy: f64, vz: f64) -> (f64, f64, f64) {
        // q * v * q^-1  (using expanded formula)
        let (w, x, y, z) = (self.w, self.x, self.y, self.z);
        let rx = (1.0 - 2.0 * (y * y + z * z)) * vx
               + 2.0 * (x * y - z * w) * vy
               + 2.0 * (x * z + y * w) * vz;
        let ry = 2.0 * (x * y + z * w) * vx
               + (1.0 - 2.0 * (x * x + z * z)) * vy
               + 2.0 * (y * z - x * w) * vz;
        let rz = 2.0 * (x * z - y * w) * vx
               + 2.0 * (y * z + x * w) * vy
               + (1.0 - 2.0 * (x * x + y * y)) * vz;
        (rx, ry, rz)
    }

    /// Quaternion multiplication: self ∘ other
    pub fn mul(self, o: Quat) -> Self {
        Quat {
            w: self.w * o.w - self.x * o.x - self.y * o.y - self.z * o.z,
            x: self.w * o.x + self.x * o.w + self.y * o.z - self.z * o.y,
            y: self.w * o.y - self.x * o.z + self.y * o.w + self.z * o.x,
            z: self.w * o.z + self.x * o.y - self.y * o.x + self.z * o.w,
        }
    }
}

// ── Pose ──────────────────────────────────────────────────────────────────────

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct Pose {
    pub translation: [f64; 3],
    pub orientation: Quat,
    pub torsions: Vec<f64>,
}

#[allow(dead_code)]
impl Pose {
    pub fn identity(n_tors: usize) -> Self {
        Pose {
            translation: [0.0, 0.0, 0.0],
            orientation: Quat::identity(),
            torsions: vec![0.0; n_tors],
        }
    }

    /// Degrees of freedom: 3 (translation) + 4 (quaternion) + n_tors
    pub fn n_dof(&self) -> usize {
        7 + self.torsions.len()
    }
}

// ── Apply pose to reference conformer ────────────────────────────────────────
// ref_coords: reference ligand coords with centre of mass at origin.
// Returns the rotated+translated coordinates.
// Torsion application is omitted in this phase (rigid docking).

pub fn apply_pose(ref_coords: &[[f64; 3]], pose: &Pose) -> Vec<[f64; 3]> {
    let (tx, ty, tz) = (pose.translation[0], pose.translation[1], pose.translation[2]);
    ref_coords
        .iter()
        .map(|&[rx, ry, rz]| {
            let (x, y, z) = pose.orientation.rotate(rx, ry, rz);
            [x + tx, y + ty, z + tz]
        })
        .collect()
}

/// Apply a flexible pose: first bend the reference conformation by the pose's
/// torsion angles (via the torsion tree), re-centre the result on the origin,
/// then apply the rigid rotation + translation. Falls back to the rigid path
/// when the tree has no rotatable bonds.
pub fn apply_pose_flex(
    ref_coords: &[[f64; 3]],
    pose: &Pose,
    tree: &crate::docking::torsion::TorsionTree,
) -> Vec<[f64; 3]> {
    if tree.n_tors() == 0 {
        return apply_pose(ref_coords, pose);
    }
    let mut internal = crate::docking::torsion::apply_torsions(ref_coords, tree, &pose.torsions);
    centre_coords(&mut internal);
    let (tx, ty, tz) = (pose.translation[0], pose.translation[1], pose.translation[2]);
    internal
        .iter()
        .map(|&[rx, ry, rz]| {
            let (x, y, z) = pose.orientation.rotate(rx, ry, rz);
            [x + tx, y + ty, z + tz]
        })
        .collect()
}

/// Centre a set of coordinates by subtracting the centroid.
/// Returns the centroid.
pub fn centre_coords(coords: &mut Vec<[f64; 3]>) -> [f64; 3] {
    let n = coords.len() as f64;
    let cx = coords.iter().map(|c| c[0]).sum::<f64>() / n;
    let cy = coords.iter().map(|c| c[1]).sum::<f64>() / n;
    let cz = coords.iter().map(|c| c[2]).sum::<f64>() / n;
    for c in coords.iter_mut() {
        c[0] -= cx;
        c[1] -= cy;
        c[2] -= cz;
    }
    [cx, cy, cz]
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quat_identity_no_rotation() {
        let q = Quat::identity();
        let (rx, ry, rz) = q.rotate(1.0, 2.0, 3.0);
        assert!((rx - 1.0).abs() < 1e-9);
        assert!((ry - 2.0).abs() < 1e-9);
        assert!((rz - 3.0).abs() < 1e-9);
    }

    #[test]
    fn quat_rotate_90_around_z() {
        // 90° around z-axis: (1,0,0) → (0,1,0)
        let q = Quat::from_axis_angle(0.0, 0.0, 1.0, std::f64::consts::FRAC_PI_2);
        let (rx, ry, rz) = q.rotate(1.0, 0.0, 0.0);
        assert!((rx - 0.0).abs() < 1e-9, "rx={rx}");
        assert!((ry - 1.0).abs() < 1e-9, "ry={ry}");
        assert!((rz - 0.0).abs() < 1e-9, "rz={rz}");
    }

    #[test]
    fn quat_normalise_unit() {
        let q = Quat { w: 2.0, x: 0.0, y: 0.0, z: 0.0 }.normalise();
        assert!((q.w - 1.0).abs() < 1e-9);
    }

    #[test]
    fn quat_mul_identity() {
        let q = Quat::from_axis_angle(0.0, 1.0, 0.0, 0.3);
        let id = Quat::identity();
        let qr = q.mul(id);
        assert!((qr.w - q.w).abs() < 1e-9);
        assert!((qr.x - q.x).abs() < 1e-9);
    }

    #[test]
    fn apply_pose_translation_only() {
        let refs: Vec<[f64; 3]> = vec![[1.0, 0.0, 0.0], [-1.0, 0.0, 0.0]];
        let pose = Pose {
            translation: [5.0, 0.0, 0.0],
            orientation: Quat::identity(),
            torsions: vec![],
        };
        let out = apply_pose(&refs, &pose);
        assert!((out[0][0] - 6.0).abs() < 1e-9);
        assert!((out[1][0] - 4.0).abs() < 1e-9);
    }

    #[test]
    fn centre_coords_returns_centroid() {
        let mut coords = vec![[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]];
        let cen = centre_coords(&mut coords);
        assert!((cen[0] - 2.0).abs() < 1e-9);
        assert!((coords[0][0] - (-1.0)).abs() < 1e-9);
        assert!((coords[1][0] - 1.0).abs() < 1e-9);
    }

    #[test]
    fn apply_pose_preserves_distances() {
        let refs: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.54, 0.0, 0.0]];
        let pose = Pose {
            translation: [3.0, 7.0, -2.0],
            orientation: Quat::from_axis_angle(1.0, 0.0, 0.0, 1.234),
            torsions: vec![],
        };
        let out = apply_pose(&refs, &pose);
        let dx = out[0][0] - out[1][0];
        let dy = out[0][1] - out[1][1];
        let dz = out[0][2] - out[1][2];
        let d = (dx * dx + dy * dy + dz * dz).sqrt();
        assert!((d - 1.54).abs() < 1e-9, "distance not preserved: {d}");
    }
}
