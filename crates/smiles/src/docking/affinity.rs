use crate::docking::atomtype::VinaType;
use crate::docking::score::{pairwise, vdw};

// ── 3D Affinity (interaction) Maps ───────────────────────────────────────────
// Pre-computes the interaction energy between each grid point and all protein
// atoms, for a fixed set of probe atom types. During docking search, each
// ligand atom's energy is looked up via trilinear interpolation.
//
// Grid layout: (type, ix, iy, iz), flat row-major: ix*(ny*nz) + iy*nz + iz.

// Probe atom types for which maps are built.
// Order must match PROBE_TYPES below for index lookups.
pub const PROBE_TYPES: &[VinaType] = &[
    VinaType::C,
    VinaType::A,
    VinaType::N,
    VinaType::NA,
    VinaType::NS,
    VinaType::OA,
    VinaType::SA,
    VinaType::HD,
    VinaType::P,
    VinaType::F,
    VinaType::CL,
    VinaType::BR,
    VinaType::I,
    VinaType::S,
];

pub const N_PROBE: usize = PROBE_TYPES.len();

pub fn probe_index(t: VinaType) -> Option<usize> {
    PROBE_TYPES.iter().position(|&p| p == t)
}

pub struct AffinityMap {
    pub center: [f64; 3],
    pub origin: [f64; 3], // actual lower corner of grid = center - size
    pub n: [usize; 3],
    pub spacing: f64,
    /// Flat storage: data[probe_idx][ix * ny * nz + iy * nz + iz]
    data: Vec<Vec<f64>>,
}

impl AffinityMap {
    /// Build affinity maps from protein atoms.
    /// center: binding-site center (Å)
    /// size:   box half-widths in each dimension (Å)
    /// spacing: grid spacing (Å), default 0.375 Å in Vina
    pub fn build(
        prot_coords: &[[f64; 3]],
        prot_types: &[VinaType],
        center: [f64; 3],
        size: [f64; 3],
        spacing: f64,
    ) -> Self {
        let origin = [
            center[0] - size[0],
            center[1] - size[1],
            center[2] - size[2],
        ];
        let n = [
            ((2.0 * size[0]) / spacing).ceil() as usize + 1,
            ((2.0 * size[1]) / spacing).ceil() as usize + 1,
            ((2.0 * size[2]) / spacing).ceil() as usize + 1,
        ];
        let n_pts = n[0] * n[1] * n[2];

        let mut data: Vec<Vec<f64>> = (0..N_PROBE).map(|_| vec![0.0; n_pts]).collect();

        for ix in 0..n[0] {
            let gx = origin[0] + ix as f64 * spacing;
            for iy in 0..n[1] {
                let gy = origin[1] + iy as f64 * spacing;
                for iz in 0..n[2] {
                    let gz = origin[2] + iz as f64 * spacing;
                    let idx = ix * n[1] * n[2] + iy * n[2] + iz;

                    for (pidx, &probe) in PROBE_TYPES.iter().enumerate() {
                        let mut e = 0.0;
                        for (pc, &pt) in prot_coords.iter().zip(prot_types.iter()) {
                            let dx = gx - pc[0];
                            let dy = gy - pc[1];
                            let dz = gz - pc[2];
                            let r = (dx * dx + dy * dy + dz * dz).sqrt();
                            // Include only interactions within 8 + max_vdw cutoff
                            let d = r - vdw(probe) - vdw(pt);
                            if d <= 8.0 {
                                e += pairwise(probe, pt, r);
                            }
                        }
                        data[pidx][idx] = e;
                    }
                }
            }
        }

        AffinityMap { center, origin, n, spacing, data }
    }

    /// Trilinear interpolation of the affinity map for atom type `t` at
    /// position `pos`. Returns `None` if pos is outside the grid.
    pub fn interpolate(&self, t: VinaType, pos: [f64; 3]) -> Option<f64> {
        let pidx = probe_index(t)?;

        // Convert pos to fractional grid coordinates
        let fx = (pos[0] - self.origin[0]) / self.spacing;
        let fy = (pos[1] - self.origin[1]) / self.spacing;
        let fz = (pos[2] - self.origin[2]) / self.spacing;

        let ix = fx as usize;
        let iy = fy as usize;
        let iz = fz as usize;

        if ix + 1 >= self.n[0] || iy + 1 >= self.n[1] || iz + 1 >= self.n[2] {
            return None; // outside grid
        }

        let tx = fx - ix as f64;
        let ty = fy - iy as f64;
        let tz = fz - iz as f64;

        let v = |x: usize, y: usize, z: usize| -> f64 {
            self.data[pidx][x * self.n[1] * self.n[2] + y * self.n[2] + z]
        };

        // Trilinear interpolation
        let e = (1.0 - tx) * (1.0 - ty) * (1.0 - tz) * v(ix, iy, iz)
            + tx       * (1.0 - ty) * (1.0 - tz) * v(ix+1, iy, iz)
            + (1.0 - tx) * ty       * (1.0 - tz) * v(ix, iy+1, iz)
            + (1.0 - tx) * (1.0 - ty) * tz       * v(ix, iy, iz+1)
            + tx       * ty       * (1.0 - tz) * v(ix+1, iy+1, iz)
            + tx       * (1.0 - ty) * tz       * v(ix+1, iy, iz+1)
            + (1.0 - tx) * ty       * tz       * v(ix, iy+1, iz+1)
            + tx       * ty       * tz       * v(ix+1, iy+1, iz+1);

        Some(e)
    }

    /// Score all ligand atoms against the pre-computed maps.
    ///
    /// Non-polar hydrogens (and any atom type without a probe map, e.g. plain
    /// `H`) are skipped: like AutoDock Vina's united-atom treatment, their
    /// contribution is folded into the parent heavy atom and they are not
    /// scored individually. Returns `None` only if a *scored* (heavy / polar)
    /// atom falls outside the grid.
    pub fn score_ligand(
        &self,
        lig_coords: &[[f64; 3]],
        lig_types: &[VinaType],
        n_rot: u32,
    ) -> Option<f64> {
        let mut e = 0.0;
        for (c, &t) in lig_coords.iter().zip(lig_types.iter()) {
            if probe_index(t).is_none() {
                continue; // non-polar H / unmapped type → not scored
            }
            e += self.interpolate(t, *c)?;
        }
        Some(e / (1.0 + crate::docking::score::W_ROT * n_rot as f64))
    }

    #[allow(dead_code)]
    pub fn n_pts(&self) -> usize {
        self.n[0] * self.n[1] * self.n[2]
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::docking::score::total_score;

    fn single_carbon_map() -> AffinityMap {
        let prot_coords = vec![[0.0, 0.0, 0.0_f64]];
        let prot_types = vec![VinaType::C];
        AffinityMap::build(
            &prot_coords,
            &prot_types,
            [5.0, 5.0, 5.0],
            [6.0, 6.0, 6.0],
            0.5,
        )
    }

    #[test]
    fn grid_has_correct_dimensions() {
        let m = single_carbon_map();
        // size 6.0 / spacing 0.5 = 12 intervals → 13 points per axis (+1)
        assert_eq!(m.n[0], 25);
        assert_eq!(m.n[1], 25);
        assert_eq!(m.n[2], 25);
    }

    #[test]
    fn interpolate_at_grid_centre() {
        // Grid centred at (5,5,5). Protein C at (0,0,0).
        // The probe C at (5,5,5) is ~8.66 Å from origin — within cutoff range.
        let m = single_carbon_map();
        let val = m.interpolate(VinaType::C, [5.0, 5.0, 5.0]);
        assert!(val.is_some(), "grid centre should be within bounds");
        assert!(val.unwrap().is_finite());
    }

    #[test]
    fn outside_grid_returns_none() {
        let m = single_carbon_map();
        let val = m.interpolate(VinaType::C, [100.0, 100.0, 100.0]);
        assert!(val.is_none(), "far-away point should be None");
    }

    #[test]
    fn probe_index_roundtrip() {
        for (i, &t) in PROBE_TYPES.iter().enumerate() {
            assert_eq!(probe_index(t), Some(i));
        }
        assert_eq!(probe_index(VinaType::Unknown), None);
    }

    #[test]
    fn map_vs_direct_at_grid_point() {
        // Verify that the map value at an exact grid point matches direct calculation.
        let prot_coords = vec![[0.0, 0.0, 0.0_f64]];
        let prot_types = vec![VinaType::C];
        let center = [4.0, 4.0, 4.0];
        let size = [5.0, 5.0, 5.0];
        let spacing = 1.0;
        let map = AffinityMap::build(&prot_coords, &prot_types, center, size, spacing);

        // Probe C at exactly grid point (ix=1, iy=1, iz=1)
        // origin is stored in map.origin
        let gx = map.origin[0] + 1.0 * spacing;
        let gy = map.origin[1] + 1.0 * spacing;
        let gz = map.origin[2] + 1.0 * spacing;

        let map_val = map.interpolate(VinaType::C, [gx, gy, gz]).unwrap();
        let direct = total_score(
            &[[gx, gy, gz]],
            &[VinaType::C],
            &prot_coords,
            &prot_types,
            0,
        );

        assert!(
            (map_val - direct).abs() < 1e-9,
            "map={map_val} vs direct={direct} at ({gx},{gy},{gz})"
        );
    }

    #[test]
    fn score_ligand_returns_none_outside_grid() {
        let m = single_carbon_map();
        let result = m.score_ligand(&[[100.0, 0.0, 0.0]], &[VinaType::C], 0);
        assert!(result.is_none());
    }

    #[test]
    fn score_ligand_inside_grid() {
        let m = single_carbon_map();
        let result = m.score_ligand(&[[5.0, 5.0, 5.0]], &[VinaType::C], 0);
        assert!(result.is_some());
        assert!(result.unwrap().is_finite());
    }
}
