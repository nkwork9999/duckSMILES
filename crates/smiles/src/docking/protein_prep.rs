use crate::docking::atomtype::VinaType;
use crate::docking::pdb::{parse_pdb, PdbAtom};

// ── Protein preparation (protonation + polar-H addition) ─────────────────────
//
// The raw PDB→PDBQT path in pdbqt.rs types heavy atoms but never adds the polar
// hydrogens a protein needs to *donate* H-bonds. Without HD atoms a docked
// ligand can only accept H-bonds from the protein backbone/side chains, never
// receive them — a large scoring gap. This module is the "Reduce + PROPKA-lite"
// step: it assigns standard protonation states at a given pH and places polar
// hydrogens geometrically so the affinity map carries real donor points.
//
// Scope / honest limits:
//  - Protonation states are the standard physiological assignments keyed on
//    residue + a pKa cutoff, NOT a full Poisson-Boltzmann / hydrogen-bond
//    network optimisation. Asn/Gln/His "flips" are not searched (HIE tautomer
//    is assumed for neutral His).
//  - Hydrogen placement uses idealised cone geometry from the donor's heavy
//    neighbours; azimuthal (rotatable O-H/N-H) orientation is not optimised
//    against acceptors. This is sufficient for a grid-based Vina map but is not
//    a force-field H-bond network minimisation.
//  - Termini: backbone N is protonated with a single amide H uniformly; the
//    extra N-terminal NH3+ hydrogens are not added.

const H_BOND_N: f64 = 1.01;
const H_BOND_O: f64 = 0.96;
const H_BOND_S: f64 = 1.34;

/// A prepared protein atom: 3D position + assigned Vina type. Added polar
/// hydrogens appear as extra entries with `VinaType::HD`.
#[derive(Clone, Debug)]
pub struct PreparedAtom {
    pub pos: [f64; 3],
    pub vtype: VinaType,
    pub name: String,
    pub res_name: String,
}

/// Build a prepared protein from PDB text at the given pH.
///
/// Steps: drop solvent (HOH/WAT), assign pH-aware heavy-atom Vina types, then
/// add polar hydrogens (HD) for every donor heavy atom that does not already
/// carry one in the input.
pub fn prepare_protein(pdb_text: &str, ph: f64) -> Vec<PreparedAtom> {
    let raw = parse_pdb(pdb_text);

    // Heavy atoms (skip solvent and existing H — we re-derive H ourselves).
    let heavy: Vec<&PdbAtom> = raw
        .iter()
        .filter(|a| !is_solvent(&a.res_name))
        .filter(|a| a.element != "H" && a.element != "D")
        .collect();

    let mut out: Vec<PreparedAtom> = Vec::with_capacity(heavy.len());

    // Type every heavy atom (pH-aware).
    for a in &heavy {
        let vtype = heavy_type(&a.element, a.name.trim(), a.res_name.trim(), ph, a.is_hetatm);
        out.push(PreparedAtom {
            pos: [a.x, a.y, a.z],
            vtype,
            name: a.name.trim().to_string(),
            res_name: a.res_name.trim().to_string(),
        });
    }

    // Add polar hydrogens for donors.
    let heavy_positions: Vec<[f64; 3]> = heavy.iter().map(|a| [a.x, a.y, a.z]).collect();
    let mut added: Vec<PreparedAtom> = Vec::new();

    for (i, a) in heavy.iter().enumerate() {
        if let Some((n_h, h_elem)) = donor_spec(a.res_name.trim(), a.name.trim(), &a.element, ph) {
            let donor = [a.x, a.y, a.z];
            let neighbours = heavy_neighbours(&heavy_positions, i, donor);
            let bond_len = match h_elem {
                "N" => H_BOND_N,
                "O" => H_BOND_O,
                "S" => H_BOND_S,
                _ => H_BOND_N,
            };
            for hpos in place_polar_hydrogens(donor, &neighbours, n_h, bond_len) {
                added.push(PreparedAtom {
                    pos: hpos,
                    vtype: VinaType::HD,
                    name: "H".to_string(),
                    res_name: a.res_name.trim().to_string(),
                });
            }
        }
    }

    out.extend(added);
    out
}

/// Convenience: prepared protein as parallel coordinate / type vectors for the
/// affinity-map builder.
pub fn prepared_coords_types(prepared: &[PreparedAtom]) -> (Vec<[f64; 3]>, Vec<VinaType>) {
    let coords = prepared.iter().map(|a| a.pos).collect();
    let types = prepared.iter().map(|a| a.vtype).collect();
    (coords, types)
}

/// Serialise a prepared protein to PDBQT text (one ATOM record per atom, with
/// the Vina type in the trailing column). Added polar hydrogens appear as HD.
pub fn prepared_to_pdbqt(prepared: &[PreparedAtom]) -> String {
    let mut s = String::with_capacity(prepared.len() * 80);
    for (i, a) in prepared.iter().enumerate() {
        let serial = (i + 1) % 100000;
        let name = if a.name.len() <= 4 { a.name.clone() } else { a.name[..4].to_string() };
        let res = if a.res_name.len() <= 3 { a.res_name.clone() } else { a.res_name[..3].to_string() };
        // Columns follow the ligand PDBQT writer in pdbqt.rs.
        s.push_str(&format!(
            "ATOM  {serial:>5} {name:<4} {res:<3} A   1    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00    {charge:>6.3} {vt}\n",
            serial = serial, name = name, res = res,
            x = a.pos[0], y = a.pos[1], z = a.pos[2],
            charge = 0.0, vt = a.vtype.as_str(),
        ));
    }
    s
}

// ── solvent / residue helpers ───────────────────────────────────────────────

fn is_solvent(res: &str) -> bool {
    matches!(res.trim(), "HOH" | "WAT" | "DOD" | "H2O")
}

// ── heavy-atom Vina typing (pH-aware) ───────────────────────────────────────

fn heavy_type(element: &str, name: &str, res: &str, ph: f64, is_hetatm: bool) -> VinaType {
    match element {
        "C" => carbon_type(name, res),
        "O" => VinaType::OA, // backbone/side-chain O are acceptors (and OH also donate)
        "S" => VinaType::SA,
        "N" => nitrogen_type(name, res, ph),
        "P" => VinaType::P,
        "F" => VinaType::F,
        "CL" | "Cl" => VinaType::CL,
        "BR" | "Br" => VinaType::BR,
        "I" => VinaType::I,
        "MG" | "Mg" => VinaType::MG,
        "ZN" | "Zn" => VinaType::ZN,
        "FE" | "Fe" => VinaType::FE,
        "MN" | "Mn" => VinaType::MN,
        "CU" | "Cu" => VinaType::CU,
        "CA" | "Ca" if is_hetatm => VinaType::CA,
        _ => VinaType::Unknown,
    }
}

fn carbon_type(name: &str, res: &str) -> VinaType {
    let arom: &[(&str, &[&str])] = &[
        ("PHE", &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
        ("TYR", &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
        ("TRP", &["CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]),
        ("HIS", &["CG", "CD2", "CE1"]),
        ("HID", &["CG", "CD2", "CE1"]),
        ("HIE", &["CG", "CD2", "CE1"]),
        ("HIP", &["CG", "CD2", "CE1"]),
    ];
    for (r, atoms) in arom {
        if res == *r && atoms.contains(&name) {
            return VinaType::A;
        }
    }
    VinaType::C
}

/// Nitrogen type after protonation. Donor-only N is generic `N` (its donation
/// is carried by the added HD); only genuine lone-pair acceptors get `NA`.
fn nitrogen_type(name: &str, res: &str, ph: f64) -> VinaType {
    // Neutral His (pH above ~6) in the HIE tautomer: ND1 is a lone-pair acceptor.
    if (res == "HIS" || res == "HIE") && name == "ND1" && ph >= 6.0 {
        return VinaType::NA;
    }
    if res == "HID" && name == "NE2" && ph >= 6.0 {
        return VinaType::NA;
    }
    // Everything else that is a recognised donor (backbone, Lys, Arg, Asn,
    // Gln, Trp, protonated His) becomes generic N; the HD does the donating.
    VinaType::N
}

// ── donor specification (which heavy atoms get H, how many) ──────────────────
//
// Returns Some((n_hydrogens, h_bond_element)) for donor heavy atoms at the
// given pH; None for acceptors / non-donors.

fn donor_spec(res: &str, name: &str, element: &str, ph: f64) -> Option<(usize, &'static str)> {
    // Backbone amide N (every residue except proline) donates one H.
    if name == "N" && res != "PRO" && element == "N" {
        return Some((1, "N"));
    }

    match (res, name) {
        // Hydroxyls: donor + acceptor, one H.
        ("SER", "OG") => Some((1, "O")),
        ("THR", "OG1") => Some((1, "O")),
        ("TYR", "OH") => Some((1, "O")),
        // Thiol: weak donor, one H.
        ("CYS", "SG") => Some((1, "S")),
        // Primary amides: NH2 (two H, planar).
        ("ASN", "ND2") => Some((2, "N")),
        ("GLN", "NE2") => Some((2, "N")),
        // Indole NH (Trp): one H.
        ("TRP", "NE1") => Some((1, "N")),
        // Lysine ammonium NH3+ (protonated below pKa ~10.5): three H.
        ("LYS", "NZ") if ph < 10.5 => Some((3, "N")),
        // Arginine guanidinium (protonated below pKa ~12.5).
        ("ARG", "NE") if ph < 12.5 => Some((1, "N")),
        ("ARG", "NH1") if ph < 12.5 => Some((2, "N")),
        ("ARG", "NH2") if ph < 12.5 => Some((2, "N")),
        // Histidine: neutral HIE tautomer donates from NE2; if protonated
        // (HIP, pH below ~6) both ring N's donate.
        ("HIS", "NE2") | ("HIE", "NE2") if ph >= 6.0 => Some((1, "N")),
        ("HIS", "ND1") | ("HID", "ND1") if ph >= 6.0 => Some((1, "N")),
        ("HIS", "ND1") | ("HIS", "NE2") | ("HIP", "ND1") | ("HIP", "NE2") if ph < 6.0 => {
            Some((1, "N"))
        }
        // Asp/Glu carboxylates: protonated (donor) only below pKa ~4.
        ("ASP", "OD1") | ("ASP", "OD2") if ph < 4.0 => Some((1, "O")),
        ("GLU", "OE1") | ("GLU", "OE2") if ph < 4.0 => Some((1, "O")),
        _ => None,
    }
}

// ── geometry: find heavy neighbours of a donor ──────────────────────────────

const BOND_CUTOFF: f64 = 1.85; // Å; covalent neighbour cutoff for heavy atoms

fn heavy_neighbours(positions: &[[f64; 3]], donor_idx: usize, donor: [f64; 3]) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for (j, p) in positions.iter().enumerate() {
        if j == donor_idx {
            continue;
        }
        let d = dist(donor, *p);
        if d > 0.1 && d < BOND_CUTOFF {
            out.push(*p);
        }
    }
    out
}

// ── geometry: place polar hydrogens on a cone around the open direction ──────

/// Place `n_h` hydrogens at distance `bond_len` from `donor`. The "open"
/// direction is opposite the average bond direction to the heavy neighbours;
/// hydrogens are spread on a cone about it (half-angle by H count) and evenly
/// in azimuth. When there are no heavy neighbours a default axis is used.
pub fn place_polar_hydrogens(
    donor: [f64; 3],
    neighbours: &[[f64; 3]],
    n_h: usize,
    bond_len: f64,
) -> Vec<[f64; 3]> {
    if n_h == 0 {
        return Vec::new();
    }

    // Open axis = -normalize(sum of unit vectors toward neighbours).
    let mut axis = [0.0_f64; 3];
    for nb in neighbours {
        let v = norm3(sub(*nb, donor));
        axis = [axis[0] + v[0], axis[1] + v[1], axis[2] + v[2]];
    }
    axis = if mag3(axis) > 1e-6 {
        let a = norm3(axis);
        [-a[0], -a[1], -a[2]]
    } else {
        [0.0, 0.0, 1.0] // no neighbours: arbitrary up
    };

    // Cone half-angle by hydrogen count (sp2 planar vs sp3 tetrahedral).
    let half_angle_deg: f64 = match n_h {
        1 => 0.0,    // straight along the open axis
        2 => 60.0,   // planar NH2 (≈120° at N)
        _ => 70.5,   // tetrahedral NH3+ (≈109.5° at N)
    };
    let half = half_angle_deg.to_radians();

    // Perpendicular basis (u, w) spanning the plane normal to `axis`.
    let (u, w) = perpendicular_basis(axis);

    let mut hs = Vec::with_capacity(n_h);
    for k in 0..n_h {
        let phi = if n_h == 1 {
            0.0
        } else {
            2.0 * std::f64::consts::PI * (k as f64) / (n_h as f64)
        };
        // direction = cos(half)*axis + sin(half)*(cosφ*u + sinφ*w)
        let (sh, ch) = half.sin_cos();
        let (sp, cp) = phi.sin_cos();
        let dir = [
            ch * axis[0] + sh * (cp * u[0] + sp * w[0]),
            ch * axis[1] + sh * (cp * u[1] + sp * w[1]),
            ch * axis[2] + sh * (cp * u[2] + sp * w[2]),
        ];
        let dir = norm3(dir);
        hs.push([
            donor[0] + bond_len * dir[0],
            donor[1] + bond_len * dir[1],
            donor[2] + bond_len * dir[2],
        ]);
    }
    hs
}

fn perpendicular_basis(axis: [f64; 3]) -> ([f64; 3], [f64; 3]) {
    // Pick a reference not parallel to axis.
    let reference = if axis[0].abs() < 0.9 {
        [1.0, 0.0, 0.0]
    } else {
        [0.0, 1.0, 0.0]
    };
    let u = norm3(cross(axis, reference));
    let w = norm3(cross(axis, u));
    (u, w)
}

// ── small vector helpers ────────────────────────────────────────────────────

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}
fn mag3(a: [f64; 3]) -> f64 {
    (a[0] * a[0] + a[1] * a[1] + a[2] * a[2]).sqrt()
}
fn norm3(a: [f64; 3]) -> [f64; 3] {
    let m = mag3(a);
    if m < 1e-12 {
        [0.0, 0.0, 0.0]
    } else {
        [a[0] / m, a[1] / m, a[2] / m]
    }
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn dist(a: [f64; 3], b: [f64; 3]) -> f64 {
    mag3(sub(a, b))
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // Tiny PDB: one serine residue backbone + OG, plus an aspartate side chain.
    const SER_ASP_PDB: &str = "\
ATOM      1  N   SER A   1       0.000   0.000   0.000  1.00  0.00           N  \n\
ATOM      2  CA  SER A   1       1.450   0.000   0.000  1.00  0.00           C  \n\
ATOM      3  C   SER A   1       2.000   1.400   0.000  1.00  0.00           C  \n\
ATOM      4  O   SER A   1       1.300   2.400   0.000  1.00  0.00           O  \n\
ATOM      5  CB  SER A   1       2.000  -1.000   1.000  1.00  0.00           C  \n\
ATOM      6  OG  SER A   1       3.400  -1.000   1.000  1.00  0.00           O  \n\
ATOM      7  CG  ASP B   2      10.000   0.000   0.000  1.00  0.00           C  \n\
ATOM      8  OD1 ASP B   2      11.000   0.700   0.000  1.00  0.00           O  \n\
ATOM      9  OD2 ASP B   2      11.000  -0.700   0.000  1.00  0.00           O  \n\
";

    #[test]
    fn adds_polar_hydrogens_as_hd() {
        let prepared = prepare_protein(SER_ASP_PDB, 7.4);
        let n_hd = prepared.iter().filter(|a| a.vtype == VinaType::HD).count();
        // Backbone N (1 H) + Ser OG (1 H) = 2 donor H at pH 7.4.
        // Asp carboxylate is deprotonated → no H.
        assert_eq!(n_hd, 2, "expected 2 polar H (backbone N + Ser OG)");
    }

    #[test]
    fn aspartate_deprotonated_at_physiological_ph() {
        let prepared = prepare_protein(SER_ASP_PDB, 7.4);
        // Both Asp oxygens present as OA acceptors, no HD attached to them.
        let asp_o = prepared
            .iter()
            .filter(|a| a.res_name == "ASP" && a.vtype == VinaType::OA)
            .count();
        assert_eq!(asp_o, 2, "both Asp carboxylate O should be OA acceptors");
    }

    #[test]
    fn aspartate_protonated_at_low_ph() {
        let prepared = prepare_protein(SER_ASP_PDB, 2.0);
        // At pH 2 the carboxylate is protonated → at least one extra HD appears.
        let n_hd = prepared.iter().filter(|a| a.vtype == VinaType::HD).count();
        assert!(n_hd >= 3, "low pH should protonate Asp (more HD): {n_hd}");
    }

    #[test]
    fn solvent_is_removed() {
        let pdb = "\
ATOM      1  N   SER A   1       0.000   0.000   0.000  1.00  0.00           N  \n\
HETATM    2  O   HOH A 100       5.000   5.000   5.000  1.00  0.00           O  \n\
";
        let prepared = prepare_protein(pdb, 7.4);
        assert!(
            prepared.iter().all(|a| a.res_name != "HOH"),
            "water must be stripped"
        );
    }

    #[test]
    fn hydrogen_at_correct_bond_length() {
        // O-H placed at ~0.96 Å from the Ser OG.
        let prepared = prepare_protein(SER_ASP_PDB, 7.4);
        let og = prepared.iter().find(|a| a.name == "OG").unwrap().pos;
        let hd: Vec<_> = prepared.iter().filter(|a| a.vtype == VinaType::HD).collect();
        // The H nearest to OG should be ~0.96 Å away.
        let nearest = hd
            .iter()
            .map(|h| dist(h.pos, og))
            .fold(f64::INFINITY, f64::min);
        assert!((nearest - 0.96).abs() < 0.05, "O-H length wrong: {nearest}");
    }

    #[test]
    fn lysine_gets_three_hydrogens() {
        let pdb = "\
ATOM      1  CE  LYS A   1       0.000   0.000   0.000  1.00  0.00           C  \n\
ATOM      2  NZ  LYS A   1       1.500   0.000   0.000  1.00  0.00           N  \n\
";
        let prepared = prepare_protein(pdb, 7.4);
        let n_hd = prepared.iter().filter(|a| a.vtype == VinaType::HD).count();
        assert_eq!(n_hd, 3, "Lys NZ should get 3 ammonium hydrogens");
        // Each H ~1.01 Å from NZ.
        let nz = prepared.iter().find(|a| a.name == "NZ").unwrap().pos;
        for h in prepared.iter().filter(|a| a.vtype == VinaType::HD) {
            assert!((dist(h.pos, nz) - 1.01).abs() < 0.05, "N-H length");
        }
    }

    #[test]
    fn placed_hydrogen_points_away_from_neighbour() {
        // Single neighbour at +x; the 1-H donor should place H on the -x side.
        let donor = [0.0, 0.0, 0.0];
        let neigh = vec![[1.5, 0.0, 0.0]];
        let hs = place_polar_hydrogens(donor, &neigh, 1, 1.0);
        assert_eq!(hs.len(), 1);
        assert!(hs[0][0] < 0.0, "H should point away from the neighbour: {:?}", hs[0]);
    }

    #[test]
    fn two_hydrogens_are_distinct_and_correct_length() {
        let donor = [0.0, 0.0, 0.0];
        let neigh = vec![[1.5, 0.0, 0.0]];
        let hs = place_polar_hydrogens(donor, &neigh, 2, 1.01);
        assert_eq!(hs.len(), 2);
        assert!(dist(hs[0], hs[1]) > 0.5, "the two H should be separated");
        for h in &hs {
            assert!((dist(*h, donor) - 1.01).abs() < 1e-6, "N-H length");
        }
    }
}
