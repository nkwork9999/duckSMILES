pub mod pdb_parser;
pub mod cif_parser;
pub mod xyz_parser;

#[derive(Debug, Clone)]
pub struct Atom {
    pub model_num: i32,
    pub chain: String,
    pub resid: i32,
    pub resname: String,
    pub atom_name: String,
    pub altloc: char,
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub occupancy: f32,
    pub b_factor: f32,
    pub element: String,
    pub is_hetatm: bool,
}

// =============================================================================
// C FFI exports
// =============================================================================

fn as_str(ptr: *const u8, len: usize) -> &'static str {
    unsafe { std::str::from_utf8_unchecked(std::slice::from_raw_parts(ptr, len)) }
}

fn parse_structure(text: &str, format: u8, include_hetatm: bool) -> Vec<Atom> {
    match format {
        1 => pdb_parser::parse_pdb(text, include_hetatm),
        2 => cif_parser::parse_cif(text, include_hetatm),
        3 => xyz_parser::parse_xyz(text),
        _ => {
            if text.starts_with("data_") || text.contains("_atom_site.") {
                cif_parser::parse_cif(text, include_hetatm)
            } else if text.lines().next().map(|l| l.trim().parse::<usize>().is_ok()).unwrap_or(false) {
                xyz_parser::parse_xyz(text)
            } else {
                pdb_parser::parse_pdb(text, include_hetatm)
            }
        }
    }
}

fn centroid(atoms: &[Atom]) -> Option<(f64, f64, f64)> {
    if atoms.is_empty() {
        return None;
    }
    let n = atoms.len() as f64;
    let (sx, sy, sz) = atoms.iter().fold((0.0, 0.0, 0.0), |acc, atom| {
        (
            acc.0 + atom.x as f64,
            acc.1 + atom.y as f64,
            acc.2 + atom.z as f64,
        )
    });
    Some((sx / n, sy / n, sz / n))
}

fn radius_of_gyration(atoms: &[Atom]) -> Option<f64> {
    let (cx, cy, cz) = centroid(atoms)?;
    let n = atoms.len() as f64;
    let sum_sq = atoms.iter().map(|atom| {
        let dx = atom.x as f64 - cx;
        let dy = atom.y as f64 - cy;
        let dz = atom.z as f64 - cz;
        dx * dx + dy * dy + dz * dz
    }).sum::<f64>();
    Some((sum_sq / n).sqrt())
}

fn coordinate_bounds(atoms: &[Atom]) -> Option<[(f64, f64); 3]> {
    let first = atoms.first()?;
    let mut min_x = first.x as f64;
    let mut min_y = first.y as f64;
    let mut min_z = first.z as f64;
    let mut max_x = min_x;
    let mut max_y = min_y;
    let mut max_z = min_z;
    for atom in atoms.iter().skip(1) {
        let x = atom.x as f64;
        let y = atom.y as f64;
        let z = atom.z as f64;
        min_x = min_x.min(x);
        min_y = min_y.min(y);
        min_z = min_z.min(z);
        max_x = max_x.max(x);
        max_y = max_y.max(y);
        max_z = max_z.max(z);
    }
    Some([(min_x, max_x), (min_y, max_y), (min_z, max_z)])
}

/// Auto-detect format and count atoms.
/// format: 0=auto, 1=pdb, 2=cif, 3=xyz
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_atom_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    let text = as_str(data, len);
    let atoms = parse_structure(text, format, true);
    atoms.len() as i32
}

/// Count unique chains in PDB/CIF. Returns count or -1.
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_chain_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    let text = as_str(data, len);
    let atoms = parse_structure(text, format, true);
    let mut chains: Vec<&str> = atoms.iter().map(|a| a.chain.as_str()).collect();
    chains.sort_unstable();
    chains.dedup();
    chains.len() as i32
}

/// Count unique residues. Returns count or -1.
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_residue_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    let text = as_str(data, len);
    let atoms = parse_structure(text, format, false);
    let mut residues: Vec<(String, i32)> = atoms.iter().map(|a| (a.chain.clone(), a.resid)).collect();
    residues.sort();
    residues.dedup();
    residues.len() as i32
}

/// Count models/frames. Returns count.
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_model_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return 0; }
    let text = as_str(data, len);
    let atoms = parse_structure(text, format, true);
    let max_model = atoms.iter().map(|a| a.model_num).max().unwrap_or(0);
    max_model
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_centroid_x(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    centroid(&atoms).map(|c| c.0).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_centroid_y(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    centroid(&atoms).map(|c| c.1).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_centroid_z(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    centroid(&atoms).map(|c| c.2).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_radius_of_gyration(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    radius_of_gyration(&atoms).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_min_x(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[0].0).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_max_x(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[0].1).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_min_y(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[1].0).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_max_y(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[1].1).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_min_z(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[2].0).unwrap_or(f64::NAN)
}

#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_max_z(data: *const u8, len: usize, format: u8) -> f64 {
    if data.is_null() || len == 0 { return f64::NAN; }
    let atoms = parse_structure(as_str(data, len), format, true);
    coordinate_bounds(&atoms).map(|b| b[2].1).unwrap_or(f64::NAN)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_pdb_atom_count() {
        let pdb = b"ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
        let count = ds_structure_atom_count(pdb.as_ptr(), pdb.len(), 1);
        assert_eq!(count, 1);
    }

    #[test]
    fn test_xyz_atom_count() {
        let xyz = b"3\nwater\nO  0.000  0.000  0.117\nH  0.000  0.756 -0.469\nH  0.000 -0.756 -0.469\n";
        let count = ds_structure_atom_count(xyz.as_ptr(), xyz.len(), 3);
        assert_eq!(count, 3);
    }

    #[test]
    fn test_xyz_model_count() {
        let xyz = b"2\nframe1\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n2\nframe2\nC 0.0 0.0 0.5\nO 1.0 0.0 0.5\n";
        let count = ds_structure_model_count(xyz.as_ptr(), xyz.len(), 3);
        assert_eq!(count, 2);
    }

    #[test]
    fn test_xyz_centroid_and_radius() {
        let xyz = b"2\nco\nC 0.0 0.0 0.0\nO 2.0 0.0 0.0\n";
        assert_eq!(ds_structure_centroid_x(xyz.as_ptr(), xyz.len(), 3), 1.0);
        assert_eq!(ds_structure_centroid_y(xyz.as_ptr(), xyz.len(), 3), 0.0);
        assert_eq!(ds_structure_centroid_z(xyz.as_ptr(), xyz.len(), 3), 0.0);
        assert_eq!(ds_structure_radius_of_gyration(xyz.as_ptr(), xyz.len(), 3), 1.0);
    }

    #[test]
    fn test_xyz_bounds() {
        let xyz = b"3\ntri\nC -1.0 2.0 0.5\nO 3.0 -2.0 1.5\nN 0.0 1.0 -4.0\n";
        assert_eq!(ds_structure_min_x(xyz.as_ptr(), xyz.len(), 3), -1.0);
        assert_eq!(ds_structure_max_x(xyz.as_ptr(), xyz.len(), 3), 3.0);
        assert_eq!(ds_structure_min_y(xyz.as_ptr(), xyz.len(), 3), -2.0);
        assert_eq!(ds_structure_max_y(xyz.as_ptr(), xyz.len(), 3), 2.0);
        assert_eq!(ds_structure_min_z(xyz.as_ptr(), xyz.len(), 3), -4.0);
        assert_eq!(ds_structure_max_z(xyz.as_ptr(), xyz.len(), 3), 1.5);
    }
}
