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

/// Auto-detect format and count atoms.
/// format: 0=auto, 1=pdb, 2=cif, 3=xyz
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_atom_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    let text = as_str(data, len);
    let atoms = match format {
        1 => pdb_parser::parse_pdb(text, true),
        2 => cif_parser::parse_cif(text, true),
        3 => xyz_parser::parse_xyz(text),
        _ => {
            // Auto-detect
            if text.starts_with("data_") || text.contains("_atom_site.") {
                cif_parser::parse_cif(text, true)
            } else if text.lines().next().map(|l| l.trim().parse::<usize>().is_ok()).unwrap_or(false) {
                xyz_parser::parse_xyz(text)
            } else {
                pdb_parser::parse_pdb(text, true)
            }
        }
    };
    atoms.len() as i32
}

/// Count unique chains in PDB/CIF. Returns count or -1.
#[unsafe(no_mangle)]
pub extern "C" fn ds_structure_chain_count(data: *const u8, len: usize, format: u8) -> i32 {
    if data.is_null() || len == 0 { return -1; }
    let text = as_str(data, len);
    let atoms = match format {
        1 => pdb_parser::parse_pdb(text, true),
        2 => cif_parser::parse_cif(text, true),
        _ => {
            if text.starts_with("data_") || text.contains("_atom_site.") {
                cif_parser::parse_cif(text, true)
            } else {
                pdb_parser::parse_pdb(text, true)
            }
        }
    };
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
    let atoms = match format {
        1 => pdb_parser::parse_pdb(text, false),
        2 => cif_parser::parse_cif(text, false),
        _ => {
            if text.starts_with("data_") || text.contains("_atom_site.") {
                cif_parser::parse_cif(text, false)
            } else {
                pdb_parser::parse_pdb(text, false)
            }
        }
    };
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
    let atoms = match format {
        3 => xyz_parser::parse_xyz(text),
        1 => pdb_parser::parse_pdb(text, true),
        2 => cif_parser::parse_cif(text, true),
        _ => {
            if text.lines().next().map(|l| l.trim().parse::<usize>().is_ok()).unwrap_or(false) {
                xyz_parser::parse_xyz(text)
            } else if text.starts_with("data_") || text.contains("_atom_site.") {
                cif_parser::parse_cif(text, true)
            } else {
                pdb_parser::parse_pdb(text, true)
            }
        }
    };
    let max_model = atoms.iter().map(|a| a.model_num).max().unwrap_or(0);
    max_model
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
}
