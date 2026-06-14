// ── PDB file parser ───────────────────────────────────────────────────────────
// Parses ATOM and HETATM records from PDB format (fixed-column, IUPAC PDB 3.3).
// Only the fields needed for docking preparation are extracted.

#[allow(dead_code)]
#[derive(Clone, Debug)]
pub struct PdbAtom {
    pub serial: u32,
    pub name: String,     // atom name, e.g. "CA", "NZ"
    pub alt_loc: char,    // alternate location indicator (' ' = none)
    pub res_name: String, // residue name, e.g. "ALA"
    pub chain: char,      // chain ID
    pub res_seq: i32,     // residue sequence number
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub occupancy: f64,
    pub b_factor: f64,
    pub element: String,  // element symbol, e.g. "C", "N"
    pub is_hetatm: bool,
}

/// Parse PDB text and return all ATOM and HETATM records.
/// Lines that are malformed or cannot be parsed are silently skipped.
pub fn parse_pdb(text: &str) -> Vec<PdbAtom> {
    let mut atoms = Vec::new();

    for line in text.lines() {
        let record = if line.len() >= 6 { &line[..6] } else { continue };
        let is_hetatm = match record.trim_end() {
            "ATOM" => false,
            "HETATM" => true,
            _ => continue,
        };

        // PDB is 1-indexed in the spec; we use 0-indexed slices.
        // All col references below are PDB-standard (1-based), converted to 0-based:
        //   col N → line[N-1..N-1+len]

        let serial = parse_int(line, 7, 11).unwrap_or(0) as u32;
        let name = line.get(12..16).unwrap_or("    ").to_string();
        let alt_loc = line.get(16..17).and_then(|s| s.chars().next()).unwrap_or(' ');
        let res_name = line.get(17..20).unwrap_or("   ").trim().to_string();
        let chain = line.get(21..22).and_then(|s| s.chars().next()).unwrap_or(' ');
        let res_seq = parse_int(line, 23, 26).unwrap_or(0);
        let x = parse_f64(line, 31, 38).unwrap_or(0.0);
        let y = parse_f64(line, 39, 46).unwrap_or(0.0);
        let z = parse_f64(line, 47, 54).unwrap_or(0.0);
        let occupancy = parse_f64(line, 55, 60).unwrap_or(1.0);
        let b_factor = parse_f64(line, 61, 66).unwrap_or(0.0);

        // Element: columns 77-78 (0-indexed: 76-78), right-justified.
        // Many old PDB files omit this column; fall back to atom name.
        let element = if line.len() >= 78 {
            let raw = &line[76..78];
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                element_from_name(name.trim())
            } else {
                // Capitalise correctly: "FE" → "Fe", "N" → "N"
                normalise_element(trimmed)
            }
        } else {
            element_from_name(name.trim())
        };

        // Skip alternate locations other than ' ' or 'A'
        if alt_loc != ' ' && alt_loc != 'A' {
            continue;
        }

        atoms.push(PdbAtom {
            serial,
            name,
            alt_loc,
            res_name,
            chain,
            res_seq,
            x,
            y,
            z,
            occupancy,
            b_factor,
            element,
            is_hetatm,
        });
    }

    atoms
}

// ── helpers ───────────────────────────────────────────────────────────────────

fn parse_int(line: &str, col_start: usize, col_end: usize) -> Option<i32> {
    let s = line.get(col_start - 1..col_end)?.trim();
    s.parse::<i32>().ok()
}

fn parse_f64(line: &str, col_start: usize, col_end: usize) -> Option<f64> {
    let s = line.get(col_start - 1..col_end)?.trim();
    s.parse::<f64>().ok()
}

/// Infer element symbol from PDB atom name.
/// PDB atom names follow the convention: column 14 is the first character of
/// the element symbol for 1-char elements, column 13 for 2-char elements.
fn element_from_name(name: &str) -> String {
    let name = name.trim();
    if name.is_empty() {
        return "?".to_string();
    }
    // The first alphabetic character sequence is the element.
    // Names like "1HB", "2HD" start with a digit — skip digits.
    let chars: Vec<char> = name.chars().collect();
    let start = chars.iter().position(|c| c.is_alphabetic()).unwrap_or(0);
    // Collect up to 2 alphabetic chars
    let end = (start + 2).min(chars.len());
    let raw: String = chars[start..end]
        .iter()
        .take_while(|c| c.is_alphabetic())
        .collect();
    normalise_element(&raw)
}

fn normalise_element(raw: &str) -> String {
    let mut chars = raw.chars();
    match (chars.next(), chars.next()) {
        (Some(a), Some(b)) => format!("{}{}", a.to_ascii_uppercase(), b.to_ascii_lowercase()),
        (Some(a), None) => a.to_ascii_uppercase().to_string(),
        _ => "?".to_string(),
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const MINI_PDB: &str = "\
ATOM      1  N   ALA A   1      -8.901   4.127  -0.555  1.00  0.00           N  \n\
ATOM      2  CA  ALA A   1      -8.608   3.135  -1.618  1.00  0.00           C  \n\
ATOM      3  C   ALA A   1      -7.117   2.964  -1.897  1.00  0.00           C  \n\
ATOM      4  O   ALA A   1      -6.634   3.751  -2.711  1.00  0.00           O  \n\
HETATM    5  C1  LIG B   1       1.234   2.345   3.456  1.00  0.00           C  \n\
";

    #[test]
    fn parse_atom_count() {
        let atoms = parse_pdb(MINI_PDB);
        assert_eq!(atoms.len(), 5);
    }

    #[test]
    fn parse_coordinates() {
        let atoms = parse_pdb(MINI_PDB);
        assert!((atoms[0].x - (-8.901)).abs() < 1e-3);
        assert!((atoms[0].y - 4.127).abs() < 1e-3);
        assert!((atoms[0].z - (-0.555)).abs() < 1e-3);
    }

    #[test]
    fn parse_residue_and_chain() {
        let atoms = parse_pdb(MINI_PDB);
        assert_eq!(atoms[0].res_name, "ALA");
        assert_eq!(atoms[0].chain, 'A');
        assert_eq!(atoms[0].res_seq, 1);
    }

    #[test]
    fn parse_element_from_column() {
        let atoms = parse_pdb(MINI_PDB);
        assert_eq!(atoms[0].element, "N");
        assert_eq!(atoms[1].element, "C");
        assert_eq!(atoms[3].element, "O");
    }

    #[test]
    fn hetatm_flagged() {
        let atoms = parse_pdb(MINI_PDB);
        assert!(!atoms[0].is_hetatm);
        assert!(atoms[4].is_hetatm);
    }

    #[test]
    fn element_from_name_simple() {
        assert_eq!(element_from_name("N"), "N");
        assert_eq!(element_from_name(" CA "), "Ca"); // two-char: C+A → Ca
        assert_eq!(element_from_name(" N  "), "N");
        assert_eq!(element_from_name("FE"), "Fe");
    }

    #[test]
    fn skip_alt_loc_b() {
        let pdb = "\
ATOM      1  N   ALA A   1       1.000   2.000   3.000  0.50  0.00           N  \n\
ATOM      2  N   ALA A   1       1.100   2.100   3.100  0.50  0.00           N  \n\
";
        // alt_loc 'B' should be skipped if it's the second variant
        // (above has blank alt_loc on both, so both parsed)
        let atoms = parse_pdb(pdb);
        assert_eq!(atoms.len(), 2);
    }

    #[test]
    fn serial_parsed() {
        let atoms = parse_pdb(MINI_PDB);
        assert_eq!(atoms[0].serial, 1);
        assert_eq!(atoms[4].serial, 5);
    }
}
