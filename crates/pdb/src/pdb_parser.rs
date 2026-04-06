use crate::Atom;

/// Parse PDB fixed-column format
pub fn parse_pdb(text: &str, include_hetatm: bool) -> Vec<Atom> {
    let mut atoms = Vec::new();
    let mut model_num = 1i32;

    for line in text.lines() {
        if line.starts_with("MODEL") {
            model_num = line[6..].trim().parse().unwrap_or(model_num);
        } else if line.starts_with("ATOM") || (include_hetatm && line.starts_with("HETATM")) {
            if line.len() < 54 {
                continue;
            }
            let is_hetatm = line.starts_with("HETATM");
            let atom_name = line.get(12..16).unwrap_or("    ").trim().to_string();
            let altloc = line.as_bytes().get(16).map(|&b| b as char).unwrap_or(' ');
            let resname = line.get(17..20).unwrap_or("   ").trim().to_string();
            let chain = line.get(21..22).unwrap_or(" ").trim().to_string();
            let resid = line.get(22..26).unwrap_or("    ").trim().parse().unwrap_or(0);
            let x = line.get(30..38).unwrap_or("").trim().parse().unwrap_or(0.0);
            let y = line.get(38..46).unwrap_or("").trim().parse().unwrap_or(0.0);
            let z = line.get(46..54).unwrap_or("").trim().parse().unwrap_or(0.0);
            let occupancy = line.get(54..60).unwrap_or("").trim().parse().unwrap_or(1.0);
            let b_factor = line.get(60..66).unwrap_or("").trim().parse().unwrap_or(0.0);
            let element = line.get(76..78).unwrap_or("").trim().to_string();

            atoms.push(Atom {
                model_num, chain, resid, resname, atom_name, altloc,
                x, y, z, occupancy, b_factor, element, is_hetatm,
            });
        }
    }
    atoms
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_atom_line() {
        let pdb = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000  1.00 20.00           C  \n";
        let atoms = parse_pdb(pdb, false);
        assert_eq!(atoms.len(), 1);
        assert_eq!(atoms[0].atom_name, "CA");
        assert_eq!(atoms[0].resname, "ALA");
        assert_eq!(atoms[0].chain, "A");
        assert_eq!(atoms[0].resid, 1);
        assert!((atoms[0].x - 1.0).abs() < 0.01);
        assert!((atoms[0].y - 2.0).abs() < 0.01);
        assert!((atoms[0].z - 3.0).abs() < 0.01);
        assert_eq!(atoms[0].element, "C");
    }

    #[test]
    fn test_hetatm_filter() {
        let pdb = "ATOM      1  CA  ALA A   1       0.0     0.0     0.0   1.00  0.00           C  \nHETATM    2  O   HOH A 100       1.0     1.0     1.0   1.00  0.00           O  \n";
        assert_eq!(parse_pdb(pdb, false).len(), 1);
        assert_eq!(parse_pdb(pdb, true).len(), 2);
    }
}
