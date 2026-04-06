use crate::Atom;

/// Parse XYZ format: line 1 = count, line 2 = comment, then N lines of "Element x y z"
pub fn parse_xyz(text: &str) -> Vec<Atom> {
    let mut atoms = Vec::new();
    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;
    let mut frame = 1i32;

    while i < lines.len() {
        let count_line = lines[i].trim();
        let num_atoms: usize = match count_line.parse() {
            Ok(n) => n,
            Err(_) => { i += 1; continue; }
        };
        i += 1; // skip count
        if i >= lines.len() { break; }
        i += 1; // skip comment

        for _ in 0..num_atoms {
            if i >= lines.len() { break; }
            let parts: Vec<&str> = lines[i].split_whitespace().collect();
            if parts.len() >= 4 {
                atoms.push(Atom {
                    model_num: frame,
                    chain: String::new(),
                    resid: 0,
                    resname: String::new(),
                    atom_name: parts[0].to_string(),
                    altloc: ' ',
                    x: parts[1].parse().unwrap_or(0.0),
                    y: parts[2].parse().unwrap_or(0.0),
                    z: parts[3].parse().unwrap_or(0.0),
                    occupancy: 1.0,
                    b_factor: 0.0,
                    element: parts[0].to_string(),
                    is_hetatm: false,
                });
            }
            i += 1;
        }
        frame += 1;
    }
    atoms
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_xyz() {
        let xyz = "3\nwater\nO  0.000  0.000  0.117\nH  0.000  0.756 -0.469\nH  0.000 -0.756 -0.469\n";
        let atoms = parse_xyz(xyz);
        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].element, "O");
        assert!((atoms[0].z - 0.117).abs() < 0.001);
    }

    #[test]
    fn test_multi_frame() {
        let xyz = "2\nframe1\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n2\nframe2\nC 0.0 0.0 0.5\nO 1.0 0.0 0.5\n";
        let atoms = parse_xyz(xyz);
        assert_eq!(atoms.len(), 4);
        assert_eq!(atoms[0].model_num, 1);
        assert_eq!(atoms[2].model_num, 2);
    }
}
