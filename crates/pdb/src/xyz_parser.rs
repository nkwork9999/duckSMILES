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
            Err(_) => {
                i += 1;
                continue;
            }
        };
        i += 1; // skip count
        if i >= lines.len() {
            break;
        }
        i += 1; // skip comment

        for _ in 0..num_atoms {
            if i >= lines.len() {
                break;
            }
            let parts: Vec<&str> = lines[i].split_whitespace().collect();
            if parts.len() >= 4 {
                let coordinates: Option<Vec<f32>> = parts[1..4]
                    .iter()
                    .map(|part| part.parse::<f32>().ok().filter(|v| v.is_finite()))
                    .collect();
                let Some(coordinates) = coordinates else {
                    i += 1;
                    continue;
                };
                atoms.push(Atom {
                    model_num: frame,
                    chain: String::new(),
                    resid: 0,
                    resname: String::new(),
                    atom_name: parts[0].to_string(),
                    altloc: ' ',
                    x: coordinates[0],
                    y: coordinates[1],
                    z: coordinates[2],
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
    fn invalid_coordinates_are_not_silently_zeroed() {
        for bad in ["invalid", "NaN", "inf", "-inf", "1e999"] {
            for axis in 0..3 {
                let mut coordinates = ["1", "2", "3"];
                coordinates[axis] = bad;
                let xyz = format!("2\ncomment\nC {}\nO 4 5 6\n", coordinates.join(" "));
                let atoms = parse_xyz(&xyz);
                assert_eq!(atoms.len(), 1, "{bad} at axis {axis}");
                assert_eq!(atoms[0].element, "O");
                assert_eq!(atoms[0].x, 4.0);
            }
        }
    }

    #[test]
    fn invalid_atom_keeps_frame_numbering() {
        let atoms = parse_xyz("1\nbad\nC bad 0 0\n1\ngood\nO 1 2 3\n");
        assert_eq!(atoms.len(), 1);
        assert_eq!(atoms[0].model_num, 2);
    }

    #[test]
    fn test_parse_xyz() {
        let xyz =
            "3\nwater\nO  0.000  0.000  0.117\nH  0.000  0.756 -0.469\nH  0.000 -0.756 -0.469\n";
        let atoms = parse_xyz(xyz);
        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].element, "O");
        assert!((atoms[0].z - 0.117).abs() < 0.001);
    }

    #[test]
    fn test_multi_frame() {
        let xyz =
            "2\nframe1\nC 0.0 0.0 0.0\nO 1.0 0.0 0.0\n2\nframe2\nC 0.0 0.0 0.5\nO 1.0 0.0 0.5\n";
        let atoms = parse_xyz(xyz);
        assert_eq!(atoms.len(), 4);
        assert_eq!(atoms[0].model_num, 1);
        assert_eq!(atoms[2].model_num, 2);
    }
}
