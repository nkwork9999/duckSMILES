/// SDF/MOL V2000/V3000 parser

#[derive(Debug, Clone)]
pub struct Atom {
    pub x: f32,
    pub y: f32,
    pub z: f32,
    pub symbol: String,
}

#[derive(Debug, Clone)]
pub struct Bond {
    pub atom1: u32,
    pub atom2: u32,
    pub bond_type: u8,
}

#[derive(Debug, Clone)]
pub struct Molecule {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub bonds: Vec<Bond>,
    pub properties: Vec<(String, String)>,
}

// Atomic weights for common elements
fn atomic_weight(sym: &str) -> f64 {
    match sym {
        "H" => 1.008,
        "He" => 4.003,
        "C" => 12.011,
        "N" => 14.007,
        "O" => 15.999,
        "F" => 18.998,
        "P" => 30.974,
        "S" => 32.06,
        "Cl" => 35.45,
        "Br" => 79.904,
        "I" => 126.904,
        "Na" => 22.990,
        "K" => 39.098,
        "Ca" => 40.078,
        "Fe" => 55.845,
        "Si" => 28.086,
        "Se" => 78.971,
        "B" => 10.81,
        _ => 0.0,
    }
}

impl Molecule {
    pub fn formula(&self) -> String {
        let mut counts = std::collections::BTreeMap::new();
        for a in &self.atoms {
            *counts.entry(a.symbol.clone()).or_insert(0u32) += 1;
        }
        // Hill system: C first, then H, then alphabetical
        let mut result = String::new();
        if let Some(&c) = counts.get("C") {
            result.push('C');
            if c > 1 {
                result.push_str(&c.to_string());
            }
            counts.remove("C");
            if let Some(&h) = counts.get("H") {
                result.push('H');
                if h > 1 {
                    result.push_str(&h.to_string());
                }
                counts.remove("H");
            }
        }
        for (sym, count) in &counts {
            result.push_str(sym);
            if *count > 1 {
                result.push_str(&count.to_string());
            }
        }
        result
    }

    pub fn weight(&self) -> f64 {
        self.atoms.iter().map(|a| atomic_weight(&a.symbol)).sum()
    }

    pub fn has_3d(&self) -> bool {
        self.atoms.iter().any(|a| a.z.abs() > 1e-4)
    }

    pub fn centroid(&self) -> Option<(f64, f64, f64)> {
        if self.atoms.is_empty() {
            return None;
        }
        let n = self.atoms.len() as f64;
        let (sx, sy, sz) = self.atoms.iter().fold((0.0, 0.0, 0.0), |acc, atom| {
            (
                acc.0 + atom.x as f64,
                acc.1 + atom.y as f64,
                acc.2 + atom.z as f64,
            )
        });
        Some((sx / n, sy / n, sz / n))
    }

    pub fn radius_of_gyration(&self) -> Option<f64> {
        let (cx, cy, cz) = self.centroid()?;
        let n = self.atoms.len() as f64;
        let sum_sq = self
            .atoms
            .iter()
            .map(|atom| {
                let dx = atom.x as f64 - cx;
                let dy = atom.y as f64 - cy;
                let dz = atom.z as f64 - cz;
                dx * dx + dy * dy + dz * dz
            })
            .sum::<f64>();
        Some((sum_sq / n).sqrt())
    }

    pub fn coordinate_bounds(&self) -> Option<[(f64, f64); 3]> {
        let first = self.atoms.first()?;
        let mut min_x = first.x as f64;
        let mut min_y = first.y as f64;
        let mut min_z = first.z as f64;
        let mut max_x = min_x;
        let mut max_y = min_y;
        let mut max_z = min_z;
        for atom in self.atoms.iter().skip(1) {
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

    pub fn property(&self, key: &str) -> Option<&str> {
        self.properties
            .iter()
            .find(|(name, _)| name == key)
            .map(|(_, value)| value.as_str())
    }

    pub fn properties_json(&self) -> String {
        properties_to_json(&self.properties)
    }

    pub fn atoms_json(&self) -> String {
        atoms_to_json(&self.atoms)
    }

    pub fn bonds_json(&self) -> String {
        bonds_to_json(&self.bonds)
    }

    pub fn mol_json(&self) -> String {
        let mut out = String::from("{\"name\":\"");
        out.push_str(&json_escape(&self.name));
        out.push_str("\",\"formula\":\"");
        out.push_str(&json_escape(&self.formula()));
        out.push_str("\",\"weight\":");
        push_json_f64(&mut out, self.weight());
        out.push_str(",\"num_atoms\":");
        out.push_str(&self.atoms.len().to_string());
        out.push_str(",\"num_bonds\":");
        out.push_str(&self.bonds.len().to_string());
        out.push_str(",\"has_3d\":");
        out.push_str(if self.has_3d() { "true" } else { "false" });
        out.push_str(",\"atoms\":");
        out.push_str(&self.atoms_json());
        out.push_str(",\"bonds\":");
        out.push_str(&self.bonds_json());
        out.push_str(",\"properties\":");
        out.push_str(&self.properties_json());
        out.push('}');
        out
    }
}

fn parse_f32(s: &str) -> f32 {
    s.trim().parse::<f32>().unwrap_or(0.0)
}

fn parse_u32(s: &str) -> u32 {
    s.trim().parse::<u32>().unwrap_or(0)
}

fn parse_property_name(line: &str) -> Option<String> {
    if !line.starts_with('>') {
        return None;
    }
    let start = line.find('<')?;
    let rest = &line[start + 1..];
    let end = rest.find('>')?;
    Some(rest[..end].to_string())
}

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for ch in s.chars() {
        match ch {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            '\u{08}' => out.push_str("\\b"),
            '\u{0c}' => out.push_str("\\f"),
            ch if ch <= '\u{1f}' => {
                out.push_str(&format!("\\u{:04x}", ch as u32));
            }
            _ => out.push(ch),
        }
    }
    out
}

fn push_json_f32(out: &mut String, value: f32) {
    if value == 0.0 {
        out.push('0');
    } else {
        out.push_str(&value.to_string());
    }
}

fn push_json_f64(out: &mut String, value: f64) {
    if value == 0.0 {
        out.push('0');
    } else {
        out.push_str(&value.to_string());
    }
}

pub fn atoms_to_json(atoms: &[Atom]) -> String {
    let mut out = String::from("[");
    for (i, atom) in atoms.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"index\":");
        out.push_str(&(i + 1).to_string());
        out.push_str(",\"symbol\":\"");
        out.push_str(&json_escape(&atom.symbol));
        out.push_str("\",\"x\":");
        push_json_f32(&mut out, atom.x);
        out.push_str(",\"y\":");
        push_json_f32(&mut out, atom.y);
        out.push_str(",\"z\":");
        push_json_f32(&mut out, atom.z);
        out.push('}');
    }
    out.push(']');
    out
}

pub fn bonds_to_json(bonds: &[Bond]) -> String {
    let mut out = String::from("[");
    for (i, bond) in bonds.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"index\":");
        out.push_str(&(i + 1).to_string());
        out.push_str(",\"atom1\":");
        out.push_str(&bond.atom1.to_string());
        out.push_str(",\"atom2\":");
        out.push_str(&bond.atom2.to_string());
        out.push_str(",\"bond_type\":");
        out.push_str(&bond.bond_type.to_string());
        out.push('}');
    }
    out.push(']');
    out
}

pub fn properties_to_json(properties: &[(String, String)]) -> String {
    let mut out = String::from("[");
    for (i, (name, value)) in properties.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"name\":\"");
        out.push_str(&json_escape(name));
        out.push_str("\",\"value\":\"");
        out.push_str(&json_escape(value));
        out.push_str("\"}");
    }
    out.push(']');
    out
}

pub fn sdf_properties_json(mols: &[Molecule]) -> String {
    let mut out = String::from("[");
    for (i, mol) in mols.iter().enumerate() {
        if i > 0 {
            out.push(',');
        }
        out.push_str("{\"record\":");
        out.push_str(&(i + 1).to_string());
        out.push_str(",\"name\":\"");
        out.push_str(&json_escape(&mol.name));
        out.push_str("\",\"properties\":");
        out.push_str(&mol.properties_json());
        out.push('}');
    }
    out.push(']');
    out
}

fn parse_sdf_properties(lines: &[&str], start: usize) -> Vec<(String, String)> {
    let mut properties = Vec::new();
    let mut i = start;
    while i < lines.len() {
        let line = lines[i];
        if line.starts_with("$$$$") {
            break;
        }
        let Some(prop_name) = parse_property_name(line) else {
            i += 1;
            continue;
        };

        i += 1;
        let mut value_lines = Vec::new();
        while i < lines.len() {
            let value_line = lines[i];
            if value_line.starts_with("$$$$") || value_line.trim().is_empty() {
                break;
            }
            value_lines.push(value_line);
            i += 1;
        }
        properties.push((prop_name, value_lines.join("\n")));

        while i < lines.len() && lines[i].trim().is_empty() {
            i += 1;
        }
    }
    properties
}

fn parse_mol_v2000(lines: &[&str]) -> Option<Molecule> {
    let name = lines[0].trim().to_string();
    let counts = lines[3];
    if counts.len() < 6 {
        return None;
    }
    let num_atoms = counts[0..3].trim().parse::<usize>().ok()?;
    let num_bonds = counts[3..6].trim().parse::<usize>().ok()?;

    if lines.len() < 4 + num_atoms + num_bonds {
        return None;
    }

    let mut atoms = Vec::with_capacity(num_atoms);
    for i in 0..num_atoms {
        let line = lines[4 + i];
        if line.len() < 34 {
            continue;
        }
        let x = parse_f32(&line[0..10]);
        let y = parse_f32(&line[10..20]);
        let z = parse_f32(&line[20..30]);
        let symbol = line[30..34].trim().to_string();
        atoms.push(Atom { x, y, z, symbol });
    }

    let mut bonds = Vec::with_capacity(num_bonds);
    for i in 0..num_bonds {
        let line = lines[4 + num_atoms + i];
        if line.len() < 9 {
            continue;
        }
        let a1 = parse_u32(&line[0..3]);
        let a2 = parse_u32(&line[3..6]);
        let bt = line[6..9].trim().parse::<u8>().unwrap_or(1);
        bonds.push(Bond {
            atom1: a1,
            atom2: a2,
            bond_type: bt,
        });
    }

    let properties = parse_sdf_properties(lines, 4 + num_atoms + num_bonds);

    Some(Molecule {
        name,
        atoms,
        bonds,
        properties,
    })
}

fn v3000_records(lines: &[&str]) -> Vec<String> {
    let mut records = Vec::new();
    let mut current = String::new();
    for line in lines {
        if !line.starts_with("M  V30") {
            continue;
        }
        let mut payload = line.get(6..).unwrap_or("").trim_start().to_string();
        let continued = payload.trim_end().ends_with('-');
        if continued {
            payload = payload.trim_end_matches('-').trim_end().to_string();
        }

        if !current.is_empty() && !payload.is_empty() {
            current.push(' ');
        }
        current.push_str(&payload);

        if !continued {
            records.push(current.trim().to_string());
            current.clear();
        }
    }
    if !current.trim().is_empty() {
        records.push(current.trim().to_string());
    }
    records
}

fn parse_mol_v3000(lines: &[&str]) -> Option<Molecule> {
    let name = lines[0].trim().to_string();
    let records = v3000_records(lines);
    let counts = records
        .iter()
        .find(|record| record.starts_with("COUNTS "))?;
    let counts_parts: Vec<&str> = counts.split_whitespace().collect();
    if counts_parts.len() < 3 {
        return None;
    }
    let num_atoms = counts_parts[1].parse::<usize>().ok()?;
    let num_bonds = counts_parts[2].parse::<usize>().ok()?;

    let mut atoms = Vec::with_capacity(num_atoms);
    let mut bonds = Vec::with_capacity(num_bonds);
    let mut in_atom_block = false;
    let mut in_bond_block = false;

    for record in &records {
        match record.as_str() {
            "BEGIN ATOM" => {
                in_atom_block = true;
                in_bond_block = false;
                continue;
            }
            "END ATOM" => {
                in_atom_block = false;
                continue;
            }
            "BEGIN BOND" => {
                in_bond_block = true;
                in_atom_block = false;
                continue;
            }
            "END BOND" => {
                in_bond_block = false;
                continue;
            }
            _ => {}
        }

        if in_atom_block {
            let parts: Vec<&str> = record.split_whitespace().collect();
            if parts.len() < 6 {
                continue;
            }
            let symbol = parts[1].to_string();
            let x = parse_f32(parts[2]);
            let y = parse_f32(parts[3]);
            let z = parse_f32(parts[4]);
            atoms.push(Atom { x, y, z, symbol });
        } else if in_bond_block {
            let parts: Vec<&str> = record.split_whitespace().collect();
            if parts.len() < 4 {
                continue;
            }
            let bt = parts[1].parse::<u8>().unwrap_or(1);
            let a1 = parse_u32(parts[2]);
            let a2 = parse_u32(parts[3]);
            bonds.push(Bond {
                atom1: a1,
                atom2: a2,
                bond_type: bt,
            });
        }
    }

    if atoms.len() != num_atoms || bonds.len() != num_bonds {
        return None;
    }

    let property_start = lines
        .iter()
        .position(|line| line.trim() == "M  END")
        .map(|idx| idx + 1)
        .unwrap_or(lines.len());
    let properties = parse_sdf_properties(lines, property_start);

    Some(Molecule {
        name,
        atoms,
        bonds,
        properties,
    })
}

/// Parse a single MOL block (V2000 or V3000)
pub fn parse_mol(text: &str) -> Option<Molecule> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 4 {
        return None;
    }
    if lines[3].contains("V3000") || lines.iter().any(|line| line.starts_with("M  V30")) {
        parse_mol_v3000(&lines)
    } else {
        parse_mol_v2000(&lines)
    }
}

/// Parse SDF file (multiple molecules separated by $$$$)
pub fn parse_sdf(text: &str) -> Vec<Molecule> {
    text.split("$$$$")
        .map(|block| block.trim_matches(|c| c == '\n' || c == '\r'))
        .filter(|block| !block.trim().is_empty())
        .filter_map(parse_mol)
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_MOL: &str = "ethanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0\n    2.3100    1.3300    0.0000 O   0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\nM  END\n";
    const SAMPLE_MOL_V3000: &str = "ethanol_v3000\n  test\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.0000 0.0000 0.0000 0\nM  V30 2 C 1.5400 0.0000 0.0000 0 CFG=0\nM  V30 3 O 2.3100 1.3300 1.0000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2 CFG=0\nM  V30 2 2 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";
    const SAMPLE_MOL_V3000_CONTINUED: &str = "continued_v3000\n  test\n\n  0  0  0     0  0            999 V3000\nM  V30 BEGIN CTAB\nM  V30 COUNTS 3 2 0 0 0\nM  V30 BEGIN ATOM\nM  V30 1 C 0.0000 0.0000 0.0000 0\nM  V30 2 C 1.5400 0.0000 0.0000 0 ATTCHPT=(2 1 -\nM  V30 2)\nM  V30 3 O 2.3100 1.3300 1.0000 0\nM  V30 END ATOM\nM  V30 BEGIN BOND\nM  V30 1 1 1 2\nM  V30 2 2 2 3\nM  V30 END BOND\nM  V30 END CTAB\nM  END\n";

    #[test]
    fn test_parse_mol() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.atoms[0].symbol, "C");
        assert_eq!(mol.atoms[2].symbol, "O");
    }

    #[test]
    fn test_parse_mol_v3000() {
        let mol = parse_mol(SAMPLE_MOL_V3000).unwrap();
        assert_eq!(mol.name, "ethanol_v3000");
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.atoms[0].symbol, "C");
        assert_eq!(mol.atoms[2].symbol, "O");
        assert_eq!(mol.bonds[0].atom1, 1);
        assert_eq!(mol.bonds[0].atom2, 2);
        assert_eq!(mol.bonds[0].bond_type, 1);
        assert_eq!(mol.bonds[1].bond_type, 2);
        assert!(mol.has_3d());
        assert!((mol.centroid().unwrap().2 - 0.333333).abs() < 0.01);
    }

    #[test]
    fn test_parse_mol_v3000_line_continuation() {
        let mol = parse_mol(SAMPLE_MOL_V3000_CONTINUED).unwrap();
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.formula(), "C2O");
    }

    #[test]
    fn test_formula() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        assert_eq!(mol.formula(), "C2O");
        let mol_v3000 = parse_mol(SAMPLE_MOL_V3000).unwrap();
        assert_eq!(mol_v3000.formula(), "C2O");
    }

    #[test]
    fn test_weight() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        let w = mol.weight();
        assert!((w - (12.011 * 2.0 + 15.999)).abs() < 0.01);
        let mol_v3000 = parse_mol(SAMPLE_MOL_V3000).unwrap();
        let w_v3000 = mol_v3000.weight();
        assert!((w_v3000 - (12.011 * 2.0 + 15.999)).abs() < 0.01);
    }

    #[test]
    fn test_atoms_bonds_and_mol_json() {
        let mol = parse_mol(SAMPLE_MOL_V3000).unwrap();
        assert_eq!(
            mol.atoms_json(),
            "[{\"index\":1,\"symbol\":\"C\",\"x\":0,\"y\":0,\"z\":0},{\"index\":2,\"symbol\":\"C\",\"x\":1.54,\"y\":0,\"z\":0},{\"index\":3,\"symbol\":\"O\",\"x\":2.31,\"y\":1.33,\"z\":1}]"
        );
        assert_eq!(
            mol.bonds_json(),
            "[{\"index\":1,\"atom1\":1,\"atom2\":2,\"bond_type\":1},{\"index\":2,\"atom1\":2,\"atom2\":3,\"bond_type\":2}]"
        );

        let mol_with_props = parse_sdf(&format!(
            "{}> <ID>\nV3000-1\n\n> <NOTE>\nline one\nline two\n\n$$$$\n",
            SAMPLE_MOL_V3000
        ))
        .pop()
        .unwrap();
        assert_eq!(
            mol_with_props.mol_json(),
            "{\"name\":\"ethanol_v3000\",\"formula\":\"C2O\",\"weight\":40.021,\"num_atoms\":3,\"num_bonds\":2,\"has_3d\":true,\"atoms\":[{\"index\":1,\"symbol\":\"C\",\"x\":0,\"y\":0,\"z\":0},{\"index\":2,\"symbol\":\"C\",\"x\":1.54,\"y\":0,\"z\":0},{\"index\":3,\"symbol\":\"O\",\"x\":2.31,\"y\":1.33,\"z\":1}],\"bonds\":[{\"index\":1,\"atom1\":1,\"atom2\":2,\"bond_type\":1},{\"index\":2,\"atom1\":2,\"atom2\":3,\"bond_type\":2}],\"properties\":[{\"name\":\"ID\",\"value\":\"V3000-1\"},{\"name\":\"NOTE\",\"value\":\"line one\\nline two\"}]}"
        );
    }

    #[test]
    fn test_sdf() {
        let sdf = format!("{}$$$$\n{}$$$$\n", SAMPLE_MOL, SAMPLE_MOL);
        let mols = parse_sdf(&sdf);
        assert_eq!(mols.len(), 2);
        let sdf_v3000 = format!("{}$$$$\n{}$$$$\n", SAMPLE_MOL, SAMPLE_MOL_V3000);
        let mols_v3000 = parse_sdf(&sdf_v3000);
        assert_eq!(mols_v3000.len(), 2);
        assert_eq!(mols_v3000[1].name, "ethanol_v3000");
        assert_eq!(mols_v3000[1].formula(), "C2O");
    }

    #[test]
    fn test_sdf_v3000_properties() {
        let sdf = format!(
            "{}> <ID>\nV3000-1\n\n> <NOTE>\nline one\nline two\n\n$$$$\n",
            SAMPLE_MOL_V3000
        );
        let mol = parse_sdf(&sdf).pop().unwrap();
        assert_eq!(mol.property("ID"), Some("V3000-1"));
        assert_eq!(mol.property("NOTE"), Some("line one\nline two"));
        assert_eq!(
            mol.properties_json(),
            "[{\"name\":\"ID\",\"value\":\"V3000-1\"},{\"name\":\"NOTE\",\"value\":\"line one\\nline two\"}]"
        );
    }

    #[test]
    fn test_sdf_properties_preserve_multiline_duplicates_and_empty_values() {
        let sdf = format!(
            "{}> <ID>\n123\n\n> <NOTE>\nline one\nline two\n\n> <ID>\n456\n\n> <EMPTY>\n\n$$$$\n",
            SAMPLE_MOL
        );
        let mol = parse_sdf(&sdf).pop().unwrap();
        assert_eq!(mol.property("ID"), Some("123"));
        assert_eq!(mol.property("NOTE"), Some("line one\nline two"));
        assert_eq!(mol.property("EMPTY"), Some(""));
        assert_eq!(mol.properties.len(), 4);
        assert_eq!(
            mol.properties_json(),
            "[{\"name\":\"ID\",\"value\":\"123\"},{\"name\":\"NOTE\",\"value\":\"line one\\nline two\"},{\"name\":\"ID\",\"value\":\"456\"},{\"name\":\"EMPTY\",\"value\":\"\"}]"
        );
    }

    #[test]
    fn test_sdf_properties_json_includes_record_names() {
        let sdf = format!(
            "{}> <ID>\n1\n\n$$$$\n{}> <ID>\n2\n\n$$$$\n",
            SAMPLE_MOL, SAMPLE_MOL
        );
        let mols = parse_sdf(&sdf);
        assert_eq!(
            sdf_properties_json(&mols),
            "[{\"record\":1,\"name\":\"ethanol\",\"properties\":[{\"name\":\"ID\",\"value\":\"1\"}]},{\"record\":2,\"name\":\"ethanol\",\"properties\":[{\"name\":\"ID\",\"value\":\"2\"}]}]"
        );
    }
}
