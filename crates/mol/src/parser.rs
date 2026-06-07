/// SDF/MOL V2000 parser

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
        "H" => 1.008, "He" => 4.003,
        "C" => 12.011, "N" => 14.007, "O" => 15.999, "F" => 18.998,
        "P" => 30.974, "S" => 32.06, "Cl" => 35.45, "Br" => 79.904,
        "I" => 126.904, "Na" => 22.990, "K" => 39.098, "Ca" => 40.078,
        "Fe" => 55.845, "Si" => 28.086, "Se" => 78.971, "B" => 10.81,
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
            if c > 1 { result.push_str(&c.to_string()); }
            counts.remove("C");
            if let Some(&h) = counts.get("H") {
                result.push('H');
                if h > 1 { result.push_str(&h.to_string()); }
                counts.remove("H");
            }
        }
        for (sym, count) in &counts {
            result.push_str(sym);
            if *count > 1 { result.push_str(&count.to_string()); }
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
        let sum_sq = self.atoms.iter().map(|atom| {
            let dx = atom.x as f64 - cx;
            let dy = atom.y as f64 - cy;
            let dz = atom.z as f64 - cz;
            dx * dx + dy * dy + dz * dz
        }).sum::<f64>();
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
}

fn parse_f32(s: &str) -> f32 {
    s.trim().parse::<f32>().unwrap_or(0.0)
}

fn parse_u32(s: &str) -> u32 {
    s.trim().parse::<u32>().unwrap_or(0)
}

/// Parse a single MOL block (V2000)
pub fn parse_mol(text: &str) -> Option<Molecule> {
    let lines: Vec<&str> = text.lines().collect();
    if lines.len() < 4 {
        return None;
    }

    let name = lines[0].trim().to_string();
    // lines[1]: program/timestamp, lines[2]: comment
    // lines[3]: counts line
    let counts = lines[3];
    if counts.len() < 6 {
        return None;
    }
    let num_atoms = counts[0..3].trim().parse::<usize>().ok()?;
    let num_bonds = counts[3..6].trim().parse::<usize>().ok()?;

    if lines.len() < 4 + num_atoms + num_bonds {
        return None;
    }

    // Atom block
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

    // Bond block
    let mut bonds = Vec::with_capacity(num_bonds);
    for i in 0..num_bonds {
        let line = lines[4 + num_atoms + i];
        if line.len() < 9 {
            continue;
        }
        let a1 = parse_u32(&line[0..3]);
        let a2 = parse_u32(&line[3..6]);
        let bt = line[6..9].trim().parse::<u8>().unwrap_or(1);
        bonds.push(Bond { atom1: a1, atom2: a2, bond_type: bt });
    }

    // Properties (SDF style)
    let mut properties = Vec::new();
    let mut prop_name = String::new();
    let mut in_prop = false;
    for &line in &lines[4 + num_atoms + num_bonds..] {
        if line.starts_with("M  END") {
            continue;
        }
        if line.starts_with("> <") {
            let end = line.find('>').map(|_| {
                line[3..].find('>').map(|p| p + 3).unwrap_or(line.len())
            }).unwrap_or(line.len());
            prop_name = line[3..end].to_string();
            in_prop = true;
        } else if in_prop {
            let val = line.trim();
            if val.is_empty() {
                in_prop = false;
            } else {
                properties.push((prop_name.clone(), val.to_string()));
                in_prop = false;
            }
        }
    }

    Some(Molecule { name, atoms, bonds, properties })
}

/// Parse SDF file (multiple molecules separated by $$$$)
pub fn parse_sdf(text: &str) -> Vec<Molecule> {
    text.split("$$$$")
        .filter(|block| !block.trim().is_empty())
        .filter_map(|block| parse_mol(block.trim()))
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    const SAMPLE_MOL: &str = "ethanol\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0\n    2.3100    1.3300    0.0000 O   0  0  0  0  0\n  1  2  1  0\n  2  3  2  0\nM  END\n";

    #[test]
    fn test_parse_mol() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.bonds.len(), 2);
        assert_eq!(mol.atoms[0].symbol, "C");
        assert_eq!(mol.atoms[2].symbol, "O");
    }

    #[test]
    fn test_formula() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        assert_eq!(mol.formula(), "C2O");
    }

    #[test]
    fn test_weight() {
        let mol = parse_mol(SAMPLE_MOL).unwrap();
        let w = mol.weight();
        assert!((w - (12.011 * 2.0 + 15.999)).abs() < 0.01);
    }

    #[test]
    fn test_sdf() {
        let sdf = format!("{}$$$$\n{}$$$$\n", SAMPLE_MOL, SAMPLE_MOL);
        let mols = parse_sdf(&sdf);
        assert_eq!(mols.len(), 2);
    }
}
