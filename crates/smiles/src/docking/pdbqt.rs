use crate::conformer::generate_conformer;
use crate::docking::atomtype::{assign_types, pdb_atom_type, VinaType};
use crate::docking::pdb::parse_pdb;
use crate::parser::parse;

// ── PDBQT format ─────────────────────────────────────────────────────────────
// PDBQT extends PDB with two extra columns at the end of each ATOM line:
//   col 71-76: partial charge (6.3f, prefixed with '+' or '-')
//   col 78-79: AutoDock/Vina atom type (e.g. "C", "A", "OA", "HD")
//
// Example:
//   ATOM      1  C   LIG A   1       1.234   2.345   3.456  1.00  0.00    +0.000 C
//   ATOM      2  C   LIG A   1       2.345   3.456   4.567  1.00  0.00    +0.000 A
//   ATOM      3  OA  LIG A   1       3.456   4.567   5.678  1.00  0.00    -0.300 OA
//
// For ligands: ROOT / ENDROOT / TORSDOF markers added.
// For proteins: plain ATOM records only.

// ── PDBQT line formatter ──────────────────────────────────────────────────────

fn pdbqt_atom_line(
    serial: u32,
    atom_name: &str,
    res_name: &str,
    chain: char,
    res_seq: i32,
    x: f64,
    y: f64,
    z: f64,
    charge: f64,
    vtype: VinaType,
    is_hetatm: bool,
) -> String {
    let record = if is_hetatm { "HETATM" } else { "ATOM  " };
    let type_str = vtype.as_str();
    let charge_str = if charge >= 0.0 {
        format!("+{:.3}", charge)
    } else {
        format!("{:.3}", charge)
    };
    // PDBQT: cols 1-6 record, 7-11 serial, 13-16 name, 18-20 resname,
    //        22 chain, 23-26 resseq, 31-38 x, 39-46 y, 47-54 z,
    //        55-60 occ, 61-66 bfac, 71-76 charge, 78-79 type
    format!(
        "{record}{serial:>5}  {name:<4}{resname:<4}{chain}{resseq:>4}    {x:>8.3}{y:>8.3}{z:>8.3}  1.00  0.00    {charge:>6} {vtype}\n",
        record = record,
        serial = serial,
        name = format!("{:<4}", atom_name),
        resname = format!("{:<3}", res_name),
        chain = chain,
        resseq = res_seq,
        x = x,
        y = y,
        z = z,
        charge = charge_str,
        vtype = type_str,
    )
}

// ── Ligand PDBQT from SMILES ──────────────────────────────────────────────────
// Generates a rigid PDBQT for a small molecule from its SMILES string.
// Coordinates come from Phase 1 (distance geometry + UFF-lite minimisation).
// Charges are set to 0.0 (Vina does not use partial charges).
// Returns None if the SMILES is invalid or embedding fails.

pub fn smiles_to_pdbqt(smiles: &str, seed: u64) -> Option<String> {
    let mol_bare = parse(smiles)?;
    let mol = mol_bare.with_explicit_hydrogens();
    let coords = generate_conformer(smiles, seed)?;
    let types = assign_types(&mol);

    let n = mol.atoms.len();
    let mut out = String::new();

    out.push_str(&format!("REMARK  SMILES = {}\n", smiles));
    out.push_str("ROOT\n");

    for i in 0..n {
        let atom = &mol.atoms[i];
        // Skip non-polar hydrogens for Vina (optional, but standard practice)
        // We keep them here so the output is complete; Vina handles both.
        let x = coords[3 * i];
        let y = coords[3 * i + 1];
        let z = coords[3 * i + 2];
        let type_str = types[i].as_str();
        // atom_name: element + serial within element (simple scheme)
        let atom_name = format!("{}{}", atom.symbol, i + 1);
        let line = pdbqt_atom_line(
            (i + 1) as u32,
            &atom_name,
            "LIG",
            'A',
            1,
            x, y, z,
            0.0,
            types[i],
            false,
        );
        let _ = type_str;
        out.push_str(&line);
    }

    out.push_str("ENDROOT\n");
    out.push_str(&format!("TORSDOF 0\n"));

    Some(out)
}

// ── Protein PDBQT from PDB text ───────────────────────────────────────────────
// Parses a PDB file and writes a PDBQT with Vina atom types.
// Hydrogens from the PDB are included as HD/H.
// If no H are present (standard protein PDB), Vina will use polar protein H
// from its internal knowledge — this is fine for rigid docking.

pub fn pdb_to_pdbqt(pdb_text: &str) -> String {
    let atoms = parse_pdb(pdb_text);
    let mut out = String::new();

    for atom in &atoms {
        let elem = atom.element.as_str();
        let vtype = pdb_atom_type(elem, &atom.name, &atom.res_name, atom.is_hetatm);
        if vtype == VinaType::Unknown {
            continue; // skip atoms we cannot type (metals not in our table, etc.)
        }
        let line = pdbqt_atom_line(
            atom.serial,
            atom.name.trim(),
            &atom.res_name,
            atom.chain,
            atom.res_seq,
            atom.x, atom.y, atom.z,
            0.0,
            vtype,
            atom.is_hetatm,
        );
        out.push_str(&line);
    }

    out
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn smiles_to_pdbqt_methane() {
        let pdbqt = smiles_to_pdbqt("C", 0).unwrap();
        assert!(pdbqt.contains("ROOT"), "missing ROOT");
        assert!(pdbqt.contains("ENDROOT"), "missing ENDROOT");
        assert!(pdbqt.contains("TORSDOF"), "missing TORSDOF");
        assert!(pdbqt.contains("LIG"), "missing residue name");
        // methane: 1 C + 4 H = 5 atoms
        let atom_lines = pdbqt.lines().filter(|l| l.starts_with("ATOM")).count();
        assert_eq!(atom_lines, 5, "wrong atom count: {atom_lines}");
    }

    #[test]
    fn smiles_to_pdbqt_benzene_has_aromatic_type() {
        let pdbqt = smiles_to_pdbqt("c1ccccc1", 0).unwrap();
        // Should contain type "A" for aromatic carbons
        assert!(pdbqt.contains(" A\n") || pdbqt.contains(" A\r"), "no aromatic type A in\n{pdbqt}");
    }

    #[test]
    fn smiles_to_pdbqt_invalid_smiles() {
        assert!(smiles_to_pdbqt("not-a-smiles", 0).is_none());
    }

    #[test]
    fn pdb_to_pdbqt_basic() {
        let pdb = "\
ATOM      1  N   ALA A   1      -8.901   4.127  -0.555  1.00  0.00           N  \n\
ATOM      2  CA  ALA A   1      -8.608   3.135  -1.618  1.00  0.00           C  \n\
ATOM      3  CG  PHE B   2       1.000   2.000   3.000  1.00  0.00           C  \n\
";
        let pdbqt = pdb_to_pdbqt(pdb);
        // N backbone → NS
        assert!(pdbqt.contains(" NS\n") || pdbqt.contains(" NS\r"), "no NS in\n{pdbqt}");
        // ALA CA → C
        assert!(pdbqt.contains(" C\n") || pdbqt.contains(" C\r"), "no C in\n{pdbqt}");
        // PHE CG → A
        assert!(pdbqt.contains(" A\n") || pdbqt.contains(" A\r"), "no A in\n{pdbqt}");
    }

    #[test]
    fn pdbqt_line_coord_format() {
        let pdbqt = smiles_to_pdbqt("C", 0).unwrap();
        // Every ATOM line should be parseable with x,y,z in cols 31-54
        for line in pdbqt.lines() {
            if !line.starts_with("ATOM") { continue; }
            let x: f64 = line[30..38].trim().parse().expect("x not parseable");
            let y: f64 = line[38..46].trim().parse().expect("y not parseable");
            let z: f64 = line[46..54].trim().parse().expect("z not parseable");
            assert!(x.is_finite() && y.is_finite() && z.is_finite());
        }
    }

    #[test]
    fn pdbqt_charge_format() {
        let pdbqt = smiles_to_pdbqt("CO", 0).unwrap();
        // Every ATOM line should contain a charge in "+N.NNN" format
        for line in pdbqt.lines() {
            if !line.starts_with("ATOM") { continue; }
            // charge is at position 70..76 (0-indexed)
            if line.len() >= 76 {
                let charge_str = &line[69..75];
                let charge: f64 = charge_str.trim().parse()
                    .expect(&format!("charge not parseable: '{charge_str}' in '{line}'"));
                assert!((charge - 0.0).abs() < 1e-6);
            }
        }
    }
}
