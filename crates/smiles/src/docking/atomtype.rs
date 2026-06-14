use crate::parser::Molecule;
use crate::conformer::params::Hybridization;

// ── Vina atom types ───────────────────────────────────────────────────────────
// Used by AutoDock Vina's scoring function.  Each type maps to a set of
// interaction parameters (steric, H-bond, hydrophobic, …).
// Reference: Trott & Olson 2010, J. Comput. Chem., Vina source atom_type.h

#[allow(dead_code)]
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum VinaType {
    C,   // aliphatic carbon
    A,   // aromatic carbon
    N,   // generic nitrogen (no lone pair available / no H)
    NA,  // H-bond acceptor nitrogen (lone pair, e.g. pyridine-N)
    NS,  // H-bond donor nitrogen (N–H)
    O,   // generic oxygen
    OA,  // H-bond acceptor oxygen (C=O, ether, …)
    S,   // generic sulfur
    SA,  // H-bond acceptor sulfur (thioether, thiol)
    P,   // phosphorus
    F,   // fluorine
    CL,  // chlorine
    BR,  // bromine
    I,   // iodine
    H,   // non-polar hydrogen (bonded to C)
    HD,  // polar hydrogen / H-bond donor (bonded to N or O)
    MG,  // magnesium
    ZN,  // zinc
    FE,  // iron
    CA,  // calcium
    MN,  // manganese
    CU,  // copper
    Unknown,
}

impl VinaType {
    /// Vina PDBQT type string (as it appears in the last column of a PDBQT line).
    pub fn as_str(self) -> &'static str {
        match self {
            VinaType::C  => "C",
            VinaType::A  => "A",
            VinaType::N  => "N",
            VinaType::NA => "NA",
            VinaType::NS => "NS",
            VinaType::O  => "O",
            VinaType::OA => "OA",
            VinaType::S  => "S",
            VinaType::SA => "SA",
            VinaType::P  => "P",
            VinaType::F  => "F",
            VinaType::CL => "Cl",
            VinaType::BR => "Br",
            VinaType::I  => "I",
            VinaType::H  => "H",
            VinaType::HD => "HD",
            VinaType::MG => "Mg",
            VinaType::ZN => "Zn",
            VinaType::FE => "Fe",
            VinaType::CA => "Ca",
            VinaType::MN => "Mn",
            VinaType::CU => "Cu",
            VinaType::Unknown => "?",
        }
    }
}

// ── Assign Vina types to all atoms in a Molecule ─────────────────────────────
// The molecule should have explicit hydrogens already added.

pub fn assign_types(mol: &Molecule) -> Vec<VinaType> {
    let n = mol.atoms.len();
    let mut types = vec![VinaType::Unknown; n];

    for i in 0..n {
        let atom = &mol.atoms[i];
        let nbrs = mol.neighbors(i);
        let sym = atom.symbol.as_str();

        types[i] = match sym {
            "C" => {
                if atom.aromatic {
                    VinaType::A
                } else {
                    let hyb = mol.hybridization(i);
                    if hyb == Hybridization::SP2 {
                        // sp2 non-aromatic carbon (C=O, C=N…) stays C
                        VinaType::C
                    } else {
                        VinaType::C
                    }
                }
            }
            "N" => classify_nitrogen(mol, i, &nbrs),
            "O" => {
                // Oxygen is almost always an H-bond acceptor in organic mols.
                // The only exception would be positively charged O (very rare).
                VinaType::OA
            }
            "S" => VinaType::SA,
            "P" => VinaType::P,
            "F" => VinaType::F,
            "Cl" => VinaType::CL,
            "Br" => VinaType::BR,
            "I" => VinaType::I,
            "H" => classify_hydrogen(mol, i, &nbrs),
            "Mg" => VinaType::MG,
            "Zn" => VinaType::ZN,
            "Fe" => VinaType::FE,
            "Ca" => VinaType::CA,
            "Mn" => VinaType::MN,
            "Cu" => VinaType::CU,
            _ => VinaType::Unknown,
        };
    }

    types
}

fn classify_nitrogen(
    mol: &Molecule,
    idx: usize,
    nbrs: &[(usize, crate::parser::BondOrder)],
) -> VinaType {
    use crate::parser::BondOrder;

    let atom = &mol.atoms[idx];

    // Has an H neighbour → H-bond donor (NS)
    let has_h = nbrs.iter().any(|(j, _)| mol.atoms[*j].symbol == "H");
    if has_h {
        return VinaType::NS;
    }

    // Aromatic N with lone pair available (pyridine-like) → NA
    if atom.aromatic {
        return VinaType::NA;
    }

    // sp2 N with no H and non-aromatic: e.g. amide N, imine N
    // Amide N (bonded to C=O) is not a good H-bond acceptor → N
    // Imine N (C=N-) is acceptor → NA
    let has_double = nbrs.iter().any(|(_, o)| *o == BondOrder::Double);
    if has_double {
        // If the double bond is C=N → this N is an acceptor
        return VinaType::NA;
    }

    // Amine N with no H (tertiary amine) → NA (lone pair available)
    // But only if not positively charged
    if atom.charge == 0 {
        VinaType::NA
    } else {
        VinaType::N
    }
}

fn classify_hydrogen(
    mol: &Molecule,
    _idx: usize,
    nbrs: &[(usize, crate::parser::BondOrder)],
) -> VinaType {
    // Polar H: bonded to N or O → HD (H-bond donor)
    for (j, _) in nbrs {
        let sym = mol.atoms[*j].symbol.as_str();
        if sym == "N" || sym == "O" || sym == "S" {
            return VinaType::HD;
        }
    }
    VinaType::H
}

// ── Assign Vina type from PDB element + atom name ─────────────────────────────
// Used for protein atoms that come from a PDB file (no bond graph available).
// `element`: element symbol from PDB columns 77-78 (e.g. "C", "N", "O")
// `atom_name`: PDB atom name (e.g. "CA", "NZ", "OG1")
// `res_name`: 3-letter residue name (e.g. "ALA", "PHE")
// `is_hetatm`: true if the record is HETATM

pub fn pdb_atom_type(element: &str, atom_name: &str, res_name: &str, is_hetatm: bool) -> VinaType {
    match element {
        "C" => carbon_type_from_pdb(atom_name, res_name),
        "N" => nitrogen_type_from_pdb(atom_name, res_name),
        "O" => VinaType::OA,
        "S" => VinaType::SA,
        "P" => VinaType::P,
        "F" => VinaType::F,
        "CL" | "Cl" => VinaType::CL,
        "BR" | "Br" => VinaType::BR,
        "I" => VinaType::I,
        "H" | "D" => hydrogen_type_from_pdb(atom_name),
        "MG" | "Mg" => VinaType::MG,
        "ZN" | "Zn" => VinaType::ZN,
        "FE" | "Fe" => VinaType::FE,
        "CA" | "Ca" if res_name == "CA " || is_hetatm => VinaType::CA,
        "MN" | "Mn" => VinaType::MN,
        "CU" | "Cu" => VinaType::CU,
        _ => VinaType::Unknown,
    }
}

fn carbon_type_from_pdb(atom_name: &str, res_name: &str) -> VinaType {
    // Aromatic carbons in standard residues
    let arom_res_atoms: &[(&str, &[&str])] = &[
        ("PHE", &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
        ("TYR", &["CG", "CD1", "CD2", "CE1", "CE2", "CZ"]),
        ("TRP", &["CG", "CD1", "CD2", "CE2", "CE3", "CZ2", "CZ3", "CH2"]),
        ("HIS", &["CG", "CD2", "CE1"]),
        ("HID", &["CG", "CD2", "CE1"]),
        ("HIE", &["CG", "CD2", "CE1"]),
        ("HIP", &["CG", "CD2", "CE1"]),
    ];
    let name = atom_name.trim();
    for (res, atoms) in arom_res_atoms {
        if res_name.trim() == *res && atoms.contains(&name) {
            return VinaType::A;
        }
    }
    VinaType::C
}

fn nitrogen_type_from_pdb(atom_name: &str, res_name: &str) -> VinaType {
    let name = atom_name.trim();
    let res = res_name.trim();
    // Nitrogen H-bond donors (have H in protein context)
    // backbone NH, Lys NZ, Arg NH/NE, Trp NE1, His ND1/NE2, Asn ND2, Gln NE2
    match (res, name) {
        ("LYS", "NZ")                => VinaType::NS,
        ("ARG", "NH1" | "NH2" | "NE") => VinaType::NS,
        ("TRP", "NE1")               => VinaType::NS,
        ("HIS" | "HID", "ND1")         => VinaType::NS,
        ("HIE", "NE2")                 => VinaType::NS,
        ("HIP", "ND1" | "NE2")         => VinaType::NS,
        ("ASN", "ND2")                 => VinaType::NS,
        ("GLN", "NE2")                 => VinaType::NS,
        // backbone N is a donor
        (_, "N")                       => VinaType::NS,
        // H-bond acceptor nitrogens
        ("HIE", "ND1")                 => VinaType::NA,
        ("HID", "NE2")                 => VinaType::NA,
        (_, _)                         => VinaType::N,
    }
}

fn hydrogen_type_from_pdb(atom_name: &str) -> VinaType {
    // Polar H on N or O → HD.  Convention: H names starting with H are
    // usually non-polar in carbon context; OH/NH hydrogens have specific names.
    // For protein PDBQT preparation, all H attached to N/O are HD.
    // Since we don't have the bond graph, use atom name heuristics.
    let name = atom_name.trim();
    // Names like HN, HG1, HO, H1, H2, HE… — need context, use H by default.
    // Many workflows omit H from protein PDBQT (non-polar protein prep).
    if name.starts_with("HN") || name == "HO" || name.starts_with("HO") {
        VinaType::HD
    } else {
        VinaType::H
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn types_for(smi: &str) -> Vec<VinaType> {
        let mol = parse(smi).unwrap().with_explicit_hydrogens();
        assign_types(&mol)
    }

    #[test]
    fn benzene_aromatic_carbons() {
        // c1ccccc1: atoms 0-5 aromatic C, 6-11 H
        let t = types_for("c1ccccc1");
        for i in 0..6 {
            assert_eq!(t[i], VinaType::A, "atom {i} should be A");
        }
        for i in 6..12 {
            assert_eq!(t[i], VinaType::H, "atom {i} should be H");
        }
    }

    #[test]
    fn ethane_aliphatic_carbon() {
        let t = types_for("CC");
        assert_eq!(t[0], VinaType::C);
        assert_eq!(t[1], VinaType::C);
    }

    #[test]
    fn ethanol_oxygen_and_polar_h() {
        // CCO: C0, C1, O2, H3-H5 (on C0), H6-H7 (on C1), H8 (on O)
        let t = types_for("CCO");
        // O should be OA
        let o_idx = t.iter().position(|&x| x == VinaType::OA).unwrap();
        assert_eq!(o_idx, 2);
        // H on O should be HD
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let n = mol.atoms.len();
        for i in 0..n {
            if mol.atoms[i].symbol == "H" {
                let bond_to_o = mol.neighbors(i).iter().any(|(j, _)| mol.atoms[*j].symbol == "O");
                if bond_to_o {
                    assert_eq!(t[i], VinaType::HD, "H on O should be HD at idx {i}");
                } else {
                    assert_eq!(t[i], VinaType::H, "H on C should be H at idx {i}");
                }
            }
        }
    }

    #[test]
    fn pyridine_na_nitrogen() {
        // c1ccncc1 — aromatic N should be NA
        let t = types_for("c1ccncc1");
        let mol = parse("c1ccncc1").unwrap().with_explicit_hydrogens();
        let n_idx = mol.atoms.iter().position(|a| a.symbol == "N").unwrap();
        assert_eq!(t[n_idx], VinaType::NA);
    }

    #[test]
    fn amine_ns_nitrogen() {
        // CN: methylamine — N has H → NS
        let t = types_for("CN");
        let mol = parse("CN").unwrap().with_explicit_hydrogens();
        let n_idx = mol.atoms.iter().position(|a| a.symbol == "N").unwrap();
        assert_eq!(t[n_idx], VinaType::NS);
    }

    #[test]
    fn pdb_type_backbone_carbon() {
        assert_eq!(pdb_atom_type("C", "CA", "ALA", false), VinaType::C);
    }

    #[test]
    fn pdb_type_phe_aromatic() {
        assert_eq!(pdb_atom_type("C", "CG", "PHE", false), VinaType::A);
        assert_eq!(pdb_atom_type("C", "CE1", "PHE", false), VinaType::A);
    }

    #[test]
    fn pdb_type_backbone_n() {
        assert_eq!(pdb_atom_type("N", "N", "ALA", false), VinaType::NS);
    }

    #[test]
    fn pdb_type_lys_nz() {
        assert_eq!(pdb_atom_type("N", "NZ", "LYS", false), VinaType::NS);
    }

    #[test]
    fn pdb_type_oxygen_always_oa() {
        assert_eq!(pdb_atom_type("O", "O", "ALA", false), VinaType::OA);
        assert_eq!(pdb_atom_type("O", "OG", "SER", false), VinaType::OA);
    }

    #[test]
    fn as_str_roundtrip() {
        assert_eq!(VinaType::A.as_str(), "A");
        assert_eq!(VinaType::NA.as_str(), "NA");
        assert_eq!(VinaType::OA.as_str(), "OA");
        assert_eq!(VinaType::CL.as_str(), "Cl");
        assert_eq!(VinaType::HD.as_str(), "HD");
    }
}
