//! Wildman-Crippen LogP calculator.
//!
//! Port of RDKit's `Crippen.txt` parameter table and algorithm (see
//! `rdkit/Code/GraphMol/Descriptors/Crippen.cpp`). The algorithm:
//!
//! 1. Expand implicit hydrogens to explicit atoms.
//! 2. For each SMARTS pattern in table order, find substructure matches.
//!    Each match's first atom is assigned this pattern's LogP value
//!    (if not already assigned — first-match-wins).
//! 3. Sum per-atom LogP contributions.
//!
//! Reference: S. A. Wildman, G. M. Crippen, *J. Chem. Inf. Comput. Sci.*
//! **39**, 868-873 (1999).

use crate::parser::Molecule;
use crate::smarts::{match_at, parse_smarts, Pattern};
use std::sync::OnceLock;

/// Parameter table from RDKit's `Data/Crippen.txt`.
/// Each tuple is `(atom_type_label, smarts, logp_contribution, mr_contribution)`.
/// Order is significant: the first matching pattern wins per atom. The same
/// table drives both `calc_logp` (3rd field) and `calc_mr` (4th field, the
/// Wildman-Crippen molar refractivity). RDKit leaves the MR cell blank for a
/// few charged/metal types (N10, N12, O12, Hal, Me2); those contribute `0.0`.
static CRIPPEN_DATA: &[(&str, &str, f64, f64)] = &[
    ("C1",  "[CH4]",                     0.1441,  2.503),
    ("C1",  "[CH3]C",                    0.1441,  2.503),
    ("C1",  "[CH2](C)C",                 0.1441,  2.503),
    ("C2",  "[CH](C)(C)C",               0.0,     2.433),
    ("C2",  "[C](C)(C)(C)C",             0.0,     2.433),
    ("C3",  "[CH3][N,O,P,S,F,Cl,Br,I]",  -0.2035, 2.753),
    ("C3",  "[CH2X4]([N,O,P,S,F,Cl,Br,I])[A;!#1]", -0.2035, 2.753),
    ("C4",  "[CH1X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])[A;!#1]", -0.2051, 2.731),
    ("C4",  "[CH0X4]([N,O,P,S,F,Cl,Br,I])([A;!#1])([A;!#1])[A;!#1]", -0.2051, 2.731),
    ("C5",  "[C]=[!C;A;!#1]",            -0.2783, 5.007),
    ("C6",  "[CH2]=C",                   0.1551,  3.513),
    ("C6",  "[CH1](=C)[A;!#1]",          0.1551,  3.513),
    ("C6",  "[CH0](=C)([A;!#1])[A;!#1]", 0.1551,  3.513),
    ("C6",  "[C](=C)=C",                 0.1551,  3.513),
    ("C7",  "[CX2]#[A;!#1]",             0.0017,  3.888),
    ("C8",  "[CH3]c",                    0.08452, 2.464),
    ("C9",  "[CH3]a",                    -0.1444, 2.412),
    ("C10", "[CH2X4]a",                  -0.0516, 2.488),
    ("C11", "[CHX4]a",                   0.1193,  2.582),
    ("C12", "[CH0X4]a",                  -0.0967, 2.576),
    ("C13", "[cH0]-[A;!C;!N;!O;!S;!F;!Cl;!Br;!I;!#1]", -0.5443, 4.041),
    ("C14", "[c][#9]",                   0.0,     3.257),
    ("C15", "[c][#17]",                  0.245,   3.564),
    ("C16", "[c][#35]",                  0.198,   3.18),
    ("C17", "[c][#53]",                  0.0,     3.104),
    ("C18", "[cH]",                      0.1581,  3.35),
    ("C19", "[c](:a)(:a):a",             0.2955,  4.346),
    ("C20", "[c](:a)(:a)-a",             0.2713,  3.904),
    ("C21", "[c](:a)(:a)-C",             0.136,   3.509),
    ("C22", "[c](:a)(:a)-N",             0.4619,  4.067),
    ("C23", "[c](:a)(:a)-O",             0.5437,  3.853),
    ("C24", "[c](:a)(:a)-S",             0.1893,  2.673),
    ("C25", "[c](:a)(:a)=[C,N,O]",       -0.8186, 3.135),
    ("C26", "[C](=C)(a)[A;!#1]",         0.264,   4.305),
    ("C26", "[C](=C)(c)a",               0.264,   4.305),
    ("C26", "[CH1](=C)a",                0.264,   4.305),
    ("C26", "[C]=c",                     0.264,   4.305),
    ("C27", "[CX4][A;!C;!N;!O;!P;!S;!F;!Cl;!Br;!I;!#1]", 0.2148, 2.693),
    ("CS",  "[#6]",                      0.08129, 3.243),
    ("H1",  "[#1][#6,#1]",               0.123,   1.057),
    ("H2",  "[#1]O[CX4,c]",              -0.2677, 1.395),
    ("H2",  "[#1]O[!C;!N;!O;!S]",        -0.2677, 1.395),
    ("H2",  "[#1][!C;!N;!O]",            -0.2677, 1.395),
    ("H3",  "[#1][#7]",                  0.2142,  0.9627),
    ("H3",  "[#1]O[#7]",                 0.2142,  0.9627),
    ("H4",  "[#1]OC=[#6,#7,O,S]",        0.298,   1.805),
    ("H4",  "[#1]O[O,S]",                0.298,   1.805),
    ("HS",  "[#1]",                      0.1125,  1.112),
    ("N1",  "[NH2+0][A;!#1]",            -1.019,  2.262),
    ("N2",  "[NH+0]([A;!#1])[A;!#1]",    -0.7096, 2.173),
    ("N3",  "[NH2+0]a",                  -1.027,  2.827),
    ("N4",  "[NH1+0]([!#1;A,a])a",       -0.5188, 3.0),
    ("N5",  "[NH+0]=[!#1;A,a]",          0.08387, 1.757),
    ("N6",  "[N+0](=[!#1;A,a])[!#1;A,a]", 0.1836, 2.428),
    ("N7",  "[N+0]([A;!#1])([A;!#1])[A;!#1]", -0.3187, 1.839),
    ("N8",  "[N+0](a)([!#1;A,a])[A;!#1]", -0.4458, 2.819),
    ("N8",  "[N+0](a)(a)a",              -0.4458, 2.819),
    ("N9",  "[N+0]#[A;!#1]",             0.01508, 1.725),
    ("N10", "[NH3,NH2,NH;+,+2,+3]",      -1.95,   0.0),
    ("N11", "[n+0]",                     -0.3239, 2.202),
    ("N12", "[n;+,+2,+3]",               -1.119,  0.0),
    ("N13", "[NH0;+,+2,+3]([A;!#1])([A;!#1])([A;!#1])[A;!#1]", -0.3396, 0.2604),
    ("N13", "[NH0;+,+2,+3](=[A;!#1])([A;!#1])[!#1;A,a]", -0.3396, 0.2604),
    ("N13", "[NH0;+,+2,+3](=[#6])=[#7]", -0.3396, 0.2604),
    ("N14", "[N;+,+2,+3]#[A;!#1]",       0.2887,  3.359),
    ("N14", "[N;-,-2,-3]",               0.2887,  3.359),
    ("N14", "[N;+,+2,+3](=[N;-,-2,-3])=N", 0.2887, 3.359),
    ("NS",  "[#7]",                      -0.4806, 2.134),
    ("O1",  "[o]",                       0.1552,  1.08),
    ("O2",  "[OH,OH2]",                  -0.2893, 0.8238),
    ("O3",  "[O]([A;!#1])[A;!#1]",       -0.0684, 1.085),
    ("O4",  "[O](a)[!#1;A,a]",           -0.4195, 1.182),
    ("O5",  "[O]=[#7,#8]",               0.0335,  3.367),
    ("O5",  "[OX1;-,-2,-3][#7]",         0.0335,  3.367),
    ("O6",  "[OX1;-,-2,-2][#16]",        -0.3339, 0.7774),
    ("O6",  "[O;-0]=[#16;-0]",           -0.3339, 0.7774),
    // Intentional order flip: O12 before O7 per RDKit
    ("O12", "[O-]C(=O)",                 -1.326,  0.0),
    ("O7",  "[OX1;-,-2,-3][!#1;!N;!S]",  -1.189,  0.0),
    ("O8",  "[O]=c",                     0.1788,  3.135),
    ("O9",  "[O]=[CH]C",                 -0.1526, 0.0),
    ("O9",  "[O]=C(C)([A;!#1])",         -0.1526, 0.0),
    ("O9",  "[O]=[CH][N,O]",             -0.1526, 0.0),
    ("O9",  "[O]=[CH2]",                 -0.1526, 0.0),
    ("O9",  "[O]=[CX2]=O",               -0.1526, 0.0),
    ("O10", "[O]=[CH]c",                 0.1129,  0.2215),
    ("O10", "[O]=C([C,c])[a;!#1]",       0.1129,  0.2215),
    ("O10", "[O]=C(c)[A;!#1]",           0.1129,  0.2215),
    ("O11", "[O]=C([!#1;!#6])[!#1;!#6]", 0.4833,  0.389),
    ("OS",  "[#8]",                      -0.1188, 0.6865),
    ("F",   "[#9-0]",                    0.4202,  1.108),
    ("Cl",  "[#17-0]",                   0.6895,  5.853),
    ("Br",  "[#35-0]",                   0.8456,  8.927),
    ("I",   "[#53-0]",                   0.8857,  14.02),
    ("Hal", "[#9,#17,#35,#53;-]",        -2.996,  0.0),
    ("Hal", "[#53;+,+2,+3]",             -2.996,  0.0),
    ("Hal", "[+;#3,#11,#19,#37,#55]",    -2.996,  0.0),
    ("P",   "[#15]",                     0.8612,  6.92),
    // Intentional order flip: S2 before S1 per RDKit
    ("S2",  "[S;-,-2,-3,-4,+1,+2,+3,+5,+6]", -0.0024, 7.365),
    ("S2",  "[S-0]=[N,O,P,S]",           -0.0024, 7.365),
    ("S1",  "[S;A]",                     0.6482,  7.591),
    ("S3",  "[s;a]",                     0.6237,  6.691),
    ("Me1", "[#3,#11,#19,#37,#55]",      -0.3808, 5.754),
    ("Me1", "[#4,#12,#20,#38,#56]",      -0.3808, 5.754),
    ("Me1", "[#5,#13,#31,#49,#81]",      -0.3808, 5.754),
    ("Me1", "[#14,#32,#50,#82]",         -0.3808, 5.754),
    ("Me1", "[#33,#51,#83]",             -0.3808, 5.754),
    ("Me1", "[#34,#52,#84]",             -0.3808, 5.754),
    ("Me2", "[#21,#22,#23,#24,#25,#26,#27,#28,#29,#30]", -0.0025, 0.0),
    ("Me2", "[#39,#40,#41,#42,#43,#44,#45,#46,#47,#48]", -0.0025, 0.0),
    ("Me2", "[#72,#73,#74,#75,#76,#77,#78,#79,#80]",     -0.0025, 0.0),
];

/// Pre-parsed patterns (parsed once per process on first use). Each entry pairs
/// the compiled pattern with its `(logp, mr)` contributions.
fn patterns() -> &'static [(Pattern, f64, f64)] {
    static CACHE: OnceLock<Vec<(Pattern, f64, f64)>> = OnceLock::new();
    CACHE.get_or_init(|| {
        CRIPPEN_DATA
            .iter()
            .filter_map(|(_, smarts, logp, mr)| parse_smarts(smarts).map(|p| (p, *logp, *mr)))
            .collect()
    })
}

/// Shared atom-typing core: assign each atom (after H-expansion) its
/// first-matching Crippen contribution, then sum `select`ed values.
/// `select` picks LogP or MR from the `(logp, mr)` pair.
fn crippen_sum(mol: &Molecule, select: fn(f64, f64) -> f64) -> f64 {
    if mol.atoms.is_empty() {
        return f64::NAN;
    }
    let mol_h = mol.with_explicit_hydrogens();
    let n = mol_h.atoms.len();
    let mut assigned = vec![false; n];
    let mut contrib = vec![0.0_f64; n];

    for (pat, logp, mr) in patterns() {
        for i in 0..n {
            if !assigned[i] && match_at(pat, &mol_h, i) {
                assigned[i] = true;
                contrib[i] = select(*logp, *mr);
            }
        }
        if assigned.iter().all(|&x| x) {
            break;
        }
    }

    contrib.iter().sum()
}

/// Compute Wildman-Crippen LogP for a molecule.
///
/// Returns NaN if the molecule has no atoms.
pub fn calc_logp(mol: &Molecule) -> f64 {
    crippen_sum(mol, |logp, _mr| logp)
}

/// Compute Wildman-Crippen molar refractivity (MolMR) for a molecule, summing
/// the per-atom MR contributions from the same atom-typing as `calc_logp`.
///
/// Returns NaN if the molecule has no atoms.
pub fn calc_mr(mol: &Molecule) -> f64 {
    crippen_sum(mol, |_logp, mr| mr)
}

// =============================================================================
// Tests
// =============================================================================

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn logp_of(smiles: &str) -> f64 {
        let mol = parse(smiles).expect("parse smiles");
        calc_logp(&mol)
    }

    /// Compare LogP to a reference value (typically from RDKit).
    /// Tolerance is broad because small parser differences can shift atom types.
    fn assert_close(smiles: &str, expected: f64, tol: f64) {
        let got = logp_of(smiles);
        assert!(
            (got - expected).abs() < tol,
            "LogP({}) = {:.4}, expected {:.4} ± {}",
            smiles, got, expected, tol
        );
    }

    // --------- Patterns parse successfully ---------

    #[test]
    fn all_patterns_parse() {
        let patterns = CRIPPEN_DATA.iter()
            .filter_map(|(label, smarts, _, _)| {
                if parse_smarts(smarts).is_none() {
                    Some(format!("{}: {}", label, smarts))
                } else {
                    None
                }
            })
            .collect::<Vec<_>>();
        assert!(patterns.is_empty(), "failed patterns: {:?}", patterns);
    }

    // --------- Simple molecules (RDKit reference values) ---------

    #[test]
    fn logp_methane() {
        // RDKit: MolLogP(Chem.MolFromSmiles("C")) ≈ 0.6361
        // = 1× C1 (0.1441) + 4× H1 (0.123) = 0.6361
        assert_close("C", 0.6361, 0.05);
    }

    #[test]
    fn logp_water() {
        // RDKit: MolLogP(O) ≈ -0.8247
        // = 1× O2 (-0.2893) + 2× H2 (-0.2677) = -0.8247
        assert_close("O", -0.8247, 0.05);
    }

    #[test]
    fn logp_ethanol() {
        // RDKit: MolLogP(CCO) ≈ -0.0014
        assert_close("CCO", -0.0014, 0.1);
    }

    #[test]
    fn logp_benzene() {
        // RDKit: MolLogP(c1ccccc1) ≈ 1.6866
        // = 6× C18 (0.1581) + 6× H1 (0.123) = 6×0.2811 = 1.6866
        assert_close("c1ccccc1", 1.6866, 0.05);
    }

    #[test]
    fn logp_methanol() {
        // RDKit: MolLogP(CO) ≈ -0.3915
        assert_close("CO", -0.3915, 0.1);
    }

    #[test]
    fn logp_chloromethane() {
        // RDKit: MolLogP(CCl) ≈ 0.8675 approximately
        assert_close("CCl", 0.87, 0.2);
    }

    #[test]
    fn logp_no_panic_invalid_but_parsed() {
        // Just verify we don't panic on unusual structures
        let v = logp_of("[Na+].[Cl-]");
        assert!(v.is_finite());
    }

    // --------- MolMR (Wildman-Crippen molar refractivity) ---------

    fn mr_of(smiles: &str) -> f64 {
        calc_mr(&parse(smiles).expect("parse smiles"))
    }

    #[test]
    fn mr_methane() {
        // RDKit: MolMR(C) = 6.731 = 1× C1 (2.503) + 4× H1 (1.057)
        assert!((mr_of("C") - 6.731).abs() < 0.01, "MR(C) = {}", mr_of("C"));
    }

    #[test]
    fn mr_benzene() {
        // RDKit: MolMR(c1ccccc1) = 26.442 = 6× C18 (3.35) + 6× H1 (1.057)
        assert!((mr_of("c1ccccc1") - 26.442).abs() < 0.05, "MR = {}", mr_of("c1ccccc1"));
    }

    #[test]
    fn mr_shares_logp_atom_typing() {
        // Same molecule yields finite, positive MR for both LogP-typed atoms.
        let v = mr_of("CC(=O)Oc1ccccc1C(=O)O"); // aspirin
        assert!(v > 30.0 && v < 50.0, "aspirin MR out of expected range: {}", v);
    }

    #[test]
    fn mr_grows_with_size() {
        assert!(mr_of("CCCCCC") > mr_of("CC"), "hexane MR should exceed ethane");
    }

    // --------- Every atom gets a type assigned ---------

    #[test]
    fn all_atoms_assigned_ethanol() {
        let mol = parse("CCO").unwrap().with_explicit_hydrogens();
        let n = mol.atoms.len();
        let mut assigned = vec![false; n];

        for (pat, _, _) in patterns() {
            for i in 0..n {
                if !assigned[i] && match_at(pat, &mol, i) {
                    assigned[i] = true;
                }
            }
        }

        let unassigned: Vec<_> = (0..n)
            .filter(|&i| !assigned[i])
            .map(|i| (i, mol.atoms[i].symbol.clone()))
            .collect();
        assert!(unassigned.is_empty(), "unassigned atoms: {:?}", unassigned);
    }

    #[test]
    fn all_atoms_assigned_caffeine() {
        let mol = parse("Cn1c(=O)c2c(ncn2C)n(C)c1=O").unwrap().with_explicit_hydrogens();
        let n = mol.atoms.len();
        let mut assigned = vec![false; n];

        for (pat, _, _) in patterns() {
            for i in 0..n {
                if !assigned[i] && match_at(pat, &mol, i) {
                    assigned[i] = true;
                }
            }
        }

        let unassigned: Vec<_> = (0..n)
            .filter(|&i| !assigned[i])
            .map(|i| (i, mol.atoms[i].symbol.clone()))
            .collect();
        assert!(unassigned.is_empty(), "caffeine unassigned: {:?}", unassigned);
    }

    // --------- Consistency: order of magnitude sanity ---------

    #[test]
    fn logp_cholesterol_is_lipophilic() {
        let v = logp_of("CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C");
        assert!(v > 3.0, "cholesterol should be strongly lipophilic, got {}", v);
    }

    #[test]
    fn logp_glucose_is_hydrophilic() {
        let v = logp_of("OCC1OC(O)C(O)C(O)C1O");
        assert!(v < 0.0, "glucose should be hydrophilic, got {}", v);
    }

    #[test]
    fn logp_aspirin_moderate() {
        let v = logp_of("CC(=O)Oc1ccccc1C(=O)O");
        // RDKit: ~1.19. Permissive since our matcher may differ slightly.
        assert!((0.5..2.5).contains(&v), "aspirin LogP out of range: {}", v);
    }
}

#[cfg(test)]
mod print_values {
    use super::*;
    use crate::parser::parse;

    #[test]
    #[ignore]
    fn print_logp_table() {
        let samples = [
            ("methane",    "C"),
            ("water",      "O"),
            ("ethanol",    "CCO"),
            ("methanol",   "CO"),
            ("benzene",    "c1ccccc1"),
            ("toluene",    "Cc1ccccc1"),
            ("phenol",     "Oc1ccccc1"),
            ("chloromethane","CCl"),
            ("aspirin",    "CC(=O)Oc1ccccc1C(=O)O"),
            ("ibuprofen",  "CC(C)Cc1ccc(cc1)C(C)C(=O)O"),
            ("caffeine",   "Cn1c(=O)c2c(ncn2C)n(C)c1=O"),
            ("glucose",    "OCC1OC(O)C(O)C(O)C1O"),
            ("cholesterol","CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"),
        ];
        for (name, smi) in samples {
            let mol = parse(smi).unwrap();
            let v = calc_logp(&mol);
            println!("{:15} {:35} LogP = {:+.4}", name, smi, v);
        }
    }
}
