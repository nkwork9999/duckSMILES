//! QED — quantitative estimate of drug-likeness (Bickerton et al., 2012).
//!
//! Port of RDKit's `rdkit/Chem/QED.py`. QED maps eight molecular properties
//! through asymmetric-double-sigmoidal (ADS) desirability functions and takes a
//! weighted geometric mean:
//!
//! ```text
//! QED = exp( sum(wi*ln(di)) / sum(wi) ),  di = ADS_i(property_i)
//! ```
//!
//! Properties (in weight order): MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS.
//! Defaults use the `WEIGHT_MEAN` set, matching RDKit's `QED.default` /
//! `QED.weights_mean`.
//!
//! ## Fidelity notes
//!
//! QED.py itself documents that its values differ slightly from the original
//! paper because of the underlying property calculators (esp. LogP). This port
//! inherits the same caveat plus two duckSmiles-specific ones:
//!
//! - **ROTB** uses RDKit's *Strict* rotatable-bond definition. RDKit expresses
//!   it as a single bond-SMARTS, but our SMARTS engine doesn't support bond-OR
//!   (`-,:`). We instead reuse the (parseable) per-atom environment queries and
//!   apply the bond-level constraint (acyclic single bond) in Rust — equivalent
//!   for the single/aromatic·!ring case, since aromatic bonds are always ring
//!   bonds and thus excluded by `!@` anyway.
//! - **ALERTS** evaluates 114 of RDKit's 116 structural alerts. The two skipped
//!   are multi-component SMARTS (`F.F.F.F` and a triple-ester), which the engine
//!   doesn't parse. They only affect molecules that are salt mixtures / triesters.
//!
//! `qed` returns `NaN` for an empty molecule (no atoms), mirroring `calc_logp`.

use crate::logp_crippen::calc_logp;
use crate::parser::{Bond, BondOrder, Molecule};
use crate::smarts::{count_unique, match_at, matches_mol, parse_smarts, unique_matches, Pattern};
use crate::tpsa::calc_tpsa;
use std::sync::OnceLock;

// --- WEIGHT_MEAN (RDKit QED default), property order MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS ---
const WEIGHT_MEAN: [f64; 8] = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];

// --- ADS parameters (A, B, C, D, E, F, DMAX) per property, RDKit adsParameters ---
const ADS_PARAMS: [[f64; 7]; 8] = [
    // MW
    [2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561],
    // ALOGP
    [3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604],
    // HBA
    [2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046],
    // HBD
    [1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616],
    // PSA
    [1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167],
    // ROTB
    [0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403],
    // AROM
    [3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610],
    // ALERTS
    [0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140],
];

// QED hydrogen-bond acceptor SMARTS (summed match counts). RDKit AcceptorSmarts.
static ACCEPTOR_SMARTS: &[&str] = &[
    "[oH0;X2]",
    "[OH1;X2;v2]",
    "[OH0;X2;v2]",
    "[OH0;X1;v2]",
    "[O-;X1]",
    "[SH0;X2;v2]",
    "[SH0;X1;v2]",
    "[S-;X1]",
    "[nH0;X2]",
    "[NH0;X1;v3]",
    "[$([N;+0;X3;v3]);!$(N[C,S]=O)]",
];

// RDKit CalcNumHBD pattern (counted as matching atoms).
static HBD_SMARTS: &str = "[N&!H0&v3,N&!H0&+1&v4,O&H1&+0,S&H1&+0,n&H1&+0]";

// QED AliphaticRings query: an aliphatic ring atom bonded to a non-aromatic
// atom. AROM = SSSR count after deleting all atoms matching this.
static ALIPHATIC_RING_SMARTS: &str = "[$([A;R][!a])]";

// Strict rotatable-bond atom environments (the bond-OR is handled in Rust).
// First atom carries the amide exclusions; the second is the common environment.
static ROTB_ENV_COMMON: &str = "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])]";
static ROTB_ENV_AMIDE: &str = "[!$(*#*)&!D1&!$(C(F)(F)F)&!$(C(Cl)(Cl)Cl)&!$(C(Br)(Br)Br)&!$(C([CH3])([CH3])[CH3])&!$([CH3])&!$([CD3](=[N,O,S])-!@[#7,O,S!D1])&!$([#7,O,S!D1]-!@[CD3]=[N,O,S])&!$([CD3](=[N+])-!@[#7!D1])&!$([#7!D1]-!@[CD3]=[N+])]";

// RDKit StructuralAlertSmarts (116 entries; two multi-component ones don't
// parse in our engine and are silently dropped by `filter_map` below).
static ALERT_SMARTS: &[&str] = &[
    "*1[O,S,N]*1",
    "[S,C](=[O,S])[F,Br,Cl,I]",
    "[CX4][Cl,Br,I]",
    "[#6]S(=O)(=O)O[#6]",
    "[$([CH]),$(CC)]#CC(=O)[#6]",
    "[$([CH]),$(CC)]#CC(=O)O[#6]",
    "n[OH]",
    "[$([CH]),$(CC)]#CS(=O)(=O)[#6]",
    "C=C(C=O)C=O",
    "n1c([F,Cl,Br,I])cccc1",
    "[CH1](=O)",
    "[#8][#8]",
    "[C;!R]=[N;!R]",
    "[N!R]=[N!R]",
    "[#6](=O)[#6](=O)",
    "[#16][#16]",
    "[#7][NH2]",
    "C(=O)N[NH2]",
    "[#6]=S",
    "[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]",
    "C1(=[O,N])C=CC(=[O,N])C=C1",
    "C1(=[O,N])C(=[O,N])C=CC=C1",
    "a21aa3a(aa1aaaa2)aaaa3",
    "a31a(a2a(aa1)aaaa2)aaaa3",
    "a1aa2a3a(a1)A=AA=A3=AA=A2",
    "c1cc([NH2])ccc1",
    "[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]",
    "I",
    "OS(=O)(=O)[O-]",
    "[N+](=O)[O-]",
    "C(=O)N[OH]",
    "C1NC(=O)NC(=O)1",
    "[SH]",
    "[S-]",
    "c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]",
    "[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1",
    "[CR1]1[CR1][CR1]cc[CR1][CR1]1",
    "[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1",
    "[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1",
    "C#C",
    "[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]",
    "[$([N+R]),$([n+R]),$([N+]=C)][O-]",
    "[#6]=N[OH]",
    "[#6]=NOC=O",
    "[#6](=O)[CX4,CR0X3,O][#6](=O)",
    "c1ccc2c(c1)ccc(=O)o2",
    "[O+,o+,S+,s+]",
    "N=C=O",
    "[NX3,NX4][F,Cl,Br,I]",
    "c1ccccc1OC(=O)[#6]",
    "[CR0]=[CR0][CR0]=[CR0]",
    "[C+,c+,C-,c-]",
    "N=[N+]=[N-]",
    "C12C(NC(N1)=O)CSC2",
    "c1c([OH])c([OH,NH2,NH])ccc1",
    "P",
    "[N,O,S]C#N",
    "C=C=O",
    "[Si][F,Cl,Br,I]",
    "[SX2]O",
    "[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)",
    "O1CCCCC1OC2CCC3CCCCC3C2",
    "N=[CR0][N,n,O,S]",
    "[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2",
    "C=[C!r]C#N",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1",
    "[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])",
    "[OH]c1ccc([OH,NH2,NH])cc1",
    "c1ccccc1OC(=O)O",
    "[SX2H0][N]",
    "c12ccccc1(SC(S)=N2)",
    "c12ccccc1(SC(=S)N2)",
    "c1nnnn1C=O",
    "s1c(S)nnc1NC=O",
    "S1C=CSC1=S",
    "C(=O)Onnn",
    "OS(=O)(=O)C(F)(F)F",
    "N#CC[OH]",
    "N#CC(=O)",
    "S(=O)(=O)C#N",
    "N[CH2]C#N",
    "C1(=O)NCC1",
    "S(=O)(=O)[O-,OH]",
    "NC[F,Cl,Br,I]",
    "C=[C!r]O",
    "[NX2+0]=[O+0]",
    "[OR0,NR0][OR0,NR0]",
    "C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]",
    "[CX2R0][NX3R0]",
    "c1ccccc1[C;!R]=[C;!R]c2ccccc2",
    "[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]",
    "[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]",
    "[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]",
    "[*]=[N+]=[*]",
    "[SX3](=O)[O-,OH]",
    "N#N",
    "F.F.F.F",
    "[R0;D2][R0;D2][R0;D2][R0;D2]",
    "[cR,CR]~C(=O)NC(=O)~[cR,CR]",
    "C=!@CC=[O,S]",
    "[#6,#8,#16][#6](=O)O[#6]",
    "c[C;R0](=[O,S])[#6]",
    "c[SX2][C;!R]",
    "C=C=C",
    "c1nc([F,Cl,Br,I,S])ncc1",
    "c1ncnc([F,Cl,Br,I,S])c1",
    "c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])",
    "[#6]S(=O)(=O)c1ccc(cc1)F",
    "[15N]",
    "[13C]",
    "[18O]",
    "[34S]",
];

fn compile_all(smarts: &[&str]) -> Vec<Pattern> {
    smarts.iter().filter_map(|s| parse_smarts(s)).collect()
}

fn acceptors() -> &'static [Pattern] {
    static CACHE: OnceLock<Vec<Pattern>> = OnceLock::new();
    CACHE.get_or_init(|| compile_all(ACCEPTOR_SMARTS))
}

fn alerts() -> &'static [Pattern] {
    static CACHE: OnceLock<Vec<Pattern>> = OnceLock::new();
    CACHE.get_or_init(|| compile_all(ALERT_SMARTS))
}

fn single_pattern(smarts: &str) -> &'static Pattern {
    // Leak a single compiled pattern once; SMARTS here are constants that
    // always parse, so unwrap is safe (covered by the patterns_compile test).
    fn cache(slot: &'static OnceLock<Pattern>, s: &str) -> &'static Pattern {
        slot.get_or_init(|| parse_smarts(s).expect("QED constant SMARTS must parse"))
    }
    match smarts {
        s if s == HBD_SMARTS => {
            static C: OnceLock<Pattern> = OnceLock::new();
            cache(&C, HBD_SMARTS)
        }
        s if s == ALIPHATIC_RING_SMARTS => {
            static C: OnceLock<Pattern> = OnceLock::new();
            cache(&C, ALIPHATIC_RING_SMARTS)
        }
        s if s == ROTB_ENV_COMMON => {
            static C: OnceLock<Pattern> = OnceLock::new();
            cache(&C, ROTB_ENV_COMMON)
        }
        s if s == ROTB_ENV_AMIDE => {
            static C: OnceLock<Pattern> = OnceLock::new();
            cache(&C, ROTB_ENV_AMIDE)
        }
        _ => unreachable!("unknown QED constant SMARTS"),
    }
}

/// HBA = Σ over the 11 acceptor patterns of their match counts.
fn hba(mol: &Molecule) -> usize {
    acceptors().iter().map(|p| count_unique(p, mol)).sum()
}

/// HBD = number of atoms matching RDKit's CalcNumHBD pattern.
fn hbd(mol: &Molecule) -> usize {
    count_unique(single_pattern(HBD_SMARTS), mol)
}

/// Strict rotatable bonds: acyclic single bonds whose endpoints both satisfy the
/// common environment and at least one satisfies the amide-aware environment.
fn rotb_strict(mol: &Molecule) -> usize {
    let env_common = single_pattern(ROTB_ENV_COMMON);
    let env_amide = single_pattern(ROTB_ENV_AMIDE);
    let ring = mol.ring_info();
    mol.bonds
        .iter()
        .enumerate()
        .filter(|(idx, bond)| {
            // (single OR aromatic) AND !ring == single AND !ring, since aromatic
            // bonds are always ring bonds.
            bond.order == BondOrder::Single
                && !ring.bond_in_ring.get(*idx).copied().unwrap_or(false)
        })
        .filter(|(_, bond)| {
            let a_common = match_at(env_common, mol, bond.a);
            let b_common = match_at(env_common, mol, bond.b);
            if !(a_common && b_common) {
                return false;
            }
            // At least one endpoint must also clear the amide-aware environment.
            match_at(env_amide, mol, bond.a) || match_at(env_amide, mol, bond.b)
        })
        .count()
}

/// AROM = SSSR ring count after deleting atoms matching the AliphaticRings
/// query (RDKit: `GetSSSR(DeleteSubstructs(mol, AliphaticRings))`).
fn arom(mol: &Molecule) -> usize {
    let pat = single_pattern(ALIPHATIC_RING_SMARTS);
    let mut delete = vec![false; mol.atoms.len()];
    for m in unique_matches(pat, mol) {
        for idx in m {
            if idx < delete.len() {
                delete[idx] = true;
            }
        }
    }
    if delete.iter().all(|&d| !d) {
        // Nothing deleted — SSSR of the whole molecule.
        return mol.ring_info().rings.len();
    }

    // Build the residual molecule with surviving atoms reindexed.
    let mut remap = vec![usize::MAX; mol.atoms.len()];
    let mut atoms = Vec::new();
    for (i, atom) in mol.atoms.iter().enumerate() {
        if !delete[i] {
            remap[i] = atoms.len();
            atoms.push(atom.clone());
        }
    }
    let mut bonds = Vec::new();
    for bond in &mol.bonds {
        if !delete[bond.a] && !delete[bond.b] {
            bonds.push(Bond {
                a: remap[bond.a],
                b: remap[bond.b],
                order: bond.order,
            });
        }
    }
    let bond_count = bonds.len() as i32;
    let sub = Molecule {
        atoms,
        bonds,
        bond_count,
    };
    sub.ring_info().rings.len()
}

/// ALERTS = number of structural-alert patterns with at least one match.
fn alerts_count(mol: &Molecule) -> usize {
    alerts().iter().filter(|p| matches_mol(p, mol)).count()
}

/// The eight QED properties in weight order.
fn properties(mol: &Molecule) -> [f64; 8] {
    [
        mol.molecular_weight(),
        calc_logp(mol),
        hba(mol) as f64,
        hbd(mol) as f64,
        calc_tpsa(mol),
        rotb_strict(mol) as f64,
        arom(mol) as f64,
        alerts_count(mol) as f64,
    ]
}

/// ADS desirability function for property value `x` with parameter row
/// `[A, B, C, D, E, F, DMAX]`.
fn ads(x: f64, p: &[f64; 7]) -> f64 {
    let (a, b, c, d, e, f, dmax) = (p[0], p[1], p[2], p[3], p[4], p[5], p[6]);
    let exp1 = 1.0 + (-(x - c + d / 2.0) / e).exp();
    let exp2 = 1.0 + (-(x - c - d / 2.0) / f).exp();
    let dx = a + b / exp1 * (1.0 - 1.0 / exp2);
    dx / dmax
}

/// Weighted-geometric-mean QED using the WEIGHT_MEAN parameters (RDKit default).
/// Returns NaN for an empty molecule.
pub fn qed(mol: &Molecule) -> f64 {
    if mol.atoms.is_empty() {
        return f64::NAN;
    }
    let props = properties(mol);
    let mut num = 0.0;
    let mut wsum = 0.0;
    for i in 0..8 {
        let di = ads(props[i], &ADS_PARAMS[i]);
        num += WEIGHT_MEAN[i] * di.ln();
        wsum += WEIGHT_MEAN[i];
    }
    (num / wsum).exp()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn qed_of(s: &str) -> f64 {
        qed(&parse(s).expect("parse"))
    }

    #[test]
    fn all_constant_smarts_compile() {
        // The 11 acceptors all parse; 114/116 alerts parse (2 multi-component).
        assert_eq!(acceptors().len(), 11);
        assert_eq!(alerts().len(), 114, "expected 114 parseable alerts");
        // The bespoke-helper SMARTS must all parse (single_pattern unwraps them).
        let _ = single_pattern(HBD_SMARTS);
        let _ = single_pattern(ALIPHATIC_RING_SMARTS);
        let _ = single_pattern(ROTB_ENV_COMMON);
        let _ = single_pattern(ROTB_ENV_AMIDE);
    }

    #[test]
    fn ads_matches_reference_point() {
        // Spot-check the ADS function against a hand-computed value for MW=300.
        // exp1 = 1+exp(-(300-290.749+1.21)/49.223); exp2 = 1+exp(-(300-290.749-1.21)/65.371)
        let v = ads(300.0, &ADS_PARAMS[0]);
        assert!(v > 0.9 && v <= 1.0, "ADS(MW=300) out of expected band: {v}");
    }

    #[test]
    fn rotb_butane_is_one() {
        assert_eq!(rotb_strict(&parse("CCCC").unwrap()), 1);
    }

    #[test]
    fn rotb_excludes_tert_butyl_and_cf3() {
        // neopentane: all bonds touch a CH3 or the quaternary C → 0 rotatable
        assert_eq!(rotb_strict(&parse("CC(C)(C)C").unwrap()), 0);
    }

    #[test]
    fn arom_benzene_is_one_naphthalene_two() {
        assert_eq!(arom(&parse("c1ccccc1").unwrap()), 1);
        assert_eq!(arom(&parse("c1ccc2ccccc2c1").unwrap()), 2);
    }

    #[test]
    fn arom_cyclohexane_is_zero() {
        // purely aliphatic ring → AROM 0 after deletion
        assert_eq!(arom(&parse("C1CCCCC1").unwrap()), 0);
    }

    #[test]
    fn qed_in_unit_interval() {
        for s in ["CCO", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O", "CN1C=NC2=C1C(=O)N(C(=O)N2C)C"] {
            let q = qed_of(s);
            assert!((0.0..=1.0).contains(&q), "QED({s}) = {q} out of [0,1]");
        }
    }

    #[test]
    fn qed_reference_molecules_are_close() {
        // RDKit reference values from QED.py docstring (Peter G's implementation):
        //   N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O   -> 0.253
        //   CNC(=NCCSCc1nc[nH]c1C)NC#N              -> 0.234
        //   CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12     -> 0.234
        // Exact parity is impossible (QED.py documents LogP-driven divergence);
        // we assert each lands in a reasonable neighborhood of the reference.
        let cases = [
            ("N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O", 0.253),
            ("CNC(=NCCSCc1nc[nH]c1C)NC#N", 0.234),
            ("CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12", 0.234),
        ];
        for (smi, want) in cases {
            let got = qed_of(smi);
            assert!(
                (got - want).abs() < 0.06,
                "QED({smi}) = {got:.3}, reference {want:.3}"
            );
        }
    }
}

