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
//! The pinned registered definition removes ordinary explicit hydrogen vertices
//! with the default RemoveHs rules before computing the eight properties.
//! QED's acceptor queries differ from the standard NumHAcceptors descriptor.
//! Structural alerts are pattern-presence counts, not occurrence counts.
//! Empty inputs retain the existing native NaN contract.

use crate::logp_crippen::calc_logp;
use crate::parser::{BondOrder, Molecule, RingInfo};
use crate::smarts::{
    multi_matches_mol, parse_smarts, parse_smarts_multi, MultiPattern, Pattern, PreparedMatcher,
};
use std::sync::OnceLock;

// --- WEIGHT_MEAN (RDKit QED default), property order MW,ALOGP,HBA,HBD,PSA,ROTB,AROM,ALERTS ---
const WEIGHT_MEAN: [f64; 8] = [0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95];

// --- ADS parameters (A, B, C, D, E, F, DMAX) per property, RDKit adsParameters ---
const ADS_PARAMS: [[f64; 7]; 8] = [
    // MW
    [
        2.817065973,
        392.5754953,
        290.7489764,
        2.419764353,
        49.22325677,
        65.37051707,
        104.9805561,
    ],
    // ALOGP
    [
        3.172690585,
        137.8624751,
        2.534937431,
        4.581497897,
        0.822739154,
        0.576295591,
        131.3186604,
    ],
    // HBA
    [
        2.948620388,
        160.4605972,
        3.615294657,
        4.435986202,
        0.290141953,
        1.300669958,
        148.7763046,
    ],
    // HBD
    [
        1.618662227,
        1010.051101,
        0.985094388,
        0.000000001,
        0.713820843,
        0.920922555,
        258.1632616,
    ],
    // PSA
    [
        1.876861559,
        125.2232657,
        62.90773554,
        87.83366614,
        12.01999824,
        28.51324732,
        104.5686167,
    ],
    // ROTB
    [
        0.010000000,
        272.4121427,
        2.558379970,
        1.565547684,
        1.271567166,
        2.758063707,
        105.4420403,
    ],
    // AROM
    [
        3.217788970,
        957.7374108,
        2.274627939,
        0.000000001,
        1.317690384,
        0.375760881,
        312.3372610,
    ],
    // ALERTS
    [
        0.010000000,
        1199.094025,
        -0.09002883,
        0.000000001,
        0.185904477,
        0.875193782,
        417.7253140,
    ],
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

// AROM is the ordinary SSSR dimension after deleting this query's atoms.
static ALIPHATIC_RING_SMARTS: &str = "[$([A;R][!a])]";

// RDKit StructuralAlertSmarts (116 entries, two of them dot-separated).
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

/// Alerts include two dot-separated patterns (`F.F.F.F` and a triple ester),
/// so they are compiled as multi-component SMARTS rather than dropped.
fn compile_all_multi(smarts: &[&str]) -> Vec<MultiPattern> {
    smarts
        .iter()
        .filter_map(|s| parse_smarts_multi(s))
        .collect()
}

fn acceptors() -> &'static [Pattern] {
    static CACHE: OnceLock<Vec<Pattern>> = OnceLock::new();
    CACHE.get_or_init(|| compile_all(ACCEPTOR_SMARTS))
}

fn alerts() -> &'static [MultiPattern] {
    static CACHE: OnceLock<Vec<MultiPattern>> = OnceLock::new();
    CACHE.get_or_init(|| compile_all_multi(ALERT_SMARTS))
}

fn aliphatic_pattern() -> &'static Pattern {
    static CACHE: OnceLock<Pattern> = OnceLock::new();
    CACHE.get_or_init(|| parse_smarts(ALIPHATIC_RING_SMARTS).expect("constant QED query"))
}

/// Default RemoveHs for the supported concrete SMILES graph. Isotopic,
/// isolated, bridging, hydride, H-only-neighbor and defining-stereo Hs remain.
/// The original graph is untouched. Attached-H totals transfer to the neighbor.
fn hydrogen_suppressed(mol: &Molecule) -> Option<Molecule> {
    let adj = mol.adjacency();
    let removed: Vec<_> = mol
        .atoms
        .iter()
        .enumerate()
        .map(|(i, atom)| {
            if atom.symbol != "H"
                || atom.isotope.is_some_and(|v| v != 0)
                || atom.charge == -1
                || adj[i].len() != 1
            {
                return false;
            }
            let (neighbor, bond) = adj[i][0];
            if mol.atoms[neighbor].symbol == "H" {
                return false;
            }
            // A directed single H bond at a terminal double-bond endpoint defines
            // its stereochemistry and is retained by default RemoveHs.
            !(adj[neighbor].len() == 2
                && mol.bonds[bond].direction != 0
                && adj[neighbor]
                    .iter()
                    .any(|&(_, b)| mol.bonds[b].order == BondOrder::Double))
        })
        .collect();
    if !removed.iter().any(|&r| r) {
        return None;
    }
    let mut remap = vec![usize::MAX; mol.atoms.len()];
    let mut atoms = Vec::new();
    for (i, atom) in mol.atoms.iter().enumerate() {
        if removed[i] {
            continue;
        }
        remap[i] = atoms.len();
        let mut atom = atom.clone();
        let transferred = adj[i].iter().filter(|&&(j, _)| removed[j]).count() as i32;
        atom.hydrogen += transferred;
        atom.nbr_order.clear();
        atoms.push(atom);
    }
    let bonds: Vec<_> = mol
        .bonds
        .iter()
        .filter(|b| !removed[b.a] && !removed[b.b])
        .map(|b| {
            let mut b = *b;
            b.a = remap[b.a];
            b.b = remap[b.b];
            b
        })
        .collect();
    Some(Molecule {
        atoms,
        bond_count: bonds.len() as i32,
        bonds,
    })
}

fn acceptor_counts(mol: &Molecule, matcher: &PreparedMatcher<'_>) -> Vec<usize> {
    // All eleven acceptor patterns have one top-level query atom. Recursive
    // environments constrain that atom; they do not add mapped atoms.
    acceptors()
        .iter()
        .map(|p| {
            (0..mol.atoms.len())
                .filter(|&i| matcher.matches(p, i))
                .take(1000)
                .count()
        })
        .collect()
}

fn alert_flags(mol: &Molecule, matcher: &PreparedMatcher<'_>) -> Vec<bool> {
    alerts()
        .iter()
        .map(|p| {
            if p.components.len() == 1 {
                (0..mol.atoms.len()).any(|i| matcher.matches(&p.components[0], i))
            } else {
                // Only two queries need disjoint multi-component embeddings.
                multi_matches_mol(p, mol)
            }
        })
        .collect()
}

fn arom_with_context(mol: &Molecule, matcher: &PreparedMatcher<'_>) -> usize {
    let deleted: Vec<_> = (0..mol.atoms.len())
        .map(|i| matcher.matches(aliphatic_pattern(), i))
        .collect();
    let mut parent: Vec<_> = (0..mol.atoms.len()).collect();
    fn root(parent: &mut [usize], mut i: usize) -> usize {
        while parent[i] != i {
            parent[i] = parent[parent[i]];
            i = parent[i];
        }
        i
    }
    let mut cycles = 0;
    for b in &mol.bonds {
        if deleted[b.a] || deleted[b.b] || b.order == BondOrder::Dative {
            continue;
        }
        let a = root(&mut parent, b.a);
        let c = root(&mut parent, b.b);
        if a == c {
            cycles += 1;
        } else {
            parent[a] = c;
        }
    }
    // E - V + components is the ordinary SSSR dimension. Symmetric SSSR can
    // contain extra rings and is deliberately not used for this property.
    cycles
}

fn properties_on_graph(
    mol: &Molecule,
    ring: &RingInfo,
    base: Option<&crate::workbench::BasicDescriptors>,
) -> [f64; 8] {
    let matcher = PreparedMatcher::with_ring_info(mol, ring);
    [
        base.map_or_else(|| mol.molecular_weight(), |b| b.mw),
        base.map_or_else(|| calc_logp(mol), |b| b.logp),
        acceptor_counts(mol, &matcher).iter().sum::<usize>() as f64,
        base.map_or_else(|| crate::num_h_donors(mol) as f64, |b| b.hbd as f64),
        base.map_or_else(
            || crate::tpsa::calc_tpsa_with_ring_info(mol, ring),
            |b| b.tpsa,
        ),
        base.map_or_else(
            || crate::num_rotatable_bonds_with_ring_bonds(mol, &ring.bond_in_ring) as f64,
            |b| b.rotatable_bonds as f64,
        ),
        arom_with_context(mol, &matcher) as f64,
        alert_flags(mol, &matcher).iter().filter(|&&v| v).count() as f64,
    ]
}

/// QED has its own default-H-removal graph; the input graph stays intact.
fn properties(mol: &Molecule) -> [f64; 8] {
    let suppressed = hydrogen_suppressed(mol);
    let graph = suppressed.as_ref().unwrap_or(mol);
    properties_on_graph(graph, &graph.ring_info(), None)
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
    score(properties(mol))
}

fn score(props: [f64; 8]) -> f64 {
    let mut num = 0.0;
    let mut wsum = 0.0;
    for i in 0..8 {
        let di = ads(props[i], &ADS_PARAMS[i]);
        num += WEIGHT_MEAN[i] * di.ln();
        wsum += WEIGHT_MEAN[i];
    }
    (num / wsum).exp()
}

/// Registered workbench path reuses existing properties and ring topology
/// unless QED's default explicit-H removal needs a separate graph.
pub(crate) fn calculate(
    mol: &Molecule,
    ring: &RingInfo,
    base: &crate::workbench::BasicDescriptors,
) -> f64 {
    match hydrogen_suppressed(mol) {
        Some(graph) => score(properties_on_graph(&graph, &graph.ring_info(), None)),
        None => score(properties_on_graph(mol, ring, Some(base))),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn prepared(source: &str) -> Molecule {
        crate::workbench::prepare(source).unwrap().0
    }

    #[test]
    fn element_queries_match_independent_token_truth_table() {
        let mut errors = Vec::new();
        let mut checked = 0;
        for line in include_str!("../tests/fixtures/smarts_element_tokens.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let query = parse_smarts(fields[0]).unwrap();
            let mol = prepared(fields[1]);
            let rings = mol.ring_info();
            let matcher = PreparedMatcher::with_ring_info(&mol, &rings);
            let actual = (0..mol.atoms.len())
                .filter(|&i| matcher.matches(&query, i))
                .map(|i| i.to_string())
                .collect::<Vec<_>>()
                .join(",");
            if actual != fields[2] {
                errors.push(format!(
                    "{} on {}: {actual:?} vs {:?}",
                    fields[0], fields[1], fields[2]
                ));
            }
            checked += 1;
        }
        assert_eq!(checked, 638);
        assert!(
            errors.is_empty(),
            "{} token mismatches:\n{}",
            errors.len(),
            errors.join("\n")
        );
    }

    #[test]
    fn default_hydrogen_removal_matches_independent_retained_graphs() {
        let mut checked = 0;
        for line in include_str!("../tests/fixtures/qed_hydrogen_contexts.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let mol = prepared(fields[0]);
            let suppressed = hydrogen_suppressed(&mol);
            let graph = suppressed.as_ref().unwrap_or(&mol);
            let atoms = graph
                .atoms
                .iter()
                .map(|a| {
                    format!(
                        "{},{},{},{},{}",
                        a.symbol,
                        a.isotope.unwrap_or(0),
                        a.hydrogen,
                        a.charge,
                        u8::from(a.aromatic)
                    )
                })
                .collect::<Vec<_>>()
                .join(";");
            assert_eq!(atoms, fields[1], "{} retained atoms", fields[0]);
            let bonds = graph
                .bonds
                .iter()
                .map(|b| {
                    let order = match b.order {
                        BondOrder::Single => 1,
                        BondOrder::Double => 2,
                        BondOrder::Triple => 3,
                        BondOrder::Aromatic => 12,
                        BondOrder::Dative => 17,
                    };
                    (b.a, b.b, order)
                })
                .collect::<Vec<_>>();
            let reference = fields[2]
                .split(';')
                .filter(|b| !b.is_empty())
                .map(|b| {
                    let f: Vec<usize> = b.split(',').map(|v| v.parse().unwrap()).collect();
                    (f[0], f[1], f[2])
                })
                .collect::<Vec<_>>();
            let canonical = |mut edges: Vec<(usize, usize, usize)>| {
                for (a, b, order) in &mut edges {
                    // Ordinary edges are undirected; dative endpoints retain
                    // donor/acceptor direction. Bond array order is immaterial.
                    if *order != 17 && *a > *b {
                        std::mem::swap(a, b);
                    }
                }
                edges.sort_unstable();
                edges
            };
            assert_eq!(
                canonical(bonds),
                canonical(reference),
                "{} retained bonds",
                fields[0]
            );
            checked += 1;
        }
        assert_eq!(checked, 51);
    }

    #[test]
    fn registered_reference_properties_and_scores() {
        let mut differences = Vec::new();
        let mut checked = 0;
        for line in include_str!("../tests/fixtures/descriptor_qed.tsv").lines() {
            if line.starts_with('#') {
                continue;
            }
            let columns: Vec<_> = line.split('\t').collect();
            let source = columns[0];
            let mol = prepared(source);
            let suppressed = hydrogen_suppressed(&mol);
            let graph = suppressed.as_ref().unwrap_or(&mol);
            let ring = graph.ring_info();
            let matcher = PreparedMatcher::with_ring_info(graph, &ring);
            let expected_acceptors: Vec<usize> =
                columns[3].split(',').map(|v| v.parse().unwrap()).collect();
            let expected_alerts: Vec<bool> = columns[4].split(',').map(|v| v == "1").collect();
            if acceptor_counts(graph, &matcher) != expected_acceptors {
                differences.push(format!("{source}: QED acceptor query vector differs"));
            }
            if alert_flags(graph, &matcher) != expected_alerts {
                differences.push(format!("{source}: QED structural-alert vector differs"));
            }
            let expected: Vec<f64> = columns[2].split(',').map(|v| v.parse().unwrap()).collect();
            let actual = properties(&mol);
            for (index, (&a, &b)) in actual.iter().zip(&expected).enumerate() {
                let tolerance = if index == 0 { 1e-8 } else { 1e-9 };
                if !a.is_finite() || (a - b).abs() > tolerance {
                    differences.push(format!(
                        "{source}: property {index}: native={a}, reference={b}"
                    ));
                }
            }
            let reference: f64 = columns[1].parse().unwrap();
            let score = qed(&mol);
            if !score.is_finite() || (score - reference).abs() > 1e-9 {
                differences.push(format!(
                    "{source}: QED: native={score}, reference={reference}"
                ));
            }
            checked += 1;
        }
        assert_eq!(checked, 413);
        assert!(
            differences.is_empty(),
            "{} differences; first40:\n{}",
            differences.len(),
            differences
                .iter()
                .take(40)
                .cloned()
                .collect::<Vec<_>>()
                .join("\n")
        );
    }

    fn qed_of(s: &str) -> f64 {
        qed(&parse(s).expect("parse"))
    }

    #[test]
    fn all_constant_smarts_compile() {
        // The 11 acceptors and all 116 alerts must compile.
        assert_eq!(acceptors().len(), 11);
        assert_eq!(alerts().len(), 116, "every alert must compile");
        let _ = aliphatic_pattern();
        assert!(acceptors().iter().all(|p| p.atoms.len() == 1));
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
        assert_eq!(properties(&parse("CCCC").unwrap())[5], 1.);
    }

    #[test]
    fn rotb_excludes_tert_butyl_and_cf3() {
        // neopentane: all bonds touch a CH3 or the quaternary C → 0 rotatable
        assert_eq!(properties(&parse("CC(C)(C)C").unwrap())[5], 0.);
    }

    #[test]
    fn arom_benzene_is_one_naphthalene_two() {
        assert_eq!(properties(&parse("c1ccccc1").unwrap())[6], 1.);
        assert_eq!(properties(&parse("c1ccc2ccccc2c1").unwrap())[6], 2.);
    }

    #[test]
    fn arom_cyclohexane_is_zero() {
        // purely aliphatic ring → AROM 0 after deletion
        assert_eq!(properties(&parse("C1CCCCC1").unwrap())[6], 0.);
    }

    #[test]
    fn qed_in_unit_interval() {
        for s in [
            "CCO",
            "c1ccccc1",
            "CC(=O)Oc1ccccc1C(=O)O",
            "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
        ] {
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
        // These rounded original-paper values are only a historical smoke test;
        // the independent pinned fixture above supplies the precise checks.
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
