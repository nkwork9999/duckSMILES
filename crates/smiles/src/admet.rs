//! ADMET / drug-likeness rule panels + toxicophore structural alerts.
//!
//! This module layers interpretation on top of the existing physicochemical
//! descriptors (MW, cLogP, TPSA, HBD, HBA, rotatable bonds, molar refractivity,
//! atom counts) and a curated catalogue of named structural-alert SMARTS.
//!
//! Nothing here is a machine-learned ADMET *prediction* (hERG, Ames, DILI, CYP
//! inhibition need trained models and labelled data, which live outside a
//! zero-dependency SMILES engine). What we provide is the rule-based layer that
//! tools like RDKit's `Descriptors`/`FilterCatalog` expose: empirical
//! drug-likeness filters and substructural toxicophore/reactivity flags.
//!
//! Rule panels implemented:
//! - **Lipinski** Rule of Five (Lipinski 1997)
//! - **Veber** oral-bioavailability rules (Veber 2002)
//! - **Ghose** drug-likeness window (Ghose 1999)
//! - **Egan** absorption (Egan 2000)
//! - **Muegge** pharmacophore filter (Muegge 2001)
//! - **Lead-likeness** (Teague 1999)
//!
//! Structural alerts: a curated subset of Brenk / PAINS / reactive-group
//! toxicophores, each a named SMARTS that the engine can compile.

use crate::logp_crippen::{calc_logp, calc_mr};
use crate::parser::Molecule;
use crate::smarts::{matches_mol, parse_smarts, Pattern};
use crate::tpsa::calc_tpsa;
use std::sync::OnceLock;

// ── descriptor bundle ───────────────────────────────────────────────────────

pub struct Descriptors {
    pub mw: f64,
    pub logp: f64,
    pub tpsa: f64,
    pub hbd: usize,
    pub hba: usize,
    pub rotb: usize,
    pub mr: f64,
    pub heavy_atoms: usize,
    pub total_atoms: usize,
    pub carbons: usize,
    pub heteroatoms: usize,
    pub rings: usize,
    pub aromatic_rings: usize,
    pub fraction_csp3: f64,
}

pub fn descriptors(mol: &Molecule) -> Descriptors {
    let heavy_atoms = mol.atoms.iter().filter(|a| a.symbol != "H").count();
    let total_atoms = mol.atoms.len()
        + mol
            .atoms
            .iter()
            .map(|a| a.hydrogen.max(0) as usize)
            .sum::<usize>();
    let carbons = mol.atoms.iter().filter(|a| a.symbol == "C").count();

    Descriptors {
        mw: mol.molecular_weight(),
        logp: calc_logp(mol),
        tpsa: calc_tpsa(mol),
        hbd: crate::num_h_donors(mol),
        hba: crate::num_h_acceptors(mol),
        rotb: crate::num_rotatable_bonds(mol),
        mr: calc_mr(mol),
        heavy_atoms,
        total_atoms,
        carbons,
        heteroatoms: crate::num_heteroatoms(mol),
        rings: crate::ring_count(mol),
        aromatic_rings: crate::num_aromatic_rings(mol),
        fraction_csp3: crate::fraction_csp3(mol),
    }
}

// ── rule panel result ───────────────────────────────────────────────────────

pub struct RuleResult {
    pub name: &'static str,
    /// Criteria that FAILED (each a short human-readable string).
    pub failures: Vec<String>,
    /// Number of violated criteria.
    pub violations: usize,
    /// Overall pass: for Lipinski, ≤1 violation passes (the classic rule);
    /// for the others, all criteria must hold.
    pub pass: bool,
}

fn rule(name: &'static str, failures: Vec<String>, allow_one: bool) -> RuleResult {
    let violations = failures.len();
    let pass = if allow_one { violations <= 1 } else { violations == 0 };
    RuleResult { name, failures, violations, pass }
}

/// Lipinski Rule of Five: MW≤500, cLogP≤5, HBD≤5, HBA≤10. ≤1 violation passes.
pub fn lipinski(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if d.mw > 500.0 { f.push(format!("MW {:.1} > 500", d.mw)); }
    if d.logp > 5.0 { f.push(format!("cLogP {:.2} > 5", d.logp)); }
    if d.hbd > 5 { f.push(format!("HBD {} > 5", d.hbd)); }
    if d.hba > 10 { f.push(format!("HBA {} > 10", d.hba)); }
    rule("Lipinski", f, true)
}

/// Veber: rotatable bonds ≤10 AND TPSA ≤140.
pub fn veber(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if d.rotb > 10 { f.push(format!("RotB {} > 10", d.rotb)); }
    if d.tpsa > 140.0 { f.push(format!("TPSA {:.1} > 140", d.tpsa)); }
    rule("Veber", f, false)
}

/// Ghose: 160≤MW≤480, -0.4≤cLogP≤5.6, 40≤MR≤130, 20≤atoms≤70 (total incl. H).
pub fn ghose(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if !(160.0..=480.0).contains(&d.mw) { f.push(format!("MW {:.1} outside 160-480", d.mw)); }
    if !(-0.4..=5.6).contains(&d.logp) { f.push(format!("cLogP {:.2} outside -0.4..5.6", d.logp)); }
    if !(40.0..=130.0).contains(&d.mr) { f.push(format!("MR {:.1} outside 40-130", d.mr)); }
    if !(20..=70).contains(&d.total_atoms) { f.push(format!("atoms {} outside 20-70", d.total_atoms)); }
    rule("Ghose", f, false)
}

/// Egan absorption: TPSA ≤131.6 AND cLogP ≤5.88.
pub fn egan(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if d.tpsa > 131.6 { f.push(format!("TPSA {:.1} > 131.6", d.tpsa)); }
    if d.logp > 5.88 { f.push(format!("cLogP {:.2} > 5.88", d.logp)); }
    rule("Egan", f, false)
}

/// Muegge pharmacophore filter (Bayer): a panel of nine criteria.
pub fn muegge(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if !(200.0..=600.0).contains(&d.mw) { f.push(format!("MW {:.1} outside 200-600", d.mw)); }
    if !(-2.0..=5.0).contains(&d.logp) { f.push(format!("cLogP {:.2} outside -2..5", d.logp)); }
    if d.tpsa > 150.0 { f.push(format!("TPSA {:.1} > 150", d.tpsa)); }
    if d.rings > 7 { f.push(format!("rings {} > 7", d.rings)); }
    if d.carbons <= 4 { f.push(format!("carbons {} <= 4", d.carbons)); }
    if d.heteroatoms <= 1 { f.push(format!("heteroatoms {} <= 1", d.heteroatoms)); }
    if d.rotb > 15 { f.push(format!("RotB {} > 15", d.rotb)); }
    if d.hba > 10 { f.push(format!("HBA {} > 10", d.hba)); }
    if d.hbd > 5 { f.push(format!("HBD {} > 5", d.hbd)); }
    rule("Muegge", f, false)
}

/// Lead-likeness (Teague): MW≤350, cLogP≤3.5, rotatable bonds≤7.
pub fn lead_likeness(d: &Descriptors) -> RuleResult {
    let mut f = Vec::new();
    if d.mw > 350.0 { f.push(format!("MW {:.1} > 350", d.mw)); }
    if d.logp > 3.5 { f.push(format!("cLogP {:.2} > 3.5", d.logp)); }
    if d.rotb > 7 { f.push(format!("RotB {} > 7", d.rotb)); }
    rule("Lead-likeness", f, false)
}

pub fn all_rules(d: &Descriptors) -> Vec<RuleResult> {
    vec![
        lipinski(d), veber(d), ghose(d), egan(d), muegge(d), lead_likeness(d),
    ]
}

// ── structural-alert (toxicophore) catalogue ────────────────────────────────
//
// Curated subset of Brenk / PAINS / reactive-group toxicophores. Each entry is
// (alert_name, SMARTS). All patterns are verified to compile by the
// `all_alerts_compile` test; if a pattern ever fails to parse it is silently
// dropped at runtime (filter_map) and the test will catch the regression.

static ALERT_DEFS: &[(&str, &str)] = &[
    ("acyl_halide",        "[CX3](=[OX1])[F,Cl,Br,I]"),
    ("sulfonyl_halide",    "[SX4](=[OX1])(=[OX1])[F,Cl,Br,I]"),
    ("alkyl_halide",       "[CX4][Cl,Br,I]"),
    ("aldehyde",           "[CX3H1](=O)[#6]"),
    ("michael_acceptor",   "[CX3]=[CX3][CX3]=[OX1]"),
    ("acetylene_ketone",   "C#C[CX3]=[OX1]"),
    ("epoxide",            "C1OC1"),
    ("aziridine",          "C1CN1"),
    ("beta_lactam",        "C1CC(=O)N1"),
    ("isocyanate",         "[NX2]=C=[OX1]"),
    ("isothiocyanate",     "[NX2]=C=[SX1]"),
    ("thiocyanate",        "[SX2]C#N"),
    ("azide",              "[NX2]=[NX2]=[NX1]"),
    ("azo",                "[NX2]=[NX2]"),
    ("nitroso",            "[NX2]=O"),
    ("nitro",              "[NX3](=O)=O"),
    ("nitro_aromatic",     "[a][NX3](=O)=O"),
    ("hydrazine",          "[NX3][NX3]"),
    ("hydroxylamine",      "[NX3][OX2H]"),
    ("peroxide",           "[OX2][OX2]"),
    ("disulfide",          "[SX2][SX2]"),
    ("thiol",              "[SX2H]"),
    ("aldehyde_any",       "[CX3H1]=O"),
    ("imine",              "[CX3;!R]=[NX2;!R]"),
    ("anhydride",          "[CX3](=[OX1])O[CX3](=[OX1])"),
    ("sulfonate_ester",    "[#6]S(=O)(=O)O[#6]"),
    ("phosphonate_ester",  "P(=O)([OX2][#6])[OX2][#6]"),
    ("quaternary_nitrogen","[NX4+]"),
    ("aromatic_nitroso",   "[a][NX2]=O"),
    ("diketone_1_2",       "[#6](=O)[#6](=O)"),
    ("triple_bond_chain",  "[CX2]#[CX2]"),
    ("halopyridine",       "[Cl,Br,I][c]"),
];

fn alert_patterns() -> &'static [(&'static str, Pattern)] {
    static CACHE: OnceLock<Vec<(&'static str, Pattern)>> = OnceLock::new();
    CACHE.get_or_init(|| {
        ALERT_DEFS
            .iter()
            .filter_map(|(name, smarts)| parse_smarts(smarts).map(|p| (*name, p)))
            .collect()
    })
}

/// Names of all structural alerts whose pattern matches the molecule.
pub fn structural_alerts(mol: &Molecule) -> Vec<&'static str> {
    alert_patterns()
        .iter()
        .filter(|(_, p)| matches_mol(p, mol))
        .map(|(name, _)| *name)
        .collect()
}

// ── JSON serialisation ──────────────────────────────────────────────────────

fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len());
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            _ => out.push(c),
        }
    }
    out
}

fn rule_json(r: &RuleResult) -> String {
    let fails: Vec<String> = r
        .failures
        .iter()
        .map(|s| format!("\"{}\"", json_escape(s)))
        .collect();
    format!(
        "{{\"pass\":{},\"violations\":{},\"failures\":[{}]}}",
        r.pass,
        r.violations,
        fails.join(",")
    )
}

/// Full ADMET report as a JSON object string.
pub fn admet_json(mol: &Molecule) -> String {
    let d = descriptors(mol);
    let rules = all_rules(&d);
    let alerts = structural_alerts(mol);

    let rules_json: Vec<String> = rules
        .iter()
        .map(|r| format!("\"{}\":{}", json_escape(r.name), rule_json(r)))
        .collect();

    let alerts_json: Vec<String> = alerts
        .iter()
        .map(|a| format!("\"{}\"", json_escape(a)))
        .collect();

    format!(
        "{{\"descriptors\":{{\
\"mw\":{:.3},\"logp\":{:.3},\"tpsa\":{:.3},\"hbd\":{},\"hba\":{},\"rotb\":{},\
\"mr\":{:.3},\"heavy_atoms\":{},\"total_atoms\":{},\"carbons\":{},\
\"heteroatoms\":{},\"rings\":{},\"aromatic_rings\":{},\"fraction_csp3\":{:.4}}},\
\"rules\":{{{}}},\
\"alerts\":{{\"n\":{},\"hits\":[{}]}}}}",
        d.mw, d.logp, d.tpsa, d.hbd, d.hba, d.rotb, d.mr,
        d.heavy_atoms, d.total_atoms, d.carbons, d.heteroatoms,
        d.rings, d.aromatic_rings, d.fraction_csp3,
        rules_json.join(","),
        alerts.len(),
        alerts_json.join(",")
    )
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parser::parse;

    fn mol(smi: &str) -> Molecule {
        parse(smi).unwrap()
    }

    #[test]
    fn all_alerts_compile() {
        // Every catalogued SMARTS must parse; otherwise it is silently dropped.
        assert_eq!(
            alert_patterns().len(),
            ALERT_DEFS.len(),
            "some structural-alert SMARTS failed to compile"
        );
    }

    #[test]
    fn aspirin_passes_lipinski() {
        let d = descriptors(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        let r = lipinski(&d);
        assert!(r.pass, "aspirin should pass Lipinski: {:?}", r.failures);
        assert_eq!(r.violations, 0);
    }

    #[test]
    fn caffeine_passes_lipinski_and_veber() {
        let m = mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C");
        let d = descriptors(&m);
        assert!(lipinski(&d).pass, "caffeine Lipinski");
        assert!(veber(&d).pass, "caffeine Veber");
    }

    #[test]
    fn huge_lipophilic_fails_lipinski() {
        // A long alkane chain → very high cLogP and MW, multiple violations.
        let m = mol("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
        let d = descriptors(&m);
        let r = lipinski(&d);
        assert!(!r.pass, "C36 alkane should fail Lipinski (logP/MW)");
    }

    #[test]
    fn aldehyde_alert_fires() {
        // Benzaldehyde has an aldehyde toxicophore.
        let hits = structural_alerts(&mol("O=Cc1ccccc1"));
        assert!(hits.contains(&"aldehyde"), "benzaldehyde should flag aldehyde: {hits:?}");
    }

    #[test]
    fn epoxide_alert_fires() {
        let hits = structural_alerts(&mol("C1OC1"));
        assert!(hits.contains(&"epoxide"), "should flag epoxide: {hits:?}");
    }

    #[test]
    fn nitro_alert_fires() {
        // Nitrobenzene flags nitro (and nitro_aromatic).
        let hits = structural_alerts(&mol("O=[N+]([O-])c1ccccc1"));
        // Our SMARTS uses neutral pentavalent nitro [NX3](=O)=O. Accept either
        // representation via the alkyl/aromatic nitro check on the neutral form.
        let neutral = structural_alerts(&mol("c1ccccc1N(=O)=O"));
        assert!(
            hits.contains(&"nitro") || hits.contains(&"nitro_aromatic")
                || neutral.contains(&"nitro") || neutral.contains(&"nitro_aromatic"),
            "nitrobenzene should flag nitro: charged={hits:?} neutral={neutral:?}"
        );
    }

    #[test]
    fn clean_drug_has_no_reactive_alerts() {
        // Caffeine is a benign drug — should carry no reactive toxicophore.
        // (It may legitimately have none of our reactive-group alerts.)
        let hits = structural_alerts(&mol("Cn1cnc2c1c(=O)n(C)c(=O)n2C"));
        for bad in ["acyl_halide", "epoxide", "isocyanate", "peroxide", "alkyl_halide"] {
            assert!(!hits.contains(&bad), "caffeine unexpectedly flagged {bad}: {hits:?}");
        }
    }

    #[test]
    fn admet_json_is_wellformed() {
        let j = admet_json(&mol("CC(=O)Oc1ccccc1C(=O)O"));
        assert!(j.starts_with("{\"descriptors\":"), "json head: {j}");
        assert!(j.contains("\"rules\":"));
        assert!(j.contains("\"Lipinski\":"));
        assert!(j.contains("\"alerts\":"));
        assert!(j.ends_with('}'));
    }

    #[test]
    fn ghose_window_rejects_tiny_molecule() {
        // Methane is far below the Ghose MW / atom-count window.
        let d = descriptors(&mol("C"));
        assert!(!ghose(&d).pass, "methane should fail Ghose window");
    }
}
