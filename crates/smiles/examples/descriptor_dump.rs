//! Dump scalar descriptors for each SMILES read from stdin, one per line, as
//! TSV. Used to cross-check against Python RDKit.
use std::io::{self, BufRead, Write};

use ducksmiles_smiles::*;

fn s(smi: &str) -> (*const u8, usize) {
    (smi.as_ptr(), smi.len())
}

fn str_fn(f: unsafe extern "C" fn(*const u8, usize, *mut u8, usize) -> i32, smi: &str) -> String {
    let mut buf = vec![0u8; 8192];
    let (p, l) = s(smi);
    let n = unsafe { f(p, l, buf.as_mut_ptr(), buf.len()) };
    if n < 0 {
        return "ERR".to_string();
    }
    let n = n as usize;
    if n > buf.len() {
        return "TOOBIG".to_string();
    }
    String::from_utf8_lossy(&buf[..n]).to_string()
}

const HEADER: &[&str] = &[
    "smiles",
    "num_atoms",
    "num_bonds",
    "formula",
    "mol_weight",
    "exact_mass",
    "tpsa",
    "logp",
    "mol_mr",
    "qed",
    "hbd",
    "hba",
    "rotb",
    "heteroatoms",
    "fcsp3",
    "ring_count",
    "arom_rings",
    "aliph_rings",
    "satur_rings",
    "arom_carbo",
    "arom_hetero",
    "satur_carbo",
    "satur_hetero",
    "aliph_carbo",
    "aliph_hetero",
    "lipinski_viol",
    "canonical",
    "murcko",
];

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    writeln!(out, "{}", HEADER.join("\t")).ok();
    for line in stdin.lock().lines() {
        let smi = match line {
            Ok(v) => v,
            Err(_) => break,
        };
        let smi = smi.trim().to_string();
        if smi.is_empty() {
            continue;
        }
        let (p, l) = s(&smi);
        let mut cols: Vec<String> = vec![smi.clone()];
        cols.push(ds_mol_num_atoms(p, l).to_string());
        cols.push(ds_mol_num_bonds(p, l).to_string());
        cols.push(str_fn(ds_mol_formula, &smi));
        {
            for v in [
                ds_mol_weight(p, l),
                ds_mol_exact_mass(p, l),
                ds_tpsa(p, l),
                ds_logp_crippen(p, l),
                ds_mol_mr(p, l),
                ds_qed(p, l),
            ] {
                cols.push(format!("{:.4}", v));
            }
            for v in [
                ds_num_h_donors(p, l),
                ds_num_h_acceptors(p, l),
                ds_num_rotatable_bonds(p, l),
                ds_num_heteroatoms(p, l),
            ] {
                cols.push(v.to_string());
            }
            cols.push(format!("{:.4}", ds_fraction_csp3(p, l)));
            for v in [
                ds_ring_count(p, l),
                ds_num_aromatic_rings(p, l),
                ds_num_aliphatic_rings(p, l),
                ds_num_saturated_rings(p, l),
                ds_num_aromatic_carbocycles(p, l),
                ds_num_aromatic_heterocycles(p, l),
                ds_num_saturated_carbocycles(p, l),
                ds_num_saturated_heterocycles(p, l),
                ds_num_aliphatic_carbocycles(p, l),
                ds_num_aliphatic_heterocycles(p, l),
                ds_lipinski_violations(p, l),
            ] {
                cols.push(v.to_string());
            }
        }
        cols.push(str_fn(ds_canonical_smiles, &smi));
        cols.push(str_fn(ds_murcko_scaffold, &smi));
        writeln!(out, "{}", cols.join("\t")).ok();
    }
}
