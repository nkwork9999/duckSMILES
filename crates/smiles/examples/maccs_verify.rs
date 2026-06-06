//! Emit MACCS on-bits (1..=166) for each SMILES read from stdin, one per line,
//! as `SMILES\tbit,bit,bit`. Used to cross-check against Python RDKit.
use std::io::{self, BufRead, Write};
use ducksmiles_smiles::verify::{maccs_bits, on_bits, parse};

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for line in stdin.lock().lines() {
        let smi = match line {
            Ok(s) => s,
            Err(_) => break,
        };
        let smi = smi.trim();
        if smi.is_empty() {
            continue;
        }
        match parse(smi) {
            Some(mol) => {
                let bits = on_bits(&maccs_bits(&mol));
                let joined: Vec<String> = bits.iter().map(|b| b.to_string()).collect();
                writeln!(out, "{}\t{}", smi, joined.join(",")).ok();
            }
            None => {
                writeln!(out, "{}\tPARSE_FAIL", smi).ok();
            }
        }
    }
}
