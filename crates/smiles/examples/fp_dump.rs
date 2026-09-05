//! Emit MACCS on-bits and Morgan(r=2,2048) on-bits for each SMILES from stdin,
//! as `SMILES\tmaccs_bits\tmorgan_bits`. Cross-check against Python RDKit.
use std::io::{self, BufRead, Write};

use ducksmiles_smiles::{ds_maccs_keys, ds_morgan_fp_bits};

fn on_bits(buf: &[u8], n: usize, start: usize) -> Vec<usize> {
    let mut v = Vec::new();
    for b in start..n {
        if buf[b / 8] >> (b % 8) & 1 == 1 {
            v.push(b);
        }
    }
    v
}

fn join(v: &[usize]) -> String {
    v.iter()
        .map(|x| x.to_string())
        .collect::<Vec<_>>()
        .join(",")
}

fn main() {
    let stdin = io::stdin();
    let stdout = io::stdout();
    let mut out = stdout.lock();
    for line in stdin.lock().lines() {
        let smi = match line {
            Ok(v) => v,
            Err(_) => break,
        };
        let smi = smi.trim().to_string();
        if smi.is_empty() {
            continue;
        }
        let (p, l) = (smi.as_ptr(), smi.len());
        let mut mac = vec![0u8; 32];
        let nm = ds_maccs_keys(p, l, mac.as_mut_ptr(), mac.len());
        let mut mor = vec![0u8; 256];
        let nn = ds_morgan_fp_bits(p, l, 2, 2048, mor.as_mut_ptr(), mor.len());
        if nm < 0 || nn < 0 {
            writeln!(out, "{}\tERR\tERR", smi).ok();
            continue;
        }
        writeln!(
            out,
            "{}\t{}\t{}",
            smi,
            join(&on_bits(&mac, 167, 1)),
            join(&on_bits(&mor, 2048, 0))
        )
        .ok();
    }
}
