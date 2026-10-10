//! One SMILES per line; one JSON array of six sparse count vectors per line.
//! Errors produce null, preserving row alignment. No dependencies or C ABI changes.
use ducksmiles_smiles::morgan_experimental::{fingerprints, fingerprints_v2};
use ducksmiles_smiles::verify::{morgan_bits, parse};
use std::io::{self, BufRead, Write};
fn main() -> io::Result<()> {
    let mode = std::env::args().nth(1).unwrap_or_default();
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let smi = line?;
        let result = if mode == "--v2" {
            fingerprints_v2(smi.trim(), 2048)
        } else if mode == "--legacy" {
            parse(smi.trim())
                .ok_or("parse failed".to_string())
                .map(|m| {
                    vec![morgan_bits(&m, 2, 2048)
                        .iter()
                        .enumerate()
                        .flat_map(|(i, b)| {
                            (0..8).filter_map(move |j| {
                                if b & (1 << j) != 0 {
                                    Some((i * 8 + j, 1))
                                } else {
                                    None
                                }
                            })
                        })
                        .collect()]
                })
        } else {
            fingerprints(smi.trim(), 2048)
        };
        match result {
            Ok(fps) => {
                let rows: Vec<_> = fps
                    .iter()
                    .map(|fp| {
                        format!(
                            "[{}]",
                            fp.iter()
                                .map(|(i, c)| format!("[{i},{c}]"))
                                .collect::<Vec<_>>()
                                .join(",")
                        )
                    })
                    .collect();
                writeln!(out, "[{}]", rows.join(","))?;
            }
            Err(e) => {
                eprintln!("{e}");
                writeln!(out, "null")?;
            }
        }
    }
    out.flush()
}
