//! Streaming token encoder. Input: comma-separated u32 set per line.
//! Output: comma-separated MinHash coordinates. Arguments: dimensions [seed].
use ducksmiles_smiles::minhash::MinHash;
use std::io::{self, BufRead, Write};
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().collect();
    let dimensions = args.get(1).map(|s| s.parse()).transpose()?.unwrap_or(256);
    let seed = args.get(2).map(|s| s.parse()).transpose()?.unwrap_or(42);
    let encoder = MinHash::new(dimensions, seed)?;
    if args.get(3).map(String::as_str) == Some("--coefficients") {
        for (a, b) in encoder.coefficients() {
            println!("{a},{b}");
        }
        return Ok(());
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        let s = line?;
        let tokens: Vec<u32> = if s.trim().is_empty() {
            vec![]
        } else {
            s.trim()
                .split(',')
                .map(str::parse)
                .collect::<Result<_, _>>()?
        };
        let values: Vec<_> = encoder.encode(&tokens).iter().map(u32::to_string).collect();
        writeln!(out, "{}", values.join(","))?;
    }
    out.flush()?;
    Ok(())
}
