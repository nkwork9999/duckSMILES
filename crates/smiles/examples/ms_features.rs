//! SMILES input; dependency-free JSON lines. Error objects retain row alignment.
use ducksmiles_smiles::ms_features::calculate;
use std::io::{self, BufRead, Write};
fn main() -> io::Result<()> {
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        match calculate(line?.trim(), 2048) {
            Ok(a) => {
                write!(out, "{{\"channels\":[")?;
                for (i, c) in a.channels.iter().enumerate() {
                    if i > 0 {
                        write!(out, ",")?;
                    }
                    write!(out, "[")?;
                    for (j, (b, v)) in c.iter().enumerate() {
                        if j > 0 {
                            write!(out, ",")?;
                        }
                        write!(out, "[{b},{v}]")?;
                    }
                    write!(out, "]")?;
                }
                write!(
                    out,
                    "],\"properties\":{:?},\"formula\":{:?},\"mass\":{},\"fragments\":[",
                    a.properties, a.formula, a.mass
                )?;
                for (i, f) in a.fragments.iter().enumerate() {
                    if i > 0 {
                        write!(out, ",")?;
                    }
                    write!(
                        out,
                        "[{:.8},{},{},{},{},{},{},{}]",
                        f.mass,
                        f.formula[0],
                        f.cuts,
                        u8::from(f.ring),
                        f.basic,
                        f.acidic,
                        f.environment,
                        f.atoms.len()
                    )?;
                }
                writeln!(out, "]}}")?;
            }
            Err(e) => writeln!(out, "{{\"error\":{e:?}}}")?,
        }
    }
    out.flush()
}
