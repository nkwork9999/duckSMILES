use ducksmiles_smiles::paper_features::calculate;
use std::io::{self, BufRead, Write};
fn main() -> io::Result<()> {
    let channel = std::env::args()
        .nth(1)
        .map(|s| s.parse::<usize>())
        .transpose()
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidInput, e))?;
    if channel.is_some_and(|i| i >= 50) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "channel must be0..49",
        ));
    }
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        match calculate(line?.trim(), 2048) {
            Ok(x) => {
                let x = if let Some(i) = channel {
                    vec![x[i].clone()]
                } else {
                    x
                };
                write!(out, "[")?;
                for (i, c) in x.iter().enumerate() {
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
                writeln!(out, "]")?;
            }
            Err(e) => writeln!(out, "{{\"error\":{e:?}}}")?,
        }
    }
    out.flush()
}
