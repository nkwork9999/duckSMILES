//! One SMILES per line. --positive returns five channels; --audit adds ion receipts.
use ducksmiles_smiles::hr_features::{calculate, SPEC};
use std::io::{self, BufRead, Write};

fn main() -> io::Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args
        .iter()
        .any(|arg| arg != "--positive" && arg != "--audit")
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected --positive or --audit",
        ));
    }
    let positive = args.iter().any(|arg| arg == "--positive");
    let audit = args.iter().any(|arg| arg == "--audit");
    let mut out = io::BufWriter::new(io::stdout().lock());
    for line in io::stdin().lock().lines() {
        match calculate(line?.trim(), 2048) {
            Ok(features) => {
                write!(out, "{{\"spec\":{SPEC:?},\"channels\":[")?;
                for (i, channel) in features
                    .channels
                    .iter()
                    .take(if positive { 5 } else { 10 })
                    .enumerate()
                {
                    if i != 0 {
                        write!(out, ",")?;
                    }
                    write!(out, "[")?;
                    for (j, (bin, value)) in channel.iter().enumerate() {
                        if j != 0 {
                            write!(out, ",")?;
                        }
                        write!(out, "[{bin},{value}]")?;
                    }
                    write!(out, "]")?;
                }
                write!(
                    out,
                    "],\"eligible_fragments\":{},\"unsupported_boundary_fragments\":{}",
                    features.eligible_fragments, features.unsupported_boundary_fragments
                )?;
                if audit {
                    write!(out, ",\"ions\":[")?;
                    for (i, ion) in features.ions.iter().enumerate() {
                        if i != 0 {
                            write!(out, ",")?;
                        }
                        let fragment = &features.fragments[ion.fragment];
                        write!(out, "{{\"atoms\":{:?},\"original_formula\":{:?},\"formula\":{:?},\"polarity\":{},\"hydrogen_shift\":{},\"mz\":{:.12},\"cuts\":{}}}",
                               fragment.atoms, fragment.formula, ion.formula, ion.polarity, ion.hydrogen_shift, ion.mz, fragment.cuts)?;
                    }
                    write!(out, "]")?;
                }
                writeln!(out, "}}")?;
            }
            Err(error) => writeln!(out, "{{\"error\":{error:?}}}")?,
        }
    }
    out.flush()
}
