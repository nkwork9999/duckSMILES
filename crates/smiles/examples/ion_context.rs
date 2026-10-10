//! Default input: one SMILES per line. --positive:24 channels. --audit:receipts.
//! --match input: polarity<TAB>precursor<TAB>mzCSV<TAB>intensityCSV<TAB>SMILES.
use ducksmiles_smiles::ion_context;
use std::io::{self, BufRead, Write};

fn number_list(text: &str) -> Result<Vec<f64>, String> {
    if text.is_empty() {
        return Ok(Vec::new());
    }
    text.split(',')
        .map(|x| x.parse::<f64>().map_err(|_| "invalid numeric array".into()))
        .collect()
}
fn match_line(line: &str) -> Result<[f64; 8], String> {
    let fields: Vec<_> = line.split('\t').collect();
    if fields.len() != 5 {
        return Err("expected five tab-separated fields".into());
    }
    let mode = fields[0].parse::<i8>().map_err(|_| "invalid polarity")?;
    let precursor = fields[1].parse::<f64>().map_err(|_| "invalid precursor")?;
    ion_context::match_smiles(
        fields[4],
        mode,
        precursor,
        &number_list(fields[2])?,
        &number_list(fields[3])?,
    )
}
fn main() -> io::Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args
        .iter()
        .any(|a| a != "--positive" && a != "--audit" && a != "--match")
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "expected --positive, --audit or --match",
        ));
    }
    let positive = args.iter().any(|a| a == "--positive");
    let audit = args.iter().any(|a| a == "--audit");
    let matching = args.iter().any(|a| a == "--match");
    let mut out = io::BufWriter::new(io::stdout().lock());
    for input in io::stdin().lock().lines() {
        let line = input?;
        if matching {
            match match_line(&line) {
                Ok(scores) => writeln!(
                    out,
                    "{{\"spec\":{:?},\"scores\":{:?}}}",
                    ion_context::MATCH_SPEC,
                    scores
                )?,
                Err(error) => writeln!(out, "{{\"error\":{error:?}}}")?,
            }
            continue;
        }
        match ion_context::calculate(line.trim(), 2048) {
            Err(error) => writeln!(out, "{{\"error\":{error:?}}}")?,
            Ok(x) => {
                write!(out, "{{\"spec\":{:?},\"channels\":[", ion_context::SPEC)?;
                for (i, ch) in x
                    .channels
                    .iter()
                    .take(if positive { 24 } else { 48 })
                    .enumerate()
                {
                    if i > 0 {
                        write!(out, ",")?;
                    }
                    write!(out, "[")?;
                    for (j, (b, v)) in ch.iter().enumerate() {
                        if j > 0 {
                            write!(out, ",")?;
                        }
                        write!(out, "[{b},{v}]")?;
                    }
                    write!(out, "]")?;
                }
                write!(
                    out,
                    "],\"eligible_fragments\":{},\"unsupported_boundary_fragments\":{}",
                    x.hr.eligible_fragments, x.hr.unsupported_boundary_fragments
                )?;
                if audit {
                    write!(out, ",\"ions\":[")?;
                    for (i, ion) in x.hr.ions.iter().enumerate() {
                        if i > 0 {
                            write!(out, ",")?;
                        }
                        let f = &x.hr.fragments[ion.fragment];
                        write!(out,"{{\"atoms\":{:?},\"formula\":{:?},\"polarity\":{},\"hydrogen_shift\":{},\"mz\":{:.12},\"cuts\":{},\"ring\":{}}}",f.atoms,ion.formula,ion.polarity,ion.hydrogen_shift,ion.mz,f.cuts,f.ring)?;
                    }
                    write!(out, "],\"paths\":[")?;
                    for (i, (a, b)) in x.paths.iter().enumerate() {
                        if i > 0 {
                            write!(out, ",")?;
                        }
                        write!(out, "[{a},{b}]")?;
                    }
                    write!(out, "]")?;
                }
                writeln!(out, "}}")?;
            }
        }
    }
    out.flush()
}
