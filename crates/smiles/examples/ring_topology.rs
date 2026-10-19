//! Raw-parser topology diagnostic, not a prepared-chemistry or accuracy oracle.
//! Input: name<TAB>SMILES. Output retains every row, including raw parse failures.
use ducksmiles_smiles::verify::parse;
use std::io::{self, BufRead};

fn main() {
    for line in io::stdin().lock().lines() {
        let line = line.expect("read source row");
        let (name, source) = line.split_once('\t').expect("name<TAB>SMILES");
        match parse(source) {
            Some(mol) => println!(
                "{}\t{}\t{:?}\t{:?}",
                name,
                mol.atoms.len(),
                mol.bonds.iter().map(|b| (b.a, b.b)).collect::<Vec<_>>(),
                mol.ring_info().rings
            ),
            None => println!("{name}\tunsupported\t[]\t[]"),
        }
    }
}
