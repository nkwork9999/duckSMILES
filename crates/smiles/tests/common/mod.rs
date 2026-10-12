#![allow(dead_code)]

use std::collections::{BTreeSet, HashMap};
use std::fs;
use std::path::PathBuf;

pub type StringFn = extern "C" fn(*const u8, usize, *mut u8, usize) -> i32;

pub fn data_path(name: &str) -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests/data")
        .join(name)
}

pub fn tsv(name: &str) -> Vec<HashMap<String, String>> {
    let text = fs::read_to_string(data_path(name)).unwrap_or_else(|e| panic!("read {name}: {e}"));
    let mut lines = text.lines();
    let header: Vec<&str> = lines.next().expect("TSV header").split('\t').collect();
    lines
        .enumerate()
        .map(|(line_no, line)| {
            let fields: Vec<&str> = line.split('\t').collect();
            assert_eq!(
                fields.len(),
                header.len(),
                "{name}: malformed line {}",
                line_no + 2
            );
            header
                .iter()
                .zip(fields)
                .map(|(k, v)| ((*k).to_owned(), v.to_owned()))
                .collect()
        })
        .collect()
}

pub fn ffi_string(f: StringFn, smiles: &str) -> Result<String, i32> {
    let mut buf = vec![0_u8; 16_384];
    let n = f(smiles.as_ptr(), smiles.len(), buf.as_mut_ptr(), buf.len());
    if n < 0 {
        return Err(n);
    }
    let n = n as usize;
    if n > buf.len() {
        return Err(-2);
    }
    Ok(String::from_utf8(buf[..n].to_vec()).expect("FFI returned non-UTF-8"))
}

pub fn on_bits(buf: &[u8], first: usize, end: usize) -> BTreeSet<usize> {
    (first..end)
        .filter(|&n| (buf[n / 8] >> (n % 8)) & 1 == 1)
        .collect()
}

pub fn fixture_bits(value: &str) -> BTreeSet<usize> {
    if value.is_empty() {
        return BTreeSet::new();
    }
    value
        .split(',')
        .map(|x| x.parse::<usize>().expect("bit index"))
        .collect()
}

pub fn mismatch_report(name: &str, mismatches: &[String], total: usize) -> String {
    let shown = mismatches
        .iter()
        .take(20)
        .cloned()
        .collect::<Vec<_>>()
        .join("\n");
    format!(
        "{name}: {} / {total} mismatches ({:.2}%)\n{shown}{}",
        mismatches.len(),
        100.0 * mismatches.len() as f64 / total.max(1) as f64,
        if mismatches.len() > 20 {
            "\n... (report capped at 20 lines)"
        } else {
            ""
        }
    )
}
