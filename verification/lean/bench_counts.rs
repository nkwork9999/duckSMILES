// Same prepared element labels and query order as Bench.lean.
use std::{env, fs, time::Instant};

fn main() {
    let args: Vec<_> = env::args().collect();
    let text = fs::read_to_string(&args[1]).unwrap();
    let n: usize = args[2].parse().unwrap();
    let rows: Vec<Vec<&str>> = text
        .lines()
        .map(|s| s.split_whitespace().collect())
        .collect();
    let probes = ["C", "N", "O", "Cl", "F"];
    let start = Instant::now();
    let mut checksum = 0;
    for i in 0..n {
        let key = probes[i % probes.len()];
        for atoms in &rows {
            checksum += atoms.iter().filter(|atom| **atom == key).count();
        }
    }
    println!("{},{}", start.elapsed().as_nanos(), checksum);
}
