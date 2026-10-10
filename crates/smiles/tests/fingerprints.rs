mod common;

use common::{fixture_bits, on_bits, tsv};
use ducksmiles_smiles::*;
use std::collections::BTreeSet;

fn maccs(s: &str) -> BTreeSet<usize> {
    let mut b = [0u8; 21];
    assert_eq!(
        ds_maccs_keys(s.as_ptr(), s.len(), b.as_mut_ptr(), b.len()),
        21,
        "{s}"
    );
    on_bits(&b, 1, 167)
}
fn morgan(s: &str) -> BTreeSet<usize> {
    let mut b = [0u8; 256];
    assert_eq!(
        ds_morgan_fp_bits(s.as_ptr(), s.len(), 2, 2048, b.as_mut_ptr(), b.len()),
        256,
        "{s}"
    );
    on_bits(&b, 0, 2048)
}
fn tanimoto(a: &BTreeSet<usize>, b: &BTreeSet<usize>) -> f64 {
    let i = a.intersection(b).count();
    let u = a.len() + b.len() - i;
    if u == 0 { 1.0 } else { i as f64 / u as f64 }
}

#[test]
fn maccs_matches_rdkit_exactly() {
    let rows = tsv("rdkit_maccs.tsv");
    let mut over = [0usize; 167];
    let mut under = [0usize; 167];
    let mut molecules = 0;
    for r in &rows {
        let got = maccs(&r["smiles"]);
        let want = fixture_bits(&r["on_bits"]);
        if got != want {
            molecules += 1;
            for &x in got.difference(&want) {
                over[x] += 1
            }
            for &x in want.difference(&got) {
                under[x] += 1
            }
        }
    }
    let details = (1..167)
        .filter(|&i| over[i] + under[i] > 0)
        .map(|i| format!("key {i}: over={}, under={}", over[i], under[i]))
        .collect::<Vec<_>>();
    assert_eq!(
        molecules,
        0,
        "{molecules}/{} molecules mismatch\n{}",
        rows.len(),
        details.join("\n")
    );
}

fn ranks(values: &[f64]) -> Vec<f64> {
    let mut order: Vec<usize> = (0..values.len()).collect();
    order.sort_by(|&a, &b| values[a].total_cmp(&values[b]));
    let mut out = vec![0.0; values.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i + 1;
        while j < order.len() && values[order[j]] == values[order[i]] {
            j += 1
        }
        let rank = (i + j - 1) as f64 / 2.0 + 1.0;
        for k in i..j {
            out[order[k]] = rank
        }
        i = j
    }
    out
}
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let ma = a.iter().sum::<f64>() / a.len() as f64;
    let mb = b.iter().sum::<f64>() / b.len() as f64;
    let mut xy = 0.0;
    let mut xx = 0.0;
    let mut yy = 0.0;
    for i in 0..a.len() {
        let x = a[i] - ma;
        let y = b[i] - mb;
        xy += x * y;
        xx += x * x;
        yy += y * y
    }
    xy / (xx * yy).sqrt()
}

#[test]
fn morgan_density_and_similarity_track_rdkit() {
    let rows = tsv("rdkit_morgan.tsv");
    let mut duck = Vec::new();
    let mut rdkit = Vec::new();
    let mut density_bad = Vec::new();
    for r in &rows {
        let d = morgan(&r["smiles"]);
        let q = fixture_bits(&r["on_bits"]);
        let low = (q.len() as f64 * 0.5).floor() as usize;
        let high = (q.len() * 2) + 3;
        if d.len() < low || d.len() > high {
            density_bad.push(format!(
                "{}: duck={}, rdkit={}",
                r["smiles"],
                d.len(),
                q.len()
            ));
        }
        duck.push(d);
        rdkit.push(q);
    }
    assert!(
        density_bad.is_empty(),
        "Morgan density outside sane band:\n{}",
        density_bad.join("\n")
    );
    let mut ds = Vec::new();
    let mut rs = Vec::new();
    for i in 0..rows.len() {
        for j in i + 1..rows.len() {
            ds.push(tanimoto(&duck[i], &duck[j]));
            rs.push(tanimoto(&rdkit[i], &rdkit[j]));
        }
    }
    let rho = pearson(&ranks(&ds), &ranks(&rs));
    assert!(
        rho >= 0.95,
        "Spearman pairwise Tanimoto rho={rho:.6}, required >= 0.95"
    );
}

#[test]
fn morgan_is_spelling_invariant() {
    let mut groups = std::collections::BTreeMap::<String, Vec<String>>::new();
    for r in tsv("smiles_variants.tsv") {
        groups
            .entry(r["id"].clone())
            .or_default()
            .push(r["smiles"].clone());
    }
    for (id, v) in groups {
        let first = morgan(&v[0]);
        for s in &v[1..] {
            assert_eq!(morgan(s), first, "id={id}: {} / {s}", v[0]);
        }
    }
}
