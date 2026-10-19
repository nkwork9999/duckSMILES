//! Reproducible engine-only baseline; not an iPhone or RDKit benchmark.
//! cargo run --release -p ducksmiles_smiles --example workbench_benchmark
use ducksmiles_smiles::{self as ds, workbench};
use std::{hint::black_box, time::Instant};

const PANEL: [&str; 8] = [
    "CCO",
    "CC(=O)Oc1ccccc1C(=O)O",
    "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
    "c1ccccc1",
    "[Na+].[Cl-]",
    "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
    "CC(=O)Nc1ccc(O)cc1",
    "N[C@@H](C)C(=O)O",
];
fn scalar(s: &str) -> ([f64; 12], String) {
    let p = s.as_ptr();
    let n = s.len();
    let mut formula = [0u8; 256];
    let len = ds::ds_mol_formula(p, n, formula.as_mut_ptr(), formula.len());
    assert!(len > 0 && len as usize <= formula.len());
    (
        [
            ds::ds_mol_weight(p, n),
            ds::ds_mol_exact_mass(p, n),
            ds::ds_logp_crippen(p, n),
            ds::ds_tpsa(p, n),
            ds::ds_num_h_donors(p, n) as f64,
            ds::ds_num_h_acceptors(p, n) as f64,
            ds::ds_mol_num_atoms(p, n) as f64,
            ds::ds_mol_num_bonds(p, n) as f64,
            ds::ds_mol_formal_charge(p, n) as f64,
            ds::ds_mol_num_fragments(p, n) as f64,
            ds::ds_num_rotatable_bonds(p, n) as f64,
            ds::ds_fraction_csp3(p, n),
        ],
        String::from_utf8(formula[..len as usize].to_vec()).unwrap(),
    )
}
fn fused(s: &str) -> ([f64; 12], String) {
    let d = workbench::analyze(s).unwrap();
    (
        [
            d.mw,
            d.exact_mass,
            d.logp,
            d.tpsa,
            d.hbd as f64,
            d.hba as f64,
            d.heavy_atoms as f64,
            d.bonds as f64,
            d.formal_charge as f64,
            d.fragments as f64,
            d.rotatable_bonds as f64,
            d.fraction_csp3,
        ],
        d.formula,
    )
}
fn timed(mut f: impl FnMut()) -> f64 {
    let t = Instant::now();
    f();
    t.elapsed().as_secs_f64() * 1000.0
}
fn main() {
    // A parity failure aborts before any speed result is reported.
    for s in PANEL {
        let (a, fa) = scalar(s);
        let (b, fb) = fused(s);
        assert_eq!(fa, fb, "formula: {s}");
        for (i, (x, y)) in a.iter().zip(b).enumerate() {
            assert!(
                x.is_finite() && y.is_finite() && (x - y).abs() < 1e-9,
                "{s} descriptor {i}: {x} vs {y}"
            );
        }
    }
    let inputs: Vec<_> = PANEL.into_iter().cycle().take(1_000).collect();
    let batches: Vec<_> = inputs.chunks(64).map(|c| c.join("\n")).collect();
    let mut output = vec![0u8; 64 * 1024];
    // Warm both paths and include JSON serialization in the native ABI timing.
    for s in &inputs {
        black_box(scalar(s));
        black_box(fused(s));
    }
    let mut scalar_ms = Vec::new();
    let mut fused_ms = Vec::new();
    let mut batch_ms = Vec::new();
    for trial in 0..5 {
        let measure_scalar = || {
            timed(|| {
                for s in &inputs {
                    black_box(scalar(black_box(s)));
                }
            })
        };
        let measure_fused = || {
            timed(|| {
                for s in &inputs {
                    black_box(fused(black_box(s)));
                }
            })
        };
        if trial % 2 == 0 {
            scalar_ms.push(measure_scalar());
            fused_ms.push(measure_fused());
        } else {
            fused_ms.push(measure_fused());
            scalar_ms.push(measure_scalar());
        }
        batch_ms.push(timed(|| {
            for input in &batches {
                // SAFETY: disjoint, owned input and writable output buffers.
                let n = unsafe {
                    workbench::ds_workbench_analyze(
                        input.as_ptr(),
                        input.len(),
                        output.as_mut_ptr(),
                        output.len(),
                    )
                };
                assert!(n > 0 && n as usize <= output.len());
                let json = std::str::from_utf8(&output[..n as usize]).unwrap();
                assert!(json.starts_with("{\"schema\":1,\"rows\":["));
                black_box(json);
            }
        }));
    }
    println!("{{\"scope\":\"mac_engine_only\",\"arch\":\"{}\",\"rows\":1000,\"unique_molecules\":8,\"batch_size\":64,\"trials\":5,\"parity_molecules\":8,\"scalar_ms\":{:?},\"fused_ms\":{:?},\"batch_json_ms\":{:?}}}",std::env::consts::ARCH,scalar_ms,fused_ms,batch_ms);
}
