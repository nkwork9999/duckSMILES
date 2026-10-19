use ducksmiles_smiles::{descriptor_profile::SPECS, workbench::analyze};

#[test]
fn large_input_path257_matches_all_independent_descriptor_references() {
    let source = include_str!("fixtures/descriptor_path257.tsv");
    let lines: Vec<_> = source.lines().filter(|l| !l.starts_with('#')).collect();
    assert_eq!(lines.len(), 1);
    let f: Vec<_> = lines[0].split('\t').collect();
    assert_eq!(f.len(), SPECS.len() + 2);
    let actual = analyze(f[0]).expect("the 257-atom path must be supported");
    assert_eq!(actual.formula, f[1]);
    let mut differences = Vec::new();
    for ((spec, actual), expected) in SPECS.iter().zip(actual.descriptor_values).zip(&f[2..]) {
        let expected = if *expected == "null" {
            None
        } else {
            Some(expected.parse::<f64>().unwrap())
        };
        let same = match (actual, expected) {
            (None, None) => true,
            (Some(a), Some(e)) => {
                (a - e).abs() <= spec.tolerance.max(spec.relative_tolerance * e.abs())
            }
            _ => false,
        };
        if !same {
            differences.push(format!("{}: {:?} vs {:?}", spec.id, actual, expected));
        }
    }
    assert!(differences.is_empty(), "{}", differences.join("\n"));
}
