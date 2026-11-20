use ducksmiles_analytics::execute;
fn number(j: &str, key: &str) -> f64 {
    let p = j.find(&format!("\"{key}\":")).unwrap() + key.len() + 3;
    j[p..]
        .split([',', '}', ']'])
        .next()
        .unwrap()
        .parse()
        .unwrap()
}
#[test]
fn ideal_perovskite_and_charge_are_independent() {
    let a = 2_f64.sqrt() * 2. - 1.;
    let out = execute(
        "perovskite_descriptors",
        &[vec![a, 1., 1., 2., 4., -2.]],
        &[],
    )
    .unwrap();
    assert!((number(&out, "goldschmidt_t") - 1.).abs() < 1e-14);
    assert_eq!(number(&out, "octahedral_mu"), 1.);
    assert!(out.contains("\"charge_balanced\":true"));
    let out = execute(
        "perovskite_descriptors",
        &[vec![a, 1., 1., 1., 4., -2.]],
        &[],
    )
    .unwrap();
    assert_eq!(number(&out, "charge_sum"), -1.);
    assert!(out.contains("\"charge_balanced\":false"));
}
#[test]
fn viscosity_fit_roundtrip_and_training_residuals() {
    let data: Vec<Vec<f64>> = [280., 300., 330., 360.]
        .iter()
        .map(|t| vec![*t, (-12_f64 + 1500. / t).exp()])
        .collect();
    let out = execute("viscosity_fit", &data, &[]).unwrap();
    assert!((number(&out, "activation_energy_J_mol") - 1500. * 8.31446261815324).abs() < 1e-8);
    assert!(number(&out, "rmse_log") < 1e-12);
    assert!(out.contains("\"temperature_range_K\":[280,360]"));
    let out = execute("viscosity_arrhenius", &[vec![300.]], &[-12., 1500.]).unwrap();
    assert!(out.contains(&(-7_f64).exp().to_string()));
}
#[test]
fn mixtures_obey_endpoints_order_and_interaction() {
    for rows in [
        vec![vec![0.5, 0.001], vec![0.5, 0.004]],
        vec![vec![0.5, 0.004], vec![0.5, 0.001]],
    ] {
        let ideal = execute("viscosity_mix", &rows, &[]).unwrap();
        assert!((number(&ideal, "viscosity_Pa_s") - 0.002).abs() < 1e-15);
        let interaction = execute("viscosity_mix", &rows, &[2.]).unwrap();
        assert!((number(&interaction, "viscosity_Pa_s") - 0.002 * 0.5_f64.exp()).abs() < 1e-15);
    }
    let out = execute(
        "viscosity_mix",
        &[vec![1., 0.001], vec![0., 0.004]],
        &[100.],
    )
    .unwrap();
    assert!((number(&out, "viscosity_Pa_s") - 0.001).abs() < 1e-15);
    let out = execute(
        "mixture_descriptors",
        &[vec![0.25, 2.], vec![0.75, 6.]],
        &[],
    )
    .unwrap();
    assert!(out.contains("\"means\":[5]"));
    assert!(out.contains(&3_f64.sqrt().to_string()));
}
#[test]
fn invalid_inputs_and_unrepresentable_results_are_errors() {
    for (op, d, p) in [
        (
            "perovskite_descriptors",
            vec![vec![1., 0., 1., 1., 2., -1.]],
            vec![],
        ),
        (
            "perovskite_descriptors",
            vec![vec![1., 1., 1., -1., 2., -1.]],
            vec![],
        ),
        ("viscosity_arrhenius", vec![vec![0.]], vec![1., 2.]),
        ("viscosity_arrhenius", vec![vec![300.]], vec![1000., 0.]),
        ("viscosity_arrhenius", vec![vec![300.]], vec![-1000., 0.]),
        ("viscosity_fit", vec![vec![300., 1.]; 3], vec![]),
        (
            "viscosity_fit",
            vec![vec![300., 1.], vec![310., 0.], vec![320., 2.]],
            vec![],
        ),
        ("viscosity_mix", vec![vec![0.4, 1.], vec![0.4, 2.]], vec![]),
        ("viscosity_mix", vec![vec![1., -1.]], vec![]),
        ("viscosity_mix", vec![vec![1., 1.]], vec![0.]),
        ("mixture_descriptors", vec![vec![1., f64::INFINITY]], vec![]),
        (
            "mixture_descriptors",
            vec![vec![0.5, 1.], vec![0.5]],
            vec![],
        ),
    ] {
        assert!(execute(op, &d, &p).is_err(), "{op}: {d:?}");
    }
}
