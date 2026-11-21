use ducksmiles_mol::parser::{parse_sdf, project_sdf_properties_json, sdf_properties_json};

const WATER: &str = "water\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 O   0  0\nM  END";

#[test]
fn multiproperty_projection_matches_full_reader_on_valid_records() {
    let wanted = format!("{WATER}\n>  <NAME>\nwater\nwell\n\n>  <EMPTY>\n");
    let records = format!("{wanted}\n\n$$$$\n{wanted}");
    let raw = format!("{wanted}\n\n>  <IGNORE>\nother\n\n$$$$\n{wanted}\n\n>  <EXTRA>\ntext");
    assert_eq!(
        project_sdf_properties_json(&raw.replace('\n', "\r\n"), &["NAME", "EMPTY"]),
        sdf_properties_json(&parse_sdf(&records))
    );
}

#[test]
fn tags_can_be_read_without_parsing_coordinates() {
    let without_coordinates = WATER.replace("0.0000", "   NaN");
    let source = format!("{without_coordinates}\n>  <ID>\nmy-id");
    assert!(parse_sdf(&source).is_empty());
    assert!(project_sdf_properties_json(&source, &["ID"]).contains("\"value\":\"my-id\""));
}
