use crate::Atom;

/// Parse mmCIF _atom_site loop
pub fn parse_cif(text: &str, include_hetatm: bool) -> Vec<Atom> {
    let mut atoms = Vec::new();
    let lines: Vec<&str> = text.lines().collect();
    let mut i = 0;

    // Find _atom_site loop
    while i < lines.len() {
        if lines[i].trim() == "loop_" {
            i += 1;
            // Collect column names
            let mut columns: Vec<String> = Vec::new();
            while i < lines.len() && lines[i].trim().starts_with("_atom_site.") {
                columns.push(lines[i].trim().to_string());
                i += 1;
            }
            if columns.is_empty() {
                continue;
            }

            // Map column names to indices
            let col_idx = |name: &str| columns.iter().position(|c| c == name);
            let i_group = col_idx("_atom_site.group_PDB");
            let i_name = col_idx("_atom_site.label_atom_id").or_else(|| col_idx("_atom_site.auth_atom_id"));
            let i_resname = col_idx("_atom_site.label_comp_id").or_else(|| col_idx("_atom_site.auth_comp_id"));
            let i_chain = col_idx("_atom_site.label_asym_id").or_else(|| col_idx("_atom_site.auth_asym_id"));
            let i_resid = col_idx("_atom_site.label_seq_id").or_else(|| col_idx("_atom_site.auth_seq_id"));
            let i_x = col_idx("_atom_site.Cartn_x");
            let i_y = col_idx("_atom_site.Cartn_y");
            let i_z = col_idx("_atom_site.Cartn_z");
            let i_occ = col_idx("_atom_site.occupancy");
            let i_bfac = col_idx("_atom_site.B_iso_or_equiv");
            let i_elem = col_idx("_atom_site.type_symbol");
            let i_model = col_idx("_atom_site.pdbx_PDB_model_num");

            // Parse data rows
            while i < lines.len() {
                let line = lines[i].trim();
                if line.is_empty() || line.starts_with('#') || line.starts_with("loop_") || line.starts_with('_') {
                    break;
                }
                let fields = split_cif_fields(line);
                if fields.len() < columns.len() {
                    i += 1;
                    continue;
                }

                let group = i_group.map(|idx| &fields[idx][..]).unwrap_or("ATOM");
                let is_hetatm = group == "HETATM";
                if !include_hetatm && is_hetatm {
                    i += 1;
                    continue;
                }

                let get = |idx: Option<usize>| idx.map(|j| fields[j].as_str()).unwrap_or("");

                atoms.push(Atom {
                    model_num: get(i_model).parse().unwrap_or(1),
                    chain: get(i_chain).to_string(),
                    resid: get(i_resid).parse().unwrap_or(0),
                    resname: get(i_resname).to_string(),
                    atom_name: get(i_name).to_string(),
                    altloc: ' ',
                    x: get(i_x).parse().unwrap_or(0.0),
                    y: get(i_y).parse().unwrap_or(0.0),
                    z: get(i_z).parse().unwrap_or(0.0),
                    occupancy: get(i_occ).parse().unwrap_or(1.0),
                    b_factor: get(i_bfac).parse().unwrap_or(0.0),
                    element: get(i_elem).to_string(),
                    is_hetatm,
                });
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    atoms
}

/// Split CIF data line into fields, handling single-quoted strings
fn split_cif_fields(line: &str) -> Vec<String> {
    let mut fields = Vec::new();
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    while i < chars.len() {
        while i < chars.len() && chars[i].is_whitespace() {
            i += 1;
        }
        if i >= chars.len() { break; }
        if chars[i] == '\'' {
            i += 1;
            let start = i;
            while i < chars.len() && !(chars[i] == '\'' && (i + 1 >= chars.len() || chars[i + 1].is_whitespace())) {
                i += 1;
            }
            fields.push(chars[start..i].iter().collect());
            if i < chars.len() { i += 1; } // skip closing quote
        } else {
            let start = i;
            while i < chars.len() && !chars[i].is_whitespace() {
                i += 1;
            }
            fields.push(chars[start..i].iter().collect());
        }
    }
    fields
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_cif() {
        let cif = "data_test\nloop_\n_atom_site.group_PDB\n_atom_site.label_atom_id\n_atom_site.label_comp_id\n_atom_site.label_asym_id\n_atom_site.label_seq_id\n_atom_site.Cartn_x\n_atom_site.Cartn_y\n_atom_site.Cartn_z\n_atom_site.type_symbol\nATOM CA ALA A 1 1.0 2.0 3.0 C\nATOM N   ALA A 1 4.0 5.0 6.0 N\n";
        let atoms = parse_cif(cif, false);
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].atom_name, "CA");
        assert!((atoms[0].x - 1.0).abs() < 0.01);
    }
}
