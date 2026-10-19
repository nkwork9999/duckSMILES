//! Drawing graph prepared by the same pipeline as descriptor calculations.
//! Separate from the legacy raw graph API. No descriptor arithmetic changes.
use crate::parser::BondOrder;
use std::fmt::Write;

const SCHEMA: &str = "ducksmiles-workbench-drawing-v1";

fn json_string(s: &str) -> String {
    let mut out = String::from("\"");
    for c in s.chars() {
        match c {
            '\\' => out.push_str("\\\\"),
            '"' => out.push_str("\\\""),
            c if c.is_control() => {
                let _ = write!(out, "\\u{:04x}", c as u32);
            }
            c => out.push(c),
        }
    }
    out.push('"');
    out
}

fn number_or_null<T: std::fmt::Display>(value: Option<T>) -> String {
    value
        .map(|n| n.to_string())
        .unwrap_or_else(|| "null".to_string())
}

/// Prepared input-order graph. Unsupported stereo is disclosed, not invented.
pub(crate) fn calculate(source: &str) -> Option<String> {
    let (mol, _) = crate::workbench::prepare(source)?;
    if mol.atoms.len() > crate::native_graph::MAX_ATOMS
        || mol.bonds.len() > crate::native_graph::MAX_BONDS
    {
        return None;
    }
    let atomic_numbers = crate::electronic::atomic_numbers(&mol)?;
    let unrepresented_stereo = source.contains('|')
        || ["@SP", "@TB", "@OH", "@AL"]
            .iter()
            .any(|tag| source.contains(tag));
    let mut out = format!(
        "{{\"schema\":\"{SCHEMA}\",\"input_smiles\":{},\"preparation_revision\":\"{}\",\"atom_order\":\"input_parse_order\",\"unrepresented_stereo\":{},\"atoms\":[",
        json_string(source), crate::workbench::REVISION, unrepresented_stereo
    );
    for (i, atom) in mol.atoms.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let explicit_neighbors = mol
            .neighbors(i)
            .iter()
            .filter(|(j, _)| mol.atoms[*j].symbol == "H")
            .count();
        let bracket_h = if atom.in_bracket && !atom.prepared_h_is_implicit {
            atom.hydrogen
        } else {
            0
        };
        let implicit_h = atom.hydrogen - bracket_h;
        let chirality = atom
            .chirality
            .as_ref()
            .map(|s| json_string(s))
            .unwrap_or_else(|| "null".into());
        let neighbors = if atom.chirality.is_some() {
            atom.nbr_order
                .iter()
                .map(|n| n.to_string())
                .collect::<Vec<_>>()
                .join(",")
        } else {
            String::new()
        };
        let _ = write!(out,
            "{{\"index\":{i},\"symbol\":{},\"atomic_number\":{},\"isotope\":{},\"formal_charge\":{},\"aromatic\":{},\"in_bracket\":{},\"bracket_h_count\":{bracket_h},\"implicit_h_count\":{implicit_h},\"explicit_h_neighbor_count\":{explicit_neighbors},\"atom_map\":{},\"chirality_token\":{chirality},\"stereo_neighbor_order\":[{neighbors}]}}",
            json_string(&atom.symbol), atomic_numbers[i], number_or_null(atom.isotope),
            atom.charge, atom.aromatic, atom.in_bracket, number_or_null(atom.atom_map));
    }
    out.push_str("],\"bonds\":[");
    for (i, bond) in mol.bonds.iter().enumerate() {
        if i != 0 {
            out.push(',');
        }
        let kind = match bond.order {
            BondOrder::Single => "single",
            BondOrder::Double => "double",
            BondOrder::Triple => "triple",
            BondOrder::Aromatic => "aromatic",
            BondOrder::Dative => "dative",
        };
        let _ = write!(
            out,
            "{{\"index\":{i},\"begin\":{},\"end\":{},\"bond_type\":\"{kind}\",\"direction\":{}}}",
            bond.a, bond.b, bond.direction
        );
    }
    out.push_str("]}");
    Some(out)
}

/// Prepared drawing JSON. Returns required bytes (without NUL), or -1.
/// Short buffers are not written. Null output is valid only with zero capacity.
/// # Safety
/// Input/output must address the specified disjoint readable/writable ranges.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn ds_workbench_drawing_graph(
    input: *const u8,
    length: usize,
    output: *mut u8,
    capacity: usize,
) -> i32 {
    std::panic::catch_unwind(|| {
        if input.is_null()
            || length == 0
            || length > crate::native_graph::MAX_INPUT_BYTES
            || (output.is_null() && capacity != 0)
        {
            return -1;
        }
        // SAFETY: the caller guarantees readable input; the byte limit is checked.
        let bytes = unsafe { std::slice::from_raw_parts(input, length) };
        let Ok(source) = std::str::from_utf8(bytes) else {
            return -1;
        };
        match calculate(source) {
            Some(json) => crate::write_required(json.as_bytes(), output, capacity),
            None => -1,
        }
    })
    .unwrap_or(-1)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn drawing_uses_prepared_nitro_charges_and_bond_orders() {
        let json = calculate("CN(=O)=O").expect("prepared nitro graph");
        assert!(json.contains("\"formal_charge\":1"));
        assert!(json.contains("\"formal_charge\":-1"));
    }

    #[test]
    fn drawing_keeps_directed_coordination_and_unrepresented_stereo_notice() {
        let json = calculate("N->[Cu]").expect("dative drawing graph");
        assert!(json.contains("\"begin\":0,\"end\":1,\"bond_type\":\"dative\""));
        let stereo =
            calculate("[Pt@SP1](Cl)(F)(Br)I").expect("non-tetrahedral input remains available");
        assert!(stereo.contains("\"unrepresented_stereo\":true"));
    }

    #[test]
    fn drawing_abi_queries_capacity_preserves_short_buffers_and_rejects_bad_inputs() {
        let source = b"CCO";
        let n = unsafe {
            ds_workbench_drawing_graph(source.as_ptr(), source.len(), std::ptr::null_mut(), 0)
        };
        assert!(n > 0, "drawing ABI must return a required byte count");
        let mut short = [0x5a_u8; 8];
        assert_eq!(
            unsafe {
                ds_workbench_drawing_graph(
                    source.as_ptr(),
                    source.len(),
                    short.as_mut_ptr(),
                    short.len(),
                )
            },
            n
        );
        assert_eq!(short, [0x5a; 8], "size query must not partially write JSON");
        let mut out = vec![0x5a_u8; n as usize + 4];
        assert_eq!(
            unsafe {
                ds_workbench_drawing_graph(
                    source.as_ptr(),
                    source.len(),
                    out.as_mut_ptr(),
                    n as usize,
                )
            },
            n
        );
        assert_eq!(&out[n as usize..], &[0x5a; 4]);
        assert!(std::str::from_utf8(&out[..n as usize])
            .unwrap()
            .contains("\"input_smiles\":\"CCO\""));
        assert_eq!(
            unsafe { ds_workbench_drawing_graph(std::ptr::null(), 3, std::ptr::null_mut(), 0) },
            -1
        );
        assert_eq!(
            unsafe {
                ds_workbench_drawing_graph(source.as_ptr(), source.len(), std::ptr::null_mut(), 1)
            },
            -1
        );
        let bad = [0xff];
        assert_eq!(
            unsafe { ds_workbench_drawing_graph(bad.as_ptr(), bad.len(), std::ptr::null_mut(), 0) },
            -1
        );
        assert!(
            calculate(&"C".repeat(257)).is_none(),
            "drawing budget remains256 atoms"
        );
    }
}
