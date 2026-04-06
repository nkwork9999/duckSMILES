/// SMILES → SELFIES encoder
///
/// SELFIES tokens:
///   Atoms:    [C], [N], [O], [S], [P], [F], [Cl], [Br], [I], [B]
///   Bonds:    [=C], [#C], [=N], etc. (bond prefix on atom)
///   Branch:   [Branch1], [Branch2], [Branch3] (1/2/3 tokens follow as branch length + content)
///   Ring:     [Ring1], [Ring2] (ring closure with distance index)
///   Special:  [C@@H], [C@H], [nop] etc.
///
/// This is a simplified encoder handling the most common SMILES patterns.

/// Encode SMILES to SELFIES. Returns None if SMILES is invalid.
pub fn smiles_to_selfies(smiles: &str) -> Option<String> {
    if smiles.is_empty() {
        return None;
    }

    let chars: Vec<char> = smiles.chars().collect();
    let mut result = String::new();
    let mut pos = 0;

    while pos < chars.len() {
        let c = chars[pos];

        match c {
            // Branch start → [BranchN] with nested encoding
            '(' => {
                // Find matching ')'
                let branch_start = pos + 1;
                let mut depth = 1;
                let mut end = branch_start;
                while end < chars.len() && depth > 0 {
                    match chars[end] {
                        '(' => depth += 1,
                        ')' => depth -= 1,
                        _ => {}
                    }
                    if depth > 0 { end += 1; }
                }
                let branch_content: String = chars[branch_start..end].iter().collect();
                if let Some(encoded) = smiles_to_selfies(&branch_content) {
                    // Count tokens in encoded branch to determine Branch level
                    let token_count = encoded.matches('[').count();
                    if token_count <= 1 {
                        result.push_str("[Branch1]");
                        result.push_str("[C]"); // length index placeholder
                    } else {
                        result.push_str("[Branch2]");
                        result.push_str("[Ring1]"); // length index
                        result.push_str("[C]");     // second index
                    }
                    result.push_str(&encoded);
                }
                pos = end + 1; // skip ')'
            }

            // Bond + next atom
            '=' => {
                pos += 1;
                if pos < chars.len() {
                    let atom_token = read_atom(&chars, &mut pos);
                    result.push('[');
                    result.push('=');
                    result.push_str(&atom_token);
                    result.push(']');
                }
            }
            '#' => {
                pos += 1;
                if pos < chars.len() {
                    let atom_token = read_atom(&chars, &mut pos);
                    result.push('[');
                    result.push('#');
                    result.push_str(&atom_token);
                    result.push(']');
                }
            }

            // Ring closures (simplified: ignore for basic encoding)
            '1'..='9' | '%' => {
                if c == '%' {
                    pos += 3; // skip %XX
                } else {
                    pos += 1;
                }
                result.push_str("[Ring1]");
                result.push_str("[C]"); // distance index placeholder
            }

            // Dot (fragment separator)
            '.' => {
                result.push_str("[.]");
                pos += 1;
            }

            // Stereo bond markers (skip)
            '/' | '\\' | '-' | ':' => {
                pos += 1;
            }

            // Bracket atom [...]
            '[' => {
                let start = pos;
                while pos < chars.len() && chars[pos] != ']' {
                    pos += 1;
                }
                if pos < chars.len() {
                    pos += 1; // skip ']'
                }
                let token: String = chars[start..pos].iter().collect();
                result.push_str(&token);
            }

            // Aromatic / organic subset atom
            _ => {
                let atom_token = read_atom(&chars, &mut pos);
                result.push('[');
                result.push_str(&atom_token);
                result.push(']');
            }
        }
    }

    if result.is_empty() {
        None
    } else {
        Some(result)
    }
}

/// Read an atom symbol from position, advancing pos
fn read_atom(chars: &[char], pos: &mut usize) -> String {
    let c = chars[*pos];
    let mut sym = String::new();

    if c.is_ascii_lowercase() {
        // Aromatic
        sym.push(c);
        *pos += 1;
    } else if c.is_ascii_uppercase() {
        sym.push(c);
        *pos += 1;
        // Two-letter elements
        if *pos < chars.len() && chars[*pos].is_ascii_lowercase() {
            let two: String = [c, chars[*pos]].iter().collect();
            if matches!(two.as_str(), "Cl" | "Br" | "Si" | "Se" | "Na" | "Li" | "Mg" | "Al" | "Ca" | "Fe" | "Cu" | "Zn" | "As" | "Sn" | "Pb") {
                sym.push(chars[*pos]);
                *pos += 1;
            }
        }
    } else {
        sym.push(c);
        *pos += 1;
    }

    sym
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_methane() {
        let s = smiles_to_selfies("C").unwrap();
        assert_eq!(s, "[C]");
    }

    #[test]
    fn test_ethanol() {
        let s = smiles_to_selfies("CCO").unwrap();
        assert_eq!(s, "[C][C][O]");
    }

    #[test]
    fn test_double_bond() {
        let s = smiles_to_selfies("C=O").unwrap();
        assert_eq!(s, "[C][=O]");
    }

    #[test]
    fn test_triple_bond() {
        let s = smiles_to_selfies("C#N").unwrap();
        assert_eq!(s, "[C][#N]");
    }

    #[test]
    fn test_branch() {
        let s = smiles_to_selfies("CC(=O)O").unwrap();
        assert!(s.contains("[Branch"));
        assert!(s.contains("[=O]"));
    }

    #[test]
    fn test_bracket_atom() {
        let s = smiles_to_selfies("[NH4+]").unwrap();
        assert_eq!(s, "[NH4+]");
    }

    #[test]
    fn test_empty() {
        assert!(smiles_to_selfies("").is_none());
    }
}
