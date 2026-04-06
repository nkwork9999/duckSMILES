/// SELFIES → SMILES decoder
///
/// Parses SELFIES tokens and reconstructs SMILES string.
/// Every valid SELFIES string produces a valid molecule (by design).

/// Decode SELFIES to SMILES. Returns None if input is empty.
pub fn selfies_to_smiles(selfies: &str) -> Option<String> {
    if selfies.is_empty() {
        return None;
    }

    let tokens = tokenize(selfies);
    if tokens.is_empty() {
        return None;
    }

    let mut result = String::new();

    let mut i = 0;
    while i < tokens.len() {
        let token = &tokens[i];

        if token == "." {
            result.push('.');
            i += 1;
            continue;
        }

        if token.starts_with("Branch") {
            // BranchN: next N tokens describe branch length index(es), then content
            let level: usize = token.chars().last()
                .and_then(|c| c.to_digit(10))
                .unwrap_or(1) as usize;

            // Skip N index tokens
            let skip = level;
            i += 1; // past BranchN

            // Collect branch body tokens (after index tokens)
            let index_end = (i + skip).min(tokens.len());

            // The branch content is everything after the index tokens until
            // the next token at the same level. For simplicity, take one token
            // as branch content per Branch1, or up to the count for Branch2+.
            i = index_end;
            if i < tokens.len() {
                let atom = decode_atom_token(&tokens[i]);
                result.push('(');
                result.push_str(&atom);
                result.push(')');
                i += 1;
            }
            continue;
        }

        if token.starts_with("Ring") {
            let level: usize = token.chars().last()
                .and_then(|c| c.to_digit(10))
                .unwrap_or(1) as usize;
            // Skip index tokens
            i += 1 + level;
            // Ring closures would need tracking; simplified: skip
            continue;
        }

        // Regular atom token
        let atom = decode_atom_token(token);
        result.push_str(&atom);
        i += 1;
    }

    if result.is_empty() {
        None
    } else {
        Some(result)
    }
}

/// Tokenize SELFIES string: split by brackets
fn tokenize(selfies: &str) -> Vec<String> {
    let mut tokens = Vec::new();
    let chars: Vec<char> = selfies.chars().collect();
    let mut i = 0;

    while i < chars.len() {
        if chars[i] == '[' {
            let start = i + 1;
            i += 1;
            while i < chars.len() && chars[i] != ']' {
                i += 1;
            }
            if i < chars.len() {
                let token: String = chars[start..i].iter().collect();
                tokens.push(token);
                i += 1;
            }
        } else {
            i += 1;
        }
    }
    tokens
}

/// Convert a single SELFIES atom token to SMILES
fn decode_atom_token(token: &str) -> String {
    if token == "nop" {
        return String::new();
    }
    if token == "." {
        return ".".to_string();
    }

    let mut result = String::new();

    // Bond prefix
    let rest = if token.starts_with('=') {
        result.push('=');
        &token[1..]
    } else if token.starts_with('#') {
        result.push('#');
        &token[1..]
    } else {
        token
    };

    // Atom symbol
    let needs_bracket = rest.contains('+') || rest.contains('-')
        || rest.contains('H') && rest.len() > 1 && !matches!(rest, "H")
        || rest.contains('@')
        || !is_organic_subset(rest);

    if needs_bracket && !rest.starts_with('[') {
        // Check if it's a simple aromatic or organic atom
        if is_simple_atom(rest) {
            result.push_str(rest);
        } else {
            result.push('[');
            result.push_str(rest);
            result.push(']');
        }
    } else {
        result.push_str(rest);
    }

    result
}

fn is_organic_subset(s: &str) -> bool {
    matches!(s, "B" | "C" | "N" | "O" | "P" | "S" | "F" | "Cl" | "Br" | "I"
             | "c" | "n" | "o" | "s" | "p" | "b")
}

fn is_simple_atom(s: &str) -> bool {
    is_organic_subset(s) || matches!(s, "H")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tokenize() {
        let tokens = tokenize("[C][=O][N]");
        assert_eq!(tokens, vec!["C", "=O", "N"]);
    }

    #[test]
    fn test_decode_simple() {
        let s = selfies_to_smiles("[C][C][O]").unwrap();
        assert_eq!(s, "CCO");
    }

    #[test]
    fn test_decode_double_bond() {
        let s = selfies_to_smiles("[C][=O]").unwrap();
        assert_eq!(s, "C=O");
    }

    #[test]
    fn test_decode_triple() {
        let s = selfies_to_smiles("[C][#N]").unwrap();
        assert_eq!(s, "C#N");
    }

    #[test]
    fn test_decode_bracket_atom() {
        let s = selfies_to_smiles("[NH4+]").unwrap();
        assert_eq!(s, "[NH4+]");
    }

    #[test]
    fn test_decode_branch() {
        let s = selfies_to_smiles("[C][Branch1][C][=O][O]").unwrap();
        // Should produce C(=O)O
        assert!(s.contains("(=O)"));
    }

    #[test]
    fn test_empty() {
        assert!(selfies_to_smiles("").is_none());
    }

    #[test]
    fn test_nop() {
        let s = selfies_to_smiles("[C][nop][C]").unwrap();
        assert_eq!(s, "CC");
    }
}
