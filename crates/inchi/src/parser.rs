/// Parse InChI string layers. InChI format:
/// `InChI=1S/<formula>/c<connections>/h<hydrogens>/q<charge>/p<protons>/b<stereo_bond>/t<stereo_tet>/m<type>/s<inv>/i<isotope>`

const INCHI_PREFIX: &[u8] = b"InChI=";
const INCHI_STD_PREFIX: &[u8] = b"InChI=1S/";

pub fn is_valid(s: &[u8]) -> bool {
    s.starts_with(INCHI_PREFIX) && s.len() > INCHI_PREFIX.len() + 2
}

pub fn is_standard(s: &[u8]) -> bool {
    s.starts_with(INCHI_STD_PREFIX)
}

/// Extract version: "1S" or "1"
pub fn version(s: &[u8]) -> Option<&[u8]> {
    if !s.starts_with(INCHI_PREFIX) {
        return None;
    }
    let rest = &s[INCHI_PREFIX.len()..];
    let end = rest.iter().position(|&b| b == b'/')?;
    Some(&rest[..end])
}

/// Extract molecular formula (first layer after version)
pub fn formula(s: &[u8]) -> Option<&[u8]> {
    let prefix_len = if s.starts_with(INCHI_STD_PREFIX) {
        INCHI_STD_PREFIX.len()
    } else if s.starts_with(INCHI_PREFIX) {
        // Find first '/' after "InChI=X"
        let rest = &s[INCHI_PREFIX.len()..];
        let slash = rest.iter().position(|&b| b == b'/')?;
        INCHI_PREFIX.len() + slash + 1
    } else {
        return None;
    };
    let rest = &s[prefix_len..];
    let end = rest.iter().position(|&b| b == b'/').unwrap_or(rest.len());
    if end == 0 {
        return None;
    }
    Some(&rest[..end])
}

/// Extract a specific layer by its prefix character (c, h, q, p, b, t, m, s, i, f, r)
pub fn layer(s: &[u8], prefix: u8) -> Option<&[u8]> {
    // Find "/<prefix>" in the string
    let needle = [b'/', prefix];
    let mut pos = 0;
    while pos + 1 < s.len() {
        if s[pos] == needle[0] && s[pos + 1] == needle[1] {
            let start = pos + 2;
            let end = s[start..]
                .iter()
                .position(|&b| b == b'/')
                .map(|p| start + p)
                .unwrap_or(s.len());
            return Some(&s[start..end]);
        }
        pos += 1;
    }
    None
}

pub fn has_stereo(s: &[u8]) -> bool {
    layer(s, b'b').is_some() || layer(s, b't').is_some()
}

pub fn num_stereo_centers(s: &[u8]) -> i32 {
    match layer(s, b't') {
        Some(t_layer) => {
            if t_layer.is_empty() {
                0
            } else {
                t_layer.iter().filter(|&&b| b == b',').count() as i32 + 1
            }
        }
        None => 0,
    }
}

/// Compare skeleton: formula + connections + hydrogens must match
pub fn skeleton_match(a: &[u8], b: &[u8]) -> bool {
    let fa = formula(a);
    let fb = formula(b);
    if fa != fb {
        return false;
    }
    let ca = layer(a, b'c');
    let cb = layer(b, b'c');
    if ca != cb {
        return false;
    }
    let ha = layer(a, b'h');
    let hb = layer(b, b'h');
    ha == hb
}

// --- InChIKey ---

pub fn inchikey_is_valid(s: &[u8]) -> bool {
    // Format: XXXXXXXXXXXXXX-XXXXXXXXXX-X (14-10-1 with dashes)
    if s.len() != 27 {
        return false;
    }
    if s[14] != b'-' || s[25] != b'-' {
        return false;
    }
    s.iter().enumerate().all(|(i, &b)| {
        if i == 14 || i == 25 {
            b == b'-'
        } else {
            b.is_ascii_uppercase()
        }
    })
}

pub fn inchikey_connectivity(s: &[u8]) -> Option<&[u8]> {
    if s.len() >= 14 {
        Some(&s[..14])
    } else {
        None
    }
}

pub fn inchikey_stereo(s: &[u8]) -> Option<&[u8]> {
    if s.len() >= 25 {
        Some(&s[15..25])
    } else {
        None
    }
}

pub fn inchikey_protonation(s: &[u8]) -> Option<&[u8]> {
    if s.len() == 27 {
        Some(&s[26..27])
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const ACETIC: &[u8] = b"InChI=1S/C2H4O2/c1-2(3)4/h1H3,(H,3,4)";
    const KEY: &[u8] = b"QTBSBXVTEAMEQO-UHFFFAOYSA-N";

    #[test]
    fn test_valid() {
        assert!(is_valid(ACETIC));
        assert!(!is_valid(b"garbage"));
    }

    #[test]
    fn test_standard() {
        assert!(is_standard(ACETIC));
        assert!(!is_standard(b"InChI=1/C2H4O2/c"));
    }

    #[test]
    fn test_version() {
        assert_eq!(version(ACETIC), Some(b"1S" as &[u8]));
    }

    #[test]
    fn test_formula() {
        assert_eq!(formula(ACETIC), Some(b"C2H4O2" as &[u8]));
    }

    #[test]
    fn test_layers() {
        assert_eq!(layer(ACETIC, b'c'), Some(b"1-2(3)4" as &[u8]));
        assert_eq!(layer(ACETIC, b'h'), Some(b"1H3,(H,3,4)" as &[u8]));
        assert_eq!(layer(ACETIC, b'q'), None);
    }

    #[test]
    fn test_stereo() {
        assert!(!has_stereo(ACETIC));
        let stereo = b"InChI=1S/C4H8/c1-3-4-2/h3-4H,1-2H3/b3-4+";
        assert!(has_stereo(stereo));
    }

    #[test]
    fn test_stereo_centers() {
        let tet = b"InChI=1S/C21H30O4/c1-2/t13-,14+,15+";
        assert_eq!(num_stereo_centers(tet), 3);
    }

    #[test]
    fn test_skeleton_match() {
        assert!(skeleton_match(ACETIC, ACETIC));
    }

    #[test]
    fn test_inchikey() {
        assert!(inchikey_is_valid(KEY));
        assert!(!inchikey_is_valid(b"short"));
        assert_eq!(inchikey_connectivity(KEY), Some(b"QTBSBXVTEAMEQO" as &[u8]));
        assert_eq!(inchikey_stereo(KEY), Some(b"UHFFFAOYSA" as &[u8]));
        assert_eq!(inchikey_protonation(KEY), Some(b"N" as &[u8]));
    }
}
