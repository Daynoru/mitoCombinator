use anyhow::Result;
use std::collections::HashMap;

/// Restriction enzyme with recognition site
#[derive(Debug, Clone)]
pub struct Enzyme {
    pub name: String,
    pub recognition: String,
    pub cut_top: i32,
    pub cut_bottom: i32,
}

/// Enzyme name -> list of site positions in genome
pub type RestrictionMap = HashMap<String, Vec<usize>>;

/// Parse CSV file: Enzyme,Recognition
pub fn parse_enzymes(path: &str) -> Result<Vec<Enzyme>> {
    let content = std::fs::read_to_string(path)?;
    let mut enzymes = Vec::new();

    for line in content.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with("Enzyme") {
            continue;
        }
        let parts: Vec<&str> = line.splitn(2, ',').collect();
        if parts.len() < 2 {
            continue;
        }
        let name = parts[0].trim().to_string();
        let raw = parts[1].trim().to_string();
        if name.is_empty() || raw.is_empty() {
            continue;
        }
        let (cut_top, cut_bottom) = parse_cut_positions(&raw);
        enzymes.push(Enzyme { name, recognition: raw, cut_top, cut_bottom });
    }
    eprintln!("Loaded {} enzymes from {}", enzymes.len(), path);
    Ok(enzymes)
}

/// Extract cut positions from notation like "G/AATTC" or "CCGC(-3/-1)"
fn parse_cut_positions(raw: &str) -> (i32, i32) {
    if raw.contains('/') && !raw.contains('(') {
        let pos = raw.find('/').unwrap_or(0) as i32;
        return (pos, pos);
    }
    if let Some(paren) = raw.find('(') {
        let inner = &raw[paren + 1..].trim_end_matches(')');
        let parts: Vec<&str> = inner.split('/').collect();
        if parts.len() == 2 {
            let top: i32 = parts[0].replace("none", "0").parse().unwrap_or(0);
            let bottom: i32 = parts[1].replace("none", "0").parse().unwrap_or(0);
            return (top, bottom);
        }
    }
    (0, 0)
}

/// Convert recognition site (with IUPAC ambiguity codes) to plain DNA alphabet
pub fn recognition_to_pattern(recognition: &str) -> String {
    recognition.chars().filter(|c| c.is_alphabetic()).collect()
}

/// Expand IUPAC ambiguity code to possible nucleotides
fn expand_iupac(c: char) -> &'static [u8] {
    match c.to_ascii_uppercase() {
        'A' => b"A",   'C' => b"C",   'G' => b"G",   'T' | 'U' => b"T",
        'R' => b"AG",  'Y' => b"CT",  'M' => b"AC",  'K' => b"GT",
        'S' => b"CG",  'W' => b"AT",  'H' => b"ACT", 'B' => b"CGT",
        'V' => b"ACG", 'D' => b"AGT", 'N' => b"ACGT",
        _ => b"N",
    }
}

/// Find all positions of a restriction site in sequence
pub fn find_sites(seq: &[u8], enzyme: &Enzyme) -> Vec<usize> {
    let pattern = recognition_to_pattern(&enzyme.recognition);
    if pattern.is_empty() {
        return Vec::new();
    }
    let mut positions = Vec::new();
    for i in 0..seq.len().saturating_sub(pattern.len()) {
        if matches_iupac(&seq[i..i + pattern.len()], &pattern) {
            positions.push(i);
        }
    }
    positions
}

/// Check if a DNA segment matches an IUPAC pattern
fn matches_iupac(seq: &[u8], pattern: &str) -> bool {
    seq.iter()
        .zip(pattern.chars())
        .all(|(&b, p)| expand_iupac(p).contains(&b.to_ascii_uppercase()))
}

/// Build full restriction map: enzyme -> all site positions in genome
pub fn build_map(seq: &[u8], enzymes: &[Enzyme]) -> RestrictionMap {
    let mut map = RestrictionMap::new();
    for enzyme in enzymes {
        let sites = find_sites(seq, enzyme);
        if !sites.is_empty() {
            map.insert(enzyme.name.clone(), sites);
        }
    }
    map
}