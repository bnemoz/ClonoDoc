//! Minimal GenBank parser (`docs/03_ARCHITECTURE.md` §4).
//!
//! Enough to populate a [`crate::model::Vector`]: the `LOCUS` line (name, length,
//! topology), the `FEATURES` table (type, location, qualifiers), and `ORIGIN`
//! (sequence). Feature locations are kept 1-based inclusive (GenBank convention)
//! and converted to 0-based only at the model boundary.

use crate::model::Feature;
use crate::seq;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct GbFeature {
    pub kind: String,
    /// 1-based inclusive start/end (min start / max end for join()).
    pub start: usize,
    pub end: usize,
    pub complement: bool,
    pub qualifiers: BTreeMap<String, String>,
}

#[derive(Debug, Clone)]
pub struct GbRecord {
    pub name: String,
    pub length: usize,
    pub circular: bool,
    pub sequence: String,
    pub features: Vec<GbFeature>,
}

pub fn read_path(path: &Path) -> Result<GbRecord> {
    let text = std::fs::read_to_string(path)?;
    parse(&text).map_err(|m| Error::parse(path.display().to_string(), m))
}

/// Parse a location string like `1..57`, `complement(58..85)`, `join(1..3,10..20)`.
/// Returns `(start, end, complement)` 1-based inclusive (min start / max end).
fn parse_location(loc: &str) -> Option<(usize, usize, bool)> {
    let loc = loc.trim();
    let complement = loc.contains("complement");
    // Strip wrappers, keep only digits and `..` and commas.
    let digits: String = loc
        .chars()
        .map(|c| {
            if c.is_ascii_digit() || c == '.' || c == ',' {
                c
            } else {
                ' '
            }
        })
        .collect();
    let mut starts = Vec::new();
    let mut ends = Vec::new();
    for part in digits.split([',', ' ']).filter(|p| !p.is_empty()) {
        let nums: Vec<usize> = part
            .split("..")
            .filter(|s| !s.is_empty())
            .filter_map(|s| s.trim_matches('.').parse::<usize>().ok())
            .collect();
        if let (Some(&a), Some(&b)) = (nums.first(), nums.last()) {
            starts.push(a.min(b));
            ends.push(a.max(b));
        }
    }
    match (starts.iter().min(), ends.iter().max()) {
        (Some(&s), Some(&e)) => Some((s, e, complement)),
        _ => None,
    }
}

pub fn parse(text: &str) -> std::result::Result<GbRecord, String> {
    let mut name = String::new();
    let mut length = 0usize;
    let mut circular = false;
    let mut sequence = String::new();
    let mut features: Vec<GbFeature> = Vec::new();

    let mut section = Section::Header;
    let mut cur: Option<GbFeature> = None;
    let mut cur_qual_key: Option<String> = None;

    for raw in text.lines() {
        if raw.starts_with("LOCUS") {
            let toks: Vec<&str> = raw.split_whitespace().collect();
            if toks.len() >= 2 {
                name = toks[1].to_string();
            }
            // length is the token followed by "bp"
            for w in toks.windows(2) {
                if w[1] == "bp" {
                    length = w[0].parse().unwrap_or(0);
                }
            }
            circular = raw.to_ascii_lowercase().contains("circular");
            continue;
        }
        if raw.starts_with("FEATURES") {
            section = Section::Features;
            continue;
        }
        if raw.starts_with("ORIGIN") {
            if let Some(f) = cur.take() {
                features.push(f);
            }
            section = Section::Origin;
            continue;
        }
        if raw.starts_with("//") {
            break;
        }

        match section {
            Section::Header => {}
            Section::Features => {
                // A feature key line starts at column 5 (after 5 spaces) with a
                // non-space in column 6; qualifier lines start with `/` after 21 spaces.
                let trimmed = raw.trim_start();
                let indent = raw.len() - trimmed.len();
                if indent <= 7 && !trimmed.is_empty() && !trimmed.starts_with('/') {
                    // New feature line: "<key>   <location>"
                    if let Some(f) = cur.take() {
                        features.push(f);
                    }
                    cur_qual_key = None;
                    let mut parts = trimmed.splitn(2, char::is_whitespace);
                    let kind = parts.next().unwrap_or("").to_string();
                    let loc = parts.next().unwrap_or("").trim();
                    if let Some((s, e, comp)) = parse_location(loc) {
                        cur = Some(GbFeature {
                            kind,
                            start: s,
                            end: e,
                            complement: comp,
                            qualifiers: BTreeMap::new(),
                        });
                    } else {
                        cur = Some(GbFeature {
                            kind,
                            start: 0,
                            end: 0,
                            complement: false,
                            qualifiers: BTreeMap::new(),
                        });
                    }
                } else if let Some(body) = trimmed.strip_prefix('/') {
                    // Qualifier: /key="value" or /key=value
                    if let Some(f) = cur.as_mut() {
                        if let Some(eq) = body.find('=') {
                            let key = body[..eq].to_string();
                            let val = body[eq + 1..].trim_matches('"').to_string();
                            cur_qual_key = Some(key.clone());
                            f.qualifiers.insert(key, val);
                        } else {
                            cur_qual_key = Some(body.to_string());
                            f.qualifiers.insert(body.to_string(), String::new());
                        }
                    }
                } else if let (Some(f), Some(k)) = (cur.as_mut(), cur_qual_key.as_ref()) {
                    // Continuation of a multi-line qualifier value.
                    let cont = trimmed.trim_matches('"');
                    if let Some(v) = f.qualifiers.get_mut(k) {
                        if !v.is_empty() {
                            v.push(' ');
                        }
                        v.push_str(cont);
                    }
                }
            }
            Section::Origin => {
                sequence.push_str(&seq::clean(raw));
            }
        }
    }

    if sequence.is_empty() {
        return Err("no ORIGIN sequence found".into());
    }
    if length == 0 {
        length = sequence.len();
    }
    Ok(GbRecord {
        name,
        length,
        circular,
        sequence,
        features,
    })
}

enum Section {
    Header,
    Features,
    Origin,
}

/// Map the parsed feature table to the canonical role spans the model expects,
/// by keyword (`docs/03_ARCHITECTURE.md` §4). 1-based inclusive spans preserved.
pub fn map_roles(rec: &GbRecord) -> BTreeMap<String, Feature> {
    let mut roles: BTreeMap<String, Feature> = BTreeMap::new();
    let name_of = |f: &GbFeature| -> Option<String> {
        f.qualifiers
            .get("standard_name")
            .or_else(|| f.qualifiers.get("gene"))
            .or_else(|| f.qualifiers.get("note"))
            .cloned()
    };
    for f in &rec.features {
        if f.start == 0 {
            continue;
        }
        let hay = format!(
            "{} {} {}",
            f.kind,
            name_of(f).unwrap_or_default(),
            f.qualifiers.get("note").cloned().unwrap_or_default()
        )
        .to_ascii_lowercase();

        let role = if f.kind.eq_ignore_ascii_case("sig_peptide")
            || hay.contains("signal peptide")
            || hay.contains("sig_peptide")
        {
            Some("signal_peptide")
        } else if hay.contains("c region")
            || hay.contains("ighg")
            || hay.contains("igha")
            || hay.contains("igkc")
            || hay.contains("iglc")
            || hay.contains("constant")
        {
            Some("constant_region")
        } else if (hay.contains("ampr") || hay.contains("kanr") || hay.contains("resistance"))
            && !hay.contains("promoter")
        {
            // The resistance *gene* (CDS), not its promoter.
            Some("resistance")
        } else if hay.contains("ori") || hay.contains("origin") || hay.contains("puc") {
            Some("origin")
        } else if hay.contains("cmv") || (hay.contains("promoter") && !hay.contains("ampr")) {
            Some("promoter")
        } else {
            None
        };

        if let Some(role) = role {
            // Keep the first (usually broadest) match per role; don't overwrite
            // a constant_region with a sub-feature like enhancer.
            roles.entry(role.to_string()).or_insert(Feature {
                start: f.start,
                end: f.end,
                name: name_of(f),
            });
        }
    }
    roles
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_location_forms() {
        assert_eq!(parse_location("1..57"), Some((1, 57, false)));
        assert_eq!(parse_location("complement(58..85)"), Some((58, 85, true)));
        assert_eq!(parse_location("join(1..3,10..20)"), Some((1, 20, false)));
    }

    #[test]
    fn parses_tiny_record() {
        let gb = "LOCUS       test        12 bp    DNA     circular UNA 22-MAY-2026\n\
                  FEATURES             Location/Qualifiers\n     \
                  sig_peptide     1..6\n                     /standard_name=\"Signal peptide\"\n\
                  ORIGIN\n        1 atggaa ctgggg\n//\n";
        let rec = parse(gb).unwrap();
        assert_eq!(rec.length, 12);
        assert!(rec.circular);
        assert_eq!(rec.sequence, "ATGGAACTGGGG");
        assert_eq!(rec.features.len(), 1);
        assert_eq!(rec.features[0].start, 1);
        let roles = map_roles(&rec);
        assert!(roles.contains_key("signal_peptide"));
    }
}
