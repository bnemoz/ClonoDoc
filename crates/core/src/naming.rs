//! Sample-name parsing & chain pairing (`docs/04_NAMING.md`).
//!
//! Extracts `(ab_id, chain_class)` from an arbitrary record id using a
//! configurable [`NamingProfile`]. Forgiving by design; when it cannot decide
//! confidently it returns `needs_confirmation = true` rather than guessing.

use crate::model::ChainClass;
pub use crate::model::NamingProfile;
use std::collections::BTreeMap;

const BOUNDARY: char = '\u{1f}'; // unit separator used as an internal token boundary

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct NameParse {
    pub ab_id: String,
    pub chain_class: ChainClass,
    pub confidence: Confidence,
    pub needs_confirmation: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Confidence {
    High,
    Low,
    None,
}

/// Map a synonym string to its chain class via the profile's synonym table key.
fn class_of(key: &str) -> ChainClass {
    match key {
        "heavy" => ChainClass::Heavy,
        "kappa" => ChainClass::Kappa,
        "lambda" => ChainClass::Lambda,
        "light" => ChainClass::Light,
        _ => ChainClass::Unknown,
    }
}

/// Build a flat list of `(synonym_lower, class_key, is_single_letter)`, sorted
/// longest-first so `heavychain` matches before `heavy` before `hc` before `h`.
fn synonym_index(profile: &NamingProfile) -> Vec<(String, String, bool)> {
    let mut out: Vec<(String, String, bool)> = Vec::new();
    for (class_key, syns) in &profile.chain_synonyms {
        for s in syns {
            let lower = s.to_ascii_lowercase();
            let single = lower.chars().count() == 1;
            out.push((lower, class_key.clone(), single));
        }
    }
    // Longest-match-first; stable on ties by class key for determinism.
    out.sort_by(|a, b| b.0.len().cmp(&a.0.len()).then(a.1.cmp(&b.1)));
    out
}

/// Parse a record id with the given profile.
pub fn parse_name(record_id: &str, profile: &NamingProfile) -> NameParse {
    // 1. Normalize separators into a single boundary char; keep original tokens too.
    let mut normalized = record_id.to_string();
    for sep in &profile.separators {
        if !sep.is_empty() {
            normalized = normalized.replace(sep.as_str(), &BOUNDARY.to_string());
        }
    }
    // 2. Tokenize on the boundary (drop empties).
    let tokens: Vec<&str> = normalized
        .split(BOUNDARY)
        .filter(|t| !t.is_empty())
        .collect();
    if tokens.is_empty() {
        return NameParse {
            ab_id: record_id.to_string(),
            chain_class: ChainClass::Unknown,
            confidence: Confidence::None,
            needs_confirmation: true,
        };
    }

    let index = synonym_index(profile);

    // 3. Find the chain token: prefer the rightmost token that exactly equals a
    //    synonym, longest-match-first. The chain suffix is the common convention.
    for (ti, token) in tokens.iter().enumerate().rev() {
        let tl = token.to_ascii_lowercase();
        if let Some((_, class_key, single)) = index.iter().find(|(syn, _, _)| *syn == tl) {
            let chain_class = class_of(class_key);
            let ab_id = reconstruct_ab_id(record_id, &tokens, ti, profile);
            let (confidence, needs) = if *single {
                // Single-letter tokens are dangerous; always route to confirmation.
                (Confidence::Low, true)
            } else {
                (Confidence::High, false)
            };
            return NameParse {
                ab_id,
                chain_class,
                confidence,
                needs_confirmation: needs,
            };
        }
    }

    // 5. Fallback to the profile regex (lightweight: anchored suffix capture).
    if let Some(rx) = &profile.id_regex {
        if let Some(parse) = try_regex(record_id, rx) {
            return parse;
        }
    }

    // 6. Nothing matched.
    NameParse {
        ab_id: record_id.to_string(),
        chain_class: ChainClass::Unknown,
        confidence: Confidence::None,
        needs_confirmation: true,
    }
}

/// Reconstruct `ab_id` by removing the matched chain token (and the separators it
/// consumed) from the **original** id, preserving internal separators and case.
fn reconstruct_ab_id(
    original: &str,
    tokens: &[&str],
    matched_idx: usize,
    _profile: &NamingProfile,
) -> String {
    // If the matched token is the final token (the usual case), strip it as a
    // suffix from the original, then trim trailing separators.
    if matched_idx == tokens.len() - 1 {
        let token = tokens[matched_idx];
        if let Some(pos) = original
            .to_ascii_lowercase()
            .rfind(&token.to_ascii_lowercase())
        {
            // Only treat as a suffix strip if the token is genuinely at the end.
            if pos + token.len() == original.len() {
                let head = &original[..pos];
                return head.trim_end_matches(['_', '-', ' ', ':']).to_string();
            }
        }
    }
    // Otherwise rebuild from the surviving tokens (rare; internal chain token).
    let kept: Vec<&str> = tokens
        .iter()
        .enumerate()
        .filter(|(i, _)| *i != matched_idx)
        .map(|(_, t)| *t)
        .collect();
    kept.join("_")
}

/// Minimal regex fallback. Rather than pull a regex engine into `core`, we honor
/// the *intent* of the default `id_regex`: an `ab_id` followed by a trailing
/// chain token among a fixed set. This keeps `core` dependency-light while still
/// catching ids whose chain token was not in the synonym table.
fn try_regex(record_id: &str, _rx: &str) -> Option<NameParse> {
    // Recognized trailing chain tokens (longest-first).
    const SUFFIXES: &[(&str, ChainClass)] = &[
        ("heavychain", ChainClass::Heavy),
        ("lightchain", ChainClass::Light),
        ("lambdachain", ChainClass::Lambda),
        ("kappachain", ChainClass::Kappa),
        ("heavy", ChainClass::Heavy),
        ("light", ChainClass::Light),
        ("lambda", ChainClass::Lambda),
        ("kappa", ChainClass::Kappa),
        ("hchain", ChainClass::Heavy),
        ("lchain", ChainClass::Light),
        ("hc", ChainClass::Heavy),
        ("lc", ChainClass::Light),
        ("vh", ChainClass::Heavy),
        ("vk", ChainClass::Kappa),
        ("vl", ChainClass::Lambda),
    ];
    let lower = record_id.to_ascii_lowercase();
    for (suffix, class) in SUFFIXES {
        if lower.ends_with(suffix) {
            let head = &record_id[..record_id.len() - suffix.len()];
            let ab_id = head.trim_end_matches(['_', '-', ' ', ':']).to_string();
            if !ab_id.is_empty() {
                return Some(NameParse {
                    ab_id,
                    chain_class: *class,
                    confidence: Confidence::High,
                    needs_confirmation: false,
                });
            }
        }
    }
    None
}

/// The built-in default naming profile (mirrors `reference/example_library.json5`).
/// Used as a fallback when a project references a profile the library lacks.
pub fn default_profile() -> NamingProfile {
    let mut chain_synonyms = BTreeMap::new();
    chain_synonyms.insert(
        "heavy".into(),
        ["heavychain", "heavy", "hchain", "hc", "vh", "igh", "h"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    chain_synonyms.insert(
        "kappa".into(),
        ["kappachain", "kappa", "igk", "vk", "k"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    chain_synonyms.insert(
        "lambda".into(),
        ["lambdachain", "lambda", "igl", "vl", "l"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    chain_synonyms.insert(
        "light".into(),
        ["lightchain", "light", "lchain", "lc"]
            .iter()
            .map(|s| s.to_string())
            .collect(),
    );
    NamingProfile {
        id: "default".into(),
        display_name: "Default".into(),
        chain_synonyms,
        id_regex: Some(
            r"^(?<ab_id>.+?)[_\-\s]*(?<chain>heavy|light|hc|lc|kappa|lambda|[hkl])$".into(),
        ),
        separators: vec!["_".into(), "-".into(), " ".into(), ":".into()],
    }
}

/// Group parsed records by `ab_id` for per-antibody rollup (`docs/04` §4).
pub fn group_by_ab(parses: &[(String, NameParse)]) -> BTreeMap<String, Vec<String>> {
    let mut groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for (record_id, p) in parses {
        groups
            .entry(p.ab_id.clone())
            .or_default()
            .push(record_id.clone());
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    fn default_profile() -> NamingProfile {
        let mut chain_synonyms = BTreeMap::new();
        chain_synonyms.insert(
            "heavy".to_string(),
            vec!["heavychain", "heavy", "hchain", "hc", "vh", "igh", "h"]
                .into_iter()
                .map(String::from)
                .collect(),
        );
        chain_synonyms.insert(
            "kappa".to_string(),
            vec!["kappachain", "kappa", "igk", "vk", "k"]
                .into_iter()
                .map(String::from)
                .collect(),
        );
        chain_synonyms.insert(
            "lambda".to_string(),
            vec!["lambdachain", "lambda", "igl", "vl", "l"]
                .into_iter()
                .map(String::from)
                .collect(),
        );
        chain_synonyms.insert(
            "light".to_string(),
            vec!["lightchain", "light", "lchain", "lc"]
                .into_iter()
                .map(String::from)
                .collect(),
        );
        NamingProfile {
            id: "default".into(),
            display_name: "Default".into(),
            chain_synonyms,
            id_regex: Some(
                r"^(?<ab_id>.+?)[_\-\s]*(?<chain>heavy|light|hc|lc|kappa|lambda|[hkl])$".into(),
            ),
            separators: vec!["_".into(), "-".into(), " ".into(), ":".into()],
        }
    }

    fn p(id: &str) -> NameParse {
        parse_name(id, &default_profile())
    }

    #[test]
    fn worked_examples_from_doc_04() {
        let a = p("HVA-0195-r3-d02_heavy");
        assert_eq!(a.ab_id, "HVA-0195-r3-d02");
        assert_eq!(a.chain_class, ChainClass::Heavy);
        assert_eq!(a.confidence, Confidence::High);

        let b = p("HVA-0195-r3-d02_light");
        assert_eq!(b.ab_id, "HVA-0195-r3-d02");
        assert_eq!(b.chain_class, ChainClass::Light);

        // Colon + double underscore must not corrupt ab_id.
        let c = p("UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4__heavy");
        assert_eq!(c.ab_id, "UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4");
        assert_eq!(c.chain_class, ChainClass::Heavy);

        let d = p("mab1_HC");
        assert_eq!(d.ab_id, "mab1");
        assert_eq!(d.chain_class, ChainClass::Heavy);

        let e = p("Ab_007_lambda");
        assert_eq!(e.ab_id, "Ab_007");
        assert_eq!(e.chain_class, ChainClass::Lambda);

        let f = p("clone3-heavychain");
        assert_eq!(f.ab_id, "clone3");
        assert_eq!(f.chain_class, ChainClass::Heavy);

        let g = p("7G12_kappa");
        assert_eq!(g.ab_id, "7G12");
        assert_eq!(g.chain_class, ChainClass::Kappa);
    }

    #[test]
    fn single_letter_routes_to_confirmation() {
        let s = p("sample12_H");
        assert_eq!(s.ab_id, "sample12");
        assert_eq!(s.chain_class, ChainClass::Heavy);
        assert_eq!(s.confidence, Confidence::Low);
        assert!(s.needs_confirmation);
    }

    #[test]
    fn no_token_is_unknown_and_needs_confirmation() {
        let w = p("weird_name_no_token");
        assert_eq!(w.ab_id, "weird_name_no_token");
        assert_eq!(w.chain_class, ChainClass::Unknown);
        assert!(w.needs_confirmation);
    }

    #[test]
    fn deterministic() {
        assert_eq!(p("HVA-0195-r3-d02_heavy"), p("HVA-0195-r3-d02_heavy"));
    }
}
