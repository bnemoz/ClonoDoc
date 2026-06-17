//! In-silico assembly, overhang handling, locus detection, and circular
//! topology normalization (`docs/01_DESIGN.md` §1.5, §3.1, §3.2).

use crate::align::best_ungapped_identity;
use crate::model::{Locus, OverhangSet, Vector};
use crate::seq;

/// Result of stripping overhangs from an ordered fragment.
#[derive(Debug, Clone)]
pub struct Stripped {
    pub core: String,
    pub oh5_mismatches: usize,
    pub oh3_mismatches: usize,
    pub oh5_present: bool,
    pub oh3_present: bool,
}

/// Count mismatches between two equal-length byte slices.
fn mismatches(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b.iter()).filter(|(x, y)| x != y).count()
}

/// Strip the 5′ and 3′ overhangs from an ordered fragment, tolerating up to
/// `max_mismatch` typos per overhang (but reporting them). Returns the core
/// (the insert) plus which overhangs were found.
pub fn strip_overhangs(ordered: &str, oh5: &str, oh3: &str, max_mismatch: usize) -> Stripped {
    let s = seq::clean(ordered);
    let bytes = s.as_bytes();
    let o5 = oh5.as_bytes();
    let o3 = oh3.as_bytes();

    let (oh5_present, oh5_mm, start) = if bytes.len() >= o5.len() {
        let mm = mismatches(&bytes[..o5.len()], o5);
        (mm <= max_mismatch, mm, o5.len())
    } else {
        (false, o5.len(), 0)
    };

    let (oh3_present, oh3_mm, end) = if bytes.len() >= o3.len() {
        let tail = &bytes[bytes.len() - o3.len()..];
        let mm = mismatches(tail, o3);
        (mm <= max_mismatch, mm, bytes.len() - o3.len())
    } else {
        (false, o3.len(), bytes.len())
    };

    // Only carve the core out of the region actually flanked by present overhangs.
    let lo = if oh5_present { start } else { 0 };
    let hi = if oh3_present { end } else { bytes.len() };
    let core = if lo <= hi {
        s[lo..hi].to_string()
    } else {
        String::new()
    };

    Stripped {
        core,
        oh5_mismatches: oh5_mm,
        oh3_mismatches: oh3_mm,
        oh5_present,
        oh3_present,
    }
}

/// Build the assembled construct (`docs/01_DESIGN.md` §3.2):
/// `vector[0..oh5_end] ++ core ++ vector[oh3_start..]`.
pub fn assemble(vector: &Vector, core: &str) -> String {
    let v = seq::clean(&vector.sequence);
    let oh5_end = vector.insertion_site.oh5_end.min(v.len());
    let oh3_start = vector.insertion_site.oh3_start.min(v.len());
    let core = seq::clean(core);
    let mut out = String::with_capacity(v.len() + core.len());
    out.push_str(&v[..oh5_end]);
    out.push_str(&core);
    out.push_str(&v[oh3_start..]);
    out
}

/// The half-open insert window `[oh5_end, oh5_end + core_len)` within the
/// assembled construct.
pub fn insert_window(vector: &Vector, core_len: usize) -> std::ops::Range<usize> {
    let start = vector.insertion_site.oh5_end;
    start..(start + core_len)
}

/// Detect the locus a fragment carries from its overhangs.
///
/// Returns `(locus, ambiguous)`. Because `oh_5[IGK] == oh_5[IGL]`, the 5′
/// overhang alone cannot split κ from λ — the 3′ overhang is the discriminator
/// (`docs/01_DESIGN.md` §1.5). Heavy is unambiguous from either end.
pub fn detect_locus_from_overhangs(
    fragment: &str,
    set: &OverhangSet,
    max_mismatch: usize,
) -> Option<Locus> {
    let s = seq::clean(fragment);
    let bytes = s.as_bytes();
    let mut candidates = Vec::new();
    for locus in [Locus::Igh, Locus::Igk, Locus::Igl] {
        let oh5 = set.oh5(locus);
        let oh3 = set.oh3(locus);
        let oh5_ok = oh5
            .map(|o| {
                let o = o.as_bytes();
                bytes.len() >= o.len() && mismatches(&bytes[..o.len()], o) <= max_mismatch
            })
            .unwrap_or(false);
        let oh3_ok = oh3
            .map(|o| {
                let o = o.as_bytes();
                bytes.len() >= o.len()
                    && mismatches(&bytes[bytes.len() - o.len()..], o) <= max_mismatch
            })
            .unwrap_or(false);
        // The 3′ overhang is the authoritative discriminator; require it when
        // available, falling back to the 5′ match only for a partial fragment.
        if oh3_ok {
            candidates.push((locus, 2));
        } else if oh5_ok {
            candidates.push((locus, 1));
        }
    }
    // Prefer the strongest (3′-confirmed) candidate; ties resolved deterministically.
    candidates.sort_by(|a, b| b.1.cmp(&a.1).then(locus_rank(a.0).cmp(&locus_rank(b.0))));
    candidates.first().map(|(l, _)| *l)
}

fn locus_rank(l: Locus) -> u8 {
    match l {
        Locus::Igh => 0,
        Locus::Igk => 1,
        Locus::Igl => 2,
    }
}

/// Detect locus from a translated protein's constant-region motif
/// (`docs/01_DESIGN.md` §1.5). Backbone-independent confirmation for partial reads.
pub fn detect_locus_from_protein(protein: &str) -> Option<Locus> {
    // Heavy: V ends `WG.GT[LT]VTVSS` then constant starts `ASTKGP`.
    if protein.contains("ASTKGP") || protein.contains("VTVSS") {
        return Some(Locus::Igh);
    }
    // Kappa: V ends `FG.GTKVEIK`, constant starts `RTVAAP`.
    if protein.contains("TKVEIK") || protein.contains("RTVAAP") {
        return Some(Locus::Igk);
    }
    // Lambda: constant starts `GQPKAAP`; V often ends `...TVL`.
    if protein.contains("GQPKAAP") || protein.contains("QPKAAPSV") {
        return Some(Locus::Igl);
    }
    None
}

/// Orientation/rotation outcome of normalizing a circular read.
#[derive(Debug, Clone)]
pub struct Normalized {
    pub sequence: String,
    pub reverse_complemented: bool,
    pub rotation: usize,
}

/// Normalize a circular read against a reference (`docs/01_DESIGN.md` §3.1):
/// pick the strand whose best ungapped match to a `landmark` is higher, then
/// rotate so the landmark sits near the start. `landmark` is a codon-opt-invariant
/// backbone anchor (e.g. the start of the constant region or the leader).
pub fn normalize_circular(read: &str, landmark: &str) -> Normalized {
    let fwd = seq::clean(read);
    let rev = seq::revcomp(&fwd);

    let (chosen, rced) = {
        let f = best_landmark_identity(&fwd, landmark);
        let r = best_landmark_identity(&rev, landmark);
        if r > f {
            (rev, true)
        } else {
            (fwd, false)
        }
    };

    // Rotate so the landmark's best match is at index 0 (treat as circular).
    let rotation = best_landmark_offset(&chosen, landmark).unwrap_or(0);
    let rotated = rotate(&chosen, rotation);

    Normalized {
        sequence: rotated,
        reverse_complemented: rced,
        rotation,
    }
}

/// Best ungapped identity of `landmark` against any circular offset of `text`.
fn best_landmark_identity(text: &str, landmark: &str) -> f64 {
    if landmark.is_empty() || text.is_empty() {
        return 0.0;
    }
    // Double the text to emulate circularity, then slide the (shorter) landmark.
    let doubled = format!("{text}{text}");
    let cap = text.len(); // only need offsets within one full turn
    let t = doubled.as_bytes();
    let l = landmark.as_bytes();
    if l.len() > text.len() {
        return best_ungapped_identity(text, landmark);
    }
    let mut best = 0usize;
    for off in 0..cap {
        let mut hit = 0;
        for k in 0..l.len() {
            if t[off + k] == l[k] {
                hit += 1;
            }
        }
        best = best.max(hit);
    }
    best as f64 / l.len() as f64
}

/// Circular offset in `text` where `landmark` matches best.
fn best_landmark_offset(text: &str, landmark: &str) -> Option<usize> {
    if landmark.is_empty() || text.is_empty() || landmark.len() > text.len() {
        return None;
    }
    let doubled = format!("{text}{text}");
    let t = doubled.as_bytes();
    let l = landmark.as_bytes();
    let mut best = (0usize, 0usize); // (hits, offset)
    for off in 0..text.len() {
        let mut hit = 0;
        for k in 0..l.len() {
            if t[off + k] == l[k] {
                hit += 1;
            }
        }
        if hit > best.0 {
            best = (hit, off);
        }
    }
    Some(best.1)
}

/// Rotate a string left by `n` (circular).
fn rotate(s: &str, n: usize) -> String {
    if s.is_empty() {
        return String::new();
    }
    let n = n % s.len();
    let b = s.as_bytes();
    let mut out = Vec::with_capacity(s.len());
    out.extend_from_slice(&b[n..]);
    out.extend_from_slice(&b[..n]);
    String::from_utf8(out).unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::{ChainClass, Feature, InsertionSite, Provenance};
    use std::collections::BTreeMap;

    fn french_overhangs() -> OverhangSet {
        let mut oh_5 = BTreeMap::new();
        oh_5.insert(
            "IGH".into(),
            "GCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGT".into(),
        );
        oh_5.insert("IGK".into(), "ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC".into());
        oh_5.insert("IGL".into(), "ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC".into());
        let mut oh_3 = BTreeMap::new();
        oh_3.insert("IGH".into(), "GCTAGCACCAAGGGCCCATCGGTCTTCC".into());
        oh_3.insert("IGK".into(), "CGTACGGTGGCTGCACCATCTGTCTTCATC".into());
        oh_3.insert(
            "IGL".into(),
            "GGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCGAGTGAGGAGCTTCAAGCCAACAAGGCC".into(),
        );
        OverhangSet {
            id: "french_default".into(),
            display_name: "French".into(),
            oh_5,
            oh_3,
            provenance: Provenance::default(),
        }
    }

    #[test]
    fn igh_overhang_vector_mapping() {
        // From verified_facts §2: oh_5[IGH] == vector 17..57 (0-based 16..57);
        // oh_3[IGH] == vector 58..85 (0-based 57..85). These are checked against
        // the real vector in the integration tests; here we check the lengths.
        let set = french_overhangs();
        assert_eq!(set.oh5(Locus::Igh).unwrap().len(), 41);
        assert_eq!(set.oh3(Locus::Igh).unwrap().len(), 28);
        // 5' shared between light loci; 3' distinct.
        assert_eq!(set.oh5(Locus::Igk), set.oh5(Locus::Igl));
        assert_ne!(set.oh3(Locus::Igk), set.oh3(Locus::Igl));
    }

    #[test]
    fn locus_from_overhang_uses_3prime_for_light() {
        let set = french_overhangs();
        // Build a fake kappa fragment: oh5(IGK) + filler + oh3(IGK).
        let frag = format!(
            "{}AAAAAA{}",
            set.oh5(Locus::Igk).unwrap(),
            set.oh3(Locus::Igk).unwrap()
        );
        assert_eq!(
            detect_locus_from_overhangs(&frag, &set, 1),
            Some(Locus::Igk)
        );
        let frag_l = format!(
            "{}AAAAAA{}",
            set.oh5(Locus::Igl).unwrap(),
            set.oh3(Locus::Igl).unwrap()
        );
        assert_eq!(
            detect_locus_from_overhangs(&frag_l, &set, 1),
            Some(Locus::Igl)
        );
    }

    fn tiny_vector() -> Vector {
        Vector {
            id: "t".into(),
            display_name: "t".into(),
            chain_class: ChainClass::Heavy,
            isotype: "IgG1".into(),
            topology: "circular".into(),
            length: 8,
            sequence: "AAAACCCC".into(),
            sequence_sha256: String::new(),
            features: BTreeMap::<String, Feature>::new(),
            insertion_site: InsertionSite {
                oh5_end: 4,
                oh3_start: 4,
            },
            constant_anchor_aa: "P".into(),
            overhang_set: "x".into(),
            overhang_locus: Locus::Igh,
            provenance: Provenance::default(),
        }
    }

    #[test]
    fn assembly_merges_core_between_flanks() {
        let v = tiny_vector();
        assert_eq!(assemble(&v, "GG"), "AAAAGGCCCC");
        assert_eq!(insert_window(&v, 2), 4..6);
    }

    #[test]
    fn circular_rotation_and_strand() {
        // A circular text; the read is a rotation of it, reverse-complemented.
        let text = "ATGCATGCAAAATTTTGGGGCCCC";
        let landmark = "AAAATTTTGGGG";
        let read = rotate(text, 7);
        let rc = seq::revcomp(&read);
        let norm = normalize_circular(&rc, landmark);
        assert!(norm.reverse_complemented);
        assert!(norm.sequence.starts_with(landmark));
    }
}
