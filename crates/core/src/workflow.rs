//! High-level orchestration helpers shared by the CLI and GUI: building an
//! ad-hoc project from a library and extracting expected insert cores from an
//! order file. Keeps both front-ends thin and consistent.

use crate::assemble;
use crate::model::{
    sha256_hex, ChainClass, InsertionSite, Library, Locus, OverhangSet, Project, Provenance, Vector,
};
use crate::naming;
use crate::seq;
use crate::seqio::genbank::GbRecord;
use crate::seqio::SeqRecord;
use std::collections::BTreeMap;

/// Build a [`Vector`] from a parsed GenBank record for the guided library
/// importer (first-use library construction). The insertion site is located by
/// finding the chosen overhang set's 5′/3′ arms in the vector sequence; if they
/// are not present (different overhangs), it falls back to the signal-peptide /
/// constant-region feature boundaries. The constant-region protein anchor is
/// read off the in-frame translation just past the insertion site.
#[allow(clippy::too_many_arguments)]
pub fn vector_from_genbank(
    gb: &GbRecord,
    id: &str,
    display_name: &str,
    chain_class: ChainClass,
    isotype: &str,
    set: &OverhangSet,
    locus: Locus,
    added_by: &str,
) -> Vector {
    let sequence = seq::clean(&gb.sequence);
    let features = crate::seqio::genbank::map_roles(gb);

    // Locate overhangs in the vector to compute the insertion site.
    let oh5 = set.oh5(locus).unwrap_or("");
    let oh3 = set.oh3(locus).unwrap_or("");
    let oh5_end = match (!oh5.is_empty()).then(|| sequence.find(oh5)).flatten() {
        Some(p) => p + oh5.len(),
        None => features
            .get("signal_peptide")
            .map(|f| f.range0().end)
            .unwrap_or(0),
    };
    let oh3_start = match (!oh3.is_empty()).then(|| sequence.find(oh3)).flatten() {
        Some(p) => p,
        None => features
            .get("constant_region")
            .map(|f| f.range0().start)
            .unwrap_or(oh5_end),
    };

    // Constant-region anchor: first residues of the in-frame translation
    // beginning at the insertion site (codon-invariant backbone).
    let anchor: String = seq::translate_to_stop(&sequence[oh3_start.min(sequence.len())..])
        .chars()
        .take(8)
        .collect();

    let length = sequence.len();
    let sha = sha256_hex(&sequence);
    Vector {
        id: id.to_string(),
        display_name: display_name.to_string(),
        chain_class,
        isotype: isotype.to_string(),
        topology: if gb.circular { "circular" } else { "linear" }.to_string(),
        length,
        sequence,
        sequence_sha256: sha,
        features,
        insertion_site: InsertionSite { oh5_end, oh3_start },
        constant_anchor_aa: anchor,
        overhang_set: set.id.clone(),
        overhang_locus: locus,
        provenance: Provenance {
            added_by: added_by.to_string(),
            added: String::new(),
            source_file: Some(gb.name.clone()),
        },
    }
}

/// Build an in-memory project: assign the chosen (or first) heavy/κ/λ vectors,
/// the chosen (or first) overhang set, and the first naming profile.
pub fn ad_hoc_project(
    lib: &Library,
    heavy_vector: Option<&str>,
    overhang_set: Option<&str>,
) -> Project {
    let overhang_set = overhang_set
        .map(|s| s.to_string())
        .or_else(|| lib.overhang_sets.first().map(|o| o.id.clone()))
        .unwrap_or_default();
    let naming_profile = lib
        .naming_profiles
        .first()
        .map(|p| p.id.clone())
        .unwrap_or_else(|| "default".to_string());

    let mut vector_assignments = BTreeMap::new();
    let heavy = heavy_vector.map(|s| s.to_string()).or_else(|| {
        lib.vectors
            .iter()
            .find(|v| v.chain_class == ChainClass::Heavy)
            .map(|v| v.id.clone())
    });
    if let Some(h) = heavy {
        vector_assignments.insert("heavy".to_string(), h);
    }
    if let Some(k) = lib
        .vectors
        .iter()
        .find(|v| v.chain_class == ChainClass::Kappa)
    {
        vector_assignments.insert("kappa".to_string(), k.id.clone());
    }
    if let Some(l) = lib
        .vectors
        .iter()
        .find(|v| v.chain_class == ChainClass::Lambda)
    {
        vector_assignments.insert("lambda".to_string(), l.id.clone());
    }

    Project {
        schema_version: 1,
        name: "project".into(),
        created: String::new(),
        description: String::new(),
        naming_profile,
        overhang_set,
        vector_assignments,
        ground_truth: None,
        order: None,
        sequencing: None,
        name_overrides: BTreeMap::new(),
    }
}

/// Extract expected insert cores from an order file, keyed by `(AB_ID_UPPER, class)`.
/// A `light` record's κ/λ locus is resolved from the fragment's 3′ overhang.
pub fn order_cores(
    records: &[SeqRecord],
    lib: &Library,
    project: &Project,
    has_overhangs: bool,
) -> BTreeMap<(String, ChainClass), String> {
    let Some(set) = lib.overhang_set(&project.overhang_set) else {
        return BTreeMap::new();
    };
    let prof = lib
        .naming_profile(&project.naming_profile)
        .cloned()
        .unwrap_or_else(naming::default_profile);
    let mm = lib.alignment.overhang_max_mismatch;
    let mut map = BTreeMap::new();
    for r in records {
        let np = naming::parse_name(&r.id, &prof);
        let locus = match np.chain_class {
            ChainClass::Heavy => Some(Locus::Igh),
            ChainClass::Kappa => Some(Locus::Igk),
            ChainClass::Lambda => Some(Locus::Igl),
            ChainClass::Light => assemble::detect_locus_from_overhangs(&r.sequence, set, mm),
            ChainClass::Unknown => None,
        };
        let Some(locus) = locus else { continue };
        let core = if has_overhangs {
            let s = assemble::strip_overhangs(
                &r.sequence,
                set.oh5(locus).unwrap_or(""),
                set.oh3(locus).unwrap_or(""),
                mm,
            );
            if s.oh5_present && s.oh3_present {
                s.core
            } else {
                continue;
            }
        } else {
            seq::clean(&r.sequence)
        };
        map.insert((np.ab_id.to_ascii_uppercase(), locus.chain_class()), core);
    }
    map
}
