//! High-level orchestration helpers shared by the CLI and GUI: building an
//! ad-hoc project from a library and extracting expected insert cores from an
//! order file. Keeps both front-ends thin and consistent.

use crate::assemble;
use crate::model::{ChainClass, Library, Locus, Project};
use crate::naming;
use crate::seq;
use crate::seqio::SeqRecord;
use std::collections::BTreeMap;

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
