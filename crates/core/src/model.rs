//! Serde data model for the lab-global **Library** and per-campaign **Project**
//! (`docs/02_CONFIG_SCHEMA.md`). Stored as `json5` so the shared library can carry
//! comments. Coordinates in the on-disk feature table are 1-based inclusive
//! (GenBank convention); [`Feature::range0`] converts to 0-based half-open once.

use crate::{Error, Result};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::Path;

/// Chain class as parsed from a sample name. `Light` is deliberately distinct
/// from κ/λ — many panels label light chains only as `light` and leave the locus
/// to downstream autodetect (`docs/04_NAMING.md` §2).
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum ChainClass {
    Heavy,
    Kappa,
    Lambda,
    Light,
    Unknown,
}

impl ChainClass {
    pub fn as_str(self) -> &'static str {
        match self {
            ChainClass::Heavy => "heavy",
            ChainClass::Kappa => "kappa",
            ChainClass::Lambda => "lambda",
            ChainClass::Light => "light",
            ChainClass::Unknown => "unknown",
        }
    }

    /// The immunoglobulin loci this class is compatible with.
    pub fn loci(self) -> &'static [Locus] {
        match self {
            ChainClass::Heavy => &[Locus::Igh],
            ChainClass::Kappa => &[Locus::Igk],
            ChainClass::Lambda => &[Locus::Igl],
            // A `light` token may resolve to κ or λ.
            ChainClass::Light => &[Locus::Igk, Locus::Igl],
            ChainClass::Unknown => &[],
        }
    }
}

/// An immunoglobulin locus / overhang key.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum Locus {
    #[serde(rename = "IGH")]
    Igh,
    #[serde(rename = "IGK")]
    Igk,
    #[serde(rename = "IGL")]
    Igl,
}

impl Locus {
    pub fn key(self) -> &'static str {
        match self {
            Locus::Igh => "IGH",
            Locus::Igk => "IGK",
            Locus::Igl => "IGL",
        }
    }
    pub fn chain_class(self) -> ChainClass {
        match self {
            Locus::Igh => ChainClass::Heavy,
            Locus::Igk => ChainClass::Kappa,
            Locus::Igl => ChainClass::Lambda,
        }
    }
}

/// A 1-based inclusive feature span (GenBank convention), with optional name.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Feature {
    pub start: usize,
    pub end: usize,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub name: Option<String>,
}

impl Feature {
    /// Convert to a 0-based half-open `[start, end)` Rust slice range.
    pub fn range0(&self) -> std::ops::Range<usize> {
        (self.start.saturating_sub(1))..self.end
    }
}

/// The insertion site: 0-based indices in the vector where the linearized ends
/// meet the insert. For the French IGH vector both equal 57 (overhangs adjacent).
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct InsertionSite {
    pub oh5_end: usize,
    pub oh3_start: usize,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct Provenance {
    #[serde(default)]
    pub added_by: String,
    #[serde(default)]
    pub added: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source_file: Option<String>,
}

/// A vector backbone.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Vector {
    pub id: String,
    pub display_name: String,
    pub chain_class: ChainClass,
    pub isotype: String,
    #[serde(default = "default_topology")]
    pub topology: String,
    #[serde(default)]
    pub length: usize,
    pub sequence: String,
    #[serde(default)]
    pub sequence_sha256: String,
    #[serde(default)]
    pub features: BTreeMap<String, Feature>,
    pub insertion_site: InsertionSite,
    pub constant_anchor_aa: String,
    pub overhang_set: String,
    pub overhang_locus: Locus,
    #[serde(default)]
    pub provenance: Provenance,
}

fn default_topology() -> String {
    "circular".to_string()
}

impl Vector {
    pub fn is_circular(&self) -> bool {
        self.topology.eq_ignore_ascii_case("circular")
    }
    /// Recompute and return the sha256 of the (cleaned) sequence.
    pub fn compute_sha256(&self) -> String {
        sha256_hex(&crate::seq::clean(&self.sequence))
    }
}

/// An overhang set: 5′ and 3′ homology arms keyed by locus.
/// `oh_5[IGK] == oh_5[IGL]` by design; `oh_3` distinguishes κ from λ.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OverhangSet {
    pub id: String,
    #[serde(default)]
    pub display_name: String,
    pub oh_5: BTreeMap<String, String>,
    pub oh_3: BTreeMap<String, String>,
    #[serde(default)]
    pub provenance: Provenance,
}

impl OverhangSet {
    pub fn oh5(&self, locus: Locus) -> Option<&str> {
        self.oh_5.get(locus.key()).map(|s| s.as_str())
    }
    pub fn oh3(&self, locus: Locus) -> Option<&str> {
        self.oh_3.get(locus.key()).map(|s| s.as_str())
    }
}

/// A naming profile: chain synonyms (per class), a fallback regex, and separators.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NamingProfile {
    pub id: String,
    #[serde(default)]
    pub display_name: String,
    pub chain_synonyms: BTreeMap<String, Vec<String>>,
    #[serde(default)]
    pub id_regex: Option<String>,
    #[serde(default = "default_separators")]
    pub separators: Vec<String>,
}

fn default_separators() -> Vec<String> {
    vec!["_".into(), "-".into(), " ".into(), ":".into()]
}

/// Alignment + threshold settings. No magic numbers in code — everything is here.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AlignmentSettings {
    #[serde(default = "default_matrix")]
    pub protein_matrix: String,
    pub protein_gap_open: i32,
    pub protein_gap_extend: i32,
    pub nt_match: i32,
    pub nt_mismatch: i32,
    pub nt_gap_open: i32,
    pub nt_gap_extend: i32,
    pub backbone_identity_min: f64,
    pub overhang_max_mismatch: usize,
}

fn default_matrix() -> String {
    "BLOSUM62".to_string()
}

impl Default for AlignmentSettings {
    fn default() -> Self {
        AlignmentSettings {
            protein_matrix: "BLOSUM62".into(),
            protein_gap_open: -11,
            protein_gap_extend: -1,
            nt_match: 2,
            nt_mismatch: -3,
            nt_gap_open: -5,
            nt_gap_extend: -1,
            backbone_identity_min: 0.97,
            overhang_max_mismatch: 1,
        }
    }
}

/// Advisory optimization checks (mirror dnachisel). Never fail an order.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OptimizationChecks {
    pub gc_max_global: f64,
    pub gc_max_window: f64,
    pub gc_window: usize,
    pub rare_codon_min_frequency: f64,
    #[serde(default)]
    pub species: String,
    #[serde(default)]
    pub enabled: bool,
}

impl Default for OptimizationChecks {
    fn default() -> Self {
        OptimizationChecks {
            gc_max_global: 0.56,
            gc_max_window: 0.64,
            gc_window: 60,
            rare_codon_min_frequency: 0.05,
            species: "h_sapiens".into(),
            enabled: true,
        }
    }
}

/// The lab-global library file.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Library {
    #[serde(default = "default_schema_version")]
    pub schema_version: u32,
    #[serde(default)]
    pub vectors: Vec<Vector>,
    #[serde(default)]
    pub overhang_sets: Vec<OverhangSet>,
    #[serde(default)]
    pub naming_profiles: Vec<NamingProfile>,
    #[serde(default)]
    pub alignment: AlignmentSettings,
    #[serde(default)]
    pub optimization_checks: OptimizationChecks,
}

fn default_schema_version() -> u32 {
    1
}

impl Library {
    pub fn load(path: impl AsRef<Path>) -> Result<Library> {
        let text = std::fs::read_to_string(&path)?;
        Library::from_json5(&text)
    }

    pub fn from_json5(text: &str) -> Result<Library> {
        json5::from_str(text).map_err(|e| Error::parse("library.json5", e.to_string()))
    }

    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        // Serialize as pretty JSON (valid json5) so it round-trips cleanly.
        let text = serde_json::to_string_pretty(self).map_err(|e| Error::Config(e.to_string()))?;
        std::fs::write(path, text)?;
        Ok(())
    }

    pub fn vector(&self, id: &str) -> Option<&Vector> {
        self.vectors.iter().find(|v| v.id == id)
    }
    pub fn overhang_set(&self, id: &str) -> Option<&OverhangSet> {
        self.overhang_sets.iter().find(|o| o.id == id)
    }
    pub fn naming_profile(&self, id: &str) -> Option<&NamingProfile> {
        self.naming_profiles.iter().find(|p| p.id == id)
    }

    /// Additive merge of another library (`docs/02` §1.4): new vectors / overhang
    /// sets are added; vectors are deduped by `sequence_sha256`. Returns a list of
    /// conflict messages for the caller to surface (never drops data silently).
    pub fn merge_from(&mut self, other: &Library) -> Vec<String> {
        let mut conflicts = Vec::new();
        let hash_of = |v: &Vector| {
            if v.sequence_sha256.is_empty() {
                v.compute_sha256()
            } else {
                v.sequence_sha256.clone()
            }
        };
        for v in &other.vectors {
            let hash = hash_of(v);
            match self.vectors.iter().find(|e| hash_of(e) == hash) {
                // Identical entry (same id + sequence): nothing to do.
                Some(e) if e.id == v.id => {}
                // Same sequence under a different id: keep both, surface a conflict.
                Some(e) => {
                    conflicts.push(format!(
                        "vector '{}' has the same sequence as existing '{}'; kept both",
                        v.id, e.id
                    ));
                    self.vectors.push(v.clone());
                }
                // New sequence: add it.
                None => self.vectors.push(v.clone()),
            }
        }
        for o in &other.overhang_sets {
            if !self.overhang_sets.iter().any(|e| e.id == o.id) {
                self.overhang_sets.push(o.clone());
            }
        }
        for p in &other.naming_profiles {
            if !self.naming_profiles.iter().any(|e| e.id == p.id) {
                self.naming_profiles.push(p.clone());
            }
        }
        conflicts
    }
}

// ---------------------------------------------------------------------------
// Project
// ---------------------------------------------------------------------------

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GroundTruthSpec {
    pub file: String,
    #[serde(default)]
    pub format: String,
    #[serde(default)]
    pub columns: BTreeMap<String, String>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrderSpec {
    pub file: String,
    #[serde(default)]
    pub format: String,
    #[serde(default = "default_true")]
    pub has_overhangs: bool,
}

fn default_true() -> bool {
    true
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SequencingInput {
    pub file: String,
    #[serde(default)]
    pub format: String,
    #[serde(default)]
    pub mode: Option<String>,
}

#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SequencingSpec {
    #[serde(default)]
    pub mode: String,
    #[serde(default)]
    pub inputs: Vec<SequencingInput>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NameOverride {
    pub ab_id: String,
    pub chain_class: ChainClass,
}

/// A per-campaign project.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Project {
    #[serde(default = "default_schema_version")]
    pub schema_version: u32,
    pub name: String,
    #[serde(default)]
    pub created: String,
    #[serde(default)]
    pub description: String,
    #[serde(default = "default_profile")]
    pub naming_profile: String,
    pub overhang_set: String,
    pub vector_assignments: BTreeMap<String, String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub ground_truth: Option<GroundTruthSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub order: Option<OrderSpec>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub sequencing: Option<SequencingSpec>,
    #[serde(default)]
    pub name_overrides: BTreeMap<String, NameOverride>,
}

fn default_profile() -> String {
    "default".to_string()
}

impl Project {
    pub fn load(path: impl AsRef<Path>) -> Result<Project> {
        let text = std::fs::read_to_string(&path)?;
        json5::from_str(&text).map_err(|e| Error::parse("project.json5", e.to_string()))
    }
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let text = serde_json::to_string_pretty(self).map_err(|e| Error::Config(e.to_string()))?;
        std::fs::write(path, text)?;
        Ok(())
    }

    /// The vector id assigned to a given chain class for this project.
    pub fn vector_for(&self, class: ChainClass) -> Option<&str> {
        self.vector_assignments
            .get(class.as_str())
            .map(|s| s.as_str())
    }
}

/// sha256 of a byte string as lowercase hex.
pub fn sha256_hex(data: &str) -> String {
    use sha2::{Digest, Sha256};
    let mut hasher = Sha256::new();
    hasher.update(data.as_bytes());
    let digest = hasher.finalize();
    digest.iter().map(|b| format!("{:02x}", b)).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn feature_range_is_zero_based_half_open() {
        let f = Feature {
            start: 58,
            end: 1047,
            name: None,
        };
        assert_eq!(f.range0(), 57..1047);
    }

    #[test]
    fn locus_class_mapping() {
        assert_eq!(Locus::Igh.chain_class(), ChainClass::Heavy);
        assert_eq!(Locus::Igk.chain_class(), ChainClass::Kappa);
        assert_eq!(ChainClass::Light.loci(), &[Locus::Igk, Locus::Igl]);
    }
}
