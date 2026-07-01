//! Verdict taxonomy shared by both gates (`docs/01_DESIGN.md` §2, §3).
//!
//! Every verdict carries a human-readable `reason` a bench scientist can act on
//! without reading code, plus structured fields for the report tables.

use serde::{Deserialize, Serialize};

/// Severity ordering so the "worst wins" when several signals apply
/// (`docs/01_DESIGN.md` §3.4). Higher discriminant = worse.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub enum Severity {
    /// Informational, never blocks a PASS.
    Advisory,
    /// Clean pass.
    Pass,
    /// A real problem the scientist must look at.
    Fail,
}

/// Gate-1 (pre-cloning order QC) verdict kinds.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Gate1Kind {
    Pass,
    TranslationDrift,
    FrameshiftAtJunction,
    OverhangMissing,
    OverhangWrongLocus,
    LocusMismatch,
    NoGroundTruth,
    GcWarning,
    RareCodonWarning,
}

impl Gate1Kind {
    pub fn severity(self) -> Severity {
        use Gate1Kind::*;
        match self {
            Pass => Severity::Pass,
            NoGroundTruth | GcWarning | RareCodonWarning => Severity::Advisory,
            _ => Severity::Fail,
        }
    }

    pub fn label(self) -> &'static str {
        use Gate1Kind::*;
        match self {
            Pass => "PASS",
            TranslationDrift => "TRANSLATION_DRIFT",
            FrameshiftAtJunction => "FRAMESHIFT_AT_JUNCTION",
            OverhangMissing => "OVERHANG_MISSING",
            OverhangWrongLocus => "OVERHANG_WRONG_LOCUS",
            LocusMismatch => "LOCUS_MISMATCH",
            NoGroundTruth => "NO_GROUND_TRUTH",
            GcWarning => "GC_WARNING",
            RareCodonWarning => "RARE_CODON_WARNING",
        }
    }
}

/// A single ordered record's Gate-1 result.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gate1Verdict {
    pub record_id: String,
    pub ab_id: String,
    pub chain_class: String,
    /// The primary (worst) verdict kind.
    pub kind: Gate1Kind,
    /// Any advisory flags raised alongside the primary verdict.
    pub advisories: Vec<Gate1Kind>,
    pub reason: String,
    /// Core (insert) length in nt, when extracted.
    pub core_len: Option<usize>,
    /// Index of the first premature stop codon (aa), when relevant.
    pub premature_stop_aa: Option<usize>,
    /// Whether the assembled ORF reads through into the constant region.
    pub reads_through: Option<bool>,
}

impl Gate1Verdict {
    pub fn severity(&self) -> Severity {
        self.kind.severity()
    }
    pub fn passed(&self) -> bool {
        matches!(self.kind, Gate1Kind::Pass | Gate1Kind::NoGroundTruth)
    }
}

/// Gate-2 (post-cloning sequencing QC) verdict kinds, in worst-wins order.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum Gate2Kind {
    WrongVector,
    EmptyVector,
    JunctionFrameshift,
    WrongInsertSwap,
    PointMutation,
    SilentVariant,
    InsufficientRead,
    Pass,
}

impl Gate2Kind {
    /// Worst-wins rank (higher = worse), matching `docs/01_DESIGN.md` §3.4:
    /// `WRONG_VECTOR > EMPTY_VECTOR > JUNCTION_FRAMESHIFT > WRONG_INSERT_SWAP >
    ///  POINT_MUTATION > SILENT_VARIANT > PASS`. `InsufficientRead` is a
    /// non-determination ranked just above PASS.
    pub fn rank(self) -> u8 {
        use Gate2Kind::*;
        match self {
            WrongVector => 7,
            EmptyVector => 6,
            JunctionFrameshift => 5,
            WrongInsertSwap => 4,
            PointMutation => 3,
            SilentVariant => 2,
            InsufficientRead => 1,
            Pass => 0,
        }
    }

    pub fn severity(self) -> Severity {
        use Gate2Kind::*;
        match self {
            Pass => Severity::Pass,
            SilentVariant => Severity::Advisory,
            _ => Severity::Fail,
        }
    }

    pub fn label(self) -> &'static str {
        use Gate2Kind::*;
        match self {
            WrongVector => "WRONG_VECTOR",
            EmptyVector => "EMPTY_VECTOR",
            JunctionFrameshift => "JUNCTION_FRAMESHIFT",
            WrongInsertSwap => "WRONG_INSERT_SWAP",
            PointMutation => "POINT_MUTATION",
            SilentVariant => "SILENT_VARIANT",
            InsufficientRead => "INSUFFICIENT_READ",
            Pass => "PASS",
        }
    }
}

/// A point substitution located in the protein alignment.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mutation {
    pub position_aa: usize,
    pub wt: char,
    pub mut_aa: char,
}

/// A single sequencing read's Gate-2 result for one chain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Gate2Verdict {
    pub record_id: String,
    pub ab_id: String,
    pub chain_class: String,
    pub kind: Gate2Kind,
    pub reason: String,
    /// Autodetected backbone (vector id) and its percent identity, if observed.
    pub backbone_vector: Option<String>,
    pub backbone_identity: Option<f64>,
    /// "full" or "flanks_only" (partial Sanger).
    pub backbone_observed: String,
    pub mutations: Vec<Mutation>,
    pub premature_stop_aa: Option<usize>,
    pub reads_through: Option<bool>,
    /// If a sample swap is suspected, the better-matching panel member.
    pub suspected_identity: Option<String>,
    /// The gapped protein alignment (observed ORF vs expected), when computed —
    /// drives the Wetlab alignment view. `aligned_observed`/`aligned_expected`
    /// are equal-length strings using `-` for gaps.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned_observed: Option<String>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub aligned_expected: Option<String>,
    /// Residue index (0-based, into the aligned strings) of the leader→V boundary,
    /// so the UI can label variable-region positions.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub leader_aa_len: Option<usize>,
}

impl Gate2Verdict {
    pub fn severity(&self) -> Severity {
        self.kind.severity()
    }
    pub fn passed(&self) -> bool {
        matches!(self.kind, Gate2Kind::Pass | Gate2Kind::SilentVariant)
    }
}

/// Per-antibody rollup: an antibody PASSES only if every expected chain passes.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct AntibodyRollup {
    pub ab_id: String,
    pub passed: bool,
    pub note: String,
    /// record_ids of the chains that contributed.
    pub chains: Vec<String>,
}
