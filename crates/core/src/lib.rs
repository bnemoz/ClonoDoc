//! abclone-core — the headless verification engine for `abclone-verify`.
//!
//! This crate contains **zero** UI dependencies. Every verdict the tool can
//! produce is computed here and exercised by the integration tests against the
//! fixtures in `reference/` and `test_data/`. The GUI (`abclone-app`) and the
//! CLI (`abclone-cli`) are thin orchestration layers over this engine.
//!
//! Reading order mirrors `docs/`:
//! * [`seq`]    — sequence primitives (alphabet detection, translation).
//! * [`seqio`]  — parsers (FASTA, CSV, XLSX, GenBank, AB1).
//! * [`model`]  — serde data model (Library, Vector, OverhangSet, Project).
//! * [`naming`] — sample-name parsing & chain pairing (`docs/04_NAMING.md`).
//! * [`align`]  — affine-gap pairwise alignment (nt + BLOSUM62 protein).
//! * [`assemble`] — in-silico construct assembly (`docs/01_DESIGN.md` §3.2).
//! * [`gate1`]  — pre-cloning order QC.
//! * [`gate2`]  — post-cloning sequencing QC.
//! * [`report`] — verdict export (CSV / XLSX-ready / HTML).

pub mod align;
pub mod assemble;
pub mod gate1;
pub mod gate2;
pub mod model;
pub mod naming;
pub mod report;
pub mod seq;
pub mod seqio;
pub mod verdict;

pub use verdict::{Gate1Verdict, Gate2Verdict, Severity};

/// Crate-wide error type. Parsers and engines return this; the apps wrap it.
#[derive(Debug, thiserror::Error)]
pub enum Error {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),
    #[error("parse error in {context}: {message}")]
    Parse { context: String, message: String },
    #[error("config error: {0}")]
    Config(String),
    #[error("{0}")]
    Other(String),
}

impl Error {
    pub fn parse(context: impl Into<String>, message: impl Into<String>) -> Self {
        Error::Parse {
            context: context.into(),
            message: message.into(),
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;
