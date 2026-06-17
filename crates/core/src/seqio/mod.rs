//! Sequence & table parsers (`docs/01_DESIGN.md` §1.1, `docs/03_ARCHITECTURE.md` §4–5).
//!
//! Liberal in input by design: tolerate wrapped lines, lowercase, trailing
//! whitespace, gaps, and messy headers. Strict, explicit verdicts come later.

pub mod ab1;
pub mod fasta;
pub mod genbank;
pub mod tabular;

use crate::seq::{detect_type, SeqType};

/// A generic id + sequence record (a FASTA entry, an order row, a read).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqRecord {
    pub id: String,
    pub sequence: String,
}

impl SeqRecord {
    pub fn seq_type(&self) -> SeqType {
        detect_type(&self.sequence)
    }
}

/// A ground-truth panel row: an antibody id with its heavy and/or light chain.
#[derive(Debug, Clone, Default)]
pub struct GroundTruthRow {
    pub ab_id: String,
    pub heavy: Option<String>,
    pub light: Option<String>,
}

/// Dispatch a load by declared format string (case-insensitive). Falls back to
/// the file extension when `format` is empty.
pub fn read_id_seq_auto(path: &std::path::Path, format: &str) -> crate::Result<Vec<SeqRecord>> {
    let fmt = resolve_format(path, format);
    match fmt.as_str() {
        "fasta" => fasta::read_path(path),
        "xlsx" => tabular::read_xlsx_id_seq(path),
        "csv" => tabular::read_csv_id_seq(path),
        other => Err(crate::Error::parse(
            path.display().to_string(),
            format!("unsupported order format '{other}'"),
        )),
    }
}

pub(crate) fn resolve_format(path: &std::path::Path, declared: &str) -> String {
    if !declared.is_empty() {
        return declared.to_ascii_lowercase();
    }
    match path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .as_deref()
    {
        Some("fasta") | Some("fa") | Some("fna") | Some("fasta.txt") => "fasta".into(),
        Some("xlsx") | Some("xls") => "xlsx".into(),
        Some("csv") => "csv".into(),
        Some("gb") | Some("gbk") | Some("genbank") => "genbank".into(),
        Some("ab1") => "ab1".into(),
        _ => "fasta".into(),
    }
}
