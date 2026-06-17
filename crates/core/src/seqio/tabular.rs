//! CSV / XLSX readers with liberal (case-insensitive, trimmed) header matching.

use super::{GroundTruthRow, SeqRecord};
use crate::seq;
use crate::{Error, Result};
use std::collections::BTreeMap;
use std::path::Path;

/// Normalize a header cell for matching: lowercase, trimmed, non-alphanumerics dropped.
fn norm_header(h: &str) -> String {
    h.trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|c| c.is_ascii_alphanumeric())
        .collect()
}

fn find_col(headers: &[String], aliases: &[&str]) -> Option<usize> {
    let normed: Vec<String> = headers.iter().map(|h| norm_header(h)).collect();
    for alias in aliases {
        let a = norm_header(alias);
        if let Some(i) = normed.iter().position(|h| *h == a) {
            return Some(i);
        }
    }
    None
}

const ID_ALIASES: &[&str] = &["sequence_id", "seq_id", "id", "name", "sequenceid", "seqid"];
const SEQ_ALIASES: &[&str] = &["sequence", "seq"];

// ---------------------------------------------------------------------------
// CSV
// ---------------------------------------------------------------------------

fn csv_rows(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    let mut rdr = csv::ReaderBuilder::new()
        .flexible(true)
        .from_path(path)
        .map_err(|e| Error::parse(path.display().to_string(), e.to_string()))?;
    let headers: Vec<String> = rdr
        .headers()
        .map_err(|e| Error::parse(path.display().to_string(), e.to_string()))?
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mut rows = Vec::new();
    for rec in rdr.records() {
        let rec = rec.map_err(|e| Error::parse(path.display().to_string(), e.to_string()))?;
        rows.push(rec.iter().map(|s| s.to_string()).collect());
    }
    Ok((headers, rows))
}

pub fn read_csv_id_seq(path: &Path) -> Result<Vec<SeqRecord>> {
    let (headers, rows) = csv_rows(path)?;
    let id_col = find_col(&headers, ID_ALIASES)
        .ok_or_else(|| Error::parse(path.display().to_string(), "no id column found"))?;
    let seq_col = find_col(&headers, SEQ_ALIASES)
        .ok_or_else(|| Error::parse(path.display().to_string(), "no sequence column found"))?;
    Ok(rows
        .into_iter()
        .filter_map(|r| {
            let id = r.get(id_col)?.trim().to_string();
            let s = seq::clean(r.get(seq_col)?);
            if id.is_empty() {
                None
            } else {
                Some(SeqRecord { id, sequence: s })
            }
        })
        .collect())
}

// ---------------------------------------------------------------------------
// XLSX (calamine)
// ---------------------------------------------------------------------------

fn xlsx_rows(path: &Path) -> Result<(Vec<String>, Vec<Vec<String>>)> {
    use calamine::{open_workbook_auto, Data, Reader};
    let mut wb = open_workbook_auto(path)
        .map_err(|e| Error::parse(path.display().to_string(), e.to_string()))?;
    let sheet_name = wb
        .sheet_names()
        .first()
        .cloned()
        .ok_or_else(|| Error::parse(path.display().to_string(), "workbook has no sheets"))?;
    let range = wb
        .worksheet_range(&sheet_name)
        .map_err(|e| Error::parse(path.display().to_string(), e.to_string()))?;
    let cell = |d: &Data| -> String {
        match d {
            Data::String(s) => s.clone(),
            Data::Float(f) => {
                if f.fract() == 0.0 {
                    format!("{}", *f as i64)
                } else {
                    format!("{f}")
                }
            }
            Data::Int(i) => format!("{i}"),
            Data::Bool(b) => format!("{b}"),
            _ => String::new(),
        }
    };
    let mut iter = range.rows();
    let headers: Vec<String> = match iter.next() {
        Some(h) => h.iter().map(cell).collect(),
        None => return Ok((Vec::new(), Vec::new())),
    };
    let rows: Vec<Vec<String>> = iter.map(|r| r.iter().map(cell).collect()).collect();
    Ok((headers, rows))
}

pub fn read_xlsx_id_seq(path: &Path) -> Result<Vec<SeqRecord>> {
    let (headers, rows) = xlsx_rows(path)?;
    let id_col = find_col(&headers, ID_ALIASES)
        .ok_or_else(|| Error::parse(path.display().to_string(), "no id column found"))?;
    let seq_col = find_col(&headers, SEQ_ALIASES)
        .ok_or_else(|| Error::parse(path.display().to_string(), "no sequence column found"))?;
    Ok(rows
        .into_iter()
        .filter_map(|r| {
            let id = r.get(id_col)?.trim().to_string();
            let s = seq::clean(r.get(seq_col).map(|x| x.as_str()).unwrap_or(""));
            if id.is_empty() {
                None
            } else {
                Some(SeqRecord { id, sequence: s })
            }
        })
        .collect())
}

// ---------------------------------------------------------------------------
// Ground-truth panel (ab_id + H_seq + L_seq), CSV / XLSX / FASTA
// ---------------------------------------------------------------------------

const AB_ALIASES: &[&str] = &["ab_id", "abid", "antibody", "id", "name"];
const HEAVY_ALIASES: &[&str] = &["h_seq", "hseq", "heavy", "heavyseq", "hc", "vh", "h"];
const LIGHT_ALIASES: &[&str] = &["l_seq", "lseq", "light", "lightseq", "lc", "vl", "l"];

/// Read a ground-truth panel from a table (CSV/XLSX), honoring optional column
/// aliases supplied by the project (`columns: { ab_id, heavy, light }`).
pub fn read_ground_truth_table(
    path: &Path,
    overrides: &BTreeMap<String, String>,
) -> Result<Vec<GroundTruthRow>> {
    let fmt = super::resolve_format(path, "");
    let (headers, rows) = if fmt == "csv" {
        csv_rows(path)?
    } else {
        xlsx_rows(path)?
    };

    let col = |key: &str, defaults: &[&str]| -> Option<usize> {
        if let Some(alias) = overrides.get(key) {
            if let Some(i) = find_col(&headers, &[alias.as_str()]) {
                return Some(i);
            }
        }
        find_col(&headers, defaults)
    };

    let ab_col = col("ab_id", AB_ALIASES)
        .ok_or_else(|| Error::parse(path.display().to_string(), "no ab_id column found"))?;
    let h_col = col("heavy", HEAVY_ALIASES);
    let l_col = col("light", LIGHT_ALIASES);

    let mut out = Vec::new();
    for r in rows {
        let ab_id = match r.get(ab_col) {
            Some(s) if !s.trim().is_empty() => s.trim().to_string(),
            _ => continue,
        };
        let pick = |c: Option<usize>| -> Option<String> {
            c.and_then(|i| r.get(i))
                .map(|s| seq::clean(s))
                .filter(|s| !s.is_empty())
        };
        out.push(GroundTruthRow {
            ab_id,
            heavy: pick(h_col),
            light: pick(l_col),
        });
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn header_normalization() {
        assert_eq!(norm_header(" Sequence_ID "), "sequenceid");
        assert_eq!(norm_header("H_seq"), "hseq");
    }

    #[test]
    fn finds_columns_liberally() {
        let h = vec!["Sequence_ID".to_string(), "Sequence".to_string()];
        assert_eq!(find_col(&h, ID_ALIASES), Some(0));
        assert_eq!(find_col(&h, SEQ_ALIASES), Some(1));
    }
}
