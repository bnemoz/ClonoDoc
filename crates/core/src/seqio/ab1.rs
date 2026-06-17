//! AB1 / ABIF parser (`docs/03_ARCHITECTURE.md` §5).
//!
//! ABIF = header + a directory of 28-byte tagged records (big-endian). We extract
//! base calls from `PBAS` tag 2 (fallback tag 1). Quality (`PCON`) is parsed only
//! as **advisory** and is off by default — Plasmidsaurus Nanopore-derived AB1
//! quality is not meaningful, so verdicts never gate on it.

use crate::{Error, Result};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct Ab1 {
    pub bases: String,
    /// Per-base quality, advisory only (default consumers ignore this).
    pub quality: Option<Vec<u8>>,
}

struct DirEntry {
    name: String,
    number: i32,
    element_type: i16,
    #[allow(dead_code)]
    element_size: i16,
    num_elements: i32,
    data_size: i32,
    data_offset: i32,
}

fn be_i16(b: &[u8], o: usize) -> i16 {
    i16::from_be_bytes([b[o], b[o + 1]])
}
fn be_i32(b: &[u8], o: usize) -> i32 {
    i32::from_be_bytes([b[o], b[o + 1], b[o + 2], b[o + 3]])
}

fn read_dir_entry(b: &[u8], off: usize) -> Option<DirEntry> {
    if off + 28 > b.len() {
        return None;
    }
    Some(DirEntry {
        name: String::from_utf8_lossy(&b[off..off + 4]).to_string(),
        number: be_i32(b, off + 4),
        element_type: be_i16(b, off + 8),
        element_size: be_i16(b, off + 10),
        num_elements: be_i32(b, off + 12),
        data_size: be_i32(b, off + 16),
        data_offset: be_i32(b, off + 20),
    })
}

impl DirEntry {
    /// Slice of the raw data this entry points at. Values ≤ 4 bytes are stored
    /// inline in the `data_offset` field itself.
    fn data<'a>(&self, b: &'a [u8], raw_offset_field: &'a [u8]) -> Option<&'a [u8]> {
        let size = self.data_size.max(0) as usize;
        if size <= 4 {
            Some(&raw_offset_field[..size.min(4)])
        } else {
            let off = self.data_offset as usize;
            if off + size <= b.len() {
                Some(&b[off..off + size])
            } else {
                None
            }
        }
    }
}

pub fn read_path(path: &Path) -> Result<Ab1> {
    let bytes = std::fs::read(path)?;
    parse(&bytes).map_err(|m| Error::parse(path.display().to_string(), m))
}

pub fn parse(b: &[u8]) -> std::result::Result<Ab1, String> {
    if b.len() < 34 || &b[0..4] != b"ABIF" {
        return Err("not an ABIF/AB1 file (bad magic)".into());
    }
    // The root directory entry sits at byte offset 6.
    let root = read_dir_entry(b, 6).ok_or("truncated ABIF header")?;
    let dir_start = root.data_offset as usize;
    let n = root.num_elements.max(0) as usize;

    let mut bases: Option<String> = None;
    let mut bases_rank = i32::MAX; // prefer tag number 2, then 1
    let mut quality: Option<Vec<u8>> = None;

    for i in 0..n {
        let off = dir_start + i * 28;
        let Some(entry) = read_dir_entry(b, off) else {
            break;
        };
        // The inline-data field is the 4 bytes at off+20 within the entry.
        let raw_offset_field = &b[off + 20..off + 24];
        match entry.name.as_str() {
            "PBAS" => {
                if let Some(data) = entry.data(b, raw_offset_field) {
                    let s: String = data
                        .iter()
                        .map(|&c| c as char)
                        .filter(|c| c.is_ascii_alphabetic())
                        .collect();
                    // Prefer PBAS2 over PBAS1.
                    let rank = if entry.number == 2 { 0 } else { 1 };
                    if rank < bases_rank {
                        bases = Some(s);
                        bases_rank = rank;
                    }
                }
            }
            "PCON" if (entry.element_type == 1 || entry.element_type == 2) => {
                if let Some(data) = entry.data(b, raw_offset_field) {
                    quality = Some(data.to_vec());
                }
            }
            _ => {}
        }
    }

    match bases {
        Some(bases) => Ok(Ab1 { bases, quality }),
        None => Err("no PBAS base-call tag found".into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a minimal synthetic ABIF with one PBAS2 record to exercise the parser.
    fn synth_ab1(bases: &str) -> Vec<u8> {
        let mut out = Vec::new();
        out.extend_from_slice(b"ABIF");
        out.extend_from_slice(&1i16.to_be_bytes()); // version
                                                    // Root dir entry (28 bytes) at offset 6. It describes the directory.
                                                    // We'll place the directory array right after a small gap.
        let dir_offset: i32 = 6 + 28; // directory starts immediately after root entry
        let num_entries: i32 = 1;
        // root entry
        out.extend_from_slice(b"tdir");
        out.extend_from_slice(&1i32.to_be_bytes()); // number
        out.extend_from_slice(&1023i16.to_be_bytes()); // element type
        out.extend_from_slice(&28i16.to_be_bytes()); // element size
        out.extend_from_slice(&num_entries.to_be_bytes());
        out.extend_from_slice(&(num_entries * 28).to_be_bytes()); // data size
        out.extend_from_slice(&dir_offset.to_be_bytes());
        out.extend_from_slice(&0i32.to_be_bytes()); // handle

        // The PBAS data will live after the directory array.
        let data_offset: i32 = dir_offset + 28;
        // PBAS2 entry
        out.extend_from_slice(b"PBAS");
        out.extend_from_slice(&2i32.to_be_bytes()); // number = 2
        out.extend_from_slice(&2i16.to_be_bytes()); // element type char
        out.extend_from_slice(&1i16.to_be_bytes()); // element size
        out.extend_from_slice(&(bases.len() as i32).to_be_bytes());
        out.extend_from_slice(&(bases.len() as i32).to_be_bytes()); // data size
        out.extend_from_slice(&data_offset.to_be_bytes());
        out.extend_from_slice(&0i32.to_be_bytes());

        // base-call payload
        out.extend_from_slice(bases.as_bytes());
        out
    }

    #[test]
    fn parses_synthetic_pbas2() {
        let ab1 = synth_ab1("ACGTACGTNN");
        let parsed = parse(&ab1).unwrap();
        assert_eq!(parsed.bases, "ACGTACGTNN");
    }

    #[test]
    fn rejects_non_abif() {
        assert!(parse(b"NOPEnotanabiffile___________________").is_err());
    }
}
