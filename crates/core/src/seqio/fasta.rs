//! FASTA reader. Tolerates wrapped lines, blank lines, lowercase, gaps and `*`.

use super::SeqRecord;
use crate::seq;
use crate::Result;
use std::path::Path;

pub fn read_path(path: &Path) -> Result<Vec<SeqRecord>> {
    let text = std::fs::read_to_string(path)?;
    Ok(parse(&text))
}

/// Parse FASTA text into records. The sequence is cleaned (uppercased, gaps and
/// non-alphabetic characters dropped).
pub fn parse(text: &str) -> Vec<SeqRecord> {
    let mut records = Vec::new();
    let mut cur_id: Option<String> = None;
    let mut cur_seq = String::new();
    for line in text.lines() {
        let line = line.trim_end();
        if let Some(rest) = line.strip_prefix('>') {
            if let Some(id) = cur_id.take() {
                records.push(SeqRecord {
                    id,
                    sequence: seq::clean(&cur_seq),
                });
                cur_seq.clear();
            }
            // The id is the first whitespace-delimited token after '>'.
            let id = rest.split_whitespace().next().unwrap_or("").to_string();
            cur_id = Some(id);
        } else if cur_id.is_some() {
            cur_seq.push_str(line.trim());
        }
    }
    if let Some(id) = cur_id.take() {
        records.push(SeqRecord {
            id,
            sequence: seq::clean(&cur_seq),
        });
    }
    records
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_wrapped_and_lowercase() {
        let text = ">rec1 some description\nACGT\nacgt\n\n>rec2\nTTTT";
        let recs = parse(text);
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].id, "rec1");
        assert_eq!(recs[0].sequence, "ACGTACGT");
        assert_eq!(recs[1].id, "rec2");
    }
}
