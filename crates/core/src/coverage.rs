//! Per-base read-depth coverage of an expected construct.
//!
//! Each read is anchored to the expected reference by a seed match (tolerating
//! the arbitrary rotation/strand of a circular Plasmidsaurus read), and the span
//! it covers increments the depth. The result drives the Wetlab coverage view.
//!
//! This is deliberately a coarse, fast mapping (seed + span), not a full
//! alignment: it answers "which positions of the construct are supported by
//! reads, and how deeply", which is what the coverage chart shows.

use crate::seq;

/// A coverage profile across a reference of length `len`.
#[derive(Debug, Clone)]
pub struct Coverage {
    /// Per-base depth, length == reference length.
    pub depth: Vec<u32>,
}

impl Coverage {
    pub fn len(&self) -> usize {
        self.depth.len()
    }
    pub fn is_empty(&self) -> bool {
        self.depth.is_empty()
    }
    pub fn mean(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        self.depth.iter().map(|&d| d as f64).sum::<f64>() / self.depth.len() as f64
    }
    pub fn max(&self) -> u32 {
        self.depth.iter().copied().max().unwrap_or(0)
    }
    /// Fraction of positions with depth ≥ 1.
    pub fn breadth(&self) -> f64 {
        if self.depth.is_empty() {
            return 0.0;
        }
        let covered = self.depth.iter().filter(|&&d| d > 0).count();
        covered as f64 / self.depth.len() as f64
    }
}

const SEED: usize = 18;

/// Build a coverage profile for `reference` from a set of reads.
///
/// For each read, the best seed offset is found on the reference (and its
/// reverse complement / circular wrap); the covered span then increments depth.
/// Reads that don't anchor confidently are skipped (they contribute nothing
/// rather than smearing false coverage).
pub fn coverage_profile(reference: &str, reads: &[String]) -> Coverage {
    let reference = seq::clean(reference);
    let len = reference.len();
    let mut depth = vec![0u32; len];
    if len == 0 {
        return Coverage { depth };
    }
    // Double the reference to emulate circularity for anchoring.
    let doubled = format!("{reference}{reference}");

    for raw in reads {
        let read = seq::clean(raw);
        if read.len() < SEED {
            continue;
        }
        // Try both strands; keep the better-anchored orientation.
        let candidates = [read.clone(), seq::revcomp(&read)];
        let mut best: Option<(usize, usize)> = None; // (start_in_ref, read_len)
        for cand in &candidates {
            if let Some(start) = anchor(&doubled, cand, len) {
                best = Some((start, cand.len()));
                break;
            }
        }
        if let Some((start, rlen)) = best {
            // Mark the covered span, wrapping around the circular reference.
            let span = rlen.min(len);
            for k in 0..span {
                let idx = (start + k) % len;
                depth[idx] = depth[idx].saturating_add(1);
            }
        }
    }
    Coverage { depth }
}

/// Find where a read anchors on the (doubled) reference using its 5′ seed.
/// Returns a 0-based start within one reference length, or None.
fn anchor(doubled: &str, read: &str, ref_len: usize) -> Option<usize> {
    let seed = &read.as_bytes()[..SEED.min(read.len())];
    let hay = doubled.as_bytes();
    // Scan one full turn of the reference for an exact seed match.
    for start in 0..ref_len {
        if start + seed.len() <= hay.len() && &hay[start..start + seed.len()] == seed {
            return Some(start);
        }
    }
    None
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_read_covers_everything_once() {
        let reference = "ATGCGTACGTTAGCCGATCGATCGGATCAGCTAGCTAGCAT";
        let cov = coverage_profile(reference, &[reference.to_string()]);
        assert_eq!(cov.len(), reference.len());
        assert!(cov.depth.iter().all(|&d| d == 1));
        assert!((cov.breadth() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn two_reads_stack_depth() {
        let reference = "ATGCGTACGTTAGCCGATCGATCGGATCAGCTAGCTAGCAT";
        let cov = coverage_profile(reference, &[reference.to_string(), reference.to_string()]);
        assert!(cov.depth.iter().all(|&d| d == 2));
        assert!((cov.mean() - 2.0).abs() < 1e-9);
    }

    #[test]
    fn partial_read_leaves_a_gap() {
        let reference = "ATGCGTACGTTAGCCGATCGATCGGATCAGCTAGCTAGCATTTGGCCAA";
        // A read covering only the first ~half.
        let half = &reference[..24];
        let cov = coverage_profile(reference, &[half.to_string()]);
        assert!(cov.depth[0] >= 1);
        assert_eq!(cov.depth[reference.len() - 1], 0, "tail uncovered");
        assert!(cov.breadth() < 1.0);
    }

    #[test]
    fn reverse_complement_read_still_maps() {
        let reference = "ATGCGTACGTTAGCCGATCGATCGGATCAGCTAGCTAGCAT";
        let rc = seq::revcomp(reference);
        let cov = coverage_profile(reference, &[rc]);
        assert!(cov.breadth() > 0.9, "RC read should cover the reference");
    }
}
