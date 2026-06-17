//! Pairwise alignment (`docs/01_DESIGN.md` §1.4).
//!
//! A self-contained affine-gap (Gotoh) global aligner plus BLOSUM62 protein
//! scoring and configurable nucleotide scoring. All scores/penalties come from
//! the library `alignment` block — no magic numbers in the engine.
//!
//! We deliberately implement the aligner here rather than pulling a heavyweight
//! bioinformatics crate: the sequences are short (≤ ~6 kb), the algorithm is
//! textbook, and a local implementation is deterministic, dependency-light, and
//! exactly testable against the golden fixtures.

use crate::model::AlignmentSettings;

/// Outcome of a global pairwise alignment.
#[derive(Debug, Clone)]
pub struct Alignment {
    pub score: i32,
    pub aligned_a: String,
    pub aligned_b: String,
    /// Number of aligned columns where both residues are present and equal.
    pub matches: usize,
    /// Aligned columns where both present but differ.
    pub mismatches: usize,
    /// Columns where exactly one side is a gap.
    pub gaps: usize,
}

impl Alignment {
    /// Percent identity over aligned (non-gap) columns.
    pub fn identity(&self) -> f64 {
        let aligned = self.matches + self.mismatches;
        if aligned == 0 {
            0.0
        } else {
            self.matches as f64 / aligned as f64
        }
    }
    /// Percent identity over the full alignment length (gaps count against).
    pub fn identity_with_gaps(&self) -> f64 {
        let total = self.matches + self.mismatches + self.gaps;
        if total == 0 {
            0.0
        } else {
            self.matches as f64 / total as f64
        }
    }
}

/// A scoring scheme: substitution score for two residues + affine gap penalties.
pub struct Scoring<'a> {
    pub sub: Box<dyn Fn(u8, u8) -> i32 + 'a>,
    pub gap_open: i32,
    pub gap_extend: i32,
}

impl<'a> Scoring<'a> {
    /// Protein BLOSUM62 scoring from the library settings.
    pub fn protein(s: &AlignmentSettings) -> Scoring<'static> {
        Scoring {
            sub: Box::new(blosum62),
            gap_open: s.protein_gap_open,
            gap_extend: s.protein_gap_extend,
        }
    }
    /// Nucleotide match/mismatch scoring from the library settings.
    pub fn nucleotide(s: &AlignmentSettings) -> Scoring<'static> {
        let m = s.nt_match;
        let mm = s.nt_mismatch;
        Scoring {
            sub: Box::new(move |a, b| if a == b { m } else { mm }),
            gap_open: s.nt_gap_open,
            gap_extend: s.nt_gap_extend,
        }
    }
}

const NEG: i32 = i32::MIN / 4; // sentinel "very negative" that won't overflow on add

/// Global (Needleman–Wunsch) alignment with affine gaps (Gotoh, three matrices).
///
/// `gap_open` is charged on the *first* gap residue, `gap_extend` on each
/// subsequent one. Both are expected to be negative.
pub fn global_affine(a: &str, b: &str, sc: &Scoring) -> Alignment {
    let a = a.as_bytes();
    let b = b.as_bytes();
    let n = a.len();
    let m = b.len();

    // M: best score ending with a (mis)match; X: ending with gap in b (a aligned to '-');
    // Y: ending with gap in a ('-' aligned to b).
    let idx = |i: usize, j: usize| i * (m + 1) + j;
    let mut mm = vec![NEG; (n + 1) * (m + 1)];
    let mut gx = vec![NEG; (n + 1) * (m + 1)];
    let mut gy = vec![NEG; (n + 1) * (m + 1)];
    // Traceback: 0=diag(M), 1=up(X, gap in b), 2=left(Y, gap in a).
    let mut tb = vec![0u8; (n + 1) * (m + 1)];

    mm[idx(0, 0)] = 0;
    for i in 1..=n {
        gx[idx(i, 0)] = sc.gap_open + (i as i32 - 1) * sc.gap_extend;
        mm[idx(i, 0)] = gx[idx(i, 0)];
        tb[idx(i, 0)] = 1;
    }
    for j in 1..=m {
        gy[idx(0, j)] = sc.gap_open + (j as i32 - 1) * sc.gap_extend;
        mm[idx(0, j)] = gy[idx(0, j)];
        tb[idx(0, j)] = 2;
    }

    for i in 1..=n {
        for j in 1..=m {
            let s = (sc.sub)(a[i - 1], b[j - 1]);
            let diag = best3(
                mm[idx(i - 1, j - 1)],
                gx[idx(i - 1, j - 1)],
                gy[idx(i - 1, j - 1)],
            );
            let m_here = diag.saturating_add(s);

            // Gap in b (consume a): either open from M/Y or extend X.
            let open_x =
                best3(mm[idx(i - 1, j)], gy[idx(i - 1, j)], NEG).saturating_add(sc.gap_open);
            let ext_x = gx[idx(i - 1, j)].saturating_add(sc.gap_extend);
            gx[idx(i, j)] = open_x.max(ext_x);

            // Gap in a (consume b).
            let open_y =
                best3(mm[idx(i, j - 1)], gx[idx(i, j - 1)], NEG).saturating_add(sc.gap_open);
            let ext_y = gy[idx(i, j - 1)].saturating_add(sc.gap_extend);
            gy[idx(i, j)] = open_y.max(ext_y);

            mm[idx(i, j)] = m_here;
            // Decide which state wins this cell for traceback.
            let best = m_here.max(gx[idx(i, j)]).max(gy[idx(i, j)]);
            tb[idx(i, j)] = if best == m_here {
                0
            } else if best == gx[idx(i, j)] {
                1
            } else {
                2
            };
        }
    }

    // Traceback from (n, m), following the winning state per cell.
    let mut i = n;
    let mut j = m;
    let mut ra = Vec::new();
    let mut rb = Vec::new();
    while i > 0 || j > 0 {
        let dir = if i == 0 {
            2
        } else if j == 0 {
            1
        } else {
            tb[idx(i, j)]
        };
        match dir {
            0 => {
                ra.push(a[i - 1]);
                rb.push(b[j - 1]);
                i -= 1;
                j -= 1;
            }
            1 => {
                ra.push(a[i - 1]);
                rb.push(b'-');
                i -= 1;
            }
            _ => {
                ra.push(b'-');
                rb.push(b[j - 1]);
                j -= 1;
            }
        }
    }
    ra.reverse();
    rb.reverse();

    let mut matches = 0;
    let mut mismatches = 0;
    let mut gaps = 0;
    for (x, y) in ra.iter().zip(rb.iter()) {
        if *x == b'-' || *y == b'-' {
            gaps += 1;
        } else if x == y {
            matches += 1;
        } else {
            mismatches += 1;
        }
    }

    let score = best3(mm[idx(n, m)], gx[idx(n, m)], gy[idx(n, m)]);
    Alignment {
        score,
        aligned_a: String::from_utf8(ra).unwrap(),
        aligned_b: String::from_utf8(rb).unwrap(),
        matches,
        mismatches,
        gaps,
    }
}

#[inline]
fn best3(a: i32, b: i32, c: i32) -> i32 {
    a.max(b).max(c)
}

/// Convenience: protein global alignment using library settings.
pub fn align_protein(a: &str, b: &str, settings: &AlignmentSettings) -> Alignment {
    global_affine(a, b, &Scoring::protein(settings))
}

/// Convenience: nucleotide global alignment using library settings.
pub fn align_nt(a: &str, b: &str, settings: &AlignmentSettings) -> Alignment {
    global_affine(a, b, &Scoring::nucleotide(settings))
}

/// Ungapped percent identity of `query` against the best offset within `subject`
/// (used for fast backbone identity over a long landmark). Returns the best
/// fraction of matching positions over `query.len()`.
pub fn best_ungapped_identity(query: &str, subject: &str) -> f64 {
    let q = query.as_bytes();
    let s = subject.as_bytes();
    if q.is_empty() || s.len() < q.len() {
        return 0.0;
    }
    let mut best = 0usize;
    for off in 0..=(s.len() - q.len()) {
        let mut hit = 0;
        for k in 0..q.len() {
            if q[k] == s[off + k] {
                hit += 1;
            }
        }
        if hit > best {
            best = hit;
        }
    }
    best as f64 / q.len() as f64
}

// ---------------------------------------------------------------------------
// BLOSUM62
// ---------------------------------------------------------------------------

/// BLOSUM62 substitution score for two amino-acid bytes (uppercase). Unknown or
/// `X` residues score 0 vs anything (neutral), gaps should never reach here.
pub fn blosum62(a: u8, b: u8) -> i32 {
    let ia = aa_index(a);
    let ib = aa_index(b);
    match (ia, ib) {
        (Some(x), Some(y)) => BLOSUM62_MATRIX[x][y] as i32,
        _ => 0,
    }
}

const AA_ORDER: &[u8; 24] = b"ARNDCQEGHILKMFPSTWYVBZX*";

fn aa_index(c: u8) -> Option<usize> {
    AA_ORDER.iter().position(|&x| x == c.to_ascii_uppercase())
}

// Standard BLOSUM62, ordered as AA_ORDER (A R N D C Q E G H I L K M F P S T W Y V B Z X *).
#[rustfmt::skip]
const BLOSUM62_MATRIX: [[i8; 24]; 24] = [
    [ 4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,-2,-1, 0,-4],
    [-1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,-1, 0,-1,-4],
    [-2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3, 3, 0,-1,-4],
    [-2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [ 0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,-3,-3,-2,-4],
    [-1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2, 0, 3,-1,-4],
    [-1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,-1,-2,-1,-4],
    [-2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3, 0, 0,-1,-4],
    [-1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,-3,-3,-1,-4],
    [-1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,-4,-3,-1,-4],
    [-1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2, 0, 1,-1,-4],
    [-1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,-3,-1,-1,-4],
    [-2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,-3,-3,-1,-4],
    [-1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,-2,-1,-2,-4],
    [ 1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2, 0, 0, 0,-4],
    [ 0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,-1,-1, 0,-4],
    [-3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,-4,-3,-2,-4],
    [-2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,-3,-2,-1,-4],
    [ 0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4,-3,-2,-1,-4],
    [-2,-1, 3, 4,-3, 0, 1,-1, 0,-3,-4, 0,-3,-3,-2, 0,-1,-4,-3,-3, 4, 1,-1,-4],
    [-1, 0, 0, 1,-3, 3, 4,-2, 0,-3,-3, 1,-1,-3,-1, 0,-1,-3,-2,-2, 1, 4,-1,-4],
    [ 0,-1,-1,-1,-2,-1,-1,-1,-1,-1,-1,-1,-1,-2,-2, 0, 0,-2,-1,-1,-1,-1,-1,-4],
    [-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4,-4, 1],
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn blosum_self_scores() {
        assert_eq!(blosum62(b'A', b'A'), 4);
        assert_eq!(blosum62(b'W', b'W'), 11);
        assert_eq!(blosum62(b'C', b'C'), 9);
        // Conservative substitution scores positive.
        assert!(blosum62(b'I', b'L') > 0);
    }

    #[test]
    fn identical_protein_aligns_perfectly() {
        let s = AlignmentSettings::default();
        let al = align_protein("MELGLRWVFL", "MELGLRWVFL", &s);
        assert_eq!(al.mismatches, 0);
        assert_eq!(al.gaps, 0);
        assert!((al.identity() - 1.0).abs() < 1e-9);
    }

    #[test]
    fn single_substitution_detected() {
        let s = AlignmentSettings::default();
        let al = align_protein("MELGLRWVFL", "MELGKRWVFL", &s);
        assert_eq!(al.mismatches, 1);
        assert_eq!(al.gaps, 0);
    }

    #[test]
    fn deletion_creates_a_gap() {
        let s = AlignmentSettings::default();
        let al = align_protein("MELGLRWVFL", "MELGRWVFL", &s);
        assert_eq!(al.gaps, 1);
    }

    #[test]
    fn nt_identity() {
        let s = AlignmentSettings::default();
        let al = align_nt("ACGTACGT", "ACGTACGT", &s);
        assert!((al.identity() - 1.0).abs() < 1e-9);
    }
}
