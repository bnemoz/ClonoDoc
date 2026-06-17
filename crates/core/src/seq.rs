//! Sequence primitives: nt/AA classification and translation.
//!
//! The foundational principle (`docs/01_DESIGN.md` §0): codon optimization means
//! the variable region is verified at the **amino-acid** level, while backbone and
//! overhangs are verified at the **nucleotide** level. Everything here is the
//! shared machinery both gates rely on.

/// The molecule kind a string represents.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum SeqType {
    Nt,
    Aa,
}

/// Normalize a raw sequence string: uppercase, drop whitespace, `*`, `-` and digits.
///
/// Liberal in input (the design mantra): tolerate wrapped lines, lowercase, gaps,
/// and FASTA stop markers.
pub fn clean(raw: &str) -> String {
    raw.chars()
        .filter(|c| c.is_ascii_alphabetic())
        .map(|c| c.to_ascii_uppercase())
        .collect()
}

/// Autodetect whether a cleaned sequence is nucleotide or amino acid.
///
/// Per `docs/01_DESIGN.md` §1.2: treat as **nt** iff the alphabet ⊆ {A,C,G,T,U,N}
/// and (length % 3 == 0 OR ≥90% of characters are in the DNA alphabet). Otherwise **AA**.
pub fn detect_type(seq: &str) -> SeqType {
    let s = clean(seq);
    if s.is_empty() {
        return SeqType::Aa;
    }
    let dna = |c: char| matches!(c, 'A' | 'C' | 'G' | 'T' | 'U' | 'N');
    let dna_count = s.chars().filter(|&c| dna(c)).count();
    let frac = dna_count as f64 / s.len() as f64;
    let all_dna = dna_count == s.len();
    // nt iff: a clean codon-length DNA string, OR overwhelmingly DNA alphabet
    // (the ≥90% rule covers partial/ambiguous-base reads not a multiple of 3).
    if (all_dna && s.len().is_multiple_of(3)) || frac >= 0.90 {
        SeqType::Nt
    } else {
        SeqType::Aa
    }
}

/// Reverse complement of a nucleotide string (cleaned first).
pub fn revcomp(seq: &str) -> String {
    clean(seq)
        .chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'U' => 'A',
            'G' => 'C',
            'C' => 'G',
            'N' => 'N',
            other => other,
        })
        .collect()
}

/// Translate a nucleotide sequence in frame 0 using the standard genetic code.
///
/// `*` marks a stop. Incomplete trailing codons are dropped. Unknown codons
/// (containing `N` or non-ACGT) translate to `X`.
pub fn translate(nt: &str) -> String {
    let s = clean(nt);
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() / 3);
    let mut i = 0;
    while i + 3 <= bytes.len() {
        out.push(codon(&bytes[i..i + 3]));
        i += 3;
    }
    out
}

/// Translate until (and excluding) the first stop codon.
pub fn translate_to_stop(nt: &str) -> String {
    let full = translate(nt);
    match full.find('*') {
        Some(idx) => full[..idx].to_string(),
        None => full,
    }
}

/// Map a single codon (3 bytes) to its amino acid.
fn codon(c: &[u8]) -> char {
    let key = [
        c[0].to_ascii_uppercase(),
        c[1].to_ascii_uppercase(),
        c[2].to_ascii_uppercase(),
    ];
    // U -> T normalization for RNA input.
    let n = |b: u8| if b == b'U' { b'T' } else { b };
    match [n(key[0]), n(key[1]), n(key[2])] {
        [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => 'F',
        [b'T', b'T', b'A'] | [b'T', b'T', b'G'] => 'L',
        [b'C', b'T', _] => 'L',
        [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => 'I',
        [b'A', b'T', b'G'] => 'M',
        [b'G', b'T', _] => 'V',
        [b'T', b'C', _] => 'S',
        [b'C', b'C', _] => 'P',
        [b'A', b'C', _] => 'T',
        [b'G', b'C', _] => 'A',
        [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => 'Y',
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] => '*',
        [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => 'H',
        [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => 'Q',
        [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => 'N',
        [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => 'K',
        [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => 'D',
        [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => 'E',
        [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => 'C',
        [b'T', b'G', b'A'] => '*',
        [b'T', b'G', b'G'] => 'W',
        [b'C', b'G', _] => 'R',
        [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => 'S',
        [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => 'R',
        [b'G', b'G', _] => 'G',
        _ => 'X',
    }
}

/// GC fraction over a cleaned nucleotide string. Returns 0.0 for empty input.
pub fn gc_fraction(nt: &str) -> f64 {
    let s = clean(nt);
    if s.is_empty() {
        return 0.0;
    }
    let gc = s.chars().filter(|&c| c == 'G' || c == 'C').count();
    gc as f64 / s.len() as f64
}

/// Maximum GC fraction over any window of the given size. Returns the global
/// fraction if the sequence is shorter than the window.
pub fn max_windowed_gc(nt: &str, window: usize) -> f64 {
    let s = clean(nt);
    if s.len() <= window || window == 0 {
        return gc_fraction(&s);
    }
    let bytes = s.as_bytes();
    let is_gc = |b: u8| b == b'G' || b == b'C';
    let mut count = bytes[..window].iter().filter(|&&b| is_gc(b)).count();
    let mut best = count;
    for i in window..bytes.len() {
        if is_gc(bytes[i]) {
            count += 1;
        }
        if is_gc(bytes[i - window]) {
            count -= 1;
        }
        best = best.max(count);
    }
    best as f64 / window as f64
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn detects_nt_and_aa() {
        assert_eq!(detect_type("ATGCATGCA"), SeqType::Nt);
        assert_eq!(detect_type("atg cat gca"), SeqType::Nt);
        // A protein with characters outside the DNA alphabet.
        assert_eq!(detect_type("MELGLRWVFLVAILEGVQC"), SeqType::Aa);
        // Ambiguous-base read, length not %3 but overwhelmingly DNA.
        assert_eq!(detect_type("ATGCATGCATG"), SeqType::Nt);
    }

    #[test]
    fn translates_standard_code() {
        // The French vector leader.
        assert_eq!(
            translate("ATGGAACTGGGGCTCCGCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGT"),
            "MELGLRWVFLVAILEGVQC"
        );
    }

    #[test]
    fn stop_truncation_works() {
        assert_eq!(translate_to_stop("ATGTAAATG"), "M");
        assert_eq!(translate("ATGTAAATG"), "M*M");
    }

    #[test]
    fn revcomp_roundtrip() {
        assert_eq!(revcomp("ATGC"), "GCAT");
        assert_eq!(revcomp(&revcomp("ATGCGGTA")), "ATGCGGTA");
    }

    #[test]
    fn gc_math() {
        assert!((gc_fraction("GGGGCCCC") - 1.0).abs() < 1e-9);
        assert!((gc_fraction("ATATATAT") - 0.0).abs() < 1e-9);
        assert!((gc_fraction("ATGC") - 0.5).abs() < 1e-9);
    }
}
