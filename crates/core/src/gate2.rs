//! Gate 2 — sequencing verification, post-cloning (`docs/01_DESIGN.md` §3).
//!
//! Pipeline per read: topology-normalize (strand + rotation for full plasmids) →
//! backbone identity → locate & translate the ORF → global protein-align observed
//! vs expected → classify into the verdict taxonomy. Optional silent-SNP layer
//! when the IDT nt order is loaded, and sample-swap detection across the panel.

use crate::align::{self, best_ungapped_identity};
use crate::assemble;
use crate::model::{ChainClass, Library, Locus, OverhangSet, Project, Vector};
use crate::naming::{self, NameParse};
use crate::seq;
use crate::seqio::{GroundTruthRow, SeqRecord};
use crate::verdict::{Gate2Kind, Gate2Verdict, Mutation};
use std::collections::BTreeMap;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SeqMode {
    FullPlasmid,
    PartialSanger,
}

impl SeqMode {
    pub fn parse(s: &str) -> SeqMode {
        if s.eq_ignore_ascii_case("partial_sanger") {
            SeqMode::PartialSanger
        } else {
            SeqMode::FullPlasmid
        }
    }
}

pub struct Gate2Context<'a> {
    pub library: &'a Library,
    pub project: &'a Project,
    pub overhang_set: &'a OverhangSet,
    pub ground_truth: BTreeMap<String, GroundTruthRow>,
    /// (AB_ID_UPPER, chain_class) → expected insert core nt (from the IDT order).
    pub order_cores: BTreeMap<(String, ChainClass), String>,
    pub mode: SeqMode,
}

impl<'a> Gate2Context<'a> {
    pub fn new(
        library: &'a Library,
        project: &'a Project,
        overhang_set: &'a OverhangSet,
        ground_truth: &[GroundTruthRow],
        order_cores: BTreeMap<(String, ChainClass), String>,
        mode: SeqMode,
    ) -> Self {
        let map = ground_truth
            .iter()
            .map(|r| (r.ab_id.to_ascii_uppercase(), r.clone()))
            .collect();
        Gate2Context {
            library,
            project,
            overhang_set,
            ground_truth: map,
            order_cores,
            mode,
        }
    }

    fn resolve_parse(&self, record_id: &str, profile: &naming::NamingProfile) -> NameParse {
        if let Some(ov) = self.project.name_overrides.get(record_id) {
            return NameParse {
                ab_id: ov.ab_id.clone(),
                chain_class: ov.chain_class,
                confidence: naming::Confidence::High,
                needs_confirmation: false,
            };
        }
        naming::parse_name(record_id, profile)
    }

    /// Candidate vectors for backbone identification: the project's assigned
    /// vectors first, then any others in the library.
    fn candidate_vectors(&self) -> Vec<&'a Vector> {
        let mut out: Vec<&Vector> = Vec::new();
        for id in self.project.vector_assignments.values() {
            if let Some(v) = self.library.vector(id) {
                if !out.iter().any(|x| x.id == v.id) {
                    out.push(v);
                }
            }
        }
        for v in &self.library.vectors {
            if !out.iter().any(|x| x.id == v.id) {
                out.push(v);
            }
        }
        out
    }

    /// A codon-invariant backbone landmark for a vector: the start of the constant
    /// region (not codon-optimized), or the 3′ flank if no constant feature.
    fn backbone_landmark(v: &Vector) -> String {
        let s = seq::clean(&v.sequence);
        let start = v
            .features
            .get("constant_region")
            .map(|f| f.range0().start)
            .unwrap_or(v.insertion_site.oh3_start);
        let end = (start + 200).min(s.len());
        s[start..end].to_string()
    }

    pub fn run(&self, reads: &[SeqRecord]) -> Vec<Gate2Verdict> {
        let profile = self
            .library
            .naming_profile(&self.project.naming_profile)
            .cloned()
            .unwrap_or_else(naming::default_profile);
        reads.iter().map(|r| self.run_one(r, &profile)).collect()
    }

    fn run_one(&self, rec: &SeqRecord, profile: &naming::NamingProfile) -> Gate2Verdict {
        let parse = self.resolve_parse(&rec.id, profile);
        let mut v = Gate2Verdict {
            record_id: rec.id.clone(),
            ab_id: parse.ab_id.clone(),
            chain_class: parse.chain_class.as_str().to_string(),
            kind: Gate2Kind::Pass,
            reason: String::new(),
            backbone_vector: None,
            backbone_identity: None,
            backbone_observed: match self.mode {
                SeqMode::FullPlasmid => "full".into(),
                SeqMode::PartialSanger => "flanks_only".into(),
            },
            mutations: Vec::new(),
            premature_stop_aa: None,
            reads_through: None,
            suspected_identity: None,
        };

        // Resolve locus/vector. For light, detect κ/λ from the read.
        let locus = match parse.chain_class {
            ChainClass::Heavy => Some(Locus::Igh),
            ChainClass::Kappa => Some(Locus::Igk),
            ChainClass::Lambda => Some(Locus::Igl),
            ChainClass::Light => {
                assemble::detect_locus_from_protein(&best_frame_protein(&rec.sequence))
            }
            ChainClass::Unknown => None,
        };
        let expected_vector = locus
            .map(|l| l.chain_class())
            .and_then(|c| self.project.vector_for(c))
            .and_then(|id| self.library.vector(id));

        // Read adequacy.
        if seq::clean(&rec.sequence).len() < 60 {
            v.kind = Gate2Kind::InsufficientRead;
            v.reason = "read too short to cover the insert/junction".into();
            return v;
        }

        // Identify the backbone (full-plasmid only) and orient the read.
        let mut working = seq::clean(&rec.sequence);
        if self.mode == SeqMode::FullPlasmid {
            let candidates = self.candidate_vectors();
            let mut best: Option<(&Vector, f64)> = None;
            for cand in &candidates {
                let landmark = Self::backbone_landmark(cand);
                let fwd = best_ungapped_identity_circular(&working, &landmark);
                let rev = best_ungapped_identity_circular(&seq::revcomp(&working), &landmark);
                let id = fwd.max(rev);
                if best.map(|(_, b)| id > b).unwrap_or(true) {
                    best = Some((cand, id));
                }
            }
            if let Some((bv, id)) = best {
                v.backbone_vector = Some(bv.id.clone());
                v.backbone_identity = Some(id);
                // Normalize using the matched vector's landmark (strand + rotation),
                // anchoring on the leader so the ORF starts near index 0.
                let leader = {
                    let s = seq::clean(&bv.sequence);
                    s[..bv.insertion_site.oh5_end.min(s.len())].to_string()
                };
                let norm = assemble::normalize_circular(&working, &leader);
                working = norm.sequence;

                // WRONG_VECTOR: backbone class disagrees with the name's locus.
                if id >= self.library.alignment.backbone_identity_min {
                    let class_ok = match parse.chain_class {
                        ChainClass::Heavy => bv.chain_class == ChainClass::Heavy,
                        ChainClass::Kappa => bv.chain_class == ChainClass::Kappa,
                        ChainClass::Lambda => bv.chain_class == ChainClass::Lambda,
                        // A light token may legitimately hit κ or λ.
                        ChainClass::Light => {
                            matches!(bv.chain_class, ChainClass::Kappa | ChainClass::Lambda)
                        }
                        ChainClass::Unknown => true,
                    };
                    if !class_ok {
                        v.kind = Gate2Kind::WrongVector;
                        v.reason = format!(
                            "backbone matches {} ({:.1}% id) but the name says {}",
                            bv.display_name,
                            id * 100.0,
                            parse.chain_class.as_str()
                        );
                        return v;
                    }
                }
            }
        }

        let vector = match expected_vector.or_else(|| {
            v.backbone_vector
                .as_ref()
                .and_then(|id| self.library.vector(id))
        }) {
            Some(vec) => vec,
            None => {
                v.kind = Gate2Kind::InsufficientRead;
                v.reason = "no reference vector available to verify against".into();
                return v;
            }
        };
        v.reads_through = None;

        // Build the expected protein (leader + V + constant).
        let expected_core = self
            .order_cores
            .get(&(
                parse.ab_id.to_ascii_uppercase(),
                effective_class(parse.chain_class, locus),
            ))
            .cloned();
        let expected_aa = if let Some(core) = &expected_core {
            seq::translate_to_stop(&assemble::assemble(vector, core))
        } else if let Some(row) = self.ground_truth.get(&parse.ab_id.to_ascii_uppercase()) {
            let truth = match parse.chain_class {
                ChainClass::Heavy => row.heavy.as_deref(),
                _ => row.light.as_deref(),
            };
            match truth {
                Some(t) => {
                    let v_aa = if seq::detect_type(t) == seq::SeqType::Nt {
                        seq::translate_to_stop(t)
                    } else {
                        seq::clean(t)
                    };
                    let s = seq::clean(&vector.sequence);
                    let leader_aa =
                        seq::translate(&s[..vector.insertion_site.oh5_end.min(s.len())]);
                    let const_aa =
                        seq::translate_to_stop(&s[vector.insertion_site.oh3_start.min(s.len())..]);
                    format!("{leader_aa}{v_aa}{const_aa}")
                }
                None => {
                    v.kind = Gate2Kind::InsufficientRead;
                    v.reason =
                        "no expected sequence (order or ground truth) for this antibody".into();
                    return v;
                }
            }
        } else {
            v.kind = Gate2Kind::InsufficientRead;
            v.reason = "no expected sequence (order or ground truth) for this antibody".into();
            return v;
        };

        // Locate & translate the observed ORF. Anchor on the leader.
        let leader_aa = {
            let s = seq::clean(&vector.sequence);
            seq::translate(&s[..vector.insertion_site.oh5_end.min(s.len())])
        };
        let observed_aa = locate_orf_protein(&working, &leader_aa);

        let anchor = &vector.constant_anchor_aa;
        let observed_reads_through = observed_aa.contains(anchor.as_str());
        v.reads_through = Some(observed_reads_through);

        // Align observed vs expected.
        let al = align::align_protein(&observed_aa, &expected_aa, &self.library.alignment);

        if !observed_reads_through {
            // No read-through → junction frameshift (premature stop).
            let orf_full = locate_orf_full(&working, &leader_aa);
            let stop = orf_full.find('*');
            v.premature_stop_aa = stop;
            v.kind = Gate2Kind::JunctionFrameshift;
            v.reason = format!(
                "ORF diverges and hits a premature stop (aa {}); no read-through into the constant region ({} absent)",
                stop.map(|s| s.to_string()).unwrap_or_else(|| "?".into()),
                anchor
            );
            return v;
        }

        // Reads through. Empty-vector signature: the observed protein is the
        // leader fused straight to the constant region (no variable domain). We
        // detect it directly by comparing to the empty-vector ORF.
        let empty_aa = seq::translate_to_stop(&assemble::assemble(vector, ""));
        let v_len_aa = expected_aa.len().saturating_sub(empty_aa.len());
        if v_len_aa > 10 {
            let empty_al = align::align_protein(&observed_aa, &empty_aa, &self.library.alignment);
            // Near-identical to the empty ORF (few gaps) ⇒ no insert present.
            if empty_al.identity_with_gaps() > 0.97 && empty_al.gaps < observed_aa.len() / 10 + 2 {
                v.kind = Gate2Kind::EmptyVector;
                v.reason =
                    "leader reads straight into the constant region with no variable domain (empty/religated vector)"
                        .into();
                return v;
            }
        }

        // Collect substitutions in the V region.
        let muts = collect_mutations(&al.aligned_a, &al.aligned_b, leader_aa.len());
        if !muts.is_empty() {
            // Sample-swap check: does another panel member explain the read better?
            if let Some(swap) = self.detect_swap(&observed_aa, &parse, vector) {
                v.kind = Gate2Kind::WrongInsertSwap;
                v.suspected_identity = Some(swap);
                v.reason = format!(
                    "observed insert does not match {} but matches another panel member ({})",
                    parse.ab_id,
                    v.suspected_identity.as_deref().unwrap_or("?")
                );
                return v;
            }
            v.kind = Gate2Kind::PointMutation;
            let m = &muts[0];
            v.reason = format!(
                "{} point mutation(s) in the variable region; first at aa {} ({}→{})",
                muts.len(),
                m.position_aa,
                m.wt,
                m.mut_aa
            );
            v.mutations = muts;
            return v;
        }

        // Clean at the protein level. Optional silent-SNP layer (nt) if order loaded.
        if let Some(core) = &expected_core {
            if let Some(obs_core) = extract_insert_nt(&working, vector, core.len()) {
                let nt = align::align_nt(&obs_core, core, &self.library.alignment);
                if nt.mismatches > 0 || nt.gaps > 0 {
                    v.kind = Gate2Kind::SilentVariant;
                    v.reason = format!(
                        "protein identical to expected; {} synonymous nt change(s) in the insert (silent)",
                        nt.mismatches + nt.gaps
                    );
                    return v;
                }
            }
        }

        v.kind = Gate2Kind::Pass;
        v.reason =
            "insert matches the expected antibody and reads in-frame through the constant region"
                .into();
        v
    }

    /// Does `observed_aa` match a *different* panel member's V better than the
    /// expected one? Returns the suspected ab_id.
    fn detect_swap(&self, observed_aa: &str, parse: &NameParse, vector: &Vector) -> Option<String> {
        let s = seq::clean(&vector.sequence);
        let leader_len = seq::translate(&s[..vector.insertion_site.oh5_end.min(s.len())]).len();
        let obs_v: String = observed_aa.chars().skip(leader_len).collect();

        let mut best: Option<(String, f64)> = None;
        for (ab_id, row) in &self.ground_truth {
            let truth = match parse.chain_class {
                ChainClass::Heavy => row.heavy.as_deref(),
                _ => row.light.as_deref(),
            };
            let Some(t) = truth else { continue };
            let truth_aa = if seq::detect_type(t) == seq::SeqType::Nt {
                seq::translate_to_stop(t)
            } else {
                seq::clean(t)
            };
            let al = align::align_protein(&obs_v, &truth_aa, &self.library.alignment);
            let id = al.identity();
            if best.as_ref().map(|(_, b)| id > *b).unwrap_or(true) {
                best = Some((ab_id.clone(), id));
            }
        }
        match best {
            Some((ab_id, id)) if id >= 0.95 && !ab_id.eq_ignore_ascii_case(&parse.ab_id) => {
                Some(ab_id)
            }
            _ => None,
        }
    }
}

/// The chain class to key order cores by: a `light` parse resolves to κ/λ.
fn effective_class(class: ChainClass, locus: Option<Locus>) -> ChainClass {
    match class {
        ChainClass::Light => match locus {
            Some(Locus::Igk) => ChainClass::Kappa,
            Some(Locus::Igl) => ChainClass::Lambda,
            _ => ChainClass::Light,
        },
        other => other,
    }
}

/// Best ungapped identity of `landmark` against any circular offset of `text`.
fn best_ungapped_identity_circular(text: &str, landmark: &str) -> f64 {
    if landmark.is_empty() || text.is_empty() {
        return 0.0;
    }
    if landmark.len() > text.len() {
        return best_ungapped_identity(text, landmark);
    }
    let doubled = format!("{text}{text}");
    let t = doubled.as_bytes();
    let l = landmark.as_bytes();
    let mut best = 0usize;
    for off in 0..text.len() {
        let mut hit = 0;
        for k in 0..l.len() {
            if t[off + k] == l[k] {
                hit += 1;
            }
        }
        best = best.max(hit);
    }
    best as f64 / l.len() as f64
}

/// Translate the read in the frame that maximizes the run before the first stop,
/// used as a quick protein view for light-locus motif detection.
fn best_frame_protein(read: &str) -> String {
    let s = seq::clean(read);
    let mut best = String::new();
    for frame in 0..3 {
        let p = seq::translate_to_stop(&s[frame..]);
        if p.len() > best.len() {
            best = p;
        }
    }
    best
}

/// Locate the ORF starting at the leader and translate to the first stop.
fn locate_orf_protein(read: &str, leader_aa: &str) -> String {
    locate_orf_full(read, leader_aa)
        .split('*')
        .next()
        .unwrap_or("")
        .to_string()
}

/// Locate the ORF starting at the leader and translate full length (keeping `*`).
fn locate_orf_full(read: &str, leader_aa: &str) -> String {
    let s = seq::clean(read);
    // Try all three frames and both the read and its start; pick the frame whose
    // translation contains the leader peptide (or its prefix) earliest.
    let probe = &leader_aa[..leader_aa.len().min(8)];
    let mut best: Option<(usize, String)> = None;
    for frame in 0..3 {
        if frame >= s.len() {
            continue;
        }
        let prot = seq::translate(&s[frame..]);
        if let Some(pos) = prot.find(probe) {
            // Translate from the leader's start codon in this frame.
            let nt_start = frame + pos * 3;
            let from_leader = seq::translate(&s[nt_start..]);
            let score = pos;
            if best.as_ref().map(|(b, _)| score < *b).unwrap_or(true) {
                best = Some((score, from_leader));
            }
        }
    }
    match best {
        Some((_, p)) => p,
        // Fallback: frame 0 from the start.
        None => seq::translate(&s),
    }
}

/// Extract the observed insert nt window from a normalized read, given the
/// expected core length, by locating the leader and skipping to oh5_end.
fn extract_insert_nt(read: &str, vector: &Vector, core_len: usize) -> Option<String> {
    let s = seq::clean(read);
    let vs = seq::clean(&vector.sequence);
    let leader = &vs[..vector.insertion_site.oh5_end.min(vs.len())];
    // Find the leader (nt, codon-invariant) in the read.
    let probe = &leader[..leader.len().min(30)];
    let pos = s.find(probe)?;
    let insert_start = pos + vector.insertion_site.oh5_end;
    if insert_start + core_len <= s.len() {
        Some(s[insert_start..insert_start + core_len].to_string())
    } else {
        None
    }
}

/// Collect substitutions from an alignment, reporting positions relative to the
/// start of the variable region (i.e. after the leader).
fn collect_mutations(aligned_obs: &str, aligned_exp: &str, leader_len: usize) -> Vec<Mutation> {
    let mut muts = Vec::new();
    let mut exp_pos = 0usize; // residue index into the expected protein
    for (o, e) in aligned_obs.bytes().zip(aligned_exp.bytes()) {
        let exp_present = e != b'-';
        if exp_present {
            exp_pos += 1;
        }
        if o != b'-' && e != b'-' && o != e {
            // Report position within the V region (1-based).
            let pos_in_v = exp_pos.saturating_sub(leader_len);
            if pos_in_v >= 1 {
                muts.push(Mutation {
                    position_aa: pos_in_v,
                    wt: e as char,
                    mut_aa: o as char,
                });
            }
        }
    }
    muts
}

#[cfg(test)]
mod tests {
    // Gate-2 behavior is exercised end-to-end against assembled references in
    // the crate integration tests (tests/), where the real vector is available.
    use super::*;

    #[test]
    fn seqmode_parse() {
        assert_eq!(SeqMode::parse("partial_sanger"), SeqMode::PartialSanger);
        assert_eq!(SeqMode::parse("full_plasmid"), SeqMode::FullPlasmid);
        assert_eq!(SeqMode::parse(""), SeqMode::FullPlasmid);
    }
}
