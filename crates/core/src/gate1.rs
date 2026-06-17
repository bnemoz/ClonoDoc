//! Gate 1 — in-silico order QC, pre-cloning (`docs/01_DESIGN.md` §2).
//!
//! For each ordered record: parse the name → locus; strip & locus-check the
//! overhangs; translate the core and (if a panel is loaded) compare to ground
//! truth at the AA level; assemble the full construct and confirm it is
//! **productive** (stop-free read-through into the constant region). The
//! read-through check — not a V-only AA check — is the linchpin that catches the
//! 379-nt junction frameshift whose V domain looks perfect.

use crate::align;
use crate::assemble;
use crate::model::{ChainClass, Library, Locus, OverhangSet, Project, Vector};
use crate::naming::{self, NameParse};
use crate::seq;
use crate::seqio::{GroundTruthRow, SeqRecord};
use crate::verdict::{Gate1Kind, Gate1Verdict};
use std::collections::BTreeMap;

/// Everything Gate 1 needs, resolved from the library + project.
pub struct Gate1Context<'a> {
    pub library: &'a Library,
    pub project: &'a Project,
    pub overhang_set: &'a OverhangSet,
    /// ab_id (uppercased) → panel row.
    pub ground_truth: BTreeMap<String, GroundTruthRow>,
    pub has_overhangs: bool,
}

impl<'a> Gate1Context<'a> {
    pub fn new(
        library: &'a Library,
        project: &'a Project,
        overhang_set: &'a OverhangSet,
        ground_truth: &[GroundTruthRow],
        has_overhangs: bool,
    ) -> Self {
        let map = ground_truth
            .iter()
            .map(|r| (r.ab_id.to_ascii_uppercase(), r.clone()))
            .collect();
        Gate1Context {
            library,
            project,
            overhang_set,
            ground_truth: map,
            has_overhangs,
        }
    }

    /// Resolve the chain class for a record, honoring per-project name overrides.
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

    /// Resolve which library vector and locus apply to a record. For a `light`
    /// record the κ/λ locus is detected from the ordered fragment's 3′ overhang.
    fn resolve_vector(
        &self,
        class: ChainClass,
        fragment: &str,
    ) -> (Option<Locus>, Option<&'a Vector>) {
        let locus = match class {
            ChainClass::Heavy => Some(Locus::Igh),
            ChainClass::Kappa => Some(Locus::Igk),
            ChainClass::Lambda => Some(Locus::Igl),
            ChainClass::Light => assemble::detect_locus_from_overhangs(
                fragment,
                self.overhang_set,
                self.library.alignment.overhang_max_mismatch,
            )
            .or_else(|| assemble::detect_locus_from_protein(&seq::translate(fragment))),
            ChainClass::Unknown => None,
        };
        let vector = locus
            .map(|l| l.chain_class())
            .and_then(|c| self.project.vector_for(c))
            .and_then(|id| self.library.vector(id));
        (locus, vector)
    }

    /// Run Gate 1 over a set of ordered records.
    pub fn run(&self, records: &[SeqRecord]) -> Vec<Gate1Verdict> {
        let profile = self
            .library
            .naming_profile(&self.project.naming_profile)
            .cloned()
            .unwrap_or_else(naming::default_profile);
        records.iter().map(|r| self.run_one(r, &profile)).collect()
    }

    fn run_one(&self, rec: &SeqRecord, profile: &naming::NamingProfile) -> Gate1Verdict {
        let parse = self.resolve_parse(&rec.id, profile);
        let mut advisories = Vec::new();
        let mut verdict = Gate1Verdict {
            record_id: rec.id.clone(),
            ab_id: parse.ab_id.clone(),
            chain_class: parse.chain_class.as_str().to_string(),
            kind: Gate1Kind::Pass,
            advisories: Vec::new(),
            reason: String::new(),
            core_len: None,
            premature_stop_aa: None,
            reads_through: None,
        };

        // 1–2. Locus & overhang stripping.
        let (locus, vector) = self.resolve_vector(parse.chain_class, &rec.sequence);
        let core;
        if self.has_overhangs {
            let Some(locus) = locus else {
                verdict.kind = Gate1Kind::OverhangMissing;
                verdict.reason =
                    "could not determine locus, so overhangs could not be located".into();
                return verdict;
            };
            let oh5 = self.overhang_set.oh5(locus).unwrap_or("");
            let oh3 = self.overhang_set.oh3(locus).unwrap_or("");
            let stripped = assemble::strip_overhangs(
                &rec.sequence,
                oh5,
                oh3,
                self.library.alignment.overhang_max_mismatch,
            );
            if !stripped.oh5_present || !stripped.oh3_present {
                // Maybe the fragment carries a *different* locus' overhangs.
                let other = assemble::detect_locus_from_overhangs(
                    &rec.sequence,
                    self.overhang_set,
                    self.library.alignment.overhang_max_mismatch,
                );
                verdict.kind = match other {
                    Some(o)
                        if o.chain_class() != parse.chain_class
                            && parse.chain_class != ChainClass::Light =>
                    {
                        verdict.reason = format!(
                            "name says {} but the fragment carries {} overhangs",
                            parse.chain_class.as_str(),
                            o.key()
                        );
                        Gate1Kind::OverhangWrongLocus
                    }
                    _ => {
                        verdict.reason = format!(
                            "expected {} overhangs not found at both ends (5′ {}, 3′ {})",
                            locus.key(),
                            yn(stripped.oh5_present),
                            yn(stripped.oh3_present)
                        );
                        Gate1Kind::OverhangMissing
                    }
                };
                return verdict;
            }
            if stripped.oh5_mismatches > 0 || stripped.oh3_mismatches > 0 {
                verdict.reason = format!(
                    "overhang typo tolerated (5′ {} mm, 3′ {} mm); ",
                    stripped.oh5_mismatches, stripped.oh3_mismatches
                );
            }
            core = stripped.core;
        } else {
            // The fragment is already the bare core (e.g. an untrimmed xlsx order).
            core = seq::clean(&rec.sequence);
        }
        verdict.core_len = Some(core.len());

        // 3. Translate the core.
        let core_aa = seq::translate(&core);
        let core_aa_to_stop = seq::translate_to_stop(&core);

        // 4. Compare to ground truth, when available.
        let truth = self.ground_truth.get(&parse.ab_id.to_ascii_uppercase());
        let mut had_truth = false;
        if let Some(row) = truth {
            let truth_seq = match parse.chain_class {
                ChainClass::Heavy => row.heavy.as_deref(),
                _ => row.light.as_deref(),
            };
            if let Some(t) = truth_seq {
                had_truth = true;
                let truth_aa = if seq::detect_type(t) == seq::SeqType::Nt {
                    seq::translate_to_stop(t)
                } else {
                    seq::clean(t)
                };
                let al = align::align_protein(&core_aa_to_stop, &truth_aa, &self.library.alignment);
                if al.mismatches > 0 || al.gaps > 0 {
                    verdict.kind = Gate1Kind::TranslationDrift;
                    let pos = first_difference(&al.aligned_a, &al.aligned_b);
                    verdict.reason = format!(
                        "{}translated core differs from ground truth ({} mismatch(es), {} gap(s); first at aa {})",
                        verdict.reason,
                        al.mismatches,
                        al.gaps,
                        pos.map(|p| p.to_string()).unwrap_or_else(|| "?".into())
                    );
                    return verdict;
                }
            }
        }

        // 5. Assemble & test productivity (the read-through linchpin).
        if let Some(vector) = vector {
            let assembled = assemble::assemble(vector, &core);
            let orf = seq::translate(&assembled);
            let first_stop = orf.find('*');
            let prefix = match first_stop {
                Some(i) => &orf[..i],
                None => &orf[..],
            };
            let reads_through = prefix.contains(&vector.constant_anchor_aa);
            verdict.reads_through = Some(reads_through);
            verdict.premature_stop_aa = first_stop;
            if !reads_through {
                verdict.kind = Gate1Kind::FrameshiftAtJunction;
                verdict.reason = format!(
                    "{}V-domain translates but the assembled ORF hits a premature stop at aa {} and does not read through into the constant region ({} absent)",
                    verdict.reason,
                    first_stop.map(|i| i.to_string()).unwrap_or_else(|| "?".into()),
                    vector.constant_anchor_aa
                );
                return verdict;
            }
        } else {
            // No vector for this locus: fall back to a structural productivity check.
            verdict.reads_through = None;
            if core.len() % 3 != 0 {
                verdict.kind = Gate1Kind::FrameshiftAtJunction;
                verdict.reason = format!(
                    "{}core length {} is not a multiple of 3 (frameshift at junction)",
                    verdict.reason,
                    core.len()
                );
                return verdict;
            }
            if core_aa.contains('*') {
                verdict.kind = Gate1Kind::FrameshiftAtJunction;
                let pos = core_aa.find('*').unwrap();
                verdict.premature_stop_aa = Some(pos);
                verdict.reason = format!(
                    "{}premature stop at aa {} within the translated core",
                    verdict.reason, pos
                );
                return verdict;
            }
            verdict.reason = format!(
                "{}no vector assigned for this locus; productivity checked structurally only",
                verdict.reason
            );
        }

        // 6. Advisory optimization checks (never fail the order).
        if self.library.optimization_checks.enabled {
            let oc = &self.library.optimization_checks;
            let gc_too_high = seq::gc_fraction(&core) > oc.gc_max_global
                || seq::max_windowed_gc(&core, oc.gc_window) > oc.gc_max_window;
            if gc_too_high {
                advisories.push(Gate1Kind::GcWarning);
            }
        }

        // Final verdict shaping.
        verdict.advisories = advisories;
        if had_truth {
            verdict.kind = Gate1Kind::Pass;
            if verdict.reason.is_empty() {
                verdict.reason = "core translates identically to ground truth; assembled construct is productive (reads through into the constant region)".into();
            } else {
                verdict.reason = format!(
                    "{}core matches ground truth; construct productive",
                    verdict.reason
                );
            }
        } else {
            verdict.kind = Gate1Kind::NoGroundTruth;
            if verdict.reason.is_empty() {
                verdict.reason = "no ground-truth panel entry for this antibody; assembled construct is productive (advisory pass)".into();
            } else {
                verdict.reason = format!(
                    "{}construct productive; no ground-truth entry (advisory)",
                    verdict.reason
                );
            }
        }
        verdict
    }
}

fn yn(b: bool) -> &'static str {
    if b {
        "found"
    } else {
        "absent"
    }
}

/// First aligned column where the two aligned strings differ (1-based, counting
/// only residues of the first sequence).
fn first_difference(a: &str, b: &str) -> Option<usize> {
    let mut pos = 0usize;
    for (x, y) in a.bytes().zip(b.bytes()) {
        if x != b'-' {
            pos += 1;
        }
        if x != y {
            return Some(pos.max(1));
        }
    }
    None
}
