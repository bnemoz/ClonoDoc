//! Verdict reporting (`docs/01_DESIGN.md` §5): CSV and self-contained HTML, plus
//! the per-antibody rollup (an antibody PASSES only if all its chains pass).

use crate::naming;
use crate::verdict::{AntibodyRollup, Gate1Verdict, Gate2Verdict};
use std::collections::BTreeMap;

/// Roll Gate-1 verdicts up per antibody.
pub fn rollup_gate1(verdicts: &[Gate1Verdict]) -> Vec<AntibodyRollup> {
    let mut groups: BTreeMap<String, Vec<&Gate1Verdict>> = BTreeMap::new();
    for v in verdicts {
        groups.entry(v.ab_id.clone()).or_default().push(v);
    }
    groups
        .into_iter()
        .map(|(ab_id, vs)| {
            let passed = vs.iter().all(|v| v.passed());
            let chains = vs.iter().map(|v| v.record_id.clone()).collect();
            let note = if passed {
                format!("{} chain(s) pass", vs.len())
            } else {
                let bad: Vec<&str> = vs
                    .iter()
                    .filter(|v| !v.passed())
                    .map(|v| v.kind.label())
                    .collect();
                format!("flagged: {}", bad.join(", "))
            };
            AntibodyRollup {
                ab_id,
                passed,
                note,
                chains,
            }
        })
        .collect()
}

/// Roll Gate-2 verdicts up per antibody, flagging incomplete pairs.
pub fn rollup_gate2(
    verdicts: &[Gate2Verdict],
    expected_chains: &BTreeMap<String, usize>,
) -> Vec<AntibodyRollup> {
    let mut groups: BTreeMap<String, Vec<&Gate2Verdict>> = BTreeMap::new();
    for v in verdicts {
        groups.entry(v.ab_id.clone()).or_default().push(v);
    }
    groups
        .into_iter()
        .map(|(ab_id, vs)| {
            let all_pass = vs.iter().all(|v| v.passed());
            let expected = expected_chains.get(&ab_id).copied().unwrap_or(vs.len());
            let complete = vs.len() >= expected;
            let passed = all_pass && complete;
            let chains = vs.iter().map(|v| v.record_id.clone()).collect();
            let note = if !complete {
                format!(
                    "INCOMPLETE_PAIR ({} of {} chains present)",
                    vs.len(),
                    expected
                )
            } else if passed {
                format!("{} chain(s) pass", vs.len())
            } else {
                let bad: Vec<&str> = vs
                    .iter()
                    .filter(|v| !v.passed())
                    .map(|v| v.kind.label())
                    .collect();
                format!("flagged: {}", bad.join(", "))
            };
            AntibodyRollup {
                ab_id,
                passed,
                note,
                chains,
            }
        })
        .collect()
}

fn csv_escape(s: &str) -> String {
    if s.contains([',', '"', '\n']) {
        format!("\"{}\"", s.replace('"', "\"\""))
    } else {
        s.to_string()
    }
}

/// Gate-1 verdicts → CSV.
pub fn gate1_csv(verdicts: &[Gate1Verdict]) -> String {
    let mut out = String::from(
        "record_id,ab_id,chain_class,verdict,severity,core_len,premature_stop_aa,reads_through,reason\n",
    );
    for v in verdicts {
        out.push_str(&format!(
            "{},{},{},{},{:?},{},{},{},{}\n",
            csv_escape(&v.record_id),
            csv_escape(&v.ab_id),
            csv_escape(&v.chain_class),
            v.kind.label(),
            v.severity(),
            v.core_len.map(|x| x.to_string()).unwrap_or_default(),
            v.premature_stop_aa
                .map(|x| x.to_string())
                .unwrap_or_default(),
            v.reads_through.map(|x| x.to_string()).unwrap_or_default(),
            csv_escape(&v.reason),
        ));
    }
    out
}

/// Gate-2 verdicts → CSV.
pub fn gate2_csv(verdicts: &[Gate2Verdict]) -> String {
    let mut out = String::from(
        "record_id,ab_id,chain_class,verdict,severity,backbone_vector,backbone_identity,backbone_observed,reads_through,suspected_identity,reason\n",
    );
    for v in verdicts {
        out.push_str(&format!(
            "{},{},{},{},{:?},{},{},{},{},{},{}\n",
            csv_escape(&v.record_id),
            csv_escape(&v.ab_id),
            csv_escape(&v.chain_class),
            v.kind.label(),
            v.severity(),
            csv_escape(v.backbone_vector.as_deref().unwrap_or("")),
            v.backbone_identity
                .map(|x| format!("{:.4}", x))
                .unwrap_or_default(),
            csv_escape(&v.backbone_observed),
            v.reads_through.map(|x| x.to_string()).unwrap_or_default(),
            csv_escape(v.suspected_identity.as_deref().unwrap_or("")),
            csv_escape(&v.reason),
        ));
    }
    out
}

fn html_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn verdict_color(label: &str) -> &'static str {
    match label {
        "PASS" => "#1a7f37",
        "NO_GROUND_TRUTH" | "SILENT_VARIANT" | "GC_WARNING" | "RARE_CODON_WARNING" => "#9a6700",
        _ => "#cf222e",
    }
}

/// A self-contained HTML report for a project (`docs/01_DESIGN.md` §5).
pub fn html_report(project_name: &str, gate1: &[Gate1Verdict], gate2: &[Gate2Verdict]) -> String {
    let mut s = String::new();
    s.push_str(&format!(
        "<!doctype html><html><head><meta charset=\"utf-8\"><title>abclone-verify — {}</title>",
        html_escape(project_name)
    ));
    s.push_str(
        "<style>body{font-family:-apple-system,Segoe UI,Roboto,sans-serif;margin:2rem;color:#1f2328}\
         h1{font-size:1.4rem}h2{margin-top:2rem;border-bottom:1px solid #d0d7de;padding-bottom:.3rem}\
         table{border-collapse:collapse;width:100%;font-size:.85rem;margin-top:.5rem}\
         th,td{border:1px solid #d0d7de;padding:.35rem .5rem;text-align:left;vertical-align:top}\
         th{background:#f6f8fa}.v{font-weight:600}.reason{color:#57606a}\
         .summary{display:flex;gap:1rem;flex-wrap:wrap;margin:.5rem 0}\
         .chip{padding:.2rem .6rem;border-radius:1rem;color:#fff;font-size:.8rem}</style></head><body>",
    );
    s.push_str(&format!(
        "<h1>abclone-verify report — {}</h1>",
        html_escape(project_name)
    ));

    if !gate1.is_empty() {
        s.push_str("<h2>Gate 1 — Order QC (pre-cloning)</h2>");
        s.push_str(&summary_chips_g1(gate1));
        s.push_str("<table><tr><th>Record</th><th>Antibody</th><th>Chain</th><th>Verdict</th><th>Core len</th><th>Reads through</th><th>Reason</th></tr>");
        for v in gate1 {
            s.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td class=\"v\" style=\"color:{}\">{}</td><td>{}</td><td>{}</td><td class=\"reason\">{}</td></tr>",
                html_escape(&v.record_id),
                html_escape(&v.ab_id),
                html_escape(&v.chain_class),
                verdict_color(v.kind.label()),
                v.kind.label(),
                v.core_len.map(|x| x.to_string()).unwrap_or_default(),
                v.reads_through.map(|x| x.to_string()).unwrap_or_default(),
                html_escape(&v.reason),
            ));
        }
        s.push_str("</table>");
    }

    if !gate2.is_empty() {
        s.push_str("<h2>Gate 2 — Sequencing QC (post-cloning)</h2>");
        s.push_str("<table><tr><th>Record</th><th>Antibody</th><th>Chain</th><th>Verdict</th><th>Backbone</th><th>Identity</th><th>Observed</th><th>Reason</th></tr>");
        for v in gate2 {
            s.push_str(&format!(
                "<tr><td>{}</td><td>{}</td><td>{}</td><td class=\"v\" style=\"color:{}\">{}</td><td>{}</td><td>{}</td><td>{}</td><td class=\"reason\">{}</td></tr>",
                html_escape(&v.record_id),
                html_escape(&v.ab_id),
                html_escape(&v.chain_class),
                verdict_color(v.kind.label()),
                v.kind.label(),
                html_escape(v.backbone_vector.as_deref().unwrap_or("—")),
                v.backbone_identity.map(|x| format!("{:.1}%", x * 100.0)).unwrap_or_default(),
                html_escape(&v.backbone_observed),
                html_escape(&v.reason),
            ));
        }
        s.push_str("</table>");
    }

    s.push_str("<p class=\"reason\" style=\"margin-top:2rem\">Generated by abclone-verify.</p>");
    s.push_str("</body></html>");
    s
}

fn summary_chips_g1(verdicts: &[Gate1Verdict]) -> String {
    let mut counts: BTreeMap<&str, usize> = BTreeMap::new();
    for v in verdicts {
        *counts.entry(v.kind.label()).or_default() += 1;
    }
    let mut s = String::from("<div class=\"summary\">");
    for (label, n) in counts {
        s.push_str(&format!(
            "<span class=\"chip\" style=\"background:{}\">{}: {}</span>",
            verdict_color(label),
            label,
            n
        ));
    }
    s.push_str("</div>");
    s
}

/// Expected chain counts per antibody from the parsed records (default: count of
/// distinct classes seen). Helper for [`rollup_gate2`].
pub fn expected_chain_counts(
    record_ids: &[String],
    profile: &naming::NamingProfile,
) -> BTreeMap<String, usize> {
    let mut map: BTreeMap<String, usize> = BTreeMap::new();
    for id in record_ids {
        let p = naming::parse_name(id, profile);
        *map.entry(p.ab_id).or_default() += 1;
    }
    map
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::verdict::Gate1Kind;

    fn g1(id: &str, ab: &str, kind: Gate1Kind) -> Gate1Verdict {
        Gate1Verdict {
            record_id: id.into(),
            ab_id: ab.into(),
            chain_class: "heavy".into(),
            kind,
            advisories: vec![],
            reason: "r".into(),
            core_len: Some(378),
            premature_stop_aa: None,
            reads_through: Some(true),
        }
    }

    #[test]
    fn rollup_requires_all_chains_pass() {
        let vs = vec![
            g1("a_heavy", "a", Gate1Kind::Pass),
            g1("a_light", "a", Gate1Kind::FrameshiftAtJunction),
            g1("b_heavy", "b", Gate1Kind::Pass),
        ];
        let roll = rollup_gate1(&vs);
        let a = roll.iter().find(|r| r.ab_id == "a").unwrap();
        let b = roll.iter().find(|r| r.ab_id == "b").unwrap();
        assert!(!a.passed);
        assert!(b.passed);
    }

    #[test]
    fn csv_has_header_and_rows() {
        let vs = vec![g1("a_heavy", "a", Gate1Kind::Pass)];
        let csv = gate1_csv(&vs);
        assert!(csv.starts_with("record_id,"));
        assert!(csv.contains("a_heavy"));
        assert!(csv.contains("PASS"));
    }

    #[test]
    fn html_is_self_contained() {
        let vs = vec![g1("a_heavy", "a", Gate1Kind::Pass)];
        let html = html_report("demo", &vs, &[]);
        assert!(html.contains("<!doctype html>"));
        assert!(html.contains("a_heavy"));
    }
}
