//! `clonodoc-cli` — a thin headless runner over `clonodoc-core`.
//!
//! Usage:
//!   clonodoc-cli gate1 --library <lib.json5> --order <order.fasta|xlsx>
//!                     [--ground-truth <panel.csv|xlsx>] [--heavy-vector <id>]
//!                     [--overhang-set <id>] [--no-overhangs]
//!                     [--csv <out.csv>] [--html <out.html>]
//!   clonodoc-cli gate2 --library <lib.json5> --order <order.fasta> --reads <reads.fasta>
//!                     [--ground-truth <panel>] [--partial-sanger]
//!                     [--csv <out.csv>] [--html <out.html>]
//!   clonodoc-cli inspect-vector <vector.gb>
//!
//! It builds an in-memory project (heavy vector assigned, first overhang set /
//! naming profile) so a lab member can run a QC gate from one command.

use anyhow::{anyhow, bail, Context, Result};
use clonodoc_core::gate1::Gate1Context;
use clonodoc_core::gate2::{Gate2Context, SeqMode};
use clonodoc_core::model::{ChainClass, Library, Project};
use clonodoc_core::report;
use clonodoc_core::seqio::{self, fasta, genbank, GroundTruthRow, SeqRecord};
use clonodoc_core::{naming, workflow};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

struct Args {
    map: BTreeMap<String, String>,
    flags: Vec<String>,
    positional: Vec<String>,
}

fn parse_args(raw: &[String]) -> Args {
    let mut map = BTreeMap::new();
    let mut flags = Vec::new();
    let mut positional = Vec::new();
    let mut i = 0;
    while i < raw.len() {
        let a = &raw[i];
        if let Some(key) = a.strip_prefix("--") {
            // boolean flags
            if matches!(key, "no-overhangs" | "partial-sanger") {
                flags.push(key.to_string());
                i += 1;
            } else if i + 1 < raw.len() {
                map.insert(key.to_string(), raw[i + 1].clone());
                i += 2;
            } else {
                flags.push(key.to_string());
                i += 1;
            }
        } else {
            positional.push(a.clone());
            i += 1;
        }
    }
    Args {
        map,
        flags,
        positional,
    }
}

impl Args {
    fn get(&self, k: &str) -> Option<&str> {
        self.map.get(k).map(|s| s.as_str())
    }
    fn require(&self, k: &str) -> Result<&str> {
        self.get(k).ok_or_else(|| anyhow!("missing required --{k}"))
    }
    fn flag(&self, k: &str) -> bool {
        self.flags.iter().any(|f| f == k)
    }
}

fn main() {
    if let Err(e) = run() {
        eprintln!("error: {e:#}");
        std::process::exit(1);
    }
}

fn run() -> Result<()> {
    let argv: Vec<String> = std::env::args().skip(1).collect();
    let Some((cmd, rest)) = argv.split_first() else {
        print_usage();
        return Ok(());
    };
    let args = parse_args(rest);
    match cmd.as_str() {
        "gate1" => cmd_gate1(&args),
        "gate2" => cmd_gate2(&args),
        "inspect-vector" => cmd_inspect_vector(&args),
        "-h" | "--help" | "help" => {
            print_usage();
            Ok(())
        }
        other => {
            print_usage();
            bail!("unknown command '{other}'");
        }
    }
}

fn print_usage() {
    eprintln!(
        "clonodoc-cli — ClonoDoc antibody cloning verifier\n\n\
         COMMANDS\n  \
         gate1  --library L --order O [--ground-truth G] [--heavy-vector ID] [--overhang-set ID] [--no-overhangs] [--csv F] [--html F]\n  \
         gate2  --library L --order O --reads R [--ground-truth G] [--partial-sanger] [--csv F] [--html F]\n  \
         inspect-vector <vector.gb>\n"
    );
}

/// Build an in-memory project from the library and CLI flags.
fn build_project(lib: &Library, args: &Args) -> Result<Project> {
    if lib.overhang_sets.is_empty() {
        bail!("library has no overhang sets");
    }
    Ok(workflow::ad_hoc_project(
        lib,
        args.get("heavy-vector"),
        args.get("overhang-set"),
    ))
}

fn load_library(path: &str) -> Result<Library> {
    Library::load(path).with_context(|| format!("loading library {path}"))
}

fn load_ground_truth(args: &Args) -> Result<Vec<GroundTruthRow>> {
    match args.get("ground-truth") {
        None => Ok(Vec::new()),
        Some(p) => {
            let path = Path::new(p);
            let fmt = seqio_format(path);
            if fmt == "fasta" {
                // FASTA panel: id carries ab_id + chain token.
                let recs = fasta::read_path(path)?;
                let prof = naming::default_profile();
                let mut map: BTreeMap<String, GroundTruthRow> = BTreeMap::new();
                for r in recs {
                    let np = naming::parse_name(&r.id, &prof);
                    let row = map
                        .entry(np.ab_id.clone())
                        .or_insert_with(|| GroundTruthRow {
                            ab_id: np.ab_id.clone(),
                            heavy: None,
                            light: None,
                        });
                    match np.chain_class {
                        ChainClass::Heavy => row.heavy = Some(r.sequence),
                        _ => row.light = Some(r.sequence),
                    }
                }
                Ok(map.into_values().collect())
            } else {
                seqio::tabular::read_ground_truth_table(path, &BTreeMap::new())
                    .map_err(|e| anyhow!("{e}"))
            }
        }
    }
}

fn seqio_format(path: &Path) -> String {
    match path
        .extension()
        .and_then(|e| e.to_str())
        .map(|e| e.to_ascii_lowercase())
        .as_deref()
    {
        Some("fasta") | Some("fa") | Some("fna") => "fasta".into(),
        Some("xlsx") | Some("xls") => "xlsx".into(),
        Some("csv") => "csv".into(),
        _ => "fasta".into(),
    }
}

fn load_order(path: &str) -> Result<Vec<SeqRecord>> {
    let p = Path::new(path);
    seqio::read_id_seq_auto(p, "").map_err(|e| anyhow!("{e}"))
}

fn cmd_gate1(args: &Args) -> Result<()> {
    let lib = load_library(args.require("library")?)?;
    let project = build_project(&lib, args)?;
    let set = lib
        .overhang_set(&project.overhang_set)
        .ok_or_else(|| anyhow!("overhang set '{}' not found", project.overhang_set))?;
    let order = load_order(args.require("order")?)?;
    let ground_truth = load_ground_truth(args)?;
    let has_overhangs = !args.flag("no-overhangs");

    let ctx = Gate1Context::new(&lib, &project, set, &ground_truth, has_overhangs);
    let verdicts = ctx.run(&order);

    println!("Gate 1 — Order QC ({} records)\n", verdicts.len());
    print_gate1_table(&verdicts);
    let rollup = report::rollup_gate1(&verdicts);
    let passing = rollup.iter().filter(|r| r.passed).count();
    println!(
        "\nAntibodies: {} of {} pass all chains",
        passing,
        rollup.len()
    );

    if let Some(csv) = args.get("csv") {
        std::fs::write(csv, report::gate1_csv(&verdicts))?;
        println!("wrote {csv}");
    }
    if let Some(html) = args.get("html") {
        std::fs::write(html, report::html_report(&project.name, &verdicts, &[]))?;
        println!("wrote {html}");
    }
    Ok(())
}

fn cmd_gate2(args: &Args) -> Result<()> {
    let lib = load_library(args.require("library")?)?;
    let project = build_project(&lib, args)?;
    let set = lib
        .overhang_set(&project.overhang_set)
        .ok_or_else(|| anyhow!("overhang set '{}' not found", project.overhang_set))?;
    let ground_truth = load_ground_truth(args)?;

    let cores = match args.get("order") {
        Some(o) => {
            let order = load_order(o)?;
            workflow::order_cores(&order, &lib, &project, !args.flag("no-overhangs"))
        }
        None => BTreeMap::new(),
    };
    let reads = load_order(args.require("reads")?)?;
    let mode = if args.flag("partial-sanger") {
        SeqMode::PartialSanger
    } else {
        SeqMode::FullPlasmid
    };

    let ctx = Gate2Context::new(&lib, &project, set, &ground_truth, cores, mode);
    let verdicts = ctx.run(&reads);

    println!("Gate 2 — Sequencing QC ({} reads)\n", verdicts.len());
    print_gate2_table(&verdicts);

    if let Some(csv) = args.get("csv") {
        std::fs::write(csv, report::gate2_csv(&verdicts))?;
        println!("wrote {csv}");
    }
    if let Some(html) = args.get("html") {
        std::fs::write(html, report::html_report(&project.name, &[], &verdicts))?;
        println!("wrote {html}");
    }
    Ok(())
}

fn cmd_inspect_vector(args: &Args) -> Result<()> {
    let path = args
        .positional
        .first()
        .map(PathBuf::from)
        .ok_or_else(|| anyhow!("usage: inspect-vector <vector.gb>"))?;
    let rec = genbank::read_path(&path).map_err(|e| anyhow!("{e}"))?;
    println!(
        "{}: {} bp, {}",
        rec.name,
        rec.length,
        if rec.circular { "circular" } else { "linear" }
    );
    let roles = genbank::map_roles(&rec);
    for (role, f) in &roles {
        println!(
            "  {:16} {}..{}  {}",
            role,
            f.start,
            f.end,
            f.name.as_deref().unwrap_or("")
        );
    }
    Ok(())
}

fn print_gate1_table(verdicts: &[clonodoc_core::verdict::Gate1Verdict]) {
    println!(
        "{:<48} {:<8} {:<22} {:<8}",
        "RECORD", "CHAIN", "VERDICT", "READTHRU"
    );
    for v in verdicts {
        println!(
            "{:<48} {:<8} {:<22} {:<8}",
            truncate(&v.record_id, 48),
            v.chain_class,
            v.kind.label(),
            v.reads_through
                .map(|b| b.to_string())
                .unwrap_or_else(|| "-".into()),
        );
    }
}

fn print_gate2_table(verdicts: &[clonodoc_core::verdict::Gate2Verdict]) {
    println!(
        "{:<40} {:<8} {:<20} {:<14} {:<8}",
        "RECORD", "CHAIN", "VERDICT", "BACKBONE", "ID%"
    );
    for v in verdicts {
        println!(
            "{:<40} {:<8} {:<20} {:<14} {:<8}",
            truncate(&v.record_id, 40),
            v.chain_class,
            v.kind.label(),
            truncate(v.backbone_vector.as_deref().unwrap_or("-"), 14),
            v.backbone_identity
                .map(|x| format!("{:.1}", x * 100.0))
                .unwrap_or_else(|| "-".into()),
        );
    }
}

fn truncate(s: &str, n: usize) -> String {
    if s.len() <= n {
        s.to_string()
    } else {
        format!("{}…", &s[..n - 1])
    }
}
