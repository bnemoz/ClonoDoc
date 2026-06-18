//! Golden integration tests against the real fixtures in `reference/` and
//! `test_data/`. These encode the empirically-verified facts in
//! `reference/verified_facts.md` and must never regress.

use clonodoc_core::assemble;
use clonodoc_core::gate1::Gate1Context;
use clonodoc_core::gate2::{Gate2Context, SeqMode};
use clonodoc_core::model::Library;
use clonodoc_core::model::{ChainClass, Locus, Project};
use clonodoc_core::seq;
use clonodoc_core::seqio::{self, fasta, genbank, SeqRecord};
use clonodoc_core::verdict::{Gate1Kind, Gate2Kind};
use std::collections::BTreeMap;
use std::path::PathBuf;

fn root() -> PathBuf {
    // crates/core -> workspace root is two levels up.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
}

fn load_library() -> Library {
    let text = std::fs::read_to_string(root().join("reference/example_library.json5")).unwrap();
    Library::from_json5(&text).unwrap()
}

fn order_records() -> Vec<SeqRecord> {
    fasta::read_path(&root().join("test_data/IDT_ordered_sequences_correct.fasta")).unwrap()
}

fn test_project(lib: &Library) -> Project {
    let mut vector_assignments = BTreeMap::new();
    vector_assignments.insert("heavy".to_string(), "french_igg1_heavy".to_string());
    Project {
        schema_version: 1,
        name: "test".into(),
        created: "2026-06-17".into(),
        description: String::new(),
        naming_profile: lib.naming_profiles[0].id.clone(),
        overhang_set: lib.overhang_sets[0].id.clone(),
        vector_assignments,
        ground_truth: None,
        order: None,
        sequencing: None,
        name_overrides: BTreeMap::new(),
    }
}

// --- Phase 1: parsing -------------------------------------------------------

#[test]
fn genbank_vector_parses_to_verified_facts() {
    let rec =
        genbank::read_path(&root().join("test_data/IgG1_heavy_chain_french_vector.gb")).unwrap();
    assert_eq!(rec.length, 5738, "vector length");
    assert!(rec.circular, "topology circular");
    assert_eq!(rec.sequence.len(), 5738);

    let roles = genbank::map_roles(&rec);
    // sig_peptide 1..57; constant 58..1047; AmpR 2774..3634 (verified_facts §1).
    assert_eq!(roles["signal_peptide"].start, 1);
    assert_eq!(roles["signal_peptide"].end, 57);
    assert_eq!(roles["constant_region"].start, 58);
    assert_eq!(roles["constant_region"].end, 1047);
    assert_eq!(roles["resistance"].start, 2774);
    assert_eq!(roles["resistance"].end, 3634);
}

#[test]
fn order_fasta_has_42_records() {
    let recs = order_records();
    assert_eq!(recs.len(), 42);
    assert!(recs.iter().any(|r| r.id == "HVA-0195-r3-d02_heavy"));
    assert!(recs
        .iter()
        .any(|r| r.id == "UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4__heavy"));
}

#[test]
fn library_vector_sequence_hash_matches() {
    let lib = load_library();
    let v = lib.vector("french_igg1_heavy").unwrap();
    assert_eq!(v.length, 5738);
    assert_eq!(v.compute_sha256(), v.sequence_sha256);
}

#[test]
fn overhangs_map_to_vector_positions() {
    // verified_facts §2: oh_5[IGH] == vector 0-based 16..57; oh_3[IGH] == 57..85.
    let lib = load_library();
    let v = lib.vector("french_igg1_heavy").unwrap();
    let set = lib.overhang_set("french_default").unwrap();
    let vs = seq::clean(&v.sequence);
    assert_eq!(&vs[16..57], set.oh5(Locus::Igh).unwrap());
    assert_eq!(&vs[57..85], set.oh3(Locus::Igh).unwrap());
    // 5' shared between light loci; 3' distinct (the κ/λ discriminator).
    assert_eq!(set.oh5(Locus::Igk), set.oh5(Locus::Igl));
    assert_ne!(set.oh3(Locus::Igk), set.oh3(Locus::Igl));
}

#[test]
fn guided_genbank_import_builds_correct_vector() {
    // The first-use library builder: import the .gb + the bundled overhang set,
    // and it must reproduce the verified insertion site and constant anchor.
    let lib = load_library();
    let set = lib.overhang_set("french_default").unwrap();
    let gb =
        genbank::read_path(&root().join("test_data/IgG1_heavy_chain_french_vector.gb")).unwrap();
    let v = clonodoc_core::workflow::vector_from_genbank(
        &gb,
        "my_heavy",
        "My Heavy Vector",
        ChainClass::Heavy,
        "IgG1",
        set,
        Locus::Igh,
        "tester",
    );
    assert_eq!(v.length, 5738);
    assert!(v.topology.eq_ignore_ascii_case("circular"));
    // Overhangs located in the vector → insertion site 57/57 (verified_facts §3).
    assert_eq!(v.insertion_site.oh5_end, 57);
    assert_eq!(v.insertion_site.oh3_start, 57);
    assert_eq!(v.constant_anchor_aa, "ASTKGPSV");
    assert_eq!(v.compute_sha256(), v.sequence_sha256);

    // And the built vector is usable: a 378-core assembles to a productive ORF.
    let core = heavy_core_378();
    let orf = seq::translate(&assemble::assemble(&v, &core));
    let stop = orf.find('*').unwrap();
    assert_eq!(stop, 475);
    assert!(orf[..stop].contains(&v.constant_anchor_aa));
}

// --- Phase 3: the two golden assembly fixtures ------------------------------

fn heavy_core_378() -> String {
    let recs = order_records();
    let r = recs
        .iter()
        .find(|r| r.id == "HVA-0195-r3-d02_heavy")
        .unwrap();
    let lib = load_library();
    let set = lib.overhang_set("french_default").unwrap();
    let stripped = assemble::strip_overhangs(
        &r.sequence,
        set.oh5(Locus::Igh).unwrap(),
        set.oh3(Locus::Igh).unwrap(),
        1,
    );
    assert!(stripped.oh5_present && stripped.oh3_present);
    assert_eq!(stripped.core.len(), 378);
    stripped.core
}

#[test]
fn golden_pass_fixture_reads_through() {
    let lib = load_library();
    let v = lib.vector("french_igg1_heavy").unwrap();
    let core = heavy_core_378();
    let assembled = assemble::assemble(v, &core);
    let orf = seq::translate(&assembled);
    let first_stop = orf.find('*').unwrap();
    assert_eq!(first_stop, 475, "natural stop at aa 475");
    assert!(
        orf[..first_stop].contains("ASTKGPSV"),
        "reads through into constant"
    );
    assert!(orf.starts_with("MELGLRWVFLVAILEGVQCQVQLVESGGG"));
}

#[test]
fn golden_fail_fixture_frameshifts_at_146() {
    let lib = load_library();
    let v = lib.vector("french_igg1_heavy").unwrap();
    // One stray trailing base => 379-nt core (the verified Gate-1 catch).
    let core_bad = format!("{}G", heavy_core_378());
    assert_eq!(core_bad.len(), 379);
    let assembled = assemble::assemble(v, &core_bad);
    let orf = seq::translate(&assembled);
    let first_stop = orf.find('*').unwrap();
    assert_eq!(first_stop, 146, "premature stop at aa 146");
    assert!(!orf[..first_stop].contains("ASTKGPSV"), "no read-through");
    // The V domain protein is identical to the PASS case up to the junction.
    let good = seq::translate(&assemble::assemble(v, &heavy_core_378()));
    assert!(good.starts_with(&orf[..140]));
}

// --- Phase 4: Gate 1 --------------------------------------------------------

#[test]
fn gate1_correct_order_all_productive() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    let ctx = Gate1Context::new(&lib, &project, set, &[], true);
    let verdicts = ctx.run(&order_records());
    assert_eq!(verdicts.len(), 42);
    // Every record passes (advisory NO_GROUND_TRUTH; no FAILs).
    for v in &verdicts {
        assert!(
            v.passed(),
            "record {} should pass but got {:?}: {}",
            v.record_id,
            v.kind,
            v.reason
        );
    }
    // Heavy records (with the assigned vector) must read through into the constant.
    let heavy = verdicts
        .iter()
        .find(|v| v.record_id == "HVA-0195-r3-d02_heavy")
        .unwrap();
    assert_eq!(heavy.reads_through, Some(true));
    assert_eq!(heavy.core_len, Some(378));
}

#[test]
fn gate1_flags_379_frameshift() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // Reconstruct a wrong order: oh5 + 379core + oh3.
    let core_bad = format!("{}G", heavy_core_378());
    let frag = format!(
        "{}{}{}",
        set.oh5(Locus::Igh).unwrap(),
        core_bad,
        set.oh3(Locus::Igh).unwrap()
    );
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: frag,
    };
    let ctx = Gate1Context::new(&lib, &project, set, &[], true);
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate1Kind::FrameshiftAtJunction);
    assert_eq!(v.reads_through, Some(false));
    assert_eq!(v.premature_stop_aa, Some(146));
}

#[test]
fn gate1_wrong_locus_overhang() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // A heavy-named record carrying κ overhangs.
    let core = heavy_core_378();
    let frag = format!(
        "{}{}{}",
        set.oh5(Locus::Igk).unwrap(),
        core,
        set.oh3(Locus::Igk).unwrap()
    );
    let rec = SeqRecord {
        id: "mab9_heavy".into(),
        sequence: frag,
    };
    let ctx = Gate1Context::new(&lib, &project, set, &[], true);
    let v = &ctx.run(&[rec])[0];
    assert!(
        matches!(
            v.kind,
            Gate1Kind::OverhangWrongLocus | Gate1Kind::OverhangMissing
        ),
        "got {:?}: {}",
        v.kind,
        v.reason
    );
}

// --- Phase 5: Gate 2 (synthetic reads from the assembled reference) ---------

fn full_plasmid_read(core: &str) -> String {
    let lib = load_library();
    let v = lib.vector("french_igg1_heavy").unwrap();
    assemble::assemble(v, core)
}

fn gate2_ctx<'a>(
    lib: &'a Library,
    project: &'a Project,
    set: &'a clonodoc_core::model::OverhangSet,
    cores: BTreeMap<(String, ChainClass), String>,
    mode: SeqMode,
) -> Gate2Context<'a> {
    Gate2Context::new(lib, project, set, &[], cores, mode)
}

fn heavy_cores_map() -> BTreeMap<(String, ChainClass), String> {
    let mut m = BTreeMap::new();
    m.insert(
        ("HVA-0195-R3-D02".to_string(), ChainClass::Heavy),
        heavy_core_378(),
    );
    m
}

#[test]
fn gate2_perfect_read_passes() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    let read = full_plasmid_read(&heavy_core_378());
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::Pass, "reason: {}", v.reason);
    assert_eq!(v.reads_through, Some(true));
}

#[test]
fn gate2_reverse_complement_rotated_still_passes() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // Rotate the circular read and reverse-complement it.
    let assembled = full_plasmid_read(&heavy_core_378());
    let rotated = format!("{}{}", &assembled[1500..], &assembled[..1500]);
    let rc = seq::revcomp(&rotated);
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: rc,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::Pass, "reason: {}", v.reason);
}

#[test]
fn gate2_point_mutation_detected() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // Substitute one codon inside the V region (change codon ~30 to GCT=Ala).
    let mut core: Vec<u8> = heavy_core_378().into_bytes();
    let pos = 30 * 3; // codon 30
    core[pos] = b'G';
    core[pos + 1] = b'C';
    core[pos + 2] = b'T';
    let core = String::from_utf8(core).unwrap();
    let read = full_plasmid_read(&core);
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::PointMutation, "reason: {}", v.reason);
    assert!(!v.mutations.is_empty());
}

#[test]
fn gate2_junction_frameshift_detected() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // Delete one base near the 3' junction => frameshift, no read-through.
    let mut core = heavy_core_378();
    core.pop(); // 377 nt
    let read = full_plasmid_read(&core);
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(
        v.kind,
        Gate2Kind::JunctionFrameshift,
        "reason: {}",
        v.reason
    );
    assert_eq!(v.reads_through, Some(false));
}

#[test]
fn gate2_empty_vector_detected() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // No insert at all: leader fused straight into the constant region.
    let read = full_plasmid_read("");
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::EmptyVector, "reason: {}", v.reason);
}

#[test]
fn gate2_silent_variant_detected() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();
    // Make a synonymous change: find a Leu codon CTG and change to CTC (still Leu).
    let core = heavy_core_378();
    let idx = core.find("CTG").expect("a CTG codon exists");
    // ensure codon-aligned
    let idx = idx - (idx % 3);
    let mut bytes = core.clone().into_bytes();
    if &bytes[idx..idx + 3] == b"CTG" {
        bytes[idx + 2] = b'C'; // CTG -> CTC, both Leu
    }
    let mutated = String::from_utf8(bytes).unwrap();
    // The expected core (in the order map) is the original; the read carries the
    // synonymous change.
    let read = full_plasmid_read(&mutated);
    let ctx = gate2_ctx(&lib, &project, set, heavy_cores_map(), SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert!(
        matches!(v.kind, Gate2Kind::SilentVariant | Gate2Kind::Pass),
        "expected silent/pass, got {:?}: {}",
        v.kind,
        v.reason
    );
}

fn heavy_core_for(id: &str) -> String {
    let recs = order_records();
    let r = recs.iter().find(|r| r.id == id).unwrap();
    let lib = load_library();
    let set = lib.overhang_set("french_default").unwrap();
    let s = assemble::strip_overhangs(
        &r.sequence,
        set.oh5(Locus::Igh).unwrap(),
        set.oh3(Locus::Igh).unwrap(),
        1,
    );
    assert!(s.oh5_present && s.oh3_present);
    s.core
}

#[test]
fn gate2_wrong_vector_detected() {
    // Build a distinct κ "vector" by cloning the heavy backbone and mutating its
    // constant-region landmark, so backbone identification can tell them apart.
    let mut lib = load_library();
    let mut kappa = lib.vector("french_igg1_heavy").unwrap().clone();
    kappa.id = "fake_kappa".into();
    kappa.display_name = "Fake Kappa".into();
    kappa.chain_class = ChainClass::Kappa;
    kappa.overhang_locus = Locus::Igk;
    // Mutate 200 nt at the start of the constant region to a distinct landmark.
    let cstart = kappa.features["constant_region"].range0().start;
    let mut bytes = seq::clean(&kappa.sequence).into_bytes();
    for i in cstart..(cstart + 200).min(bytes.len()) {
        bytes[i] = if (i % 2) == 0 { b'A' } else { b'T' };
    }
    kappa.sequence = String::from_utf8(bytes).unwrap();
    kappa.sequence_sha256 = kappa.compute_sha256();
    lib.vectors.push(kappa.clone());

    let mut project = test_project(&lib);
    project
        .vector_assignments
        .insert("kappa".into(), "fake_kappa".into());
    let set = lib.overhang_set("french_default").unwrap().clone();

    // A heavy insert placed in the κ backbone, but named heavy.
    let read = assemble::assemble(&kappa, &heavy_core_378());
    let ctx = Gate2Context::new(
        &lib,
        &project,
        &set,
        &[],
        heavy_cores_map(),
        SeqMode::FullPlasmid,
    );
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::WrongVector, "reason: {}", v.reason);
    assert_eq!(v.backbone_vector.as_deref(), Some("fake_kappa"));
}

#[test]
fn gate2_wrong_insert_swap_detected() {
    let lib = load_library();
    let project = test_project(&lib);
    let set = lib.overhang_set("french_default").unwrap();

    let core_a = heavy_core_378(); // HVA-0195-r3-d02
    let core_b = heavy_core_for("HVA-0194-r3-d09_heavy");

    // Panel knows both antibodies' heavy V (AA).
    let gt = vec![
        GroundTruthRow {
            ab_id: "HVA-0195-r3-d02".into(),
            heavy: Some(seq::translate_to_stop(&core_a)),
            light: None,
        },
        GroundTruthRow {
            ab_id: "HVA-0194-r3-d09".into(),
            heavy: Some(seq::translate_to_stop(&core_b)),
            light: None,
        },
    ];

    // Order map expects A's core for the A record, but the read carries B's insert.
    let mut cores = BTreeMap::new();
    cores.insert(("HVA-0195-R3-D02".to_string(), ChainClass::Heavy), core_a);

    let read = full_plasmid_read(&core_b);
    let ctx = Gate2Context::new(&lib, &project, set, &gt, cores, SeqMode::FullPlasmid);
    let rec = SeqRecord {
        id: "HVA-0195-r3-d02_heavy".into(),
        sequence: read,
    };
    let v = &ctx.run(&[rec])[0];
    assert_eq!(v.kind, Gate2Kind::WrongInsertSwap, "reason: {}", v.reason);
    assert_eq!(v.suspected_identity.as_deref(), Some("HVA-0194-R3-D09"));
}

// keep a reference to seqio so unused-import lints stay quiet if helpers change
use seqio::GroundTruthRow;
