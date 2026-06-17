# 02 — Configuration schema

Two scopes, matching the natural structure of the work:

- **Library (lab-global):** vectors + overhang sets + naming profiles + alignment
  settings. One file, shared and enriched by the whole lab.
- **Project (per-project):** the inputs and results for one cloning campaign.

Format is **`json5`** (JSON with comments and trailing commas) so the shared file can be
self-documenting. All examples below are valid json5.

---

## 1. Library file — `library.json5`

Default location (per-OS, via the `directories` crate):
- Linux: `~/.config/abclone-verify/library.json5`
- macOS: `~/Library/Application Support/abclone-verify/library.json5`
- Windows: `%APPDATA%\abclone-verify\library.json5`

A lab can instead point the app at a **shared** path (network drive, synced folder) so
everyone reads/writes the same library. The app must merge-import additively and dedupe
by sequence hash (§1.4).

```json5
{
  schema_version: 1,

  // ---- Vector backbones -------------------------------------------------
  vectors: [
    {
      id: "french_igg1_heavy",          // stable key, referenced by projects
      display_name: "French IgG1 Heavy Chain",
      chain_class: "heavy",             // heavy | kappa | lambda
      isotype: "IgG1",                  // free text shown in reports (IgG1/IgG4/IgA/kappa/lambda…)
      topology: "circular",             // circular | linear
      length: 5738,

      // Full vector nucleotide sequence (the ORIGIN). Stored uppercase, no spaces.
      sequence: "ATGGAACTGGGG…",        // (full 5738 nt)
      sequence_sha256: "…",             // for dedupe + change detection

      // Feature spans are 1-based inclusive (GenBank convention). The parser
      // converts to 0-based internally. Auto-populated from the .gb feature table;
      // user confirms/edits in the UI.
      features: {
        signal_peptide:   { start: 1,    end: 57   },
        constant_region:  { start: 58,   end: 1047, name: "IGHG1" },
        resistance:       { start: 2774, end: 3634, name: "AmpR"  },
        origin:           { start: 3779, end: 4452, name: "pUC"   },
        promoter:         { start: 4811, end: 5428, name: "CMV"   },
      },

      // Insertion site = where the linearized vector ends meet the insert.
      // 0-based; for the French IGH vector both equal 57 (overhangs adjacent).
      insertion_site: { oh5_end: 57, oh3_start: 57 },

      // The constant-region protein anchor used to confirm in-frame read-through.
      constant_anchor_aa: "ASTKGPSV",

      // Which overhang set this vector expects (see below). A vector and its
      // overhangs must be locus-consistent.
      overhang_set: "french_default",
      overhang_locus: "IGH",

      provenance: { added_by: "benjamin", added: "2026-05-22", source_file: "IgG1_heavy_chain_french_vector.gb" },
    },
    // … kappa and lambda light-chain vectors as separate entries …
  ],

  // ---- Overhang sets ----------------------------------------------------
  // oh_5 is shared between IGK and IGL; oh_3 differs (κ vs λ discriminator).
  overhang_sets: [
    {
      id: "french_default",
      display_name: "French vector overhangs",
      oh_5: {
        IGH: "GCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGT",
        IGK: "ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC",
        IGL: "ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC",   // == IGK on purpose
      },
      oh_3: {
        IGH: "GCTAGCACCAAGGGCCCATCGGTCTTCC",
        IGK: "CGTACGGTGGCTGCACCATCTGTCTTCATC",
        IGL: "GGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCGAGTGAGGAGCTTCAAGCCAACAAGGCC",
      },
      provenance: { added_by: "benjamin", added: "2026-05-22" },
    },
    // Other labs add their own sets here; projects pick which to use.
  ],

  // ---- Naming profiles (see docs/04_NAMING.md) --------------------------
  naming_profiles: [
    {
      id: "default",
      display_name: "Default (broad synonyms)",
      // Longest-match-first; matched case-insensitively as a delimited token.
      chain_synonyms: {
        heavy:  ["heavychain", "heavy", "hchain", "hc", "vh", "igh", "h"],
        kappa:  ["kappachain", "kappa", "igk", "vk", "k"],
        lambda: ["lambdachain", "lambda", "igl", "vl", "l"],
        light:  ["lightchain", "light", "lchain", "lc"],   // resolves to κ/λ by autodetect
      },
      // Fallback regex with named groups, tried if no synonym token is found.
      id_regex: "^(?<ab_id>.+?)[_\\-\\s]*(?<chain>heavy|light|hc|lc|kappa|lambda|[hkl])$",
      // Separators normalized before tokenizing.
      separators: ["_", "-", " ", ":"],
    },
  ],

  // ---- Alignment + threshold settings -----------------------------------
  alignment: {
    protein_matrix: "BLOSUM62",
    protein_gap_open: -11,
    protein_gap_extend: -1,
    nt_match: 2, nt_mismatch: -3, nt_gap_open: -5, nt_gap_extend: -1,
    backbone_identity_min: 0.97,   // below this, backbone match is "uncertain"
    overhang_max_mismatch: 1,      // tolerate 1 nt typo in an overhang, but report it
  },

  // ---- Optimization advisory checks (mirror dnachisel) ------------------
  optimization_checks: {
    gc_max_global: 0.56,
    gc_max_window: 0.64, gc_window: 60,
    rare_codon_min_frequency: 0.05, species: "h_sapiens",
    enabled: true,   // advisory only; never fails an order
  },
}
```

### 1.4 Sharing & enrichment rules
- **Additive import:** importing another lab member's `library.json5` adds new vectors /
  overhang sets; it never silently overwrites.
- **Dedupe by `sequence_sha256`** for vectors and by exact-sequence for overhang sets.
- On a hash collision with differing metadata, surface a conflict for the user to resolve
  (keep both with suffixed ids, or pick one). Never drop data silently.
- `provenance.added_by` lets the lab see who contributed each entry.

---

## 2. Project file — `<project>/project.json5`

Each project is a folder (shown in the sidebar, §`03_ARCHITECTURE.md`). The folder holds
`project.json5`, the loaded input files (copied in or referenced by path), and generated
reports.

```json5
{
  schema_version: 1,
  name: "round4_panel",
  created: "2026-06-17",
  description: "Round-4 wetlab panel, French IgG1 + κ/λ",

  naming_profile: "default",      // ref into library.naming_profiles
  overhang_set: "french_default", // ref into library.overhang_sets

  // Which library vector each chain class maps to for this project.
  vector_assignments: {
    heavy:  "french_igg1_heavy",
    kappa:  "french_igk_light",
    lambda: "french_igl_light",
  },

  // ---- Inputs (paths relative to the project folder, or absolute) -------
  ground_truth: {
    file: "round4_wetlab_panel.csv",
    format: "csv",                // csv | xlsx | fasta
    columns: { ab_id: "ab_id", heavy: "H_seq", light: "L_seq" }, // alias overrides
    // seq_type per chain auto-detected (nt vs AA); user can pin it.
  },

  order: {
    file: "IDT_ordered_sequences_correct.fasta",
    format: "fasta",              // fasta | xlsx
    has_overhangs: true,
  },

  sequencing: {
    mode: "full_plasmid",         // full_plasmid | partial_sanger | mixed
    inputs: [
      { file: "seq/HVA-0195-r3-d02_heavy.fasta", format: "fasta" },
      { file: "seq/junctions_plate1.ab1",        format: "ab1", mode: "partial_sanger" },
      // a whole folder may be globbed; one circular FASTA == one sample.
    ],
  },

  // ---- Results (written by the app; safe to delete & regenerate) --------
  results: {
    gate1: "results/order_qc.json",
    gate2: "results/sequencing_qc.json",
    report_html: "results/report.html",
  },
}
```

### Notes
- **Library is global; the other three inputs are per-project.** This mirrors how
  overhangs are reused across rounds while each panel/order/sequencing set is campaign-specific.
- All references between project and library are **by `id`**, so renaming a display name
  never breaks a project.
- The app must tolerate partially-filled projects (e.g. ground truth + order loaded but no
  sequencing yet → Gate 1 runnable, Gate 2 greyed out).
