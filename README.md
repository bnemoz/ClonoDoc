# Antibody Cloning Verifier (`abclone-verify`)

A small, portable, cross-platform (Windows / Linux / macOS) desktop application that
verifies antibody cloning at two independent stages:

1. **Pre-cloning (in-silico order QC)** — Before you send anything to IDT, confirm that
   the codon-optimized, overhang-appended order faithfully encodes your ground-truth
   antibody and will assemble in-frame into the chosen vector.
2. **Post-cloning (sequencing verification)** — After transformation, miniprep, and
   whole-plasmid or Sanger sequencing, confirm that the recovered construct is correct:
   right backbone (isotype, resistance, locus consistent with the sample name), right
   insert (no mutations, no frameshifts, no indels), in-frame read-through into the
   constant region.

It is the standalone replacement for the manual Geneious workflow, designed to run
**batch** (e.g. 20 chains for 10 antibodies arriving from sequencing at once) and to be
**lab-shareable**: the vector/overhang library lives in a single human-readable file that
any lab member can enrich.

---

## Why this tool exists (the core insight)

Inserts are **codon-optimized**, so the ordered/cloned nucleotide sequence deliberately
*does not* match the original source at the nucleotide level — only the amino-acid
sequence is preserved. We confirmed this empirically on real data: the longest exact
nucleotide match between a codon-optimized insert and its source vector region was 6 bp
(pure coincidence). **Therefore the variable region must be verified at the protein
level.** The vector backbone and the cloning overhangs, by contrast, are *not*
codon-changed and are verified at the nucleotide level.

This single fact drives the entire verification strategy (see `docs/01_DESIGN.md`).

---

## The two QC gates at a glance

| | Gate 1 — Order QC (pre-cloning) | Gate 2 — Sequencing QC (post-cloning) |
|---|---|---|
| **Question** | "Did I order the right molecule?" | "Did the cloning work?" |
| **Inputs** | Ground truth + IDT order file | IDT order (or ground truth) + sequencing reads |
| **Failure class caught** | Optimization bug, wrong/missing overhang, off-by-one frame error, locus mismatch | Point mutation, junction frameshift/indel, empty/religated vector, wrong vector, sample swap |
| **Cost of skipping** | Wasted synthesis + weeks of cloning on a dead design | False confidence in a bad clone |
| **Runs** | Once, before ordering | Once, after sequencing returns |

A real example of a Gate-1 catch is included as a regression fixture: an order whose
insert core was 379 nt instead of 378 (one stray trailing base) — the V-domain
translated perfectly but the construct frameshifted at the 3′ junction and never read
into the constant region. Gate 1 flags this in milliseconds, for free, before the order
is placed. See `reference/verified_facts.md`.

---

## Technology stack (and why)

- **Language: Rust.** No runtime, no interpreter; compiles to a single self-contained
  native binary per OS. The user explicitly does not need source-readability for end
  users, which frees the choice toward the most robust portable option.
- **GUI: `egui` / `eframe`.** Immediate-mode GUI that produces one statically-linked
  binary with **no webview and no installer** — drop the executable in a shared folder
  and it runs. Tables, forms, a project sidebar, and inline alignment views are all
  well within egui's wheelhouse.
- **Alignment: `rust-bio`** (pairwise, affine gaps, custom/BLOSUM62 scoring).
- **Spreadsheets: `calamine`** (read `.xlsx`) + the `csv` crate.
- **GenBank parsing: `gb-io`** if available on crates.io at build time (verify), else a
  small hand-rolled feature-table parser (the format is simple; see `docs/03_ARCHITECTURE.md`).
- **AB1 / ABIF parsing: hand-rolled** (~150 lines; no mature crate exists, and the
  format is a straightforward tagged binary directory). **Quality scores are ignored by
  default** — Plasmidsaurus Nanopore-derived AB1 quality is not meaningful (see DESIGN).
- **Config: `json5`** so the lab-shared library file can carry inline comments.
- **File dialogs: `rfd`** (native open/save dialogs on all three OSes).

Rationale for *not* choosing alternatives (Tauri/web UI, Go+Fyne, Python+packaging) is
in `docs/03_ARCHITECTURE.md`.

---

## Repository layout (target)

```
abclone-verify/
├── README.md                     ← this file
├── docs/
│   ├── 01_DESIGN.md              ← biology, algorithms, failure taxonomy (the WHAT & WHY)
│   ├── 02_CONFIG_SCHEMA.md       ← annotated library + project config (json5)
│   ├── 03_ARCHITECTURE.md        ← Rust crates, modules, data flow, UI, threading
│   ├── 04_NAMING.md              ← configurable sample-name parser + chain pairing
│   └── 05_IMPLEMENTATION_PLAN.md ← phased build for Claude Code + test fixtures
├── reference/
│   ├── verified_facts.md         ← empirically confirmed data used as fixtures
│   └── example_library.json5     ← pre-populated French IgG1 vector + overhang sets
└── (src/, Cargo.toml, tests/ …   ← created by Claude Code during implementation)
```

## Build (target behavior, for Claude Code to realize)

```bash
cargo build --release            # → target/release/abclone-verify[.exe]
```

Single binary, no external runtime. Cross-compilation targets: `x86_64-pc-windows-gnu`,
`x86_64-unknown-linux-gnu`, `aarch64-apple-darwin` / `x86_64-apple-darwin`.

---

## Reading order for the implementer

1. `reference/verified_facts.md` — ground-truth data and fixtures (start here).
2. `docs/01_DESIGN.md` — the verification logic.
3. `docs/02_CONFIG_SCHEMA.md` + `reference/example_library.json5` — the data model.
4. `docs/04_NAMING.md` — name parsing/pairing.
5. `docs/03_ARCHITECTURE.md` — how to structure the Rust code and UI.
6. `docs/05_IMPLEMENTATION_PLAN.md` — phased plan + acceptance tests.
