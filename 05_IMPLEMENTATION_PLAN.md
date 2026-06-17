# 05 — Implementation plan (for Claude Code)

Phased so each milestone is independently testable and the science (in `core`) is proven
before any UI is built. Build `core` first, headless, against the fixtures.

---

## Phase 0 — Scaffold
- `cargo new` workspace; `crates/core` (lib) + `crates/app` (eframe bin).
- Add dependencies (`03_ARCHITECTURE.md` §1). Verify `gb-io` availability; if absent,
  plan the custom parser (§4 there).
- CI: `cargo test`, `cargo clippy`, `cargo fmt --check`, release builds for the three OS
  targets.
- **Done when:** empty `egui` window opens and `cargo test` runs (0 tests).

## Phase 1 — `seqio` (parsers)
- FASTA, CSV, XLSX (`calamine`) readers with liberal header matching.
- GenBank parser → `Vector` (features + sequence + topology).
- ABIF/AB1 parser → base calls (+ advisory quality).
- nt-vs-AA autodetector.
- **Tests:** parse the real files; assert vector length 5738, topology circular, the
  feature spans in `reference/verified_facts.md` §1, and that the order FASTA yields 42
  records with the expected ids.

## Phase 2 — `model` + `naming`
- serde structs for `Library`, `Vector`, `OverhangSet`, `Project`, records; json5 load/save.
- Name parser + pairing (`04_NAMING.md`), including overrides.
- **Tests:** every row of `04_NAMING.md` §3 (incl. the `UNREG:…__heavy` and single-letter
  confirm cases); round-trip a library + project through json5.

## Phase 3 — `align` + `assemble`
- rust-bio wrappers: nt (configurable scoring) and protein (BLOSUM62).
- Circular normalization (strand + rotation / reference-doubling).
- In-silico assembly (`01_DESIGN.md` §3.2); insert-window computation.
- Locus detection from overhang pair (κ/λ via 3′) and from constant-region motif.
- **Tests (golden, from `reference/verified_facts.md`):**
  - Assemble the **correct** `HVA-0195-r3-d02_heavy`: ORF stop-free until aa 475,
    `ASTKGPSV` present in-frame ⇒ productive.
  - Assemble the **wrong** (379 nt) core: premature stop at aa 146, no read-through.
  - `oh_5[IGK] == oh_5[IGL]`, `oh_3[IGK] != oh_3[IGL]`; κ detected from 3′ only.
  - `oh_5[IGH]` maps to vector 17..57; `oh_3[IGH]` to 58..85.

## Phase 4 — `gate1` (order QC)
- Implement the six steps of `01_DESIGN.md` §2; emit per-record verdicts + reasons.
- **Acceptance:**
  - Correct order file → all `PASS` (given ground truth; `NO_GROUND_TRUTH` advisory if the
    panel isn't loaded).
  - First/wrong xlsx core (379) → `FRAMESHIFT_AT_JUNCTION` with reason citing the
    premature stop and no read-through.
  - A heavy record given κ overhangs → `OVERHANG_WRONG_LOCUS` / `LOCUS_MISMATCH`.
  - A core whose translation differs from ground truth → `TRANSLATION_DRIFT` with positions.

## Phase 5 — `gate2` (sequencing QC)
- Topology normalize → backbone identity (full-plasmid) → ORF translate → protein align →
  verdict taxonomy (`01_DESIGN.md` §3). Silent-SNP layer when IDT nt loaded. Sample-swap
  detection across the panel. Partial-Sanger path reports `backbone_observed=flanks_only`.
- Per-antibody rollup + batch summary.
- **Acceptance (construct synthetic reads from the assembled references):**
  - Perfect read → `PASS`.
  - Single-codon substitution in V → `POINT_MUTATION` at the right position.
  - 1-nt deletion at the 3′ junction → `JUNCTION_FRAMESHIFT`, no read-through.
  - Empty vector (leader→constant, no insert) → `EMPTY_VECTOR`.
  - Heavy insert in a κ backbone → `WRONG_VECTOR`.
  - Read matching a *different* panel member → `WRONG_INSERT_SWAP`.
  - Synonymous nt change vs IDT, identical AA → `SILENT_VARIANT` (not a fail).
  - Reverse-complement / rotated input → same verdict as forward (normalization works).

## Phase 6 — `report`
- Verdict structs → CSV/XLSX (`rust_xlsxwriter`) + self-contained HTML.
- Human-readable reason strings on every verdict.

## Phase 7 — `app` (egui)
- Sidebar project explorer; global Library settings panel; the four load buttons;
  Setup / Gate 1 / Gate 2 tabs; sortable colored verdict table; detail/alignment pane;
  confirm controls (vector, name, nt/AA); off-thread runs with progress (`rayon` + channel).
- **Acceptance:** end-to-end on the real files — create a project, load library + ground
  truth + correct order → Gate 1 all PASS; swap in the 379 order → the offending row flags
  FRAMESHIFT; (with synthetic or real reads) Gate 2 batch table populates and exports.

## Phase 8 — Packaging
- Release binaries for Windows/Linux/macOS; bundle a starter `library.json5`
  (`reference/example_library.json5`) on first run; document the shared-library path option.

---

## Cross-cutting acceptance: the two golden fixtures must never regress
1. **PASS:** correct `HVA-0195-r3-d02_heavy` assembles to a stop-free ORF reading through
   into `ASTKGPSV…SLSLSPGK`.
2. **FAIL:** the 379-nt core frameshifts (`…WGQGTTVTVSS·G·*`, stop at aa 146, no
   read-through) and is reported as a junction frameshift — **with the V-domain protein
   still identical to the PASS case**, proving the read-through check (not a V-only AA
   check) is what catches it.

Keep both wired as integration tests from Phase 3 onward.

---

## Notes for the implementer
- Honor the **core/UI split**: no `egui` types in `crates/core`.
- Coordinates: GenBank is 1-based inclusive; convert once at parse time and work 0-based
  internally. Document the convention at every boundary.
- Never gate verdicts on AB1 quality.
- All thresholds/scores come from the library `alignment` block — no magic numbers in code.
- Be liberal in input parsing, strict and explicit in verdicts; every non-PASS carries a
  reason a bench scientist can act on without reading code.
