# 03 — Architecture

Implementation guidance for the Rust + `egui` build. This is a recommended structure;
deviate where a cleaner design emerges, but preserve the **core/UI separation** (§3) so
the verification engine is testable headlessly.

---

## 1. Stack decision & rationale

| Concern | Choice | Why / why not the alternative |
|---|---|---|
| Language | **Rust** | Single static binary, no runtime. End users never read source, so readability isn't a constraint. |
| GUI | **`eframe`/`egui`** | One binary, **no webview, no installer**, cross-platform. Tauri reintroduces a webview + JS build chain; Go+Fyne has a weaker alignment ecosystem. |
| Alignment | **`rust-bio`** | Pairwise affine-gap DP, custom & BLOSUM62 scoring. Sequences here are short (≤~6 kb), so full/banded DP is fine. |
| XLSX read | **`calamine`** | Pure-Rust, no Excel/LibreOffice dependency. |
| XLSX write | **`rust_xlsxwriter`** | Pure-Rust report export. |
| CSV | **`csv`** | Standard. |
| GenBank | **`gb-io`** *(verify on crates.io)* or custom | If unavailable/limited, hand-roll (§4). |
| AB1/ABIF | **custom** | No mature crate; format is a simple tagged binary directory (§5). |
| Config | **`json5`** + **`serde`** | Comments in the lab-shared file. |
| Paths | **`directories`** | Per-OS config/data dirs. |
| Dialogs | **`rfd`** | Native open/save. |
| Hashing | **`sha2`** | Vector dedupe. |
| Errors | **`anyhow`** (app) / **`thiserror`** (lib) | Ergonomic. |
| Parallel batch | **`rayon`** | Per-construct verification is embarrassingly parallel. |

---

## 2. Crate layout (workspace)

```
clonodoc/
├── Cargo.toml                  # [workspace]
├── crates/
│   ├── core/                   # headless verification engine (no UI deps) — lib
│   │   ├── seqio/              # FASTA, XLSX, CSV, GenBank, AB1 readers
│   │   ├── model/              # Vector, OverhangSet, Project, Library, Records (serde)
│   │   ├── naming/             # name parser + pairing (docs/04)
│   │   ├── align/              # rust-bio wrappers: nt, protein, circular-normalize
│   │   ├── assemble/           # in-silico assembly (docs/01 §3.2)
│   │   ├── gate1/              # order QC
│   │   ├── gate2/              # sequencing QC
│   │   └── report/             # verdict structs, CSV/XLSX/HTML export
│   └── app/                    # eframe binary (depends on core)
│       ├── main.rs
│       └── ui/                 # sidebar, panels, tables, detail/alignment view
└── tests/                      # integration tests using reference/ fixtures
```

**Rule:** `core` has **zero** UI dependencies and is fully exercised by `tests/` against
the fixtures in `reference/verified_facts.md`. The CLI-less engine must produce every
verdict without egui. This keeps the science testable and the GUI thin.

---

## 3. Data flow

```
        ┌─────────────── Library (global json5) ───────────────┐
        │  vectors · overhang_sets · naming_profiles · settings │
        └───────────────────────────────────────────────────────┘
                              │ (by id)
        ┌─────────────────────▼─────────────────────┐
        │            Project (per-folder json5)       │
        │  ground_truth · order · sequencing · refs   │
        └───────┬───────────────────────────┬─────────┘
                │                            │
        ┌───────▼────────┐          ┌────────▼─────────┐
        │  Gate 1 engine │          │  Gate 2 engine    │
        │  (pre-cloning) │          │  (post-cloning)   │
        └───────┬────────┘          └────────┬─────────┘
                │   assemble + align (shared core::align/assemble)
                └──────────────┬─────────────┘
                       ┌────────▼────────┐
                       │  Verdicts +      │
                       │  report export   │
                       └─────────────────┘
```

Shared between gates: `seqio`, `model`, `naming`, `align`, `assemble`. Gate 1 and Gate 2
are thin orchestration layers over those.

---

## 4. Minimal GenBank parser (if not using `gb-io`)

Enough to populate a `Vector`:
- `LOCUS` line → name, length, `circular|linear`.
- `FEATURES` block → for each feature: type (col 6+), location string, qualifiers
  (`/standard_name=`, `/gene=`, `/note=`). Parse simple `start..end`,
  `complement(start..end)`, and `join(...)` (take min start / max end for our coarse
  spans; we only need approximate feature windows). Ignore `<`/`>` partial markers.
- `ORIGIN` → concatenate the base lines, strip digits/whitespace, uppercase.
- Map feature `standard_name`/`gene` to roles by keyword:
  `signal peptide|sig_peptide → signal_peptide`; `C Region|IGHG|IGHA|constant|Fc →
  constant_region`; `AmpR|KanR|CmR|resistance → resistance`; `ori|origin → origin`;
  `CMV|promoter → promoter`. Show the mapping in the UI for confirmation (annotations are
  Geneious-derived and occasionally messy).

---

## 5. AB1 / ABIF parser (custom)

ABIF = header + a directory of tagged records.
- Bytes 0..4: `"ABIF"` magic; next 2 bytes version.
- A directory entry is 28 bytes: tag name (4 ASCII), tag number (i32), element type
  (i16), element size (i16), num elements (i32), data size (i32), data offset (i32),
  handle (i32). Read big-endian.
- The directory's own location is given by the entry at offset 26.
- Extract **`PBAS` tag number 2** (or `1`) → ASCII base calls (the sequence).
- Extract `PCON` (quality) **only as advisory**; default off.
- Return `{ bases: String, quality: Option<Vec<u8>> }`. Do not fail if `PCON` is absent.

Keep this ~150 lines and unit-test against a couple of real AB1s during implementation.

---

## 6. UI structure (egui)

```
┌──────────────┬───────────────────────────────────────────────┐
│  SIDEBAR     │  MAIN PANEL                                    │
│              │                                                │
│  [Library]   │  Tabs: [Setup] [Gate 1: Order] [Gate 2: Seq]   │
│  ──────────  │                                                │
│  Projects:   │  Setup tab:                                    │
│   ▸ round4   │    • 4 load buttons (see below)                │
│   ▸ round5   │    • parsed-record table w/ confirm dropdowns  │
│   ▸ …        │                                                │
│  [+ New]     │  Gate tabs:                                    │
│              │    • batch verdict table (sortable, colored)   │
│              │    • row click → detail pane (alignment view)  │
│              │    • [Run] and [Export] buttons                │
└──────────────┴───────────────────────────────────────────────┘
```

**The four load buttons** (user's spec), all on the Setup tab of the selected project:
1. **Load Library** → vectors + overhang sets (opens global settings; shared file).
2. **Load Ground Truth** → project VDJ panel (CSV/XLSX/FASTA; AA or nt).
3. **Load IDT Order** → ordered sequences with overhangs (FASTA/XLSX) → enables Gate 1.
4. **Load Sequencing Results** → Plasmidsaurus FASTA / Sanger AB1 → enables Gate 2.

- **Sidebar = project explorer.** Projects are folders; `[+ New]` creates one. Library
  sits above projects (global). Selecting a project loads its `project.json5`.
- **Confirmations** (user-in-the-loop): autodetected vector, low-confidence name parses,
  and nt-vs-AA calls each surface an inline confirm control before a gate runs.
- **Detail/alignment view:** monospace, V region vs expected with mismatches highlighted;
  for frameshift, show where the ORF diverges and the premature stop; for nt silent
  variants, a separate nt track over the insert window.

---

## 7. Threading
- Verification runs off the UI thread. Use `rayon` to verify constructs in parallel inside
  a gate; send progress + results back to the UI via a channel (`std::sync::mpsc` or
  `crossbeam`). egui repaints on message receipt.
- File parsing for large batches likewise off-thread; show a progress bar.

---

## 8. Testing hooks
- `core` exposes pure functions: `parse_name`, `detect_locus`, `assemble`, `translate`,
  `align_protein`, `run_gate1`, `run_gate2`. Each is unit-tested.
- Integration tests load the fixtures (`reference/`) and assert the golden PASS/FAIL
  verdicts (especially the 378-vs-379 frameshift case). See `05_IMPLEMENTATION_PLAN.md`.
