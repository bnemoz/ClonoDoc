# 01 — Design: verification logic, algorithms, failure taxonomy

This document specifies *what* the tool computes and *why*. Implementation structure is
in `03_ARCHITECTURE.md`; configuration shapes are in `02_CONFIG_SCHEMA.md`.

---

## 0. Foundational principle

Codon optimization means the cloned/ordered **nucleotide** sequence of the variable
region intentionally diverges from the source while preserving the **protein**. So:

- **Variable region** → verify at the **amino-acid** level.
- **Backbone + overhangs** (not codon-changed) → verify at the **nucleotide** level.
- **Frame / read-through** → verify by translating the whole assembled ORF and locating
  the stop codon.

Three failure modes therefore produce three distinct, separable signatures:

| Failure | AA of V region | Frame at 3′ junction | Constant region |
|---|---|---|---|
| Point mutation | ≥1 mismatch | intact | reads through |
| Junction indel / frameshift | may look perfect, then diverges | broken | premature stop, no read-through |
| Empty / religated vector | absent (no V at all) | n/a | leader fused straight to constant |

All three are caught by one operation: translate the full ORF, globally align the
protein to the expected protein, inspect where it breaks.

---

## 1. Shared primitives

### 1.1 Sequence ingestion
- **FASTA**: standard; tolerate wrapped lines, lowercase, trailing whitespace, `*`/`-`.
- **XLSX**: via `calamine`; expect a header row, columns `sequence_id|seq_id|id` and
  `sequence|seq`. Be liberal in header matching (case-insensitive, trimmed).
- **CSV**: `csv` crate; same liberal header matching. Ground-truth panel columns
  `ab_id`, `H_seq`, `L_seq` (configurable aliases).
- **GenBank** (`.gb`/`.gbk`): parse `LOCUS` (length, topology), `FEATURES` (type, span,
  `/standard_name`, `/gene`, `/note`), and `ORIGIN` (sequence). See `03_ARCHITECTURE.md`
  for the minimal parser spec.
- **AB1 / ABIF**: hand-rolled. Read the tag directory; extract base calls from `PBAS2`
  (fallback `PBAS1`). Quality (`PCON`) parsed but **advisory only**, default off.

### 1.2 nt vs AA autodetection
A string is treated as **nt** iff its alphabet ⊆ `{A,C,G,T,U,N}` (case-insensitive) and
its length is a multiple of 3 OR ≥ 90% of characters are in the DNA alphabet. Otherwise
**AA**. Always expose the decision in the UI for confirmation; never silently mis-call.

### 1.3 Translation
Standard genetic code, frame 0 unless an explicit ORF start is located. Stop = `*`.
Provide a "translate until first stop" and a "translate full length" mode.

### 1.4 Alignment
- **Nucleotide**: `rust-bio` pairwise, semiglobal or global as appropriate, match `+2` /
  mismatch `-3` / gap-open `-5` / gap-extend `-1` (tunable in config). Used for backbone
  identity, overhang matching, and (when an IDT nt reference exists) silent-SNP detection.
- **Protein**: `rust-bio` pairwise, **BLOSUM62**, gap-open `-11` / gap-extend `-1`
  (tunable). Used for the V-region and full-ORF comparisons.
- **Circular topology**: before any whole-plasmid alignment, **normalize** (§3.1).

### 1.5 Locus / chain-class detection
Two independent signals, cross-checked:
1. **Overhang match**: which `(oh_5, oh_3)` pair the sequence carries. Remember
   `oh_5` is shared between κ and λ → **use `oh_3` (or the constant region) to split κ/λ.**
2. **Constant-region motif** (protein): heavy ends V at `…WGxGT[LT]VTVSS` then `ASTKGP…`;
   κ ends `…FG.GTKVEIK`; λ ends `…TVL`. Used as a backbone-independent confirmation and
   for partial reads lacking a full overhang.

The pair must agree; disagreement is itself a reportable flag.

---

## 2. Gate 1 — In-silico order QC (pre-cloning)

**Inputs:** project ground truth (panel), the IDT order file (with overhangs), the
selected vector library + overhang set. **No sequencing data.**

For each ordered record:

1. **Parse name** → `ab_id` + chain token → locus class (`docs/04_NAMING.md`).
2. **Strip overhangs.** Find `oh_5[locus]` as a prefix and `oh_3[locus]` as a suffix
   (allow ≤1 mismatch to tolerate a typo, but report any mismatch). The remainder is the
   `core`. Failure here ⇒ **OVERHANG_MISSING / OVERHANG_WRONG_LOCUS**.
   - Check the overhang's locus matches the name's locus: a `_heavy` carrying κ overhangs
     ⇒ **LOCUS_MISMATCH**.
3. **Translate `core`** (frame 0) → `core_aa`.
4. **Compare to ground truth.** Look up `ab_id` + chain in the panel; obtain
   `truth_aa` (translate if the panel cell is nt). Global protein align `core_aa` vs
   `truth_aa`.
   - Identical ⇒ translation faithful (EnforceTranslation invariant held).
   - Mismatch ⇒ **TRANSLATION_DRIFT** (report positions; this means the optimizer, the
     reverse-translation, or a manual edit changed the protein).
   - No matching panel entry ⇒ **NO_GROUND_TRUTH** (advisory, not fail).
5. **Assemble & test productivity.** Build `assembled = vector[..oh5_end] + core +
   vector[oh3_start..]` (§3.2), translate the full ORF from the leader ATG, and confirm:
   - no premature stop before the natural constant-region stop, AND
   - the constant-region anchor (e.g. `ASTKGPSV` for IgG1 H, locus-appropriate for L) is
     present in-frame.
   - Premature stop / no read-through ⇒ **FRAMESHIFT_AT_JUNCTION** (this is the 379-bug
     fixture; `core_len % 3 != 0` is a fast pre-check but the authoritative test is the
     ORF translation).
6. **Advisory optimization checks** (do not fail the order, just warn): overall GC ≤
   configured ceiling; windowed GC; presence of any disallowed rare codons. These mirror
   the dnachisel constraints so a hand-edited order is flagged.

**Gate-1 verdicts:** `PASS`, `TRANSLATION_DRIFT`, `FRAMESHIFT_AT_JUNCTION`,
`OVERHANG_MISSING`, `OVERHANG_WRONG_LOCUS`, `LOCUS_MISMATCH`, `NO_GROUND_TRUTH`,
plus advisory `GC_WARNING` / `RARE_CODON_WARNING`.

---

## 3. Gate 2 — Sequencing verification (post-cloning)

**Inputs:** sequencing reads (full-plasmid FASTA or partial Sanger AB1/FASTA) + the
expected reference (IDT order file and/or ground-truth panel) + the vector library.

### 3.0 Reference reconstruction
For each expected construct, build the **expected full plasmid** in silico (§3.2) and the
**expected protein** = `leader_aa + V_aa + constant_aa`, where `V_aa` comes from the
ground truth (always available at AA level) and, if the IDT order is loaded, the
**expected insert nt** (`core`) is also available for silent-SNP detection.

### 3.1 Topology normalization (full-plasmid mode)
A Plasmidsaurus consensus is circular, arbitrary rotation, arbitrary strand.
1. **Strand**: detect by matching a long invariant backbone landmark (e.g. AmpR or the
   constant region) on both strands; pick the strand with the better match; if the
   reverse-complement wins, RC the read.
2. **Rotation**: locate a canonical anchor (the leader ATG region, or the start of the
   constant region — both are codon-opt-invariant) and rotate the linear representation
   so the construct ORF starts near the beginning. Equivalent and simpler for alignment:
   **double the reference** (`ref+ref`) and align the read against it, or double the read.
3. Only after normalization do backbone and ORF alignments run.

### 3.2 In-silico assembly (shared by Gate 1 and Gate 2 reference build)
```
oh5_end   = index in vector just past oh_5[locus]      (French IGH: 57)
oh3_start = index in vector at the start of oh_3[locus] (French IGH: 57)
assembled = vector[0 .. oh5_end] ++ core ++ vector[oh3_start ..]
```
Store `assembled` (nt) and its translated ORF as the expected references.
The insert window within `assembled` is `[oh5_end, oh5_end + core.len())`.

### 3.3 Check #1 — Backbone identity (full-plasmid only)
nt-align the read (minus the insert window) against each candidate vector backbone.
Best match yields:
- **isotype** (from the matched vector's constant-region feature: IgG1 / IgG4 / IgA / κ / λ…),
- **resistance** (from the matched vector's resistance feature: AmpR / KanR …),
- **locus class** of the backbone,
- a confidence (percent identity over the backbone).
Cross-check against the name's locus (heavy token must hit a heavy backbone; light token
may hit κ **or** λ). Mismatch ⇒ **WRONG_VECTOR**. User confirms the autodetected vector.

For **partial Sanger**, the backbone is not in the read: identity is inferred only from
the immediate flanks/overhangs and reported as `backbone_observed = flanks_only` —
never a full PASS on backbone.

### 3.4 Check #2 — Insert + frame (both modes)
1. Locate the ORF in the read: anchor on the leader (codon-opt-invariant nt) for
   full-plasmid; for partial Sanger, anchor on the 5′ overhang / first in-frame codon.
2. Translate the ORF (full length) → `observed_aa`.
3. **Global protein-align** `observed_aa` vs `expected_aa` (leader+V+constant).
   - Clean match through the constant region ⇒ insert correct **and** Fc in-frame
     (Check #1 corroborated at protein level) ⇒ contributes to `PASS`.
   - Isolated substitutions in the V span ⇒ **POINT_MUTATION** (report position, wt→mut,
     and whether it lies in a CDR vs framework if IMGT/Kabat numbering is enabled;
     otherwise just linear position).
   - Divergence + premature stop + no constant read-through ⇒ **JUNCTION_FRAMESHIFT**;
     localize the indel by the alignment gap.
   - No V at all, leader→constant directly ⇒ **EMPTY_VECTOR**.
4. **Silent-SNP layer (optional, when IDT nt loaded):** nt-align the read's insert window
   against the expected `core`. nt mismatches with **no** AA change ⇒ **SILENT_VARIANT**
   (reported, not failed by default — these are synthesis SNPs that don't change protein).
5. **Sample-swap detection:** if `observed_aa` mismatches the *expected* antibody but
   matches a *different* panel member better, flag **WRONG_INSERT_SWAP** with the
   suspected correct identity.

### 3.5 Read adequacy
If a read is too short to cover the insert, has a coverage gap across the junction, or is
dominated by ambiguous bases, return **INSUFFICIENT_READ** rather than a false PASS/FAIL.

**Gate-2 verdicts (per chain):** `PASS`, `POINT_MUTATION`, `SILENT_VARIANT`,
`JUNCTION_FRAMESHIFT`, `EMPTY_VECTOR`, `WRONG_VECTOR`, `WRONG_INSERT_SWAP`,
`INSUFFICIENT_READ`. Severity ordering (worst wins when multiple apply):
`WRONG_VECTOR > EMPTY_VECTOR > JUNCTION_FRAMESHIFT > WRONG_INSERT_SWAP > POINT_MUTATION >
SILENT_VARIANT > PASS`.

---

## 4. Pairing and per-antibody rollup
- Group chains by `ab_id` (parser output). An antibody **PASSES** only if **all** its
  expected chains (typically one heavy + one light) individually pass.
- Missing expected chain ⇒ **INCOMPLETE_PAIR**.
- Extra/unexpected chain for an `ab_id` ⇒ flag for review.
- Batch summary: counts by verdict, a per-antibody table, and a drill-down per chain
  showing the protein alignment (and nt alignment for the insert window).

---

## 5. Output / reporting
- **In-app**: sortable batch table (antibody × chain × verdict × confidence), color-coded;
  click a row → detail pane with the alignment(s), the located mutation/indel, and the
  autodetected vector to confirm.
- **Export**: CSV/XLSX summary (calamine/`rust_xlsxwriter`), plus an optional
  self-contained HTML report per project. PDF is out of scope for v1.
- Every verdict carries a human-readable reason string (e.g. *"V-region identical to
  ground truth; +1 frameshift at 3′ junction → premature stop at codon 146; no
  constant-region read-through"*).

---

## 6. Explicit non-goals (v1)
- No de novo antibody numbering required for a PASS/FAIL (IMGT/Kabat CDR localization is
  a *nice-to-have* enhancement to mutation reporting, behind a flag).
- No dependency on Python/abstar/abutils/dnachisel — the relevant invariants are
  re-implemented natively.
- No trace (chromatogram) rendering for AB1 in v1; base calls only.
- Quality-score gating is disabled (Plasmidsaurus quality is not meaningful).
