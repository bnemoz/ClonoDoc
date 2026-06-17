# Verified reference facts & fixtures

Everything below was confirmed by direct computation on the user's real files
(`IgG1_heavy_chain_french_vector.gb`, `IDT_order.xlsx`, `IDT_order_correct.xlsx`,
`IDT_ordered_sequences_correct.fasta`, `optimize_clone.ipynb`). Use these as
**golden fixtures** for tests. Coordinates are stated in both 1-based (GenBank
convention) and 0-based (Rust slice convention) where it matters.

---

## 1. The French IgG1 heavy-chain vector

- File: `IgG1_heavy_chain_french_vector.gb`, **5738 bp, circular**.
- Reading frame opens at the leader ATG at position 1.

Key features parsed from the GenBank feature table:

| Feature | 1-based span | Notes |
|---|---|---|
| `sig_peptide` | 1..57 | Leader; 19 codons; translates `MELGLRWVFLVAILEGVQC` |
| 5′ overhang (`IgG1 seamless overhang 5' Heavy chain`, primer_bind) | 17..57 | exact homology arm |
| C-region (`IGHG1 C Region`, gene) | 58..1047 | starts `ASTKGPSVF…`, ends `…SLSLSPGK` |
| 3′ overhang (`IgG1 seamless overhang 3' Heavy chain`, primer_bind, complement) | 58..85 | exact homology arm |
| `Fc` | 352..1047 | hinge+CH2+CH3 |
| `AmpR` (CDS) | 2774..3634 | **resistance gene** |
| `AmpR promoter` | 2669..2773 | |
| pUC `ori` | 3779..4452 | |
| `CMV` promoter | 4811..5428 | expression |
| `CMV enhancer` | 4824..5203 | |

**Empty-vector ORF** (no insert; leader fused straight to constant region):
`MELGLRWVFLVAILEGVQC` + `ASTKGPSVF…SLSLSPGK`. This is the signature of an
empty/religated vector at Gate 2.

---

## 2. Overhang sets (from the notebook, French vector)

```
oh_5 = {
  'IGH': 'GCTGGGTTTTCCTTGTTGCTATTCTCGAGGGTGTCCAGTGT'   # 41 nt
  'IGK': 'ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC'        # 36 nt
  'IGL': 'ATCCTTTTTCTAGTAGCAACTGCAACCGGTGTACAC'        # 36 nt  (== IGK)
}
oh_3 = {
  'IGH': 'GCTAGCACCAAGGGCCCATCGGTCTTCC'                # 28 nt
  'IGK': 'CGTACGGTGGCTGCACCATCTGTCTTCATC'              # 30 nt
  'IGL': 'GGTCAGCCCAAGGCTGCCCCCTCGGTCACTCTGTTCCCGCCCTCGAGTGAGGAGCTTCAAGCCAACAAGGCC'  # 72 nt
}
```

**Critical detail for locus autodetect:**
- `oh_5['IGK'] == oh_5['IGL']` → the 5′ overhang is **shared** between light loci and
  therefore **cannot** distinguish κ from λ.
- `oh_3['IGK'] != oh_3['IGL']` → the 3′ overhang **does** distinguish them.
- ⇒ **Light-chain locus (κ vs λ) must be called from the 3′ end / constant region, never
  from the 5′ overhang.** Heavy is unambiguous from either end.

Overhang ↔ vector mapping (verified exact substring matches):
- `oh_5['IGH']` == vector 17..57 (0-based 16..57).
- `oh_3['IGH']` == vector 58..85 (0-based 57..85).
- The two IGH overhangs are **directly adjacent** in the empty vector (oh_5 ends at 57,
  oh_3 begins at 58). The insert's job is to land its V-region core between them.

---

## 3. The assembly model (verified)

The ordered fragment is `oh_5[locus] + optimized_core + oh_3[locus]`. The overhangs are
exact copies of the vector flanks, so in-silico assembly **merges** them with the vector
(they are not duplicated):

```
assembled = vector[0 .. oh5_end]  +  core  +  vector[oh3_start ..]
```

where `oh5_end` is the 0-based index just past the 5′ overhang in the vector and
`oh3_start` is the 0-based start of the 3′ overhang. For the French IGH vector,
`oh5_end == oh3_start == 57`.

### Golden PASS fixture (correct order)
- `HVA-0195-r3-d02_heavy` from `IDT_ordered_sequences_correct.fasta`.
- Ordered length 447 = 41 (oh5) + **378** (core, `378 % 3 == 0`) + 28 (oh3).
- Assembled ORF (leader→V→C) is stop-free until the **natural** stop at aa 475:
  `MELGLRWVFLVAILEGVQC·QVQLVESGGG…WGQGTTVTVSS·ASTKGPSVFPLAPSS…SLSLSPGK*`
- Reads through into the constant region (`ASTKGPSV…` present). **PASS.**

### Golden FAIL fixture (frameshift — Gate 1 catch)
- `HVA-0195-r3-d02_heavy` from the *first* `IDT_order.xlsx` (the "wrong" file).
- Core length **379** (`379 % 3 == 1`) — one stray trailing base.
- Assembled ORF: perfect V-domain, then `…WGQGTTVTVSS·G·*` — **premature stop at aa 146**,
  no constant-region read-through. **FAIL: junction frameshift / non-codon-boundary trim.**
- The V-domain *protein* is identical to the PASS fixture; only the frame at the 3′
  junction differs. This is exactly why the read-through check (translate full ORF, look
  for the stop) is the linchpin — an AA-only check of the V region alone would call this
  PASS incorrectly.

---

## 4. The optimization pipeline being verified (`optimize_clone.ipynb`)

This is *the user's* convention; others vary (see `docs/04_NAMING.md`). It is documented
here so Gate 1 knows what invariants to check, **not** so the tool reproduces it.

```
round4_wetlab_panel.csv            # columns: ab_id, H_seq, L_seq  (AA designed, or nt patient)
        │
        ├─ abverse.reverse_translate(AA → candidate nt)
        │
        ├─ abstar.run(...)         # annotate → locus (IGH/IGK/IGL), sequence_aa, productive
        │
        ├─ DnaOptimizationProblem  # dnachisel
        │     constraints:
        │        EnforceTranslation(sequence_aa)      ← THE INVARIANT Gate 1 verifies
        │        EnforceGCContent(maxi=0.56)
        │        AvoidRareCodons(min_frequency=0.05, species=h_sapiens)
        │        UniquifyAllKmers(8)
        │     objectives:
        │        CodonOptimize(h_sapiens, boost=10)
        │        EnforceGCContent(maxi=0.64, window=60)
        │        UniquifyAllKmers(6)
        │     (input trimmed to a multiple of 3 first: offset = len%3, drop trailing)
        │
        ├─ cloned = oh_5[locus] + optimized + oh_3[locus]
        │
        └─ write IDT_ordered_sequences_correct.fasta
```

The notebook's own QC (cells 18–23) reduces to two abstar fields: **`productive`**
(in-frame, no premature stop) and **`locus`** (correct chain class). Gate 1 reproduces
the *spirit* of these natively (no abstar dependency) plus the explicit
ground-truth-AA comparison.

**Gate-1 invariants (what to check), all passed by the correct file, failed by the wrong one:**
1. Overhangs present and **locus-matched** (κ overhang on a κ chain, etc.).
2. `translate(core) == ground_truth_AA` (the EnforceTranslation invariant, after the fact).
3. Assembled construct is **productive** (stop-free read-through into the constant region).
4. (Advisory) GC ≤ ceiling; no reintroduced rare codons.

---

## 5. Ground-truth input schema (from the notebook)

A table with one row per antibody:

| column | meaning | type |
|---|---|---|
| `ab_id` | antibody identifier (stem) | string |
| `H_seq` | heavy chain | **AA** (computer-designed) or **nt** (patient-isolated) |
| `L_seq` | light chain | AA or nt |

Comparison is **always at AA level** (codon-opt preserves AA, not nt). If a cell is nt,
translate it first; if AA, use directly. Auto-detect nt vs AA by alphabet
(`{A,C,G,T,U,N}` ⊂ alphabet and length % 3 == 0 ⇒ treat as nt candidate; otherwise AA).
Accept FASTA as an alternative to the table (record id → `ab_id` + chain token).

---

## 6. Sample / sequence naming observed

- Standard: `HVA-0195-r3-d02_heavy`, `HVA-0195-r3-d02_light` → `ab_id = HVA-0195-r3-d02`,
  chain token `heavy`/`light`.
- De novo design: `UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4__heavy`
  → contains a colon and a **double** underscore; `ab_id` is everything before the final
  chain token. The parser must survive this (see `docs/04_NAMING.md`).
- Light chains here carry only `light` (no κ/λ token); κ/λ is resolved by autodetect from
  the 3′ end. The κ example (`B10`) translates to a Vκ ending `…WTFGGGTKVEIK`.

---

## 7. Sequencing input modes

- **Plasmidsaurus full-plasmid**: one circular consensus FASTA **per sample**. Arbitrary
  rotation and strand. Both Gate-2 checks possible (backbone + insert + read-through).
- **Partial Sanger (other vendors)**: junction + insert only, AB1 or FASTA. Backbone not
  fully observed → backbone identity is *inferred from flanks only* and must be reported
  as such, not silently passed.
- **AB1 quality**: Plasmidsaurus AB1 is assembled from Oxford Nanopore passes; its Phred
  values are not meaningful. **Do not gate on quality.** Parse base calls (`PBAS`/`PBAS2`
  tag); treat quality (`PCON`) as advisory/off by default.
