# 04 — Sample-name parsing & chain pairing

Naming conventions vary between people and projects. The parser must extract, from an
arbitrary record id, two things:

1. `ab_id` — the antibody identifier (the grouping key for pairing).
2. `chain_class` — one of `heavy | kappa | lambda | light` (`light` = κ/λ to be resolved
   later by autodetect).

It must be **configurable** (per-project `naming_profile`) and **forgiving**, and when it
cannot decide confidently it must **ask the user** rather than guess.

---

## 1. Algorithm

```
parse(record_id, profile) -> { ab_id, chain_class, confidence, needs_confirmation }
```

1. **Normalize separators.** Replace every char in `profile.separators`
   (default `_ - <space> :`) with a single `\x1f` token boundary. Keep the original
   string too (for `ab_id` reconstruction).
2. **Tokenize** on the boundary. Lowercase a copy for matching.
3. **Find the chain token, longest-match-first.** Iterate `profile.chain_synonyms` with
   synonyms sorted by descending length so `heavychain` matches before `heavy` before
   `hc` before `h`. Prefer a token at the **end** of the id (chain suffix is the common
   convention); if multiple candidates, the rightmost wins.
   - Map the matched synonym to its class via the synonym→class table.
4. **Reconstruct `ab_id`** = the original id with the matched chain token (and any
   adjacent separators it consumed) removed. Trim leftover trailing separators.
5. **Fallback to `id_regex`** if step 3 found nothing: apply the profile's regex with
   named groups `ab_id` and `chain`; map `chain` through the synonym table.
6. **Confidence & confirmation.**
   - High: a multi-letter token matched (`heavy`, `kappa`, `lc`, …).
   - Low / `needs_confirmation = true`: only a **single-letter** token matched (`h`/`k`/`l`)
     — these are dangerous (they collide with arbitrary text), so flag for user review.
   - None: no token and no regex match → `needs_confirmation = true`, `chain_class` unknown.

Single-letter synonyms are included for completeness but **always** route to confirmation
unless the user has pinned a profile that trusts them.

---

## 2. Default synonym table

| class | synonyms (longest-first) |
|---|---|
| `heavy` | `heavychain`, `heavy`, `hchain`, `hc`, `vh`, `igh`, `h` |
| `kappa` | `kappachain`, `kappa`, `igk`, `vk`, `k` |
| `lambda` | `lambdachain`, `lambda`, `igl`, `vl`, `l` |
| `light` | `lightchain`, `light`, `lchain`, `lc` |

`light` is deliberately distinct from κ/λ: many panels (including the user's) label light
chains only as `light` and leave κ/λ to be **autodetected from the 3′ overhang / constant
region** (see `01_DESIGN.md` §1.5). When a record parses to `light`, the locus is resolved
downstream, not at name-parse time.

---

## 3. Worked examples (must pass as tests)

| record id | ab_id | chain_class | confidence |
|---|---|---|---|
| `HVA-0195-r3-d02_heavy` | `HVA-0195-r3-d02` | heavy | high |
| `HVA-0195-r3-d02_light` | `HVA-0195-r3-d02` | light | high |
| `UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4__heavy` | `UNREG:GTTCATTGTCATGCCG_d02_w74_esmfold_bb42m4` | heavy | high |
| `mab1_HC` | `mab1` | heavy | high |
| `Ab_007_lambda` | `Ab_007` | lambda | high |
| `clone3-heavychain` | `clone3` | heavy | high |
| `7G12_kappa` | `7G12` | kappa | high |
| `sample12_H` | `sample12` | heavy | **low → confirm** |
| `weird_name_no_token` | `weird_name_no_token` | unknown | **confirm** |

Note the double-underscore / colon case must not corrupt `ab_id`: only the **final**
chain token is stripped; internal `_d02_`, `:`, etc. are preserved.

---

## 4. Pairing & rollup
- Group parsed records by `ab_id`.
- Expected membership per antibody is configurable (default: exactly one heavy + one
  light). Use the project's `vector_assignments` keys to know which classes are expected.
- A `light` record satisfies the "one light" expectation regardless of whether it later
  resolves to κ or λ.
- **Missing** expected chain → `INCOMPLETE_PAIR`. **Extra** chain for an `ab_id` → review flag.
- The pairing is presented in the batch table grouped by antibody, each row expandable to
  its chains.

---

## 5. Manual override
- The UI exposes the parse result per record with an editable `ab_id` and a `chain_class`
  dropdown. Overrides are stored in the project file:
  ```json5
  name_overrides: { "weird_name_no_token": { ab_id: "WN-01", chain_class: "heavy" } }
  ```
- Overrides take precedence over the parser. They are per-project (different campaigns may
  reuse an id differently).

---

## 6. Robustness requirements
- Case-insensitive matching; never mutate the stored original id.
- Tolerant of mixed separators within one id.
- Deterministic: the same id + profile always yields the same parse.
- Never throw on a weird id — return `unknown` + `needs_confirmation` instead.
