# Handoff: Cloning Workbench — Desktop GUI (Direction A, "Bench")

## Overview
A cross-platform desktop application for molecular-cloning verification. Scientists use it to (1) curate a reusable library of reference parts, (2) load the files for a specific cloning project, (3) **verify an in-silico cloning design before ordering inserts**, and (4) **verify the finished clone against sequencing data after the wetlab work**. This document specifies the redesigned GUI so it can be rebuilt in the app's real codebase.

Single persistent window, no OS chrome. Dark theme, spacious/calm density, designed for a desktop viewport of **~1400px wide and up** (degrades acceptably down to ~1100px). Left sidebar (file navigation) is always present; the four workflows are top tabs in the main pane.

---

## About the Design Files
The files in this bundle are **design references created in HTML** — a working prototype showing the intended look, layout, and behavior. **They are not production code to copy.** The task is to **recreate this design in the target codebase's existing environment** (the framework the app already uses — e.g. Tauri/Electron + React/Svelte, Qt/QML, etc.), using its established components, state, and styling patterns. If no UI framework is established yet, pick the one most appropriate for a cross-platform desktop scientific tool and implement the design there.

The prototype is built as a "Design Component" (`.dc.html`) and needs the bundled `support.js` runtime + an internet connection (Google Fonts) to render. Open `Cloning Workbench - A.dc.html` in a browser to view it. The tabs are clickable — navigate all four screens.

## Fidelity
**High-fidelity (hifi).** Colors, typography, spacing, and component treatments are final and intended to be matched closely. The sample data (sequences, plasmid names, check results) is illustrative — wire it to real data. Exact values are in **Design Tokens** below.

---

## Design Tokens

### Color — surfaces & text (dark theme)
| Token | Hex | Use |
|---|---|---|
| `bg/app` | `#0E1217` | App background, deepest panels, sidebar search field |
| `bg/panel` | `#11161c` | Sidebar, cards, list rows |
| `bg/panel-alt` | `#10151b` | Drop zone / inset panels |
| `surface/chip` | `#161d24` | Inactive filter chips, format badges |
| `surface/muted` | `#1a222b` / `#1f2730` | Annotation bars, backbone segments |
| `border/strong` | `#2a333e` | Ghost-button borders, scrollbar thumb |
| `border/card` | `#1f272f` | Card & list-row borders |
| `border/soft` | `#1d242c` / `#232b34` | Sidebar dividers, search field border |
| `divider/hair` | `#161d24` | Header-strip bottom borders |
| `text/primary` | `#E7EDF3` | Headings, primary values |
| `text/secondary` | `#cdd6df` | Ghost button labels, active tree items |
| `text/body` | `#9aa7b4` / `#aeb9c4` | Body copy, secondary list text |
| `text/muted` | `#7e8c99` | Captions, detail metadata |
| `text/faint` | `#5e6b78` / `#566270` | Section labels, hints, placeholders |

### Color — brand & semantic (Scripps Research palette)
| Token | Hex | Use |
|---|---|---|
| `accent` (Nitrogen Sky) | `#5BC2D2` | **Primary accent**: active tab, primary buttons, selection, insert/EGFP highlight, junction divider. Themeable (see Tweaks). |
| `gold` (Infinity Gold) | `#FFC951` | Secondary/attention: AmpR feature, warnings, wetlab "review" verdict + button |
| `green/pass` | `#57C98A` | Pass states, ✓ check icons, status pills |
| `purple` (Mobius Blue) | `#8291C6` | `ori` feature, primer items/roles |
| `feature/green` | `#6EB744` | `lacZα` feature, enzyme items |
| `red/mismatch` | `#EF634F` / `#EF6A53` | Mismatch base highlight, error |
| Spectrum strip | `linear-gradient(90deg,#FFC951,#EF634F,#A54472,#193D66,#5BC2D2)` | 3px brand bar across the very top of the window |

### Color — DNA base coding (used in every sequence/chromatogram view)
`A = #6EB744` (green) · `T = #EF6A53` (red) · `C = #5BC2D2` (cyan) · `G = #FFC951` (gold) · unknown/`N = #9aa7b4`.

### Typography
- **UI font:** `IBM Plex Sans` (400/500/600/700).
- **Mono font:** `IBM Plex Mono` (400/500/600) — **all** sequence text, lengths/bp, enzyme sites, positions, numeric stats.
- Scale: page title `21px/600`, letter-spacing `-0.015em`; section/card title `13px/600`; card subtitle & body `12.5px/400`; list item name `14px/500`; tab label `13px`; section labels (e.g. "REFERENCE LIBRARY") `10.5px/600`, letter-spacing `0.08em`, uppercase; captions `11–11.5px`; sequence bases `15px/500`; big stat numbers `18px/600` mono.

### Spacing, radius, shadow
- Window content padding: `20–26px`. Card padding: `18–20px`. Gaps between cards/rows: `8px` (lists), `14–20px` (sections).
- Radius: cards/panels `14px`; list rows / inputs / buttons `9–11px`; chips/pills `18–20px`; small badges `6–7px`; status circles `50%`.
- Sidebar width `262px`. Right detail/preview panels `316px`. Tab-bar item padding `16px 14px`.
- No heavy shadows in-app (flat dark surfaces separated by 1px borders). Primary buttons: solid `accent` fill, text `#0E1217`.

---

## Global Shell (present on all screens)

**Top spectrum strip:** full-width 3px gradient bar (token above), flush to top edge.

**Left sidebar (`262px`, `#11161c`, right border `#1d242c`):**
- Brand lockup: infinity-loop mark (two overlapping rings, 30×18 — left ring stroke `accent`, right ring stroke `#FFC951`, stroke-width 2.4) + wordmark "Helix" `14.5px/600` and tagline "CLONING WORKBENCH" `10.5px` faint. *(Product name is a placeholder — rename freely.)*
- Search field: `#0E1217` bg, `#232b34` border, radius 9, magnifier icon + placeholder "Search sequences…".
- Two tree sections with `10.5px/600` uppercase labels: **REFERENCE LIBRARY** (Vectors, Enzymes — collapsed, chevron-right) and **PROJECTS** (GFP-reporter — expanded, with children `pUC19-EGFP.gb` [selected: `accent` text, bg `#13202a`, 2px left accent border] and `insert_EGFP.fasta`; CRISPR-library collapsed). Folder icons; selected file uses a file-glyph in `accent`.
- Footer: 28px avatar circle (`#263a77` bg, initials `#bcc8e6`) + user name `12px` muted, top border `#1d242c`.

**Main pane = tab bar + screen body:**
- **Tab bar** (border-bottom `#1d242c`): four tabs — `Reference Library`, `Project Files`, `In-silico Check`, `Wetlab Verify` — each an icon + label. Active tab: `accent` text, `600` weight, 2px `accent` bottom-border (overlapping the bar's border). Inactive: `#8b98a5`, transparent bottom-border. Clicking a tab swaps the body.
- **Header strip** (per screen): breadcrumb (optional) + page title left; action buttons + status pill right; bottom border `#161d24`.
- **Body:** scrollable, padding ~24px.

Tab icons (16px, stroke=currentColor 1.4): Library = rounded rect with a top divider line; Project Files = folder; In-silico = double-helix (two mirrored S-curves); Wetlab = flask.

---

## Screens / Views

### 1. Reference Library
**Purpose:** browse/curate reusable reference parts (vectors, enzymes, primers, features) reused across projects; import new ones; send one into a project.
**Layout:** header strip → body split into a **list (flex:1)** + **detail panel (`316px`, right, self-aligned top)**, `20px` gap.
- **Header:** title "Reference Library", subtitle "Curated backbones, enzymes & primers · reused across projects". Right: ghost button **Import…** + primary **+ Add reference**.
- **Filter row:** chips `All` (active: `accent` bg, `#0E1217` text), `Vectors`, `Enzymes`, `Primers`, `Features` (inactive: `#161d24` bg, `#aeb9c4`); right-aligned "9 references" count, faint.
- **List rows** (`#11161c` bg, `#1f272f` border, radius 11, padding `14px 16px`, gap 8 between): colored 10px square type-dot · name (`14px/500`) + type label (`11px` uppercase faint) in a 150px column · `meta` (mono, `9aa7b4`, ellipsized) · `tags` (`11.5px` muted, right). Selected row: bg `#13202a`, border `rgba(91,194,210,0.30)`.
  Sample items: pUC19 [selected], pET-28a(+), pcDNA3.1(+) (vectors, gold dot); EGFP (feature, accent dot); EcoRI, BamHI, BsaI (enzymes, green dot); M13 Fwd, T7 promoter (primers, purple dot).
- **Detail panel:** "VECTOR · SELECTED" label → "pUC19" `18px/600` + "2,686 bp" mono → **circular plasmid mini-map** (SVG: thin `#222b34` backbone ring + 3 thick rounded colored arcs for AmpR/ori/lacZα, center label) → feature list rows (color square + name + mono length right) → source line ("Source · Addgene · GenBank L09137") above a top border → full-width primary **Use in project →** button.

### 2. Project Files
**Purpose:** load the sequence files for the active project and assign each a role (backbone / insert / primer); preview a file.
**Layout:** header → body split **list (flex:1)** + **preview panel (`316px`)**.
- **Header:** breadcrumb "Projects / GFP-reporter", title "Project Files". Right: ghost **New project** + primary **Load files…**.
- **Drop zone:** `1.5px dashed #2e3a45`, radius 13, `#10151b` bg, centered upload-arrow icon + "Drop sequence files here — or **browse**" (browse in accent) + mono hint ".gb · .gbk · .fasta · .dna · .seq · .ape".
- **Loaded files** (`#11161c` rows, border `#1f272f`): file-glyph (colored per role) · file name (mono `13.5px`, 190px, ellipsized) · format badge (`#161d24` pill) · length (mono) · **role pill** right (tinted by role — Backbone=gold tint, Insert=accent tint, Primer=purple tint — with a dropdown caret) · green ✓ status. Sample: pUC19-backbone.gb (Backbone), insert_EGFP.fasta (Insert), primer_F.seq / primer_R.seq (Primer).
- **Preview panel:** "PREVIEW" label → file name (mono) + "GenBank · 2,686 bp · circular" → **linear feature bar** (14px tall flex row of colored segments) with `1 … 2,686` mono end-labels → "First bases" mono block (`#0E1217` inset, base-colored characters, wraps).

### 3. In-silico Check  *(the hero screen — already approved)*
**Purpose:** verify a simulated assembly (Gibson / restriction-ligation / Golden Gate) is correct before ordering the insert.
**Layout:** header → top region `grid 392px / 1fr` (plasmid map card + checks/stats column) → full-width junction viewer card.
- **Header:** breadcrumb to `pUC19-EGFP.gb`; title "pUC19–EGFP" + mono meta "4,238 bp · circular · Gibson". Right: green status pill "5 / 6 passed" + ghost **Export .gb** + primary **Order insert →**.
- **Plasmid map card:** "Simulated construct map" → large circular SVG (backbone ring + 4 colored feature arcs: EGFP insert=accent, lacZα=green, AmpR=gold, ori=purple; tick-marks + labels for EcoRI & BamHI cut sites; center name + bp) → wrapping legend (color square + name [+ bp]).
- **Verification checks:** title + "in-silico assembly · Gibson". Each check row (`#11161c`, border `#1f272f`, radius 10): status circle (✓ green on `rgba(87,201,138,.12)`, or `!` gold on `rgba(255,201,81,.12)`) · label (`13px/500`) · mono detail right. Six checks: Insert orientation, Reading frame, EcoRI junction, BamHI junction, Premature stop codons (all pass); Internal cut sites (warn — "1 internal EcoRI").
- **Stat tiles:** 4-col grid — 4,238 total bp · 717 insert bp · 52% GC · 238 ORF aa (big mono number + faint label).
- **Junction viewer:** "5′ junction — EcoRI · vector ↔ insert" + A/T/C/G base-color legend. Inside a horizontal-scroll, fixed-width block: an **annotation track** (left "pUC19 MCS" on `#1a222b`, right "EGFP CDS →" on accent tint) above the **sequence row** — base characters (mono 15px, colored per base, 13px wide cells), split by a 2px dashed `accent` vertical divider at the junction.

### 4. Wetlab Verify
**Purpose:** confirm the finished clone matches the expected construct using Sanger + whole-plasmid sequencing; surface variants.
**Layout:** header → verdict banner → full-width coverage card → `grid 1fr / 360px` (alignment card + [chromatogram card over per-feature card]).
- **Header:** title "Wetlab Verification", subtitle "Sanger + whole-plasmid reads vs. expected pUC19–EGFP". Right: gold pill "98.7% identity" (mono) + ghost **Upload reads…**.
- **Verdict banner:** gold-tinted gradient bg (`rgba(255,201,81,.12)`→`.03`), gold border, radius 14. 42px gold check-circle + "Construct confirmed — 1 silent variant" `16px/600` + sub "EGFP CDS intact · c.345C>T is synonymous (no aa change) · review recommended" + gold **Approve clone** button. *(Pass-clean state would use the green palette instead.)*
- **Coverage card:** "Read coverage across construct" + mono "mean depth 142×". SVG (viewBox 1000×150, preserveAspectRatio none): filled accent coverage area with a dip; dashed-gold "low" callout rect over the dip; below it a feature bar of colored segments (AmpR gold / ori purple / EGFP accent) with dark labels.
- **Alignment card:** "Read alignment · EGFP CDS region" + mono "pos 331–357". Horizontal-scroll block: **Expected ref** row (120px label + base-colored mono chars, 15px cells) over **Read F · Sanger** row (same), where mismatched bases render red on `rgba(239,106,83,.22)` with a rounded highlight; a gold caret note "↑ c.345C>T · silent (Gly)" beneath.
- **Chromatogram card:** "Chromatogram · variant call" → SVG (340×110) baseline + dashed-gold highlight rect over the variant base + one colored bell-curve per base (`M cx-16 92 Q cx 92-h cx+16 92`, stroke = base color, width 2.4) → base letters row (mono, base-colored, 40px cells) → caption "Clean single peak at variant — high-confidence call".
- **Per-feature result card:** "Per-feature result" → compact rows (status circle + feature name + mono detail): AmpR ✓, ori ✓, EGFP CDS ! (silent c.345C>T), lac promoter ✓.

---

## Interactions & Behavior
- **Tab navigation:** clicking a top tab sets the active screen and re-renders the body. Active styling per the Tab bar spec. (Cursor: pointer on tabs.)
- **Selection:** clicking a Reference Library row selects it and populates the detail panel (prototype shows the selected state on pUC19).
- **Buttons:** primary = solid `accent`, dark text; ghost = transparent, `#2a333e` border, `#cdd6df` text. All buttons `white-space:nowrap`. Hover (to implement in code): primary → slightly lighten/raise; ghost → border `accent`, text `#E7EDF3`.
- **Role pill (Project Files):** opens a menu to set Backbone / Insert / Primer (caret indicates a dropdown).
- **Horizontal scroll:** sequence/alignment blocks scroll horizontally when wider than their card; rest of layout is fixed-pane (no page-level horizontal scroll on desktop widths).
- **Status semantics:** green = pass/confirmed; gold = warning/needs-review; red = mismatch/error. Keep these consistent everywhere.
- No animation is required for v1 beyond standard hover/focus transitions (~120–160ms ease).

## State Management
- `activeTab`: `'library' | 'files' | 'check' | 'wetlab'` (default `'library'`).
- `selectedReference`: id of the highlighted library item → drives the detail panel.
- `activeProject` + `projectFiles[]`: each `{ name, format, length, role, parsed }`; `role` is editable.
- `selectedFile`: drives the Project Files preview.
- In-silico result model: `{ construct, lengthBp, insertBp, gc, orfAa, method, checks[] }` where each check is `{ label, state: 'pass'|'warn'|'fail', detail }`.
- Wetlab result model: `{ identity, meanDepth, verdict, variants[], coverage[], alignment: { ref, reads[] }, features[] }`.
- `accentTheme`: `'Cyan' | 'Gold' | 'Green'` → sets the CSS variable `--accent` (`#5BC2D2 / #FFC951 / #57C98A`). All accent usage reads this variable.
- Data fetching: parse uploaded sequence files (GenBank/FASTA/.seq/.dna/.ape) into the feature/sequence models; run the in-silico assembly simulation; align uploaded reads to the expected construct for the wetlab view.

## Design Tokens — quick reference (CSS variables to seed)
```
--bg-app:#0E1217;  --bg-panel:#11161c;  --surface-chip:#161d24;
--border-card:#1f272f;  --border-strong:#2a333e;  --border-soft:#1d242c;  --divider-hair:#161d24;
--text-primary:#E7EDF3;  --text-secondary:#cdd6df;  --text-body:#9aa7b4;  --text-muted:#7e8c99;  --text-faint:#5e6b78;
--accent:#5BC2D2;  --gold:#FFC951;  --green:#57C98A;  --purple:#8291C6;  --feature-green:#6EB744;  --red:#EF6A53;
--base-A:#6EB744;  --base-T:#EF6A53;  --base-C:#5BC2D2;  --base-G:#FFC951;
--radius-card:14px; --radius-row:11px; --radius-pill:20px;
font-ui:'IBM Plex Sans'; font-mono:'IBM Plex Mono';
```

## Assets
- **Fonts:** IBM Plex Sans + IBM Plex Mono (Google Fonts, or bundle locally for an offline desktop app).
- **Icons:** all inline SVG line icons (sidebar chevrons/folders/files, search, tab icons, upload arrow, checkmark) — recreate with the codebase's icon set or keep as inline SVG. No raster assets.
- **Plasmid maps, coverage chart, chromatogram:** all drawn as inline SVG from data — reimplement as data-driven SVG/canvas components.
- **Brand:** infinity-loop mark + top spectrum strip are nods to the Scripps Research brand. Swap in the official logo lockup if/when desired.

## Files
- `Cloning Workbench - A.dc.html` — the high-fidelity prototype (open in a browser; click the tabs to see all four screens). Needs `support.js` + internet (fonts) to render.
- `support.js` — runtime required by the prototype. Not part of the design; do not ship.
