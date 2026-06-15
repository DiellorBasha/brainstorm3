# Wavelet Designer — convert the FilterDesigner into a localized wavelet tool

**Date:** 2026-06-15
**Status:** Design approved; implementation pending
**Branch:** `feat/filterbank-designer` (brainstorm3) + paired `bst-java` fork branch
**Author:** Diellor Basha

## Problem

The interactive designer we built conflates two distinct concepts:

- A **filter** has *no localization*. Its input is a source map (a cortical vector
  field); a spectral kernel (low/high/band-pass, heat) reshapes that field across the
  Dirac eigenvalue spectrum to isolate a scale. No vertex, no direction, no tiling.
  This is exactly `bst_dirac_eigenmodes_filter` / `process_dirac_filter`.
- A **wavelet** is *localized*. A vertex (delta) + a direction + kernel parameters
  define an atom. Like a temporal CWT it has two uses: **synthesis** (reconstruct the
  atom to see it) and **analysis** (inner-product the atom against a source map to get
  a coefficient measuring how well the field matches the wavelet there). Tiling the
  spectrum yields the multi-scale bank — the wavelet-transform analogue.

The current panel is a wavelet synthesizer with a filter's "source map" input bolted
on, and its seed direction is chosen in the *world* frame (azimuth/elevation in MRI
X/Y/Z). Two things are wrong: the filter/wavelet concepts are mixed, and a world-frame
direction is not physically meaningful for a cortical source (it can point into the
skull). The Dirac operator does act on full 3-D ambient vector fields, but the seed
*direction* must be specified in the **local cortical frame** so the chosen orientation
is meaningful; the field it generates is still ambient 3-D.

## Goal

Convert the existing panel into a correct, self-contained **Wavelet Designer**
(synthesis + design + tiling). The pure Filter designer is a *separate future tool* and
is out of scope here. The wavelet *analysis* path (CWT coefficients from a source map)
is deferred to its own later round.

This round delivers:

1. Rename the tool and its DB artifact from "filter(bank)" to "wavelet".
2. Make the input **localization-only** (click a vertex); remove the filter-style
   "Active source map" input.
3. Specify the seed direction in the **local `manifold_` frame** (U,V tangent, N normal)
   instead of world az/el; the embedded seed is a full ambient 3-vector the Dirac diffuses.
4. Keep kernel + scale + chirality + spectrum tiling + save exactly as they are.

## Non-goals

- The pure **Filter designer** (kernel-only, no localization) — separate future tool.
- The **analysis / CWT-coefficients** path (apply a wavelet/bank to a source map →
  coefficient map) — deferred to its own round; the source-map input returns there as
  an *analysis target*, not a filter input.
- Changing the Dirac math (`bst_dirac_filter`, `bst_dirac_eigenmodes_filter`,
  `bst_filterbank_tiles`) — unchanged; they already operate on ambient directions.

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Disposition of the current panel | It **is** the Wavelet Designer; rename + refine, don't rebuild |
| This round's scope | **Correctness first**: rename + local-frame direction + drop filter-style source input + keep synthesis/tiling. Analysis-coefficients deferred |
| Seed direction parameterization | **Two angles in the local frame** {U,V,N}: in-plane φ∈[0,360°) within (U,V), tilt θ∈[−90,90°] toward N; default θ=+90° (outward normal) |
| DB node rename | **Rename** `filterbank_` → `wavelet_` (full schema/DB/tree churn, done now while unmerged) |

## Architecture

### Renames (filter → wavelet)

| Old | New |
|---|---|
| `toolbox/gui/panel_filter_designer.m` | `toolbox/gui/panel_wavelet_designer.m` |
| `toolbox/gui/view_filter_designer.m` | `toolbox/gui/view_wavelet_designer.m` |
| panel name `'FilterDesigner'` | `'WaveletDesigner'` |
| tree menu "Design filterbank…" | "Design wavelet…" |
| `db_template('filterbankmat')` | `db_template('waveletmat')` |
| `db_add_filterbank.m` | `db_add_wavelet.m` |
| `bst_get('FilterbankFile', …)` | `bst_get('WaveletFile', …)` |
| file type `filterbank` (`filterbank_*.mat`) | `wavelet` (`wavelet_*.mat`) |
| `Surface.Filterbank` DB list / `ParentEigen` | `Surface.Wavelet` DB list / `ParentEigen` |
| bst-java node type `filterbank` | `wavelet` |
| tests `test_db_add_filterbank`, `test_filter_designer_session` | `test_db_add_wavelet`, `test_wavelet_designer_session` |

The saved struct keeps its shape (`Tiles`, `Tiling = struct(Wavelet, Opts)`,
`ParentEigen`, `Provenance`); only the type name changes. `bst_filterbank_tiles` keeps
its name (it is still the spectrum-tiling module that consumes one wavelet).

### Input: localization-only (panel section 1)

Remove the `jInputDelta` / `jInputSource` radios and the `SetSeedSource` path. The Input
section becomes: a one-line hint ("click a vertex on the cortex to place the wavelet")
plus the direction and chirality controls. Seeding happens only through the
`figure_3d` click hook (`SetSeedVertex`). `process_dirac_filter` already provides the
no-localization filtering, so nothing is lost.

### Seed direction in the local manifold frame (panel section 1, the real change)

- **Dependency:** the parent surface's `manifold_` node. On session open,
  find-or-create it via `tess_manifold(SurfaceFile)` (mirrors how the eigenbasis is
  found-or-created), using the manifold's default gauge. Store the loaded manifold in
  the panel session state.
- **Per-vertex frame:** derive the combed per-vertex frame `(U,V,N)` from the manifold's
  `Embedded`/`Gauge` via `view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)`
  (already exists). `U,V` span the tangent plane; `N` is the outward normal.
- **Direction control:** the two existing sliders, reinterpreted:
  - In-plane angle φ ∈ [0,360°) (relabel "In-plane angle"),
  - Tilt θ ∈ [−90,90°] toward N (relabel "Tilt → normal"), default θ=+90°.
- **Embedding** at the clicked vertex `v`:
  `d = cosθ·(cosφ·U_v + sinφ·V_v) + sinθ·N_v`, normalized — a full ambient 3-vector.
  This `d` is passed to the existing seed pipeline (`SetSeedVertex` builds the quaternion
  delta `(0,d)` and projects it); the Dirac then diffuses it as a full 3-D field.
- **Gauge caveat:** φ's zero-reference depends on the manifold's tangent comb (gauge);
  θ (toward N) is gauge-invariant. Documented in the panel help; acceptable for an
  exploration tool where φ is swept.

### Unchanged

`BuildWavelet` (kernel/params/chirality + the new local-frame direction), `BuildTiling`,
`Refresh`, `ComputeField`, `OnSave`, the spectrum strip, the live-preview controller,
the projection/filter split — all unchanged except `SeedDirection` (now frame-based) and
the removed source-map input. `figure_3d` click-to-seed hook unchanged.

## Components / files

**brainstorm3 (rename + edit):**
- `panel_wavelet_designer.m` — renamed; `SeedDirection` uses the local frame; input
  radios removed; sliders relabelled; session state holds the manifold frame.
- `view_wavelet_designer.m` — renamed orchestrator; find-or-create the `manifold_` node
  and pass the derived frame to the panel; everything else (preview results node, tab
  activation, source-vector display, teardown) unchanged.
- `db_template.m` — `waveletmat` template; `Surface.Wavelet` list (replaces
  `filterbankmat` / `Surface.Filterbank`).
- `db_update.m` / `bst_startup.m` — DB version bump (5.06) backfilling `Surface.Wavelet`
  via `NormalizeSurfaceArray` (same path as the 5.05 filterbank add). The unmerged 5.05
  `Filterbank` field is superseded; the migration adds `Wavelet`.
- `bst_get.m` — `WaveletFile` accessor (replaces `FilterbankFile`).
- `db_add_wavelet.m` — renamed from `db_add_filterbank`.
- `file_gettype.m` / `file_fullpath.m` — register `wavelet` (`_wavelet_`) in place of
  `filterbank`.
- `tess_manifold.m` — reused (find-or-create the manifold); a thin helper
  (`bst_manifold_frame` or reuse `view_manifold('DeriveVertexFrame')`) returns `(U,V,N)`.
- `node_create_subject.m` / `tree_callbacks.m` — nest `wavelet_` under the eigen node;
  "Design wavelet…" menu; wavelet popup (Edit/Delete); cascade delete.

**bst-java fork:** `wavelet` node type (renders with the default icon, like the others).

**Tests (rename + extend):**
- `test_wavelet_schema` (was filterbank schema) — `waveletmat` + `Surface.Wavelet`.
- `test_db_add_wavelet` — round-trip under the eigen node.
- `test_filterbank_tiles` — unchanged (the tiling module keeps its name).
- `test_dirac_filter_coeffs` — unchanged.
- `test_wavelet_designer_session` — open → seed (local-frame direction) → save → nested
  node → delete; plus a check that the seed direction equals
  `cosθ(cosφ·U+sinφ·V)+sinθ·N` at the clicked vertex.

## Data flow (seed)

1. User clicks vertex `v` → `figure_3d` hook → `panel_wavelet_designer('SetSeedVertex', v)`.
2. `SeedDirection`: read φ,θ from the sliders; fetch `(U_v,V_v,N_v)` from the cached
   manifold frame; `d = normalize(cosθ(cosφ·U_v+sinφ·V_v)+sinθ·N_v)`.
3. Build the quaternion delta `(0,d)` at `v`; project once (`ReturnCoeffs`); cache.
4. `Refresh` reconstructs the atom (single wavelet) or, if tiling is on, the selected
   tile, and pushes it to the preview (unchanged).

## Error handling

- **No manifold / manifold compute fails:** the find-or-create reports via
  `bst_progress`/`bst_error`; if it cannot be produced, the session aborts cleanly with a
  message (the designer requires a frame to place a directional seed).
- **Eigenbasis vs manifold vertex mismatch:** both derive from the same parent surface;
  assert vertex counts match before embedding, else abort with a clear message.
- **Click outside the eigen support:** unchanged (existing `vertexNotFound` guard).
- **Teardown:** unchanged idempotent teardown (removes the transient preview node, closes
  the figure and spectrum strip).

## Testing

- **Schema/DB round-trip** (`test_wavelet_schema`, `test_db_add_wavelet`): the renamed
  template/accessor/add path save and resolve a `wavelet_` node nested under the eigen
  node; reload is identical.
- **Direction embedding** (headless, in `test_wavelet_designer_session` or a focused
  test): for known `(U,V,N)` and φ,θ, the seed equals
  `cosθ(cosφ·U+sinφ·V)+sinθ·N` (unit), and θ=+90° gives exactly `N`.
- **Tiling module** (`test_filterbank_tiles`): unchanged — still consumes one wavelet +
  opts.
- **Session lifecycle** (`test_wavelet_designer_session`): open → seed → save → nested
  `wavelet_` node → delete; full teardown, no orphan preview node.

## Build order

1. **Rename pass** (mechanical): files, identifiers, panel/tab names, menu text, tests.
   Keep behaviour identical; run the renamed suite green.
2. **DB rename + migration** (`waveletmat`, `Surface.Wavelet`, `WaveletFile`,
   `db_add_wavelet`, `file_gettype`/`file_fullpath`, v5.06) + bst-java `wavelet` node.
3. **Input localization-only**: remove the source-map radios + `SetSeedSource`.
4. **Local-frame direction**: manifold find-or-create on open, frame derivation,
   `SeedDirection` rewrite, slider relabel + tilt default, the embedding test.
5. End-to-end session test on the renamed, frame-seeded designer.
