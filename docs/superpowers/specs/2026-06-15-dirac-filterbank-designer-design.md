# Spatial filterbank designer — interactive design session

**Date:** 2026-06-15
**Status:** Design approved; implementation pending
**Branch:** `feat/filterbank-designer` (brainstorm3) + a paired branch on the `bst-java` fork
**Author:** Diellor Basha

## Problem

The Dirac (and scalar Laplace–Beltrami / connection-Laplacian) spatial filterbank
has a complete *compute* layer — `bst_dirac_eigenmodes_filter`, `bst_eigenmodes_filter`,
the `bst_eigfilter_*` kernel registry, the synthesis wavelet `bst_dirac_filter`, and
the batch `process_dirac_filter` / `process_eigenmodes_filter`. What is missing is an
**interactive designer**: a way to pick a kernel from a dropdown, click a delta on the
cortical surface, choose operator-specific parameters (for Dirac: direction and
chirality), see the resulting field live on the cortex, tile the eigenmode spectrum
into a small bank of wavelets, and save that bank as a reusable artifact.

Today the only interactive feedback is `view_eigfilter_response`, a 1-D plot of
`g(λ)` with no cortical surface. There is no delta-click, no live preview loop, no
direction/chirality picker, and no saved filterbank object.

## Goal

A **transient interactive design session** that:

1. opens on an existing eigenbasis (`eigen_` node), reusing `figure_3d` for the cortex;
2. lets the user pick a kernel and parameters, a delta seed (or a loaded source map),
   and — for Dirac — a direction and chirality;
3. previews the filtered field live on the cortex as parameters change;
4. tiles the eigenmode spectrum into a small bank of wavelets (× chiralities);
5. saves the bank as a `filterbank_` node nested under the `eigen_` node, reusable on
   any data via the existing compute cores;
6. tears itself down completely on Save or Cancel (no permanent docked panel).

The architecture is **operator-general** (scale `g(λ)` is shared across all three
operators) with **Dirac wired end-to-end first**; LBO and connection-Laplacian are
registered-but-deferred render modes (the panel simply hides the direction/chirality
controls when the eigen node's `Variant` is not `Dirac`).

## Non-goals

- Small-multiples montage of all tiles at once (the chosen layout is one selected
  tile + a spectrum strip).
- A permanent always-docked designer panel.
- Materializing/storing filtered fields in the saved node (the node stores the recipe
  bank only; fields are regenerated on demand).
- Designing temporal filters (this is the *spatial* eigenmode domain).
- Resolving the connection-Laplacian PSD/phase rendering (deferred).

## Decisions (from brainstorming)

| Question | Decision |
|---|---|
| Operator scope | General framework, **Dirac wired first**; LBO/connection stubbed |
| Input field | **Delta OR loaded source map** (toggle) |
| Saved artifact | A **filterbank**: the user designs one filter, then a tiling generator adds more filters spanning the eigenmode spectrum (e.g. 4 wavelets) × optional chiralities; stored as a **recipe bank** (no materialized fields) |
| Bank display | **Selected tile on a `figure_3d` surface + a `g(λ)` spectrum strip** showing all tiles; click a curve to switch the displayed tile |
| Render architecture | **Reuse `figure_3d`** via a temporary results file + a docked control panel (Approach A) |
| Panel lifetime | **Transient session**, not permanent: docked only while designing; on Save/Cancel the panel undocks and the figures close |
| Node nesting | `filterbank_` is a **true child of the `eigen_` node** (new derived-anatomy nesting level → requires bst-java fork work) |

## Architecture

Three brainstorm3 GUI components with a clean split of responsibility, plus the
schema/DB/tree plumbing and the bst-java fork node type.

### Components and responsibilities

**`toolbox/gui/view_filter_designer.m` — session orchestrator (glue).**
A `view_*` entry point, called once to start a session. Owns lifecycle, not widgets
or math:
- finds-or-loads the eigenbasis via `tess_eigen` (now cached → fast);
- opens a `figure_3d` of the parent cortex (a temporary results file holds the
  preview field);
- docks the control panel with
  `gui_show('panel_filter_designer', 'BrainstormTab', 'Filter designer', …)`;
- links panel ↔ figure (passes the figure handle to the panel; sets the figure's
  `CloseRequestFcn` so closing the figure ends the session);
- runs teardown on Save/Cancel: `gui_hide('FilterDesigner')`, close the preview
  figure, `view_eigfilter_response('close')`, delete the temporary results file;
- enforces a singleton session (one designer at a time, keyed by panel tag).

**`toolbox/gui/panel_filter_designer.m` — control panel (state + interaction + live controller).**
A `panel_*` module (`CreatePanel()` builds Java widgets and registers callbacks). The
only component that persists for the whole session and holds session state:
- widgets: input radio (Delta / Source map), kernel dropdown
  (`bst_eigfilter_kernel('list')`), auto-generated parameter controls
  (`bst_eigfilter_kernel('info',name).params` → one slider/value-box per param using
  its `default` and `range`), direction (normal / tangent / custom xyz) and chirality
  (None / + / −, axis = normal default) — both shown only when `Variant=='Dirac'`,
  tiling (Tiles N + geometric/dyadic spacing, optional × chiralities), Save / Cancel;
- state: the current design, the active tile index, and the **cached projection**
  `c = Phiᵀ·B·ψ`;
- live-preview controller: on any control change, regenerate the tile gains `hⱼ(λ)`,
  reconstruct the **selected** tile only (`diag(hⱼ)·c` + optional `P±` for Dirac), push
  the field into the preview results file, and refresh `figure_3d`.

If the panel grows unwieldy (widgets + state + controller in one file), the planned
split is to extract the live-preview controller (projection cache, tile
reconstruction, figure refresh) into a `bst_filter_designer_state.m` helper, leaving
`panel_filter_designer` as widgets + callbacks. Start unified; split only if earned.

**`toolbox/gui/view_eigfilter_response.m` — spectrum strip (passive display, EXTENDED).**
Already exists as a 1-curve `g(λ)` plot with no state and no cortex knowledge. Extend
from "plot one curve" to "plot all N tiles, highlight the active one, and select a tile
on curve-click" (the click calls back into the panel). Remains driveable standalone
with any kernel handle, exactly as today.

### Component relationships

```
view_filter_designer  ── opens ──▶  figure_3d (cortex preview, temp results file)
        │                                 ▲
        │ docks                           │ pushes filtered field
        ▼                                 │
panel_filter_designer ───────────────────┘   (state + live controller)
        │  draws / reads curve-clicks
        ▼
view_eigfilter_response  (g(λ) tiles strip)
```

## Data flow — live preview (projection / filter split)

The interactive loop is cheap because projection is separated from filtering:

1. **Seed** (delta-click or source-load, once per seed): embed the input as a
   pure-imaginary quaternion field `ψ` and compute `c = Phiᵀ·B·ψ`. Cache `c` in the
   panel state. (For a delta: `ψ` is the direction quaternion at the clicked vertex;
   for a source map: `ψ` is the selected time frame.)
2. **Filter** (on every kernel/param/direction/chirality/tile change): regenerate the
   per-tile gains `hⱼ(λ)` from the eigfilter registry; reconstruct only the selected
   tile, `Jfilt = reconstruct(diag(hⱼ)·c)` with the optional helicity projector `P±`
   for Dirac. This is a matrix multiply on the already-cached basis → interactive on
   slider drags.
3. **Display:** write `Jfilt` into the preview's temporary results file and refresh the
   `figure_3d` (data + colormap + the source-vector quiver override). The spectrum strip
   redraws all N tile curves with the active one highlighted.

The seed projection is the only O(K·nVert) step and runs once per seed; all live
parameter changes are O(K) gain evaluation + one reconstruction.

## Saved artifact — `filterbank_*.mat`

A recipe bank (no materialized fields), saved as a true child of the `eigen_` node:

- `Variant` — operator variant inherited from the eigen node (`Dirac`, etc.);
- `EigenFile` — reference to the parent `eigen_` node (the basis the recipes apply in);
- `ParentEigen` — DB path of the parent eigen node (for nesting);
- `Tiles(j)` — one entry per tile: `{Kernel, Params, Direction, Chirality, Axis}`;
- `Tiling` — generator metadata (kernel family, N tiles, spacing rule, λ-range,
  chirality set) so the bank can be regenerated/edited;
- `Comment`, `Provenance` (compute date, nxr version, design vertex if delta-seeded).

The bank is tiny and reapplicable to any data via `bst_dirac_eigenmodes_filter` /
`process_dirac_filter`. Double-clicking the node re-opens the designer pre-loaded from
the saved recipe.

## Tree integration and schema

### brainstorm3
- `toolbox/db/db_template.m` — add a `filterbankmat` template (the schema above) and a
  `Filterbank` child list on the eigen node entry; bump the `db_update` DB version with
  a migration that backfills the empty `Filterbank` field on existing eigen nodes.
- `toolbox/db/db_add_filterbank.m` — save a `filterbank_*.mat` alongside the anatomy
  folder and register it as a child of the parent `eigen_` node (mirrors
  `db_add_eigen`, but nests under the eigen node rather than the surface).
- `toolbox/tree/tree_callbacks.m` — "Design filterbank…" context item on the `eigen_`
  node (launches `view_filter_designer`); populate the new `filterbank_` child level;
  double-click a `filterbank_` node → re-open the designer pre-loaded.

### bst-java (fork)
- New `filterbank` tree node type + icon in the Java tree renderer.
- Render `filterbank_` nodes as children of `eigen_` nodes (the new nesting level).
- Dev-only on a feature branch of the fork; never pushed/merged upstream.

## Reused, not rebuilt

`figure_3d` (surface + source-vector quiver + true-size toggle + Snapshot),
`bst_dirac_eigenmodes_filter` (project → scale → chirality → reconstruct core),
`bst_eigfilter_kernel` registry (`list`/`info` introspection drives the auto-form),
`tess_eigen` (cached find-or-load basis access), `gui_show`/`gui_hide` (transient
docking), and the derived-anatomy node pattern (`manifold_`/`eigen_`/`operator_`).

## Error handling

- **No eigenbasis / wrong variant:** the entry point requires a resolved `eigen_` node;
  if absent it offers to create one via `tess_eigen` (with progress), else aborts.
- **Source-map input mismatch:** a loaded source map must be unconstrained
  (3 components/vertex) and span the same vertices as the basis; otherwise the input
  toggle disables "Source map" with a tooltip explaining why (mirrors
  `process_dirac_filter`'s checks).
- **Delta outside eigen support:** clicking a vertex not in the basis support warns and
  is ignored (mirrors `bst_dirac_filter`'s `vertexNotFound`).
- **Out-of-range params:** the auto-generated controls clamp to each param's `range`
  from the registry metadata; bandpass/non-prior-admissible kernels are flagged from
  the `bandpass`/`priorAdmissible` metadata.
- **Teardown safety:** teardown is idempotent and runs from a single path used by Save,
  Cancel, and the figure `CloseRequestFcn`, so closing any window cleans up the rest and
  deletes the temporary results file (no orphaned docked panel or temp file).

## Testing

- **`bst_eigfilter_kernel` introspection:** every registered kernel returns
  well-formed `info` metadata (`display`, `params{default,range}` for each param) so the
  auto-form has no gaps. Pure, no GUI.
- **Projection/filter split (headless):** drive the live-preview controller without a
  figure — seed once, vary params, assert each tile equals
  `bst_dirac_eigenmodes_filter` applied directly (single source of truth), and that
  re-seeding is the only projection performed.
- **Tiling generator:** N tiles span the requested λ-range with the expected spacing;
  `× chiralities` doubles the bank with matched `±` entries; flat/identity tile round-
  trips the seed.
- **Save / reload round-trip:** a saved `filterbank_` node reloads to an identical
  recipe bank, nests under the correct `eigen_` node, and reproduces the same fields when
  reapplied to the same data.
- **Session lifecycle:** open → Save and open → Cancel both fully tear down (panel
  hidden, figures closed, temp results file deleted, no orphan handles), verified by
  handle/file scans before and after.

## Build order

1. **Schema + DB + save/reload** (`db_template`, `db_add_filterbank`, migration) with the
   headless round-trip test — establishes the artifact independently of any GUI.
2. **Live-preview controller** (projection/filter split, tiling generator) with the
   headless controller test — the design math, testable without a figure.
3. **Spectrum strip** extension of `view_eigfilter_response` (multi-tile + clickable).
4. **Control panel** `panel_filter_designer` wiring controllers + widgets to the figure.
5. **Orchestrator** `view_filter_designer` (open/link/teardown) + `tree_callbacks` entry.
6. **bst-java fork** node type/icon + eigen-child rendering.
7. End-to-end session lifecycle test.
