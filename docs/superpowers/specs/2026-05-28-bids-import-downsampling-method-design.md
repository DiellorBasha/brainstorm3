# BIDS Import — Cortex Downsampling Method Design

**Date:** 2026-05-28
**Status:** Approved (brainstorming)
**Author:** Diellor Basha

## Goal

Let the Brainstorm BIDS importer (`process_import_bids`) choose the cortex
downsampling method — `reducepatch` (legacy) or `icosphere` (manifold,
FreeSurfer/MNE-style) — and default it to **icosphere / ico5** so a fresh import
produces the **20484-vertex, 2-component, 2-manifold** cortex the eigenmode
pipeline expects. This is the prerequisite for rebuilding the OMEGA_Tutorial
protocol from BIDS with correct ico-downsampled anatomy.

## Background

- The icosphere machinery already exists and is exercised today: `tess_downsize`
  has an `'icosphere'` method, and `import_anatomy_fs` accepts a `Method` argument
  (arg #9) plus an interactive dialog (method radio → ico-level radio). See
  `dev/ico-downsize.md`.
- **The gap:** `process_import_bids` calls `import_anatomy_fs` with only 6 args
  (`process_import_bids.m:486`), so `Method` defaults to `'reducepatch'`. The BIDS
  importer therefore *cannot* produce icosphere anatomy regardless of options.
- The BIDS importer's Run path is **non-interactive** (`OPTIONS.isInteractive = 0`
  at line 143; `isInteractiveAnat = 0` at line 458), so `import_anatomy_fs`'s own
  dialog never fires through the pipeline. The user-facing surface is the **process
  options panel** (the `nvertices` / `mni` / `anatregister` radios at lines 94–115),
  not a popup.

## Verified contracts (read from source)

- `import_anatomy_fs` signature:
  `import_anatomy_fs(iSubject, FsDir, nVertices, isInteractive, sFid, isExtraMaps=0, isVolumeAtlas=1, isKeepMri=0, Method=[])`.
  The current BIDS call passes 6 args, so the effective defaults in play are
  `isExtraMaps=0, isVolumeAtlas=1, isKeepMri=0, Method='reducepatch'`.
- Non-interactive icosphere path (`import_anatomy_fs.m:141–156`): with
  `Method='icosphere'` and a non-empty `nVertices`, it computes
  `nVertHemi = round(nVertices/2)` and `tess_downsize` snaps to the nearest ico
  count. So passing `nVertices = 20484` with `Method='icosphere'` yields per-hemi
  10242 → **ico5 exactly**.
- `process_import_bids` uses `eval(macro_method)` dispatch, so new subfunctions are
  callable by string name (e.g. `process_import_bids('GetIcoVertexCount','ico5')`).

## Scope

**In scope:** changes to `toolbox/process/functions/process_import_bids.m` only,
plus tests. No changes to `import_anatomy_fs` / `tess_downsize` (already done).

**Out of scope (deferred, flagged):**
- Extending icosphere to the non-FreeSurfer per-hemisphere importers
  (`import_anatomy_cat/bs/bv/civet`). They share the same `tess_downsize`
  reducepatch pattern; generalization is future work.
- The full OMEGA_Tutorial protocol rebuild (recordings → covariances → head model
  → eigenmodes). That is the **next** increment, which *consumes* this change.

## Design

### Component 1 — Options panel (`GetDescription`)

Add two `radio_linelabel` controls that mirror `import_anatomy_fs`'s own dialog:

- `downsamplemethod`: `Reducepatch` / `Icosphere` — **default `'icosphere'`**.
- `icolevel`: `ico3` / `ico4` / `ico5` / `ico6` — **default `'ico5'`**.

Relabel the existing `nvertices` control to
`'Number of vertices (cortex, reducepatch): '` to signal it is consulted only on
the reducepatch path. Its default (15000) and type are unchanged, preserving the
reducepatch fallback count.

### Component 2 — Pure helper `GetIcoVertexCount(level)`

New macro-dispatched subfunction (`%#ok<DEFNU>`), pure and unit-testable:

| level | total vertices |
|-------|----------------|
| ico3  | 1284  |
| ico4  | 5124  |
| ico5  | 20484 |
| ico6  | 81924 |

Errors on an unrecognized level. Used to translate the `icolevel` choice into the
`nVertices` value handed to `import_anatomy_fs` (which then snaps via `round/2`).

### Component 3 — `Run` option wiring

Read the two new options into `OPTIONS.DownsampleMethod` and `OPTIONS.IcoLevel`
alongside the existing `OPTIONS.nVertices`.

### Component 4 — `ImportBidsDataset` defaults + interactive path

- `Def_OPTIONS` gains `DownsampleMethod='icosphere'` and `IcoLevel='ico5'`.
- Interactive path (`isInteractive=1`, the `java_dialog` at line ~469): ask a
  **method** radio first; if icosphere → ask an **ico-level** radio; if reducepatch
  → keep the existing vertex-count input. Asked once for the whole import (as today).

### Component 5 — Anatomy import call site

At the format `switch` (lines 484–497):

- **FreeSurfer:** if `DownsampleMethod=='icosphere'`, set
  `nVertArg = GetIcoVertexCount(IcoLevel)` and `methodArg='icosphere'`; else
  `nVertArg = OPTIONS.nVertices`, `methodArg='reducepatch'`. Extend the call to nine
  args:
  `import_anatomy_fs(iSubject, dir, nVertArg, isInteractiveAnat, sMriFid, 0, 1, 0, methodArg)`.
  The `0,1,0` are the verified current effective defaults, so reducepatch behavior is
  preserved exactly.
- **CAT12 / BrainSuite / BrainVISA / CIVET:** if icosphere was requested, append a
  `Warning` to `Messages` ("Icosphere downsampling is currently FreeSurfer-only;
  subject `<X>` (`<format>`) imported with reducepatch at `<nVertices>` vertices.")
  and call the existing importer unchanged with `OPTIONS.nVertices`. A code comment
  marks the deferred extension.

## Error handling

- `GetIcoVertexCount` errors on an unrecognized level (guards typos / future levels).
- Existing `nVertices` validation (`process_import_bids.m:138`, `< 50` → error) is
  retained for the reducepatch path.
- Non-FreeSurfer + icosphere is a **warning, not an error** — mixed-format datasets
  still import fully (FreeSurfer subjects get icosphere; others get reducepatch).

## Testing

1. **`dev/tests/test_process_import_bids_options.m`** (new; mirrors the existing
   `test_process_eigenmodes_*_options.m` pattern): load `GetDescription`; assert
   `downsamplemethod` is a radio defaulting to `'icosphere'`, `icolevel` is a radio
   defaulting to `'ico5'`, and `nvertices` still exists.
2. **Pure-helper assertions** (same test file): `GetIcoVertexCount` returns the four
   level→count mappings (ico5 → 20484, etc.) and errors on an unknown level.
3. **Static analysis:** `check_matlab_code` on the modified file — no new findings
   beyond pre-existing lint.
4. **Live check (one FreeSurfer subject):** run the BIDS import with default options
   on OMEGA sub-0002's anatomy and confirm the cortex is **20484 vertices, 2
   components, manifold** — proof the default reaches
   `tess_downsize('icosphere', 10242)`. Doubles as the first slice of the eventual
   OMEGA rebuild.

## Consequences

- `tutorial_omega.m` calls this importer via `CallProcess` without setting a method,
  so `bst_process` fills the `GetDescription` defaults → it now produces ico5
  anatomy automatically. The OMEGA rebuild increment can reuse it directly.
- Existing reducepatch-based imports are unaffected only if the user explicitly
  selects `reducepatch`; the **default behavior of the fork's BIDS importer changes
  to icosphere/ico5** (intended).
