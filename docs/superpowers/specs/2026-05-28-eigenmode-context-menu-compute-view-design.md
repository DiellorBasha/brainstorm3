# Eigenmode Context Menus: Compute + Interactive Viewer — Design

**Date:** 2026-05-28
**Status:** Approved design, ready for implementation plan

## Background

The scalar LBO eigenmode pipeline exists end-to-end, but every feature is
reachable **only** through the Process tab (drag files → Run → pick process).
There is zero context-menu (right-click) integration: a `grep` across
`toolbox/gui/` and `tree_callbacks.m` finds no eigenmode references. The only
process that produces an on-surface visualization is `process_eigenmodes_view`,
which writes a persistent Results node whose "time" axis is overloaded as the
mode index.

This is the first of several incremental additions that make the eigenmode
methods explorable by clicking. It covers **two coupled additions shipped as
one increment** (compute produces no visible artifact on its own — it can only
be validated by viewing the modes):

1. Right-click a **cortex** surface → **Compute eigenmodes** (brief dialog).
2. Right-click a **cortex** surface → **View eigenmodes** — a dedicated,
   transient viewer modeled on `view_leadfield_sensitivity.m`: **←/→ steps
   through modes**, the cortex colormap updates live, and the legend shows the
   mode index, eigenvalue, and spatial wavelength.

### Reference implementation (the visualization template)

`toolbox/gui/view_leadfield_sensitivity.m` (`DisplayMode='Surface'`) is the
proven pattern we mirror:

- Opens a surface figure with the source-overlay pipeline.
- Backs up the figure `KeyPressFcn` and installs a custom handler; arrow keys
  step the displayed channel (we step the mode index instead).
- Each step recomputes a per-vertex scalar, writes it into the figure's
  `TessInfo(1).Data` / `DataLimitValue` / `DataMinMax` appdata, and refreshes
  via `panel_surface('UpdateSurfaceColormap', hFig)`.
- A bottom-left `uicontrol` legend reports the current selection; an `H` key
  shows a help popup; unhandled keys fall through to the backed-up callback
  (preserving rotate/zoom).

### Existing building blocks reused

- `process_eigenmodes.m` — `Run` already loads the surface, validates manifold
  (`tess_manifold`), computes `tess_eigenmodes(V, F, ...)`, and stores via
  `out_tess_eigenmodes`. Indexed under *Import > Import anatomy* (Index 25).
- `in_tess_eigenmodes(SurfaceFile)` → `[Eigenmodes, isComputed]`; `Eigenmodes`
  has `.Vectors [nVert × nModes]`, `.Values [nModes × 1]`, `.nModes`,
  `.MassType`.
- `tree_callbacks.m` surface-node menu block: `case {'scalp','cortex',...}`
  (~line 1162). `Display` item at ~1171; cortex-specific items (`Extract
  envelope`, …) inside `if ~bst_get('ReadOnly')` at ~1263.

## Goals

- Two context-menu items on the cortex node: *Compute eigenmodes*,
  *View eigenmodes*.
- A leadfield-style transient viewer with arrow-key mode stepping and live
  colormap update.
- Logic split into pure, headlessly testable helpers; thin GUI wrappers.
- No regression to the existing `process_eigenmodes` Run path.

## Non-goals (YAGNI)

- Surface types other than cortex (the process supports them; the menu does not
  expose them yet).
- Saving anything to the database from the viewer (it is transient).
- Multi-mode montage / side-by-side modes.
- A repair-if-nonmanifold toggle in the dialog (repair stays off).

## Design

### Component 1 — Menu hooks (`toolbox/tree/tree_callbacks.m`)

In the surface-node block `case {'scalp','cortex','outerskull','innerskull','other'}`,
**scoped to `strcmpi(nodeType,'cortex')`**, add a separator and two items
grouped together immediately after the existing `Display` item (~line 1171):

- **Compute eigenmodes** — gated by `~bst_get('ReadOnly')`:
  `gui_component('MenuItem', jPopup, [], 'Compute eigenmodes', IconLoader.ICON_SURFACE_CORTEX, [], @(h,ev)bst_call(@process_eigenmodes, 'ComputeInteractive', iSubject, filenameRelative))`
- **View eigenmodes** — always shown (works in read-only):
  `gui_component('MenuItem', jPopup, [], 'View eigenmodes', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_eigenmodes, filenameRelative))`

*View eigenmodes* is always enabled; if no modes exist it shows a friendly
error pointing to Compute (so the menu build does not read the surface file on
every right-click).

### Component 2 — `process_eigenmodes('ComputeInteractive', iSubject, SurfaceFile)`

New macro-dispatched method following the `ComputeInteractive` idiom
(cf. `process_fem_mesh`):

1. `java_dialog` collects **number of eigenmodes** (default 300) and **mass
   type** (barycentric / Voronoi / Galerkin). Cancel → return (no-op).
2. If the surface already has eigenmodes (`in_tess_eigenmodes`), confirm
   overwrite; decline → no-op.
3. Call the shared compute core with RemoveDC=true, Repair=false, under
   `bst_progress`. Errors (e.g. non-manifold) surface via `bst_error`.

**Refactor:** extract the existing `Run` body into a shared core
`Compute(SurfaceFile, nModes, MassType, RemoveDC, Repair)` returning the stored
`Eigenmodes` struct (or erroring). Both `Run` and `ComputeInteractive` call it.
This must be behavior-preserving for the existing process — verified by the
existing `process_eigenmodes` tests plus the round-trip test below.

### Component 3 — `view_eigenmodes.m` (new, `toolbox/gui/`)

`hFig = view_eigenmodes(SurfaceFile, iMode)` (iMode default 1):

1. `[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile)`; if `~isComputed` →
   `bst_error('Run "Compute eigenmodes" first.')` and return `[]`.
2. Open a cortex surface figure with the source-overlay pipeline (mechanism per
   the leadfield template; the exact bootstrap call — `view_surface` +
   manual overlay vs. an empty-overlay `view_surface_data` — is resolved in the
   plan, with the in-memory mode matrix that `process_eigenmodes_view` builds as
   the proven fallback substrate).
3. Render the current mode via the pure helper `GetModeDisplay` (below):
   write `.Data` into `TessInfo(1)`, set a **symmetric** `DataLimitValue`/
   `DataMinMax` `[-m, m]` (modes are signed +/− lobes), then
   `panel_surface('UpdateSurfaceColormap', hFig)`. Bipolar colormap.
4. Back up `KeyPressFcn`; install handler:
   - `leftarrow` → mode−1, `rightarrow` → mode+1 (clamped to `[1, K]`)
   - `pageup`/`pagedown` → ±10 (clamped)
   - `h` → help popup (mirrors leadfield viewer)
   - otherwise → call the backed-up callback (preserve rotate/zoom)
5. Bottom-left `uicontrol` legend via `GetModeDisplay().Label`:
   `Mode k / K   —   λ = <val>   —   wavelength ≈ 2π/√λ`.
   Figure name `Eigenmodes: <SurfaceFile>`; title hint
   `[←/→ next/prev mode, PgUp/PgDn ±10, H for help]`.
6. Transient — **no DB node**.

#### Pure helpers (headlessly testable)

- `GetModeDisplay(Eig, iMode)` → struct `.Data [nVert×1]`, `.CLim [-m m]`
  (`m = max(abs(Data))`, guarded against all-zero), `.Label` (formatted string),
  `.iMode` (validated). Used for both initial render and each step.
- `StepMode(iMode, delta, nModes)` → next index clamped to `[1, nModes]`.
- Wavelength: `2*pi / sqrt(lambda)` for `lambda > 0` (units follow the surface
  geometry); reported as "n/a" for `lambda <= 0`.

## Data flow

```
right-click cortex ─▶ Compute eigenmodes ─▶ java_dialog(nModes, massType)
                                          ─▶ process_eigenmodes('Compute', ...)
                                          ─▶ out_tess_eigenmodes ─▶ .Eigenmodes in surface .mat

right-click cortex ─▶ View eigenmodes ─▶ in_tess_eigenmodes(SurfaceFile)
                                       ─▶ surface figure @ mode 1 (GetModeDisplay)
                                       ─▶ ←/→ KeyPress ─▶ StepMode ─▶ GetModeDisplay
                                       ─▶ TessInfo.Data + UpdateSurfaceColormap
```

## Error handling

- View with no eigenmodes → `bst_error` directing to Compute.
- Compute on a non-manifold surface (repair off) → the core's `tess_manifold`
  check errors; `ComputeInteractive` shows it via `bst_error`.
- Dialog / overwrite-confirm cancel → silent no-op.
- Mode index clamps at `[1, K]` (no wraparound).
- All-zero mode (degenerate) → `GetModeDisplay` guards `CLim` against a
  zero-width range.

## Testing strategy

Repo idiom: `dev/tests/*.m` script-style functions printing `ALL TESTS PASSED`,
run via the MATLAB MCP (`evaluate_matlab_code` / `run_matlab_file`, **not**
`run_matlab_test_file`). GUI wrappers stay thin; logic lives in the pure helpers.

1. **Compute core round-trip** (`dev/tests/test_eigenmodes_compute_core.m`):
   build a small synthetic manifold surface (or reuse an existing test mesh),
   call the shared `Compute(...)`, then `out_`/`in_tess_eigenmodes`; assert mode
   count, ascending eigenvalues, and M-orthonormality. Confirms the `Run`
   refactor is behavior-preserving.
2. **Viewer helpers** (`dev/tests/test_view_eigenmodes_pure.m`): with a small
   `Eig` struct, assert `GetModeDisplay` returns symmetric `CLim`, correct
   `.Data` column, a non-empty `.Label`, and correct wavelength; assert
   `StepMode` clamps at 1 and K and steps ±1 / ±10 correctly.
3. **Interactive validation (user):** figure creation, arrow stepping, live
   colormap update, legend, and help popup — the parts that cannot be unit
   tested — are validated by the user clicking through in the GUI.

## Files touched

- `toolbox/tree/tree_callbacks.m` — two menu items (cortex node).
- `toolbox/process/functions/process_eigenmodes.m` — `ComputeInteractive` +
  extracted `Compute` core.
- `toolbox/gui/view_eigenmodes.m` — **new** viewer + pure helpers.
- `dev/tests/test_eigenmodes_compute_core.m` — **new**.
- `dev/tests/test_view_eigenmodes_pure.m` — **new**.
