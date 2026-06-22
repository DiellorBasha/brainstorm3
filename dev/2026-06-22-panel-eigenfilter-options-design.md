# Design — `panel_eigenfilter_options` + `panel_eigenfilter_design`

Date: 2026-06-22
Author: Diellor Basha (with Claude)
Status: Approved design — pending implementation plan

## Goal

Add a Brainstorm GUI options panel for the eigen-domain spatial filter
(`bst_eigen` `Method='filter'` → `bst_eigenfilter`), built in the same spirit as
`panel_timefreq_options` is for the Morlet time-frequency transform. The panel is the
foundation for a future responsive view of the filter's **impulse response** (the spatial
"wavelet"). This round delivers a live **spectral-response `h(λ)`** preview as the first
step toward that.

## Scope decisions (settled in brainstorming)

1. **No library changes.** The `toolbox/eigen/eigfilter/` kernel library and every
   `bst_eigfilter_*` function are kept exactly as-is (simple, working). The earlier
   consolidation idea (single `bst_eigen_*` registry, folding `evaluate`/`compose`) is
   explicitly **not** done now.
2. **Standalone panel now; process later.** Build `panel_eigenfilter_options` as a
   self-contained Brainstorm panel. Do **not** create `process_eigenfilter` yet — it
   would require the eigen_/anatomy pipeline prerequisites. The panel is shaped so a
   future `process_eigenfilter` attaches via an `'editpref'` option with no panel rewrite.
3. **Fix the panel-naming convention violation.** `bst_eigfilter_panel.m` breaks Brainstorm
   convention (panels are `panel_*`, never `bst_*…_panel`; and it is a shared *section*, not
   a full panel). Rename it to `panel_eigenfilter_design.m`.
4. **`CreatePanel(EigenFile)`** signature (Option B): the panel loads `Lambda` (and later
   `Phi`) directly from the eigen_ node.
5. **Display v1 = spectral response `h(λ)`** curve (gain vs mode index), live-updating as
   sliders move. The cortical impulse-response (the actual wavelet atom on the surface) is
   deferred to a later iteration.

## Why these names

- `panel_eigenfilter_options` mirrors `panel_timefreq_options` (the `_options` panel idiom).
- `panel_eigenfilter_design` reframes the kernel-selector widget as the **filter-design**
  component — the GUI counterpart of `morlet_design`. The design+preview lives in one
  reusable component embedded by both `panel_helmholtz` and `panel_eigenfilter_options`,
  so the future impulse-response view is written once and appears in both.
- The `bst_eigfilter_*` *library* files keep the `bst_` prefix (correct for compute/library
  functions). Only the misnamed GUI file is fixed.

## Architecture — 3 file changes, no library changes

```
toolbox/gui/panel_eigenfilter_design.m   RENAME of bst_eigfilter_panel.m (+ DrawResponse verb)
toolbox/gui/panel_helmholtz.m            EDIT: bst_eigfilter_panel(...) -> panel_eigenfilter_design(...)
toolbox/gui/panel_eigenfilter_options.m  NEW (the panel_timefreq_options analogue)
```

### Component 1 — `panel_eigenfilter_design.m` (renamed widget)

- Keeps its current `macro_method` verbs verbatim: `Kernels`, `CurrentKernel`,
  `BuildSliders`, `ParamNames`, `ReadParams`, and the internal helpers. Only the function
  name and the doc header change.
- **Adds one stateless verb**: `DrawResponse(hAxes, kernelName, params, Lambda)`.
  - Resolve the kernel handle: `g = bst_eigfilter_kernel(kernelName, params)`.
  - Evaluate the gain: `h = bst_eigfilter_evaluate(g, Lambda)`.
  - Plot `h` vs mode index `k = 1..K` into `hAxes` (eigenvalue `λ` available as context,
    e.g. axis tooltip / secondary label).
  - This is the reusable rendering core and the future home of the cortical
    impulse-response (Option B) draw verb.

### Component 2 — `panel_helmholtz.m`

- Mechanical rename of the ~6 `bst_eigfilter_panel('…')` call sites to
  `panel_eigenfilter_design('…')`. No behavior change. (Only existing caller besides the
  file's own self-references.)

### Component 3 — `panel_eigenfilter_options.m` (new)

Standard Brainstorm panel contract: `eval(macro_method)` dispatch, built with
`gui_river`/`gui_component`/`bst_mutex`, returning
`BstPanel('EigenfilterOptions', jScroll, ctrl)`.

- **`CreatePanel(EigenFile)`** — loads `EigenMat = in_bst_eigen(EigenFile)`; takes `Lambda`
  from the first non-empty hemisphere as the design axis. UI sections, top to bottom:
  1. **Eigen basis** — read-only label: Variant, mode count `K`, hemispheres present,
     source node comment.
  2. **Filter design** — embeds `panel_eigenfilter_design` (kernel dropdown + mode-index
     sliders) with `onSettle → UpdateResponse`.
  3. **Display** — a toggle (the `DisplayTimeResolution` analogue) that opens/closes a
     figure; `UpdateResponse` calls `panel_eigenfilter_design('DrawResponse', …)` and
     redraws live as sliders settle.
  4. **Comment** — output comment text field.
  5. **OK / Cancel** — `bst_mutex` release / `gui_hide`.
- **`GetPanelContents()`** — returns an OPTIONS struct that feeds `bst_eigen` directly:
  - `Method      = 'filter'`
  - `EigenFile   = <the node passed to CreatePanel>`
  - `KernelName  = panel_eigenfilter_design('CurrentKernel', …)`
  - `KernelParams= panel_eigenfilter_design('ReadParams', …)`
  - `Comment     = <comment text>`

### Arg-type-aware `CreatePanel` (future-process compatibility)

Brainstorm's `editpref` mechanism calls `CreatePanel(sProcess, sInputs)`, not
`CreatePanel(EigenFile)`. To keep Option B now and ease the future process:

- `CreatePanel` dispatches on the first argument's type:
  - `char` / file → the **EigenFile path** (implemented now).
  - `struct` (sProcess) + sInputs → reserved branch for the future
    `(sProcess, sInputs)` path that resolves the eigen node from inputs (stub/error now).
- The later `process_eigenfilter` then drops in by filling the struct branch — no rewrite.

## Data flow

```
caller → panel_eigenfilter_options('CreatePanel', EigenFile)
       → gui_show → [user designs filter, toggles live h(λ) preview] → OK
       → panel_eigenfilter_options('GetPanelContents')  → OPTIONS
       → bst_eigen(Data, OPTIONS)  → results_eigenfilter file
```

The panel configures only the filter; the input `Data` is passed to `bst_eigen` separately
by the caller. Clean separation: the panel needs only the eigenbasis (for sliders + preview).

## Error handling

- `CreatePanel`: empty/invalid `EigenFile`, or `EigenMat` missing `Lambda`/`Phi` →
  `bst_error` and return `[], []`.
- `GetPanelContents`: invalid params surface as errors (dog's `t1 < t2` is already enforced
  inside `ReadParams`).
- `DrawResponse`: guard empty `Lambda`. The sliders only ever produce single-scale params,
  so the filterbank/cell-array kernel case cannot occur here (no special handling needed).

## Validation

Brainstorm has no unit-test harness; validation is manual / script-based (run via the
MATLAB MCP against a loaded protocol):

1. Find a real `eigen_` node in the protocol.
2. `panel_eigenfilter_options('CreatePanel', EigenFile)` → `gui_show`.
3. `panel_eigenfilter_options('GetPanelContents')` → assert OPTIONS fields
   (`Method='filter'`, non-empty `EigenFile`, valid `KernelName`, struct `KernelParams`).
4. Optionally run `bst_eigen(Data, OPTIONS)` on a source file and confirm a
   `results_eigenfilter` file appears.
5. `checkcode` lint on all three files (Brainstorm idioms filtered).
6. Visual check: the `h(λ)` plot renders and updates live as sliders move.

## Out of scope (explicitly deferred)

- `process_eigenfilter` and its pipeline/anatomy prerequisites.
- Cortical impulse-response (wavelet-atom-on-surface) preview — the Option B Display.
- Any consolidation/renaming of the `bst_eigfilter_*` kernel library.
- Auto-resolving the eigen_ node from a surface (still a `TODO` in `bst_eigen`).
```
