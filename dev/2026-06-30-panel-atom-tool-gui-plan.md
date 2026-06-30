# Atom-tool panel GUI (Plan 2 of 2) Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints — this is Swing GUI work that needs the running Brainstorm session at each step; it is NOT headless-TDD-able and not suited to static subagents. The MATLAB-logic tasks (1-2) are TDD'd; the Swing tasks (3-5) are built + validated live by the controller in the GUI.

**Goal:** Replace `panel_bst_dynamics`'s dead Navigator strip with an "atom tool" — a filter selector + contextual per-kernel params + click-to-localize + live preview + Store — that authors a dynamics atom (thresholded eigfilter filter → Scout+Event) on the substrate shipped in Plan 1.

**Architecture:** The panel gains the eigfilter engine the designer already rides. A filter selector + `bst_eigfilter_controls('Sliders')`-driven contextual param rows set the kernel; a Localize toggle (the repurposed `bst_geodesic_tool` seed-picker) sets the seed; the realised field `W=bst_eigenfilter('Atom',ax,…)` previews through a new **atom-field mode** of `view_dynamics`'s overlay; Store calls `bst_dynamics('AtomFromKernel',…)` → `AddGroup` → `i_apply`. The standalone atom designer stays separate.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork (Java/Swing via `gui_component`). Engine: `bst_eigen('Axes')`, `bst_eigenfilter('Atom')`, `bst_eigfilter_controls`, `bst_dynamics('AtomFromKernel'/'Levelset'/'AddGroup')`, `view_dynamics`.

## Global Constraints

- No new dependencies. Live validation in the established Brainstorm session (preventad, sub-MTL0002 — do not clear). The designer is NOT touched (it stays the separate filter lab).
- Two design calls (approved): (1) preview = a new **atom-field mode** of `view_dynamics`'s `DynamicsOverlay` (paints a precomputed `W[nV×nT]` per frame); it does not disturb the existing Helmholtz `Op` modes. (2) Localize = the repurposed `bst_geodesic_tool` as a **seed-picker** (its `Toggle/OnClick/GetState` give the seed; its disk/`Grow` go unused) — its 3 external call sites (`figure_3d` ×2, `panel_surface` ×1) stay valid.
- The Frequency band-focus (`jFreqBand`/`OnFreqPreset` + the band↔bandpass display) is dropped this phase; frequency becomes a kernel param (a js kernel's `f0`).
- `lint` (MCP `check_matlab_code`) every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure (what changes)
- `toolbox/gui/view_dynamics.m` — add an atom-field branch to the overlay (`SetAtomField`/`ClearAtomField` verbs + `i_dynamics_overlay` paint).
- `toolbox/gui/panel_bst_dynamics.m` — replace the Navigator strip (`jNav`) with the Atom section; add `i_atom_preview`/`OnKernelChange`/`OnParamChange`/`OnLocalize`/`OnStore`; retire `i_axis_block`/`OnAxisChange`/`i_read_block`/`i_drive`/`jFreqBand`/`OnFreqPreset`/`jRegionTool`/`OnRegionTool`/`OnCaptureRegion`; update the `ctrl` struct + `i_cs`. `SyncSource` repurposed to set the atom seed.
- `dev/tests/` — `test_dynamics_atomfield_overlay.m` (Task 1), `test_atom_preview_realise.m` (Task 2). Swing tasks validated live.

---

### Task 1: Atom-field mode in the dynamics overlay

**Files:** Modify `toolbox/gui/view_dynamics.m`; Test `dev/tests/test_dynamics_atomfield_overlay.m`

**Interfaces:**
- Produces: `view_dynamics('SetAtomField', hFig, W, gv)` — stores `W[nV×nT]` (vertex×time) + global vertices `gv` into the `DynamicsOverlay` appdata (`D.AtomField=W; D.AtomGV=gv; D.Op='atom'`), then fires a repaint. `view_dynamics('ClearAtomField', hFig)` — clears it (`D.Op='none'`, `D.AtomField=[]`) and restores the native paint. The frame-paint helper `i_atom_frame_scalar(W, iT)` returns the per-vertex scalar for frame `iT` with the designer's density normalization (`i_normalize`-equivalent: one-signed → unit-mass density; signed → peak).
- Consumes: the existing `DynamicsOverlay` struct (`D.Cache/Op/srcDS/srcResult/iTess/nV`).

- [ ] **Step 1: failing test** — `dev/tests/test_dynamics_atomfield_overlay.m`: build a synthetic `W=[nV×nT]` with a localized Gaussian bump moving in time; assert `i_atom_frame_scalar(W, iT)` is the right length (`nV`), finite, and that for a one-signed `W` the scalar integrates to ~1 (unit mass) — and for a signed `W` it falls back to peak-normalized (max |.| == 1). (These are pure-function asserts; the figure paint is validated live.)
- [ ] **Step 2: run → fail** (function undefined).
- [ ] **Step 3: implement** — add the `SetAtomField`/`ClearAtomField` verbs (set `D` fields + `panel_surface('UpdateSurfaceData', hFig)` to trigger the overlay), and in `i_dynamics_overlay` add a branch: if `D.Op=='atom'` and `~isempty(D.AtomField)`, get `iT` (the same `GetTimeVector … CurrentTimeIndex` path), `scal = i_atom_frame_scalar(D.AtomField, iT)`, write `TessInfo(D.iTess).Data=scal` with the matching colormap (`source` sequential for density / `stat2` for signed), `UpdateSurfaceColormap`. Add `i_atom_frame_scalar` + the `i_normalize`-equivalent local.
- [ ] **Step 4: run → pass.**
- [ ] **Step 5: lint + commit** `feat(dynamics): atom-field mode in view_dynamics overlay (preview a realised W[nV×nT])`.

### Task 2: Panel realise/preview helper

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; Test `dev/tests/test_atom_preview_realise.m`

**Interfaces:**
- Produces: `i_atom_realise(st, kernelName, vals, seed)` — builds `ax = bst_eigen('Axes', OPTIONS)` from `st` (the linked surface + the recording's time window/rate), maps `kp = bst_eigfilter_controls('ToKernel', kernelName, vals, lmax)`, realises `[W,gv] = bst_eigenfilter('Atom', ax, kernelName, kp, seed)`, returns `(W, gv, ax)`. `i_atom_preview(st)` reads the current kernel + contextual slider values + seed from `ctrl`, calls `i_atom_realise`, and `view_dynamics('SetAtomField', st.hFig, W, gv)`.
- Consumes: Task 1's `SetAtomField`; `bst_eigen('Axes')`, `bst_eigenfilter('Atom')`, `bst_eigfilter_controls`.

- [ ] **Step 1: failing test** — `dev/tests/test_atom_preview_realise.m`: with a synthetic-axes stand-in (reuse the Task-1/substrate synthetic `ax` builder), call a thin pure form `i_atom_realise_core(ax, kernelName, kp, seed)` (the part that calls `bst_eigenfilter('Atom')`) and assert `W` is `[nV×nT]`, finite, peaks near `seed`. (The `bst_eigen('Axes')` + GUI read are validated live.)
- [ ] **Step 2: run → fail.**
- [ ] **Step 3: implement** `i_atom_realise`/`i_atom_realise_core`/`i_atom_preview` in `panel_bst_dynamics.m`.
- [ ] **Step 4: run → pass.**
- [ ] **Step 5: lint + commit** `feat(dynamics): panel realise/preview helper (eigfilter Atom -> SetAtomField)`.

### Task 3: The Atom section (Swing layout) — replaces the Navigator strip

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m` (`CreatePanel` lines ~133-150 + `ctrl` struct). **Live-validated.**

Build in place of `jNav` an "Atom" titled panel containing:
- **Filter selector** — `gui_component('combobox', …)` of kernel display names grouped (static / dynamic-ts / dynamic-js), from `bst_eigfilter_kernel('list')` + each kernel's `meta.display`. Callback → `OnKernelChange`.
- **Contextual param sub-panel** — a `JPanel` (`jAtomParams`) rebuilt by `i_build_param_rows(kernelName)` from `bst_eigfilter_controls('Sliders', kernelName, bounds)`: one labeled slider+readout row per non-empty spec row (reuse the existing slider/field idiom; skip empty-label rows). Callbacks → `OnParamChange`.
- **Localize toggle** (`jLocalize`) — replaces `jRegionTool`; callback → `OnLocalize`.
- **Threshold field** (`jThresh`, default 0.5) + **Store button** (`jStore`) — callbacks → `OnStore`.
- Update the `ctrl` struct: drop `jFreqC/jFreqW/jFreqBand/jSrcC/jSrcW/jRegionTool/jScaleC/jScaleW`; add `jKernel/jAtomParams/jLocalize/jThresh` (+ a handle list for the live param rows).

Live checkpoint: panel opens, the Atom section renders, switching the filter rebuilds the param rows. Commit `feat(gui): atom section replaces the Navigator strip in panel_bst_dynamics`.

### Task 4: Wire the interactions (kernel/param/localize/store)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`. **Live-validated.**

- `OnKernelChange` → `i_build_param_rows` + reset seed cursor + (if seed set) `i_atom_preview`.
- `OnParamChange` → `i_atom_preview` (guarded by `i_driving`).
- `OnLocalize` → `bst_geodesic_tool('Toggle', state)` (seed-pick mode); on OFF, clear the seed. `SyncSource` (already called by the tool's `Draw`) repurposed: read `GetState().seed/pos` into `st.atomSeed` (drop the disk radius), then `i_atom_preview`.
- `OnStore` → `ax`+`kp` from the current controls; `G = bst_dynamics('AtomFromKernel', ax, kernelName, kp, st.atomSeed, thr)`; `st.T = bst_dynamics('AddGroup', st.T, G)`; `i_apply(st)` + autosave. The stored atom carries its generator (Plan-1 fields).

Live checkpoint: pick a kernel → set params → toggle Localize → click cortex → see the preview update → set threshold → Store → the atom appears in the tree + a Scout+Event marker on the cortex; the `heat` filter reproduces the old geodesic-disk Scout. Commit `feat(gui): wire atom tool — kernel/param/localize/store via AtomFromKernel`.

### Task 5: Retire the dead Navigator code + reconcile

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; update `dev/tests/`. **Live + smoke.**

- Remove `i_axis_block`, `OnAxisChange`, `i_read_block`, `i_drive` (the dead axis nav), `jFreqBand`/`OnFreqPreset` (band-focus dropped), `OnRegionTool`, and `OnCaptureRegion` (superseded by Store). Reconcile `NotifySelection`/`i_sync_freq`/`i_driving` (the freq-focus echo guard) — keep only what still has a consumer.
- Confirm no orphan refs (`grep` for each removed function). Update any panel tests that referenced the removed handles (`test_nav_panel` etc.) — retarget or retire.
- Live smoke: open `view_dynamics` on a Dirac source result; author an atom end-to-end; confirm existing atom-table browse/markers still work and no console errors.

Commit `refactor(gui): retire dead Navigator nav code, atom tool is the sole authoring path`.

---

## Done criteria
- The panel's Navigator strip is replaced by the Atom section; selecting a filter shows its contextual params; click-to-localize sets a seed; the realised atom previews live and time-synced; Threshold + Store writes a Scout+Event atom carrying its generator; the `heat` filter reproduces the geodesic-disk Scout.
- Tasks 1-2 tests green (headless); Tasks 3-5 validated live in the GUI; no orphaned references; designer untouched.

## Risks / notes
- Swing layout is built live; the param-row rebuild on kernel change is the main new GUI mechanic.
- `bst_geodesic_tool` is repurposed (seed only), not retired — its `figure_3d`/`panel_surface` call sites stay valid; only the panel's use changes.
- The atom-field overlay mode coexists with the Helmholtz `Op` modes (Measure menu) — Store/preview use atom mode; the differential maps stay available.
