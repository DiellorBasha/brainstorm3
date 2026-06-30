# Atom designer — dock the controls into a temporary panel — design

**Date:** 2026-06-30
**Status:** design (approved direction; pending spec review → plan)
**Depends on:** the shared contextual-slider system — `bst_eigfilter_controls` (Sliders/ToKernel) + `panel_eigenfilter_design` atom-slider renderer (AtomKernels/BuildAtomSliders/ReadAtomVals/SetAtomVals), both shipped/committed on `development`.

---

## 1. Motivation

`view_atom_designer` currently builds its controls as MATLAB `uicontrol`s embedded on the figure_3d figure (a `uipanel` with Operator/Kernel dropdowns, Scale/Speed/Decay sliders, Connectome toggle, Save button, status text). This clutters the 3D view and predates the shared contextual-slider system. Two goals:
1. **Factor the controls off the figure** into a Brainstorm Java panel that is **temporarily docked** while a design session is open (created on launch, removed on close).
2. **Join the new developments** — the designer's params become the unified contextual-slider system (`bst_eigfilter_controls` rendered by `panel_eigenfilter_design`), so all dynamic/js kernels and the physical-units sliders appear, exactly as in the `panel_bst_dynamics` atom tool.

The figure_3d becomes **cortex-only** (surface + the `WaveletDesignerPick` seed hook); the design controls live in the docked panel.

## 2. Architecture

- **New `panel_atom_designer`** (`toolbox/gui/panel_atom_designer.m`) — a Brainstorm `BstPanel` (verb-dispatched via `macro_method`, like `panel_bst_dynamics`). `CreatePanel(callbacks, kernelKeys, bounds, currentKernel)` builds the controls and returns the `BstPanel`; it stores the design callback handles + the Java control handles in its `ctrl` struct.
- **`view_atom_designer` stays the launcher/brain.** It opens the figure (cortex + the seed hook), builds its eigen-axes / scale bounds / state, and — instead of building on-figure `uicontrol`s — creates and docks `panel_atom_designer`, passing its callback handles. All realise/normalize/overlay/fiber/save logic stays in `view_atom_designer` (closure-based, unchanged).

## 3. Controls hosted by the panel

The panel hosts what is now on the figure, rendered with `gui_component`:
- **Operator** — combobox `{connectomic, geometric}` → the designer's `OperatorChanged` (rebuilds the eigen-axes/bounds).
- **Filter** — combobox of all kernels (`panel_eigenfilter_design('AtomKernels')`) → `KernelChanged`.
- **Contextual params** — a `jParams` sub-panel filled by `panel_eigenfilter_design('BuildAtomSliders', jParams, kernel, bounds, onSettle)`; `onSettle` → the designer's `ParamChanged`. Read back via `ReadAtomVals` → `bst_eigfilter_controls('ToKernel')`.
- **Connectome** — toggle → `ToggleFibers`.
- **Save** — button → `SaveAtom`.
- **Status** — a label line (relocated from the figure's `hLabel`); the designer updates it via a panel verb `panel_atom_designer('SetStatus', text)`.

## 4. Panel ↔ designer communication

The logic stays in `view_atom_designer`; the panel relocates only the controls.
- On launch `view_atom_designer` passes a struct of handles `cb = struct('Operator',@OperatorChanged,'Kernel',@KernelChanged,'Param',@ParamChanged,'Fibers',@ToggleFibers,'Save',@SaveAtom)`; the panel's Swing callbacks call `cb.*`.
- The one read-path change: the designer's `i_phys2kernel` reads the **panel's Java sliders** — `vals = panel_atom_designer('ReadVals')` (which calls `panel_eigenfilter_design('ReadAtomVals', ctrl.jParams)`), then `bst_eigfilter_controls('ToKernel', kernel, vals, lmax)`. The old `hScale/hSpeed/hDecay` `uicontrol`s and `i_config_sliders` (on-figure rendering) are removed from `view_atom_designer`.
- The current kernel/operator come from the panel too (`panel_atom_designer('CurrentKernel')` / `'CurrentOperator'`).

## 5. Lifecycle

- Launch: after the figure opens, `view_atom_designer` builds `cb` + the kernel list + bounds, then `gui_show(panel_atom_designer('CreatePanel', cb, keys, bounds, k0), 'BrainstormTab', 'AtomDesigner', 'tools')` (docked in the tools area).
- Close: the figure's `CloseRequestFcn`/`OnClose` hides + clears the panel (`gui_hide('AtomDesigner')`), mirroring the `panel_bst_dynamics` teardown. Reopening rebuilds it.
- Only one design session at a time (the panel is a singleton tab); a second launch reuses/rebuilds it.

## 6. Preview overlay — unchanged

The designer keeps its own managed source overlay (`ImageGridAmp = W`, native time-scrub). It does NOT adopt the panel atom tool's per-frame `SetAtomField` — the full-field scrub suits the exploration lab and already works. Realise (`i_eval_atom` → `bst_eigenfilter('Atom')`), normalize (`i_normalize`), colormap, fibers, and Save (`bst_dynamics`-backed Scout+Event) are unchanged.

## 7. Scope

**In:** `panel_atom_designer` (new docked panel) + the `view_atom_designer` rewire (controls → panel; read-path → panel sliders; lifecycle).
**Out / unchanged:** the `panel_bst_dynamics` atom tool (separate); the designer's realise/normalize/overlay/fiber/save logic; the `WaveletDesignerPick` figure hook; the shared renderer + control spec (already built).

## 8. Testing

- Headless: `panel_atom_designer('CreatePanel', …)` returns a `BstPanel` whose `ctrl` has the expected handles (jOperator/jKernel/jParams/jFibers/…); `ReadVals` round-trips through the `panel_eigenfilter_design` atom sliders (covered by `test_atom_sliders_java`); `CurrentKernel`/`CurrentOperator` map the comboboxes correctly.
- Live (controller + user): launch "Design atom…" → the figure is cortex-only + the panel docks in tools; pick a vertex → atom previews; switch Operator/Filter → axes/sliders rebuild + preview; adjust params → preview updates; Connectome toggles fibers; Save writes a Scout+Event; closing the figure removes the panel.

## 9. Risks / notes

- Removing the on-figure `uicontrol`s + `i_config_sliders` from `view_atom_designer` is the main edit; the realise/save closures stay intact and keep reading `kernel`/`lmax` from the designer state, only `vals` now comes from the panel.
- `panel_atom_designer` is Java/Swing — built + validated live (the renderer itself is already headless-tested).
- The status line moves to the panel; the figure no longer shows the `hLabel` text.
- The `WaveletDesignerPick` hook stays on the figure (seed picking is a figure gesture); the panel reacts to a new seed via the designer's existing pick callback (which triggers a re-realise).
