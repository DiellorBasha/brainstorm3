# Atom designer docked-panel Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints — Swing GUI + a live designer-launch gate (the Plan-1 regression proved headless tests miss launch bugs). Task 1 (the panel) is partly headless-testable; Task 2 (the rewire) is validated by launching the designer.

**Goal:** Lift the atom designer's "Filter design" controls off the figure_3d into a standalone Brainstorm panel docked in the tools tab while a design session is open.

**Architecture:** A new `panel_atom_designer` (Java `BstPanel`, verb-dispatched) hosts Operator/Filter/contextual-params/Connectome/Save/status, reusing `panel_eigenfilter_design`'s atom-slider renderer (driven by `bst_eigfilter_controls`). `view_atom_designer` stays the launcher/brain: it opens a cortex-only figure (+ the `WaveletDesignerPick` seed hook), builds its eigen-axes/state, creates+docks the panel passing its callback handles, and reads slider params back from the panel. The realise/normalize/overlay/fiber/save logic is unchanged.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork (Java/Swing via `gui_component`, `gui_show('BrainstormTab','tools')`). Engine reuse: `bst_eigfilter_controls`, `panel_eigenfilter_design` atom sliders, `bst_eigenfilter('Atom')`, `bst_dynamics`.

## Global Constraints

- No new dependencies. Live validation in the established Brainstorm session (preventad, sub-MTL0002 — do not clear). **A live launch of "Design atom…" is a hard gate** for Task 2.
- The figure_3d becomes cortex-only (surface + the `WaveletDesignerPick` hook). All design controls live in the docked panel.
- Reuse the shared system: kernels from `panel_eigenfilter_design('AtomKernels')`; contextual sliders via `BuildAtomSliders`/`ReadAtomVals`; params via `bst_eigfilter_controls('ToKernel')`. No new per-kernel control logic.
- Designer logic stays in `view_atom_designer` (handles bridge); only the controls + the param read-path move. Overlay unchanged (`ImageGridAmp` full-field scrub).
- `lint` every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure
- Create: `toolbox/gui/panel_atom_designer.m` — the docked panel (controls + thin accessor verbs; no atom math).
- Modify: `toolbox/gui/view_atom_designer.m` — remove the on-figure `uipanel` controls; create+dock the panel; rewire `Operator/Kernel/Param/Connectome` reads to the panel; lifecycle.
- Test: `dev/tests/test_panel_atom_designer.m`.

---

### Task 1: `panel_atom_designer` — the docked controls panel

**Files:** Create `toolbox/gui/panel_atom_designer.m`; Test `dev/tests/test_panel_atom_designer.m`

**Interfaces (produces):**
- `bstPanel = panel_atom_designer('CreatePanel')` — builds the Swing panel (Operator combobox, Filter combobox = `AtomKernels`, `jParams` sub-panel, Connectome toggle, Save button, status label), returns `BstPanel('AtomDesigner', …)`.
- `panel_atom_designer('Configure', cb, bounds, kernel0, variant0)` — store the designer callbacks `cb` (appdata `AtomDesignerCB`), select the Filter=`kernel0` / Operator from `variant0`, and `BuildAtomSliders` for `kernel0`+`bounds`.
- `vals = panel_atom_designer('ReadVals')` → `panel_eigenfilter_design('ReadAtomVals', jParams)` (the `[s1 s2 s3]`).
- `k = panel_atom_designer('CurrentKernel')` ; `v = panel_atom_designer('CurrentOperator')` (returns `'LB-Connectome'`|`'Laplace-Beltrami'`).
- `panel_atom_designer('RebuildSliders', kernel, bounds)` ; `panel_atom_designer('SetStatus', text)`.
- Swing callbacks (internal): Filter→`BuildAtomSliders`+`cb.Kernel()`; Operator→`cb.Operator()`; slider settle→`cb.Param()`; Connectome→`cb.Fibers(state)`; Save→`cb.Save()`.

- [ ] **Step 1: failing test** — `dev/tests/test_panel_atom_designer.m`:

```matlab
% test_panel_atom_designer - the panel builds, configures, and reads back its controls
bstPanel = panel_atom_designer('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel returns a BstPanel');
b = struct('scaleMinMM',7,'scaleMaxMM',95,'rateMinMM2',49,'rateMaxMM2',9025);
cb = struct('Kernel',@()[], 'Operator',@()[], 'Param',@()[], 'Fibers',@(s)[], 'Save',@()[]);
panel_atom_designer('Configure', cb, b, 'gabor', 'Laplace-Beltrami');
assert(strcmp(panel_atom_designer('CurrentKernel'),'gabor'), 'Filter set to gabor');
assert(strcmp(panel_atom_designer('CurrentOperator'),'Laplace-Beltrami'), 'Operator = geometric');
panel_atom_designer('RebuildSliders', 'resonator', b);   % switch + read back
% (CurrentKernel reflects the combobox, not the rebuilt sliders; values read via ReadVals)
v = panel_atom_designer('ReadVals');
assert(numel(v)==3, 'ReadVals returns [s1 s2 s3]');
disp('OK');
```

- [ ] **Step 2: run → fail** (function undefined).
- [ ] **Step 3: implement** `panel_atom_designer.m` (verb-dispatched). `CreatePanel` mirrors `panel_bst_dynamics` layout idioms (`gui_component` MenuBar/rows; `BstPanel(panelName, jPanelNew, ctrl)`); store the Java handles + `atomKeys` in `ctrl`. Internal Swing callbacks read `cb = getappdata(0,'AtomDesignerCB')`. `CurrentOperator` maps combobox idx → `{'LB-Connectome','Laplace-Beltrami'}`.
- [ ] **Step 4: run → pass** (`OK`).
- [ ] **Step 5: lint + commit** `feat(gui): panel_atom_designer — docked controls panel for the atom designer`.

---

### Task 2: rewire `view_atom_designer` to use the docked panel

**Files:** Modify `toolbox/gui/view_atom_designer.m`. **Live-validated (launch the designer).**

**Consumes:** Task 1's `panel_atom_designer` verbs.

- [ ] **Step 1: remove the on-figure controls.** Delete the `hP = uipanel(...)` block (Operator/Kernel/Scale/Speed/Decay/Connectome/Save/status `uicontrol`s, the `i_slider` helper, and `i_config_sliders`). The figure keeps only the cortex + the `WaveletDesignerPick` hook + `SyncControls`'s non-control work.

- [ ] **Step 2: create + dock the panel + bridge.** After the figure opens and `i_build_basis` has set `lmax`/bounds, build the callbacks struct and dock the panel:

```matlab
cb = struct('Kernel',   @()bst_call(@KernelChanged), ...
            'Operator', @()bst_call(@OperatorChanged), ...
            'Param',    @()bst_call(@ParamChanged), ...
            'Fibers',   @(s)bst_call(@()OnToggleConnectome(s)), ...
            'Save',     @()bst_call(@SaveAtom));
try, gui_hide('AtomDesigner'); catch, end %#ok<CTCH>
bstPanel = panel_atom_designer('CreatePanel');
gui_show(bstPanel, 'BrainstormTab', 'tools');
b = struct('scaleMinMM',scaleMinMM,'scaleMaxMM',scaleMaxMM,'rateMinMM2',rateMinMM2,'rateMaxMM2',rateMaxMM2);
panel_atom_designer('Configure', cb, b, kernel, variant);
```

- [ ] **Step 3: rewire the callbacks to read from the panel.** `i_phys2kernel` reads the panel sliders; `KernelChanged`/`OperatorChanged`/`OnToggleConnectome` read the panel state (no `src`):

```matlab
function kp = i_phys2kernel()
    vals = panel_atom_designer('ReadVals');
    kp = bst_eigfilter_controls('ToKernel', kernel, vals, lmax);
end
function KernelChanged()
    kernel = panel_atom_designer('CurrentKernel');
    i_reset_time();  Generate();
end
function OperatorChanged()
    variant = panel_atom_designer('CurrentOperator');
    i_build_basis(variant);
    b = struct('scaleMinMM',scaleMinMM,'scaleMaxMM',scaleMaxMM,'rateMinMM2',rateMinMM2,'rateMaxMM2',rateMaxMM2);
    panel_atom_designer('RebuildSliders', kernel, b);
    Generate();
end
function ParamChanged()
    Generate();
end
function OnToggleConnectome(state)
    showFib = logical(state);
    ... (existing fiber show/hide body, using showFib)
end
```
Replace status writes `set(hLabel,'String',txt)` with `panel_atom_designer('SetStatus', txt)`.

- [ ] **Step 4: lifecycle.** In `OnClose`, hide the panel before the existing teardown: `try, gui_hide('AtomDesigner'); catch, end`. Keep the working-results-file cleanup + `CloseRequestFcn` chain. (The figure-close already routes through `OnClose`.)

- [ ] **Step 5: live launch gate.** In the GUI: right-click a cortex → "Design atom…". Verify: the figure opens **cortex-only** (no on-figure Filter panel); the **AtomDesigner** panel docks in tools with Operator/Filter/sliders/Connectome/Save/status; clicking a vertex previews an atom; switching Operator/Filter rebuilds axes/sliders + repreviews; moving a slider repreviews; Connectome toggles fibers; Save writes a Scout+Event; closing the figure removes the panel. Fix anything that misbehaves.

- [ ] **Step 6: commit** `refactor(gui): atom designer controls move to the docked panel_atom_designer (figure is cortex-only)`.

---

## Done criteria
- Launching "Design atom…" opens a cortex-only figure + the docked AtomDesigner panel; all controls work from the panel; the figure no longer carries the Filter-design `uipanel`. `test_panel_atom_designer` passes; the live launch gate passes.

## Risks / notes
- The bridge: the panel's Swing callbacks call the designer's nested-function handles (closures, valid while the session is open). `cb` lives in appdata `AtomDesignerCB`; a second launch rebuilds the panel + overwrites `cb`.
- Singleton: one design session at a time (the panel is a single 'AtomDesigner' tab) — consistent with the current single-figure designer.
- `i_reset_time` is the designer's existing time-cursor reset (used by `KernelChanged`/`i_seed`); reuse it.
