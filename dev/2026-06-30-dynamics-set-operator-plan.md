# Dynamics atom Set-operator + recording-coupled axes Implementation Plan

> **For agentic workers:** Execute INLINE (superpowers:executing-plans) with live checkpoints. Task 1 (data + axes model) is headless-TDD'd; Task 2 (the Atoms submenu + realiser verification) is built + validated live.

**Goal:** Give each Dynamics atom a per-atom operator (eigenbasis), defaulting to the launch source's operator (Dirac when Dirac-launched, else Laplace-Beltrami), selectable via an Atoms → "Set operator" submenu (Geometric/Connectomic/Tangent), with the atom realised on that basis over a 4 s window at the recording's sample rate.

**Architecture:** `G.Operator` (a `bst_eigen` Variant string) on each atom generator; `i_atom_ensure_axes` reworked into `i_atom_axes(st, variant)` with a per-variant cache + recording `Fs` + 4 s window; `i_atom_realise` uses the selected atom's operator and reduces vector/complex bases (Dirac/Connection-Laplacian) to magnitude; a static radio submenu sets the selected atom's operator.

**Tech Stack:** MATLAB R2023b, Brainstorm dev fork. `bst_eigen('Axes')`, `bst_eigenfilter('Atom')`, `bst_memory('GetTimeVector')`, `gui_component` radio menu items + `ButtonGroup`.

## Global Constraints

- No new dependencies; Brainstorm components only. Live validation in the Brainstorm session (restart `brainstorm` if it drops).
- Operator variants: Geometric=`Laplace-Beltrami`, Connectomic=`LB-Connectome`, Tangent=`Connection Laplacian`, Dirac=`Dirac` (launch default). The submenu has all four radio items; selecting one sets the **selected atom's** `G.Operator`.
- Launch default: `Dirac` if the source result behind the figure is a Dirac inverse (its `Comment` contains `dirac`, case-insensitive), else `Laplace-Beltrami`.
- Axes: 4 s window at the recording `Fs` (`nF = round(4*Fs)`, `Fs` from the source result's time vector; fallback 100 Hz); cached per variant on `st.atomAxMap`.
- Vector/complex operators → preview paints magnitude `|W|`. If a variant fails to realise, guard (warn + skip preview), never crash.
- `lint` every edited `.m`; commit after each task with `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## File structure
- Modify: `toolbox/gui/panel_bst_dynamics.m`.
- Test: `dev/tests/test_dynamics_operator.m` (Task 1).

---

### Task 1: per-atom operator + recording-coupled per-variant axes (headless)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`; Test `dev/tests/test_dynamics_operator.m`

**Interfaces (produces):**
- `i_default_atom(kernelName, kp, seed, surfaceFile, label, operator)` — now also sets `G.Operator = operator`.
- `op = i_launch_operator(st)` — `'Dirac'` if the source result's Comment contains `dirac` (case-insensitive), else `'Laplace-Beltrami'`.
- `ax = i_atom_axes(st, variant)` — find-or-build + cache (`st.atomAxMap(variant)`) the eigen-axes for `variant`, over a 4 s window at the recording `Fs`.
- `i_atom_detail(G)` — now prepends the operator: `'<op> | <kernel> . vtx N . …'`.

- [ ] **Step 1: failing test** — `dev/tests/test_dynamics_operator.m`:

```matlab
% test_dynamics_operator - per-atom operator + launch-derived default
kp = struct('lmax',40,'tau',0.3,'vals',[400 0 0]);
G  = panel_bst_dynamics('i_default_atom', 'diffusion', kp, 13, 's.mat', 'atom1', 'Dirac');
assert(strcmp(G.Operator,'Dirac'), 'atom carries its operator');
s = panel_bst_dynamics('i_atom_detail', G);  assert(contains(s,'Dirac'), 'detail shows operator');
% launch-derived default from the source-result comment
stD = struct('srcComment','MN: shared dirac kernel');  assert(strcmp(panel_bst_dynamics('i_launch_operator', stD), 'Dirac'), 'dirac comment -> Dirac');
stL = struct('srcComment','MN: 2018 (Constr) 2018');   assert(strcmp(panel_bst_dynamics('i_launch_operator', stL), 'Laplace-Beltrami'), 'else -> LBO');
disp('OK');
```
(`i_launch_operator` reads `st.srcComment` when present — a test seam; live it reads the result behind `st.hFig`.)

- [ ] **Step 2: run → fail.**
- [ ] **Step 3: implement.**
  - `i_default_atom`: add the `operator` arg → `G.Operator = operator;` (default `'Laplace-Beltrami'` if omitted).
  - `i_atom_detail`: prepend `G.Operator` (if set).
  - `i_launch_operator(st)`: if `isfield(st,'srcComment')`, `cmt = st.srcComment`; else read the source result behind `st.hFig` (`D = getappdata(st.hFig,'DynamicsOverlay'); R = GlobalData...Results(D.srcResult); cmt = R.Comment`). Return `'Dirac'` if `contains(lower(cmt),'dirac')`, else `'Laplace-Beltrami'`.
  - `i_atom_axes(st, variant)`: ensure `st.atomAxMap` (a `containers.Map`); if cached return it; else get `Fs` from the recording (`[tv] = bst_memory('GetTimeVector', srcDS, srcResult)` → `Fs = 1/median(diff(tv))`; fallback 100), `nF = round(4*Fs)`, `ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant',variant,'nModes',60,'TimeWindow',[0 (nF-1)/Fs],'SampleRate',Fs))`; cache + return. (Retain the old `i_atom_ensure_axes` as a thin wrapper that calls `i_atom_axes(st, defaultVariant)` for back-compat, or replace its callers.)
- [ ] **Step 4: run → pass.**
- [ ] **Step 5: lint + commit** `feat(dynamics): per-atom operator + launch default + per-variant recording-coupled axes`.

---

### Task 2: Set-operator submenu + realise-on-operator + verification (live)

**Files:** Modify `toolbox/gui/panel_bst_dynamics.m`. **Live-validated.**

**Consumes:** Task 1's `i_default_atom`/`i_launch_operator`/`i_atom_axes`.

- [ ] **Step 1: the submenu.** In `CreatePanel`'s Atoms menu (after "Set color", ~L64), add a static radio submenu:

```matlab
jMenuOp = gui_component('Menu', jMenuAtoms, [], 'Set operator', IconLoader.ICON_PROPERTIES, [], []);
bgOp = javax.swing.ButtonGroup();
opDefs = {'Geometric','Laplace-Beltrami'; 'Connectomic','LB-Connectome'; 'Tangent (connection Laplacian)','Connection Laplacian'; 'Dirac','Dirac'};
jOpItems = javaArray('javax.swing.JRadioButtonMenuItem', size(opDefs,1));
for io = 1:size(opDefs,1)
    jit = gui_component('radio', jMenuOp, [], opDefs{io,1}, [], [], @(h,e)bst_call(@()OnSetOperator(opDefs{io,2})));
    bgOp.add(jit);  jOpItems(io) = jit;
end
```
Add `jMenuOp`/`jOpItems` + the variant list to the `ctrl` struct (`'jOpItems',jOpItems, 'opVariants',{opDefs(:,2)'}`).

- [ ] **Step 2: realise on the atom's operator + magnitude.** Rework `i_atom_realise(st, kernel, vals, seed, variant)` to take the variant, use `ax = i_atom_axes(st, variant)`, and after `bst_eigenfilter('Atom')` reduce non-scalar output: `if ~isreal(W), W = abs(W); end` and if `size(W,1) ~= ax.nVfull` reshape/agg to per-vertex magnitude (verified live). `i_atom_preview` reads the selected atom's `G.Operator` and passes it. `i_atom_writeback` keeps `G.Operator` unchanged.

- [ ] **Step 3: `OnSetOperator`.** Set the selected atom's operator, rebuild, preview:

```matlab
function OnSetOperator(variant) %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    ia = i_field(st,'curAtom',0);  if (ia<1)||(ia>numel(st.T.Groups)), return; end
    st.T.Groups(ia).Operator = variant;  setappdata(0,'DynamicsTarget',st);
    ctrl.jAtomInfo.setText(i_atom_detail(st.T.Groups(ia)));
    i_select_op_radio(variant);  i_atom_preview();
end
```

- [ ] **Step 4: defaults + check sync.** `OnCreateAtom` sets the new atom's operator to `i_launch_operator(st)` (pass it to `i_default_atom`). `i_select_atom_load` calls `i_select_op_radio(G.Operator)` (check the matching radio; `i_select_op_radio(v)` finds `v` in `ctrl.opVariants` and `jOpItems(k).setSelected(1)`).

- [ ] **Step 5: realiser verification (the risk).** Live: in a Dirac-launched session, realise an atom on each operator. Confirm Geometric/Connectomic produce a scalar field; confirm **Dirac/Connection-Laplacian** either realise into a paintable magnitude or are guarded (try/catch in `i_atom_realise` → on failure, `bst_error`/skip preview + leave the readout noting "operator not realisable"). Adjust the magnitude reduction to the actual output shape.

- [ ] **Step 6: live gate.** Dirac-launched session → new atom defaults to **Dirac** (readout shows it); Atoms → Set operator → Geometric/Connectomic/Tangent switches the basis + re-previews (magnitude for Tangent); the readout + the checked radio track the selection. Fix anything off.

- [ ] **Step 7: commit** `feat(dynamics): Atoms > Set operator (Geometric/Connectomic/Tangent/Dirac) per atom, realise on its basis`.

---

## Done criteria
- Each atom carries its operator (default = launch source's); Atoms → Set operator switches the selected atom's basis and re-previews over a 4 s / recording-Fs window; vector/complex operators preview as magnitude or are guarded; `test_dynamics_operator` passes; the live gate passes.

## Risks / notes
- `bst_eigenfilter('Atom')` on `Dirac`/`Connection Laplacian` is unverified — Step 5 is the gate; guard rather than crash. Geometric/Connectomic must work regardless.
- A variant's eigen file may need building (slow first use) — wrap `i_atom_axes` in `bst_progress`. The Dirac file exists from the launch.
- `i_atom_axes` at high `Fs` builds a large time grid (4 s × Fs); fine for one preview.
