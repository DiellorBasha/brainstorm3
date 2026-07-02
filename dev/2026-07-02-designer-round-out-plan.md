# SP2a — Round out the atom designer to the full operator range · Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the standalone `view_atom_designer` (+ `panel_atom_designer`) support the full SP1-routable operator range (Geometric / Connectomic / Dirac / Dirac-Connectome), rendering scalar atoms as density fields and vector (Dirac) atoms as magnitude + quivers, with a seed-direction control — all data-free.

**Architecture:** `bst_eigenfilter('Atom')` becomes the single fiber-aware realiser: it already computes `[C,kind]=i_fiber(ax)` internally, so we fold the quaternion/tangent → ambient-3-vector decode (currently stranded in `panel_bst_dynamics` `i_atom_realise_core`) INTO it, returning `[W, gv, V3, isSigned]` (backward-compatible — extra outputs only). The default seed direction moves to a new app-side shared helper `bst_atom_default_dir.m`. The designer branches on `isempty(V3)` (fiber), never on operator name: scalar → existing density/peak paint; vector → `|V3|` magnitude scalar + a `QuiverVectorOverride` quiver overlay (the same figure_3d source-vector idiom the Dynamics panel uses, which bypasses the `nComponents==3` requirement when the override is set).

**Tech Stack:** MATLAB R2023b, Brainstorm (dev fork), Swing GUI via `gui_component`, `bst_eigen`/`bst_eigenfilter` eigen stack, `figure_3d` source-vector quivers.

## Global Constraints

- **Branch:** all work on `feature/designer-round-out`, created off `development` **in-place (NOT a worktree)** — the live MATLAB session needs THIS checkout on its path for the controller live pass. The controller creates the branch before dispatching Task 1.
- **Execution policy — implementer subagents do STATIC CHECKS ONLY:** `mcp__plugin_brainstorm-dev_MATLAB__check_matlab_code` (read-only lint) + `grep` + write test `.m` files + commit. **They DO NOT run/eval MATLAB or launch Brainstorm** (it hung / wiped the live session in prior sessions).
- **The CONTROLLER runs every test.** Headless unit tests run in a separate process: `matlab -batch "addpath('dev/tests'); addpath(genpath('toolbox')); BST_TEST_SURF='...'; test_X"`. Live/GUI checks run via the MATLAB-MCP `evaluate_matlab_code` in the established session (protocol `preventad`, do NOT `clear`). `check_matlab_code` misses runtime errors — the controller pass is the real gate.
- **Byte-equivalence is a hard requirement:** the relocated decode must produce output bit-identical to the panel's current `i_atom_realise_core` for the scalar AND quaternion fibers (pinned by the Task 1 unit test) so `panel_bst_dynamics` (otherwise unchanged in SP2a) keeps behaving.
- **Tests:** `dev/tests/test_*.m`, assert-based. Surface-dependent tests read the surface from `getenv('BST_TEST_SURF')` and `assumeTrue`-skip when unset. **`bst_eigen('Axes', …)` REQUIRES `SampleRate` + `TimeWindow` fields** (BuildTimeAxis).
- **Test surface:** `sub-MTL0005/tess_cortex_pial_low.mat` in protocol `preventad` (the source `sub-MTL0002`'s `R.SurfaceFile`).
- **Data-free:** SP2a touches nothing about the inverse kernel / sensor data.
- **Commits:** conventional-commit prefixes; end each commit message with the `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` trailer (internal feature branch → development, so the trailer stays; this is NOT a clean upstream-PR branch).
- **Never `clear`** in the live session (wipes GlobalData); edited `.m` auto-reload; use `rehash`.

---

## File Structure

| File | Responsibility | Task |
|---|---|---|
| `toolbox/eigen/bst_eigenfilter.m` | `Atom` verb: single fiber-aware realiser, now returns `[W, gv, V3, isSigned]`. `Fiber` unchanged. | 1 |
| `toolbox/gui/bst_atom_default_dir.m` (new) | App-side default seed direction per fiber (relocated from the panel local). | 2 |
| `toolbox/gui/panel_bst_dynamics.m` | Two byte-equivalent touches: (a) call `bst_atom_default_dir`; (b) `i_atom_realise_core` consumes `Atom`'s `V3`. | 2, 3 |
| `toolbox/gui/panel_atom_designer.m` | 4-operator selector (combobox) + seed-direction control + accessors. | 4 |
| `toolbox/gui/view_atom_designer.m` | Fiber-driven render (scalar paint vs `|V3|` magnitude + quivers); default & interactive seed direction. | 5, 6 |
| `dev/tests/test_atom_fiber_decode.m` (new) | Task 1 unit test (byte-equivalence). | 1 |
| `dev/tests/test_atom_default_dir.m` (new) | Task 2 unit test. | 2 |
| `dev/tests/test_realise_core_delegates.m` (new) | Task 3 regression test. | 3 |
| `dev/tests/test_designer_operator_map.m` (new) | Task 4 pure display↔variant map test. | 4 |

**Dependency DAG:** `1 → 3`; `1 → 5`; `2 → 4 → 5 → 6`; `2 → 5`; `2 → 6`; `4 → 6`. Linear execution order 1 → 2 → 3 → 4 → 5 → 6 satisfies it.

---

### Task 1: `bst_eigenfilter('Atom')` returns the fiber-decoded 3-vector

**Files:**
- Modify: `toolbox/eigen/bst_eigenfilter.m:66` (signature) and `:79` and after `:116` (decode)
- Test: `dev/tests/test_atom_fiber_decode.m` (create)

**Interfaces:**
- Consumes: `bst_eigen('Axes', struct(...))` → `ax` (fields `Phi/Lambda/Mass/GlobalVertices/Operator/nT/tlag/omega/NFFT/SurfaceFile`); `bst_eigenfilter('Fiber', ax)` → `[C, kind]`; `manifold_quat_imag`.
- Produces: `[W, gv, V3, isSigned] = bst_eigenfilter('Atom', ax, KernelName, KernelParams, seedVert, seedDir)`. `W` = raw modal field `[nrows x nT]` (unchanged); `gv` = block global vertices (unchanged); `V3` = full-support ambient 3-vector `[nSrc x 3]` (`[]` for the scalar fiber; frame-1 imag-3-vec for quaternion; `a·e1+b·e2` for tangent); `isSigned` = `false` for the vector fibers, `[]` for scalar.

- [ ] **Step 1: Write the failing test** — `dev/tests/test_atom_fiber_decode.m`

```matlab
function test_atom_fiber_decode
% Byte-equivalence: bst_eigenfilter('Atom') V3 output == the panel's former i_atom_realise_core decode,
% for the scalar (Laplace-Beltrami) and quaternion (Dirac) fibers. V3=[] for scalar; imag-3-vec for Dirac.
    surf = getenv('BST_TEST_SURF');
    assert(~isempty(surf), 'Set BST_TEST_SURF to a cortex surface file (skips if unset).');

    kp = struct('lmax', []);
    % ---- scalar fiber ----
    axS  = i_axes(surf, 'Laplace-Beltrami');
    seedS = axS.GlobalVertices{1}(1);
    [Ws, gvs, V3s, sgnS] = bst_eigenfilter('Atom', axS, 'heat', kp, seedS, []);
    assert(isempty(V3s), 'scalar fiber must yield V3=[]');
    assert(isempty(sgnS), 'scalar fiber must yield isSigned=[]');
    [Ws2, gvs2] = bst_eigenfilter('Atom', axS, 'heat', kp, seedS, []);   % backward-compat 2-output
    assert(isequal(Ws, Ws2) && isequal(gvs, gvs2), 'scalar [W,gv] must be unchanged');

    % ---- quaternion fiber ----
    axQ  = i_axes(surf, 'Dirac');
    seedQ = axQ.GlobalVertices{1}(1);
    dir   = [1 0 0];
    [Wq, gvq, V3q, sgnQ] = bst_eigenfilter('Atom', axQ, 'heat', kp, seedQ, dir);
    V3ref = i_ref_quat_decode(Wq, gvq, axQ);                             % former panel decode, reproduced
    assert(isequal(V3q, V3ref), 'quaternion V3 must byte-match the former decode');
    assert(size(V3q,1) == i_nsrc(axQ) && size(V3q,2) == 3, 'V3 must be [nSrc x 3]');
    assert(sgnQ == false, 'vector fiber must yield isSigned=false');
    disp('test_atom_fiber_decode PASSED');
end

function ax = i_axes(surf, variant)
    ax = bst_eigen('Axes', struct('SurfaceFile',surf, 'Variant',variant, ...
                   'nModes',60, 'TimeWindow',[0 0.5], 'SampleRate',100));
end

function nV = i_nsrc(ax)
    nV = 0; for h=1:numel(ax.GlobalVertices)
        if ~isempty(ax.GlobalVertices{h}), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    end
end

function V3 = i_ref_quat_decode(W, gv, ax)
    % Verbatim reproduction of the former panel_bst_dynamics i_atom_realise_core quaternion branch.
    nV = i_nsrc(ax);  n = numel(gv);
    im = reshape(manifold_quat_imag(W(:,1)), 3, n).';
    V3 = zeros(nV,3);  V3(gv,:) = im;
end
```

- [ ] **Step 2: (implementer) static-check the test compiles**

Run `check_matlab_code` on `dev/tests/test_atom_fiber_decode.m`. Expected: no errors (Brainstorm-idiom warnings OK). Do NOT run MATLAB.

- [ ] **Step 3: Widen the `Atom` signature and use the fiber components**

In `toolbox/eigen/bst_eigenfilter.m`, change line 66:

```matlab
function [W, gv, V3, isSigned] = Atom(ax, KernelName, KernelParams, seedVert, seedDir) %#ok<DEFNU>
```

Change line 79 (remove the unused-suppression; `C` is now consumed):

```matlab
    [C, kind] = i_fiber(ax);
```

- [ ] **Step 4: Fold the fiber decode into `Atom`**

In `toolbox/eigen/bst_eigenfilter.m`, immediately AFTER the domain `switch` block ends (after the current line 116 `end`, before the function's closing `end` at line 117), insert:

```matlab
    % --- decode the realised atom for its fiber into an ambient 3-vector V3 (frame 1) ---
    % Relocated from panel_bst_dynamics i_atom_realise_core so that realising an atom yields its
    % physical vectors directly. Byte-equivalent to the former panel decode. Scalar fiber -> V3=[].
    V3 = [];  isSigned = [];
    switch kind
        case 'quaternion'
            n  = numel(gv);
            im = reshape(manifold_quat_imag(W(:,1)), 3, n).';       % [n x 3] imag 3-vector, frame 1
            V3 = zeros(nSrc, 3);  V3(gv,:) = im;
            isSigned = false;                                       % magnitude render is one-signed
        case 'tangent'
            Fr = ax.Operator.Frame{blk};                           % operator frame for the seed's block
            a  = real(W(:,1));  b = imag(W(:,1));
            V3 = zeros(nSrc, 3);  V3(gv,:) = a.*Fr.e1 + b.*Fr.e2;   % ambient a*e1 + b*e2
            isSigned = false;
    end
    if C == 1, V3 = []; isSigned = []; end                          % scalar guard
```

(`blk`, `gv`, `nSrc`, `C`, `kind` are all already in scope from lines 73–80.)

- [ ] **Step 5: (implementer) static-check the implementation**

Run `check_matlab_code` on `toolbox/eigen/bst_eigenfilter.m`. Expected: no new errors. Then `grep -n "i_atom_realise_core\|V3" toolbox/gui/panel_bst_dynamics.m` to confirm the panel still has its own decode (untouched until Task 3) — the two coexist because `Atom`'s extra outputs are backward-compatible.

- [ ] **Step 6: Commit**

```bash
git add toolbox/eigen/bst_eigenfilter.m dev/tests/test_atom_fiber_decode.m
git commit -m "feat(eigen): bst_eigenfilter('Atom') returns fiber-decoded V3 + isSigned

Fold the quaternion/tangent -> ambient 3-vector decode into the single Atom
realiser (it already resolves the fiber). Backward-compatible extra outputs.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate (not the implementer):** run
`matlab -batch "addpath(genpath('toolbox')); addpath('dev/tests'); setenv('BST_TEST_SURF', <surf>); test_atom_fiber_decode"`
Expected: `test_atom_fiber_decode PASSED`. This pins the byte-equivalence.

---

### Task 2: Relocate the default seed direction to `bst_atom_default_dir.m`

**Files:**
- Create: `toolbox/gui/bst_atom_default_dir.m`
- Modify: `toolbox/gui/panel_bst_dynamics.m` (delete local `i_atom_default_dir` at `:726–740`; rewire call sites `:704`, `:824`, `:1588`)
- Test: `dev/tests/test_atom_default_dir.m` (create)

**Interfaces:**
- Produces: `dir = bst_atom_default_dir(ax, seedVert)` → `1` for the scalar/tangent fiber; a unit row 3-vector (seed vertex normal) for the quaternion fiber.
- Consumes: `bst_eigenfilter('Fiber', ax)`, `in_tess_bst(ax.SurfaceFile, 0)`.

- [ ] **Step 1: Write the failing test** — `dev/tests/test_atom_default_dir.m`

```matlab
function test_atom_default_dir
% bst_atom_default_dir: scalar -> 1; quaternion -> unit seed normal (byte-equiv to the former panel local).
    surf = getenv('BST_TEST_SURF');
    assert(~isempty(surf), 'Set BST_TEST_SURF (skips if unset).');

    axS = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 0.5],'SampleRate',100));
    seedS = axS.GlobalVertices{1}(1);
    assert(isequal(bst_atom_default_dir(axS, seedS), 1), 'scalar default must be 1');

    axQ = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',40,'TimeWindow',[0 0.5],'SampleRate',100));
    seedQ = axQ.GlobalVertices{1}(1);
    d = bst_atom_default_dir(axQ, seedQ);
    assert(isequal(size(d),[1 3]) && abs(norm(d)-1) < 1e-9, 'quaternion default must be a unit row 3-vector');
    S = in_tess_bst(axQ.SurfaceFile, 0);
    nref = S.VertNormals(seedQ,:);  nref = nref / norm(nref);
    assert(max(abs(d - nref)) < 1e-9, 'quaternion default must equal the unit seed normal');
    disp('test_atom_default_dir PASSED');
end
```

- [ ] **Step 2: (implementer) static-check the test** — `check_matlab_code` on the test file. Expected: no errors.

- [ ] **Step 3: Create `toolbox/gui/bst_atom_default_dir.m`** (verbatim relocation)

```matlab
function dir = bst_atom_default_dir(ax, seedVert)
% BST_ATOM_DEFAULT_DIR  Default impulse direction for an atom's operator fiber (app-side).
%   dir = bst_atom_default_dir(ax, seedVert)
%   Quaternion (Dirac) fiber -> the seed vertex surface normal (unit row 3-vector).
%   Tangent/scalar fiber     -> 1 (frame e1 / scalar amplitude).
% App-side default: orchestrators (the library realisers) carry no defaults; the GUI supplies them.
% Shared by view_atom_designer and (SP2b) panel_bst_dynamics. See atom-operator-applicability.
%
% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
% =============================================================================@

    [~, kind] = bst_eigenfilter('Fiber', ax);
    if ~strcmp(kind, 'quaternion'), dir = 1; return; end
    dir = [0 0 1];                                       % fallback if normals are missing
    try
        S = in_tess_bst(ax.SurfaceFile, 0);
        if isfield(S,'VertNormals') && ~isempty(S.VertNormals) && seedVert <= size(S.VertNormals,1)
            n = S.VertNormals(seedVert, :);  if norm(n) > 0, dir = n / norm(n); end
        end
    catch %#ok<CTCH>
    end
    dir = dir(:)';
end
```

- [ ] **Step 4: Delete the panel local and rewire the three call sites**

In `toolbox/gui/panel_bst_dynamics.m`:

Delete the local function (current lines 726–740, the `function dir = i_atom_default_dir(ax, seedVert)` block and its leading comment on 725–726).

Replace each call `i_atom_default_dir(ax, seed)` with `bst_atom_default_dir(ax, seed)` at the three sites:
- line 704 (in `i_atom_realise`): `if (nargin < 6) || isempty(seedDir), seedDir = bst_atom_default_dir(ax, seed); end`
- line 824 (in `OnPickSeedDir`, `'Normal'` branch): `dir = bst_atom_default_dir(ax, seed);`
- line 1588 (in `OnCreateAtom`): `sdir = bst_atom_default_dir(ax, seed);`

- [ ] **Step 5: (implementer) static-check + verify no stragglers**

`check_matlab_code` on `panel_bst_dynamics.m` and `bst_atom_default_dir.m`. Then:
```bash
grep -n "i_atom_default_dir" toolbox/gui/panel_bst_dynamics.m
```
Expected: **no matches** (all three rewired, local deleted).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/bst_atom_default_dir.m toolbox/gui/panel_bst_dynamics.m dev/tests/test_atom_default_dir.m
git commit -m "refactor(gui): default seed direction -> shared bst_atom_default_dir

Relocate panel_bst_dynamics's local i_atom_default_dir into an app-side shared
helper (orchestrators carry no defaults); rewire the three panel call sites.
Byte-equivalent.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate:** run `matlab -batch "... test_atom_default_dir"` → `PASSED`. Then a quick live smoke: the Dynamics panel Design preview of a scalar and a Dirac atom still previews (no regression from the relocation).

---

### Task 3: `i_atom_realise_core` consumes `Atom`'s V3 (thin)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m:677–695` (`i_atom_realise_core`)
- Test: `dev/tests/test_realise_core_delegates.m` (create)

**Interfaces:**
- Consumes: `bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir)` → `[W, gv, V3, …]` (from Task 1).
- Produces: `[W, gv, V3] = panel_bst_dynamics('i_atom_realise_core', ax, kernel, kp, seed, seedDir)` — unchanged 3-output contract; the wrapper `i_atom_realise` (`:698`) and its callers are unaffected.

- [ ] **Step 1: Write the failing/regression test** — `dev/tests/test_realise_core_delegates.m`

```matlab
function test_realise_core_delegates
% i_atom_realise_core must now be a thin passthrough of bst_eigenfilter('Atom'): identical [W,gv,V3].
    surf = getenv('BST_TEST_SURF');
    assert(~isempty(surf), 'Set BST_TEST_SURF (skips if unset).');
    kp = struct('lmax', []);

    for variant = {'Laplace-Beltrami','Dirac'}
        ax   = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant',variant{1},'nModes',50,'TimeWindow',[0 0.5],'SampleRate',100));
        seed = ax.GlobalVertices{1}(1);
        dir  = bst_atom_default_dir(ax, seed);
        [Wc, gvc, V3c]   = panel_bst_dynamics('i_atom_realise_core', ax, 'heat', kp, seed, dir);
        [Wa, gva, V3a]   = bst_eigenfilter('Atom', ax, 'heat', kp, seed, dir);
        assert(isequal(Wc,Wa) && isequal(gvc,gva) && isequal(V3c,V3a), ...
            'realise_core must equal Atom for %s', variant{1});
    end
    disp('test_realise_core_delegates PASSED');
end
```

- [ ] **Step 2: (implementer) static-check the test** — `check_matlab_code`. Expected: no errors.

- [ ] **Step 3: Replace the decode body with a delegating call**

In `toolbox/gui/panel_bst_dynamics.m`, replace the whole `i_atom_realise_core` body (lines 677–695) with:

```matlab
% Pure realise-core: run the atom on ax at (seed,dir); return the raw field W [C*n x nT], its global
% vertices gv, and the decoded full-surface ambient vectors V3 [nV x 3] ([] for scalar). The fiber
% decode now lives in bst_eigenfilter('Atom') (SP2a); this is a thin passthrough. SP2b removes it.
function [W, gv, V3] = i_atom_realise_core(ax, kernel, kp, seed, seedDir) %#ok<DEFNU>
    [W, gv, V3] = bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir);
end
```

- [ ] **Step 4: (implementer) static-check + confirm the decode is gone**

`check_matlab_code` on `panel_bst_dynamics.m`. Then:
```bash
grep -n "manifold_quat_imag\|Operator.Frame" toolbox/gui/panel_bst_dynamics.m
```
Expected: no matches inside the former `i_atom_realise_core` region (the decode now lives only in `bst_eigenfilter.m`).

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/tests/test_realise_core_delegates.m
git commit -m "refactor(gui): i_atom_realise_core delegates the fiber decode to Atom

The quaternion/tangent -> V3 decode now lives in bst_eigenfilter('Atom');
i_atom_realise_core is a thin passthrough (SP2b removes it entirely).

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate:** run `matlab -batch "... test_realise_core_delegates"` → `PASSED`. Live: Dynamics Design preview of a scalar atom AND a Dirac atom (magnitude + quivers) unchanged from before SP2a.

---

### Task 4: `panel_atom_designer` — 4-operator selector + seed-direction control

**Files:**
- Modify: `toolbox/gui/panel_atom_designer.m` (operator widget `:59–65`; `ctrl` struct `:93–94`; `Configure` `:103`; `CurrentOperator` `:129–133`; API header `:10–18`)
- Test: `dev/tests/test_designer_operator_map.m` (create)

**Interfaces:**
- Produces (panel API): `panel_atom_designer('CurrentOperator')` → one of `'Laplace-Beltrami'|'LB-Connectome'|'Dirac'|'Dirac-Connectome'`; `panel_atom_designer('CurrentDirection')` → one of `'Normal'|'+X'|'+Y'|'+Z'|'Pick-on-surface'`; `panel_atom_designer('ShowDirection', tf)`; pure maps `panel_atom_designer('i_variant_for_item', name)` and `panel_atom_designer('i_item_for_variant', variant)`. The `cb` struct handed to `Configure` gains a `Direction` field.
- Consumes: `gui_component`, `gui_river`, `java_setcb`, `bst_get('PanelControls','AtomDesigner')`, `getappdata(0,'AtomDesignerCB')`.

- [ ] **Step 1: Write the failing test** — `dev/tests/test_designer_operator_map.m`

```matlab
function test_designer_operator_map
% Pure display<->variant maps for the 4-operator selector (no Swing needed).
    pairs = { 'Geometric','Laplace-Beltrami'; 'Connectomic','LB-Connectome'; ...
              'Dirac','Dirac'; 'Dirac (connectome)','Dirac-Connectome' };
    for i = 1:size(pairs,1)
        v = panel_atom_designer('i_variant_for_item', pairs{i,1});
        assert(strcmp(v, pairs{i,2}), 'item %s -> %s (got %s)', pairs{i,1}, pairs{i,2}, v);
        nm = panel_atom_designer('i_item_for_variant', pairs{i,2});
        assert(strcmp(nm, pairs{i,1}), 'variant %s -> %s (got %s)', pairs{i,2}, pairs{i,1}, nm);
    end
    % Unknown item defaults to Geometric (Laplace-Beltrami).
    assert(strcmp(panel_atom_designer('i_variant_for_item','???'), 'Laplace-Beltrami'));
    disp('test_designer_operator_map PASSED');
end
```

- [ ] **Step 2: (implementer) static-check the test** — `check_matlab_code`. Expected: no errors.

- [ ] **Step 3: Replace the 2 operator toggles with a 4-item combobox**

In `toolbox/gui/panel_atom_designer.m`, replace lines 59–65 with:

```matlab
        % Operator: combobox over the SP1-routable operators (Connection Laplacian is out — unpersisted frame)
        gui_component('label', jDes, 'br', 'Operator:');
        opItems = {'Geometric','Connectomic','Dirac','Dirac (connectome)'};
        jOperator = gui_component('combobox', jDes, 'tab hfill', [], {opItems}, [], [], []);
        java_setcb(jOperator, 'ActionPerformedCallback', @(h,e) bst_call(@OnOperatorCb));
```

- [ ] **Step 4: Add the seed-direction control widget**

In `toolbox/gui/panel_atom_designer.m`, immediately AFTER the `jParams` row is added to `jDes` (after the current line 71 `jDes.add('br hfill', jParams);`), insert:

```matlab
        % Seed-direction control (quaternion/Dirac fibers): preset combobox, hidden for scalar operators.
        jDirRow = gui_river([0 0], [0 2 0 2]);
        gui_component('label', jDirRow, [], 'Direction: ');
        jDirCombo = gui_component('combobox', jDirRow, 'hfill', [], {{'Normal','+X','+Y','+Z','Pick-on-surface'}}, [], [], []);
        java_setcb(jDirCombo, 'ActionPerformedCallback', @(h,e) bst_call(@OnDirCb));
        jDirRow.setVisible(false);                               % no vector operator selected yet -> hidden
        jDes.add('br hfill', jDirRow);
```

- [ ] **Step 5: Update the `ctrl` struct**

In `toolbox/gui/panel_atom_designer.m`, replace the `ctrl = struct(...)` at lines 93–94 with (swap `jConn`/`jGeom` for `jOperator`; add `jDirRow`/`jDirCombo`):

```matlab
    ctrl = struct('jOperator',jOperator, 'jDirRow',jDirRow, 'jDirCombo',jDirCombo, ...
                  'jKernel',jKernel, 'jParams',jParams, ...
                  'jFibers',jFibers, 'jStatus',jStatus, 'atomKeys',{atomKeys});
```

- [ ] **Step 6: Rewrite `CurrentOperator` + add the pure maps, `CurrentDirection`, `ShowDirection`, `OnDirCb`**

In `toolbox/gui/panel_atom_designer.m`, replace `CurrentOperator` (lines 129–133) with:

```matlab
function v = CurrentOperator() %#ok<DEFNU>
    v = 'Laplace-Beltrami';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    v = i_variant_for_item(char(ctrl.jOperator.getSelectedItem()));
end

% Pure display->variant map (also used by Configure and the headless test).
function v = i_variant_for_item(name) %#ok<DEFNU>
    switch name
        case 'Connectomic',        v = 'LB-Connectome';
        case 'Dirac',              v = 'Dirac';
        case 'Dirac (connectome)', v = 'Dirac-Connectome';
        otherwise,                 v = 'Laplace-Beltrami';   % 'Geometric' / unknown
    end
end

% Pure variant->display map.
function nm = i_item_for_variant(variant) %#ok<DEFNU>
    switch variant
        case 'LB-Connectome',    nm = 'Connectomic';
        case 'Dirac',            nm = 'Dirac';
        case 'Dirac-Connectome', nm = 'Dirac (connectome)';
        otherwise,               nm = 'Geometric';
    end
end

% Selected seed-direction preset.
function v = CurrentDirection() %#ok<DEFNU>
    v = 'Normal';
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    v = char(ctrl.jDirCombo.getSelectedItem());
end

% Show/hide the seed-direction control (vector fiber -> show; scalar -> hide).
function ShowDirection(tf) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'AtomDesigner');  if isempty(ctrl), return; end
    ctrl.jDirRow.setVisible(logical(tf));
end

% Relay the direction-combo change to the designer.
function OnDirCb()
    cb = getappdata(0, 'AtomDesignerCB');
    if ~isempty(cb) && isfield(cb,'Direction'), cb.Direction(); end
end
```

- [ ] **Step 7: Preselect the operator in `Configure`**

In `toolbox/gui/panel_atom_designer.m`, replace the selection-init line 103 with:

```matlab
    ctrl.jOperator.setSelectedItem(i_item_for_variant(variant0));
```

- [ ] **Step 8: Update the API header comment**

In `toolbox/gui/panel_atom_designer.m`, update the documented API block (lines ~10–18) so `CurrentOperator` lists all four operators and add the new verbs:

```matlab
%   v    = panel_atom_designer('CurrentOperator')    % 'Laplace-Beltrami'|'LB-Connectome'|'Dirac'|'Dirac-Connectome'
%   d    = panel_atom_designer('CurrentDirection')   % 'Normal'|'+X'|'+Y'|'+Z'|'Pick-on-surface'
%   panel_atom_designer('ShowDirection', tf)         % show/hide the seed-direction control
```

Also extend the `Configure` `cb` doc to note the new `Direction` handle: `% cb = struct(Kernel/Operator/Param/Direction/Fibers/Save handles)`.

- [ ] **Step 9: (implementer) static-check + grep for stale toggle references**

`check_matlab_code` on `panel_atom_designer.m`. Then:
```bash
grep -n "jConn\|jGeom" toolbox/gui/panel_atom_designer.m
```
Expected: **no matches** (both toggles fully replaced).

- [ ] **Step 10: Commit**

```bash
git add toolbox/gui/panel_atom_designer.m dev/tests/test_designer_operator_map.m
git commit -m "feat(gui): atom designer panel - 4-operator selector + seed-direction control

Replace the 2 operator toggles with a 4-item combobox (Geometric/Connectomic/
Dirac/Dirac-connectome) and add a preset seed-direction control (hidden for
scalar fibers). New accessors CurrentDirection/ShowDirection + pure item<->variant maps.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate:** run `matlab -batch "... test_designer_operator_map"` → `PASSED`. Live: open `view_atom_designer(surf)`; the operator combobox lists all 4; the Direction row exists (initially hidden).

---

### Task 5: `view_atom_designer` — fiber-driven render (scalar paint vs magnitude + quivers)

**Files:**
- Modify: `toolbox/gui/view_atom_designer.m` (`i_eval_atom` `:302–305`; `Generate` `:167–180`; add `seedDir` state + helpers; `i_build_basis`/`OperatorChanged`/`i_seed` set the default direction)

**Interfaces:**
- Consumes: `bst_eigenfilter('Atom', …)` → `[Wloc, gv, V3, …]` (Task 1); `bst_atom_default_dir(ax, seedVtx)` (Task 2); `panel_atom_designer('CurrentOperator')` (Task 4); `figure_3d('SetShowSourceVectors', hFig, iTess, on)`; appdata `QuiverVectorOverride`.
- Produces: the designer renders any of the 4 operators — scalar → density/peak paint (unchanged); vector → `|V3|` magnitude (`source` colormap) + a quiver overlay reoriented by the default seed direction.

- [ ] **Step 1: Add `seedDir` closure state + a default-direction reset**

In `toolbox/gui/view_atom_designer.m`, add a closure variable alongside `seedVtx` (near lines 51–55) after `seedVtx` is established:

```matlab
    seedDir = bst_atom_default_dir(ax, seedVtx);     % app-side default (1 for scalar; seed normal for Dirac)
```

Add a nested helper (near the other nested helpers, e.g. after `i_seed`):

```matlab
    function i_reset_dir()
        seedDir = bst_atom_default_dir(ax, seedVtx);  % recompute the default whenever ax or the seed changes
    end
```

- [ ] **Step 2: Make `i_eval_atom` fiber-aware (return V3, reduce vector to magnitude)**

In `toolbox/gui/view_atom_designer.m`, replace `i_eval_atom` (lines 302–305) with:

```matlab
function [W, V3] = i_eval_atom(s, ax, kernel, kp, V, nV, seedDir) %#ok<INUSL>
    if nargin < 7, seedDir = []; end
    [Wloc, gv, V3] = bst_eigenfilter('Atom', ax, kernel, kp, s, seedDir);
    if isempty(V3)                                   % scalar fiber: paint the modal field directly
        W = zeros(nV, ax.nT);  W(gv,:) = Wloc;
    else                                             % vector fiber: paint per-frame |field| magnitude
        Wmag = i_field_mag(Wloc, numel(gv));         % [nGv x nT] RMS over components
        W = zeros(nV, ax.nT);  W(gv,:) = Wmag;
        if size(V3,1) < nV, V3(end+1:nV, :) = 0; end % pad to the displayed patch vertex count (quiver guard)
    end
end

% Reduce a real/complex/vector field [nc*nRows x nT] to per-row magnitude [nRows x nT] (mirrors the
% Dynamics panel's i_paintable_scalar; scalar passes through).
function s = i_field_mag(F, nRows)
    if ~isreal(F), F = abs(F); end
    if size(F,1) == nRows, s = F; return; end
    if mod(size(F,1), nRows) == 0
        nc = size(F,1) / nRows;
        s  = reshape(sqrt(sum(reshape(F, nc, nRows, []).^2, 1)), nRows, []);
    else
        s = F;
    end
end
```

- [ ] **Step 3: Add quiver set/clear helpers**

In `toolbox/gui/view_atom_designer.m`, add two nested helpers (near `i_set_cmap`):

```matlab
    function i_set_quiver(V3)
        setappdata(hFig, 'QuiverVectorOverride', V3);
        try, figure_3d('SetShowSourceVectors', hFig, iTess, 1); catch, end %#ok<CTCH>
    end
    function i_clear_quiver()
        setappdata(hFig, 'QuiverVectorOverride', []);
        try, figure_3d('SetShowSourceVectors', hFig, iTess, 0); catch, end %#ok<CTCH>
    end
```

- [ ] **Step 4: Branch `Generate` on the fiber**

In `toolbox/gui/view_atom_designer.m`, replace `Generate` (lines 167–180) with:

```matlab
    function Generate()
        try
            [W, V3] = i_eval_atom(seedVtx, ax, kernel, i_phys2kernel(), V, nV, seedDir);
            if isempty(V3)                                          % scalar fiber: density/peak by kernel class
                [W, isSigned] = i_normalize(W);
                if isSigned, i_set_cmap('stat2'); else, i_set_cmap('source'); end
                i_clear_quiver();
            else                                                    % vector fiber: magnitude + quivers
                pk = max(abs(W(:)));  if pk > 0, W = W / pk; end    % peak-normalize the magnitude
                i_set_cmap('source');
                i_set_quiver(V3);
            end
            GlobalData.DataSet(iDS).Results(iRes).ImageGridAmp = W;  % in-place update, no file I/O
            T2 = getappdata(hFig,'Surface');  T2(iTess).DataMinMax = [min(W(:)) max(W(:))];  setappdata(hFig,'Surface',T2);
            panel_surface('UpdateSurfaceData', hFig, iTess);  panel_surface('UpdateSurfaceColormap', hFig);
            if showFib, i_recolor_fibers(); end                     % keep the fiber overlay in sync with the atom
            Status();
        catch ME
            panel_atom_designer('SetStatus',['Could not propagate: ' regexprep(ME.message,'\s+',' ')]);
        end
    end
```

- [ ] **Step 5: Update the two `i_eval_atom` call sites + reset the default direction on build/seed change**

In `toolbox/gui/view_atom_designer.m`:

- The init-frame call at line 56: `W = i_eval_atom(seedVtx, ax, kernel, struct('lmax',lmax), V, nV, seedDir);`
  (If line 56 is captured only as a scalar `W` for the initial results file, keep the single-output form `[W,~] = i_eval_atom(...)`; verify at edit time and match the surrounding assignment.)
- In `OperatorChanged` (after `i_build_basis(variant)` and the seed-recenter, around line 115): add `i_reset_dir();` before `Regen();`.
- In `i_seed` (lines 271–274), after `seedVtx = vi;` add `i_reset_dir();` so the default direction follows the new seed:

```matlab
    function i_seed(vi)                          % figure_3d native pick callback: clicked cortical vertex -> seed
        if ~strcmpi(state,'design') || isempty(vi), return; end
        seedVtx = vi;  i_reset_dir();  Generate();  i_reset_time();   % new seed -> default dir + diffuse anew
    end
```

- [ ] **Step 6: (implementer) static-check + grep**

`check_matlab_code` on `view_atom_designer.m`. Then:
```bash
grep -n "i_eval_atom\|QuiverVectorOverride\|SetShowSourceVectors\|i_reset_dir" toolbox/gui/view_atom_designer.m
```
Expected: `i_eval_atom` defined with 7 args + both call sites updated; quiver helpers present; `i_reset_dir` called from init, `OperatorChanged`, and `i_seed`.

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/view_atom_designer.m
git commit -m "feat(gui): atom designer renders vector (Dirac) atoms as magnitude + quivers

i_eval_atom is now fiber-aware: scalar keeps the density/peak paint; the vector
fibers paint |V3| magnitude (source colormap) and draw a QuiverVectorOverride
quiver overlay. Default seed direction via bst_atom_default_dir, reset on
operator/seed change.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate (live — this is the deliverable's gate; no headless unit, the render is GUI):** in the live session,
`view_atom_designer(<surf>)`; select each operator in turn:
- Geometric / Connectomic → drop a seed → scalar density renders (no quivers).
- Dirac / Dirac (connectome) → drop a seed → cortex colored by `|V3|` magnitude (`source` colormap) AND black quivers appear, oriented along the seed normal (default).
Confirm switching Dirac → Geometric clears the quivers.

---

### Task 6: `view_atom_designer` — interactive seed-direction control

**Files:**
- Modify: `toolbox/gui/view_atom_designer.m` (`cb` struct `:72–76`; add `OnDirectionChanged`; direction resolver; `pickMode` state + `i_seed` dir branch; `i_sync_dir_control` on build)

**Interfaces:**
- Consumes: `panel_atom_designer('CurrentDirection')`, `panel_atom_designer('ShowDirection', tf)` (Task 4); `bst_atom_default_dir` (Task 2); `bst_eigenfilter('Fiber', ax)`.
- Produces: the Direction combobox reorients the Dirac atom — `Normal` → seed normal; `+X/+Y/+Z` → unit axes; `Pick-on-surface` → `unit(V(target) − V(seed))`; the control is shown only for the quaternion fiber.

- [ ] **Step 1: Add `pickMode` state**

In `toolbox/gui/view_atom_designer.m`, add a closure variable near `seedDir`:

```matlab
    pickMode = 'seed';       % 'seed' = clicks move the seed; 'dir' = next click sets the impulse direction
```

- [ ] **Step 2: Wire the `Direction` callback into the `cb` struct**

In `toolbox/gui/view_atom_designer.m`, extend the `cb = struct(...)` (lines 72–76) with a `Direction` handle:

```matlab
    cb = struct('Kernel',    @()bst_call(@KernelChanged), ...
                'Operator',  @()bst_call(@OperatorChanged), ...
                'Param',     @()bst_call(@ParamChanged), ...
                'Direction', @()bst_call(@OnDirectionChanged), ...
                'Fibers',    @(s)bst_call(@()OnToggleConnectome(s)), ...
                'Save',      @()bst_call(@SaveAtom));
```

- [ ] **Step 3: Add the direction handler + resolver + fiber-sync helper**

In `toolbox/gui/view_atom_designer.m`, add nested functions (near `i_reset_dir`):

```matlab
    function OnDirectionChanged()
        name = panel_atom_designer('CurrentDirection');
        if strcmpi(name, 'Pick-on-surface')
            pickMode = 'dir';                          % arm a one-shot direction pick
            panel_atom_designer('SetStatus', 'Click the cortex to point the impulse from the seed.');
            return;
        end
        seedDir = i_resolve_dir(name);
        Generate();
    end

    function d = i_resolve_dir(name)
        switch name
            case '+X', d = [1 0 0];
            case '+Y', d = [0 1 0];
            case '+Z', d = [0 0 1];
            otherwise, d = bst_atom_default_dir(ax, seedVtx);   % 'Normal'
        end
    end

    function i_sync_dir_control()                       % show the control only for the vector (quaternion) fiber
        [~, kind] = bst_eigenfilter('Fiber', ax);
        panel_atom_designer('ShowDirection', strcmp(kind, 'quaternion'));
    end
```

- [ ] **Step 4: Handle the direction-pick in `i_seed`**

In `toolbox/gui/view_atom_designer.m`, extend `i_seed` (from Task 5 Step 5) so an armed direction-pick sets `seedDir` instead of moving the seed:

```matlab
    function i_seed(vi)                          % figure_3d native pick callback
        if ~strcmpi(state,'design') || isempty(vi), return; end
        if strcmp(pickMode, 'dir')               % armed by 'Pick-on-surface': point from seed to the click
            v = V(vi,:) - V(seedVtx,:);  nv = norm(v);
            if nv > 0, seedDir = v / nv; end
            pickMode = 'seed';  Generate();  return;
        end
        seedVtx = vi;  i_reset_dir();  Generate();  i_reset_time();
    end
```

- [ ] **Step 5: Sync the control visibility on init and operator change**

In `toolbox/gui/view_atom_designer.m`:
- After the panel is configured at init (after line 80 `panel_atom_designer('Configure', …)`), add `i_sync_dir_control();`.
- In `OperatorChanged` (after `i_build_basis(variant)`), add `i_sync_dir_control();` (alongside the `i_reset_dir();` from Task 5).

- [ ] **Step 6: (implementer) static-check + grep**

`check_matlab_code` on `view_atom_designer.m`. Then:
```bash
grep -n "OnDirectionChanged\|i_resolve_dir\|i_sync_dir_control\|pickMode\|'Direction'" toolbox/gui/view_atom_designer.m
```
Expected: the `Direction` cb wired; handler/resolver/sync present; `pickMode` toggled in `i_seed`; `i_sync_dir_control` called at init and on operator change.

- [ ] **Step 7: Commit**

```bash
git add toolbox/gui/view_atom_designer.m
git commit -m "feat(gui): atom designer seed-direction control (presets + pick-on-surface)

Wire the panel Direction combobox: Normal/+X/+Y/+Z reorient the Dirac atom;
Pick-on-surface arms a one-shot cortex pick -> unit(target-seed). The control
shows only for the quaternion fiber.

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

**Controller gate (live — the deliverable is GUI):** in the live session, open `view_atom_designer(<surf>)`, select **Dirac**, drop a seed. Then:
- The Direction combobox is visible; `+X/+Y/+Z` reorient the quivers to the chosen axis (`|dot|≈1`).
- `Normal` points the quivers along the seed normal.
- `Pick-on-surface` → the status prompts; the next cortex click points the impulse from the seed toward the click (`unit(target−seed)`).
- Switch to **Geometric** → the Direction control hides and the atom renders scalar.

---

## Integration verification (controller, after Task 6)

Run the whole headless suite in one process, then the live walkthrough:

```bash
matlab -batch "addpath(genpath('toolbox')); addpath('dev/tests'); setenv('BST_TEST_SURF', <surf>); \
  test_atom_fiber_decode; test_atom_default_dir; test_realise_core_delegates; test_designer_operator_map"
```
Expected: four `... PASSED` lines.

Live (MATLAB-MCP `evaluate_matlab_code`, protocol `preventad`, do NOT `clear`):
1. `view_atom_designer(<surf>)` — all 4 operators render (scalar density / Dirac magnitude + quivers).
2. Direction presets reorient the Dirac atom; Pick-on-surface works; control hidden for scalar.
3. `panel_bst_dynamics` Dynamics Design preview of a scalar and a Dirac atom is **unchanged** (byte-equivalent relocation) — the SP2a regression guard.

Then hand off to `superpowers:requesting-code-review` (final whole-branch review = **opus**, per the execution policy — it caught a stale-node regression on SP1).

---

## Self-Review (checked against the spec)

**Spec coverage:**
- §2 Atom/Fiber taxonomy — `Atom` extended to `[W,gv,V3,isSigned]`, `Fiber` unchanged, decode relocated → **Task 1**; `i_atom_realise_core` consumes V3 → **Task 3**. ✓
- §3.1 operator selector → 4 → **Task 4**. ✓
- §3.2 fiber-driven realise + render (scalar paint / magnitude + quivers) → **Task 5**. ✓
- §3.3 seed-direction control + `bst_atom_default_dir` relocation → **Task 2** (helper + panel rewire) + **Task 4** (widget) + **Task 6** (plumbing). ✓
- §3.4 scale/mm calibration unchanged — untouched (no task needed). ✓
- §4 components table — every listed file has a task; new `bst_atom_default_dir.m` = Task 2. ✓
- §6 testing — unit byte-equivalence (Task 1), live render + direction (Tasks 5/6 controller gates). ✓
- §7 out-of-scope (SP2b) — not planned. ✓
- §8 risks — byte-equivalence pinned by Task 1 test; vector render reuses the panel quiver idiom (`QuiverVectorOverride`, verified to bypass `nComponents==3`) → Tasks 1/5. ✓

**Placeholder scan:** no TBD/"add error handling"/"similar to Task N"; every code step shows complete code.

**Type consistency:** `Atom` → `[W, gv, V3, isSigned]` used identically in Tasks 1/3/5; `bst_atom_default_dir(ax, seedVert)` signature identical in Tasks 2/5/6; panel verbs `CurrentOperator`/`CurrentDirection`/`ShowDirection` and pure maps `i_variant_for_item`/`i_item_for_variant` consistent across Tasks 4/5/6; `i_eval_atom(s, ax, kernel, kp, V, nV, seedDir)` 7-arg form consistent across Tasks 5/6.

**Note on GUI tasks (5, 6):** the deliverables are figure rendering, verified by the controller live pass rather than a headless unit (per the execution policy — subagents cannot run MATLAB/Brainstorm). The byte-equivalence that a headless test *can* pin (the decode) is fully covered by Tasks 1–4.
