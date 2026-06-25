# Dynamics "Source" differential suite — Implementation Plan (Spec 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the cortical differential maps (Divergence / Curl / Potential / Stream) a dynamics-native, ephemeral per-frame overlay driven by `panel_bst_dynamics`, and delete `view_helmholtz` / `panel_helmholtz`.

**Architecture:** `view_dynamics` opens a standard `figure_3d` source figure on the unconstrained result and installs a `CustomOverlayFcn` that runs `process_helmholtz('Compute')` each frame, painting the selected differential scalar into `TessInfo.Data` (ephemeral, not saved). `panel_bst_dynamics` provides a four-entry operator combobox and records the selected field into atoms. No new figure type; time-nav, colormap, and the raw source-vector quiver stay native to `figure_3d`.

**Tech Stack:** MATLAB, Brainstorm GUI (`view_*`/`panel_*`/`figure_3d`/`panel_surface`), `process_helmholtz('Compute')`, `tess_operators('Covariant')`.

**Spec:** `docs/superpowers/specs/2026-06-25-dynamics-source-differential-suite-design.md`

## Global Constraints

- **Port NOTHING from `view_helmholtz`'s feature model** — no irrotational/solenoidal labels, no smoothing, no `Vtot` quiver override. The figure-open + overlay is written fresh. (The raw source-vector quiver stays native to `figure_3d`; the overlay sets only the SCALAR.)
- **Four operators, one `Compute` call, cached per time index:** `Divergence→Ht.Div`, `Curl→Ht.Curl`, `Potential→Ht.Phi`, `Stream→Ht.Psi`. Switching operator is a free re-select (no recompute).
- **All four are signed scalars → `ColormapType = 'stat2'`** (diverging).
- **Curl = scalar surface vorticity `(∇×J)·n`** (already what `Compute` returns; do not change the math).
- **Ephemeral:** the overlay writes only `TessInfo.Data`; nothing is saved to disk unless a panel atom action commits it.
- **Per-figure state** lives in `getappdata(hFig,'DynamicsOverlay')` = `struct('Cov',Cov,'Op',char,'Cache',containers.Map,'srcDS',int,'srcResult',int,'iTess',int,'nV',int)`. Panel state stays in `getappdata(0,'DynamicsTarget')`.
- **Deferred, do NOT implement:** eigenmode smoothing (the Scale axis becomes inert — strip its `view_helmholtz('SetSmoothing')` calls, keep the widget), the new vortex detector, `process_poisson`, isolines.
- **MATLAB env:** Brainstorm is live (nogui, iProtocol=1 TutorialAuditory). Do NOT restart/clear Brainstorm or close figures. Edited `.m` auto-reload. Run via the MATLAB MCP. GUI 3-D figures may not fully render headless — figure-open validation is interactive (noted per task).
- **License header** on any new `.m`; authored "Diellor Basha, 2026".

## File structure

| File | Action | Responsibility |
|---|---|---|
| `toolbox/gui/view_dynamics.m` | Modify | Open the source figure + install the differential overlay; own `i_open_source_figure`, `i_dynamics_overlay`, `i_pick_scalar`, `i_find_tess`, `i_minmax`, and the `'Overlay'`/`'RefreshOverlay'` verbs |
| `toolbox/gui/panel_bst_dynamics.m` | Modify | Operator combobox + `OnMeasurement` repoint; `OnSaveCursor`/`OnRecord` read `DynamicsOverlay`; `SetTarget` drops Lambda; Scale-axis smoothing stripped; `i_op_color` for the 4 ops |
| `toolbox/gui/view_helmholtz.m` | Delete | Superseded |
| `toolbox/gui/panel_helmholtz.m` | Delete | Superseded |
| `toolbox/tree/tree_callbacks.m` | Modify (~:1885) | Redirect the "Helmholtz / vorticity (Dirac)" menu entry to `view_dynamics('FromResult', …)` |
| `dev/tests/test_dynamics_pick_scalar.m` | Create | Unit test for `i_pick_scalar` |

---

### Task 1: `i_pick_scalar` — the operator→field selector (pure, testable)

The one cleanly unit-testable unit. A pure function mapping an operator name + a `Compute` result to the per-vertex scalar.

**Files:**
- Modify: `toolbox/gui/view_dynamics.m` (add the local function + a public `'PickScalar'` verb for testing).
- Test: `dev/tests/test_dynamics_pick_scalar.m` (create).

**Interfaces:**
- Produces: `scal = view_dynamics('PickScalar', Ht, Op)` and the local `i_pick_scalar(Ht, Op)`; `Op ∈ {'Divergence','Curl','Potential','Stream'}` → `Ht.Div/Curl/Phi/Psi`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_dynamics_pick_scalar
    Ht = struct('Div',[1;2;3], 'Curl',[4;5;6], 'Phi',[7;8;9], 'Psi',[10;11;12]);
    assert(isequal(view_dynamics('PickScalar', Ht, 'Divergence'), [1;2;3]), 'Divergence->Div');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Curl'),       [4;5;6]), 'Curl->Curl');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Potential'),  [7;8;9]), 'Potential->Phi');
    assert(isequal(view_dynamics('PickScalar', Ht, 'Stream'),     [10;11;12]), 'Stream->Psi');
    ok = false; try, view_dynamics('PickScalar', Ht, 'Bogus'); catch, ok = true; end
    assert(ok, 'unknown operator must error');
    fprintf('PASS test_dynamics_pick_scalar\n');
end
```

- [ ] **Step 2: Run to verify it fails**

Run (MCP `evaluate_matlab_code`): `addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests'); rehash; test_dynamics_pick_scalar`
Expected: FAIL — `'PickScalar'` verb not handled.

- [ ] **Step 3: Add the verb + the pure helper**

In `view_dynamics.m`, add a verb branch at the top of the function (after the `'Redraw'` branch, before `'FromResult'`):
```matlab
    % --- PickScalar verb (pure, for tests): view_dynamics('PickScalar', Ht, Op) ---
    if (nargin >= 1) && ischar(varargin{1}) && strcmp(varargin{1}, 'PickScalar')
        varargout{1} = i_pick_scalar(varargin{2}, varargin{3});
        return;
    end
```
Add the local function (near the other helpers at the bottom):
```matlab
%% ===== operator name -> per-vertex scalar field (one of the Compute outputs) =====
function scal = i_pick_scalar(Ht, Op)
    switch Op
        case 'Divergence', scal = Ht.Div;
        case 'Curl',       scal = Ht.Curl;
        case 'Potential',  scal = Ht.Phi;
        case 'Stream',     scal = Ht.Psi;
        otherwise, error('view_dynamics:badOp', 'Unknown differential operator: %s', Op);
    end
end
```

- [ ] **Step 4: Run to verify it passes**

Run: `test_dynamics_pick_scalar` → `PASS test_dynamics_pick_scalar`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_dynamics.m dev/tests/test_dynamics_pick_scalar.m
git commit -m "feat(dynamics): i_pick_scalar operator->field selector (Divergence/Curl/Potential/Stream)"
```

---

### Task 2: figure-open + `CustomOverlayFcn` differential overlay in `view_dynamics`

Replace the `view_helmholtz(SrcResult)` open (view_dynamics.m:78) with a fresh source-figure open that installs the dynamics differential overlay.

**Files:**
- Modify: `toolbox/gui/view_dynamics.m` (the open path + new locals + `'Overlay'`/`'RefreshOverlay'` verbs).

**Interfaces:**
- Consumes: `i_pick_scalar` (Task 1); `process_helmholtz('Compute', J, Cov)`; `tess_operators(Surf,'Covariant')`.
- Produces: `hFig = i_open_source_figure(SrcResult)` (installs `'DynamicsOverlay'` appdata + `'CustomOverlayFcn'`); `view_dynamics('RefreshOverlay', hFig)` and `view_dynamics('Overlay', hFig)` re-run the per-frame paint. The figure's `'DynamicsOverlay'` struct is the contract `OnSaveCursor`/`OnRecord` read in Task 4.

- [ ] **Step 1: Replace the open call**

In `view_dynamics.m`, change the spatial-view open (line ~78) from:
```matlab
        hFig = view_helmholtz(SrcResult);
```
to:
```matlab
        hFig = i_open_source_figure(SrcResult);
        if isempty(hFig), return; end
```

- [ ] **Step 2: Add the `'Overlay'`/`'RefreshOverlay'` verbs**

Near the top of `view_dynamics` (after the `'PickScalar'` branch), add:
```matlab
    % --- per-frame overlay verbs: view_dynamics('Overlay'|'RefreshOverlay', hFig) ---
    if (nargin >= 2) && ischar(varargin{1}) && any(strcmp(varargin{1}, {'Overlay','RefreshOverlay'}))
        i_dynamics_overlay(varargin{2});
        return;
    end
```

- [ ] **Step 3: Add the figure-open + overlay locals**

Append to `view_dynamics.m`:
```matlab
%% ===== open a figure_3d SOURCE figure on an unconstrained result + install the overlay =====
function hFig = i_open_source_figure(SrcResult)
    global GlobalData;
    hFig = [];
    [iDS, iResult] = bst_memory('GetDataSetResult', SrcResult);
    if isempty(iResult)
        try, [iDS, iResult] = bst_memory('LoadResultsFileFull', SrcResult); catch, iResult = []; end %#ok<CTCH>
    end
    if isempty(iResult)
        bst_error('Could not load this source map.', 'Dynamics source view', 0);  return;
    end
    R = GlobalData.DataSet(iDS).Results(iResult);
    if isempty(R.nComponents) || (R.nComponents ~= 3)
        bst_error('The dynamics differential view requires an unconstrained (3-component) source.', 'Dynamics source view', 0);  return;
    end
    SurfaceFile = R.SurfaceFile;
    bst_progress('start', 'Dynamics source view', 'Loading the covariant operator...');
    Cov  = tess_operators(SurfaceFile, 'Covariant');           % find-or-create
    Surf = in_tess_bst(SurfaceFile, 0);  nV = size(Surf.Vertices, 1);
    bst_progress('stop');
    [hFig, iDSf] = view_surface_data(SurfaceFile, SrcResult, [], 'NewFigure');
    if isempty(hFig), return; end
    iTess = i_find_tess(hFig);
    % full field shown by default (the Data threshold slider drives the overlay magnitude)
    TI = getappdata(hFig,'Surface');
    if iTess <= numel(TI), TI(iTess).DataThreshold = 0; setappdata(hFig,'Surface',TI); end
    D = struct('Cov',Cov, 'Op','Divergence', ...
               'Cache',containers.Map('KeyType','double','ValueType','any'), ...
               'srcDS',iDSf, 'srcResult',iResult, 'iTess',iTess, 'nV',nV);
    setappdata(hFig, 'DynamicsOverlay', D);
    setappdata(hFig, 'CustomOverlayFcn', @(h) i_dynamics_overlay(h));   % fires per frame
    i_dynamics_overlay(hFig);                                          % paint the first frame
end

%% ===== the CustomOverlayFcn: ephemeral per-frame differential scalar =====
function i_dynamics_overlay(hFig)
    if isempty(hFig) || ~ishandle(hFig), return; end
    D = getappdata(hFig, 'DynamicsOverlay');  if isempty(D), return; end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (D.iTess > numel(TessInfo)) || ~ishandle(TessInfo(D.iTess).hPatch), return; end
    [~, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
    if isempty(iT) || iT < 1, iT = 1; end
    if ~isKey(D.Cache, iT)
        Jt = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iT, 0));   % raw 3-vector
        if size(Jt,1) ~= 3*D.nV, return; end
        D.Cache(iT) = process_helmholtz('Compute', Jt, D.Cov);        % one call, all fields
        setappdata(hFig, 'DynamicsOverlay', D);
    end
    scal = i_pick_scalar(D.Cache(iT), D.Op);
    TessInfo(D.iTess).Data         = scal;
    TessInfo(D.iTess).DataMinMax   = i_minmax(scal);
    TessInfo(D.iTess).ColormapType = 'stat2';      % all four operators are signed
    setappdata(hFig, 'Surface', TessInfo);
    panel_surface('UpdateSurfaceColormap', hFig);
end

%% ===== find the Source tess slot painted by view_surface_data =====
function iTess = i_find_tess(hFig)
    iTess = 1;  TessInfo = getappdata(hFig, 'Surface');
    for i = 1:numel(TessInfo)
        if ~isempty(TessInfo(i).DataSource) && strcmpi(TessInfo(i).DataSource.Type, 'Source'), iTess = i; return; end
    end
end

%% ===== symmetric color limits for a signed field (stat2) =====
function mm = i_minmax(s)
    m = max(abs(s(:)));  if isempty(m) || m <= 0, m = eps; end;  mm = [-m m];
end
```

- [ ] **Step 4: Static check + interactive smoke**

Static (MCP `check_matlab_code` on `view_dynamics.m`): no parse errors.
Interactive smoke (the GUI may not render headless — if so, the user runs this; document the result):
```
view_dynamics('FromResult', 'link|Subject01/S01_AEF_20131218_01_notch/results_MN_MEG_KERNEL_260605_0111.mat|Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat');
hF = bst_figures('GetCurrentFigure','3D');
D = getappdata(hF,'DynamicsOverlay'); fprintf('Op=%s nV=%d cacheN=%d\n', D.Op, D.nV, D.Cache.Count);
TI = getappdata(hF,'Surface'); assert(~isempty(TI(D.iTess).Data), 'overlay scalar not painted');
```
Expected: a cortex figure whose scalar is the Divergence map (signed/stat2), `DynamicsOverlay` installed, cache populated for the current frame. If headless-blocked, mark this step **interactive-deferred** and rely on Task 1 + Task 4 tests.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_dynamics.m
git commit -m "feat(dynamics): figure_3d-native differential overlay (CustomOverlayFcn + process_helmholtz Compute)"
```

---

### Task 3: operator combobox + `OnMeasurement` repoint in `panel_bst_dynamics`

Replace the two `Φ`/`Ψ` toggles with a four-entry combobox; repoint `OnMeasurement` onto the figure overlay.

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` — Measurement row (lines 131–137), `ctrl` struct (lines 152–156), `OnMeasurement` (lines 269–279).

**Interfaces:**
- Consumes: `view_dynamics('RefreshOverlay', hFig)` (Task 2); the figure's `'DynamicsOverlay'` struct.
- Produces: `ctrl.jMeasOp` (a combobox); `st.curOp ∈ {'Divergence','Curl','Potential','Stream'}`.

- [ ] **Step 1: Replace the Measurement-row widgets (lines 131–137)**

```matlab
    % MEASUREMENT row (differential operator selector; not an axis) + actions
    jMeas = gui_river([2 2], [0 7 2 7], 'Measurement');
    % gui_component sig: (compType, jParent, constrain, compText, compOptions, compTooltip, compCallback, fontSize)
    % combobox items go in compOptions (arg-5) as {{...}}; tooltip is arg-6, callback arg-7, NO extra [] between them.
    jMeasOp = gui_component('combobox', jMeas, '', [], {{'Divergence','Curl','Potential','Stream'}}, ...
        'Differential operator painted on the cortex (ephemeral; div/curl from process_helmholtz, potential/stream = their Poisson potentials)', ...
        @(h,e)bst_call(@OnMeasurement));
    gui_component('label', jMeas, 'tab', '  Peaks:', [], [], [], []);
    jPeaks = gui_component('text', jMeas, '', '3', {Dimension(java_scaled('value',26), BH)}, 'Extrema kept per sign', []);
    jCtrl.add(jMeas);
```

- [ ] **Step 2: Update the `ctrl` struct (line 156)**

Replace `'jMeasPot',jMeasPot, 'jMeasStr',jMeasStr,` with `'jMeasOp',jMeasOp,`:
```matlab
        'jMeasOp',jMeasOp, 'jPeaks',jPeaks, 'jPhaseItems',jPhaseItems));
```

- [ ] **Step 3: Rewrite `OnMeasurement` (lines 269–279)**

```matlab
%% ===== MEASUREMENT (differential operator descriptor; not an axis) =====
function OnMeasurement() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    name = char(ctrl.jMeasOp.getSelectedItem());          % 'Divergence'|'Curl'|'Potential'|'Stream'
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if ~isempty(D), D.Op = name; setappdata(st.hFig, 'DynamicsOverlay', D); end
    view_dynamics('RefreshOverlay', st.hFig);             % free re-select from the per-frame cache
    st.curOp = name;  setappdata(0, 'DynamicsTarget', st);
end
```

- [ ] **Step 4: Static check**

MCP `check_matlab_code` on `panel_bst_dynamics.m`: no parse errors; no remaining `jMeasPot`/`jMeasStr` references (grep).
Run: `grep -n "jMeasPot\|jMeasStr" toolbox/gui/panel_bst_dynamics.m` → no output.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): operator combobox (Divergence/Curl/Potential/Stream) drives the differential overlay"
```

---

### Task 4: repoint atom recording + `SetTarget` + Scale-axis onto the new overlay

`OnSaveCursor`/`OnRecord` read `'DynamicsOverlay'` (not `'HelmholtzState'`) and the four-operator field set; `SetTarget` drops the Lambda bootstrap; the Scale axis stops calling the deleted `view_helmholtz('SetSmoothing')`; `i_op_color` covers the four ops.

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` — `i_drive` Scale case (254–262), `OnSaveCursor` (501–521), `OnRecord` (596–614), `SetTarget` (799–805), `i_op_color` (767–773), and the `curOp` default (`'Total'`→`'Divergence'` at 501 and 609).

**Interfaces:**
- Consumes: `view_dynamics('RefreshOverlay', hFig)`; `getappdata(hFig,'DynamicsOverlay')` = `{Cache, srcDS, srcResult}`; `i_pick_scalar` operator names.
- Produces: atoms whose `Function ∈ {'divergence','curl','potential','stream'}`, `strength/charge` sampled from the matching `Compute` field.

- [ ] **Step 1: Strip the deleted-smoothing call from `i_drive` Scale (254–262)**

Replace the `case 'scale'` body with (keep the widget inert until the eigenvalue-axis work):
```matlab
        case 'scale'
            % Smoothing (eigenmode low-pass) is deferred to the eigenvalue-axis work; the Scale
            % widget is inert for now. Record the request so it round-trips, drive nothing.
            if loc.extent > 0
                st.curScale = struct('on',1,'name','heat','params',struct('t',loc.extent));
            else
                st.curScale = struct('on',0,'name','heat','params',[]);
            end
```

- [ ] **Step 2: Repoint `OnSaveCursor` (501–521)**

Replace the operator-switch + cache read:
```matlab
    op   = i_field(st, 'curOp', 'Divergence');
    switch op
        case 'Divergence', Func = 'divergence';
        case 'Curl',       Func = 'curl';
        case 'Potential',  Func = 'potential';
        case 'Stream',     Func = 'stream';
        otherwise,         Func = 'divergence';
    end
    % measured descriptor at the cursor (operator scalar at the seed, if localized + cached)
    strength = NaN;  charge = NaN;
    if isfinite(ls.center) && ~isempty(st.hFig) && ishandle(st.hFig)
        D = getappdata(st.hFig, 'DynamicsOverlay');
        if ~isempty(D) && isfield(D,'Cache')
            [~, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
            if ~isempty(D.Cache) && isKey(D.Cache, iT)
                sc = view_dynamics('PickScalar', D.Cache(iT), op);
                if ls.center>=1 && ls.center<=numel(sc), strength = sc(ls.center);  charge = sign(strength); end
            end
        end
    end
```

- [ ] **Step 3: Repoint `OnRecord` (596–614)**

Replace the `St = getappdata(...'HelmholtzState')` block + `view_helmholtz('UpdateFrame')` + the op switch:
```matlab
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D)
        java_dialog('warning', 'Record needs the linked dynamics source view (open via a Dirac result).', 'Record atoms');
        return;
    end
    view_dynamics('RefreshOverlay', st.hFig);                     % make sure the cursor frame is computed+cached
    D = getappdata(st.hFig, 'DynamicsOverlay');
    [TimeVec, iT] = bst_memory('GetTimeVector', D.srcDS, D.srcResult, 'CurrentTimeIndex');
    if isempty(D.Cache) || ~isKey(D.Cache, iT)
        java_dialog('warning', 'No decomposition at the current time.', 'Record atoms');  return;
    end
    Ht = D.Cache(iT);  tCur = TimeVec(iT);
    % scalar field + Function from the current operator (all signed)
    op   = i_field(st, 'curOp', 'Divergence');
    Scal = view_dynamics('PickScalar', Ht, op);
    switch op
        case 'Divergence', Func = 'divergence';
        case 'Curl',       Func = 'curl';
        case 'Potential',  Func = 'potential';
        case 'Stream',     Func = 'stream';
        otherwise,         Func = 'divergence';
    end
    signed = true;
```

- [ ] **Step 4: Drop the Lambda bootstrap from `SetTarget` (799–805)**

Replace lines 799–805 (the `St = getappdata(hFig,'HelmholtzState') ... Lambda` block) with nothing (delete it); `SetTarget` ends at `BuildTree();`. (Smoothing/eigenspectrum is deferred; `st.Lambda` stays `[]` from the struct initializer at line 797.)

- [ ] **Step 5: Update `i_op_color` (767–773)**

```matlab
function c = i_op_color(op)
    switch op
        case 'Divergence', c = [0.95 0.55 0.10];   % sources / sinks  -> orange
        case 'Curl',       c = [0.55 0.20 0.85];   % vorticity        -> purple
        case 'Potential',  c = [0.90 0.75 0.10];   % source potential -> amber
        case 'Stream',     c = [0.30 0.45 0.85];   % stream function  -> blue
        otherwise,         c = [0.40 0.40 0.40];   % gray
    end
end
```

- [ ] **Step 6: Validate the recording logic (extend the dynamics atom test)**

Add to `dev/tests/test_dynamics_atoms.m` a case (or a new `dev/tests/test_dynamics_source_record.m`) that primes a hidden figure's `'DynamicsOverlay'` cache with a known `Ht`, sets `st.curOp='Curl'`, drives the cursor seed, calls `OnRecord`, and asserts the newest group's `Function=='curl'` and `strength` equals the seeded `Ht.Curl` extremum. Concretely:
```matlab
function test_dynamics_source_record
    % Prime a hidden figure with a fake DynamicsOverlay cache; OnRecord must read Ht.Curl.
    hF = figure('Visible','off');  c = onCleanup(@() close(hF));
    nV = 4; Ht = struct('Div',zeros(nV,1),'Curl',[0;5;-3;0],'Phi',zeros(nV,1),'Psi',zeros(nV,1));
    D = struct('Cov',[],'Op','Curl','Cache',containers.Map('KeyType','double','ValueType','any'), ...
               'srcDS',1,'srcResult',1,'iTess',1,'nV',nV);
    D.Cache(1) = Ht;  setappdata(hF,'DynamicsOverlay',D);
    s = view_dynamics('PickScalar', D.Cache(1), 'Curl');
    assert(isequal(s, Ht.Curl), 'overlay cache wiring: PickScalar(Curl) must read Ht.Curl');
    fprintf('PASS test_dynamics_source_record (cache+pick wiring)\n');
end
```
> Note: a full `OnRecord` end-to-end needs the panel's `DynamicsTarget`/`st.T`/`Surf` context (the `test_dynamics_atoms` harness builds it). If that harness is available, extend it to drive `OnRecord` and assert `Function=='curl'`; otherwise this wiring test plus the interactive smoke (Step 7) is the gate.

Run: `test_dynamics_source_record` → PASS.

- [ ] **Step 7: Static check + interactive smoke**

`check_matlab_code` on `panel_bst_dynamics.m`: clean. `grep -n "HelmholtzState\|view_helmholtz" toolbox/gui/panel_bst_dynamics.m` → no output.
Interactive (if GUI available): open `view_dynamics('FromResult', alphaLink)`, switch the combobox Divergence→Curl→Potential→Stream (overlay repaints), "Record at cursor" → an atom with the matching `Function`.

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/tests/test_dynamics_source_record.m
git commit -m "feat(dynamics): atom recording + Scale/SetTarget repointed onto the differential overlay"
```

---

### Task 5: delete `view_helmholtz` / `panel_helmholtz`; redirect the tree menu

**Files:**
- Delete: `toolbox/gui/view_helmholtz.m`, `toolbox/gui/panel_helmholtz.m`.
- Modify: `toolbox/tree/tree_callbacks.m` (~:1885).

**Interfaces:**
- Produces: zero `view_helmholtz`/`panel_helmholtz` references in `toolbox/`.

- [ ] **Step 1: Redirect the tree menu (tree_callbacks.m ~:1885)**

Change the source-results context-menu entry from:
```matlab
                    gui_component('MenuItem', jMenuActivations, [], 'Helmholtz / vorticity (Dirac)', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_helmholtz, filenameRelative));
```
to:
```matlab
                    gui_component('MenuItem', jMenuActivations, [], 'Dynamics (differential source maps)', IconLoader.ICON_RESULTS, [], @(h,ev)bst_call(@view_dynamics, 'FromResult', filenameRelative));
```

- [ ] **Step 2: Delete the two files**

```bash
git rm toolbox/gui/view_helmholtz.m toolbox/gui/panel_helmholtz.m
```

- [ ] **Step 3: Verify the tree has no other reference + nothing else calls them**

Run: `grep -rn "view_helmholtz\|panel_helmholtz" toolbox/`
Expected: no matches (the `panel_eigenfilter_design` reference noted in exploration was to `panel_helmholtz`'s slider builder — confirm it is the reverse dependency `panel_helmholtz`→`panel_eigenfilter_design`, not `panel_eigenfilter_design`→`panel_helmholtz`; if a real caller remains, fix it). If any live reference remains, resolve it before committing.

- [ ] **Step 4: Regression — existing dynamics tests + cleanliness**

Run the dynamics test suite that exists (e.g. `test_dynamics_pick_scalar`, `test_dynamics_source_record`, and any `test_dynamics_atoms`/`test_detect_save` if runnable headless): all PASS.
Confirm `grep -rn "view_helmholtz\|panel_helmholtz" toolbox/` is empty.

- [ ] **Step 5: Commit**

```bash
git add -A toolbox/
git commit -m "refactor(dynamics): delete view_helmholtz/panel_helmholtz; tree menu opens view_dynamics"
```

---

## Self-review

**Spec coverage:**
- §5 architecture (figure_3d-native CustomOverlayFcn) → Task 2 ✅ | four-entry combobox → Task 3 ✅ | `Compute` once, cache per time index, free re-select → Task 2 (`i_dynamics_overlay`) + Task 3 (`RefreshOverlay`) ✅ | atom "Source" recording reads the selected field → Task 4 ✅ | delete view_helmholtz/panel_helmholtz + tree redirect → Task 5 ✅.
- §3 non-goals: smoothing stripped not ported (Task 4 Step 1, Task 4 Step 4 drops Lambda) ✅; no process_poisson / detector / isolines (absent) ✅; nothing ported from view_helmholtz's feature model (fresh overlay) ✅.
- §10 decisions: combobox (Task 3) ✅; labels Divergence/Curl/Potential/Stream (Tasks 3–4) ✅; `stat2` colormap (Task 2 `i_dynamics_overlay`) ✅.

**Placeholder scan:** The interactive GUI smokes (Task 2 Step 4, Task 4 Step 7) are explicitly marked interactive-deferred because 3-D figures may not render headless — they are validation *fallbacks*, with the pure `i_pick_scalar` test (Task 1) and the cache-wiring test (Task 4 Step 6) as the headless gates. No coded step is left unspecified.

**Type consistency:** the operator name set `{'Divergence','Curl','Potential','Stream'}` is identical across `i_pick_scalar` (Task 1), the combobox (Task 3), `OnMeasurement` (Task 3), `OnSaveCursor`/`OnRecord` (Task 4), and `i_op_color` (Task 4). The figure appdata key `'DynamicsOverlay'` with fields `{Cov,Op,Cache,srcDS,srcResult,iTess,nV}` is written in Task 2 and read identically in Tasks 3–4. `view_dynamics('PickScalar'|'RefreshOverlay'|'Overlay', …)` verbs are defined in Tasks 1–2 and called in Tasks 3–4.
