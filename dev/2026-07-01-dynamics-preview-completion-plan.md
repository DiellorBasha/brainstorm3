# Dynamics Preview Completion (Sub-project D) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Complete the Apply/Preview loop for the Dirac operator — lift the "scalar-only" guard so Dirac Apply produces a cortex filtered map AND a filtered-sensor overlay on the recording — plus an operator/source compatibility gate.

**Architecture:** A `i_dirac_forward` helper filters the Dirac source in its eigenbasis (`c_filt = g(λ).*project(J)`) and forwards to sensors (`D_filt = L_eig·c_filt`, `L_eig` from `bst_dirac(HeadModel)`, cached). The `i_atom_apply` Dirac branch paints the cortex (`|Phi·c_filt|`) and overlays `D_filt` vs raw on `figure_timeseries`. `i_gate_operators` greys out operators the linked source can't support (`nComponents`).

**Tech Stack:** MATLAB (Brainstorm fork). Verification via the brainstorm-dev MATLAB MCP (`check_matlab_code` static; controller-run headless + live) + `git`/`grep`.

## Global Constraints

- Files: modify `toolbox/gui/panel_bst_dynamics.m`, `toolbox/gui/view_dynamics.m`, `toolbox/gui/figure_timeseries.m`; new headless tests under `dev/`. No other files.
- **Dirac-only sensor forward:** the filtered-sensor view exists ONLY for the Dirac operator with a Dirac-dSPM source (`results_DiracEig_KERNEL_*`). Scalar operators have no sensor view (cortex-only, as shipped in B/C). Keep the scalar Apply path (B/C) unchanged.
- **Resolved math:** `D_filt = L_eig · diag(g(λ)) · K_eig · D`. `L_eig = bst_dirac(HeadModel,'nModes',K,'Tau',tau).Gain` `[nCh × 2K]` (L-then-R eigenmode leadfield). Its `Tau`/`nModes` MUST match the atom's Dirac eigenbasis (`i_atom_axes(st,'Dirac')`) so the `2K` eigenmode ordering aligns with `c_filt = [c_L; c_R]`.
- **Dirac quaternion embedding** (per hemi h): source 3-vector rows of `J[3nV]` embed into the quaternion imag slots; `C{h} = manifold_ft(Phi{h}, Mass{h}, U{h})`; reconstruct via `manifold_ift` then extract the imag 3-vector. Mirror `bst_eigenwavelet`'s `i_hemimap` Dirac layout exactly (`lIn = (0:n-1)*4 + (2:4)`, quaternion `w=0`).
- **Compat gate:** read `R.nComponents` (1 = scalar, 3 = unconstrained vector); scalar → disable Dirac + Tangent (Connection Laplacian) radio items; vector → all enabled; never auto-select a disabled operator; keep the Apply "scalar-only" message as a backstop.
- **Overlay:** filtered-sensor traces are a SECOND, tagged (`'FilteredSensorOverlay'`) line set on the recording figure's `'AxesGraph'` axes, distinct color; behind an appdata guard so ordinary recording figures are unaffected; cleared on Apply-off / operator-leaves-Dirac / session close.
- **Static-only for implementers:** implementers run `git`, `grep`, `check_matlab_code`, and WRITE headless tests but DO NOT run/evaluate MATLAB or launch Brainstorm. The controller runs the consolidated live pass.
- Commit after each task: `feat(dynamics): <summary>` with the standard session trailers. Branch `development`. Do not push.

---

### Task 1: Operator/source compatibility gate

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `i_gate_operators`; call from `SetTarget`)
- Create (test): `dev/test_operator_gate.m`

**Interfaces:**
- Consumes: `ctrl.jOpItems` (the 4 Set-operator radio items), `ctrl.opVariants` (`{'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}`), the linked source result.
- Produces: `i_gate_operators(st)` (sets `jOpItems(k).setEnabled(...)` by `nComponents`); a pure `i_gate_mask(nComponents)` → logical `[1×4]` (enabled per opVariant) that the test pins.

- [ ] **Step 1: Write the failing headless test**

Create `dev/test_operator_gate.m`:
```matlab
function tests = test_operator_gate
tests = functiontests(localfunctions);
end
function test_scalar_source_disables_vector_ops(tc)
% opVariants order: {'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}
m = panel_bst_dynamics('i_gate_mask', 1);          % scalar source
verifyEqual(tc, m, logical([1 1 0 0]));            % LB/LB-Conn enabled; Connection/Dirac disabled
end
function test_vector_source_enables_all(tc)
m = panel_bst_dynamics('i_gate_mask', 3);          % unconstrained vector source
verifyEqual(tc, m, logical([1 1 1 1]));
end
function test_unknown_defaults_permissive(tc)
m = panel_bst_dynamics('i_gate_mask', []);         % unknown -> permissive (all enabled)
verifyEqual(tc, m, logical([1 1 1 1]));
end
```

- [ ] **Step 2: Note RED state (static-only)** — `i_gate_mask` doesn't exist; do NOT run MATLAB.

- [ ] **Step 3: Add `i_gate_mask` + `i_gate_operators`**

Add to `panel_bst_dynamics.m`:
```matlab
% Which operators a source with nComponents supports (order = ctrl.opVariants):
%   {'Laplace-Beltrami','LB-Connectome','Connection Laplacian','Dirac'}.
% Scalar source (1) -> only the two scalar operators; vector (3) -> all; unknown -> permissive.
function m = i_gate_mask(nComponents) %#ok<DEFNU>
    m = true(1,4);
    if isequal(nComponents, 1), m = logical([1 1 0 0]); end
end

% Grey out the Set-operator radio items the linked source can't support (reads nComponents).
function i_gate_operators(st)
    ctrl = bst_get('PanelControls', 'Dynamics');
    if isempty(ctrl) || ~isfield(ctrl,'jOpItems') || isempty(ctrl.jOpItems), return; end
    nComp = [];
    src = i_src_resultfile_from_target(st);
    if ~isempty(src)
        try, R = in_bst_results(src, 0, 'nComponents'); nComp = R.nComponents; catch, end %#ok<CTCH>
    end
    m = i_gate_mask(nComp);
    for k = 1:min(numel(ctrl.jOpItems), numel(m))
        try, ctrl.jOpItems(k).setEnabled(logical(m(k))); catch, end %#ok<CTCH>
    end
    % if the currently-selected op is now disabled, fall back to the first enabled one
    for k = 1:numel(ctrl.jOpItems)
        if ctrl.jOpItems(k).isSelected() && ~m(k)
            f = find(m, 1);  if ~isempty(f), ctrl.jOpItems(f).setSelected(1); OnSetOperator(ctrl.opVariants{f}); end
            break;
        end
    end
end

% Resolve the launched source's results filename from the target (link|... form ok), or '' if none.
function src = i_src_resultfile_from_target(st)
    src = '';
    if isempty(st) || ~isfield(st,'hFig') || isempty(st.hFig) || ~ishandle(st.hFig), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcDS') || ~isfield(D,'srcResult') || isempty(D.srcResult), return; end
    src = i_src_resultfile(D);      % the C-era index->filename resolver
end
```

- [ ] **Step 4: Call the gate from `SetTarget`**

At the end of `SetTarget` (after `BuildTree()`), add:
```matlab
    i_gate_operators(getappdata(0, 'DynamicsTarget'));
```

- [ ] **Step 5: Static check + grep** — `check_matlab_code`; `grep -c "i_gate_mask\|i_gate_operators" toolbox/gui/panel_bst_dynamics.m` ≥ 4.

- [ ] **Step 6: Commit**
```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_operator_gate.m
git commit -m "feat(dynamics): operator/source compatibility gate (grey out by nComponents)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 2: Dirac forward core (`i_dirac_forward`) + headless forward-consistency test

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `i_dirac_forward`, `i_dirac_leadfield`)
- Create (test): `dev/test_dirac_forward.m`

**Interfaces:**
- Consumes: `ax = i_atom_axes(st,'Dirac')` (per-hemi `Phi`/`Lambda`/`Mass`/`GlobalVertices`); `bst_dirac(HeadModel,'nModes',K,'Tau',tau).Gain`; the atom's kernel gain `g(λ)`.
- Produces: `[Dfilt, Jfilt, cfilt] = i_dirac_forward(ax, Leig, J, g)` — `Dfilt [nCh × nWin]`, `Jfilt [3nV × nWin]` (filtered source 3-vector field), `cfilt [2K × nWin]`. `Leig = i_dirac_leadfield(st, ax)` (`bst_dirac` Gain, cached).

- [ ] **Step 1: Write the failing headless forward-consistency test**

Create `dev/test_dirac_forward.m`:
```matlab
function tests = test_dirac_forward
tests = functiontests(localfunctions);
end
function test_allpass_roundtrip(tc)
% Toy: 1 hemi, a FULL quaternion basis (Phi = eye(4nV), 4nV modes, Mass = I) so
% project+reconstruct is an exact identity. With all-pass g, i_dirac_forward must round-trip
% the embedded source EXACTLY (Jfilt == J), and Dfilt has the right shape. The full-basis
% round-trip is a real check (not the tautological "Dfilt==Leig*cfilt"); physics is live-verified.
nV = 4;  K = 4*nV;  nT = 3;  nCh = 5;
ax = struct('Phi',{{eye(4*nV)}}, 'Mass',{{speye(4*nV)}}, 'Lambda',{{(1:K)'}}, 'GlobalVertices',{{(1:nV)'}});
J = randn(3*nV, nT);                       % source 3-vector field
Leig = randn(nCh, K);                      % eigenbasis leadfield (1 hemi -> K cols)
g = @(l) ones(size(l));                    % all-pass
[Dfilt, Jfilt, cfilt] = panel_bst_dynamics('i_dirac_forward', ax, Leig, J, g);
verifyEqual(tc, size(Dfilt), [nCh nT]);
verifyEqual(tc, size(cfilt), [K nT]);
verifyEqual(tc, size(Jfilt), [3*nV nT]);
verifyLessThan(tc, max(abs(Jfilt(:) - J(:))), 1e-12);      % all-pass full-basis round-trip is exact
end
function test_zero_gain_zeros(tc)
% g==0 -> cfilt=0 -> Jfilt=0, Dfilt=0 (the filter removes everything).
nV = 4;  K = 4*nV;  nT = 2;  nCh = 5;
ax = struct('Phi',{{eye(4*nV)}}, 'Mass',{{speye(4*nV)}}, 'Lambda',{{(1:K)'}}, 'GlobalVertices',{{(1:nV)'}});
[Dfilt, Jfilt] = panel_bst_dynamics('i_dirac_forward', ax, randn(nCh,K), randn(3*nV,nT), @(l) zeros(size(l)));
verifyLessThan(tc, max(abs(Jfilt(:))), 1e-12);
verifyLessThan(tc, max(abs(Dfilt(:))), 1e-12);
end
```
*(Full-basis toy verifies the algebra — embed→project→gain→reconstruct; physical head model is live-verified.)*

- [ ] **Step 2: Note RED state** — `i_dirac_forward` undefined; do NOT run MATLAB.

- [ ] **Step 3: Add `i_dirac_forward` + `i_dirac_leadfield`**

```matlab
% Dirac sensor forward: filter J in the Dirac eigenbasis and forward to sensors.
%   ax   : i_atom_axes(st,'Dirac') (per-hemi Phi/Lambda/Mass/GlobalVertices; quaternion basis)
%   Leig : eigenbasis leadfield [nCh x 2K] (L-then-R), from i_dirac_leadfield
%   J    : source 3-vector field [3nV x nT]
%   g    : gain handle g(lambda)
% Returns Dfilt [nCh x nT], Jfilt [3nV x nT] (filtered 3-vector field), cfilt [2K x nT].
function [Dfilt, Jfilt, cfilt] = i_dirac_forward(ax, Leig, J, g) %#ok<DEFNU>
    nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    nT = size(J, 2);  Jfilt = zeros(3*nV, nT);  cfilt = [];
    for h = 1:numel(ax.Phi)
        Phi = ax.Phi{h};  if isempty(Phi), continue; end
        idx = ax.GlobalVertices{h}(:);  n = numel(idx);  Lam = ax.Lambda{h}(:);
        % embed the 3-vector source into the quaternion imag slots (w=0), per bst_eigenwavelet i_hemimap
        gIn = reshape([(idx-1)*3+1, (idx-1)*3+2, (idx-1)*3+3].', [], 1);            % global 3-vec rows
        lIn = reshape([(0:n-1)*4+2; (0:n-1)*4+3; (0:n-1)*4+4], [], 1);              % local quat imag slots
        U = zeros(4*n, nT);  U(lIn, :) = J(gIn, :);
        C  = manifold_ft(Phi, ax.Mass{h}, U);            % [K x nT] eigenmode coeffs
        Ch = g(Lam) .* C;                                % gain
        cfilt = [cfilt; Ch]; %#ok<AGROW>                 % stack L-then-R -> [2K x nT]
        Uf = manifold_ift(Phi, Ch);                      % [4n x nT] filtered quaternion field
        Jfilt(gIn, :) = Uf(lIn, :);                      % extract imag 3-vector
    end
    Dfilt = Leig * cfilt;                                % [nCh x nT]
end

% Dirac eigenbasis leadfield L_eig [nCh x 2K] via bst_dirac(HeadModel); cached per (headmodel,nModes,tau).
function Leig = i_dirac_leadfield(st, ax) %#ok<DEFNU>
    Leig = [];
    D = getappdata(st.hFig, 'DynamicsOverlay');  if isempty(D), return; end
    src = i_src_resultfile(D);  if isempty(src), return; end
    R = in_bst_results(src, 0, 'HeadModelFile');
    hmFile = R.HeadModelFile;
    if isempty(hmFile)
        sS = bst_get('AnyFile', src);  hm = bst_get('HeadModelForStudy', []); %#ok<NASGU>
        if isfield(sS,'HeadModel') && ~isempty(sS.HeadModel), hmFile = sS.HeadModel(sS.iHeadModel).FileName; end
    end
    if isempty(hmFile), return; end
    K = size(ax.Lambda{1},1);  tau = 0.5;   % nModes/Tau must match ax; ax built via i_atom_axes(st,'Dirac')
    key = sprintf('%s|%d|%g', hmFile, K, tau);
    M = getappdata(0, 'DynamicsDiracFwd');
    if ~isempty(M) && isstruct(M) && strcmp(M.key, key), Leig = M.Leig; return; end
    HeadModel = in_bst_headmodel(hmFile);
    CompHM = bst_dirac(HeadModel, 'nModes', K, 'Tau', tau);
    Leig = CompHM.Gain;                     % [nCh x 2K]
    setappdata(0, 'DynamicsDiracFwd', struct('key',key, 'Leig',Leig));
end
```

- [ ] **Step 4: Static check + grep** — `check_matlab_code` on both files; `grep -c "i_dirac_forward\|i_dirac_leadfield" toolbox/gui/panel_bst_dynamics.m` ≥ 3.

- [ ] **Step 5: Commit**
```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_dirac_forward.m
git commit -m "feat(dynamics): i_dirac_forward core (filter in Dirac eigenbasis -> cortex + sensor) + test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 3: Filtered-sensor overlay (`view_dynamics` + `figure_timeseries`)

**Files:**
- Modify: `toolbox/gui/view_dynamics.m` (`SetFilteredSensors`/`ClearFilteredSensors` verbs)
- Modify: `toolbox/gui/figure_timeseries.m` (a small guarded overlay-draw helper)

**Interfaces:**
- Consumes: a recording `figure_timeseries` (its `'AxesGraph'` axes).
- Produces: `view_dynamics('SetFilteredSensors', hRec, Dfilt, tWin)` (draws), `view_dynamics('ClearFilteredSensors', hRec)` (removes); `figure_timeseries('DrawFilteredOverlay', hFig)` (the draw hook, no-op unless the overlay appdata is set).

- [ ] **Step 1: Add the `figure_timeseries` overlay-draw helper**

In `figure_timeseries.m`, add a dispatchable subfunction (reachable via its `eval(macro_method)`):
```matlab
% Draw (or clear) a filtered-sensor overlay stashed on this figure by view_dynamics. No-op unless set.
function DrawFilteredOverlay(hFig) %#ok<DEFNU>
    ov = getappdata(hFig, 'FilteredSensorsOverlay');
    delete(findobj(hFig, 'Tag', 'FilteredSensorOverlay'));       % clear prior overlay lines
    if isempty(ov) || ~isfield(ov,'Dfilt') || isempty(ov.Dfilt), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesGraph');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    % scale Dfilt to the raw traces' current display scaling (butterfly): match the axes YLim span
    Dfilt = ov.Dfilt;  t = ov.tWin(:)';
    rng = max(abs(Dfilt(:)));  if rng <= 0, return; end
    yl = get(hAxes, 'YLim');  sc = 0.9 * max(abs(yl)) / rng;      % fit into the current amplitude range
    hold(hAxes, 'on');
    for c = 1:size(Dfilt,1)
        line(t, sc*Dfilt(c,:), 'Parent', hAxes, 'Color',[1 0.3 0.1 0.5], 'LineWidth',0.5, ...
             'Tag','FilteredSensorOverlay', 'HitTest','off', 'PickableParts','none');
    end
end
```
*(Butterfly-view scaling; the controller live pass tunes the scale factor / handles column view.)*

- [ ] **Step 2: Add the `view_dynamics` verbs**

In `view_dynamics.m` (dispatch region near the other verb hooks), add:
```matlab
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'SetFilteredSensors')
        SetFilteredSensors(varargin{2:end});  return;
    end
    if (nargin >= 2) && ischar(varargin{1}) && strcmp(varargin{1}, 'ClearFilteredSensors')
        ClearFilteredSensors(varargin{2});    return;
    end
```
and the functions:
```matlab
function SetFilteredSensors(hRec, Dfilt, tWin)
    if isempty(hRec) || ~ishandle(hRec), return; end
    setappdata(hRec, 'FilteredSensorsOverlay', struct('Dfilt',Dfilt, 'tWin',tWin));
    figure_timeseries('DrawFilteredOverlay', hRec);
end
function ClearFilteredSensors(hRec)
    if isempty(hRec) || ~ishandle(hRec), return; end
    setappdata(hRec, 'FilteredSensorsOverlay', []);
    figure_timeseries('DrawFilteredOverlay', hRec);
end
```

- [ ] **Step 3: Static check + grep** — `check_matlab_code` on both files; `grep -c "FilteredSensorOverlay\|DrawFilteredOverlay\|SetFilteredSensors" toolbox/gui/view_dynamics.m toolbox/gui/figure_timeseries.m`.

- [ ] **Step 4: Commit**
```bash
git add toolbox/gui/view_dynamics.m toolbox/gui/figure_timeseries.m
git commit -m "feat(dynamics): filtered-sensor overlay hooks (view_dynamics + figure_timeseries)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 4: Wire the Dirac branch into `i_atom_apply`

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (`i_atom_apply` Dirac branch)

**Interfaces:**
- Consumes: Task 2 (`i_dirac_forward`, `i_dirac_leadfield`), Task 3 (`view_dynamics('SetFilteredSensors')`), existing `i_rec_figure`-style recording lookup.
- Produces: Dirac Apply → cortex `SetFilteredField` (filtered magnitude) + sensor overlay.

- [ ] **Step 1: Replace the Dirac bail with the Dirac branch**

In `i_atom_apply`, the current guard bails for non-scalar operators. Add, BEFORE that scalar-only guard, a Dirac branch:
```matlab
    % --- Dirac operator: filter in the Dirac eigenbasis -> cortex magnitude + filtered-sensor overlay ---
    if strcmp(variant, 'Dirac')
        ax = i_atom_axes(st, 'Dirac');  if isempty(ax), return; end
        Leig = i_dirac_leadfield(st, ax);
        iWin = i_cursor_window(D.srcDS, D.srcResult, 4);  if isempty(iWin), ctrl.jAtomInfo.setText('Apply: no window'); return; end
        J = double(bst_memory('GetResultsValues', D.srcDS, D.srcResult, [], iWin, 0));   % [3nV x nWin] Dirac field
        g = bst_eigfilter_kernel(kernel, kp);
        bst_progress('start','Atom','Dirac filter -> sensors...');
        [Dfilt, Jfilt] = i_dirac_forward(ax, Leig, J, g);
        bst_progress('stop');
        % cortex: filtered magnitude
        nV = i_overlay_nv(ax);  Jmag = i_paintable_scalar(Jfilt, nV);
        pk = max(abs(Jmag(:))); if pk>0, Jmag = Jmag/pk; end
        if ~isempty(st.hFig) && ishandle(st.hFig), view_dynamics('SetFilteredField', st.hFig, Jmag, (1:nV)', iWin, false); end
        % sensor: overlay Dfilt vs raw (only if L_eig available)
        if ~isempty(Leig) && ~isempty(Dfilt)
            hRec = i_rec_figure(st);
            if isempty(hRec) || ~ishandle(hRec), try, hRec = view_timeseries(st.T.DataFile, 'MEG'); catch, hRec=[]; end, end %#ok<CTCH>
            tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  tWin = tv(iWin);
            if ~isempty(hRec) && ishandle(hRec), view_dynamics('SetFilteredSensors', hRec, Dfilt, tWin); end
            ctrl.jAtomInfo.setText(sprintf('Dirac | %s [Preview: cortex + %d-sensor overlay]', kernel, size(Dfilt,1)));
        else
            ctrl.jAtomInfo.setText('Dirac: cortex filtered (no Dirac-dSPM leadfield -> no sensor view)');
        end
        return;
    end
```
NOTE: `i_rec_figure` was deleted in sub-project A (navigator cleanup). Re-add this minimal local (used by the branch above + the Step-2 cleanup):
```matlab
% The recording DataTimeSeries figure for this session (matches st.T.DataFile), or [] if none open.
function hFig = i_rec_figure(st)
    hFig = [];  global GlobalData; %#ok<TLEV>
    if isempty(st) || ~isfield(st,'T') || isempty(st.T.DataFile), return; end
    for h = bst_figures('GetFiguresByType', {'DataTimeSeries'})'
        [~,~,iDS] = bst_figures('GetFigure', h);
        if ~isempty(iDS) && ~isempty(GlobalData.DataSet(iDS).DataFile) && file_compare(GlobalData.DataSet(iDS).DataFile, st.T.DataFile)
            hFig = h;  return;
        end
    end
end
```
Keep the scalar-operator guard + scalar branch below unchanged (LB/LB-Connectome still magnitude-filter as before).

- [ ] **Step 2: Add cleanup — clear the sensor overlay when Apply turns off / operator leaves Dirac**

In `OnApply` (and `OnSetOperator`), after the state changes, clear the overlay if not in Dirac-Apply:
```matlab
    st = getappdata(0,'DynamicsTarget');
    if ~isempty(st)
        hRec = i_rec_figure(st);
        if ~isempty(hRec) && ishandle(hRec), try, view_dynamics('ClearFilteredSensors', hRec); catch, end, end %#ok<CTCH>
    end
```
(Then `i_atom_preview` re-draws it if still in Dirac-Apply.)

- [ ] **Step 3: Static check + grep** — `check_matlab_code`; `grep -c "strcmp(variant, 'Dirac')\|SetFilteredSensors\|i_dirac_forward" toolbox/gui/panel_bst_dynamics.m`.

- [ ] **Step 4: Commit**
```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): Dirac Apply branch -> cortex filtered map + filtered-sensor overlay

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 5: controller live consolidated pass

**Files:** none (verification; fixes go to the relevant file).

- [ ] **Step 1: Headless** (controller): `test_operator_gate` (masks), `test_dirac_forward` (Dfilt=Leig·cfilt, shapes).
- [ ] **Step 2: GUI smoke** (controller, MCP): sub-MTL0002 Dirac session; set operator Dirac; add a diffusion atom; Apply → cortex shows the filtered Dirac magnitude AND the recording time series shows the `FilteredSensorOverlay` traces (raw vs filtered); edit a kernel param → both re-filter; Apply OFF → overlay cleared; load a scalar source → Dirac/Tangent greyed out. Screenshots. Tune the overlay scale/color if needed (in-place fix).
- [ ] **Step 3: Record** results in the ledger; fix regressions in-place and re-verify; clean up the session.

---

## Acceptance criteria (whole plan)
- Operator gate greys out Dirac/Tangent for a scalar source, all four for a vector source; never leaves a disabled op selected.
- `i_dirac_forward` returns `Dfilt = Leig·cfilt` (verified exactly) + a filtered 3-vector cortex field.
- Dirac Apply paints the cortex filtered magnitude AND overlays `Dfilt` vs raw on the recording; Apply-off clears the overlay.
- Scalar Apply (B/C) unchanged; `check_matlab_code` clean; controller live pass green.
