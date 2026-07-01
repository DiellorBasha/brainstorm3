# Dynamics Analyze / Reconstruct (Sub-project C) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run source maps through the designed frame → a spatial scalogram (per-scale energy over time) + synthesis residual, previewed per window in the panel and (opt-in) batched over the series by a process_, plus JTVAtoms localization into a separate dynamicsmat.

**Architecture:** A new `bst_eigenwavelet('Scalogram', ax, gCell, C)` core reuses B's cached per-hemi projection `C` to produce per-member energy `[3×nT×M]` {Global,LH,RH} + a tightness residual. The panel's Analyze reuses `i_apply_projection` (B) + `i_frame_response` (B), builds a scalogram `TimefreqMat`, and opens `view_timefreq`. An opt-in `process_source_frame` iterates contiguous windows. JTVAtoms localizes bands → a separate `dynamics_*.mat`.

**Tech Stack:** MATLAB (Brainstorm fork). Verification via the brainstorm-dev MATLAB MCP (`check_matlab_code` static; controller-run headless + live) + `git`/`grep`.

## Global Constraints

- Files: NEW `toolbox/eigen/bst_eigenwavelet.m` verb (modify), NEW `toolbox/process/functions/process_source_frame.m`, modify `toolbox/gui/panel_bst_dynamics.m`, and new headless tests under `dev/`. No other files.
- **Static, scalar-only:** the scalogram is the STATIC spatial-frame analysis on scalar operators (LB/LB-Connectome). Dynamic ts/js members are excluded from `gCell` (reuse B's `i_frame_response` gather, which already skips them). Dirac/vector keeps the "scalar-only for now" guard.
- **Window-centric:** the time window is the unit. Panel Analyze = current 4 s window (default). `process_source_frame` runs the WHOLE series only when explicitly invoked, in **contiguous non-overlapping windows**, never loading the whole series at once.
- **Reuse:** B's `i_apply_projection` (cache key `srcResult|iWin|operator`) and `i_frame_response` (returns `.gCell` static member handles + `.A`); `bst_eigenwavelet` `Analysis`/`Synthesis`/`JTVAtoms`; `view_timefreq`; `bst_dynamics`/`view_dynamics`.
- **Scalogram = `TimefreqMat`:** `TF [3 × nT × M]`, `RowNames={'Global','LH','RH'}`, `Freqs` = member scale centers (`√λ` gain-weighted centroid), `Measure='power'`, `Time` = window (or series) time vector, `Method='framescalogram'`.
- **Residual:** `‖F_modal − Frec/A‖/‖F_modal‖` in the K-mode subspace, `A=min_λ Σ_m g_m(λ)²` (frame lower bound). Tight frame ⇒ ≈0.
- **Panel Analyze preview file:** `view_timefreq` needs a saved file, so Analyze does a **find-or-replace** of a single `Frame scalogram (window)` timefreq in the source's study (not a truly transient object; keeps the DB clean). The process saves a distinct `Frame scalogram (series)`.
- **Static-only for implementers:** implementers run `git`, `grep`, `check_matlab_code`, and WRITE headless tests but DO NOT run/evaluate MATLAB or launch Brainstorm. The controller runs the consolidated live pass.
- Commit after each task: `feat(dynamics): <summary>` with the standard session trailers. Branch `development`. Do not push.

---

### Task 1: `bst_eigenwavelet('Scalogram', ax, gCell, C)` core + headless test

**Files:**
- Modify: `toolbox/eigen/bst_eigenwavelet.m` (add the `Scalogram` subfunction; it dispatches via the existing `eval(macro_method)`)
- Create (test): `dev/test_frame_scalogram.m`

**Interfaces:**
- Consumes: `ax` (from `bst_eigen('Axes')`: `ax.Phi{h}`, `ax.Lambda{h}`, `ax.Mass{h}`, `ax.GlobalVertices{h}`); `gCell` = `{@(l) g_m(l)}` static member handles; `C` = per-hemi projection cell `C{h}` `[K_h × nT]`.
- Produces: `scal = bst_eigenwavelet('Scalogram', ax, gCell, C)` → struct with `energy [3×nT×M]` (rows Global/LH/RH), `residual [1×nT]`, `resScalar` (mean), `centers [1×M]`, `A` (frame lower bound), `W [nV×nT×M]`.

- [ ] **Step 1: Write the failing headless test**

Create `dev/test_frame_scalogram.m`:
```matlab
function tests = test_frame_scalogram
tests = functiontests(localfunctions);
end

function ax = i_toy_axes()
% One-hemisphere toy: identity Phi over K modes, unit mass, lambda grid, gv=1..K.
K = 40;  I = eye(K);
ax = struct();
ax.Phi = {I};  ax.Mass = {speye(K)};  ax.Lambda = {linspace(0.1,10,K)'};  ax.GlobalVertices = {(1:K)'};
end

function test_tight_frame_residual_zero(tc)
ax = i_toy_axes();  K = numel(ax.Lambda{1});  nT = 5;
% itersine tight frame of 6 members over [0,lmax]
lmax = max(ax.Lambda{1});  Nf = 6;  gCell = cell(1,Nf);
for ii=1:Nf, gCell{ii} = bst_eigfilter_design_itersine(struct('member',ii,'Nf',Nf,'lmax',lmax)); end
C = {randn(K,nT)};
scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
verifyLessThan(tc, scal.resScalar, 1e-6);                     % tight frame reconstructs exactly
verifyEqual(tc, size(scal.energy), [3 nT Nf]);
verifyEqual(tc, squeeze(scal.energy(1,:,:)), squeeze(scal.energy(2,:,:)) + squeeze(scal.energy(3,:,:)), 'AbsTol',1e-9); % global = LH+RH (RH=0 here)
end

function test_loose_frame_residual_positive(tc)
ax = i_toy_axes();  K = numel(ax.Lambda{1});  nT = 3;
gCell = { @(l) exp(-0.5*double(l(:))) };                      % single low-pass -> not tight
C = {randn(K,nT)};
scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
verifyGreaterThan(tc, scal.resScalar, 0.1);                  % under-covers the spectrum
end
```
(Single-hemi toy: `GlobalVertices={1:K}` so LH carries all, RH empty → global=LH.)

- [ ] **Step 2: Note the RED state (static-only)**

The test references `bst_eigenwavelet('Scalogram', …)`, which does not exist yet. Do NOT run MATLAB; record the RED state as "undefined Scalogram verb / macro_method dispatch error".

- [ ] **Step 3: Add the `Scalogram` subfunction**

In `toolbox/eigen/bst_eigenwavelet.m`, add (it dispatches via the existing `eval(macro_method)` at the top; place near `Analysis`):
```matlab
%% ===== SCALOGRAM: per-member energy + tightness residual from a cached projection =====
% ax     : bst_eigen('Axes') struct (per-hemi Phi/Lambda/Mass/GlobalVertices)
% gCell  : cell of static member handles {g_m(lambda)}
% C      : per-hemi projection cell C{h} [K_h x nT] = manifold_ft(Phi{h}, Mass{h}, F(gv,:))
% Returns energy [3 x nT x M] {Global,LH,RH}, residual [1 x nT], resScalar, centers [1 x M], A, W [nV x nT x M].
function scal = Scalogram(ax, gCell, C) %#ok<DEFNU>
    M  = numel(gCell);
    nT = size(C{find(~cellfun(@isempty,C),1)}, 2);
    nV = 0; for h=1:numel(ax.GlobalVertices), if ~isempty(ax.GlobalVertices{h}), nV = max(nV, max(ax.GlobalVertices{h}(:))); end, end
    W  = zeros(nV, nT, M);
    eHemi = zeros(2, nT, M);                                  % rows: LH(1), RH(2)
    Fmod2 = 0;  Res2 = 0;                                     % modal energies for the residual
    % frame lower bound A = min over a dense grid of sum_m g_m(l)^2
    lamAll = []; for h=1:numel(ax.Lambda), lamAll=[lamAll; ax.Lambda{h}(:)]; end %#ok<AGROW>
    lg = linspace(max(min(lamAll),eps), max(lamAll), 512)';  Sg = zeros(size(lg));
    for m=1:M, v=gCell{m}(lg); Sg = Sg + real(v(:)).^2; end
    A = min(Sg);  if ~(A>0), A = 1; end
    tol = 1e-3 * max(Sg);                                     % coverage floor for the canonical dual
    for h = 1:numel(ax.Phi)
        Phi = ax.Phi{h};  if isempty(Phi) || isempty(C{h}), continue; end
        Lam = ax.Lambda{h}(:);  gv = ax.GlobalVertices{h}(:);
        Sg2 = zeros(numel(Lam),1);
        for m = 1:M
            gm = gCell{m}(Lam);  gm = real(gm(:));
            Wm = manifold_ift(Phi, gm .* C{h});               % [nGv x nT]
            W(gv,:,m) = Wm;
            eHemi(min(h,2),:,m) = sum(Wm.^2, 1);
            Sg2 = Sg2 + gm.^2;
        end
        Fmod = manifold_ift(Phi, C{h});                       % modal-space source (K-mode part)
        dual = Sg2 ./ max(Sg2, tol);                          % canonical dual: ~1 where covered, ~0 in gaps
        Frec = manifold_ift(Phi, dual .* C{h});               % exact on covered modes, 0 in frame gaps
        Fmod2 = Fmod2 + sum(Fmod.^2, 1);
        Res2  = Res2  + sum((Fmod - Frec).^2, 1);
    end
    energy = zeros(3, nT, M);
    energy(2,:,:) = eHemi(1,:,:);  energy(3,:,:) = eHemi(2,:,:);
    energy(1,:,:) = eHemi(1,:,:) + eHemi(2,:,:);
    residual = sqrt(Res2) ./ sqrt(max(Fmod2, eps));
    % per-member characteristic scale center (gain-weighted centroid of sqrt(lambda))
    centers = zeros(1, M);
    for m = 1:M, gm = abs(gCell{m}(lg)); if sum(gm)>0, centers(m) = sqrt(sum(lg.*gm)/sum(gm)); end, end
    scal = struct('energy',energy, 'residual',residual, 'resScalar',mean(residual), ...
                  'centers',centers, 'A',A, 'W',W);
end
```

- [ ] **Step 4: Static check**

`check_matlab_code` on `bst_eigenwavelet.m` and `dev/test_frame_scalogram.m`. Expected: no errors.

- [ ] **Step 5: Grep**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -c "function scal = Scalogram" toolbox/eigen/bst_eigenwavelet.m   # 1
```

- [ ] **Step 6: Commit**

```bash
git add toolbox/eigen/bst_eigenwavelet.m dev/test_frame_scalogram.m
git commit -m "feat(dynamics): bst_eigenwavelet('Scalogram') core (frame energy + residual) + test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```
*(Controller runs `test_frame_scalogram` live.)*

---

### Task 2: Panel "Analyze (window)" — scalogram TimefreqMat + view

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `OnAnalyzeWindow`, `i_scalogram_timefreq`; add an "Analyze" toolbar button)

**Interfaces:**
- Consumes: B's `i_apply_projection(st, ax, D, iWin, nV)` → `[C, ~]`; `i_frame_response(st, ax)` → `.gCell`; `bst_eigenwavelet('Scalogram', …)`; `i_atom_op`, `i_atom_axes`, `i_overlay_nv`, `i_cursor_window`, `i_cs`.
- Produces: `OnAnalyzeWindow` writes/replaces a `Frame scalogram (window)` timefreq in the source study and opens it; sets the residual in `jAtomInfo`.

- [ ] **Step 1: Add the TimefreqMat builder helper**

Add to `panel_bst_dynamics.m`:
```matlab
% Build a scalogram TimefreqMat [3 x nT x M] {Global,LH,RH} from a bst_eigenwavelet('Scalogram') result.
function FileMat = i_scalogram_timefreq(scal, timeVec, surfaceFile, dataFile, comment)
    FileMat = db_template('timefreqmat');
    FileMat.TF        = scal.energy;                          % [3 x nT x M]
    FileMat.Time      = timeVec(:)';
    FileMat.Freqs     = scal.centers(:);                      % sqrt(lambda) scale centers
    FileMat.RowNames  = {'Global','LH','RH'};
    FileMat.Measure   = 'power';
    FileMat.Method    = 'framescalogram';
    FileMat.DataType  = 'matrix';
    FileMat.SurfaceFile = surfaceFile;
    FileMat.nAvg = 1;  FileMat.Leff = 1;
    if ~isempty(dataFile), FileMat.DataFile = file_short(dataFile); end
    FileMat.Comment = comment;
end
```

- [ ] **Step 2: Add `OnAnalyzeWindow`**

```matlab
% Analyze the current 4 s window's source through the frame -> scalogram TimefreqMat + residual.
function OnAnalyzeWindow() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcResult') || isempty(D.srcResult)
        ctrl.jAtomInfo.setText('Analyze: no real source linked');  return;
    end
    variant = i_atom_op(st);
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        ctrl.jAtomInfo.setText(sprintf('%s: Analyze is scalar-only for now (use Geometric/Connectomic)', variant));  return;
    end
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    fr = i_frame_response(st, ax);
    if fr.nMembers < 1, ctrl.jAtomInfo.setText('Analyze: no static frame members'); return; end
    nV = i_overlay_nv(ax);
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);
    if isempty(iWin), ctrl.jAtomInfo.setText('Analyze: no recording window'); return; end
    bst_progress('start', 'Frame', 'Analyzing the source through the frame...');
    [C, ~] = i_apply_projection(st, ax, D, iWin, nV);          % reuse B's cache
    scal = bst_eigenwavelet('Scalogram', ax, fr.gCell, C);
    tv = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  tv = tv(iWin);
    surf = ax.SurfaceFile;  srcFile = D.srcResult;
    FileMat = i_scalogram_timefreq(scal, tv, surf, srcFile, sprintf('Frame scalogram (window) | %d members', scal_nmembers(scal)));
    TfFile = i_save_scalogram(srcFile, FileMat, 'Frame scalogram (window)');   % find-or-replace in the source study
    bst_progress('stop');
    if ~isempty(TfFile), try, view_timefreq(TfFile, 'SingleSensor'); catch, end, end %#ok<CTCH>
    ctrl.jAtomInfo.setText(sprintf('Analyze: residual %.1f%% (frame completeness)', 100*scal.resScalar));
end
function n = scal_nmembers(scal), n = size(scal.energy, 3); end
```

- [ ] **Step 3: Add the find-or-replace timefreq saver**

```matlab
% Save (find-or-replace by Comment) a timefreq FileMat into the source result's study; return its path.
function TfFile = i_save_scalogram(srcFile, FileMat, tag)
    TfFile = '';
    [sStudy, iStudy] = bst_get('AnyFile', srcFile);
    if isempty(sStudy), return; end
    % reuse an existing same-tag file for this source, else make a new path
    old = '';  iOld = [];
    if isfield(sStudy,'Timefreq') && ~isempty(sStudy.Timefreq)
        for i=1:numel(sStudy.Timefreq)
            if strncmp(sStudy.Timefreq(i).Comment, tag, numel(tag)), old = sStudy.Timefreq(i).FileName; iOld = i; break; end
        end
    end
    if ~isempty(old)
        TfFile = file_fullpath(old);
    else
        TfFile = bst_process('GetNewFilename', bst_fileparts(file_fullpath(srcFile)), 'timefreq_framescalo');
    end
    bst_save(TfFile, FileMat, 'v6');
    if ~isempty(iOld), db_add_data(iStudy, file_short(TfFile), FileMat, iOld);   % replace the same tree slot in place
    else,             db_add_data(iStudy, file_short(TfFile), FileMat); end
    panel_protocols('UpdateNode', 'Study', iStudy);
end
```

- [ ] **Step 4: Add the Analyze toolbar button in `CreatePanel`**

In the EAST `jToolbar2`, after the Apply toggle group, add:
```matlab
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_TIMEFREQ, TB_DIM}, 'Analyze: decompose the current window''s source through the frame -> spatial scalogram + residual', @(h,e)bst_call(@OnAnalyzeWindow));
```

- [ ] **Step 5: Static check + grep**

`check_matlab_code` on the file. Expected: no undefined-ref errors.
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -c "OnAnalyzeWindow\|i_scalogram_timefreq\|i_save_scalogram" toolbox/gui/panel_bst_dynamics.m   # >=4 (defs + button + calls)
```

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m
git commit -m "feat(dynamics): panel Analyze(window) -> frame scalogram TimefreqMat + residual

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 3: "Localize bands" — JTVAtoms → separate `dynamicsmat`

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add `OnLocalizeBands`; add an action)
- Create (test): `dev/test_jtvatoms_localize.m`

**Interfaces:**
- Consumes: `bst_eigenwavelet('JTVAtoms', W, ax, thr)` → `dynamicsmat` (one atom per band: `.vertices` seed, `.times` window, `.region` level set); the `W` from a fresh `Scalogram` call; `bst_dynamics('Save')`, `view_dynamics`.
- Produces: `OnLocalizeBands` saves a separate `dynamics_*.mat` and opens it.

- [ ] **Step 1: Write the failing headless test**

Create `dev/test_jtvatoms_localize.m`:
```matlab
function tests = test_jtvatoms_localize
tests = functiontests(localfunctions);
end
function test_one_atom_per_band_at_peak(tc)
% Synthetic W [nV x nT x M]: band m peaks at vertex 10*m, time m.
nV = 60; nT = 4; M = 3;  W = zeros(nV,nT,M);
for m=1:M, W(10*m, m, m) = 1; end
ax = struct('SurfaceFile','', 'TimeFile','', 'Time', (1:nT), 'tlag',(1:nT));
T = bst_eigenwavelet('JTVAtoms', W, ax, 0.5);
verifyEqual(tc, numel(T.Groups), M);
for m=1:M, verifyEqual(tc, T.Groups(m).vertices, 10*m); end
end
```
(If `JTVAtoms` needs `ax.SurfaceFile` for vertices, the test passes `''` and only checks seed indices, which come from `W`, not the surface.)

- [ ] **Step 2: Note RED state (static-only)** — `JTVAtoms` exists but the test pins its localization; do NOT run MATLAB.

- [ ] **Step 3: Add `OnLocalizeBands`**

```matlab
% Localize each frame band into a marker atom (peak vertex / time window / level set) -> separate dynamicsmat.
function OnLocalizeBands() %#ok<DEFNU>
    [ctrl, st] = i_cs();  if isempty(ctrl) || isempty(st), return; end
    D = getappdata(st.hFig, 'DynamicsOverlay');
    if isempty(D) || ~isfield(D,'srcResult') || isempty(D.srcResult)
        ctrl.jAtomInfo.setText('Localize bands: no real source linked');  return;
    end
    variant = i_atom_op(st);
    if ~any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        ctrl.jAtomInfo.setText(sprintf('%s: Localize bands is scalar-only for now', variant));  return;
    end
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    fr = i_frame_response(st, ax);  if fr.nMembers < 1, ctrl.jAtomInfo.setText('Localize: no static frame members'); return; end
    nV = i_overlay_nv(ax);
    iWin = i_cursor_window(D.srcDS, D.srcResult, 4);  if isempty(iWin), return; end
    bst_progress('start', 'Frame', 'Localizing frame bands...');
    [C, ~] = i_apply_projection(st, ax, D, iWin, nV);
    scal = bst_eigenwavelet('Scalogram', ax, fr.gCell, C);
    axL = ax;  axL.Time = bst_memory('GetTimeVector', D.srcDS, D.srcResult);  axL.Time = axL.Time(iWin);  axL.tlag = axL.Time;
    thr = i_field(st, 'atomThreshold', 0.5);
    T = bst_eigenwavelet('JTVAtoms', scal.W, axL, thr);
    T.SurfaceFile = ax.SurfaceFile;  T.DataFile = D.srcResult;
    T.Comment = sprintf('Frame bands (%s, %d members)', variant, fr.nMembers);
    bst_progress('stop');
    out = bst_fullfile(bst_fileparts(file_fullpath(D.srcResult)), sprintf('dynamics_framebands_%s.mat', datestr_safe()));
    bst_dynamics('Save', out, T);
    try, view_dynamics(out); catch, end %#ok<CTCH>
    ctrl.jAtomInfo.setText(sprintf('Localized %d frame-band atoms', numel(T.Groups)));
end
function s = datestr_safe(), s = sprintf('%09d', mod(round(now*1e5), 1e9)); end
```

- [ ] **Step 4: Wire an action** — add a toolbar button or Atoms-menu item in `CreatePanel`:
```matlab
    gui_component('ToolbarButton', jToolbar2, [], '', {IconLoader.ICON_SCOUT_SEL, TB_DIM}, 'Localize bands: localize each frame band into a marker atom -> a separate dynamics table', @(h,e)bst_call(@OnLocalizeBands));
```

- [ ] **Step 5: Static check + grep**

`check_matlab_code`. Grep `OnLocalizeBands` count ≥ 2.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m dev/test_jtvatoms_localize.m
git commit -m "feat(dynamics): Localize bands (JTVAtoms) -> separate dynamicsmat + test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 4: `process_source_frame` — opt-in whole-series scalogram

**Files:**
- Create: `toolbox/process/functions/process_source_frame.m`

**Interfaces:**
- Consumes: a source result file (input); the frame (rebuilt from a bound dynamics table OR a default itersine N); `bst_dynamics('Axes')`, `i_frame_response`-equivalent gather (inline), `bst_eigenwavelet('Scalogram')`.
- Produces: a saved whole-series scalogram `TimefreqMat` (3 rows, full time axis).

- [ ] **Step 1: Scaffold the process from the standard template**

Create `toolbox/process/functions/process_source_frame.m` following the process template (`GetDescription`/`Run`). `GetDescription`: Comment `Frame scalogram (source)`, Category `Custom`, InputTypes `{'results'}`, OutputTypes `{'timefreq'}`, options: `nframe` (integer, default 6), `variant` (combobox: Geometric/Connectomic), `winsec` (default 4). `Run(sProcess, sInputs)`:
```matlab
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    N   = sProcess.options.nframe.Value{1};
    varOpt = sProcess.options.variant.Value;
    variant = 'Laplace-Beltrami'; if iscell(varOpt) && ~isempty(varOpt) && strcmpi(varOpt{1},'Connectomic'), variant = 'LB-Connectome'; end
    winsec = sProcess.options.winsec.Value{1};
    for iIn = 1:numel(sInputs)
        R = in_bst_results(sInputs(iIn).FileName, 0, 'SurfaceFile','Time','DataFile');
        ax = bst_eigen('Axes', struct('SurfaceFile',R.SurfaceFile, 'Variant',variant, 'nModes',60, ...
                        'TimeWindow',[R.Time(1) R.Time(min(2,end))], 'SampleRate', 1/median(diff(R.Time))));
        lmax = max(ax.Lambda{1}(:));
        gCell = cell(1,N); for ii=1:N, gCell{ii} = bst_eigfilter_design_itersine(struct('member',ii,'Nf',N,'lmax',lmax)); end
        nV = 0; for h=1:numel(ax.GlobalVertices), nV=max(nV,max(ax.GlobalVertices{h}(:))); end
        tv = R.Time;  Fs = 1/median(diff(tv));  nWin = max(1, round(winsec*Fs));
        starts = 1:nWin:numel(tv);
        [iDS, iRes] = bst_memory('LoadResultsFileFull', sInputs(iIn).FileName);   % load once; page windows below
        E = [];  Res = [];  Tall = [];
        for s0 = starts
            iWin = s0:min(s0+nWin-1, numel(tv));
            F = double(bst_memory('GetResultsValues', iDS, iRes, [], iWin, 0));   % reconstruct only this window
            % reduce to scalar magnitude, project per hemi
            F = i_reduce(F, nV);
            C = cell(1,numel(ax.Phi));
            for h=1:numel(ax.Phi), gv=ax.GlobalVertices{h}(:); C{h}=manifold_ft(ax.Phi{h}, ax.Mass{h}, F(gv,:)); end
            scal = bst_eigenwavelet('Scalogram', ax, gCell, C);
            E = cat(2, E, scal.energy);  Res = [Res, scal.residual];  Tall = [Tall, tv(iWin)]; %#ok<AGROW>
        end
        FileMat = db_template('timefreqmat');
        FileMat.TF = E;  FileMat.Time = Tall;  FileMat.Freqs = scal.centers(:);
        FileMat.RowNames = {'Global','LH','RH'};  FileMat.Measure='power';  FileMat.Method='framescalogram';
        FileMat.DataType='matrix';  FileMat.SurfaceFile=R.SurfaceFile;  FileMat.nAvg=1; FileMat.Leff=1;
        FileMat.DataFile = file_short(sInputs(iIn).FileName);
        FileMat.Comment = sprintf('Frame scalogram (series) | itersine x%d, %s', N, variant);
        OutFile = bst_process('GetNewFilename', bst_fileparts(file_fullpath(sInputs(iIn).FileName)), 'timefreq_framescalo_series');
        bst_save(OutFile, FileMat, 'v6');  db_add_data(sInputs(iIn).iStudy, file_short(OutFile), FileMat);
        OutputFiles{end+1} = OutFile; %#ok<AGROW>
    end
end
% scalar magnitude reduction (k*nV -> nV), mirrors the panel's i_paintable_scalar
function s = i_reduce(F, nV)
    if ~isreal(F), F = abs(F); end
    if size(F,1)==nV, s = F; return; end
    if mod(size(F,1),nV)==0, nc=size(F,1)/nV; s = reshape(sqrt(sum(reshape(F,nc,nV,[]).^2,1)),nV,[]); else, s = F; end
end
```

- [ ] **Step 2: Static check** — `check_matlab_code` on the new file. Expected: no undefined-ref errors.

- [ ] **Step 3: Grep** — `grep -c "function OutputFiles = Run" toolbox/process/functions/process_source_frame.m` = 1.

- [ ] **Step 4: Commit**

```bash
git add toolbox/process/functions/process_source_frame.m
git commit -m "feat(dynamics): process_source_frame (opt-in whole-series frame scalogram, window-chunked)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL"
```

---

### Task 5: controller live consolidated pass

**Files:** none (verification; fixes go to the relevant file).

- [ ] **Step 1: Headless tests** (controller): `test_frame_scalogram` (tight→residual≈0, loose→>0, global=LH+RH) and `test_jtvatoms_localize` (one atom/band at peak vertex).
- [ ] **Step 2: GUI smoke** (controller, MCP): open the sub-MTL0002 Dirac session; `OnSetOperator('Laplace-Beltrami')`; `OnDesignFrame` (N=6); `OnAnalyzeWindow` → a `view_timefreq` spectrogram opens (3 rows, scale×time), residual reads ≈0% (tight); screenshot. `OnLocalizeBands` → a `dynamics_*.mat` with 6 atoms opens in `view_dynamics`; screenshot. Run `process_source_frame` on a short window range → a saved series scalogram TimefreqMat.
- [ ] **Step 3: Record** results in the ledger; fix any regression in-place and re-verify; clean up the session (no unintended saves beyond the analysis outputs).

---

## Acceptance criteria (whole plan)
- `bst_eigenwavelet('Scalogram')` returns energy `[3×nT×M]` (global=LH+RH), a residual that is ≈0 for a tight frame and >0 for a loose one, and member scale centers.
- Panel Analyze opens a scale×time spectrogram via `view_timefreq` (find-or-replace preview file) and shows the residual %.
- Localize bands writes a separate `dynamics_*.mat` with one atom per band and opens it.
- `process_source_frame` produces a whole-series scalogram TimefreqMat, computed contiguous-window by window (never the whole series at once).
- `check_matlab_code` clean; controller live pass green.
