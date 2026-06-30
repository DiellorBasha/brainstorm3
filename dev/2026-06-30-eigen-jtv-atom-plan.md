# Eigen / JTV-atom (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Give the atom designer a single, domain-aware (static / ts / js) single-atom realiser by adding `bst_eigen('Axes')` and `bst_eigenfilter('Atom')`, then retiring `dynamics/bst_atom.m`.

**Architecture:** `bst_eigen` becomes the one place axes are assembled (eigenbasis × time × temporal-frequency); `bst_dynamics('Axes')` is refactored to delegate to it; `bst_eigenfilter` gains a domain-dispatching `Atom` verb; `bst_atom`'s two real jobs (Levelset, marker-localisation accessor) move to `bst_dynamics`, its `Evaluate` is superseded by `bst_eigenfilter('Atom')`, and the file is deleted.

**Tech Stack:** MATLAB R2023b, Brainstorm (dev fork). Math primitives `manifold_ft/ift/jft/ijft` (`toolbox/math/`); kernel registry `bst_eigfilter_kernel` + `bst_eigfilter_jtv_evaluate` (`toolbox/eigen/eigfilter/`).

## Global Constraints

- MATLAB R2023b; Brainstorm dev fork. **No new dependencies.**
- Tests are plain MATLAB assertion scripts under `dev/tests/`, run with `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/<file>.m'); disp('OK')"` (or the brainstorm-dev MATLAB MCP `run_matlab_file`). A test **fails** if any `assert` throws; **passes** if the script completes and prints `OK`.
- The canonical axes struct keeps the field names `bst_dynamics('Axes')` already uses: `Phi, Lambda, Mass, GlobalVertices` (1×N cells per hemisphere), `tlag, omega, nT, NFFT, Fs, Time, Variant, SurfaceFile, EigenMat, Operator`.
- Eigenbasis resolution is **find-or-create** via `tess_eigen(SurfaceFile, Variant, 'nModes', nModes)` (NOT `GetEigenBasis`, which is find-only) — matches the designer's compute-on-demand behaviour.
- `lint` every edited `.m` (brainstorm-dev MATLAB MCP `check_matlab_code`); pre-existing globals/`try,catch`-comma/idiom warnings are acceptable, new structural warnings are not.
- Commit after each task. End every commit message with:
  `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- `toolbox/tree/tree_callbacks.m` is CRLF — not touched here, but preserve CRLF if it ever is.

---

### Task 1: `bst_eigen('Axes', OPTIONS)` + delegate `bst_dynamics('Axes')`

**Files:**
- Modify: `toolbox/eigen/bst_eigen.m` (add a leading-string `'Axes'` branch + local `BuildAxes`)
- Modify: `toolbox/dynamics/bst_dynamics.m:194-227` (`Axes` body → delegate; keep `i_load_time`)
- Test: `dev/tests/test_eigen_axes.m`

**Interfaces:**
- Produces: `ax = bst_eigen('Axes', OPTIONS)` where `OPTIONS` is a struct with required `.SurfaceFile`, optional `.Variant` (default `'Laplace-Beltrami'`), `.nModes` (default `200`), and a time axis given either as `.Time` (explicit [1×nT] seconds) or `.TimeWindow` `[t0 t1]` + `.SampleRate` (Hz). Returns the canonical axes struct (fields per Global Constraints).
- Consumes: `tess_eigen`, `in_bst_operator` (existing).

- [ ] **Step 1: Write the failing test (integration — real surface, the Brainstorm idiom)**

Create `dev/tests/test_eigen_axes.m`. Replace `SURF` with a cortex SurfaceFile from the open protocol (a low-res cortex; `tess_eigen` computes the basis on demand if absent):

```matlab
% test_eigen_axes — bst_eigen('Axes') assembles eigenbasis x time x temporal-frequency
SURF = 'Subject01/.../tess_cortex_pial_low.mat';   % <-- a real cortex SurfaceFile
ax = bst_eigen('Axes', struct('SurfaceFile',SURF, 'Variant','Laplace-Beltrami', 'nModes',20, ...
                              'TimeWindow',[0 0.99], 'SampleRate',100));
% --- temporal/spectral axis ---
assert(ax.nT == 100,                       'nT should be 100 (1 s @ 100 Hz)');
assert(abs(ax.Fs - 100) < 1e-9,            'Fs should be 100 Hz');
assert(abs(ax.tlag(2) - 0.01) < 1e-12,     'tlag step = 1/Fs');
assert(ax.NFFT == ax.nT,                   'NFFT == nT');
assert(abs(ax.omega(2) - (ax.Fs/ax.NFFT)) < 1e-9, 'omega step = Fs/NFFT');
assert(numel(ax.omega) == ax.NFFT && numel(ax.tlag) == ax.nT, 'axis lengths');
% --- spatial axis present ---
assert(iscell(ax.Phi) && ~isempty(ax.Phi{1}) && iscell(ax.Mass) && iscell(ax.GlobalVertices), 'eigenbasis fields');
assert(numel(ax.Lambda{1}) == size(ax.Phi{1},2), 'Lambda matches Phi columns');
disp('OK');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_eigen_axes.m'); disp('done')"`
Expected: FAIL — `bst_eigen` errors (the `'Axes'` verb does not exist; it tries to treat the string `'Axes'` as a data file).

- [ ] **Step 3: Add the `'Axes'` verb + local functions to `bst_eigen.m`**

At the very top of the `bst_eigen` body (immediately after the `function [OutputFiles, Messages, isError] = bst_eigen(Data, OPTIONS)` line and any `global` declarations), add a single verb branch:

```matlab
    % ===== VERB DISPATCH (string first arg) =====
    if ischar(Data) && strcmpi(Data, 'Axes')          % ax = bst_eigen('Axes', OPTIONS)
        OutputFiles = BuildAxes(OPTIONS);  Messages = '';  isError = 0;
        return;
    end
```

Add the two local functions at the end of `bst_eigen.m` (the file is script-style with subfunctions — append after the last existing one):

```matlab
function [Time, tlag, omega, nT, NFFT, Fs] = BuildTimeAxis(OPTIONS)
    % Time = actual time vector (s); tlag = lags from 0 (for ts kernels); omega = temporal-freq grid (Hz).
    if isfield(OPTIONS,'Time') && ~isempty(OPTIONS.Time)
        Time = OPTIONS.Time(:)';  Fs = 1 / (Time(2) - Time(1));
    else
        Fs = OPTIONS.SampleRate;  tw = OPTIONS.TimeWindow;
        Time = tw(1) : 1/Fs : tw(2);
    end
    nT = numel(Time);  NFFT = nT;
    omega = (0:NFFT-1) * (Fs / NFFT);            % temporal-frequency grid (Hz)
    tlag  = (0:nT-1) / Fs;                        % time-lag axis (s)
end

function ax = BuildAxes(OPTIONS)
    if ~isfield(OPTIONS,'SurfaceFile') || isempty(OPTIONS.SurfaceFile)
        error('bst_eigen(''Axes''): OPTIONS.SurfaceFile is required.');
    end
    Variant = 'Laplace-Beltrami'; if isfield(OPTIONS,'Variant') && ~isempty(OPTIONS.Variant), Variant = OPTIONS.Variant; end
    nModes  = 200;                if isfield(OPTIONS,'nModes')  && ~isempty(OPTIONS.nModes),  nModes  = OPTIONS.nModes;  end
    % --- spatial axis: eigenbasis (find-or-create) ---
    EigenMat = tess_eigen(OPTIONS.SurfaceFile, Variant, 'nModes', nModes);
    Op       = in_bst_operator(EigenMat.OperatorFile);
    % --- temporal + spectral axis ---
    [Time, tlag, omega, nT, NFFT, Fs] = BuildTimeAxis(OPTIONS);
    % --- assemble (cell fields after struct() to avoid struct-array widening) ---
    ax = struct('Variant',EigenMat.Variant, 'Time',Time, 'Fs',Fs, 'nT',nT, 'NFFT',NFFT, ...
                'omega',omega, 'tlag',tlag, 'SurfaceFile',OPTIONS.SurfaceFile);
    ax.EigenMat = EigenMat;  ax.Operator = Op;
    ax.Phi = EigenMat.Phi;   ax.Lambda = EigenMat.Lambda;  ax.Mass = Op.Mass;
    ax.GlobalVertices = EigenMat.GlobalVertices;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_eigen_axes.m'); disp('done')"`
Expected: prints `OK` then `done`, no error.

- [ ] **Step 5: Refactor `bst_dynamics('Axes')` to delegate**

Replace the body of `Axes` in `toolbox/dynamics/bst_dynamics.m` (lines ~194-227, keep the signature and `i_load_time` helper) with:

```matlab
function ax = Axes(T, variant, nModes, tWin) %#ok<DEFNU>
    % Recording-bound axes: resolve the bound recording's Time/Fs (this module's binding job),
    % then DELEGATE the canonical eigenbasis x time x freq assembly to bst_eigen('Axes').
    if (nargin < 2) || isempty(variant), variant = 'Laplace-Beltrami'; end
    if (nargin < 3) || isempty(nModes),  nModes  = 200; end
    if (nargin < 4), tWin = []; end
    if isempty(T.SurfaceFile), error('bst_dynamics(''Axes''): no SurfaceFile bound.'); end
    timeFile = T.DataFile;
    if isempty(timeFile) && ~isempty(T.Groups)
        if isfield(T.Groups,'ResultsFile') && ~isempty(T.Groups(1).ResultsFile), timeFile = T.Groups(1).ResultsFile;
        elseif isfield(T.Groups,'DataFile') && ~isempty(T.Groups(1).DataFile),   timeFile = T.Groups(1).DataFile; end
    end
    if isempty(timeFile), error('bst_dynamics(''Axes''): no DataFile/ResultsFile bound for the time axis.'); end
    [Time, ~] = i_load_time(timeFile, tWin);
    ax = bst_eigen('Axes', struct('SurfaceFile',T.SurfaceFile, 'Variant',variant, 'nModes',nModes, 'Time',Time));
    ax.TimeFile = timeFile;
end
```

- [ ] **Step 6: Lint both files**

Run brainstorm-dev MCP `check_matlab_code` on `toolbox/eigen/bst_eigen.m` and `toolbox/dynamics/bst_dynamics.m`.
Expected: no new structural warnings (pre-existing globals/idiom warnings OK).

- [ ] **Step 7: Commit**

```bash
git add toolbox/eigen/bst_eigen.m toolbox/dynamics/bst_dynamics.m dev/tests/test_eigen_axes.m
git commit -m "feat(eigen): bst_eigen('Axes') canonical assembler; bst_dynamics('Axes') delegates

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 2: `bst_eigenfilter('Atom', ax, KernelName, KernelParams, seedVert)`

**Files:**
- Modify: `toolbox/eigen/bst_eigenfilter.m` (add `Atom` local function; it is reached via the existing `eval(macro_method)` dispatch)
- Test: `dev/tests/test_eigenfilter_atom.m`

**Interfaces:**
- Produces: `[W, gv] = bst_eigenfilter('Atom', ax, KernelName, KernelParams, seedVert)` → `W` is `[nV_block × nT]` (the realised atom on the seed's eigenbasis block), `gv` the block's global vertex indices.
- Consumes: `ax` from `bst_eigen('Axes')` (Task 1); `manifold_ft/ift/jft/ijft`; `bst_eigfilter_kernel`; `bst_eigfilter_jtv_evaluate`.

- [ ] **Step 1: Write the failing test (synthetic eigenbasis; ts regression + joint-path duality)**

Create `dev/tests/test_eigenfilter_atom.m`:

```matlab
% test_eigenfilter_atom — domain-aware single-atom realiser on a SYNTHETIC eigenbasis
rng_seed = 7;  nV = 60; K = 20; nT = 64; Fs = 100;
% deterministic synthetic orthobasis (no rng dependence on Date/Math.random)
[Q,~] = qr(reshape(cos(1:(nV*K)), nV, K), 0);     % nV x K, columns orthonormal
Phi = Q;  Lam = (linspace(0.0, 5, K)').^2;  M = speye(nV);
ax = struct('nT',nT,'NFFT',nT,'Fs',Fs);
ax.Phi = {Phi}; ax.Lambda = {Lam}; ax.Mass = {M}; ax.GlobalVertices = {(1:nV)'};
ax.tlag = (0:nT-1)/Fs;  ax.omega = (0:nT-1)*(Fs/nT);
seed = 13;  loc = seed;  c0 = manifold_ft(Phi, M, full(sparse(loc,1,1,nV,1)));

% --- (a) ts atom matches the direct formula (== legacy bst_atom Evaluate path) ---
kp = struct('lmax',max(Lam),'tau',0.3);
g  = bst_eigfilter_kernel('diffusion', kp);
W_ref = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);
W_ts  = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed);
assert(isequal(size(W_ts),[nV nT]), 'ts atom shape');
assert(max(abs(W_ts(:) - W_ref(:))) < 1e-12, 'ts atom == direct manifold_ift formula');

% --- (b) static atom is constant in time ---
W_s = bst_eigenfilter('Atom', ax, 'heat', struct('lmax',max(Lam),'t',0.2), seed);
assert(max(abs(W_s(:,1) - W_s(:,end))) < 1e-12, 'static atom constant over time');

% --- (c) joint-path duality: realising a ts kernel through the JTV path matches the direct path ---
U = zeros(nV,nT); U(loc,1) = 1;
Chat = manifold_jft(Phi, M, U, nT);
Gjs  = bst_eigfilter_jtv_evaluate(g, 'ts', Lam, ax.tlag, nT);   % FFT-bridge of the ts kernel
W_joint = real(manifold_ijft(Phi, Chat .* Gjs, nT));
assert(max(abs(W_joint(:) - W_ref(:))) < 1e-9, 'JTV joint path == direct ts (validates js branch machinery)');
disp('OK');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_eigenfilter_atom.m'); disp('done')"`
Expected: FAIL — `Unrecognized ... 'Atom'` (the verb is not defined).

- [ ] **Step 3: Implement the `Atom` verb in `bst_eigenfilter.m`**

Add this local function to `toolbox/eigen/bst_eigenfilter.m` (it is dispatched by the file's existing `eval(macro_method)`):

```matlab
function [W, gv] = Atom(ax, KernelName, KernelParams, seedVert) %#ok<DEFNU>
    % Realise ONE atom: unit delta at seedVert, propagated through the kernel over the eigenbasis.
    % Domain-aware: static g(lambda) (const in time) | ts g(lambda,t) | js g(lambda,omega) (joint inverse FFT).
    if (nargin < 4) || isempty(seedVert), error('bst_eigenfilter(''Atom''): seedVert required.'); end
    blk = 0;
    for h = 1:numel(ax.GlobalVertices)
        if ~isempty(ax.GlobalVertices{h}) && any(ax.GlobalVertices{h} == seedVert), blk = h; break; end
    end
    if blk == 0, error('bst_eigenfilter(''Atom''): seed vertex %d not in the eigenbasis support.', seedVert); end
    Phi = ax.Phi{blk};  Lam = ax.Lambda{blk};  M = ax.Mass{blk};  gv = ax.GlobalVertices{blk};
    loc = find(gv == seedVert, 1);
    c0  = manifold_ft(Phi, M, full(sparse(loc,1,1,size(Phi,1),1)));    % seed in the eigenbasis [K x 1]
    kp  = KernelParams;  if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = max(Lam); end
    g    = bst_eigfilter_kernel(KernelName, kp);
    meta = bst_eigfilter_kernel('info', KernelName);
    dom  = 'static';  if isfield(meta,'domain') && ~isempty(meta.domain), dom = meta.domain; end
    switch lower(dom)
        case 'static'
            W = repmat(manifold_ift(Phi, g(Lam) .* c0), 1, ax.nT);
        case 'ts'
            W = manifold_ift(Phi, g(Lam, ax.tlag) .* c0);
        case 'js'
            U = zeros(numel(gv), ax.nT);  U(loc,1) = 1;                  % spatial delta, temporal impulse
            Chat = manifold_jft(Phi, M, U, ax.NFFT);                     % [K x NFFT]
            Gjs  = bst_eigfilter_jtv_evaluate(g, 'js', Lam, ax.omega, ax.NFFT);
            W    = real(manifold_ijft(Phi, Chat .* Gjs, ax.nT));
        otherwise
            error('bst_eigenfilter(''Atom''): unknown kernel domain ''%s''.', dom);
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_eigenfilter_atom.m'); disp('done')"`
Expected: prints `OK` then `done`.

- [ ] **Step 5: Lint**

`check_matlab_code` on `toolbox/eigen/bst_eigenfilter.m` — no new structural warnings.

- [ ] **Step 6: Commit**

```bash
git add toolbox/eigen/bst_eigenfilter.m dev/tests/test_eigenfilter_atom.m
git commit -m "feat(eigen): bst_eigenfilter('Atom') domain-aware single-atom realiser (static/ts/js)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 3: Move the marker-localisation accessor `bst_atom` → `bst_dynamics`; repoint `panel_bst_dynamics`

**Files:**
- Modify: `toolbox/dynamics/bst_dynamics.m` (add `Axes`-metadata/`NewLoc`/`Get`/`Set` + helpers `i_state`/`i_pad_cols`/`i_pad_vec`/`i_pad_pos`/`i_type`)
- Modify: `toolbox/gui/panel_bst_dynamics.m` (11 call sites: `bst_atom(` → `bst_dynamics(`)
- Modify: `toolbox/dynamics/bst_atom.m` (delete the moved functions; keep `Evaluate`/`Levelset` for now)
- Test: `dev/tests/test_dynamics_localizer.m`

**Interfaces:**
- Produces: `bst_dynamics('NewLoc', axis)`, `bst_dynamics('Get', G, axis, occ)`, `bst_dynamics('Set', G, axis, occ, loc)` — identical signatures/behaviour to the former `bst_atom` verbs. (Note: `bst_dynamics` already has an `Axes(T,...)` verb for the JOINT axes; the localisation **metadata** accessor formerly `bst_atom('Axes')` is renamed to `bst_dynamics('AxisMeta')` to avoid the name clash.)
- Consumes: the `atomgroup` struct `G` (unchanged schema).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_dynamics_localizer.m`:

```matlab
% test_dynamics_localizer — Get/Set round-trip on a group, via bst_dynamics (moved from bst_atom)
G = bst_dynamics('NewGroup','t');
loc = bst_dynamics('NewLoc','time');  loc.center = 0.4; loc.extent = 0.1;
G = bst_dynamics('Set', G, 'time', 1, loc);
out = bst_dynamics('Get', G, 'time', 1);
assert(abs(out.center - 0.4) < 1e-12 && abs(out.extent - 0.1) < 1e-12, 'time loc round-trip');
assert(strcmp(out.state,'window'), 'extent>0 => window');
M = bst_dynamics('AxisMeta');
assert(any(strcmp({M.name},'source')) && any(strcmp({M.name},'scale')), 'axis metadata');
disp('OK');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_dynamics_localizer.m'); disp('done')"`
Expected: FAIL — `Unknown command ... NewLoc` (not yet in `bst_dynamics`).

- [ ] **Step 3: Move the accessor functions into `bst_dynamics.m`**

Cut these functions from `toolbox/dynamics/bst_atom.m` and paste them into `toolbox/dynamics/bst_dynamics.m` (before its final helpers): `NewLoc`, `Get`, `Set`, and the helpers `i_state`, `i_pad_cols`, `i_pad_vec`, `i_pad_pos`, `i_type`, `i_pad_row` (verify which `i_pad_*` each uses; copy exactly, do not rewrite). Rename the former `bst_atom('Axes')` metadata function to `AxisMeta` (rename the `function A = Axes()` from bst_atom to `function A = AxisMeta()` to avoid clashing with `bst_dynamics`'s existing `Axes(T,...)`). If `bst_dynamics` already defines an identically-named helper (e.g. `i_pad_pos` exists at `bst_dynamics.m:182`), keep the existing one and drop the duplicate — verify the bodies are identical first; if they differ, keep `bst_dynamics`'s and adjust the moved code to use it.

- [ ] **Step 4: Repoint `panel_bst_dynamics.m` call sites**

In `toolbox/gui/panel_bst_dynamics.m`, replace each of these `bst_atom(` calls with `bst_dynamics(` (lines ~187, 217, 503, 509, 519, 781, 783, 791, 792, 793, 837, 854, 855, 1005). Use a single replace-all of the exact token `bst_atom('NewLoc'` → `bst_dynamics('NewLoc'`, `bst_atom('Set'` → `bst_dynamics('Set'`, `bst_atom('Get'` → `bst_dynamics('Get'`. There are no `bst_atom('Axes')` calls in this file (verified), so no `AxisMeta` repoint needed here.

- [ ] **Step 5: Run test + lint**

Run the Step-2 command. Expected: `OK`.
`check_matlab_code` on `bst_dynamics.m` and `panel_bst_dynamics.m` — no new structural warnings.

- [ ] **Step 6: Smoke-check no behavioural drift (manual note)**

Grep to confirm `panel_bst_dynamics.m` has zero remaining `bst_atom(` references:
Run: `grep -c "bst_atom(" toolbox/gui/panel_bst_dynamics.m` → expect `0`.

- [ ] **Step 7: Commit**

```bash
git add toolbox/dynamics/bst_dynamics.m toolbox/dynamics/bst_atom.m toolbox/gui/panel_bst_dynamics.m dev/tests/test_dynamics_localizer.m
git commit -m "refactor(dynamics): move atom localisation accessor bst_atom -> bst_dynamics

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

### Task 4: Move `Levelset`, repoint the designer, delete `bst_atom.m`

**Files:**
- Modify: `toolbox/dynamics/bst_dynamics.m` (add `Levelset` — moved from `bst_atom`)
- Modify: `toolbox/gui/view_atom_designer.m` (`i_eval_atom` → `bst_eigenfilter('Atom')`; `ax` build → `bst_eigen('Axes')`; `Levelset` call → `bst_dynamics`)
- Delete: `toolbox/dynamics/bst_atom.m`
- Test: `dev/tests/test_levelset.m`

**Interfaces:**
- Produces: `LS = bst_dynamics('Levelset', W, gv, thr)` → struct with `.scoutVertices`, `.eventSamples`, `.iRef` (identical to former `bst_atom('Levelset')`).
- Consumes: `bst_eigenfilter('Atom')` (Task 2), `bst_eigen('Axes')` (Task 1).

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_levelset.m`:

```matlab
% test_levelset — Levelset on a synthetic wavelet field, via bst_dynamics (moved from bst_atom)
nV = 40; nT = 20; gv = (1:nV)';
W = zeros(nV,nT);  W(10,:) = 1;  W(11,:) = 0.6;  W(:,5) = W(:,5) + 1;   % spatial + temporal peaks
LS = bst_dynamics('Levelset', W, gv, 0.5);
assert(ismember(10, LS.scoutVertices), 'peak vertex in scout');
assert(~ismember(40, LS.scoutVertices), 'far vertex excluded');
assert(~isempty(LS.eventSamples), 'event samples found');
disp('OK');
```

- [ ] **Step 2: Run test to verify it fails**

Run: `matlab -batch "addpath('<repo>/toolbox'); brainstorm setpath; run('<repo>/dev/tests/test_levelset.m'); disp('done')"`
Expected: FAIL — `Unknown command ... Levelset`.

- [ ] **Step 3: Move `Levelset` into `bst_dynamics.m`**

Cut `function LS = Levelset(W, gv, thr, iRef)` (and any unique helper it uses) from `toolbox/dynamics/bst_atom.m` into `toolbox/dynamics/bst_dynamics.m`. Body is unchanged.

- [ ] **Step 4: Run the Levelset test**

Run the Step-2 command. Expected: `OK`.

- [ ] **Step 5: Repoint `view_atom_designer.m`**

In `toolbox/gui/view_atom_designer.m`:

(a) Re-point the helper `i_eval_atom` to the new realiser (one-function change — both its call sites, the init `W = i_eval_atom(...)` and the one in `Generate`, keep working unchanged). Replace its body:

```matlab
    function W = i_eval_atom(s, ax, kernel, kp, V, nV) %#ok<INUSL>
        [Wloc, gv] = bst_eigenfilter('Atom', ax, kernel, kp, s);   % V unused now; signature kept for callers
        W = zeros(nV, ax.nT);  W(gv,:) = Wloc;
    end
```

(b) Replace the `Levelset` call at `view_atom_designer.m:344`:

```matlab
        LS = bst_dynamics('Levelset', W(ax.GlobalVertices{1},:), ax.GlobalVertices{1}, 0.5);
```

(c) Replace the eigenbasis assembly in `i_build_basis` so `ax` comes from `bst_eigen('Axes')` (carrying `omega`). Change the body of `i_build_basis(newVar)` to:

```matlab
    function i_build_basis(newVar)
        variant = newVar;
        bst_progress('start', 'Atom designer', sprintf('Building/loading %s eigenbasis...', variant));
        ax = bst_eigen('Axes', struct('SurfaceFile',SurfaceFile, 'Variant',variant, ...
                       'nModes',nModes, 'TimeWindow',[0 (nFrames-1)/100], 'SampleRate',100));
        bst_progress('stop');
        lamAll = ax.Lambda{1}(:); if numel(ax.Lambda) > 1 && ~isempty(ax.Lambda{2}), lamAll = [lamAll; ax.Lambda{2}(:)]; end
        lmax = max(lamAll);  lminPos = min(lamAll(lamAll > 1e-9));
        scaleMinMM = mm(lmax);  scaleMaxMM = mm(lminPos);
        rateMinMM2 = scaleMinMM^2;  rateMaxMM2 = scaleMaxMM^2;
    end
```

Note: `nFrames` and `mm` are already closure variables in the designer; `ax.nT`/`ax.tlag`/`ax.omega` now come from `bst_eigen('Axes')`. Remove the later standalone `nFrames = 100; ax.nT = nFrames; ax.tlag = (0:nFrames-1)/100;` line (the axes now carry these) but keep `nFrames = 100;` defined *before* the first `i_build_basis` call.

- [ ] **Step 6: Delete `bst_atom.m` and verify no callers**

```bash
git rm toolbox/dynamics/bst_atom.m
grep -rn "bst_atom(" toolbox --include='*.m' | grep -v "function varargout = bst_atom"
```
Expected: the grep returns **nothing** (no remaining callers).

- [ ] **Step 7: Regression smoke — designer realises the same field**

Run a headless realisation check (`dev/tests/test_designer_smoke.m`), exercising a real surface from the open protocol (replace `SURF` with a cortex SurfaceFile that has an eigenbasis):

```matlab
SURF = 'Subject01/.../tess_cortex_pial_low.mat';   % a real cortex with/for an eigenbasis
ax = bst_eigen('Axes', struct('SurfaceFile',SURF,'Variant','Laplace-Beltrami','nModes',50,'TimeWindow',[0 0.99],'SampleRate',100));
seed = ax.GlobalVertices{1}(1);
[W,gv] = bst_eigenfilter('Atom', ax, 'diffusion', struct('lmax',max(ax.Lambda{1}),'tau',0.3), seed);
assert(isequal(size(W),[numel(gv) ax.nT]) && all(isfinite(W(:))), 'designer realiser end-to-end');
disp('OK');
```

Run it; expected `OK`. (If no eigenbasis exists yet, `tess_eigen` computes it on demand.)

- [ ] **Step 8: Lint + Commit**

`check_matlab_code` on `bst_dynamics.m` and `view_atom_designer.m` — no new structural warnings.

```bash
git add toolbox/dynamics/bst_dynamics.m toolbox/gui/view_atom_designer.m dev/tests/test_levelset.m dev/tests/test_designer_smoke.m
git commit -m "refactor(dynamics): retire bst_atom; designer uses bst_eigen('Axes') + bst_eigenfilter('Atom')

- Levelset moved to bst_dynamics; view_atom_designer realises via bst_eigenfilter('Atom')
- bst_atom.m deleted; no remaining callers

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Done criteria (whole plan)

- `bst_atom.m` deleted; `grep bst_atom(` over `toolbox` returns no callers.
- The designer renders the existing kernels (diffusion/wave/heat/mexhat) identically (Task 2 step (a) proves ts == legacy formula; Task 4 step 7 proves end-to-end).
- A `js` kernel realises without error (Task 2 step (c) validates the joint machinery the `js` branch uses).
- `bst_dynamics('Axes')` still works (delegates); marker DB schema unchanged.
- All four `dev/tests/*.m` print `OK`.
