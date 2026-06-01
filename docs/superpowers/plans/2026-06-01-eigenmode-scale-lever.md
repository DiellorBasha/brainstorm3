# Eigenmode Scale Lever Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A stateful, surface-scoped, live-broadcast eigenmode selection lever (the spatial-scale control) that live-coarsens/smooths a displayed cortical source map via a non-destructive display-time band-limited reconstruction.

**Architecture:** A new `panel_eigenmodes` GUI panel owns all lever logic — the selection state (`GlobalData.UserModes`), the weight-shape math, and a single guarded `ApplyToColumn` entry point that filters a source column. `panel_surface('UpdateSurfaceData')` calls `ApplyToColumn` on the displayed column (the one chokepoint, so it auto-composes with time-stepping); `bst_figures('FireModesChanged')` repaints affected 3D figures. No new math — reconstruction reuses `bst_eigenmodes_filter(...,'custom','TransferFn',...)`.

**Tech Stack:** MATLAB, Brainstorm toolbox (`eval(macro_method)` dispatch, Java/Swing panels via `gui_river`/`gui_component`, `GlobalData` globals). Tests run via the MATLAB MCP tool `mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code` (function-style), NOT runtests.

**Spec:** `docs/superpowers/specs/2026-06-01-eigenmode-scale-lever-design.md`

---

## File Structure

**Create:**
- `toolbox/gui/panel_eigenmodes.m` — the lever: `CreatePanel`, state verbs (`SetCurrentMode`/`SetBand`/`SetWindowShape`/`SetActive`/`GetWeights`/`UpdatePanel`), the pure weight builder `BuildWeights`, the reconstruction entry `ApplyToColumn`, and slider/radio callbacks. All lever logic lives here.
- `dev/tests/test_eigenmode_lever_weights.m` — pure-math weight-shape + round-trip tests.
- `dev/tests/test_eigenmode_lever_state.m` — state clamp/coupling/no-op tests.
- `dev/tests/test_eigenmode_lever_integration.m` — headless end-to-end on a synthetic mesh.

**Modify:**
- `toolbox/gui/panel_surface.m` — one guarded line in `UpdateSurfaceData` (after the Results branch sets `TessInfo(iTess).Data`) calling `panel_eigenmodes('ApplyToColumn', ...)`.
- `toolbox/core/bst_figures.m` — add `FireModesChanged` subfunction; call `panel_eigenmodes('UpdatePanel', hNewFig)` in `SetCurrentFigure`'s 3D branch.
- `toolbox/gui/gui_initialize.m` — register the panel: `gui_show('panel_eigenmodes', 'BrainstormTab', 'tools')`.

## Key API facts (verified — do not re-guess)

- `[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile)` → `Eig.Vectors [nV x K]`, `Eig.Values [K x 1]`, `Eig.MassType`, `Eig.nModes`.
- `[~, MassMatrix] = tess_laplacian(Vertices, Faces, 'MassType', Eig.MassType)` → sparse `[nV x nV]`.
- `Filtered = bst_eigenmodes_filter(Eig, Data, MassMatrix, 'custom', 'TransferFn', fn)` where `fn(lambdas)` returns a `[nModes x 1]` gain. We pass `fn = @(l) W(:)` with our weight vector `W`.
- `TessInfo(iTess).Data` is the `[nV x 1]` displayed column; for source maps it is set at `panel_surface.m:~1813` via `bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex')`. `TessInfo(iTess).DataSource.Type` is `'Source'`/`'Results'` for source overlays; `TessInfo(iTess).SurfaceFile` is the surface.
- `bst_figures` and panels dispatch subfunctions via `eval(macro_method)` (call as `bst_figures('FireModesChanged')`, `panel_eigenmodes('SetBand', ...)`).
- Panels register at startup with `gui_show('panel_NAME', 'BrainstormTab', 'tools')` (see `gui_initialize.m:49-52`).

---

## Task 1: Pure weight-shape builder

The heart of the lever: build the `[1 x K]` weight vector from `(shape, kLo, kHi, iCenter)`. Fully headless, no GlobalData, no GUI. Four shapes, all defined purely in mode-index space (no eigenvalues needed):
- `single` → delta at `iCenter`
- `box` → 1 on `[kLo, kHi]`
- `tapered` → Tukey window on `[kLo, kHi]` (flat interior, cosine shoulders) — anti-ringing
- `gain` → Gaussian bell centered at `iCenter`, sigma = half the band width — smooth gain curve

**Note (spec refinement):** the spec listed `gain` as "heat-kernel via bst_eigenmodes_filter_gain". We implement `gain` as a Gaussian-in-mode-index instead — it needs no eigenvalues, keeps `BuildWeights` pure/headless-testable, and still delivers the smooth anti-ringing rolloff intended. Reconstruction still reuses `bst_eigenmodes_filter` (Task 3); only the weight *shape* is authored locally.

**Files:**
- Create: `toolbox/gui/panel_eigenmodes.m`
- Create: `dev/tests/test_eigenmode_lever_weights.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmode_lever_weights.m`:

```matlab
function test_eigenmode_lever_weights
% Pure-math tests for panel_eigenmodes('BuildWeights', ...).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

K = 200;

% single -> delta at center
w = panel_eigenmodes('BuildWeights', 'single', 1, K, 42, K);
assert(isequal(size(w), [1 K]), 'weights must be 1xK row');
assert(w(42) == 1 && sum(w) == 1, 'single must be a delta at the center');

% box -> 1 on [kLo,kHi], 0 elsewhere
w = panel_eigenmodes('BuildWeights', 'box', 30, 55, 42, K);
assert(all(w(30:55) == 1), 'box interior must be 1');
assert(sum(w) == 26, 'box must keep exactly 26 modes');
assert(w(29) == 0 && w(56) == 0, 'box must be 0 outside the band');

% tapered -> <=1 everywhere, interior 1, monotone shoulders, 0 outside
w = panel_eigenmodes('BuildWeights', 'tapered', 30, 55, 42, K);
assert(all(w >= 0 & w <= 1), 'tapered weights in [0,1]');
assert(w(29) == 0 && w(56) == 0, 'tapered 0 outside band');
mid = round((30+55)/2);
assert(abs(w(mid) - 1) < 1e-9, 'tapered interior reaches 1');
assert(w(30) <= w(31) && w(31) <= w(mid), 'tapered rising shoulder is monotone');

% gain -> Gaussian bell, peak at center, symmetric falloff, in (0,1]
w = panel_eigenmodes('BuildWeights', 'gain', 30, 55, 42, K);
assert(all(w >= 0 & w <= 1 + 1e-12), 'gain weights in [0,1]');
[~, ipk] = max(w);
assert(ipk == 42, 'gain peak at the center mode');
assert(w(42) > w(20) && w(42) > w(64), 'gain falls off away from center');

% clamping inside the builder: out-of-range band is clamped to [1,K]
w = panel_eigenmodes('BuildWeights', 'box', -5, K+100, 42, K);
assert(all(w == 1), 'box clamped to full range keeps all modes');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_weights\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `evaluate_matlab_code`:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_weights.m')
```
Expected: FAIL — `Undefined function 'panel_eigenmodes'` (file not created yet).

- [ ] **Step 3: Create panel_eigenmodes.m with the dispatch header and BuildWeights**

Create `toolbox/gui/panel_eigenmodes.m`:

```matlab
function varargout = panel_eigenmodes(varargin)
% PANEL_EIGENMODES: Eigenmode scale lever — a stateful, surface-scoped,
% live-broadcast selection in eigenmode (spatial-frequency) space. Changing the
% selection live-coarsens/smooths a displayed source map via a non-destructive
% display-time band-limited reconstruction.
%
% USAGE:  bstPanel = panel_eigenmodes('CreatePanel')
%         W  = panel_eigenmodes('BuildWeights', shape, kLo, kHi, iCenter, K)
%         panel_eigenmodes('SetBand', kLo, kHi)
%         panel_eigenmodes('SetCurrentMode', k)
%         panel_eigenmodes('SetWindowShape', shape)
%         panel_eigenmodes('SetActive', isActive)
%         W  = panel_eigenmodes('GetWeights')
%         uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u)
%         panel_eigenmodes('UpdatePanel', hFig)
%
% SEE ALSO: bst_eigenmodes_filter, in_tess_eigenmodes, panel_freq, panel_surface

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== PURE: build the [1 x K] weight vector from a window shape =====
function W = BuildWeights(shape, kLo, kHi, iCenter, K) %#ok<DEFNU>
    % Clamp all indices into [1,K]
    kLo     = min(max(round(kLo),     1), K);
    kHi     = min(max(round(kHi),     1), K);
    iCenter = min(max(round(iCenter), 1), K);
    if (kHi < kLo)
        [kLo, kHi] = deal(kHi, kLo);
    end
    W = zeros(1, K);
    switch lower(shape)
        case 'single'
            W(iCenter) = 1;
        case 'box'
            W(kLo:kHi) = 1;
        case 'tapered'
            % Tukey window over [kLo,kHi]: flat interior, cosine shoulders.
            n = kHi - kLo + 1;
            if (n <= 1)
                W(kLo:kHi) = 1;
            else
                r = 0.5;                       % shoulder fraction (each side)
                t = linspace(0, 1, n);
                wt = ones(1, n);
                edge = (r/2);
                iL = t < edge;
                iR = t > (1 - edge);
                wt(iL) = 0.5 * (1 + cos(pi * (2*t(iL)/r - 1)));
                wt(iR) = 0.5 * (1 + cos(pi * (2*t(iR)/r - 2/r + 1)));
                W(kLo:kHi) = wt;
            end
        case 'gain'
            % Gaussian bell centered at iCenter; sigma = half the band width.
            sigma = max((kHi - kLo) / 2, 1);
            k = 1:K;
            W = exp(-0.5 * ((k - iCenter) / sigma).^2);
        otherwise
            error('panel_eigenmodes:BuildWeights: unknown shape "%s".', shape);
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_weights.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_weights`

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_weights.m
git commit -m "Eigenmode lever: pure weight-shape builder (single/box/tapered/gain)"
```

---

## Task 2: State + clamping + coupled center/band

Add the `GlobalData.UserModes` state and its setter verbs. State is lazily initialized; setters clamp to `[1,K]`, keep the band center coupled to `iCurrentMode`, recompute `Weights`, and (in later tasks) broadcast. In this task the setters update state only; the broadcast call is added in Task 5.

**Files:**
- Modify: `toolbox/gui/panel_eigenmodes.m`
- Create: `dev/tests/test_eigenmode_lever_state.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmode_lever_state.m`:

```matlab
function test_eigenmode_lever_state
% State logic: init, clamp, coupled center/band, isActive default.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
global GlobalData; %#ok<GVMIS>

% Reset state for a 200-mode surface
panel_eigenmodes('ResetState', '/fake/surf.mat', 200);
st = GlobalData.UserModes;
assert(strcmp(st.SurfaceFile, '/fake/surf.mat'), 'surface stored');
assert(st.nModes == 200, 'K stored');
assert(st.isActive == 0, 'lever starts inactive');

% SetBand clamps and recenters iCurrentMode to band center
panel_eigenmodes('SetBand', 30, 55);
st = GlobalData.UserModes;
assert(isequal(st.Band, [30 55]), 'band stored');
assert(st.iCurrentMode == round((30+55)/2), 'center coupled to band midpoint');
assert(isequal(size(st.Weights), [1 200]), 'weights recomputed [1xK]');
assert(sum(st.Weights) == 26, 'box default keeps 26 modes');

% SetBand clamps out-of-range
panel_eigenmodes('SetBand', -10, 9999);
st = GlobalData.UserModes;
assert(isequal(st.Band, [1 200]), 'band clamped to [1,K]');

% SetCurrentMode (coupled): slides the band, preserving its width
panel_eigenmodes('SetBand', 30, 50);   % width 21, center 40
panel_eigenmodes('SetCurrentMode', 100);
st = GlobalData.UserModes;
assert(st.iCurrentMode == 100, 'center moved');
assert((st.Band(2) - st.Band(1)) == 20, 'band width preserved when sliding center');
assert(st.Band(1) == 90 && st.Band(2) == 110, 'band slid to recentre on 100');

% SetWindowShape 'single' collapses band to the center
panel_eigenmodes('SetWindowShape', 'single');
st = GlobalData.UserModes;
assert(strcmp(st.WindowShape, 'single'), 'shape stored');
assert(sum(st.Weights) == 1 && st.Weights(100) == 1, 'single -> delta at center');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_state\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_state.m')
```
Expected: FAIL — `Undefined function or variable 'macro_method'`-style error or `Unknown command 'ResetState'` (verbs not implemented).

- [ ] **Step 3: Add state init + setters to panel_eigenmodes.m**

Add these subfunctions to `toolbox/gui/panel_eigenmodes.m` (after `BuildWeights`):

```matlab
%% ===== STATE: lazy default =====
function st = GetState()
    global GlobalData;
    if ~isfield(GlobalData, 'UserModes') || isempty(GlobalData.UserModes) ...
            || ~isfield(GlobalData.UserModes, 'Weights')
        GlobalData.UserModes = struct(...
            'SurfaceFile',  '', ...
            'nModes',       0,  ...
            'iCurrentMode', 1,  ...
            'Weights',      [], ...
            'WindowShape',  'box', ...
            'Band',         [1 1], ...
            'isActive',     0,  ...
            'CacheSurfaceFile', '', ...
            'CacheEig',     [], ...
            'CacheMass',    []);
    end
    st = GlobalData.UserModes;
end

%% ===== STATE: reset for a surface with K modes =====
function ResetState(SurfaceFile, K) %#ok<DEFNU>
    global GlobalData;
    GetState();                          % ensure struct exists
    GlobalData.UserModes.SurfaceFile  = SurfaceFile;
    GlobalData.UserModes.nModes       = K;
    GlobalData.UserModes.Band         = [1, min(30, K)];
    GlobalData.UserModes.iCurrentMode = round(mean(GlobalData.UserModes.Band));
    GlobalData.UserModes.WindowShape  = 'box';
    GlobalData.UserModes.isActive     = 0;
    RecomputeWeights();
end

%% ===== STATE: recompute Weights from shape/band/center =====
function RecomputeWeights()
    global GlobalData;
    st = GlobalData.UserModes;
    GlobalData.UserModes.Weights = BuildWeights(st.WindowShape, ...
        st.Band(1), st.Band(2), st.iCurrentMode, st.nModes);
end

%% ===== STATE: set band (clamped); recentre iCurrentMode to band midpoint =====
function SetBand(kLo, kHi) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    kLo = min(max(round(kLo), 1), K);
    kHi = min(max(round(kHi), 1), K);
    if (kHi < kLo), [kLo, kHi] = deal(kHi, kLo); end
    GlobalData.UserModes.Band         = [kLo, kHi];
    GlobalData.UserModes.iCurrentMode = round((kLo + kHi) / 2);
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: move center (coupled — slide band, preserve width) =====
function SetCurrentMode(k) %#ok<DEFNU>
    global GlobalData;
    st = GetState();
    K  = st.nModes;
    k  = min(max(round(k), 1), K);
    halfW = round((st.Band(2) - st.Band(1)) / 2);
    kLo = min(max(k - halfW, 1), K);
    kHi = min(max(k + halfW, 1), K);
    GlobalData.UserModes.iCurrentMode = k;
    GlobalData.UserModes.Band         = [kLo, kHi];
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: window shape ('single' collapses band to center) =====
function SetWindowShape(shape) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.WindowShape = lower(shape);
    if strcmpi(shape, 'single')
        c = GlobalData.UserModes.iCurrentMode;
        GlobalData.UserModes.Band = [c, c];
    end
    RecomputeWeights();
    NotifyChanged();
end

%% ===== STATE: active toggle =====
function SetActive(isActive) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.isActive = logical(isActive);
    NotifyChanged();
end

%% ===== STATE: read canonical weights =====
function W = GetWeights() %#ok<DEFNU>
    st = GetState();
    W = st.Weights;
end

%% ===== Broadcast (real implementation added in Task 5) =====
function NotifyChanged()
    % Placeholder until Task 5 wires the figure repaint broadcast.
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_state.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_state`

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_state.m
git commit -m "Eigenmode lever: surface-scoped state, clamping, coupled center/band"
```

---

## Task 3: ApplyToColumn — guarded non-destructive reconstruction

The single entry point the surface display calls. Given a surface file and a `[nV x 1]` source column `u`, return the band-limited reconstruction when the lever is active and matches; otherwise return `u` unchanged. Caches the eigenmodes + mass matrix per surface (recomputing the mass matrix every time would be slow). Reconstruction reuses `bst_eigenmodes_filter(...,'custom','TransferFn',...)`.

**Files:**
- Modify: `toolbox/gui/panel_eigenmodes.m`
- Modify: `dev/tests/test_eigenmode_lever_integration.m` (created here)

- [ ] **Step 1: Write the failing test (synthetic mesh, no DB)**

Create `dev/tests/test_eigenmode_lever_integration.m`:

```matlab
function test_eigenmode_lever_integration
% Headless end-to-end on a synthetic mesh: ApplyToColumn must equal the
% analytic band-limited reconstruction, be a no-op when inactive, and a
% full-band weight must round-trip (uF == u).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(fullfile(repoRoot, 'toolbox', 'gui'));
addpath(fullfile(repoRoot, 'toolbox', 'math'));
addpath(fullfile(repoRoot, 'toolbox', 'anatomy'));
global GlobalData; %#ok<GVMIS>

% --- Build a small closed mesh (subdivided icosphere) ---
[V, F] = tess_sphere(642);          % ~642-vertex sphere, genus-0
K = 60;
[Vectors, Values] = tess_eigenmodes(V, F, 'nModes', K, 'MassType', 'barycentric', ...
                                    'RemoveDC', 1, 'Verbose', 0);
Eig = struct('Vectors', Vectors, 'Values', Values, 'nModes', K, 'MassType', 'barycentric');
[~, M] = tess_laplacian(V, F, 'MassType', 'barycentric');

% Inject the cache directly (ApplyToColumn uses it when SurfaceFile matches)
SurfaceFile = '/synthetic/sphere.mat';
panel_eigenmodes('ResetState', SurfaceFile, K);
panel_eigenmodes('SetCache', SurfaceFile, Eig, M);

% A random scalar field on the mesh
rng(1); u = randn(size(V,1), 1);

% --- Inactive => no-op ---
panel_eigenmodes('SetActive', 0);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
assert(isequal(uF, u), 'inactive lever must be a no-op');

% --- Active, full-band box (1..K) => round-trip identity ---
panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, K);
uF = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
assert(max(abs(uF - u)) < 1e-6 * max(abs(u)), 'full band must round-trip to u');

% --- Active, low band => smoother than the original (less high-mode energy) ---
panel_eigenmodes('SetBand', 1, 10);
uLow = panel_eigenmodes('ApplyToColumn', SurfaceFile, u);
W = panel_eigenmodes('GetWeights');
analytic = Vectors * (W(:) .* (Vectors' * (M * u)));
assert(max(abs(uLow - analytic)) < 1e-9, 'ApplyToColumn must equal analytic reconstruction');
assert(norm(uLow) < norm(u), 'low-band reconstruction must shrink overall energy');

% --- Surface mismatch => no-op (different file than cache/state) ---
uOther = panel_eigenmodes('ApplyToColumn', '/synthetic/other.mat', u);
assert(isequal(uOther, u), 'mismatched surface must be a no-op');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_integration\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
```
Expected: FAIL — `Unknown command 'SetCache'` / `ApplyToColumn` not implemented.

- [ ] **Step 3: Add SetCache + ApplyToColumn to panel_eigenmodes.m**

Add to `toolbox/gui/panel_eigenmodes.m`:

```matlab
%% ===== CACHE: store eigenmodes + mass matrix for a surface =====
function SetCache(SurfaceFile, Eig, MassMatrix) %#ok<DEFNU>
    global GlobalData;
    GetState();
    GlobalData.UserModes.CacheSurfaceFile = SurfaceFile;
    GlobalData.UserModes.CacheEig         = Eig;
    GlobalData.UserModes.CacheMass        = MassMatrix;
end

%% ===== CACHE: ensure eigenmodes + mass are loaded for a surface =====
function isOk = EnsureCache(SurfaceFile)
    global GlobalData;
    st = GetState();
    isOk = false;
    if ~isempty(st.CacheEig) && ~isempty(st.CacheMass) ...
            && file_compare(st.CacheSurfaceFile, SurfaceFile)
        isOk = true;
        return;
    end
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        return;
    end
    sSurf = in_tess_bst(SurfaceFile, 0);
    [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
    SetCache(SurfaceFile, Eig, M);
    isOk = true;
end

%% ===== APPLY: filter a displayed source column (guarded, non-destructive) =====
function uF = ApplyToColumn(SurfaceFile, u) %#ok<DEFNU>
    global GlobalData;
    uF = u;                                   % default: unchanged
    st = GetState();
    % Guards: lever off, surface mismatch, empty column
    if ~st.isActive || isempty(u) || isempty(SurfaceFile) ...
            || ~file_compare(st.SurfaceFile, SurfaceFile)
        return;
    end
    if ~EnsureCache(SurfaceFile)
        return;
    end
    Eig = GlobalData.UserModes.CacheEig;
    M   = GlobalData.UserModes.CacheMass;
    % Scalar-field guard: only filter when the column matches the mesh vertex
    % count (skip unconstrained/volume/mismatched maps rather than mis-filter).
    if (size(u,1) ~= size(Eig.Vectors,1)) || (size(u,2) ~= 1)
        return;
    end
    W = GlobalData.UserModes.Weights;
    if isempty(W) || (numel(W) ~= Eig.nModes)
        return;
    end
    % Reconstruct via the core spectral filter (custom transfer = our weights)
    uF = bst_eigenmodes_filter(Eig, u, M, 'custom', 'TransferFn', @(l) W(:));
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_integration`

If `tess_sphere` has a different signature, substitute any small closed mesh available in `toolbox/anatomy` (verify with `help tess_sphere` first); the test only needs a genus-0 mesh of a few hundred vertices.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_integration.m
git commit -m "Eigenmode lever: ApplyToColumn guarded reconstruction + per-surface cache"
```

---

## Task 4: Panel UI (CreatePanel) — band slider + readout

Build the docked panel in `panel_freq` house style: an **Active** checkbox, a **dual-handle band** rendered as two sliders (low / high) with a center readout, **window-shape** radio buttons, and a live **readout** label (kept-mode count + λ range). Callbacks fire on release into the Task 2/3 verbs. `UpdatePanel(hFig)` enables/populates from the current figure's surface or disables when ineligible.

(Brainstorm has no native dual-handle JSlider; two coupled sliders — `jSliderLo`, `jSliderHi` — are the idiomatic substitute, matching how `panel_freq` uses a single `JSlider`.)

**Files:**
- Modify: `toolbox/gui/panel_eigenmodes.m`

- [ ] **Step 1: Write the failing smoke test**

Append to `dev/tests/test_eigenmode_lever_state.m` a new function call is not possible (one function per file), so create `dev/tests/test_eigenmode_lever_panel.m`:

```matlab
function test_eigenmode_lever_panel
% Smoke test: the panel builds, and UpdatePanel enables/disables correctly.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

bstPanel = panel_eigenmodes('CreatePanel');
assert(~isempty(bstPanel), 'CreatePanel must return a BstPanel');

% Controls are registered under the panel name
ctrl = bst_get('PanelControls', 'EigenModes');
assert(~isempty(ctrl), 'panel controls must be registered');
assert(isfield(ctrl, 'jCheckActive') && isfield(ctrl, 'jSliderLo') ...
       && isfield(ctrl, 'jSliderHi') && isfield(ctrl, 'jLabelReadout'), ...
       'expected controls present');

fprintf('ALL TESTS PASSED: test_eigenmode_lever_panel\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_panel.m')
```
Expected: FAIL — `Unknown command 'CreatePanel'`.

- [ ] **Step 3: Implement CreatePanel + UpdatePanel + callbacks**

Add to `toolbox/gui/panel_eigenmodes.m`:

```matlab
%% ===== CREATE PANEL =====
function bstPanelNew = CreatePanel() %#ok<DEFNU>
    panelName = 'EigenModes';
    import java.awt.*;
    import javax.swing.*;

    jPanelNew = gui_river([2,2], [4,4,6,6], 'Spatial scale (eigenmodes)');

    % Active toggle + readout (same row)
    jCheckActive = gui_component('CheckBox', jPanelNew, '', 'Active', [], ...
        'Live-filter the displayed source map', @(h,ev)CheckActive_Callback());
    jLabelReadout = gui_component('Label', jPanelNew, 'hfill', '');
    jLabelReadout.setHorizontalAlignment(JLabel.RIGHT);

    % Band: low / high sliders (dual-handle substitute)
    gui_component('Label', jPanelNew, 'br', 'Mode band');
    jSliderLo = JSlider(1, 100, 1);
    jSliderHi = JSlider(1, 100, 30);
    java_setcb(jSliderLo, 'MouseReleasedCallback', @(h,ev)Slider_Callback());
    java_setcb(jSliderHi, 'MouseReleasedCallback', @(h,ev)Slider_Callback());
    jPanelNew.add('br hfill', jSliderLo);
    jPanelNew.add('br hfill', jSliderHi);
    jLabelBand = gui_component('Label', jPanelNew, 'br', 'lo=1  c=15  hi=30');

    % Window shape radios
    jGroup = ButtonGroup();
    gui_component('Label', jPanelNew, 'br', 'Window:');
    jRadioSingle = gui_component('Radio', jPanelNew, '',   'Single', jGroup, '', @(h,ev)Shape_Callback('single'));
    jRadioBox    = gui_component('Radio', jPanelNew, '',   'Box',    jGroup, '', @(h,ev)Shape_Callback('box'));
    jRadioTaper  = gui_component('Radio', jPanelNew, 'br', 'Taper',  jGroup, '', @(h,ev)Shape_Callback('tapered'));
    jRadioGain   = gui_component('Radio', jPanelNew, '',   'Gain',   jGroup, '', @(h,ev)Shape_Callback('gain'));
    jRadioBox.setSelected(1);

    ctrl = struct('jPanelTop',      jPanelNew, ...
                  'jCheckActive',   jCheckActive, ...
                  'jLabelReadout',  jLabelReadout, ...
                  'jSliderLo',      jSliderLo, ...
                  'jSliderHi',      jSliderHi, ...
                  'jLabelBand',     jLabelBand, ...
                  'jRadioSingle',   jRadioSingle, ...
                  'jRadioBox',      jRadioBox, ...
                  'jRadioTaper',    jRadioTaper, ...
                  'jRadioGain',     jRadioGain);
    bstPanelNew = BstPanel(panelName, jPanelNew, ctrl);
end

%% ===== CALLBACKS (read controls -> state verbs) =====
function CheckActive_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    SetActive(ctrl.jCheckActive.isSelected());
end

function Slider_Callback()
    ctrl = bst_get('PanelControls', 'EigenModes');
    SetBand(ctrl.jSliderLo.getValue(), ctrl.jSliderHi.getValue());
    RefreshControls();
end

function Shape_Callback(shape)
    SetWindowShape(shape);
    RefreshControls();
end

%% ===== UPDATE PANEL: populate/enable from the active figure's surface =====
function UpdatePanel(hFig) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', 'EigenModes');
    if isempty(ctrl)
        return;
    end
    SurfaceFile = GetFigureSurfaceWithModes(hFig);
    isEligible = ~isempty(SurfaceFile);
    SetPanelEnabled(ctrl, isEligible);
    if ~isEligible
        ctrl.jLabelReadout.setText('no eigenmodes');
        return;
    end
    [Eig, ~] = in_tess_eigenmodes(SurfaceFile);
    K = Eig.nModes;
    % (Re)initialise state if the surface changed
    st = GetState();
    if ~file_compare(st.SurfaceFile, SurfaceFile) || (st.nModes ~= K)
        ResetState(SurfaceFile, K);
    end
    ctrl.jSliderLo.setMaximum(K);  ctrl.jSliderHi.setMaximum(K);
    RefreshControls();
end

%% ===== Reflect state back into the controls + readout =====
function RefreshControls()
    ctrl = bst_get('PanelControls', 'EigenModes');
    st = GetState();
    ctrl.jSliderLo.setValue(st.Band(1));
    ctrl.jSliderHi.setValue(st.Band(2));
    ctrl.jLabelBand.setText(sprintf('lo=%d  c=%d  hi=%d', st.Band(1), st.iCurrentMode, st.Band(2)));
    nKeep = nnz(st.Weights > 1e-6);
    lamStr = '';
    if ~isempty(st.CacheEig) && ~isempty(st.CacheEig.Values)
        lam = st.CacheEig.Values;
        b = min(max(st.Band, 1), numel(lam));
        lamStr = sprintf('  lambda in [%.3g, %.3g]', lam(b(1)), lam(b(2)));
    end
    ctrl.jLabelReadout.setText(sprintf('modes %d-%d  (%d)%s', st.Band(1), st.Band(2), nKeep, lamStr));
end

%% ===== helpers =====
function SetPanelEnabled(ctrl, isOn)
    fn = fieldnames(ctrl);
    for i = 1:numel(fn)
        c = ctrl.(fn{i});
        if isa(c, 'javax.swing.JComponent')
            c.setEnabled(logical(isOn));
        end
    end
end

function SurfaceFile = GetFigureSurfaceWithModes(hFig)
    SurfaceFile = '';
    if isempty(hFig) || ~ishandle(hFig)
        return;
    end
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo)
        return;
    end
    for iTess = 1:numel(TessInfo)
        sf = TessInfo(iTess).SurfaceFile;
        if isempty(sf) || isempty(TessInfo(iTess).DataSource) ...
                || isempty(TessInfo(iTess).DataSource.FileName)
            continue;
        end
        [~, isComputed] = in_tess_eigenmodes(sf);
        if isComputed
            SurfaceFile = sf;
            return;
        end
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_panel.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_panel`

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_eigenmodes.m dev/tests/test_eigenmode_lever_panel.m
git commit -m "Eigenmode lever: panel UI (band sliders + window radios + readout)"
```

---

## Task 5: Broadcast (FireModesChanged) + figure-focus hook

Replace the `NotifyChanged` placeholder with a real broadcast, and add `FireModesChanged` to `bst_figures` (repaints visible 3D figures on the matching surface). Wire `panel_eigenmodes('UpdatePanel', hNewFig)` into `SetCurrentFigure`'s 3D branch so the panel tracks the active figure.

**Files:**
- Modify: `toolbox/gui/panel_eigenmodes.m`
- Modify: `toolbox/core/bst_figures.m`

- [ ] **Step 1: Add FireModesChanged to bst_figures.m**

Add immediately after the `FireCurrentTimeChanged` function (ends ~line 990 in `toolbox/core/bst_figures.m`):

```matlab
%% ===== FIRE MODES CHANGED =====
% Repaint visible 3D source figures after the eigenmode lever changes.
function FireModesChanged() %#ok<DEFNU>
    global GlobalData;
    if ~isfield(GlobalData, 'UserModes') || isempty(GlobalData.UserModes)
        return;
    end
    SurfaceFile = GlobalData.UserModes.SurfaceFile;
    for iDS = 1:length(GlobalData.DataSet)
        for iFig = 1:length(GlobalData.DataSet(iDS).Figure)
            sFig = GlobalData.DataSet(iDS).Figure(iFig);
            if strcmpi(get(sFig.hFigure, 'Visible'), 'off') || ~strcmpi(sFig.Id.Type, '3DViz')
                continue;
            end
            % Only repaint figures showing the lever's surface
            TessInfo = getappdata(sFig.hFigure, 'Surface');
            isMatch = false;
            for iTess = 1:numel(TessInfo)
                if file_compare(TessInfo(iTess).SurfaceFile, SurfaceFile)
                    isMatch = true; break;
                end
            end
            if ~isMatch
                continue;
            end
            % UpdateSurfaceData recomputes the column (filtered via ApplyToColumn),
            % then UpdateSurfaceColor repaints.
            panel_surface('UpdateSurfaceData', sFig.hFigure);
            figure_3d('UpdateSurfaceColor', sFig.hFigure);
        end
    end
end
```

- [ ] **Step 2: Wire NotifyChanged in panel_eigenmodes.m**

Replace the placeholder `NotifyChanged` in `toolbox/gui/panel_eigenmodes.m`:

```matlab
%% ===== Broadcast: repaint affected figures =====
function NotifyChanged()
    bst_figures('FireModesChanged');
end
```

- [ ] **Step 3: Add the focus hook in SetCurrentFigure**

In `toolbox/core/bst_figures.m`, in `SetCurrentFigure`'s 3D branch, immediately after the existing `panel_surface('UpdatePanel');` line (~1393):

```matlab
        % Update Surfaces panel
        panel_surface('UpdatePanel');
        % Update the eigenmode scale lever for the new figure's surface
        if gui_brainstorm('isTabVisible', 'EigenModes')
            panel_eigenmodes('UpdatePanel', hNewFig);
        end
```

- [ ] **Step 4: Verify no syntax regressions**

Run:
```matlab
checkcode('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/core/bst_figures.m')
checkcode('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/panel_eigenmodes.m')
```
Expected: no errors (warnings about `%#ok` unused are fine). Then re-run the state test to confirm the broadcast call does not break headless state updates (FireModesChanged returns early when no figures):
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_state.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_state`

- [ ] **Step 5: Commit**

```bash
git add toolbox/core/bst_figures.m toolbox/gui/panel_eigenmodes.m
git commit -m "Eigenmode lever: FireModesChanged broadcast + figure-focus UpdatePanel hook"
```

---

## Task 6: Surface-display hook (the one-line chokepoint)

Make `panel_surface('UpdateSurfaceData')` route the source column through the lever. This is the single integration point that makes filtering live and time-composing.

**Files:**
- Modify: `toolbox/gui/panel_surface.m`

- [ ] **Step 1: Locate the Results branch**

In `toolbox/gui/panel_surface.m`, find (≈ line 1813):
```matlab
                TessInfo(iTess).Data = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
```

- [ ] **Step 2: Add the guarded filter call directly after it**

```matlab
                TessInfo(iTess).Data = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex');
                % Eigenmode scale lever: live, non-destructive band-limited
                % reconstruction of the displayed column (no-op when inactive).
                TessInfo(iTess).Data = panel_eigenmodes('ApplyToColumn', ...
                    TessInfo(iTess).SurfaceFile, TessInfo(iTess).Data);
```

- [ ] **Step 3: Verify no syntax regressions**

Run:
```matlab
checkcode('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/panel_surface.m')
```
Expected: no new errors.

- [ ] **Step 4: Re-run the reconstruction test (the hook calls the same ApplyToColumn)**

Run:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_lever_integration.m')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_integration` (unchanged — ApplyToColumn semantics are what the hook relies on).

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/panel_surface.m
git commit -m "Eigenmode lever: route source display column through ApplyToColumn"
```

---

## Task 7: Register the panel at startup

**Files:**
- Modify: `toolbox/gui/gui_initialize.m`

- [ ] **Step 1: Add the registration line**

In `toolbox/gui/gui_initialize.m`, after the `panel_scout` registration (≈ line 52):
```matlab
gui_show('panel_scout', 'BrainstormTab', 'tools');
gui_show('panel_eigenmodes', 'BrainstormTab', 'tools');
```

- [ ] **Step 2: Verify the panel registers on a fresh start**

Run:
```matlab
brainstorm stop; brainstorm nogui
ctrl = bst_get('PanelControls', 'EigenModes');
fprintf('EigenModes registered: %d\n', ~isempty(ctrl));
```
Expected: `EigenModes registered: 1`

- [ ] **Step 3: Commit**

```bash
git add toolbox/gui/gui_initialize.m
git commit -m "Eigenmode lever: register panel_eigenmodes tool tab at startup"
```

---

## Task 8: End-to-end test on a real subject (MCP, interactive MATLAB)

Validate the full live path against a real Brainstorm protocol: compute eigenmodes, display a source map, activate the lever, drag the band, and assert the rendered surface data is the band-limited reconstruction and is smoother; toggling Active off restores the raw map.

**Files:**
- Create: `dev/tests/test_eigenmode_lever_e2e.m`

- [ ] **Step 1: Write the end-to-end test**

Create `dev/tests/test_eigenmode_lever_e2e.m`:

```matlab
function test_eigenmode_lever_e2e(SurfaceFile, ResultsFile)
% End-to-end live test. Pass a cortex SurfaceFile (with or without eigenmodes)
% and a source ResultsFile on that surface from the current protocol.
%
% USAGE: test_eigenmode_lever_e2e(SurfaceFile, ResultsFile)
if ~brainstorm('status'); brainstorm nogui; end

% Ensure eigenmodes exist on the surface
[~, isComputed] = in_tess_eigenmodes(SurfaceFile);
if ~isComputed
    process_eigenmodes('Compute', SurfaceFile, 200, 'barycentric', true, false, false);
end

% Display the source map on the cortex
hFig = view_surface_data(SurfaceFile, ResultsFile);
assert(~isempty(hFig), 'source map failed to display');

% Bring up + populate the lever, activate it
gui_brainstorm('ShowToolTab', 'EigenModes');
panel_eigenmodes('UpdatePanel', hFig);
panel_eigenmodes('SetActive', 1);

% Raw column (lever inactive) for comparison
panel_eigenmodes('SetActive', 0);
panel_surface('UpdateSurfaceData', hFig);
TessInfo = getappdata(hFig, 'Surface');
uRaw = TessInfo(1).Data;

% Low band => smoother; values must match analytic reconstruction
panel_eigenmodes('SetActive', 1);
panel_eigenmodes('SetWindowShape', 'box');
panel_eigenmodes('SetBand', 1, 15);     % fires FireModesChanged -> repaint
TessInfo = getappdata(hFig, 'Surface');
uLow = TessInfo(1).Data;
assert(norm(uLow) <= norm(uRaw) + 1e-9, 'low band must not increase energy');
assert(~isequal(uLow, uRaw), 'low band must change the displayed column');

% Toggle off restores raw
panel_eigenmodes('SetActive', 0);
panel_surface('UpdateSurfaceData', hFig);
TessInfo = getappdata(hFig, 'Surface');
assert(max(abs(TessInfo(1).Data - uRaw)) < 1e-9, 'deactivation must restore raw map');

close(hFig);
fprintf('ALL TESTS PASSED: test_eigenmode_lever_e2e\n');
end
```

- [ ] **Step 2: Run the end-to-end test**

Identify a cortex surface + a source results file in the current protocol (e.g. via `bst_get('ProtocolSubjects')` / the tree). Then run:
```matlab
test_eigenmode_lever_e2e('<SurfaceFile>', '<ResultsFile>')
```
Expected: `ALL TESTS PASSED: test_eigenmode_lever_e2e`. If no protocol data is available in the session, document this and run the headless `test_eigenmode_lever_integration` as the substitute gate.

- [ ] **Step 3: Run the full lever test suite**

```matlab
run('.../dev/tests/test_eigenmode_lever_weights.m')
run('.../dev/tests/test_eigenmode_lever_state.m')
run('.../dev/tests/test_eigenmode_lever_integration.m')
run('.../dev/tests/test_eigenmode_lever_panel.m')
```
Expected: all four print `ALL TESTS PASSED`.

- [ ] **Step 4: Commit**

```bash
git add dev/tests/test_eigenmode_lever_e2e.m
git commit -m "Eigenmode lever: end-to-end live test on a real source map"
```

---

## Done criteria

- Activating the lever on a displayed source map and dragging the band live-smooths the cortex; deactivating restores the raw map exactly.
- The stored results file is never modified (display-time only).
- All five test files pass headlessly; the e2e test passes against a real protocol (or is documented as data-gated).
- `panel_eigenmodes` owns all lever logic; the footprint in shared files is one line in `panel_surface`, one function + one hook in `bst_figures`, one registration line in `gui_initialize`.
