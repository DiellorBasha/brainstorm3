# Eigenspectrum of Source Activations Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a "View eigenspectrum" context-menu item on source-results nodes that opens a dedicated 2D viewer showing the modal power spectrum (power per LBO eigenmode, left/right hemisphere curves, vs eigenvalue or spatial wavelength) of the realized vertex source map, driven by Brainstorm's global time cursor.

**Architecture:** A single new GUI file `view_eigenmode_spectrum.m` holds pure render helpers, a thin acquisition+projection core, an IO wrapper, and a Brainstorm-managed 2D figure. Acquisition reuses the canonical on-the-fly source path (`in_bst_results` LoadFull + `bst_source_orient` for unconstrained), projection reuses `bst_eigenmodes_project` (M from `tess_laplacian`). The figure registers as a new `Id.Type='EigenSpectrum'` in the source map's dataset; two one-line `case` additions to `bst_figures.m` route figure creation and global-time-change events to it, so it animates in lockstep with the cortex display. No results-structure, file-format, or existing-process change.

**Tech Stack:** MATLAB; Brainstorm toolbox idioms (`eval(macro_method)` dispatch, `gui_component` menus, `bst_figures`/`bst_memory`/`panel_time` GUI core); tests are `dev/tests/*.m` script-functions printing `ALL TESTS PASSED`, run via the MATLAB MCP `evaluate_matlab_code`.

---

## Reference facts (confirmed against the codebase — do not re-derive)

- **Projection (reuse):** `Coeffs = bst_eigenmodes_project(Eigenmodes, Data, MassMatrix)` returns `[nModes × nTime] = Φ'·(M·Data)`. `Eigenmodes.Vectors` is `[nVert × nModes]`. (`toolbox/math/bst_eigenmodes_project.m:71`)
- **Mass matrix (reuse):** `[L, M] = tess_laplacian(Vertices, Faces, 'MassType', MassType)`. We use only `M`. (`toolbox/anatomy/tess_eigenmodes.m:157`)
- **Eigenmodes load (reuse):** `[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile)` → `.Vectors`, `.Values`, `.Component` (1=Left, 2=Right), `.CompRank`, `.nComponents`, `.MassType`. (`toolbox/io/in_tess_eigenmodes.m`)
- **Unconstrained → scalar (reuse):** `S = bst_source_orient([], 3, [], ImageGridAmp, 'rms')`. The `'rms'` case computes `sqrt(x²+y²+z²)` (the vector magnitude shown on the cortex). Precedent: `process_source_flat.m:154`. (`toolbox/inverse/bst_source_orient.m:258-260`)
- **Realized source matrix (reuse):** `ResultsMat = in_bst_results(ResultsFile, 1)` applies the kernel and returns `.ImageGridAmp [nSrc × nTime]`, `.nComponents`, `.SurfaceFile`, `.Time`. (`toolbox/io/in_bst_results.m`)
- **Dataset for time-sync (reuse):** `[iDS, iResult] = bst_memory('LoadResultsFile', ResultsFile)`. (`toolbox/core/bst_memory.m:1068`)
- **Figure creation (reuse + extend):** `[hFig, iFig] = bst_figures('CreateFigure', iDS, FigureId, 'AlwaysCreate', ResultsFile)`. `CreateFigure` has a `switch(FigureId.Type)` with `otherwise→error` (`toolbox/core/bst_figures.m:166-198`) — a new type needs a `case` there. `GetChannelsForFigure` (called at `:217`) returns empty cleanly when `FigureId.Modality` is `''` (`:` first guard) — so set `Modality=''`.
- **Global time broadcast (extend):** `panel_time('SetCurrentTime', …)` → `bst_figures('FireCurrentTimeChanged')` (`panel_time.m:251`) → `switch sFig.Id.Type` (`bst_figures.m:949-979`) dispatches to each figure type's `'CurrentTimeChangedCallback'`. The cortex view is `'3DViz'` → `panel_surface('UpdateSurfaceData')` (`:954-955`). A new `case 'EigenSpectrum'` routes time changes to our redraw. `isStatic` figures are skipped (`:945`).
- **Time→index (reuse):** `iTime = bst_closest(GlobalData.UserTimeWindow.CurrentTime, Time)`. (`toolbox/math/bst_closest.m`)
- **Select figure (reuse):** `bst_figures('SetCurrentFigure', hFig, '2D')`. (`view_spectrum.m:209`)
- **Recover (iDS,iFig) from handle (reuse):** `[hFig, iFig, iDS] = bst_figures('GetFigure', hFig)`.
- **Menu site:** `tree_callbacks.m` `case {'results','link'}` at `:1796`; "Cortical activations" submenu `jMenuActivations` built at `:1818`; "Display on cortex" gated by `ismember(sStudy.Result(iResult).HeadModelType, {'surface','mixed'})` and `~isempty(sSubject.iCortex)` at `:1823-1825`. Variable `filenameRelative` is the results file.
- **Dispatch idiom (mirror):** `view_eigenmodes.m:31-35` routes named methods to subfunctions, else calls the main GUI entry.

## File structure

- **Create `toolbox/gui/view_eigenmode_spectrum.m`** — the whole feature except the two core hooks. Contents:
  - Dispatch header (named-method router, else `ViewFigure`).
  - **Pure helpers** (headlessly tested): `ComputeModalPower`, `GetSpectrumAxis`, `GetWindowAverage`.
  - **Acquisition core** (headlessly tested): `CollapseProject` (orientation-collapse + project).
  - **IO wrapper** (interactive): `GetActivationCoeffs(ResultsFile)`.
  - **GUI** (interactive): `ViewFigure`, `CreateFigure`, `UpdateFigurePlot`, `CurrentTimeChangedCallback`, `FigureKeyPressedCallback`.
- **Modify `toolbox/core/bst_figures.m`** — one `case 'EigenSpectrum'` in `CreateFigure` (`~:196`), one in `FireCurrentTimeChanged` (`~:969`).
- **Modify `toolbox/tree/tree_callbacks.m`** — one "View eigenspectrum" item in the `{'results','link'}` activations submenu (`~:1825`).
- **Create `dev/tests/test_view_eigenmode_spectrum_pure.m`** — pure render helpers.
- **Create `dev/tests/test_eigenmode_spectrum_acquire.m`** — `CollapseProject` constrained + unconstrained.

---

### Task 1: Pure render helpers + dispatch header

**Files:**
- Create: `toolbox/gui/view_eigenmode_spectrum.m`
- Test: `dev/tests/test_view_eigenmode_spectrum_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_view_eigenmode_spectrum_pure.m`:

```matlab
function test_view_eigenmode_spectrum_pure
% Verify the pure render helpers of view_eigenmode_spectrum.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

% ---- ComputeModalPower: |theta|^2 split by hemisphere component ----
Theta     = [3+4i; 1; 0; -2];          % use complex to prove magnitude-squared
Component = [1; 1; 2; 2];
pw = view_eigenmode_spectrum('ComputeModalPower', Theta, Component);
assert(isequal(size(pw.left),  [2 1]), 'left power wrong size');
assert(isequal(size(pw.right), [2 1]), 'right power wrong size');
assert(abs(pw.left(1)  - 25) < 1e-12, '|3+4i|^2 should be 25');   % 3^2+4^2
assert(abs(pw.left(2)  - 1)  < 1e-12, 'left(2) should be 1');
assert(abs(pw.right(1) - 0)  < 1e-12, 'right(1) should be 0');
assert(abs(pw.right(2) - 4)  < 1e-12, 'right(2) should be 4');    % (-2)^2

% ---- GetSpectrumAxis: eigenvalue passthrough, wavelength = 2pi/sqrt(lambda) ----
Values = [0; 1; (2*pi)^2];             % lambda=0 -> n/a; lambda=(2pi)^2 -> wavelength 1
axE = view_eigenmode_spectrum('GetSpectrumAxis', Values, 'eigenvalue');
assert(isequal(axE.x, Values), 'eigenvalue axis must pass values through');
assert(ischar(axE.label) && ~isempty(axE.label), 'eigenvalue label missing');
axW = view_eigenmode_spectrum('GetSpectrumAxis', Values, 'wavelength');
assert(isnan(axW.x(1)), 'lambda<=0 wavelength must be NaN');
assert(abs(axW.x(2) - 2*pi) < 1e-12, 'wavelength of lambda=1 is 2pi');
assert(abs(axW.x(3) - 1)    < 1e-12, 'wavelength of lambda=(2pi)^2 is 1');

% ---- GetWindowAverage: mean power over a sample window ----
T = [1 2 3; 0 0 6];                    % [K x nTime], real
avgFull = view_eigenmode_spectrum('GetWindowAverage', T, []);
assert(abs(avgFull(1) - mean([1 4 9]))  < 1e-12, 'row1 full mean wrong');
assert(abs(avgFull(2) - mean([0 0 36])) < 1e-12, 'row2 full mean wrong');
avgWin = view_eigenmode_spectrum('GetWindowAverage', T, [1 2]);
assert(abs(avgWin(1) - mean([1 4]))  < 1e-12, 'row1 window mean wrong');
assert(abs(avgWin(2) - mean([0 0]))  < 1e-12, 'row2 window mean wrong');

fprintf('ALL TESTS PASSED: test_view_eigenmode_spectrum_pure\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `evaluate_matlab_code`:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_view_eigenmode_spectrum_pure.m')
```
Expected: FAIL — `Undefined function 'view_eigenmode_spectrum'` (file does not exist yet).

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/gui/view_eigenmode_spectrum.m` with the license header (copy the header block from `toolbox/gui/view_eigenmodes.m:11-29`, author "Diellor Basha, 2026"), then:

```matlab
function varargout = view_eigenmode_spectrum(varargin)
% VIEW_EIGENMODE_SPECTRUM: Modal power spectrum of a source map's activations.
%
% USAGE:  hFig = view_eigenmode_spectrum(ResultsFile)
%         pw   = view_eigenmode_spectrum('ComputeModalPower', ThetaCol, Component)
%         ax   = view_eigenmode_spectrum('GetSpectrumAxis', Values, mode)
%         avg  = view_eigenmode_spectrum('GetWindowAverage', Theta, iWin)
%         Th   = view_eigenmode_spectrum('CollapseProject', Eig, ImageGridAmp, nComp, M)
%
% Projects the realized vertex source map onto the surface LBO eigenmodes and
% displays power per mode (Left/Right hemisphere curves) vs eigenvalue (or
% spatial wavelength). The figure is registered in the source map's dataset and
% is driven by Brainstorm's global time cursor (see bst_figures FireCurrentTimeChanged).

methodNames = {'ComputeModalPower', 'GetSpectrumAxis', 'GetWindowAverage', ...
               'CollapseProject', 'CreateFigure', 'UpdateFigurePlot', ...
               'CurrentTimeChangedCallback'};
if (nargin >= 1) && ischar(varargin{1}) && ismember(varargin{1}, methodNames)
    [varargout{1:nargout}] = feval(varargin{:});
    return;
end
[varargout{1:nargout}] = ViewFigure(varargin{:});
end


%% ===== PURE: |theta|^2 split by hemisphere component =====
function pw = ComputeModalPower(ThetaCol, Component)
    p = abs(ThetaCol(:)) .^ 2;
    Component = Component(:);
    pw.left  = p(Component == 1);
    pw.right = p(Component == 2);
end


%% ===== PURE: spectrum x-axis (eigenvalue or spatial wavelength) =====
function ax = GetSpectrumAxis(Values, mode)
    Values = Values(:);
    switch lower(mode)
        case 'eigenvalue'
            ax.x     = Values;
            ax.label = 'Eigenvalue \lambda';
        case 'wavelength'
            x = nan(size(Values));
            pos = (Values > 0);
            x(pos) = 2 * pi ./ sqrt(Values(pos));
            ax.x     = x;
            ax.label = 'Spatial wavelength \approx 2\pi/\surd\lambda';
        otherwise
            error('Unknown spectrum axis mode: %s', mode);
    end
end


%% ===== PURE: mean modal power over a sample window =====
function avg = GetWindowAverage(Theta, iWin)
    if isempty(iWin)
        iWin = 1:size(Theta, 2);
    end
    avg = mean(abs(Theta(:, iWin)) .^ 2, 2);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run via MATLAB MCP `evaluate_matlab_code`:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_view_eigenmode_spectrum_pure.m')
```
Expected: PASS — prints `ALL TESTS PASSED: test_view_eigenmode_spectrum_pure`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmode_spectrum.m dev/tests/test_view_eigenmode_spectrum_pure.m
git commit -m "Eigenspectrum: pure render helpers (modal power, axis, window avg)"
```

---

### Task 2: Acquisition core — `CollapseProject`

**Files:**
- Modify: `toolbox/gui/view_eigenmode_spectrum.m` (add `CollapseProject` subfunction)
- Test: `dev/tests/test_eigenmode_spectrum_acquire.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmode_spectrum_acquire.m`:

```matlab
function test_eigenmode_spectrum_acquire
% Verify CollapseProject: constrained passthrough + unconstrained rms-collapse,
% both projected onto an M-orthonormal eigenbasis.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'gui'));

rng(0);
nV = 40; k = 6;
% Diagonal mass matrix and an M-orthonormal basis Phi (Phi' M Phi = I).
d = 0.5 + rand(nV, 1);
M = spdiags(d, 0, nV, nV);
A = randn(nV, k);
G = A' * (M * A);
R = chol((G + G') / 2);
Phi = A / R;
assert(norm(Phi' * (M * Phi) - eye(k), 'fro') < 1e-9, 'setup: Phi not M-orthonormal');
Em = struct('Vectors', Phi, 'Values', (1:k)', 'nModes', k, 'MassType', 'barycentric', ...
            'Component', [ones(3,1); 2*ones(3,1)], 'CompRank', [1;2;3;1;2;3], 'nComponents', 2);

% ---- Constrained (nComponents=1): a pure eigenmode field -> unit coeff at that mode ----
S = Phi(:, 3) * [1, 0, 5];                    % [nV x 3], column 1 is mode 3, scaled across time
Theta = view_eigenmode_spectrum('CollapseProject', Em, S, 1, M);
assert(isequal(size(Theta), [k 3]), 'constrained Theta wrong shape');
assert(abs(Theta(3,1) - 1) < 1e-9, 'pure mode-3 field must give unit coeff at mode 3');
assert(max(abs(Theta([1 2 4 5 6], 1))) < 1e-9, 'other modes must be ~0');
assert(max(abs(Theta(:,2))) < 1e-9, 'zero field column must give zero coeffs');
% Constrained path must equal a direct projection.
ThetaDirect = bst_eigenmodes_project(Em, S, M);
assert(max(abs(Theta(:) - ThetaDirect(:))) < 1e-9, 'constrained path must equal direct projection');

% ---- Unconstrained (nComponents=3): rms-collapse then project ----
A3 = randn(3*nV, 4);                           % [3nV x nTime]
expectedS = bst_source_orient([], 3, [], A3, 'rms');   % [nV x nTime]
expectedTheta = bst_eigenmodes_project(Em, expectedS, M);
ThetaU = view_eigenmode_spectrum('CollapseProject', Em, A3, 3, M);
assert(isequal(size(ThetaU), [k 4]), 'unconstrained Theta wrong shape');
assert(max(abs(ThetaU(:) - expectedTheta(:))) < 1e-9, 'unconstrained path must equal collapse-then-project');

fprintf('ALL TESTS PASSED: test_eigenmode_spectrum_acquire\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `evaluate_matlab_code`:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_spectrum_acquire.m')
```
Expected: FAIL — `Unknown command or method 'CollapseProject'` from the dispatcher (subfunction not defined yet).

- [ ] **Step 3: Write minimal implementation**

Add this subfunction to `toolbox/gui/view_eigenmode_spectrum.m` (after `GetWindowAverage`):

```matlab
%% ===== CORE: realized vertex field -> eigenmode coefficients =====
% ImageGridAmp: [nComp*nVert x nTime] source matrix (kernel already applied).
% nComp: 1 (constrained, signed) or 2/3 (unconstrained -> rms magnitude per vertex).
% M: [nVert x nVert] sparse mass matrix for the eigenbasis Eig.Vectors.
function Theta = CollapseProject(Eig, ImageGridAmp, nComp, M)
    if (nComp == 3) || (nComp == 2)
        S = bst_source_orient([], nComp, [], ImageGridAmp, 'rms');   % [nVert x nTime]
    else
        S = ImageGridAmp;                                            % [nVert x nTime]
    end
    Theta = bst_eigenmodes_project(Eig, S, M);                       % [nModes x nTime]
end
```

- [ ] **Step 4: Run test to verify it passes**

Run via MATLAB MCP `evaluate_matlab_code`:
```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_spectrum_acquire.m')
```
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmode_spectrum_acquire`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/view_eigenmode_spectrum.m dev/tests/test_eigenmode_spectrum_acquire.m
git commit -m "Eigenspectrum: CollapseProject core (constrained + unconstrained rms)"
```

---

### Task 3: IO wrapper + GUI (figure, plot, time callback, keys)

**Files:**
- Modify: `toolbox/gui/view_eigenmode_spectrum.m` (add `GetActivationCoeffs`, `ViewFigure`, `CreateFigure`, `UpdateFigurePlot`, `CurrentTimeChangedCallback`, `FigureKeyPressedCallback`)

This task is GUI/DB-bound and is validated by static analysis (`check_matlab_code`) plus the interactive check at the end of the plan; there is no headless unit test (matches the spec's testing strategy). Provide complete implementations.

- [ ] **Step 1: Add the IO wrapper**

Add to `toolbox/gui/view_eigenmode_spectrum.m`:

```matlab
%% ===== IO: load realized source map + eigenmodes -> coefficients =====
% Returns Theta [nModes x nTime] and Info (Values/Component/CompRank/Time/SurfaceFile).
% Returns Theta=[] (after a bst_error) when the surface has no eigenmodes.
function [Theta, Info] = GetActivationCoeffs(ResultsFile)
    Theta = [];
    Info  = [];
    % Realized vertex source map (kernel applied on the fly).
    ResultsMat = in_bst_results(ResultsFile, 1);
    if ~isfield(ResultsMat, 'SurfaceFile') || isempty(ResultsMat.SurfaceFile)
        bst_error('This source file has no associated surface.', 'Eigenspectrum', 0);
        return;
    end
    % Eigenmodes on that surface.
    [Eig, isComputed] = in_tess_eigenmodes(ResultsMat.SurfaceFile);
    if ~isComputed || isempty(Eig) || ~isfield(Eig, 'Vectors') || isempty(Eig.Vectors)
        bst_error(['No eigenmodes found on this surface.' 10 ...
                   'Right-click the cortex and run "Compute eigenmodes" first.'], 'Eigenspectrum', 0);
        return;
    end
    % Mass matrix consistent with the stored mass type.
    sSurf = in_tess_bst(ResultsMat.SurfaceFile);
    [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
    % Project.
    nComp = ResultsMat.nComponents;
    Theta = CollapseProject(Eig, ResultsMat.ImageGridAmp, nComp, M);
    Info = struct('Values',      Eig.Values(:), ...
                  'Component',   Eig.Component(:), ...
                  'CompRank',    Eig.CompRank(:), ...
                  'Time',        ResultsMat.Time, ...
                  'SurfaceFile', ResultsMat.SurfaceFile);
end
```

- [ ] **Step 2: Add the main GUI entry `ViewFigure`**

```matlab
%% ===== GUI: open the spectrum figure registered in the source map's dataset =====
function hFig = ViewFigure(ResultsFile)
    global GlobalData;
    hFig = [];
    bst_progress('start', 'Eigenspectrum', 'Loading source activations...');
    % Coefficients first (fails fast with a friendly error if no eigenmodes).
    [Theta, Info] = GetActivationCoeffs(ResultsFile);
    if isempty(Theta)
        bst_progress('stop');
        return;
    end
    % Load the results into a dataset so the global time cursor spans this file.
    [iDS, ~] = bst_memory('LoadResultsFile', ResultsFile);
    if isempty(iDS)
        bst_progress('stop');
        return;
    end
    % Create a managed figure of our type, always a fresh one.
    FigureId = db_template('FigureId');
    FigureId.Type     = 'EigenSpectrum';
    FigureId.SubType  = '';
    FigureId.Modality = '';
    [hFig, ~] = bst_figures('CreateFigure', iDS, FigureId, 'AlwaysCreate', ResultsFile);
    if isempty(hFig)
        bst_progress('stop');
        return;
    end
    % Stash everything the redraw needs.
    setappdata(hFig, 'Theta',       Theta);
    setappdata(hFig, 'SpecInfo',    Info);
    setappdata(hFig, 'ResultsFile', ResultsFile);
    setappdata(hFig, 'AxisMode',    'eigenvalue');
    setappdata(hFig, 'isStatic',    size(Theta, 2) <= 2);
    % First plot + show + select.
    UpdateFigurePlot(hFig);
    set(hFig, 'Visible', 'on');
    bst_figures('SetCurrentFigure', hFig, '2D');
    bst_progress('stop');
end
```

- [ ] **Step 3: Add `CreateFigure` (figure shell + callbacks)**

```matlab
%% ===== GUI: build the bare figure (called by bst_figures CreateFigure) =====
function hFig = CreateFigure(FigureId) %#ok<INUSD>
    hFig = figure( ...
        'Visible',       'off', ...
        'NumberTitle',   'off', ...
        'IntegerHandle', 'off', ...
        'MenuBar',       'none', ...
        'Toolbar',       'figure', ...
        'DockControls',  'on', ...
        'Units',         'pixels', ...
        'Color',         [.9 .9 .9], ...
        'Pointer',       'arrow', ...
        'Tag',           'EigenSpectrum', ...
        'Name',          'Eigenspectrum', ...
        'CloseRequestFcn', @(h,ev)bst_figures('DeleteFigure', h, ev), ...
        'KeyPressFcn',     @FigureKeyPressedCallback);
    setappdata(hFig, 'FigureId', FigureId);
    setappdata(hFig, 'isStatic', 0);
    % Single axes for the spectrum curves.
    axes('Parent', hFig, 'Tag', 'AxesEigenSpectrum', 'Units', 'normalized', ...
         'Position', [0.12 0.13 0.82 0.78]);
end
```

- [ ] **Step 4: Add `UpdateFigurePlot` (the redraw)**

```matlab
%% ===== GUI: redraw the spectrum at the current global time =====
function UpdateFigurePlot(hFig)
    global GlobalData;
    Theta    = getappdata(hFig, 'Theta');
    Info     = getappdata(hFig, 'SpecInfo');
    AxisMode = getappdata(hFig, 'AxisMode');
    if isempty(Theta)
        return;
    end
    % Current time -> column index.
    Time = Info.Time;
    if isempty(GlobalData) || isempty(GlobalData.UserTimeWindow.CurrentTime) || (numel(Time) < 2)
        iTime = 1;
    else
        iTime = bst_closest(GlobalData.UserTimeWindow.CurrentTime, Time);
    end
    iTime = max(1, min(iTime, size(Theta, 2)));
    % Power split + axis + full-window average overlay.
    pw  = ComputeModalPower(Theta(:, iTime), Info.Component);
    ax  = GetSpectrumAxis(Info.Values, AxisMode);
    avg = GetWindowAverage(Theta, []);
    Comp = Info.Component(:);
    xL = ax.x(Comp == 1);  xR = ax.x(Comp == 2);
    aL = avg(Comp == 1);   aR = avg(Comp == 2);
    % Draw.
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'AxesEigenSpectrum');
    cla(hAxes);
    hold(hAxes, 'on');
    plot(hAxes, xL, pw.left,  '-',  'Color', [0.85 0.2 0.2], 'LineWidth', 1.5);
    plot(hAxes, xR, pw.right, '-',  'Color', [0.2 0.3 0.85], 'LineWidth', 1.5);
    plot(hAxes, xL, aL, '--', 'Color', [0.85 0.2 0.2], 'LineWidth', 0.75);
    plot(hAxes, xR, aR, '--', 'Color', [0.2 0.3 0.85], 'LineWidth', 0.75);
    hold(hAxes, 'off');
    xlabel(hAxes, ax.label);
    ylabel(hAxes, 'Modal power |\theta_k|^2');
    legend(hAxes, {'Left (t)', 'Right (t)', 'Left (avg)', 'Right (avg)'}, 'Location', 'northeast');
    if (numel(Time) >= 1)
        tStr = sprintf('%1.3f s', Time(iTime));
    else
        tStr = 'n/a';
    end
    title(hAxes, sprintf('Eigenspectrum  |  t = %s  |  L=%d R=%d modes  |  [e/w axis, arrows step time]', ...
          tStr, numel(xL), numel(xR)), 'Interpreter', 'tex');
end
```

- [ ] **Step 5: Add the time callback and key handler**

```matlab
%% ===== GUI: global time cursor moved (called by bst_figures FireCurrentTimeChanged) =====
function CurrentTimeChangedCallback(hFig)
    if getappdata(hFig, 'isStatic')
        return;
    end
    UpdateFigurePlot(hFig);
end


%% ===== GUI: keys — 'e'/'w' toggle axis; arrows step the global time cursor =====
function FigureKeyPressedCallback(hFig, ev)
    global GlobalData;
    switch (ev.Key)
        case 'e'
            setappdata(hFig, 'AxisMode', 'eigenvalue');
            UpdateFigurePlot(hFig);
        case 'w'
            setappdata(hFig, 'AxisMode', 'wavelength');
            UpdateFigurePlot(hFig);
        case {'leftarrow', 'rightarrow', 'uparrow', 'downarrow', 'pageup', 'pagedown'}
            if isempty(GlobalData) || isempty(GlobalData.UserTimeWindow.SamplingRate)
                return;
            end
            sr = GlobalData.UserTimeWindow.SamplingRate;
            t  = GlobalData.UserTimeWindow.CurrentTime;
            switch (ev.Key)
                case {'leftarrow', 'downarrow'},  t = t - sr;
                case {'rightarrow', 'uparrow'},   t = t + sr;
                case 'pageup',                    t = t + 10 * sr;
                case 'pagedown',                  t = t - 10 * sr;
            end
            panel_time('SetCurrentTime', t);   % moves the one global clock -> redraws all figures
    end
end
```

- [ ] **Step 6: Static-analysis check**

Run via MATLAB MCP `check_matlab_code` on `toolbox/gui/view_eigenmode_spectrum.m`.
Expected: no errors (undefined variables/functions, syntax). Brainstorm-idiom warnings (e.g. unused `varargin`) are acceptable; there must be no **error**-level findings. Fix any error-level finding before committing.

- [ ] **Step 7: Re-run the two automated tests (no regression in helpers)**

```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_view_eigenmode_spectrum_pure.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_spectrum_acquire.m')
```
Expected: both print `ALL TESTS PASSED`.

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/view_eigenmode_spectrum.m
git commit -m "Eigenspectrum: IO wrapper + 2D viewer (figure, plot, time callback, keys)"
```

---

### Task 4: Register the figure type in `bst_figures.m`

**Files:**
- Modify: `toolbox/core/bst_figures.m` (CreateFigure switch `~:196`; FireCurrentTimeChanged switch `~:969`)

- [ ] **Step 1: Add the CreateFigure case**

In `CreateFigure`, the `switch(FigureId.Type)` block ends with:
```matlab
            case 'Video'
                hFig = figure_video('CreateFigure', FigureId);
                FigHandles = db_template('DisplayHandlesVideo');
            otherwise
                error(['Invalid figure type : ', FigureId.Type]);
```
Insert a new case immediately **before** `otherwise`:
```matlab
            case 'EigenSpectrum'
                hFig = view_eigenmode_spectrum('CreateFigure', FigureId);
                FigHandles = db_template('DisplayHandlesTimeSeries');
```

- [ ] **Step 2: Add the FireCurrentTimeChanged case**

In `FireCurrentTimeChanged`, the dispatch `switch (sFig.Id.Type)` contains:
```matlab
                case 'Spectrum'
                    figure_spectrum('CurrentTimeChangedCallback', sFig.hFigure);
```
Add immediately **after** that `case 'Spectrum'` block:
```matlab
                case 'EigenSpectrum'
                    view_eigenmode_spectrum('CurrentTimeChangedCallback', sFig.hFigure);
```

- [ ] **Step 3: Static-analysis check**

Run via MATLAB MCP `check_matlab_code` on `toolbox/core/bst_figures.m`.
Expected: no **new** error-level findings introduced by the two added cases (this is a large pre-existing file; only confirm the additions parse and reference defined names).

- [ ] **Step 4: Commit**

```bash
git add toolbox/core/bst_figures.m
git commit -m "Eigenspectrum: register EigenSpectrum figure type (create + time-sync)"
```

---

### Task 5: Context-menu hook in `tree_callbacks.m`

**Files:**
- Modify: `toolbox/tree/tree_callbacks.m` (results/link activations submenu `~:1825`)

- [ ] **Step 1: Add the menu item**

In the `case {'results', 'link'}` block, inside the `if (length(bstNodes) == 1)` body, the "Display on cortex" item is:
```matlab
                    if ismember(sStudy.Result(iResult).HeadModelType, {'surface', 'mixed'})
                        if ~isempty(sSubject) && ~isempty(sSubject.iCortex)
                            gui_component('MenuItem', jMenuActivations, [], 'Display on cortex', IconLoader.ICON_CORTEX, [], @(h,ev)view_surface_data([], filenameRelative));
                        else
                            gui_component('MenuItem', jMenuActivations, [], 'No cortex available', IconLoader.ICON_WARNING, [], []);
                        end
                    end
```
Add the eigenspectrum item immediately **after** that `end` (still inside `if (length(bstNodes) == 1)`), gated on the same surface/cortex availability:
```matlab
                    % === VIEW EIGENSPECTRUM ===
                    if ismember(sStudy.Result(iResult).HeadModelType, {'surface', 'mixed'}) && ~isempty(sSubject) && ~isempty(sSubject.iCortex)
                        gui_component('MenuItem', jMenuActivations, [], 'View eigenspectrum', IconLoader.ICON_TIMEFREQ, [], @(h,ev)view_eigenmode_spectrum(filenameRelative));
                    end
```

- [ ] **Step 2: Static-analysis check**

Run via MATLAB MCP `check_matlab_code` on `toolbox/tree/tree_callbacks.m`.
Expected: no **new** error-level findings from the added block (pre-existing file; confirm the addition parses and `filenameRelative`/`jMenuActivations`/`IconLoader` are in scope — they are, per `:1818` and `:1825`).

- [ ] **Step 3: Commit**

```bash
git add toolbox/tree/tree_callbacks.m
git commit -m "Eigenspectrum: 'View eigenspectrum' menu on source-results nodes"
```

---

### Task 6: Full regression + interactive validation handoff

**Files:** none (verification only)

- [ ] **Step 1: Run the new automated tests**

```matlab
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_view_eigenmode_spectrum_pure.m')
run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/test_eigenmode_spectrum_acquire.m')
```
Expected: both print `ALL TESTS PASSED`.

- [ ] **Step 2: Run the existing eigenmode regression suite**

Run each via MATLAB MCP `evaluate_matlab_code` (these must remain green — this increment changes none of their paths):
```matlab
run('.../dev/tests/test_eigenmodes_project_pure.m')
run('.../dev/tests/test_eigenmodes_transform_pure.m')
run('.../dev/tests/test_io_eigenmodes_roundtrip.m')
run('.../dev/tests/test_eigenmodes_perhemisphere.m')
run('.../dev/tests/test_eigenmodes_manifold_gate.m')
```
(Use the absolute repo path `/Users/diellorbasha/workspace/research/code/brainstorm3` in place of `...`.)
Expected: each prints its `ALL TESTS PASSED` line. If any fails, STOP and report — do not proceed.

- [ ] **Step 3: Report for interactive validation**

Report to the user that the feature is implemented and all automated tests pass, and list the interactive checks to perform in the Brainstorm GUI (per the spec's testing strategy item 4):
1. Right-click a **source result** (with computed surface eigenmodes) → **Cortical activations → View eigenspectrum** → the spectrum figure opens.
2. Step the **global time cursor** (slider / arrows) → the L/R curves update live, in lockstep with a cortex display of the same file if open.
3. The **window-averaged** dashed overlay is visible.
4. Press **e** / **w** → x-axis toggles between eigenvalue and spatial wavelength.
5. Right-click a source result on a surface **without** eigenmodes → friendly error pointing to "Compute eigenmodes".
6. Repeat (1)-(2) with both a **constrained** and an **unconstrained** source model.

Do not merge or open a PR; the user validates interactively first (per the established workflow).

---

## Notes & deferrals (carried from the spec)

- **Window-average follows full range** in this increment (the dashed overlay = mean over all loaded samples). Following an active time-panel sub-selection is deferred to interactive feedback (the spec allows the full-range default).
- **No persistence / no results-structure change.** The viewer is transient.
- **Out of scope:** manifold-harmonics source mapping (the `[nEigenmodes × nChannels]` kernel — separate, already implemented), other vertex-mapped scalars (PET etc.), and the 2D mode×frequency map.
