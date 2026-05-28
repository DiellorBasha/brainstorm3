# Eigenmode Coefficient Filter Library Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Apply per-mode spatial-spectral gains `h(λₖ)` to eigenmode coefficient time series (the "then filter" half of transform-first), via a shared transfer-function source reused by both the vertex filter and a new coefficient-filter process.

**Architecture:** Extract the transfer-function math into a pure `bst_eigenmodes_filter_gain` (single source, + new `tikhonov`); refactor the existing vertex-field `bst_eigenmodes_filter` to delegate to it (behavior-preserving, guarded by its existing test); add `process_eigenmodes_coeffsfilter` that applies the gain to an eigenmode-coefficient `matrix` file.

**Tech Stack:** MATLAB, Brainstorm process system (`eval(macro_method)`), `in_bst_matrix`, `in_tess_eigenmodes`, `db_template`, `bst_save`/`db_add_data`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_eigenmodes_filter_gain.m` (create) | Pure: `(lambdas, FilterType, …) → h [K×1]` gain. Single source of transfer functions. |
| `toolbox/math/bst_eigenmodes_filter.m` (modify) | Refactor: delegate gain to `bst_eigenmodes_filter_gain`; keep project→multiply→reconstruct. |
| `toolbox/process/functions/process_eigenmodes_coeffsfilter.m` (create) | Process: matrix coefficients in → `h.*θ` filtered coefficients out (+ optional vertex reconstruction). |
| `dev/tests/test_eigenmodes_filter_gain_pure.m` (create) | DB-free unit test of the gain function (each type). |
| `dev/tests/test_process_eigenmodes_coeffsfilter_options.m` (create) | DB-free options test of the process. |

**Run convention (tests):** MATLAB MCP `evaluate_matlab_code` calling the function name (script-style tests print `ALL TESTS PASSED`; NOT `run_matlab_test_file`). Lint via `check_matlab_code`.

Confirmed facts: the current `bst_eigenmodes_filter` transfer block + option defaults are `CutoffMode=50`, `ModeRange=[20 80]`, `DiffusionTime=0.01`, `MaxGain=10`, `TransferFn=[]`. Its existing guard test is `dev/tests/test_eigenmodes_filter_pure.m`. `in_bst_matrix(file)` returns `.Value [K×nTime]`, `.Time`, `.Description`, `.SurfaceFile`. The vertex filter process GUI uses options `filtertype` (radio_linelabel), `cutoffmode`, `moderange_low`, `moderange_high`, `diffusiontime`.

---

## Task 1: Shared gain function + refactor the vertex filter

**Files:**
- Create: `toolbox/math/bst_eigenmodes_filter_gain.m`
- Modify: `toolbox/math/bst_eigenmodes_filter.m`
- Test: `dev/tests/test_eigenmodes_filter_gain_pure.m` (new); `dev/tests/test_eigenmodes_filter_pure.m` (existing, must still pass)

- [ ] **Step 1: Write the failing test for the gain function**

Create `dev/tests/test_eigenmodes_filter_gain_pure.m`:

```matlab
function test_eigenmodes_filter_gain_pure
% Verify each per-mode transfer-function gain.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

lam = ((0:9)').^2;          % lambda_0=0, strictly increasing
K = numel(lam);

% lowpass: keep modes 1..4
hLP = bst_eigenmodes_filter_gain(lam, 'lowpass', 'CutoffMode', 4);
assert(isequal(hLP, [ones(4,1); zeros(K-4,1)]), 'lowpass mask wrong.');
% highpass: keep modes 4..end
hHP = bst_eigenmodes_filter_gain(lam, 'highpass', 'CutoffMode', 4);
assert(isequal(hHP, [zeros(3,1); ones(K-3,1)]), 'highpass mask wrong.');
% bandpass: keep modes 3..6
hBP = bst_eigenmodes_filter_gain(lam, 'bandpass', 'ModeRange', [3 6]);
expBP = zeros(K,1); expBP(3:6) = 1;
assert(isequal(hBP, expBP), 'bandpass mask wrong.');
% heat: t->0 is ~identity; large t suppresses high lambda monotonically
hH0 = bst_eigenmodes_filter_gain(lam, 'heat', 'DiffusionTime', 1e-12);
assert(max(abs(hH0 - 1)) < 1e-6, 'heat t->0 should be ~1.');
hH1 = bst_eigenmodes_filter_gain(lam, 'heat', 'DiffusionTime', 1);
assert(abs(hH1(1)-1) < 1e-12 && hH1(end) < 1e-6 && all(diff(hH1) <= 0), 'heat should decay with lambda.');
% inverse_heat: clamped at MaxGain
hIH = bst_eigenmodes_filter_gain(lam, 'inverse_heat', 'DiffusionTime', 1, 'MaxGain', 5);
assert(all(hIH <= 5 + 1e-9) && any(abs(hIH - 5) < 1e-9), 'inverse_heat should clamp at MaxGain.');
% tikhonov: 1 at lambda=0, decreasing, in (0,1]
hT = bst_eigenmodes_filter_gain(lam, 'tikhonov', 'RegBeta', 1);
assert(abs(hT(1)-1) < 1e-12 && all(hT > 0 & hT <= 1) && all(diff(hT) <= 0), 'tikhonov shape wrong.');
% custom: passthrough of the supplied handle
hC = bst_eigenmodes_filter_gain(lam, 'custom', 'TransferFn', @(l) 2*ones(size(l)));
assert(all(hC == 2), 'custom passthrough wrong.');
% unknown type errors
threw = false;
try, bst_eigenmodes_filter_gain(lam, 'bogus'); catch, threw = true; end
assert(threw, 'unknown filter type should error.');

fprintf('ALL TESTS PASSED: test_eigenmodes_filter_gain_pure\n');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `test_eigenmodes_filter_gain_pure`
Expected: FAIL — `Undefined function or variable 'bst_eigenmodes_filter_gain'`.

- [ ] **Step 3: Create the gain function**

Create `toolbox/math/bst_eigenmodes_filter_gain.m`:

```matlab
function h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin)
% BST_EIGENMODES_FILTER_GAIN: Per-mode transfer-function gain for eigenmode filtering.
%
% USAGE:  h = bst_eigenmodes_filter_gain(lambdas, 'lowpass',  'CutoffMode', 50)
%         h = bst_eigenmodes_filter_gain(lambdas, 'bandpass', 'ModeRange', [20 80])
%         h = bst_eigenmodes_filter_gain(lambdas, 'heat',     'DiffusionTime', 0.01)
%         h = bst_eigenmodes_filter_gain(lambdas, 'tikhonov', 'RegBeta', 1)
%
% DESCRIPTION:
%     Returns the per-mode gain vector h(lambda_k) for a spatial-spectral filter
%     in the Laplace-Beltrami eigenmode basis. This is the single source of the
%     transfer functions shared by bst_eigenmodes_filter (vertex fields) and
%     process_eigenmodes_coeffsfilter (eigenmode coefficients).
%
%     Types: 'lowpass'/'highpass'/'bandpass' (mode-index cutoffs), 'heat'
%     (exp(-t*lambda)), 'inverse_heat' (exp(+t*lambda) clamped at MaxGain),
%     'tikhonov' (1/(1+beta*lambda)), 'custom' (user TransferFn handle).
%
% INPUTS:
%     lambdas    : [K x 1] eigenvalues (K = number of modes).
%     FilterType : one of the types above.
%
% OPTIONS (name-value):
%     'CutoffMode'    (50)      mode index for low/high-pass
%     'ModeRange'     ([20 80]) [k1 k2] for band-pass
%     'DiffusionTime' (0.01)    t for heat / inverse_heat
%     'MaxGain'       (10)      clamp for inverse_heat
%     'RegBeta'       (1)       beta for tikhonov
%     'TransferFn'    ([])      function handle for custom
%
% OUTPUTS:
%     h : [K x 1] gain vector.
%
% SEE ALSO: bst_eigenmodes_filter, process_eigenmodes_coeffsfilter

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

%% ===== PARSE OPTIONS =====
CutoffMode    = 50;
ModeRange     = [20, 80];
DiffusionTime = 0.01;
MaxGain       = 10;
RegBeta       = 1;
TransferFn    = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'cutoffmode',    CutoffMode    = varargin{i+1};
        case 'moderange',     ModeRange     = varargin{i+1};
        case 'diffusiontime', DiffusionTime = varargin{i+1};
        case 'maxgain',       MaxGain       = varargin{i+1};
        case 'regbeta',       RegBeta       = varargin{i+1};
        case 'transferfn',    TransferFn    = varargin{i+1};
    end
end

lambdas = double(lambdas(:));
nModes  = numel(lambdas);
h = zeros(nModes, 1);

%% ===== TRANSFER FUNCTION =====
switch lower(FilterType)
    case 'lowpass'
        c = min(CutoffMode, nModes);
        h(1:c) = 1;
    case 'highpass'
        c = max(1, min(CutoffMode, nModes));
        h(c:end) = 1;
    case 'bandpass'
        k1 = max(1, ModeRange(1));
        k2 = min(nModes, ModeRange(2));
        h(k1:k2) = 1;
    case 'heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = exp(-DiffusionTime * lambdas);
    case 'inverse_heat'
        if DiffusionTime <= 0
            error('DiffusionTime must be positive (got %g).', DiffusionTime);
        end
        h = min(exp(DiffusionTime * lambdas), MaxGain);
    case 'tikhonov'
        if RegBeta < 0
            error('RegBeta must be non-negative (got %g).', RegBeta);
        end
        h = 1 ./ (1 + RegBeta * lambdas);
    case 'custom'
        if isempty(TransferFn) || ~isa(TransferFn, 'function_handle')
            error('Custom filter requires a TransferFn option (function handle).');
        end
        h = TransferFn(lambdas);
        if numel(h) ~= nModes
            error('TransferFn must return a vector of length %d (got %d).', nModes, numel(h));
        end
        h = h(:);
    otherwise
        error('Unknown filter type: %s. Use lowpass, highpass, bandpass, heat, inverse_heat, tikhonov, or custom.', FilterType);
end
end
```

- [ ] **Step 4: Run the gain test to verify it passes**

Run: `test_eigenmodes_filter_gain_pure`
Expected: PASS — `ALL TESTS PASSED: test_eigenmodes_filter_gain_pure`.

- [ ] **Step 5: Refactor `bst_eigenmodes_filter` to delegate to the gain function**

In `toolbox/math/bst_eigenmodes_filter.m`, replace the entire body from the line `%% ===== PARSE INPUTS =====` through the final `end` (i.e. the option-parse block, the `Phi`/`lambdas`/`nV`/`nModes` extraction, the VALIDATE block, the BUILD TRANSFER FUNCTION switch, and the APPLY block) with exactly:

```matlab
%% ===== EXTRACT EIGENMODES =====
Phi     = double(Eigenmodes.Vectors);    % [nV x nModes]
lambdas = double(Eigenmodes.Values(:));   % [nModes x 1]
nV      = size(Phi, 1);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== TRANSFER FUNCTION (shared single source) =====
h = bst_eigenmodes_filter_gain(lambdas, FilterType, varargin{:});

%% ===== APPLY: u_filtered = Phi * (h .* (Phi' * M * u)) =====
Coeffs   = Phi' * (MassMatrix * Data);
Coeffs   = bsxfun(@times, h, Coeffs);
Filtered = Phi * Coeffs;
end
```

Also, in `bst_eigenmodes_filter.m`'s header, update the `FilterType :` line (and the SEE ALSO) to mention `tikhonov` and the gain function — change:
```
%     FilterType : 'lowpass','highpass','bandpass','heat','inverse_heat','custom'
```
to:
```
%     FilterType : 'lowpass','highpass','bandpass','heat','inverse_heat','tikhonov','custom'
%                  (see bst_eigenmodes_filter_gain for the transfer functions)
```

- [ ] **Step 6: Verify BOTH math tests pass (the refactor is behavior-preserving)**

Run: `test_eigenmodes_filter_gain_pure` then `test_eigenmodes_filter_pure`
Expected: both print `ALL TESTS PASSED: …`. (The existing `test_eigenmodes_filter_pure` guards that the refactor did not change vertex-filter behavior.)

- [ ] **Step 7: Lint**

`check_matlab_code` on `toolbox/math/bst_eigenmodes_filter_gain.m` and `toolbox/math/bst_eigenmodes_filter.m` — no genuine errors (no unused-variable warnings from leftover locals in the refactored file).

- [ ] **Step 8: Commit**

```bash
git add toolbox/math/bst_eigenmodes_filter_gain.m toolbox/math/bst_eigenmodes_filter.m dev/tests/test_eigenmodes_filter_gain_pure.m
git commit -m "$(cat <<'EOF'
Extract bst_eigenmodes_filter_gain; refactor bst_eigenmodes_filter to use it

Single source of transfer functions (lowpass/highpass/bandpass/heat/
inverse_heat/tikhonov/custom); vertex filter now delegates gain to it
(behavior-preserving, guarded by test_eigenmodes_filter_pure). Adds the
tikhonov type. Covered by test_eigenmodes_filter_gain_pure.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Coefficient-filter process `process_eigenmodes_coeffsfilter`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_coeffsfilter.m`
- Test: `dev/tests/test_process_eigenmodes_coeffsfilter_options.m`

- [ ] **Step 1: Write the failing options test**

Create `dev/tests/test_process_eigenmodes_coeffsfilter_options.m`:

```matlab
function test_process_eigenmodes_coeffsfilter_options
% Verify the coefficient-filter process is a matrix-in / matrix-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_coeffsfilter('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.9) < 1e-9, 'Index must be 336.9.');
assert(isequal(sProcess.InputTypes, {'matrix'}), 'InputTypes must be {matrix}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}), 'OutputTypes must be {matrix}.');
for f = {'filtertype','cutoffmode','moderange_low','moderange_high','diffusiontime','regbeta','dorecon'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
fprintf('ALL TESTS PASSED: test_process_eigenmodes_coeffsfilter_options\n');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `test_process_eigenmodes_coeffsfilter_options`
Expected: FAIL — `Undefined function or variable 'process_eigenmodes_coeffsfilter'`.

- [ ] **Step 3: Create the process**

Create `toolbox/process/functions/process_eigenmodes_coeffsfilter.m`:

```matlab
function varargout = process_eigenmodes_coeffsfilter( varargin )
% PROCESS_EIGENMODES_COEFFSFILTER: Spatial-spectral filter of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_coeffsfilter('GetDescription')
%       OutputFiles = process_eigenmodes_coeffsfilter('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Applies a per-mode spatial-spectral gain h(lambda_k) (see
%     bst_eigenmodes_filter_gain) directly to an eigenmode-coefficient matrix
%     (matrix_eigentransform, [K x nTime]): theta_filt = h .* theta. This is the
%     "filter" step of the transform-first workflow. Optionally reconstructs the
%     filtered vertex sources Q = Phi * theta_filt.
%
% SEE ALSO: bst_eigenmodes_filter_gain, process_eigenmodes_transform

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


%% ===== GET DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Eigenmode coefficient filter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.9;   % after the dispersion analysis (336.8)
    sProcess.Description = '';
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.filtertype.Comment = {'Low-pass', 'High-pass', 'Band-pass', 'Heat (smooth)', 'Inverse heat (sharpen)', 'Tikhonov', 'Filter:'; ...
                                           'lowpass', 'highpass', 'bandpass', 'heat', 'inverse_heat', 'tikhonov', ''};
    sProcess.options.filtertype.Type    = 'radio_linelabel';
    sProcess.options.filtertype.Value   = 'heat';

    sProcess.options.cutoffmode.Comment = 'Cutoff mode index (low/high-pass): ';
    sProcess.options.cutoffmode.Type    = 'value';
    sProcess.options.cutoffmode.Value   = {50, '', 0};

    sProcess.options.moderange_low.Comment  = 'Band-pass: lower mode index: ';
    sProcess.options.moderange_low.Type     = 'value';
    sProcess.options.moderange_low.Value    = {20, '', 0};

    sProcess.options.moderange_high.Comment = 'Band-pass: upper mode index: ';
    sProcess.options.moderange_high.Type    = 'value';
    sProcess.options.moderange_high.Value   = {80, '', 0};

    sProcess.options.diffusiontime.Comment  = 'Diffusion time (heat / inverse heat): ';
    sProcess.options.diffusiontime.Type     = 'value';
    sProcess.options.diffusiontime.Value    = {0.005, '', 4};

    sProcess.options.regbeta.Comment        = 'Tikhonov beta (h = 1/(1+beta*lambda)): ';
    sProcess.options.regbeta.Type           = 'value';
    sProcess.options.regbeta.Value          = {1, '', 4};

    sProcess.options.dorecon.Comment        = 'Also reconstruct filtered vertex sources (Q = Phi * theta_filt)';
    sProcess.options.dorecon.Type           = 'checkbox';
    sProcess.options.dorecon.Value          = 0;

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Applies a per-mode gain h(&lambda;) to the ' ...
        'eigenmode coefficients (the "filter" step of transform-first).<BR>Input: an eigenmode-coefficient matrix ' ...
        '(matrix_eigentransform); requires eigenmodes on the surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = ['Eigenmode coefficient filter (' sProcess.options.filtertype.Value ')'];
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    FilterType = lower(sProcess.options.filtertype.Value);
    CutoffMode = sProcess.options.cutoffmode.Value{1};
    ModeLow    = sProcess.options.moderange_low.Value{1};
    ModeHigh   = sProcess.options.moderange_high.Value{1};
    DiffTime   = sProcess.options.diffusiontime.Value{1};
    RegBeta    = sProcess.options.regbeta.Value{1};
    DoRecon    = sProcess.options.dorecon.Value;

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        M = in_bst_matrix(sInput.FileName);
        Coeffs = double(M.Value);                       % [K x nTime]
        K = size(Coeffs, 1);

        if ~isfield(M, 'SurfaceFile') || isempty(M.SurfaceFile)
            bst_report('Error', sProcess, sInput, 'Matrix has no SurfaceFile; cannot obtain eigenmodes.');
            continue;
        end
        [Em, isComputed] = in_tess_eigenmodes(M.SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on surface: ' M.SurfaceFile]);
            continue;
        end
        if size(Em.Values, 1) < K
            bst_report('Error', sProcess, sInput, sprintf('Coefficients have %d modes but surface has only %d eigenvalues.', K, size(Em.Values,1)));
            continue;
        end
        lambdas = double(Em.Values(1:K));

        h = bst_eigenmodes_filter_gain(lambdas, FilterType, ...
            'CutoffMode', CutoffMode, 'ModeRange', [ModeLow, ModeHigh], ...
            'DiffusionTime', DiffTime, 'RegBeta', RegBeta);
        CoeffsFilt = bsxfun(@times, h, Coeffs);         % [K x nTime]

        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        % --- Filtered coefficients (matrix) ---
        Mout = db_template('matrixmat');
        Mout.Value        = CoeffsFilt;
        Mout.Time         = M.Time;
        Mout.Description  = M.Description;
        Mout.SurfaceFile  = M.SurfaceFile;
        Mout.Comment      = sprintf('EigenFilt [%s] | %s', FilterType, sInput.Comment);
        Mout.nAvg         = 1;
        Mout = bst_history('add', Mout, 'eigenmodes_coeffsfilter', Mout.Comment);
        OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenfilt');
        bst_save(OutFile, Mout, 'v6');
        db_add_data(iStudyOut, OutFile, Mout);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>

        % --- Optional filtered vertex reconstruction (results) ---
        if DoRecon
            Phi = double(Em.Vectors(:, 1:K));
            Q   = Phi * CoeffsFilt;                     % [nVert x nTime]
            ResMat = db_template('resultsmat');
            ResMat.ImageGridAmp  = Q;
            ResMat.ImagingKernel = [];
            ResMat.nComponents   = 1;
            ResMat.Time          = M.Time;
            ResMat.SurfaceFile   = M.SurfaceFile;
            ResMat.HeadModelType = 'surface';
            ResMat.Comment       = sprintf('EigenFilt recon [%s] | %s', FilterType, sInput.Comment);
            ResMat.nAvg          = 1;
            ResMat = bst_history('add', ResMat, 'eigenmodes_coeffsfilter', ResMat.Comment);
            OutFile2 = bst_process('GetNewFilename', StudyDir, 'results_eigenfilt');
            bst_save(OutFile2, ResMat, 'v6');
            db_add_data(iStudyOut, OutFile2, ResMat);
            OutputFiles{end+1} = file_short(OutFile2); %#ok<AGROW>
        end

        bst_report('Info', sProcess, sInput, sprintf('Filtered %d eigenmode coefficients with %s.', K, FilterType));
    end
end
```

- [ ] **Step 4: Run the options test to verify it passes**

Run: `test_process_eigenmodes_coeffsfilter_options`
Expected: PASS — `ALL TESTS PASSED: test_process_eigenmodes_coeffsfilter_options`.

- [ ] **Step 5: Lint**

`check_matlab_code` on the process file — only standard Brainstorm idioms (`varargout`, stale `%#ok`) acceptable.

- [ ] **Step 6: Live end-to-end validation (self-contained synthetic substrate)**

The options test does not exercise `Run`. Build a throwaway protocol with a sphere surface that has eigenmodes + a synthetic coefficient matrix, run the heat filter, confirm high modes are suppressed relative to low, then delete the protocol. In MATLAB:

```matlab
gui_brainstorm('DeleteProtocol', 'FiltSmoke');
gui_brainstorm('CreateProtocol', 'FiltSmoke', 0, 0);
[~, iSubject] = db_add_subject('FS', [], 0, 0);
sSubject = bst_get('Subject', iSubject);
anatDir = bst_fileparts(file_fullpath(sSubject.FileName));
[V, F] = tess_sphere(642);
Em = tess_eigenmodes(V, F, 'nModes', 40, 'Verbose', 0);
SurfFile = fullfile(anatDir, 'tess_filtsphere.mat');
bst_save(SurfFile, struct('Vertices', V, 'Faces', F, 'Comment', 'filtsphere', 'Eigenmodes', Em), 'v7');
db_add_surface(iSubject, SurfFile, 'filtsphere');
SurfRel = file_short(SurfFile);
K = Em.nModes; nT = 200;
Coeffs = randn(K, nT);
Mmat = db_template('matrixmat');
Mmat.Value = Coeffs; Mmat.Time = (0:nT-1)/100; Mmat.SurfaceFile = SurfRel;
Mmat.Description = arrayfun(@(k) sprintf('Mode %d', k), (1:K)', 'uni', 0);
Mmat.Comment = 'synthetic coeffs'; Mmat.nAvg = 1;
sSubject = bst_get('Subject', iSubject); iStudyFS = sSubject.iStudy(1);
StudyDir = bst_fileparts(file_fullpath(bst_get('Study', iStudyFS).FileName));
fM = bst_process('GetNewFilename', StudyDir, 'matrix_eigentransform');
bst_save(fM, Mmat, 'v6'); db_add_data(iStudyFS, fM, Mmat);
sMat = bst_process('CallProcess', 'process_select_files_matrix', [], [], 'subjectname', 'FS');
iSel = find(arrayfun(@(s) ~isempty(strfind(s.FileName,'matrix_eigentransform')), sMat), 1);
sFilt = bst_process('CallProcess', 'process_eigenmodes_coeffsfilter', sMat(iSel), [], ...
    'filtertype', 'heat', 'diffusiontime', 0.05, 'dorecon', 0);
assert(~isempty(sFilt) && ~isempty(sFilt(1).FileName), 'no output');
R = in_bst_matrix(sFilt(1).FileName);
hiIn = sum(Coeffs(end-4:end,:).^2, 'all'); loIn = sum(Coeffs(1:5,:).^2, 'all');
hiOut = sum(R.Value(end-4:end,:).^2, 'all'); loOut = sum(R.Value(1:5,:).^2, 'all');
fprintf('hi/lo energy ratio: in=%.4f out=%.4f\n', hiIn/loIn, hiOut/loOut);
assert((hiOut/loOut) < (hiIn/loIn), 'heat filter should reduce high-mode energy relative to low.');
fprintf('LIVE VALIDATION OK\n');
gui_brainstorm('DeleteProtocol', 'FiltSmoke');
```

Expected: prints `LIVE VALIDATION OK` (heat filter suppresses high modes). If substrate-registration glue needs adjustment (db_add_surface, default-study lookup), DEBUG and fix the VALIDATION SCRIPT (not Run); if Run has a bug it surfaces here — fix Run, keep the contract. ALWAYS delete `FiltSmoke` at the end. If substrate setup proves intractable, report DONE_WITH_CONCERNS (the pure tests cover the gain math; the process is thin glue) — still clean up FiltSmoke.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_coeffsfilter.m dev/tests/test_process_eigenmodes_coeffsfilter_options.m
git commit -m "$(cat <<'EOF'
Add process_eigenmodes_coeffsfilter: spatial-spectral filter of coefficients

Thin process: matrix coefficients in -> h(lambda_k) gain (via
bst_eigenmodes_filter_gain) applied as theta_filt = h.*theta -> filtered
coefficients matrix (+ optional vertex reconstruction). Validated on a
synthetic coefficient matrix (heat suppresses high modes).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Shared `bst_eigenmodes_filter_gain` (single source, + tikhonov) → Task 1. ✓
- Refactor `bst_eigenmodes_filter` to delegate (behavior-preserving, guarded) → Task 1 Steps 5–6. ✓
- `process_eigenmodes_coeffsfilter` (matrix-in → `h.*θ` → filtered coefficients + optional recon), Index 336.9 → Task 2. ✓
- Filter types lowpass/highpass/bandpass/heat/inverse_heat/tikhonov/custom; option defaults matching the vertex filter → Task 1 gain function. ✓
- λ from SurfaceFile + eigenvalue-count guard → Task 2 Run. ✓
- Wiener-from-noise-floor / time-varying / vector → NOT built (deferred), as intended. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete code; run steps have exact commands + expected output.

**Type/name consistency:** `bst_eigenmodes_filter_gain(lambdas, FilterType, …) -> h [K×1]` used identically in its test, the `bst_eigenmodes_filter` refactor (`varargin{:}` passthrough), and the process (`'CutoffMode'/'ModeRange'/'DiffusionTime'/'RegBeta'`). Process option names (`moderange_low`/`moderange_high`) match the vertex-filter GUI; mapped to `ModeRange=[low high]` for the gain call. File prefixes `matrix_eigenfilt` / `results_eigenfilt`; Index 336.9 matches the options test.
