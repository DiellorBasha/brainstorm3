# Eigenmode Spatial Transform + Spectrum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add an unregularized sensor→eigenmode spatial transform `A = pinv(L·Φ)` that produces raw eigenmode coefficient time series, then FFT them through Brainstorm's existing tools to view the joint `(λ, ω)` spectrum.

**Architecture:** A pure SVD-pseudoinverse function (`bst_eigenmodes_transform`) does the linear algebra; a process plugin (`process_eigenmodes_transform`) resolves the head model / eigenmodes / channels, applies the kernel to recordings, and saves the coefficients as a Brainstorm matrix file. That matrix file feeds the existing `process_fft`; the standard spectrum viewer renders the `(λ, ω)` plane. No regularization, no new spectral code; the existing `bst_inverse_eigenmodes` is left untouched.

**Tech Stack:** MATLAB, Brainstorm process-plugin system (`eval(macro_method)`), Brainstorm io/db (`in_bst_headmodel`, `in_tess_eigenmodes`, `db_template`, `bst_save`, `db_add_data`), `process_fft`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_eigenmodes_transform.m` (create) | Pure: `Kernel = pinv(L̃)` via SVD, `+ Info` diagnostics. No I/O, no regularization. Pairs with `bst_eigenmodes_project`. |
| `toolbox/process/functions/process_eigenmodes_transform.m` (create) | Process: resolve head model/eigenmodes/channels, build kernel, apply to data → `matrix_eigentransform`; optional `results_eigentransform` vertex reconstruction. |
| `dev/tests/test_eigenmodes_transform_pure.m` (create) | DB-free unit test of the pure function (recovery + left/right-inverse invariants, both regimes). |
| `dev/tests/test_process_eigenmodes_transform_options.m` (create) | DB-free unit test that the process exposes the right options and *omits* method/prior/SNR. |

Brainstorm auto-discovers processes by scanning `toolbox/process/functions/`, so creating the file registers it; `SubGroup='Sources'` + `Index=338` place it just before the regularized inverse (339). No manual registration.

**Run convention (all tests):** execute the function in MATLAB with the repo on the path. Via the MATLAB MCP, use `run_matlab_test_file` on the test's absolute path, or `evaluate_matlab_code` with the function name. Each test ends by printing `ALL TESTS PASSED: <name>` and errors (via `assert`) on failure.

---

## Task 1: Pure transform function `bst_eigenmodes_transform`

**Files:**
- Create: `toolbox/math/bst_eigenmodes_transform.m`
- Test: `dev/tests/test_eigenmodes_transform_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_transform_pure.m`:

```matlab
function test_eigenmodes_transform_pure
% Verify the unregularized SVD pseudoinverse transform: exact recovery in the
% noise-free case, and the left/right-inverse invariants in both regimes.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(7);

% ----- Overdetermined regime: K < nch (well-determined least squares) -----
nCh = 40; nVert = 200; K = 25;
Gain    = randn(nCh, nVert);
Phi     = orth(randn(nVert, K));        % nVert x K, full column rank
L_tilde = Gain * Phi;

c0 = randn(K, 4);                        % known coefficients
D  = L_tilde * c0;                       % noise-free data
[Kernel, Info] = bst_eigenmodes_transform(Gain, Phi);

assert(isequal(size(Kernel), [K, nCh]), 'Kernel must be [K x nch].');
Theta = Kernel * D;
assert(max(abs(Theta(:) - c0(:))) < 1e-8, 'Did not recover coefficients (overdetermined).');
assert(norm(Kernel * L_tilde - eye(K), 'fro') < 1e-8, 'Left-inverse invariant Kernel*L=I_K failed.');
assert(Info.Rank == K, 'Rank should equal K when K<nch and full rank.');
assert(isfinite(Info.ConditionNumber) && Info.ConditionNumber >= 1, 'Bad condition number.');
assert(numel(Info.SingularValues) == K, 'SingularValues length mismatch.');

% ----- Underdetermined regime: K > nch (min-norm right inverse) -----
nCh2 = 20; nVert2 = 200; K2 = 35;
Gain2    = randn(nCh2, nVert2);
Phi2     = orth(randn(nVert2, K2));
L_tilde2 = Gain2 * Phi2;
[Kernel2, Info2] = bst_eigenmodes_transform(Gain2, Phi2);

assert(isequal(size(Kernel2), [K2, nCh2]), 'Kernel must be [K x nch] (underdetermined).');
assert(norm(L_tilde2 * Kernel2 - eye(nCh2), 'fro') < 1e-8, 'Right-inverse invariant L*Kernel=I_nch failed.');
assert(Info2.Rank == nCh2, 'Rank should equal nch in the underdetermined regime.');

fprintf('ALL TESTS PASSED: test_eigenmodes_transform_pure\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_eigenmodes_transform_pure.m`
Expected: FAIL — `Undefined function or variable 'bst_eigenmodes_transform'`.

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/math/bst_eigenmodes_transform.m`:

```matlab
function [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi, varargin)
% BST_EIGENMODES_TRANSFORM: Unregularized sensor->eigenmode spatial transform.
%
% USAGE:  [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi)
%         [Kernel, Info] = bst_eigenmodes_transform(Gain, Phi, 'Tol', tol)
%
% DESCRIPTION:
%     Builds the composite transform that maps sensor recordings directly to
%     Laplace-Beltrami eigenmode coefficients, with NO regularization. This is
%     the spatial analogue of a Fourier transform: a fixed change of basis with
%     no tuning parameters. Regularization/denoising is a separate, optional
%     step applied to the coefficients afterwards (a future filter library).
%
%     The compressed lead field is L_tilde = Gain * Phi, where each column is
%     the sensor topography of one eigenmode. The transform is the Moore-Penrose
%     pseudoinverse computed via SVD (L_tilde = U*S*V'):
%
%         Kernel = pinv(L_tilde) = V * diag(1./s) * U'      [K x nch]
%         Theta  = Kernel * Data                            [K x nTime]
%
%     SVD is used (rather than the normal equations) because the correct closed
%     form depends on the K/nch ratio -- the left-inverse (K<=nch) and the
%     right-inverse (K>=nch) each have a singular normal matrix in the other
%     regime -- and because forming Gain*Phi's Gram matrix would square an
%     already-large condition number (high-lambda modes are nearly invisible to
%     the sensors). Small singular values are floored (rank guard) but not
%     otherwise weighted: the transform is unregularized but rank-safe.
%
% INPUTS:
%     Gain : [nch x nVert] constrained (fixed-orientation) lead field, already
%            restricted to the channels to use.
%     Phi  : [nVert x K]   eigenmode matrix (caller truncates to K modes).
%
% OPTIONS (name-value):
%     'Tol' : singular-value floor (rank guard). Default [] -> MATLAB pinv
%             default max(size(L_tilde))*eps(max(s)).
%
% OUTPUTS:
%     Kernel : [K x nch] transform A = pinv(L_tilde).
%     Info   : struct with fields
%              .CompressedLF    [nch x K]  L_tilde = Gain*Phi
%              .SingularValues  [r x 1]    singular values of L_tilde
%              .Rank            scalar     number of singular values > Tol
%              .ConditionNumber scalar     s(1)/s(Rank)
%              .Tol             scalar     floor used
%              .nModes          scalar     K
%
% SEE ALSO: bst_eigenmodes_project, bst_inverse_eigenmodes, in_tess_eigenmodes

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

% Parse options
Tol = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'tol', Tol = varargin{i+1};
    end
end

Gain = double(Gain);
Phi  = double(Phi);

% Compressed lead field: column k = sensor topography of eigenmode k
L_tilde = Gain * Phi;                      % [nch x K]

% SVD-based pseudoinverse (correct for K<=nch and K>=nch; no Gram squaring)
[U, S, V] = svd(L_tilde, 'econ');
s = diag(S);

% Rank-safe singular-value floor (MATLAB pinv default if Tol not supplied)
if isempty(Tol)
    if isempty(s)
        Tol = 0;
    else
        Tol = max(size(L_tilde)) * eps(max(s));
    end
end
isKeep = (s > Tol);
sinv = zeros(size(s));
sinv(isKeep) = 1 ./ s(isKeep);

Kernel = V * diag(sinv) * U';              % [K x nch]

% Diagnostics
rankEff = sum(isKeep);
Info = struct();
Info.CompressedLF   = L_tilde;
Info.SingularValues = s;
Info.Rank           = rankEff;
if rankEff >= 1
    Info.ConditionNumber = s(1) / s(rankEff);
else
    Info.ConditionNumber = Inf;
end
Info.Tol    = Tol;
Info.nModes = size(Phi, 2);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_eigenmodes_transform_pure.m`
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_transform_pure`.

- [ ] **Step 5: Lint**

Run `checkcode` on `toolbox/math/bst_eigenmodes_transform.m` (MATLAB MCP `check_matlab_code`).
Expected: no errors (Brainstorm-idiom style warnings, if any, are acceptable).

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_eigenmodes_transform.m dev/tests/test_eigenmodes_transform_pure.m
git commit -m "$(cat <<'EOF'
Add bst_eigenmodes_transform: unregularized sensor->eigenmode transform

Pure SVD pseudoinverse A = pinv(L*Phi), correct in both K<=nch and
K>=nch regimes, rank-safe but unregularized. Pairs with
bst_eigenmodes_project. Covered by test_eigenmodes_transform_pure.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Process plugin `process_eigenmodes_transform`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_transform.m`
- Test: `dev/tests/test_process_eigenmodes_transform_options.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_process_eigenmodes_transform_options.m`:

```matlab
function test_process_eigenmodes_transform_options
% Verify the transform process exposes only transform options (no regularization
% knobs), and defaults vertex reconstruction off.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_transform('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(isequal(sProcess.InputTypes, {'data','raw'}), 'InputTypes must be {data,raw}.');

assert(isfield(sProcess.options, 'nmodes'),  'Missing nmodes option.');
assert(isfield(sProcess.options, 'dorecon'), 'Missing dorecon option.');
assert(strcmp(sProcess.options.dorecon.Type, 'checkbox'), 'dorecon must be a checkbox.');
assert(sProcess.options.dorecon.Value == 0, 'Vertex reconstruction must default OFF.');

% Transform is unregularized: it must NOT expose inverse-method knobs.
assert(~isfield(sProcess.options, 'method'),     'Transform must not expose a method option.');
assert(~isfield(sProcess.options, 'prioralpha'), 'Transform must not expose a prior option.');
assert(~isfield(sProcess.options, 'snr'),        'Transform must not expose an SNR option.');

fprintf('ALL TESTS PASSED: test_process_eigenmodes_transform_options\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_process_eigenmodes_transform_options.m`
Expected: FAIL — `Undefined function or variable 'process_eigenmodes_transform'`.

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/process/functions/process_eigenmodes_transform.m`:

```matlab
function varargout = process_eigenmodes_transform( varargin )
% PROCESS_EIGENMODES_TRANSFORM: Unregularized sensor->eigenmode spatial transform.
%
% USAGE:  sProcess = process_eigenmodes_transform('GetDescription')
%       OutputFiles = process_eigenmodes_transform('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Maps sensor recordings directly to eigenmode coefficients via the
%     unregularized transform A = pinv(L*Phi) (see bst_eigenmodes_transform).
%     The output is the raw eigenmode coefficient time series Theta [K x nTime]
%     as a Brainstorm matrix file -- the spatial analogue of a time series ready
%     for FFT. Feed that matrix file into "Frequency > FFT" (process_fft) to see
%     the joint (lambda, omega) spectrum. No regularization is applied; high-mode
%     coefficients are noisy by design.
%
%     Optionally also reconstructs raw vertex sources Q = Phi * Theta.
%     Requires precomputed eigenmodes on the cortex surface.
%
% SEE ALSO: bst_eigenmodes_transform, in_tess_eigenmodes, process_eigenmodes_inverse

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
    sProcess.Comment     = 'Eigenmode transform (spatial FFT)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    % === NUMBER OF MODES ===
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};

    % === OPTIONAL VERTEX RECONSTRUCTION ===
    sProcess.options.dorecon.Comment = 'Also reconstruct raw vertex sources (Q = Phi * Theta)';
    sProcess.options.dorecon.Type    = 'checkbox';
    sProcess.options.dorecon.Value   = 0;

    % === INFO LABEL ===
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Unregularized sensor&rarr;eigenmode transform ' ...
        '(A = pinv(L&middot;&Phi;)).<BR>Coefficients are raw and noisy at high modes by design &mdash; ' ...
        'run FFT on the matrix output to see the (&lambda;,&omega;) spectrum.<BR>' ...
        'Requires precomputed eigenmodes on the cortex surface.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    nModes = sProcess.options.nmodes.Value{1};
    if nModes > 0
        Comment = sprintf('Eigenmode transform (%d modes)', nModes);
    else
        Comment = 'Eigenmode transform (auto modes)';
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};

    nModesOpt = sProcess.options.nmodes.Value{1};
    DoRecon   = sProcess.options.dorecon.Value;

    % ===== STUDY + HEAD MODEL + SURFACE =====
    [sStudy, ~, ~, ~] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.');
        return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HeadModelMat  = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HeadModelMat.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputs, 'Eigenmode transform requires a surface head model.');
        return;
    end
    SurfaceFile = HeadModelMat.SurfaceFile;

    % ===== EIGENMODES =====
    [Eigenmodes, isComputed] = in_tess_eigenmodes(SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end
    nVertEigen = size(Eigenmodes.Vectors, 1);

    % ===== CONSTRAINED GAIN (fixed orientation) =====
    HM   = in_bst_headmodel(HeadModelFile, 1);   % ApplyOrient=1 -> [nch x nVert]
    Gain = double(HM.Gain);
    nVertHM = size(Gain, 2);
    if nVertEigen ~= nVertHM
        bst_report('Error', sProcess, sInputs, ...
            sprintf(['Head model has %d vertices but eigenmodes have %d vertices.\n' ...
            'Recompute the head model (right-click study > Compute head model).'], ...
            nVertHM, nVertEigen));
        return;
    end

    % ===== CHANNELS =====
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    if isempty(ChannelFile)
        bst_report('Error', sProcess, sInputs, 'No channel file found.');
        return;
    end
    ChannelMat   = in_bst_channel(ChannelFile);
    DataMat0     = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    nAllChannels = length(ChannelMat.Channel);
    if isfield(DataMat0, 'ChannelFlag') && ~isempty(DataMat0.ChannelFlag)
        ChannelFlag = DataMat0.ChannelFlag;
    else
        ChannelFlag = ones(nAllChannels, 1);
    end
    iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iMEG)
        iMEG = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG');
    end
    if isempty(iMEG)
        bst_report('Error', sProcess, sInputs, 'No good MEG or EEG channels found.');
        return;
    end

    % ===== BUILD TRANSFORM KERNEL =====
    nCh         = numel(iMEG);
    K_available = Eigenmodes.nModes;
    if isempty(nModesOpt) || nModesOpt <= 0
        K = min(nCh, K_available);
    else
        K = min(nModesOpt, K_available);
    end
    Phi     = double(Eigenmodes.Vectors(:, 1:K));
    lambdas = double(Eigenmodes.Values(1:K));

    [Kernel, Info] = bst_eigenmodes_transform(Gain(iMEG, :), Phi);   % [K x nCh]

    bst_report('Info', sProcess, sInputs, ...
        sprintf('Eigenmode transform: %d channels -> %d modes (rank %d, condition %.1f)', ...
        nCh, K, Info.Rank, Info.ConditionNumber));

    % Row descriptions (reused across matrix outputs)
    RowNames = cell(K, 1);
    for k = 1:K
        RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k));
    end

    % ===== PER-INPUT =====
    nInputs = numel(sInputs);
    for iInput = 1:nInputs
        sInput = sInputs(iInput);
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        DataMat = in_bst_data(sInput.FileName);
        isRaw   = isstruct(DataMat.F);

        if isRaw
            % Raw: cannot precompute coefficients; save kernel-only results.
            ResMat = db_template('resultsmat');
            ResMat.ImagingKernel = Phi * Kernel;          % [nVert x nCh]
            ResMat.ImageGridAmp  = [];
            ResMat.nComponents   = 1;
            ResMat.Comment       = sprintf('EigenTransform (%d modes) | %s', K, sInput.Comment);
            ResMat.Function      = 'eigentransform';
            ResMat.Time          = [];
            ResMat.DataFile      = sInput.FileName;
            ResMat.HeadModelFile = HeadModelFile;
            ResMat.HeadModelType = 'surface';
            ResMat.SurfaceFile   = SurfaceFile;
            ResMat.GoodChannel   = iMEG;
            ResMat.ChannelFlag   = ChannelFlag;
            ResMat.nAvg          = 1;
            ResMat.Leff          = 1;
            ResMat = bst_history('add', ResMat, 'eigenmodes_transform', ...
                sprintf('Unregularized eigenmode transform: %d modes, rank %d', K, Info.Rank));

            OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigentransform');
            bst_save(OutputFile, ResMat, 'v6');
            db_add_data(iStudyOut, OutputFile, ResMat);
            OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>

            bst_report('Warning', sProcess, sInput, ...
                'Raw file: saved kernel only. Import the data to compute eigenmode coefficients and their FFT.');
        else
            F     = double(DataMat.F(iMEG, :));   % [nCh x nTime]
            Theta = Kernel * F;                   % [K x nTime]

            % --- Coefficients matrix file (this is the FFT input) ---
            MatrixMat = db_template('matrixmat');
            MatrixMat.Value        = Theta;
            MatrixMat.Time         = DataMat.Time;
            MatrixMat.nAvg         = DataMat.nAvg;
            MatrixMat.Leff         = DataMat.Leff;
            MatrixMat.SurfaceFile  = SurfaceFile;
            MatrixMat.DisplayUnits = '';
            MatrixMat.Description  = RowNames;
            MatrixMat.Comment      = sprintf('EigenTransform (%d modes) | %s', K, sInput.Comment);
            MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_transform', ...
                sprintf('Unregularized eigenmode transform: %d modes, rank %d, condition %.1f', ...
                K, Info.Rank, Info.ConditionNumber));
            MatrixMat = bst_history('add', MatrixMat, 'eigenmodes_transform', ...
                sprintf('Input: %s', sInput.FileName));

            OutputFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigentransform');
            bst_save(OutputFile, MatrixMat, 'v6');
            db_add_data(iStudyOut, OutputFile, MatrixMat);
            OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>

            % --- Optional raw vertex reconstruction ---
            if DoRecon
                ResMat = db_template('resultsmat');
                ResMat.ImagingKernel = Phi * Kernel;       % [nVert x nCh]
                ResMat.ImageGridAmp  = [];
                ResMat.nComponents   = 1;
                ResMat.Comment       = sprintf('EigenTransform recon (%d modes) | %s', K, sInput.Comment);
                ResMat.Function      = 'eigentransform';
                ResMat.Time          = DataMat.Time;
                ResMat.DataFile      = sInput.FileName;
                ResMat.HeadModelFile = HeadModelFile;
                ResMat.HeadModelType = 'surface';
                ResMat.SurfaceFile   = SurfaceFile;
                ResMat.GoodChannel   = iMEG;
                ResMat.ChannelFlag   = ChannelFlag;
                ResMat.nAvg          = DataMat.nAvg;
                ResMat.Leff          = DataMat.Leff;
                ResMat = bst_history('add', ResMat, 'eigenmodes_transform', ...
                    sprintf('Raw vertex reconstruction from %d-mode transform', K));

                OutputFile = bst_process('GetNewFilename', StudyDir, 'results_eigentransform');
                bst_save(OutputFile, ResMat, 'v6');
                db_add_data(iStudyOut, OutputFile, ResMat);
                OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>
            end
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run (MATLAB MCP `run_matlab_test_file`): `dev/tests/test_process_eigenmodes_transform_options.m`
Expected: PASS — prints `ALL TESTS PASSED: test_process_eigenmodes_transform_options`.

- [ ] **Step 5: Lint**

Run `checkcode` on `toolbox/process/functions/process_eigenmodes_transform.m` (MATLAB MCP `check_matlab_code`).
Expected: no errors (Brainstorm-idiom warnings acceptable; the `%#ok<AGROW>`/`%#ok<DEFNU>` pragmas match existing process files).

- [ ] **Step 6: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_transform.m dev/tests/test_process_eigenmodes_transform_options.m
git commit -m "$(cat <<'EOF'
Add process_eigenmodes_transform: raw eigenmode coefficients for FFT

Drives bst_eigenmodes_transform: resolves head model/eigenmodes/
channels, applies the unregularized kernel to recordings, and saves the
raw eigenmode coefficient time series as matrix_eigentransform (FFT-ready
via process_fft). Optional raw vertex reconstruction. No regularization
options. Covered by test_process_eigenmodes_transform_options.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: End-to-end smoke on OMEGA (OPTIONAL — manual, does not gate the milestone)

This exercises the full chain on real data: transform → matrix file → `process_fft` →
`(λ, ω)` image. It requires the local OMEGA dataset and a Brainstorm DB, so it is
run manually, not as an automated unit test. The spec marks it optional.

**Files:** none created (uses existing data + processes).

- [ ] **Step 1: Build the substrate** (reuse the existing harness through the head-model step)

In MATLAB:

```matlab
% Anatomy -> recordings -> covariance -> overlapping-spheres head model for sub-0002.
% (Reuses the verified setup from the icosphere regression harness.)
res = test_omega_icosphere_sourcemap('ProtocolName', 'EigenTransformSmoke');
assert(res.pass, 'OMEGA substrate setup failed; cannot run smoke.');
```

- [ ] **Step 2: Compute eigenmodes on the cortex**

```matlab
% sub-0002-test is the icosphere arm created by the harness above.
bst_process('CallProcess', 'process_eigenmodes', [], [], ...
    'subjectname', 'sub-0002-test', ...
    'surftype',    'Cortex', ...
    'nmodes',      200, ...
    'masstype',    'barycentric', ...
    'removedc',    1, ...
    'repair',      0, ...
    'overwrite',   1);
```

- [ ] **Step 3: Import a short data block, then run the transform**

```matlab
% Find the resting (raw) link for sub-0002-test, import 0-10 s, then transform.
sRaw = bst_process('CallProcess', 'process_select_files_data', [], [], ...
    'subjectname', 'sub-0002-test');
sImp = bst_process('CallProcess', 'process_import_data_time', sRaw(1), [], ...
    'timewindow', [0 10], 'split', 0, 'usectfcomp', 1, 'usessp', 1);
sTrans = bst_process('CallProcess', 'process_eigenmodes_transform', sImp, [], ...
    'nmodes', 0, 'dorecon', 0);
assert(~isempty(sTrans) && ~isempty(sTrans(1).FileName), 'Transform produced no output.');
```

- [ ] **Step 4: FFT the coefficients via the existing process**

```matlab
sFft = bst_process('CallProcess', 'process_fft', sTrans(1), [], ...
    'units', 'physical', 'avgoutput', 0);
assert(~isempty(sFft) && ~isempty(sFft(1).FileName), 'process_fft produced no output.');
disp(sFft(1).FileName);
```

- [ ] **Step 5: View the (λ, ω) spectrum**

Open `sFft(1).FileName` in Brainstorm (double-click the spectrum file). In the
spectrum viewer's image mode the rows are eigenmodes (ordered by `λₖ`, i.e.
spatial frequency) and the columns are temporal frequency `ω`. Confirm by eye:
low-`k` modes show banded structure; high-`k` modes wash into a noise floor —
the expected diagnostic picture. No assertion; this is a visual check.

---

## Self-Review

**Spec coverage:**
- Pure transform `A = pinv(L̃)` via SVD → Task 1. ✓
- Apply to recordings → `matrix_eigentransform [K × nTime]` → Task 2 (imported-data branch). ✓
- FFT via existing `process_fft` on a matrix file → Task 3 Step 4. ✓
- Visualize `(λ, ω)` in standard viewer → Task 3 Step 5. ✓
- Default `K = min(nch, available)` → Task 2 Run (`K = min(nCh, K_available)`). ✓
- Optional vertex reconstruction off by default → Task 2 (`dorecon` checkbox, default 0). ✓
- `bst_inverse_eigenmodes` untouched → no task modifies it. ✓
- Pure unit test + checkcode → Tasks 1 & 2 Steps 1–5. ✓
- Optional manual OMEGA smoke → Task 3. ✓
- Deferred items (noise floor, filter library, agreement testing, vector work) → no tasks, as intended. ✓

**Placeholder scan:** No TBD/TODO; every code step contains complete code; every run step has an exact invocation + expected result. ✓

**Type/name consistency:** `bst_eigenmodes_transform(Gain, Phi)` returns `[Kernel, Info]` with `Info.Rank`, `Info.ConditionNumber`, `Info.SingularValues`, `Info.nModes` — used identically in Task 1's test and Task 2's `Run`. File prefixes `matrix_eigentransform` / `results_eigentransform` and `Function='eigentransform'` are consistent across branches. ✓
