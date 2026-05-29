# Mode-Frequency Spectrum Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add two Brainstorm processes that compute the per-eigenmode temporal-frequency spectrum (PSD and FFT) of source maps, by feeding the **mode kernel** (`Φ'·M·ImagingKernel`) into the same `bst_psd` call the standard source PSD/FFT uses.

**Architecture:** A pure math helper builds the mode kernel; a shared engine reads the windowed sensor data + bad segments exactly as `bst_timefreq` does, then calls `bst_psd` with the mode kernel substituted for the imaging kernel (valid because `bst_psd` applies the kernel after the FFT, by linearity), and saves a standard `timefreq` file. Two thin sibling processes (PSD, FFT) supply the option panels and delegate to the engine.

**Tech Stack:** MATLAB; Brainstorm process framework (`macro_method` dispatch); `bst_psd`, `in_bst_results`, `in_tess_eigenmodes`, `tess_laplacian`, `db_template('timefreqmat')`; MATLAB MCP for running tests.

**Spec:** `docs/superpowers/specs/2026-05-29-mode-frequency-spectrum-design.md`

**Key facts verified against the codebase (do not re-derive):**
- `bst_psd(F, sfreq, WinLength, WinOverlap, BadSegments, ImagingKernel, WinFunc, PowerUnits, IsRelative)` returns `[TF, FreqVector, Nwin, Messages, TFbis]`. `TF` is `[nRows × 1 × nFreq]` **power**. It applies the kernel *after* the FFT (`bst_psd.m:157`) and skips bad-segment windows *before* the kernel (`bst_psd.m:113`). `WinLength=[]` ⇒ single whole-window transform (used for FFT).
- `bst_eigenmodes_project` computes `Phi' * (MassMatrix * Data)`.
- `Eigenmodes` struct (from `in_tess_eigenmodes(SurfaceFile)` → `[Eigenmodes, isComputed]`) has fields `.Vectors [nV×nModes]`, `.Values [nModes×1]`, `.MassMatrix [nV×nV] sparse`, `.MassType`, `.Component`, `.CompRank`.
- Save pattern: `db_template('timefreqmat')`, fill fields, `OutputFile = bst_process('GetNewFilename', StudyDir, ['timefreq_' Method])`, `bst_save(OutputFile, FileMat, 'v6')`, `db_add_data(iStudy, OutputFile, FileMat)`.
- Tests live in `dev/tests/*.m`, start Brainstorm with `if ~brainstorm('status'), brainstorm nogui; end`, and print `ALL TESTS PASSED: <name>` on success.

**Testing strategy:** The mathematical core (mode kernel + the `bst_psd` kernel-swap, including bad-segment handling) is fully covered by **pure tests that need no database** (Tasks 1–2). The DB/raw I/O path in the engine (`in_bst_results`, raw `in_fread`, `db_add_data`) is exercised by **manual GUI validation** (Task 6), per the established workflow — building a synthetic raw protocol in an automated test is out of scope for this round.

---

## File Structure

- `toolbox/math/bst_eigenmodes_modekernel.m` — **new, pure.** Builds the mode kernel `Φ(:,1:K)'·M·ImagingKernel` (or the projector `Φ(:,1:K)'·M` when no imaging kernel). No DB, no I/O.
- `toolbox/process/functions/process_eigenmodes_freq.m` — **new, shared engine** (no `GetDescription`, not in the menu). Exposes `Run(sProcess, sInputs, Method)` plus a **pure** `Compute(...)` subfunction (the spectral core, DB-free and unit-testable) and the raw/imported read helpers.
- `toolbox/process/functions/process_eigenmodes_psd.m` — **new, listed process.** PSD option panel; delegates to the engine.
- `toolbox/process/functions/process_eigenmodes_fft.m` — **new, listed process.** FFT option panel; delegates to the engine.
- `dev/tests/test_bst_eigenmodes_modekernel_pure.m` — **new.** Pure unit test for the helper.
- `dev/tests/test_eigenmodes_freq_compute.m` — **new.** Pure test for the engine's `Compute` (shape, freqs, row labels, and kernel-swap ≡ project-first parity incl. bad segments).
- `dev/tests/test_eigenmodes_freq_processes.m` — **new.** Validates both processes' `GetDescription`/`FormatComment` (no DB).

---

## Task 1: Mode kernel helper

**Files:**
- Create: `toolbox/math/bst_eigenmodes_modekernel.m`
- Test: `dev/tests/test_bst_eigenmodes_modekernel_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_bst_eigenmodes_modekernel_pure.m`:

```matlab
function test_bst_eigenmodes_modekernel
% Pure test: mode kernel equals the explicit projection, plus capping/no-kernel.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'math'));

rng(7);
nV = 40; nModesAll = 12; nChan = 25; nT = 9;
Phi = randn(nV, nModesAll);
M   = spdiags(rand(nV,1) + 0.5, 0, nV, nV);   % SPD diagonal mass matrix
Eig = struct('Vectors', Phi, 'Values', (1:nModesAll)');
K   = randn(nV, nChan);                         % imaging kernel [nV x nChan]
F   = randn(nChan, nT);                         % sensor data

% ---- Kernel case: ModeKernel*F == Phi(:,1:nModes)' * M * (K*F) ----
nModes = 7;
MK = bst_eigenmodes_modekernel(Eig, M, K, nModes);
assert(isequal(size(MK), [nModes, nChan]), 'ModeKernel size wrong');
Expected = Phi(:,1:nModes)' * (M * (K * F));
assert(max(abs(MK*F - Expected), [], 'all') < 1e-9, 'kernel-folded projection mismatch');

% ---- No-kernel case: returns the projector Phi(:,1:nModes)' * M ----
P = bst_eigenmodes_modekernel(Eig, M, [], nModes);
assert(isequal(size(P), [nModes, nV]), 'projector size wrong');
assert(max(abs(P - Phi(:,1:nModes)' * M), [], 'all') < 1e-9, 'projector mismatch');

% ---- Capping: requesting more than available uses all; first-K consistency ----
MKall = bst_eigenmodes_modekernel(Eig, M, K, 999);
assert(size(MKall,1) == nModesAll, 'cap to available modes failed');
assert(max(abs(MKall(1:nModes,:) - MK), [], 'all') < 1e-12, 'first-K rows must match');

fprintf('ALL TESTS PASSED: test_bst_eigenmodes_modekernel\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run via MATLAB MCP `run_matlab_test_file` on `dev/tests/test_bst_eigenmodes_modekernel_pure.m`.
Expected: FAIL — `Unrecognized function or variable 'bst_eigenmodes_modekernel'` (or "not found").

- [ ] **Step 3: Write the implementation**

Create `toolbox/math/bst_eigenmodes_modekernel.m`:

```matlab
function ModeKernel = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)
% BST_EIGENMODES_MODEKERNEL: Build the mode kernel mapping signals -> mode coefficients.
%
% USAGE:  ModeKernel = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)
%         Projector  = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, [], nModes)
%
% DESCRIPTION:
%     Folds the M-weighted eigenmode projection into the imaging kernel, so that
%     applying the result to sensor data yields mode coefficients directly:
%
%         ModeKernel = Phi(:,1:nModes)' * M * ImagingKernel       % [nModes x nChannels]
%         coeff(t)   = ModeKernel * F(t)                          % = Phi'*M*(K*F(t))
%
%     With an empty ImagingKernel, returns the bare projector Phi(:,1:nModes)' * M
%     ([nModes x nVertices]) to apply to a full source matrix (ImageGridAmp).
%
% INPUTS:
%     Eigenmodes    : struct with field .Vectors [nVertices x nModesAll]
%     MassMatrix    : [nVertices x nVertices] sparse mass matrix
%     ImagingKernel : [nVertices x nChannels] imaging kernel, or [] for full sources
%     nModes        : number of leading modes to keep (clamped to available count)
%
% SEE ALSO: bst_eigenmodes_project, tess_laplacian, in_tess_eigenmodes

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

Phi = double(Eigenmodes.Vectors);     % [nV x nModesAll]
nModesAll = size(Phi, 2);
if (nargin < 4) || isempty(nModes)
    nModes = nModesAll;
end
nModes = max(1, min(nModes, nModesAll));
Phi = Phi(:, 1:nModes);               % [nV x nModes]

% Projector Phi' * M  ->  [nModes x nV]
Projector = Phi' * MassMatrix;

if (nargin < 3) || isempty(ImagingKernel)
    ModeKernel = Projector;                       % [nModes x nV]
else
    ModeKernel = Projector * ImagingKernel;       % [nModes x nChannels]
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run `dev/tests/test_bst_eigenmodes_modekernel_pure.m` via MATLAB MCP.
Expected: PASS — prints `ALL TESTS PASSED: test_bst_eigenmodes_modekernel`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_modekernel.m dev/tests/test_bst_eigenmodes_modekernel_pure.m
git commit -m "Mode-frequency: add bst_eigenmodes_modekernel (Phi'*M*Kernel) + test

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Engine spectral core (`Compute`) + raw/imported readers + `Run`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_freq.m`
- Test: `dev/tests/test_eigenmodes_freq_compute.m`

This task creates the engine file with three parts: the **pure `Compute`** (tested now), the **read helpers**, and the **`Run`** orchestrator (validated manually in Task 6). All three are written in this task; only `Compute` is unit-tested here.

- [ ] **Step 1: Write the failing test (covers the pure `Compute`)**

Create `dev/tests/test_eigenmodes_freq_compute.m`:

```matlab
function test_eigenmodes_freq_compute
% Pure test of the engine's spectral core: shape, freqs, labels, and the
% kernel-swap == project-first parity (including bad-segment handling).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'process', 'functions'));

rng(11);
nChan = 18; nModes = 5; sfreq = 200; nT = 4000;
F  = randn(nChan, nT);
MK = randn(nModes, nChan);            % mode kernel [nModes x nChan]
Values = (1:nModes)';

% ---- PSD: shape, frequency vector, row labels ----
[TF, Freqs, RowNames] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, [], 'psd', 1, 50, 'mean', 'physical', Values);
assert(ndims(TF) == 3 && size(TF,1) == nModes && size(TF,2) == 1, 'PSD TF shape wrong');
assert(size(TF,3) == numel(Freqs), 'PSD freq dim mismatch');
assert(issorted(Freqs) && all(Freqs >= 0), 'freqs must be ascending and non-negative');
assert(numel(RowNames) == nModes && ischar(RowNames{1}), 'row labels wrong');
assert(all(isfinite(TF(:))) && all(TF(:) >= 0), 'PSD power must be finite and non-negative');

% ---- Parity (no bad segments): kernel-swap == project-first ----
TFref = bst_psd(MK*F, sfreq, 1, 50, [], [], 'mean', 'physical');
assert(max(abs(TF(:) - TFref(:))) < 1e-9, 'PSD kernel-swap != project-first');

% ---- Parity WITH a bad segment that drops a window ----
BadSeg = [round(1.2*sfreq); round(1.8*sfreq)];   % [start; stop] samples
[TFb] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, BadSeg, 'psd', 1, 50, 'mean', 'physical', Values);
TFbref = bst_psd(MK*F, sfreq, 1, 50, BadSeg, [], 'mean', 'physical');
assert(max(abs(TFb(:) - TFbref(:))) < 1e-9, 'bad-segment handling differs from project-first');
assert(max(abs(TFb(:) - TF(:))) > 0, 'bad segment should change the result');

% ---- FFT: single-window transform, same parity ----
[TF2, Freqs2] = process_eigenmodes_freq('Compute', ...
    F, sfreq, MK, [], 'fft', [], [], [], 'physical', Values);
assert(size(TF2,1) == nModes && size(TF2,2) == 1 && size(TF2,3) == numel(Freqs2), 'FFT shape wrong');
TF2ref = bst_psd(MK*F, sfreq, [], 0, [], [], [], 'physical');
assert(max(abs(TF2(:) - TF2ref(:))) < 1e-9, 'FFT kernel-swap != project-first');

fprintf('ALL TESTS PASSED: test_eigenmodes_freq_compute\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run `dev/tests/test_eigenmodes_freq_compute.m` via MATLAB MCP.
Expected: FAIL — `process_eigenmodes_freq` not found / no `Compute` method.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_freq.m`:

```matlab
function varargout = process_eigenmodes_freq( varargin )
% PROCESS_EIGENMODES_FREQ: Shared engine for the mode-frequency spectrum (PSD/FFT).
%
% Not listed in the process menu (no GetDescription). Driven by the sibling
% processes process_eigenmodes_psd / process_eigenmodes_fft.
%
% USAGE:  OutputFiles = process_eigenmodes_freq('Run', sProcess, sInputs, Method)
%         [TF, Freqs, RowNames, Msg] = process_eigenmodes_freq('Compute', ...
%               F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Values)
%
% DESCRIPTION:
%     Computes the per-eigenmode temporal-frequency spectrum of a source map by
%     substituting the mode kernel (Phi'*M*ImagingKernel) for the imaging kernel
%     in the standard bst_psd call. Because bst_psd applies the kernel after the
%     FFT (valid by linearity) and skips bad-segment windows before the kernel,
%     the result is identical to a source PSD/FFT with Phi'*M applied.
%
% SEE ALSO: bst_eigenmodes_modekernel, bst_psd, process_eigenmodes_psd, process_eigenmodes_fft

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


%% ===== COMPUTE (pure spectral core) =====
function [TF, Freqs, RowNames, Msg] = Compute(F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Values) %#ok<DEFNU>
    nModes = size(ModeKernel, 1);
    switch lower(Method)
        case 'psd'
            [TF, Freqs, ~, Msg] = bst_psd(F, sfreq, WinLength, WinOverlap, BadSegments, ModeKernel, WinFunc, PowerUnits);
        case 'fft'
            [TF, Freqs, ~, Msg] = bst_psd(F, sfreq, [], 0, BadSegments, ModeKernel, [], PowerUnits);
        otherwise
            error('Unsupported method: %s', Method);
    end
    % Row labels: mode index + eigenvalue
    RowNames = cell(nModes, 1);
    for k = 1:nModes
        RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, Values(k));
    end
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs, Method) %#ok<DEFNU>
    OutputFiles = {};

    % ===== OPTIONS =====
    nModesReq = sProcess.options.nmodes.Value{1};
    PowerUnits = 'physical';
    if isfield(sProcess.options, 'units') && ~isempty(sProcess.options.units.Value)
        PowerUnits = sProcess.options.units.Value;
    end
    TimeWindow = [];
    if isfield(sProcess.options, 'timewindow') && ~isempty(sProcess.options.timewindow.Value) && iscell(sProcess.options.timewindow.Value)
        TimeWindow = sProcess.options.timewindow.Value{1};
    end
    if strcmpi(Method, 'psd')
        WinLength  = sProcess.options.win_length.Value{1};
        WinOverlap = sProcess.options.win_overlap.Value{1};
        switch (sProcess.options.win_std.Value)
            case {0, 'mean'}, WinFunc = 'mean';
            case {1, 'std'},  WinFunc = 'std';
            otherwise,        WinFunc = 'mean';
        end
    else
        WinLength = []; WinOverlap = 0; WinFunc = 'mean';
    end

    bst_progress('start', 'Mode-frequency spectrum', 'Computing...');

    % ===== PER INPUT FILE =====
    for iInput = 1:length(sInputs)
        sInput = sInputs(iInput);

        % --- Load metadata only (no kernel expansion) ---
        ResultsMat = in_bst_results(sInput.FileName, 0, ...
            'ImagingKernel', 'ImageGridAmp', 'GoodChannel', 'nComponents', ...
            'DataFile', 'Time', 'SurfaceFile', 'HeadModelType', 'Atlas');

        % --- Validate ---
        if isfield(ResultsMat, 'HeadModelType') && ~isempty(ResultsMat.HeadModelType) && ~strcmpi(ResultsMat.HeadModelType, 'surface')
            bst_report('Error', sProcess, sInput, 'Only surface source models are supported.'); continue;
        end
        if isfield(ResultsMat, 'Atlas') && ~isempty(ResultsMat.Atlas)
            bst_report('Error', sProcess, sInput, 'Atlas-based source models are not supported.'); continue;
        end
        if ResultsMat.nComponents ~= 1
            bst_report('Error', sProcess, sInput, 'Constrained source orientation required (nComponents == 1).'); continue;
        end
        SurfaceFile = ResultsMat.SurfaceFile;
        if isempty(SurfaceFile)
            bst_report('Error', sProcess, sInput, 'No surface file associated with the source map.'); continue;
        end

        % --- Eigenmodes ---
        [Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on: ' SurfaceFile '. Run "Compute eigenmodes" first.']); continue;
        end
        nV_eigen = size(Eig.Vectors, 1);

        % --- Mass matrix (reuse if saved) ---
        if isfield(Eig, 'MassMatrix') && ~isempty(Eig.MassMatrix)
            M = Eig.MassMatrix;
        else
            sSurf = in_tess_bst(SurfaceFile, 0);
            [~, M] = tess_laplacian(sSurf.Vertices, sSurf.Faces, 'MassType', Eig.MassType);
        end

        % --- Mode count ---
        nModes = max(1, min(nModesReq, size(Eig.Vectors, 2)));
        if nModes < nModesReq
            bst_report('Warning', sProcess, sInput, sprintf('Requested %d modes but only %d available; using %d.', nModesReq, nModes, nModes));
        end

        % --- Read sensor/source data + bad segments (mirrors bst_timefreq) ---
        [F, TimeVector, BadSegments, errMsg] = ReadInput(sInput, ResultsMat, TimeWindow);
        if ~isempty(errMsg)
            bst_report('Error', sProcess, sInput, errMsg); continue;
        end
        sfreq = 1 ./ (TimeVector(2) - TimeVector(1));

        % --- Build the kernel passed to bst_psd ---
        if isempty(ResultsMat.ImagingKernel)
            % Full result: F are sources; check vertex match; kernel = projector
            if size(F,1) ~= nV_eigen
                bst_report('Error', sProcess, sInput, sprintf('Vertex mismatch: sources %d, eigenmodes %d.', size(F,1), nV_eigen)); continue;
            end
            ModeKernel = bst_eigenmodes_modekernel(Eig, M, [], nModes);            % [nModes x nV]
        else
            % Kernel result: F are sensors; kernel = Phi'*M*ImagingKernel
            ModeKernel = bst_eigenmodes_modekernel(Eig, M, ResultsMat.ImagingKernel, nModes);  % [nModes x nChan]
        end

        % --- Compute the spectrum (kernel-swap) ---
        [TF, Freqs, RowNames, Msg] = Compute(F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Eig.Values);
        if isempty(TF)
            bst_report('Error', sProcess, sInput, ['Spectral computation failed: ' Msg]); continue;
        end

        % --- Build timefreq file ---
        FileMat = db_template('timefreqmat');
        FileMat.TF          = TF;
        FileMat.Time        = TimeVector([1 end]);
        FileMat.Freqs       = Freqs;
        FileMat.RowNames    = RowNames;
        FileMat.Measure     = 'power';
        FileMat.Method      = Method;
        FileMat.DataType    = 'matrix';
        FileMat.SurfaceFile = SurfaceFile;
        FileMat.DataFile    = '';
        FileMat.nAvg        = 1;
        FileMat.Leff        = 1;
        switch lower(Method)
            case 'psd', FileMat.Comment = sprintf('Mode PSD [%d modes]', nModes);
            case 'fft', FileMat.Comment = sprintf('Mode FFT [%d modes]', nModes);
        end
        FileMat.Options.Method     = Method;
        FileMat.Options.Measure    = 'power';
        FileMat.Options.PowerUnits = PowerUnits;
        if strcmpi(Method, 'psd')
            FileMat.Options.WindowFunction = WinFunc;
        end
        FileMat = bst_history('add', FileMat, 'compute', sprintf('Mode-frequency spectrum (%s) of: %s', Method, sInput.FileName));
        FileMat = bst_history('add', FileMat, 'compute', sprintf('%d modes, %s mass, eigenvalue range [%.2f, %.2f]', nModes, Eig.MassType, Eig.Values(1), Eig.Values(nModes)));

        % --- Save ---
        [sStudy, iStudy] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
        OutputFile = bst_process('GetNewFilename', StudyDir, ['timefreq_' Method]);
        bst_save(OutputFile, FileMat, 'v6');
        db_add_data(iStudy, OutputFile, FileMat);
        OutputFiles{end+1} = file_short(OutputFile); %#ok<AGROW>
    end

    bst_progress('stop');
end


%% ===== READ INPUT (sensor/source data + bad segments) =====
% Mirrors the results branch of bst_timefreq (bst_timefreq.m:346-440) so that
% bad-segment handling and the time-window read are identical to a source PSD.
function [F, TimeVector, BadSegments, errMsg] = ReadInput(sInput, ResultsMat, TimeWindow)
    F = []; TimeVector = []; BadSegments = []; errMsg = '';

    if ~isempty(ResultsMat.ImagingKernel) && isempty(ResultsMat.ImageGridAmp)
        % ----- Kernel result: process the recordings file -----
        DataFile = ResultsMat.DataFile;
        if isempty(DataFile)
            errMsg = 'Kernel result requires an associated data file.'; return;
        end
        sMat = in_bst_data(DataFile);
        if isstruct(sMat.F)
            % Raw recordings: bounded windowed read
            sFile = sMat.F;
            ChannelFile = bst_get('ChannelFileForStudy', sInput.FileName);
            if isempty(ChannelFile)
                errMsg = 'No channel definition available for this file.'; return;
            end
            ChannelMat = in_bst_channel(ChannelFile);
            if (length(sFile.epochs) > 1)
                errMsg = 'Files with epochs are not supported by this process.'; return;
            end
            [F, TimeVector, BadSegments] = ReadRawRecordings(sFile, sMat.Time, ChannelMat, TimeWindow);
        else
            % Imported recordings
            F = sMat.F;
            TimeVector = sMat.Time;
            sMat.events = sMat.Events;
            sMat.prop.sfreq = 1 ./ (sMat.Time(2) - sMat.Time(1));
            isChannelEvtBad = 0;
            BadSegments = panel_record('GetBadSegments', sMat, isChannelEvtBad) - sMat.prop.sfreq * sMat.Time(1) + 1;
        end
        % Restrict to this result's good channels
        F = F(ResultsMat.GoodChannel, :);
    else
        % ----- Full result: sources themselves -----
        F = ResultsMat.ImageGridAmp;
        TimeVector = ResultsMat.Time;
        BadSegments = [];
    end

    % ----- Keep only the requested time window (mirrors bst_timefreq.m:425-440) -----
    if ~isempty(TimeWindow)
        iTime = bst_closest(TimeWindow, TimeVector);
        if (iTime(1) == iTime(2))
            errMsg = 'Selected time window is not valid for the input file.'; return;
        end
        iTime = iTime(1):iTime(2);
        TimeVector = TimeVector(iTime);
        F = F(:, iTime);
    end
end


%% ===== READ RAW RECORDINGS (bounded window) =====
% Copy of bst_timefreq's private ReadRawRecordings (bst_timefreq.m:960), with
% TimeWindow passed in explicitly. PSD bad-segment handling preserved.
function [F, TimeVector, BadSegments] = ReadRawRecordings(sFile, TimeVector, ChannelMat, TimeWindow)
    ImportOptions = db_template('ImportOptions');
    ImportOptions.ImportMode     = 'Time';
    ImportOptions.Resample       = 0;
    ImportOptions.UseCtfComp     = 1;
    ImportOptions.UseSsp         = 1;
    ImportOptions.RemoveBaseline = 'no';
    ImportOptions.DisplayMessages = 0;
    % Samples to read
    if ~isempty(TimeWindow)
        SamplesBounds = round(sFile.prop.times(1) .* sFile.prop.sfreq) + bst_closest(TimeWindow, TimeVector) - 1;
    else
        SamplesBounds = round(sFile.prop.times .* sFile.prop.sfreq);
    end
    % Read data
    [F, TimeVector] = in_fread(sFile, ChannelMat, 1, SamplesBounds, [], ImportOptions);
    % Bad segments (relative to the start of the read section)
    isChannelEvtBad = 0;
    BadSegments = panel_record('GetBadSegments', sFile, isChannelEvtBad);
    if ~isempty(BadSegments)
        BadSegments = BadSegments - SamplesBounds(1) + 1;
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run `dev/tests/test_eigenmodes_freq_compute.m` via MATLAB MCP.
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_freq_compute`.

- [ ] **Step 5: Run the Task 1 test again to confirm no regression**

Run `dev/tests/test_bst_eigenmodes_modekernel_pure.m`.
Expected: PASS.

- [ ] **Step 6: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_freq.m dev/tests/test_eigenmodes_freq_compute.m
git commit -m "Mode-frequency: add shared engine (Compute core + raw/imported read + Run)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: PSD process

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_psd.m`
- Test: `dev/tests/test_eigenmodes_freq_processes.m` (shared with Task 4; create here, extend in Task 4)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_freq_processes.m`:

```matlab
function test_eigenmodes_freq_processes
% Validate GetDescription/FormatComment of the mode-frequency processes (no DB).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
addpath(fullfile(repoRoot, 'toolbox', 'process', 'functions'));

% ---- PSD process ----
sPsd = process_eigenmodes_psd('GetDescription');
assert(~isempty(sPsd.Comment), 'PSD comment empty');
assert(strcmpi(sPsd.InputTypes{1}, 'results'), 'PSD input must be results');
assert(strcmpi(sPsd.OutputTypes{1}, 'timefreq'), 'PSD output must be timefreq');
assert(isfield(sPsd.options, 'nmodes') && isequal(sPsd.options.nmodes.Value{1}, 300), 'PSD nmodes default 300');
assert(isfield(sPsd.options, 'win_length'), 'PSD must expose window length');
assert(isfield(sPsd.options, 'win_overlap'), 'PSD must expose window overlap');
assert(isfield(sPsd.options, 'units'), 'PSD must expose units');
assert(ischar(process_eigenmodes_psd('FormatComment', sPsd)), 'PSD FormatComment must return char');

fprintf('ALL TESTS PASSED: test_eigenmodes_freq_processes\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run `dev/tests/test_eigenmodes_freq_processes.m` via MATLAB MCP.
Expected: FAIL — `process_eigenmodes_psd` not found.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_psd.m`:

```matlab
function varargout = process_eigenmodes_psd( varargin )
% PROCESS_EIGENMODES_PSD: Per-eigenmode power spectrum density (Welch) of source maps.
%
% USAGE:  sProcess = process_eigenmodes_psd('GetDescription')
%       OutputFiles = process_eigenmodes_psd('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Computes the PSD of each Laplace-Beltrami mode coefficient time series by
%     substituting the mode kernel for the imaging kernel in the standard Welch
%     PSD. Output is a timefreq file: rows = modes, x = temporal frequency.
%     Requires precomputed eigenmodes on the source surface. Constrained sources.
%
% SEE ALSO: process_eigenmodes_freq, process_eigenmodes_fft, process_psd

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
    sProcess.Comment     = 'Mode-frequency spectrum: PSD (Welch)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338.5;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Number of modes
    sProcess.options.nmodes.Comment = 'Number of modes: ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {300, '', 0};
    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Window length
    sProcess.options.win_length.Comment = 'Window length: ';
    sProcess.options.win_length.Type    = 'value';
    sProcess.options.win_length.Value   = {1, 's', []};
    % Window overlap
    sProcess.options.win_overlap.Comment = 'Window overlap ratio: ';
    sProcess.options.win_overlap.Type    = 'value';
    sProcess.options.win_overlap.Value   = {50, '%', 1};
    % Units
    sProcess.options.units.Comment = {'Physical: U<SUP>2</SUP>/Hz', 'Normalized: U<SUP>2</SUP>/Hz/s', 'Before Nov 2020', 'Units:'; ...
                                      'physical', 'normalized', 'old', ''};
    sProcess.options.units.Type    = 'radio_linelabel';
    sProcess.options.units.Value   = 'physical';
    % Save std across windows
    sProcess.options.win_std.Comment = 'Save the std across windows instead of the mean';
    sProcess.options.win_std.Type    = 'checkbox';
    sProcess.options.win_std.Value   = 0;
    % Info
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Output: timefreq file, rows = modes, x = frequency.<BR>' ...
                                           'Requires precomputed eigenmodes; constrained sources.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = process_eigenmodes_freq('Run', sProcess, sInputs, 'psd');
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run `dev/tests/test_eigenmodes_freq_processes.m` via MATLAB MCP.
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_freq_processes`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_psd.m dev/tests/test_eigenmodes_freq_processes.m
git commit -m "Mode-frequency: add PSD process (delegates to engine)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: FFT process

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_fft.m`
- Modify: `dev/tests/test_eigenmodes_freq_processes.m` (add FFT assertions)

- [ ] **Step 1: Extend the test (failing for FFT)**

In `dev/tests/test_eigenmodes_freq_processes.m`, insert before the final `fprintf`:

```matlab
% ---- FFT process ----
sFft = process_eigenmodes_fft('GetDescription');
assert(~isempty(sFft.Comment), 'FFT comment empty');
assert(strcmpi(sFft.InputTypes{1}, 'results'), 'FFT input must be results');
assert(strcmpi(sFft.OutputTypes{1}, 'timefreq'), 'FFT output must be timefreq');
assert(isfield(sFft.options, 'nmodes') && isequal(sFft.options.nmodes.Value{1}, 300), 'FFT nmodes default 300');
assert(isfield(sFft.options, 'units'), 'FFT must expose units');
assert(~isfield(sFft.options, 'win_length'), 'FFT must not expose window length');
assert(ischar(process_eigenmodes_fft('FormatComment', sFft)), 'FFT FormatComment must return char');
```

- [ ] **Step 2: Run the test to verify it fails**

Run `dev/tests/test_eigenmodes_freq_processes.m` via MATLAB MCP.
Expected: FAIL — `process_eigenmodes_fft` not found.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_fft.m`:

```matlab
function varargout = process_eigenmodes_fft( varargin )
% PROCESS_EIGENMODES_FFT: Per-eigenmode Fourier amplitude spectrum of source maps.
%
% USAGE:  sProcess = process_eigenmodes_fft('GetDescription')
%       OutputFiles = process_eigenmodes_fft('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Computes the single-window FFT of each Laplace-Beltrami mode coefficient
%     time series by substituting the mode kernel for the imaging kernel in the
%     standard FFT. Output is a timefreq file: rows = modes, x = frequency.
%     Requires precomputed eigenmodes on the source surface. Constrained sources.
%
% SEE ALSO: process_eigenmodes_freq, process_eigenmodes_psd, process_fft

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
    sProcess.Comment     = 'Mode-frequency spectrum: FFT';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338.6;
    sProcess.Description = '';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Number of modes
    sProcess.options.nmodes.Comment = 'Number of modes: ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {300, '', 0};
    % Time window
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    % Units
    sProcess.options.units.Comment = {'Physical: U<SUP>2</SUP>/Hz', 'Normalized: U<SUP>2</SUP>/Hz/s', 'Before Nov 2020', 'Units:'; ...
                                      'physical', 'normalized', 'old', ''};
    sProcess.options.units.Type    = 'radio_linelabel';
    sProcess.options.units.Value   = 'physical';
    % Info
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Output: timefreq file, rows = modes, x = frequency.<BR>' ...
                                           'Requires precomputed eigenmodes; constrained sources.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = process_eigenmodes_freq('Run', sProcess, sInputs, 'fft');
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run `dev/tests/test_eigenmodes_freq_processes.m` via MATLAB MCP.
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_freq_processes`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_fft.m dev/tests/test_eigenmodes_freq_processes.m
git commit -m "Mode-frequency: add FFT process (delegates to engine)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Full test sweep + manual GUI validation checklist

**Files:** none (verification only)

- [ ] **Step 1: Run all four mode-frequency tests**

Run each via MATLAB MCP and confirm each prints its `ALL TESTS PASSED` line:
- `dev/tests/test_bst_eigenmodes_modekernel_pure.m`
- `dev/tests/test_eigenmodes_freq_compute.m`
- `dev/tests/test_eigenmodes_freq_processes.m`

- [ ] **Step 2: Static check the new files**

Run MATLAB MCP `check_matlab_code` on each new `.m` file and confirm no errors:
- `toolbox/math/bst_eigenmodes_modekernel.m`
- `toolbox/process/functions/process_eigenmodes_freq.m`
- `toolbox/process/functions/process_eigenmodes_psd.m`
- `toolbox/process/functions/process_eigenmodes_fft.m`

- [ ] **Step 3: Write the manual GUI validation checklist to the spec folder**

This is the user's hands-on validation (raw + DB path). Present this checklist to the user (do not automate):
1. In the Process tab, select a source results file (imported and/or raw-linked, constrained).
2. Run **Sources → "Mode-frequency spectrum: PSD (Welch)"** with default 300 modes.
3. Confirm a `timefreq_psd` node appears; double-click → power spectrum opens with one line per mode, x-axis = frequency.
4. Repeat with **"Mode-frequency spectrum: FFT"**.
5. On a raw file, confirm it does not load the whole recording when a time window is set, and that a marked bad segment changes the PSD (consistency with a standard source PSD).

- [ ] **Step 4: Commit (if any cleanup was needed)**

```bash
git add -A
git commit -m "Mode-frequency: test sweep + manual validation checklist" || echo "nothing to commit"
```

> **Push:** Only when the user explicitly asks. Push **only** to `origin` (the fork `git@github.com:DiellorBasha/brainstorm3.git`), never to `upstream`.

---

## Self-Review Notes (author)

- **Spec coverage:** mode kernel (Task 1); engine `Compute` + read + save + bad-segment parity (Task 2); PSD process (Task 3); FFT process (Task 4); validation (Task 5). The joint heatmap, Morlet/Hilbert, unconstrained, and full-raw streaming are explicitly out of scope per the spec.
- **Type consistency:** `Compute(F, sfreq, ModeKernel, BadSegments, Method, WinLength, WinOverlap, WinFunc, PowerUnits, Values)` returns `[TF, Freqs, RowNames, Msg]` — used identically in the test and in `Run`. `bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)` signature matches all call sites. Option names (`nmodes`, `win_length`, `win_overlap`, `units`, `win_std`, `timewindow`) match between `GetDescription` and `Run`.
- **No placeholders:** all steps contain full code and exact run targets.
