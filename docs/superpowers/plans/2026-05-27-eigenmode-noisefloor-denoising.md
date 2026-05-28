# Eigenmode Noise-Floor Denoising Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute an empirical joint (λ,ω) noise floor from a paired empty-room recording and use it to produce an SNR-resolved spectrum + reliable-mode cutoff and a power-spectral-subtraction cleaned spectrum.

**Architecture:** A pure combine function (`bst_eigenmodes_noisefloor`) turns two PSDs (data, noise) into SNR / cleaned-PSD / Wiener-gain / K*(f). A two-input process (`process_eigenmodes_denoise`, Files A = data, Files B = empty-room) builds the data's transform kernel, applies it to *both* recordings on their common good channels, Welch-PSDs both, calls the pure function, and saves SNR + cleaned-spectrum timefreq files. Builds on M1's `bst_eigenmodes_transform`.

**Tech Stack:** MATLAB, Brainstorm process system (`eval(macro_method)`, `nInputs=2`), `bst_eigenmodes_transform`, `bst_psd` (Welch), `in_bst_data`/`in_bst_channel`/`good_channel`, `db_template('timefreqmat')`, `bst_save`/`db_add_data`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_eigenmodes_noisefloor.m` (create) | Pure: `(Pdata, Nnoise) → SNR, CleanPSD, Gain, Kstar`. No I/O, no Welch. |
| `toolbox/process/functions/process_eigenmodes_denoise.m` (create) | Process (nInputs=2): kernel → apply to data+noise on common channels → Welch PSDs → pure fn → save SNR + cleaned-spectrum timefreq. |
| `dev/tests/test_eigenmodes_noisefloor_pure.m` (create) | DB-free unit test of the pure function. |
| `dev/tests/test_process_eigenmodes_denoise_options.m` (create) | DB-free unit test of the process options/shape. |

**Run convention (tests):** execute the function in MATLAB with the repo on path via the MATLAB MCP `evaluate_matlab_code` (call the function name) or `run_matlab_file` — NOT `run_matlab_test_file` (these are script-style tests printing `ALL TESTS PASSED`, not `matlab.unittest`). Lint via `check_matlab_code`.

The disposable `EigenSmoke` protocol (subject `SmokeS`, icosphere cortex, head model, eigenmodes, a 60 s imported data block) from the M1 smoke already exists and is reused for live validation in Task 2.

---

## Task 1: Pure combine function `bst_eigenmodes_noisefloor`

**Files:**
- Create: `toolbox/math/bst_eigenmodes_noisefloor.m`
- Test: `dev/tests/test_eigenmodes_noisefloor_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_noisefloor_pure.m`:

```matlab
function test_eigenmodes_noisefloor_pure
% Verify SNR / power-subtraction / Wiener-gain / reliable-mode-cutoff math.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(3);
K = 30; nF = 40;
S = zeros(K, nF); S(1:10,:) = rand(10, nF) + 0.5;   % signal only in low modes
N = 0.2*rand(K, nF) + 0.1;                            % positive noise floor
Pdata = S + N;

Out = bst_eigenmodes_noisefloor(Pdata, N);            % defaults: Alpha=1, Floor=0, SnrThresh=1

% Power subtraction recovers the signal power exactly (Alpha=1, Floor=0).
assert(max(abs(Out.CleanPSD(:) - S(:))) < 1e-12, 'CleanPSD should equal S.');
% SNR = (S+N)/N.
SNRexp = Pdata ./ N;
assert(max(abs(Out.SNR(:) - SNRexp(:))) < 1e-12, 'SNR mismatch.');
% Wiener gain in [0,1], and exactly 0 where Pdata <= N (modes 11..30 have S=0 -> equal).
assert(all(Out.Gain(:) >= 0 & Out.Gain(:) <= 1), 'Gain must be in [0,1].');
below = (Pdata <= N);
assert(all(Out.Gain(below) == 0), 'Gain must be 0 where Pdata<=N.');

% Floor clamp respected with over-subtraction.
Out2 = bst_eigenmodes_noisefloor(Pdata, N, 'Alpha', 2, 'Floor', 0.1);
assert(all(Out2.CleanPSD(:) >= 0.1*N(:) - 1e-12), 'CleanPSD must respect the spectral floor.');

% Reliable-mode cutoff on a descending-SNR ramp at frequency 1.
Pramp = N;
Pramp(:,1) = N(:,1) .* (3 - (0:K-1)'*0.1);            % SNR(k,1) = 3 - 0.1*(k-1)
Outr = bst_eigenmodes_noisefloor(Pramp, N, 'SnrThresh', 1);
% SNR>=1  <=>  3 - 0.1*(k-1) >= 1  <=>  k <= 21
assert(Outr.Kstar(1) == 21, 'Kstar ramp wrong (got %d, expected 21).', Outr.Kstar(1));

fprintf('ALL TESTS PASSED: test_eigenmodes_noisefloor_pure\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`): `test_eigenmodes_noisefloor_pure`
Expected: FAIL — `Undefined function or variable 'bst_eigenmodes_noisefloor'`.

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/math/bst_eigenmodes_noisefloor.m`:

```matlab
function Out = bst_eigenmodes_noisefloor(Pdata, Nnoise, varargin)
% BST_EIGENMODES_NOISEFLOOR: Joint (lambda,omega) SNR + spectral-subtraction denoising.
%
% USAGE:  Out = bst_eigenmodes_noisefloor(Pdata, Nnoise)
%         Out = bst_eigenmodes_noisefloor(Pdata, Nnoise, 'Alpha',1, 'Floor',0, 'SnrThresh',1)
%
% DESCRIPTION:
%     Combines a data power spectrum and an empty-room (noise) power spectrum,
%     both in the eigenmode x frequency plane, into denoising products. Works on
%     POWER (averaged PSD), never on complex coefficients: the data and noise
%     recordings are different noise realizations, so complex subtraction would
%     add variance, whereas E[|data|^2] - E[|noise|^2] = |signal|^2.
%
%     Products:
%       SNR(k,f)      = Pdata / Nnoise
%       CleanPSD(k,f) = max(Pdata - Alpha*Nnoise, Floor*Nnoise)   (spectral subtraction)
%       Gain(k,f)     = max(Pdata - Nnoise, 0) / Pdata            (Wiener gain in [0,1])
%       Kstar(f)      = largest mode index k with SNR(k,f) >= SnrThresh (0 if none)
%
% INPUTS:
%     Pdata  : [K x nFreq] data PSD (power/Hz) per eigenmode per frequency.
%     Nnoise : [K x nFreq] empty-room PSD, same size and frequency grid.
%
% OPTIONS (name-value):
%     'Alpha'     : over-subtraction factor (>=1), default 1.
%     'Floor'     : spectral floor as a fraction of Nnoise, default 0.
%     'SnrThresh' : linear SNR threshold for the reliable-mode cutoff, default 1.
%
% OUTPUTS:
%     Out.SNR, Out.CleanPSD, Out.Gain : [K x nFreq]
%     Out.Kstar                       : [1 x nFreq]
%
% SEE ALSO: bst_eigenmodes_transform, process_eigenmodes_denoise

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
Alpha = 1; Floor = 0; SnrThresh = 1;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'alpha',     Alpha     = varargin{i+1};
        case 'floor',     Floor     = varargin{i+1};
        case 'snrthresh', SnrThresh = varargin{i+1};
    end
end

Pdata  = double(Pdata);
Nnoise = double(Nnoise);
if ~isequal(size(Pdata), size(Nnoise))
    error('bst_eigenmodes_noisefloor: Pdata and Nnoise must have the same size.');
end

%% ===== COMBINE =====
Nsafe = max(Nnoise, eps);
Out = struct();
Out.SNR      = Pdata ./ Nsafe;
Out.CleanPSD = max(Pdata - Alpha .* Nnoise, Floor .* Nnoise);
Out.Gain     = max(Pdata - Nnoise, 0) ./ max(Pdata, eps);

%% ===== RELIABLE-MODE CUTOFF =====
[K, nFreq] = size(Pdata);
Out.Kstar = zeros(1, nFreq);
for f = 1:nFreq
    idx = find(Out.SNR(:, f) >= SnrThresh, 1, 'last');
    if ~isempty(idx)
        Out.Kstar(f) = idx;
    end
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `test_eigenmodes_noisefloor_pure`
Expected: PASS — prints `ALL TESTS PASSED: test_eigenmodes_noisefloor_pure`.

- [ ] **Step 5: Lint**

`check_matlab_code` on `toolbox/math/bst_eigenmodes_noisefloor.m` — no genuine errors.

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_eigenmodes_noisefloor.m dev/tests/test_eigenmodes_noisefloor_pure.m
git commit -m "$(cat <<'EOF'
Add bst_eigenmodes_noisefloor: SNR + spectral-subtraction combine math

Pure: from data PSD + empty-room PSD compute SNR(k,f), power-subtraction
CleanPSD=max(P-aN, floor*N), Wiener gain max(P-N,0)/P in [0,1], and the
reliable-mode cutoff K*(f). Power-domain (not complex). Covered by
test_eigenmodes_noisefloor_pure.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Process `process_eigenmodes_denoise`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_denoise.m`
- Test: `dev/tests/test_process_eigenmodes_denoise_options.m`

- [ ] **Step 1: Write the failing options test**

Create `dev/tests/test_process_eigenmodes_denoise_options.m`:

```matlab
function test_process_eigenmodes_denoise_options
% Verify the denoise process is a 2-input Sources process with the right knobs.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_denoise('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(sProcess.nInputs == 2, 'Must be a 2-input process (Files A = data, Files B = empty-room).');
assert(abs(sProcess.Index - 336.6) < 1e-9, 'Index must be 336.6.');
for f = {'nmodes','noisewin','alpha','snrthresh','floorfrac'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
fprintf('ALL TESTS PASSED: test_process_eigenmodes_denoise_options\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `test_process_eigenmodes_denoise_options`
Expected: FAIL — `Undefined function or variable 'process_eigenmodes_denoise'`.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_denoise.m`:

```matlab
function varargout = process_eigenmodes_denoise( varargin )
% PROCESS_EIGENMODES_DENOISE: Joint (lambda,omega) noise floor from empty-room recordings.
%
% USAGE:  sProcess = process_eigenmodes_denoise('GetDescription')
%       OutputFiles = process_eigenmodes_denoise('Run', sProcess, sInputsA, sInputsB)
%
% DESCRIPTION:
%     Files A = data recording(s) (must have a surface head model + eigenmodes).
%     Files B = empty-room recording(s). Builds the data's eigenmode transform
%     kernel A = pinv(L*Phi), applies it to BOTH recordings on their common good
%     channels, Welch-PSDs both, then computes SNR(k,f), a power-spectral-
%     subtraction cleaned spectrum, and a reliable-mode cutoff K*(f).
%     Subtraction is on POWER (averaged PSD), never complex coefficients.
%     Both recordings must be imported (not raw). Requires precomputed eigenmodes.
%
% SEE ALSO: bst_eigenmodes_noisefloor, bst_eigenmodes_transform, in_tess_eigenmodes

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
    sProcess.Comment     = 'Eigenmode noise-floor denoising';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.6;   % just after the eigenmode transform (336.5)
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 2;       % Files A = data, Files B = empty-room
    sProcess.nMinFiles   = 1;

    sProcess.options.nmodes.Comment    = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type       = 'value';
    sProcess.options.nmodes.Value      = {0, '', 0};

    sProcess.options.noisewin.Comment  = 'Welch window length: ';
    sProcess.options.noisewin.Type     = 'value';
    sProcess.options.noisewin.Value    = {2, 's', 2};

    sProcess.options.alpha.Comment     = 'Over-subtraction factor alpha (>=1): ';
    sProcess.options.alpha.Type        = 'value';
    sProcess.options.alpha.Value       = {1, '', 2};

    sProcess.options.snrthresh.Comment = 'Reliable-mode SNR threshold (linear): ';
    sProcess.options.snrthresh.Type    = 'value';
    sProcess.options.snrthresh.Value   = {1, '', 2};

    sProcess.options.floorfrac.Comment = 'Spectral floor (fraction of noise): ';
    sProcess.options.floorfrac.Type    = 'value';
    sProcess.options.floorfrac.Value   = {0, '', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Files A = data, Files B = empty-room. ' ...
        'Subtraction is on Welch-averaged power (not complex coefficients).<BR>' ...
        'Outputs an SNR(\lambda,\omega) spectrum + cleaned power spectrum. Both inputs must be imported.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    a = sProcess.options.alpha.Value{1};
    Comment = sprintf('Eigenmode noise-floor denoising (alpha=%.1f)', a);
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputsA, sInputsB) %#ok<DEFNU>
    OutputFiles = {};
    if isempty(sInputsB)
        bst_report('Error', sProcess, sInputsA, 'Select the empty-room recording(s) as Files B.');
        return;
    end
    nModesOpt = sProcess.options.nmodes.Value{1};
    WinLen    = sProcess.options.noisewin.Value{1};
    Alpha     = sProcess.options.alpha.Value{1};
    SnrThresh = sProcess.options.snrthresh.Value{1};
    FloorFrac = sProcess.options.floorfrac.Value{1};

    % ===== DATA STUDY: head model, surface, eigenmodes, gain =====
    [sStudyA,~,~,~] = bst_get('Study', sInputsA(1).iStudy);
    if isempty(sStudyA.iHeadModel) || sStudyA.iHeadModel < 1
        bst_report('Error', sProcess, sInputsA, 'No head model for the data study (Files A).');
        return;
    end
    HeadModelFile = sStudyA.HeadModel(sStudyA.iHeadModel).FileName;
    HMmeta = in_bst_headmodel(HeadModelFile, 0, 'HeadModelType', 'SurfaceFile');
    if ~strcmpi(HMmeta.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputsA, 'Eigenmode denoise requires a surface head model.');
        return;
    end
    SurfaceFile = HMmeta.SurfaceFile;
    [Em, isC] = in_tess_eigenmodes(SurfaceFile);
    if ~isC
        bst_report('Error', sProcess, sInputsA, ['No eigenmodes on surface: ' SurfaceFile '. Run "Compute eigenmodes" first.']);
        return;
    end
    HM = in_bst_headmodel(HeadModelFile, 1);   % ApplyOrient=1 -> [nch x nVert]
    Gain = double(HM.Gain);
    if size(Gain,2) ~= size(Em.Vectors,1)
        bst_report('Error', sProcess, sInputsA, sprintf('Head model (%d) / eigenmode (%d) vertex mismatch.', size(Gain,2), size(Em.Vectors,1)));
        return;
    end

    % ===== CHANNELS: common good channels (by name) between data and noise =====
    ChA = in_bst_channel(bst_get('ChannelFileForStudy', sStudyA.FileName));
    DA  = in_bst_data(sInputsA(1).FileName);
    if isstruct(DA.F)
        bst_report('Error', sProcess, sInputsA, 'Files A must be imported data (not raw). Import a block first.');
        return;
    end
    iA = good_channel(ChA.Channel, DA.ChannelFlag, 'MEG');
    if isempty(iA), iA = good_channel(ChA.Channel, DA.ChannelFlag, 'EEG'); end

    [sStudyB,~,~,~] = bst_get('Study', sInputsB(1).iStudy);
    ChB = in_bst_channel(bst_get('ChannelFileForStudy', sStudyB.FileName));
    DB  = in_bst_data(sInputsB(1).FileName);
    if isstruct(DB.F)
        bst_report('Error', sProcess, sInputsB, 'Files B (empty-room) must be imported data (not raw). Import a block first.');
        return;
    end
    iB = good_channel(ChB.Channel, DB.ChannelFlag, 'MEG');
    if isempty(iB), iB = good_channel(ChB.Channel, DB.ChannelFlag, 'EEG'); end

    [~, ia, ib] = intersect({ChA.Channel(iA).Name}, {ChB.Channel(iB).Name}, 'stable');
    if isempty(ia)
        bst_report('Error', sProcess, sInputsA, 'No common good channels between data and empty-room.');
        return;
    end
    iCommonA = iA(ia);
    iCommonB = iB(ib);

    % ===== KERNEL ON COMMON CHANNELS =====
    nCh = numel(iCommonA);
    if isempty(nModesOpt) || nModesOpt <= 0
        K = min(nCh, Em.nModes);
    else
        K = min(nModesOpt, Em.nModes);
    end
    Phi     = double(Em.Vectors(:, 1:K));
    lambdas = double(Em.Values(1:K));
    [A, Info] = bst_eigenmodes_transform(Gain(iCommonA, :), Phi);   % [K x nCh]

    % ===== COEFFICIENTS + WELCH PSDs (density, same window) =====
    thD = A * double(DA.F(iCommonA, :));
    thN = A * double(DB.F(iCommonB, :));
    sfA = 1 / (DA.Time(2) - DA.Time(1));
    sfB = 1 / (DB.Time(2) - DB.Time(1));
    [TFd, Fv]  = bst_psd(thD, sfA, WinLen, 50, [], [], [], 'physical');
    [TFn, Fvn] = bst_psd(thN, sfB, WinLen, 50, [], [], [], 'physical');
    if numel(Fv) ~= numel(Fvn) || max(abs(Fv(:) - Fvn(:))) > 1e-6
        bst_report('Error', sProcess, sInputsA, 'Data and empty-room PSD frequency grids differ (different sampling rate or window).');
        return;
    end
    Pdata  = reshape(TFd, K, []);
    Nnoise = reshape(TFn, K, []);

    % ===== COMBINE =====
    Out = bst_eigenmodes_noisefloor(Pdata, Nnoise, 'Alpha', Alpha, 'Floor', FloorFrac, 'SnrThresh', SnrThresh);

    % ===== SAVE =====
    RowNames = cell(K,1);
    for k = 1:K, RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, lambdas(k)); end
    [sStudyOut, iStudyOut] = bst_get('Study', sInputsA(1).iStudy);
    StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

    OutputFiles{end+1} = SaveTF(Out.SNR, Fv, RowNames, DA.Time, StudyDir, iStudyOut, ...
        sprintf('EigenSNR (%d modes) | %s', K, sInputsA(1).Comment), 'timefreq_eigensnr', SurfaceFile); %#ok<AGROW>
    OutputFiles{end+1} = SaveTF(Out.CleanPSD, Fv, RowNames, DA.Time, StudyDir, iStudyOut, ...
        sprintf('EigenCleanPSD (%d modes, a=%.1f) | %s', K, Alpha, sInputsA(1).Comment), 'timefreq_eigencleanpsd', SurfaceFile); %#ok<AGROW>

    bst_report('Info', sProcess, sInputsA, sprintf('Denoise: %d modes (condition %.1f), median reliable-mode cutoff K*=%d at SNR>=%.1f.', ...
        K, Info.ConditionNumber, round(median(Out.Kstar)), SnrThresh));
end


%% ===== SAVE ONE (mode x freq) MATRIX AS A TIMEFREQ FILE =====
function OutFile = SaveTF(M2d, Freqs, RowNames, Time, StudyDir, iStudyOut, Comment, prefix, SurfaceFile)
    [K, nFreq] = size(M2d);
    TFmat = db_template('timefreqmat');
    TFmat.TF          = reshape(M2d, [K, 1, nFreq]);
    TFmat.Freqs       = Freqs(:)';
    TFmat.Time        = [Time(1), Time(end)];
    TFmat.RowNames    = RowNames;
    TFmat.Measure     = 'power';
    TFmat.Method      = 'psd';
    TFmat.DataType    = 'matrix';
    TFmat.SurfaceFile = SurfaceFile;
    TFmat.Comment     = Comment;
    TFmat.nAvg        = 1;
    TFmat = bst_history('add', TFmat, 'eigenmodes_denoise', Comment);
    FullFile = bst_process('GetNewFilename', StudyDir, prefix);
    bst_save(FullFile, TFmat, 'v6');
    db_add_data(iStudyOut, FullFile, TFmat);
    OutFile = file_short(FullFile);
end
```

- [ ] **Step 4: Run the options test to verify it passes**

Run: `test_process_eigenmodes_denoise_options`
Expected: PASS — prints `ALL TESTS PASSED: test_process_eigenmodes_denoise_options`.

- [ ] **Step 5: Lint**

`check_matlab_code` on `toolbox/process/functions/process_eigenmodes_denoise.m` — only standard Brainstorm idioms (`varargout`, stale `%#ok`) acceptable.

- [ ] **Step 6: Live end-to-end validation on the EigenSmoke substrate (verifies the Run path)**

The options test does not exercise `Run`; this step does. In MATLAB (MCP `evaluate_matlab_code`):

```matlab
% Import a 60 s empty-room block into SmokeS (subjectname REQUIRED — defaults to NewSubject otherwise)
bidsDir   = '/Users/diellorbasha/workspace/library/datasets/omega-tutorial';
noiseFile = fullfile(bidsDir,'sub-emptyroom','ses-18901014','meg','sub-emptyroom_ses-18901014_task-noise_run-01_meg.ds');
sNraw = bst_process('CallProcess','process_import_data_raw', [], [], ...
    'subjectname','SmokeS', 'datafile',{noiseFile,'CTF'}, 'channelreplace',0, 'channelalign',0, 'evtmode','value');
sNimp = bst_process('CallProcess','process_import_data_time', sNraw, [], ...
    'subjectname','SmokeS', 'timewindow',[0 60], 'split',0, 'usectfcomp',1, 'usessp',1);
% Re-select the existing data block (the M1 smoke imported one into SmokeS)
sData = bst_process('CallProcess','process_select_files_data', [], [], 'subjectname','SmokeS');
% pick the imported data block (not a raw link): the one whose study has a head model
% (the M1 smoke condition 'sub-0002_ses-01_task-rest_run-01_meg'); use the first imported block:
iData = find(arrayfun(@(s) isempty(strfind(s.FileName,'@raw')) && ~isempty(strfind(s.FileName,'data_block')), sData), 1);
assert(~isempty(iData), 'No imported data block found in SmokeS.');
sDen = bst_process('CallProcess','process_eigenmodes_denoise', sData(iData), sNimp, ...
    'nmodes',0, 'noisewin',2, 'alpha',1, 'snrthresh',1, 'floorfrac',0);
assert(numel(sDen) >= 2, 'Expected 2 timefreq outputs (SNR + cleaned).');
T = load(file_fullpath(sDen(1).FileName));
fprintf('SNR file TF size=%s Freqs=[%.2f %.1f] rows=%d\n', mat2str(size(T.TF)), T.Freqs(1), T.Freqs(end), numel(T.RowNames));
assert(all(isfinite(T.TF(:))) && ndims(T.TF)==3, 'SNR TF must be finite [K x 1 x nFreq].');
fprintf('LIVE VALIDATION OK\n');
```

Expected: prints `LIVE VALIDATION OK` with a finite `[K x 1 x nFreq]` SNR tensor. If any Brainstorm-glue detail is off (timefreq fields, `bst_psd` shape, channel matching), fix `Run`/`SaveTF` until this passes, then re-run the options test.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_denoise.m dev/tests/test_process_eigenmodes_denoise_options.m
git commit -m "$(cat <<'EOF'
Add process_eigenmodes_denoise: empty-room (lambda,omega) noise floor

Two-input process (Files A=data, Files B=empty-room): builds the data
transform kernel, applies it to both recordings on common good channels,
Welch-PSDs both, and saves an SNR(lambda,omega) spectrum + spectral-
subtraction cleaned spectrum via bst_eigenmodes_noisefloor. Validated
end-to-end on the OMEGA EigenSmoke substrate.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: (Optional, manual) SNR figure

Not required to close the milestone. After Task 2, in MATLAB, load `sDen(1)` (SNR) and `imagesc(Freqs, 1:K, 10*log10(squeeze(T.TF)))` over 0–80 Hz; confirm: the 60 Hz line shows **low** SNR (present in both data and floor → correctly not flagged as signal), the alpha band shows SNR>1 at low modes, and `Kstar` declines with frequency. Save to `dev/tests/eigentransform_denoise_snr.png`.

---

## Self-Review

**Spec coverage:**
- Empirical noise floor via same kernel on empty-room → Task 2 Run (kernel built from data study, applied to both). ✓
- Common-channel handling / rebuild A on common set → Task 2 (`intersect` by name, `bst_eigenmodes_transform(Gain(iCommonA,:),Phi)`). ✓
- Welch-averaged PSDs, density units → Task 2 (`bst_psd(...,'physical')`). ✓
- Power (not complex) subtraction + clamp; SNR; Wiener gain; K*(f) → Task 1 pure fn + test. ✓
- Diagnostic (SNR spectrum) + cleaned-spectrum outputs as timefreq → Task 2 `SaveTF`. ✓
- Reliable-mode cutoff threshold default SNR≥1 → Task 1 `SnrThresh` default 1; option in Task 2. ✓
- Wiener applied to produce cleaned *time series*: **deferred** (the spec marked it optional/off-by-default). The pure fn computes & tests `Gain`; wiring it to coefficients is left to a follow-on. Noted here so spec/plan agree.
- Index 336.6 (verified free) → Task 2 + options test. ✓
- Builds on M1 `bst_eigenmodes_transform` → reused in Task 2. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete code; run steps have exact commands + expected output.

**Type/name consistency:** `bst_eigenmodes_noisefloor(Pdata,Nnoise,...)` returns `Out.SNR/.CleanPSD/.Gain/.Kstar` — used identically in Task 1 test and Task 2 Run. `bst_eigenmodes_transform(Gain,Phi) -> [A, Info]` with `Info.ConditionNumber` matches M1. File prefixes `timefreq_eigensnr` / `timefreq_eigencleanpsd` consistent between Run and commit message.
