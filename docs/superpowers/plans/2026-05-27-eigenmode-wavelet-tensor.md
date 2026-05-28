# Eigenmode Complex Wavelet Tensor Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute and store the complex Morlet wavelet tensor `Wₖ(s,t)` of the eigenmode coefficients — the time-resolved `(λ,ω,t)` decomposition carrying amplitude and phase.

**Architecture:** A pure function (`bst_eigenmodes_wavelet`) wraps Brainstorm's validated `morlet_transform` (complex mode) and returns a `[K × nTime × nFreq]` complex tensor; a thin process (`process_eigenmodes_wavelet`) applies it to an eigenmode-coefficient `matrix` file and saves a complex, λ-labeled timefreq viewable in Brainstorm's TF viewer. All analyses (dispersion, phase coherence, etc.) are deferred.

**Tech Stack:** MATLAB, Brainstorm process system (`eval(macro_method)`), `morlet_transform`, `in_bst_matrix`, `db_template('timefreqmat')`, `bst_save`/`db_add_data`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_eigenmodes_wavelet.m` (create) | Pure: `(Coeffs, sfreq, Freqs) → complex W [K×nTime×nFreq]` via `morlet_transform('n')` + permute. No I/O. |
| `toolbox/process/functions/process_eigenmodes_wavelet.m` (create) | Process: load coefficient matrix → build freq grid → pure fn → save complex λ-labeled timefreq. |
| `dev/tests/test_eigenmodes_wavelet_pure.m` (create) | DB-free unit test of the pure function. |
| `dev/tests/test_process_eigenmodes_wavelet_options.m` (create) | DB-free unit test of the process options/shape. |

**Run convention (tests):** MATLAB MCP `evaluate_matlab_code` calling the function name (script-style tests print `ALL TESTS PASSED`; NOT `run_matlab_test_file`). Lint via `check_matlab_code`.

Reference facts (already confirmed): `morlet_transform(x, t, f, fc, FWHM_tc, squared)` returns `[nSignals × nFreq × nTime]`; `squared='n'` ⇒ complex coefficients. `in_bst_matrix(file)` returns `.Value [K×nTime]`, `.Time`, `.Description`, `.SurfaceFile`. The disposable `EigenSmoke` protocol (subject `SmokeS`) already contains a `matrix_eigentransform` coefficient file from the M1 smoke.

---

## Task 1: Pure wavelet function `bst_eigenmodes_wavelet`

**Files:**
- Create: `toolbox/math/bst_eigenmodes_wavelet.m`
- Test: `dev/tests/test_eigenmodes_wavelet_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_wavelet_pure.m`:

```matlab
function test_eigenmodes_wavelet_pure
% Verify the complex Morlet wavelet tensor: shape, complexity, amplitude
% localization, phase rate, and the default frequency grid.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

rng(5);
sfreq = 200; T = 2; nT = sfreq*T; t = (0:nT-1)/sfreq;
f0 = 10; K = 8; kSig = 3;
Coeffs = zeros(K, nT);
Coeffs(kSig,:) = cos(2*pi*f0*t);          % one mode carries a 10 Hz sinusoid
Freqs = 2:2:40;                            % grid includes 10 Hz

[W, Fout] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs);

assert(isequal(size(W), [K, nT, numel(Freqs)]), 'W must be [K x nTime x nFreq].');
assert(~isreal(W), 'W must be complex.');
assert(isequal(Fout(:)', Freqs), 'Freqs must pass through unchanged.');

% Amplitude localizes to (kSig, f0).
A = squeeze(mean(abs(W), 2));              % [K x nFreq], time-averaged amplitude
[~, imax] = max(A(:)); [kmax, fmax] = ind2sub(size(A), imax);
assert(kmax == kSig, 'Peak mode wrong (got %d, expected %d).', kmax, kSig);
[~, if0] = min(abs(Freqs - f0));
assert(fmax == if0, 'Peak frequency wrong (got bin %d).', fmax);
others = setdiff(1:K, kSig);
assert(max(A(others, if0)) < 0.1 * A(kSig, if0), 'Signal leaked to other modes.');

% Phase advances at ~2*pi*f0/sfreq per sample in the central (edge-free) region.
wsig = squeeze(W(kSig, :, if0));
ph = unwrap(angle(wsig(:)'));
mid = round(nT*0.4):round(nT*0.6);
dph = mean(diff(ph(mid)));
assert(abs(dph - 2*pi*f0/sfreq) < 0.2*(2*pi*f0/sfreq), 'Phase rate wrong (got %.4f).', dph);

% Empty Freqs -> default 40-frequency log grid.
[W2, F2] = bst_eigenmodes_wavelet(Coeffs, sfreq, []);
assert(numel(F2) == 40, 'Default grid must have 40 frequencies.');
assert(isequal(size(W2), [K, nT, 40]), 'Default tensor must be [K x nTime x 40].');

fprintf('ALL TESTS PASSED: test_eigenmodes_wavelet_pure\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`): `test_eigenmodes_wavelet_pure`
Expected: FAIL — `Undefined function or variable 'bst_eigenmodes_wavelet'`.

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/math/bst_eigenmodes_wavelet.m`:

```matlab
function [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, varargin)
% BST_EIGENMODES_WAVELET: Complex Morlet wavelet tensor of eigenmode coefficients.
%
% USAGE:  [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs)
%         [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, [], 'MorletFc',1, 'MorletFwhmTc',3)
%
% DESCRIPTION:
%     Applies the complex Morlet continuous wavelet transform to each eigenmode
%     coefficient time series, producing the time-resolved (lambda, omega, t)
%     tensor W_k(s,t). The result is COMPLEX: |W| is amplitude, arg(W) is phase.
%     This wraps Brainstorm's validated morlet_transform (squared='n', i.e. the
%     un-squared complex coefficients) and permutes its [K x nFreq x nTime]
%     output to the Brainstorm timefreq convention [K x nTime x nFreq].
%
% INPUTS:
%     Coeffs : [K x nTime] eigenmode coefficient time series.
%     sfreq  : sampling rate (Hz).
%     Freqs  : [1 x nFreq] frequencies (Hz). If empty, a default log-spaced grid
%              logspace(log10(2), log10(min(100, 0.4*sfreq)), 40) is used.
%
% OPTIONS (name-value):
%     'MorletFc'     : Morlet central frequency (default 1).
%     'MorletFwhmTc' : Morlet time resolution / FWHM in seconds (default 3).
%
% OUTPUTS:
%     W     : [K x nTime x nFreq] complex wavelet tensor.
%     Freqs : [1 x nFreq] frequencies actually used.
%
% SEE ALSO: morlet_transform, bst_eigenmodes_transform, process_eigenmodes_wavelet

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
MorletFc = 1; MorletFwhmTc = 3;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'morletfc',     MorletFc     = varargin{i+1};
        case 'morletfwhmtc', MorletFwhmTc = varargin{i+1};
    end
end

Coeffs = double(Coeffs);
[~, nTime] = size(Coeffs);

%% ===== FREQUENCY GRID =====
if nargin < 3 || isempty(Freqs)
    Freqs = logspace(log10(2), log10(min(100, 0.4*sfreq)), 40);
end
Freqs = Freqs(:)';

%% ===== MORLET CWT (complex) =====
t = (0:nTime-1) / sfreq;
% morlet_transform returns [K x nFreq x nTime]; 'n' keeps complex coefficients.
P = morlet_transform(Coeffs, t, Freqs, MorletFc, MorletFwhmTc, 'n');
% Permute to the Brainstorm timefreq convention [K x nTime x nFreq].
W = permute(P, [1 3 2]);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `test_eigenmodes_wavelet_pure`
Expected: PASS — `ALL TESTS PASSED: test_eigenmodes_wavelet_pure`.

- [ ] **Step 5: Lint**

`check_matlab_code` on `toolbox/math/bst_eigenmodes_wavelet.m` — no genuine errors.

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_eigenmodes_wavelet.m dev/tests/test_eigenmodes_wavelet_pure.m
git commit -m "$(cat <<'EOF'
Add bst_eigenmodes_wavelet: complex Morlet wavelet tensor of coefficients

Pure: wraps morlet_transform ('n' = complex), permutes [K x nFreq x nTime]
-> [K x nTime x nFreq], default log freq grid. Returns the complex
(lambda,omega,t) tensor W_k(s,t). Covered by test_eigenmodes_wavelet_pure.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Process `process_eigenmodes_wavelet`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_wavelet.m`
- Test: `dev/tests/test_process_eigenmodes_wavelet_options.m`

- [ ] **Step 1: Write the failing options test**

Create `dev/tests/test_process_eigenmodes_wavelet_options.m`:

```matlab
function test_process_eigenmodes_wavelet_options
% Verify the wavelet process is a matrix-in / timefreq-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_wavelet('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.7) < 1e-9, 'Index must be 336.7.');
assert(isequal(sProcess.InputTypes, {'matrix'}), 'InputTypes must be {matrix}.');
assert(isequal(sProcess.OutputTypes, {'timefreq'}), 'OutputTypes must be {timefreq}.');
for f = {'flo','fhi','nfreqs','morletfc','morletfwhmtc'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
fprintf('ALL TESTS PASSED: test_process_eigenmodes_wavelet_options\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `test_process_eigenmodes_wavelet_options`
Expected: FAIL — `Undefined function or variable 'process_eigenmodes_wavelet'`.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_wavelet.m`:

```matlab
function varargout = process_eigenmodes_wavelet( varargin )
% PROCESS_EIGENMODES_WAVELET: Complex Morlet wavelet tensor of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_wavelet('GetDescription')
%       OutputFiles = process_eigenmodes_wavelet('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Applies the complex Morlet CWT to an eigenmode-coefficient matrix
%     (matrix_eigentransform, [K x nTime]) and saves the time-resolved
%     (lambda, omega, t) tensor W_k(s,t) as a complex, lambda-labeled timefreq
%     file (amplitude + phase), viewable in Brainstorm's time-frequency viewer.
%
% SEE ALSO: bst_eigenmodes_wavelet, process_eigenmodes_transform

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
    sProcess.Comment     = 'Eigenmode wavelet tensor (complex)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.7;   % after the noise-floor denoise (336.6)
    sProcess.Description = '';
    sProcess.InputTypes  = {'matrix'};
    sProcess.OutputTypes = {'timefreq'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.flo.Comment    = 'Lowest frequency: ';
    sProcess.options.flo.Type       = 'value';
    sProcess.options.flo.Value      = {2, 'Hz', 2};

    sProcess.options.fhi.Comment    = 'Highest frequency (0 = auto, min(100, 0.4*sfreq)): ';
    sProcess.options.fhi.Type       = 'value';
    sProcess.options.fhi.Value      = {0, 'Hz', 2};

    sProcess.options.nfreqs.Comment = 'Number of frequencies (log-spaced): ';
    sProcess.options.nfreqs.Type    = 'value';
    sProcess.options.nfreqs.Value   = {40, '', 0};

    sProcess.options.morletfc.Comment     = 'Morlet central frequency Fc: ';
    sProcess.options.morletfc.Type        = 'value';
    sProcess.options.morletfc.Value       = {1, 'Hz', 2};

    sProcess.options.morletfwhmtc.Comment = 'Morlet time resolution (FWHM): ';
    sProcess.options.morletfwhmtc.Type    = 'value';
    sProcess.options.morletfwhmtc.Value   = {3, 's', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Complex Morlet wavelet tensor W_k(s,t) of the ' ...
        'eigenmode coefficients &mdash; the time-resolved (&lambda;,&omega;,t) decomposition (amplitude + phase).<BR>' ...
        'Input: an eigenmode-coefficient matrix (matrix_eigentransform).</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = 'Eigenmode wavelet tensor (complex)';
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    flo  = sProcess.options.flo.Value{1};
    fhi  = sProcess.options.fhi.Value{1};
    nfq  = sProcess.options.nfreqs.Value{1};
    fc   = sProcess.options.morletfc.Value{1};
    fwhm = sProcess.options.morletfwhmtc.Value{1};

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        M = in_bst_matrix(sInput.FileName);
        Coeffs = double(M.Value);                       % [K x nTime]
        if size(Coeffs, 2) < 2
            bst_report('Error', sProcess, sInput, 'Need at least 2 time samples to estimate a wavelet transform.');
            continue;
        end
        sfreq = 1 / (M.Time(2) - M.Time(1));
        fhiEff = fhi;
        if fhiEff <= 0
            fhiEff = min(100, 0.4 * sfreq);
        end
        Freqs = logspace(log10(flo), log10(fhiEff), nfq);

        [W, Freqs] = bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, 'MorletFc', fc, 'MorletFwhmTc', fwhm);

        % Row labels: reuse the matrix Description (mode labels) or build them.
        RowNames = M.Description;
        if isempty(RowNames) || numel(RowNames) ~= size(W,1)
            RowNames = cell(size(W,1), 1);
            for k = 1:size(W,1), RowNames{k} = sprintf('Mode %d', k); end
        end

        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        TFmat = db_template('timefreqmat');
        TFmat.TF       = W;                              % complex [K x nTime x nFreq]
        TFmat.Time     = M.Time;
        TFmat.Freqs    = Freqs(:)';
        TFmat.RowNames = RowNames;
        TFmat.Measure  = 'none';                         % complex (un-measured)
        TFmat.Method   = 'morlet';
        TFmat.DataType = 'matrix';
        if isfield(M, 'SurfaceFile')
            TFmat.SurfaceFile = M.SurfaceFile;
        end
        TFmat.Comment  = sprintf('EigenWavelet (%d modes, %d freqs) | %s', size(W,1), numel(Freqs), sInput.Comment);
        TFmat.nAvg     = 1;
        TFmat = bst_history('add', TFmat, 'eigenmodes_wavelet', TFmat.Comment);

        OutFile = bst_process('GetNewFilename', StudyDir, 'timefreq_eigenwavelet');
        bst_save(OutFile, TFmat, 'v6');
        db_add_data(iStudyOut, OutFile, TFmat);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>

        bst_report('Info', sProcess, sInput, sprintf('Wavelet tensor: %d modes x %d times x %d freqs (%.1f-%.1f Hz).', ...
            size(W,1), size(W,2), size(W,3), Freqs(1), Freqs(end)));
    end
end
```

- [ ] **Step 4: Run the options test to verify it passes**

Run: `test_process_eigenmodes_wavelet_options`
Expected: PASS — `ALL TESTS PASSED: test_process_eigenmodes_wavelet_options`.

- [ ] **Step 5: Lint**

`check_matlab_code` on the process file — only standard Brainstorm idioms (`varargout`, stale `%#ok`) acceptable.

- [ ] **Step 6: Live end-to-end validation on the EigenSmoke substrate (verifies the Run path)**

The options test does not exercise `Run`. The `EigenSmoke` protocol (`SmokeS`) already has a `matrix_eigentransform` coefficient file from the M1 smoke. In MATLAB (MCP `evaluate_matlab_code`):

```matlab
sMat = bst_process('CallProcess','process_select_files_matrix', [], [], 'subjectname','SmokeS');
iM = find(arrayfun(@(s) ~isempty(strfind(s.FileName,'matrix_eigentransform')), sMat), 1);
assert(~isempty(iM), 'No matrix_eigentransform found in SmokeS (run the M1 transform smoke first).');
sW = bst_process('CallProcess','process_eigenmodes_wavelet', sMat(iM), [], ...
    'flo',2, 'fhi',0, 'nfreqs',40, 'morletfc',1, 'morletfwhmtc',3);
assert(~isempty(sW) && ~isempty(sW(1).FileName), 'wavelet process produced no output');
T = load(file_fullpath(sW(1).FileName));
fprintf('W TF size=%s complex=%d Freqs=[%.1f %.1f] rows=%d\n', ...
    mat2str(size(T.TF)), ~isreal(T.TF), T.Freqs(1), T.Freqs(end), numel(T.RowNames));
assert(ndims(T.TF)==3 && ~isreal(T.TF) && all(isfinite(T.TF(:))), 'tensor must be a finite complex [K x nTime x nFreq]');
fprintf('LIVE VALIDATION OK\n');
```

Expected: prints `LIVE VALIDATION OK` with a finite **complex** `[K × nTime × nFreq]` tensor. If a Brainstorm-glue detail is off (timefreq fields rejected, complex TF not retained on load, RowNames size), DEBUG and fix `Run` until this passes, then re-run the options test.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_wavelet.m dev/tests/test_process_eigenmodes_wavelet_options.m
git commit -m "$(cat <<'EOF'
Add process_eigenmodes_wavelet: complex (lambda,omega,t) tensor

Thin process: Morlet-transforms an eigenmode-coefficient matrix via
bst_eigenmodes_wavelet and saves the complex, lambda-labeled timefreq
tensor (Measure='none', Method='morlet'). Validated end-to-end on the
OMEGA EigenSmoke substrate.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 3: (Optional, manual) amplitude view

Not required to close the increment. After Task 2, in MATLAB, load `sW(1)` and view a frequency slice over `(λ, t)`: `imagesc(T.Time, 1:size(T.TF,1), 10*log10(abs(squeeze(T.TF(:,:,fIdx))).^2))` for an alpha-band `fIdx`, or open the file in Brainstorm's TF viewer (magnitude). Confirm a time-resolved `(λ,ω)` amplitude image (e.g. an alpha band that waxes/wanes). Save to `dev/tests/eigenwavelet_smoke.png`.

---

## Self-Review

**Spec coverage:**
- Pure complex tensor via `morlet_transform('n')` + permute to `[K×nTime×nFreq]` → Task 1. ✓
- Default log-spaced grid `logspace(log10(2), log10(min(100,0.4*sfreq)), 40)` → Task 1 (and tested). ✓
- Morlet `Fc=1`/`FwhmTc=3` defaults, configurable → Task 1 options + Task 2 options. ✓
- Process: matrix-in, complex λ-labeled timefreq out, Index 336.7, `Measure='none'`, `Method='morlet'`, `DataType='matrix'` → Task 2. ✓
- Amplitude/phase recoverable (unit-tested), complex retained → Task 1 test + Task 2 live validation (`~isreal`). ✓
- Viewable in TF viewer → Task 2 (timefreq file) + Task 3 optional. ✓
- Analyses deferred (dispersion, phase coherence, etc.) → no tasks, as intended. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete code; run steps have exact commands + expected output.

**Type/name consistency:** `bst_eigenmodes_wavelet(Coeffs, sfreq, Freqs, 'MorletFc',..,'MorletFwhmTc',..) -> [W, Freqs]` used identically in Task 1 test and Task 2 Run. Tensor shape `[K×nTime×nFreq]` consistent across the pure fn, the test asserts, the process `TF`, and the live validation. File prefix `timefreq_eigenwavelet`, `Method='morlet'`, `Measure='none'` consistent.
