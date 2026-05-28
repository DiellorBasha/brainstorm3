# Eigenmode Dispersion Analysis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** From a stationary `(λ,ω)` eigenmode power spectrum, fit a wave model (peak frequency ∝ √λ → speed `c`) and a diffusion model (bandwidth ∝ λ → diffusivity `α`) and report the better-fitting regime.

**Architecture:** A pure function (`bst_eigenmodes_dispersion`) reduces each mode's spectrum to peak-frequency + bandwidth, runs two weighted through-origin fits, and picks the regime by weighted R². A thin process (`process_eigenmodes_dispersion`) loads a `(λ,ω)` timefreq, extracts power per its `Measure`, gets `λ_k` from the file's `SurfaceFile`, calls the pure function, reports the result, and saves a per-mode features matrix.

**Tech Stack:** MATLAB, Brainstorm process system (`eval(macro_method)`), `in_bst_timefreq`, `in_tess_eigenmodes`, `db_template('matrixmat')`, `bst_save`/`db_add_data`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_eigenmodes_dispersion.m` (create) | Pure: `(P, lambdas, Freqs) → {PeakFreq, Bandwidth, c, alpha, R2wave, R2diff, Regime}`. No I/O. |
| `toolbox/process/functions/process_eigenmodes_dispersion.m` (create) | Process: timefreq-in → power-per-Measure + time-avg → λ from SurfaceFile → pure fn → report + save features matrix. |
| `dev/tests/test_eigenmodes_dispersion_pure.m` (create) | DB-free unit test (synthetic wave + diffusion). |
| `dev/tests/test_process_eigenmodes_dispersion_options.m` (create) | DB-free options test. |

**Run convention (tests):** MATLAB MCP `evaluate_matlab_code` calling the function name (script-style tests print `ALL TESTS PASSED`; NOT `run_matlab_test_file`). Lint via `check_matlab_code`.

Confirmed facts: `in_bst_timefreq(file, 1)` returns `.TF [K×nTime×nFreq]`, `.Freqs`, `.Measure` (`'none'`/`'magnitude'`/`'power'`), `.SurfaceFile`. `in_tess_eigenmodes(SurfaceFile)` returns `Eigenmodes` with `.Values`. `tess_sphere(642)` makes a manifold sphere; `tess_eigenmodes(V,F,'nModes',K)` computes eigenmodes (used only by the live validation).

---

## Task 1: Pure dispersion function `bst_eigenmodes_dispersion`

**Files:**
- Create: `toolbox/math/bst_eigenmodes_dispersion.m`
- Test: `dev/tests/test_eigenmodes_dispersion_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_dispersion_pure.m`:

```matlab
function test_eigenmodes_dispersion_pure
% Verify wave-vs-diffusion discrimination + parameter recovery on synthetic spectra.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% ---- Synthetic WAVE: peak frequency f* = c0*sqrt(lambda)/(2*pi), fixed bandwidth ----
K = 30; nF = 200; Freqs = linspace(1, 50, nF);
sl = linspace(0.3, 5, K)';          % sqrt(lambda) per mode
lambdas = sl.^2;
c0 = 30;                             % m/s
fpk = c0 * sl / (2*pi);             % peak freq (Hz), ~1.4..23.9
P = zeros(K, nF);
for k = 1:K
    P(k,:) = exp(-((Freqs - fpk(k)).^2) / (2 * 2^2));   % Gaussian peak, sigma 2 Hz
end
Out = bst_eigenmodes_dispersion(P, lambdas, Freqs);
assert(strcmp(Out.Regime, 'wave'), 'Expected wave regime, got %s.', Out.Regime);
assert(Out.R2wave > Out.R2diff, 'R2wave should exceed R2diff for a wave.');
assert(abs(Out.c - c0)/c0 < 0.15, 'Wave speed should recover c0 (got %.2f vs %.2f).', Out.c, c0);

% ---- Synthetic DIFFUSION: per-mode Lorentzian half-width gamma = alpha0*lambda ----
K2 = 30; nF2 = 200; Freqs2 = linspace(1, 50, nF2);
lambdas2 = linspace(0.5, 30, K2)';
alpha0 = 0.05;
P2 = zeros(K2, nF2);
for k = 1:K2
    g = alpha0 * lambdas2(k);                       % half-width (Hz)
    P2(k,:) = 1 ./ (Freqs2.^2 + g^2);               % Lorentzian (peaks at f=0)
end
Out2 = bst_eigenmodes_dispersion(P2, lambdas2, Freqs2);
assert(strcmp(Out2.Regime, 'diffusion'), 'Expected diffusion regime, got %s.', Out2.Regime);
assert(Out2.R2diff > Out2.R2wave, 'R2diff should exceed R2wave for diffusion.');

% ---- Shape / field checks ----
assert(isequal(size(Out.PeakFreq), [K 1]) && isequal(size(Out.Bandwidth), [K 1]), 'feature shapes');
assert(isfinite(Out.c) && isfinite(Out.alpha), 'c and alpha must be finite.');

fprintf('ALL TESTS PASSED: test_eigenmodes_dispersion_pure\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`): `test_eigenmodes_dispersion_pure`
Expected: FAIL — `Undefined function or variable 'bst_eigenmodes_dispersion'`.

- [ ] **Step 3: Write minimal implementation**

Create `toolbox/math/bst_eigenmodes_dispersion.m`:

```matlab
function Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, varargin)
% BST_EIGENMODES_DISPERSION: Wave-vs-diffusion discrimination on a (lambda,omega) spectrum.
%
% USAGE:  Out = bst_eigenmodes_dispersion(P, lambdas, Freqs)
%         Out = bst_eigenmodes_dispersion(P, lambdas, Freqs, 'MinPowerFrac', 0)
%
% DESCRIPTION:
%     Discriminates the dynamical regime behind a stationary eigenmode
%     (lambda, omega) power spectrum by reducing each mode's spectrum to two
%     features and testing how they scale with the spatial frequency lambda:
%
%       Wave (dispersion curve omega = c*sqrt(lambda)): the PEAK frequency
%         f*(k) scales as sqrt(lambda_k). Weighted through-origin fit
%         f* = a*sqrt(lambda) gives speed c = 2*pi*a.
%       Diffusion (wedge; Lorentzian half-width ~ lambda): the spectral
%         BANDWIDTH w(k) scales as lambda_k. Weighted through-origin fit
%         w = b*lambda gives diffusivity alpha = 2*pi*b (proportional).
%
%     Regime = the model with the higher weighted R^2. Subtraction is on power.
%     c is in m/s when lambdas are metric LBO eigenvalues (Brainstorm surfaces
%     are in metres); the regime decision and R^2 are unit-independent.
%
% INPUTS:
%     P       : [K x nFreq] non-negative power per eigenmode per frequency.
%     lambdas : [K x 1] eigenvalues (rad^2/m^2), aligned with rows of P.
%     Freqs   : [1 x nFreq] frequencies (Hz).
%
% OPTIONS (name-value):
%     'MinPowerFrac' : drop modes whose total power is below this fraction of
%                      the max per-mode total power (default 0 = keep all;
%                      power-weighting already down-weights weak modes).
%
% OUTPUTS:
%     Out.PeakFreq  [K x 1], Out.Bandwidth [K x 1], Out.Weights [K x 1]
%     Out.c (m/s), Out.alpha, Out.R2wave, Out.R2diff
%     Out.Regime ('wave'|'diffusion'), Out.Margin (|R2wave - R2diff|)
%
% SEE ALSO: bst_eigenmodes_wavelet, process_eigenmodes_dispersion

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
MinPowerFrac = 0;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'minpowerfrac', MinPowerFrac = varargin{i+1};
    end
end

P       = max(double(P), 0);        % power is non-negative
lambdas = double(lambdas(:));        % [K x 1]
Freqs   = double(Freqs(:)');         % [1 x nFreq]

%% ===== PER-MODE FEATURES =====
p = sum(P, 2);                        % [K x 1] total power per mode (weight)
psafe = p; psafe(psafe == 0) = eps;
[~, ipk] = max(P, [], 2);            % peak-frequency index per mode
PeakFreq = Freqs(ipk)';             % [K x 1]
fbar = (P * Freqs') ./ psafe;        % power-weighted mean frequency [K x 1]
fvar = (P * (Freqs'.^2)) ./ psafe - fbar.^2;
Bandwidth = sqrt(max(fvar, 0));     % power-weighted std of frequency [K x 1]

%% ===== WEIGHTS (valid modes only) =====
valid = (lambdas > 0) & isfinite(PeakFreq) & isfinite(Bandwidth) & (p > MinPowerFrac * max(p));
w = p; w(~valid) = 0;

sl = sqrt(max(lambdas, 0));         % [K x 1]

%% ===== WAVE FIT (f* = a*sqrt(lambda)) =====
a = sum(w .* sl .* PeakFreq) / max(sum(w .* lambdas), eps);
c = 2*pi*a;
R2wave = weighted_r2(sl, PeakFreq, w);

%% ===== DIFFUSION FIT (w = b*lambda) =====
b = sum(w .* lambdas .* Bandwidth) / max(sum(w .* lambdas.^2), eps);
alpha = 2*pi*b;
R2diff = weighted_r2(lambdas, Bandwidth, w);

%% ===== REGIME =====
if R2wave >= R2diff
    Regime = 'wave';
else
    Regime = 'diffusion';
end

Out = struct('PeakFreq', PeakFreq, 'Bandwidth', Bandwidth, 'Weights', w, ...
             'c', c, 'alpha', alpha, 'R2wave', R2wave, 'R2diff', R2diff, ...
             'Regime', Regime, 'Margin', abs(R2wave - R2diff));
end


%% ===== WEIGHTED SQUARED PEARSON CORRELATION =====
function r2 = weighted_r2(x, y, w)
    sw = sum(w);
    if sw <= 0
        r2 = 0; return;
    end
    mx = sum(w .* x) / sw;
    my = sum(w .* y) / sw;
    cxy = sum(w .* (x - mx) .* (y - my));
    cxx = sum(w .* (x - mx).^2);
    cyy = sum(w .* (y - my).^2);
    denom = cxx * cyy;
    if denom <= 0
        r2 = 0;
    else
        r2 = cxy^2 / denom;
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `test_eigenmodes_dispersion_pure`
Expected: PASS — `ALL TESTS PASSED: test_eigenmodes_dispersion_pure`.

- [ ] **Step 5: Lint**

`check_matlab_code` on `toolbox/math/bst_eigenmodes_dispersion.m` — no genuine errors.

- [ ] **Step 6: Commit**

```bash
git add toolbox/math/bst_eigenmodes_dispersion.m dev/tests/test_eigenmodes_dispersion_pure.m
git commit -m "$(cat <<'EOF'
Add bst_eigenmodes_dispersion: wave-vs-diffusion discrimination

Pure: per-mode peak-frequency (wave: f*~sqrt(lambda) -> speed c) and
power-weighted bandwidth (diffusion: width~lambda -> diffusivity alpha),
regime by weighted R^2. Covered by test_eigenmodes_dispersion_pure
(synthetic wave + diffusion).

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Task 2: Process `process_eigenmodes_dispersion`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_dispersion.m`
- Test: `dev/tests/test_process_eigenmodes_dispersion_options.m`

- [ ] **Step 1: Write the failing options test**

Create `dev/tests/test_process_eigenmodes_dispersion_options.m`:

```matlab
function test_process_eigenmodes_dispersion_options
% Verify the dispersion process is a timefreq-in / matrix-out Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_dispersion('GetDescription');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.8) < 1e-9, 'Index must be 336.8.');
assert(isequal(sProcess.InputTypes, {'timefreq'}), 'InputTypes must be {timefreq}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}), 'OutputTypes must be {matrix}.');
assert(isfield(sProcess.options, 'minpowerfrac'), 'Missing minpowerfrac option.');
assert(isfield(sProcess.options, 'label_info'), 'Missing label_info option.');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_dispersion_options\n');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `test_process_eigenmodes_dispersion_options`
Expected: FAIL — `Undefined function or variable 'process_eigenmodes_dispersion'`.

- [ ] **Step 3: Write the implementation**

Create `toolbox/process/functions/process_eigenmodes_dispersion.m`:

```matlab
function varargout = process_eigenmodes_dispersion( varargin )
% PROCESS_EIGENMODES_DISPERSION: Wave-vs-diffusion regime from a (lambda,omega) spectrum.
%
% USAGE:  sProcess = process_eigenmodes_dispersion('GetDescription')
%       OutputFiles = process_eigenmodes_dispersion('Run', sProcess, sInputs)
%
% DESCRIPTION:
%     Fits a wave model (peak frequency ~ sqrt(lambda) -> speed c) and a
%     diffusion model (bandwidth ~ lambda -> diffusivity alpha) to an eigenmode
%     (lambda, omega) power spectrum (e.g. an FFT/PSD of the eigenmode
%     coefficients, or the time-averaged wavelet tensor), and reports the
%     better-fitting regime. See bst_eigenmodes_dispersion.
%
% SEE ALSO: bst_eigenmodes_dispersion, process_eigenmodes_wavelet

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
    sProcess.Comment     = 'Eigenmode dispersion (wave vs diffusion)';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.8;   % after the wavelet tensor (336.7)
    sProcess.Description = '';
    sProcess.InputTypes  = {'timefreq'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;

    sProcess.options.minpowerfrac.Comment = 'Drop modes below this fraction of max per-mode power: ';
    sProcess.options.minpowerfrac.Type    = 'value';
    sProcess.options.minpowerfrac.Value   = {0, '', 2};

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Fits a wave (peak freq &prop; &radic;&lambda;) ' ...
        'and a diffusion (bandwidth &prop; &lambda;) model to the (&lambda;,&omega;) spectrum and reports the ' ...
        'better-fitting regime + speed c / diffusivity &alpha;.<BR>Input: an eigenmode (&lambda;,&omega;) timefreq spectrum.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = 'Eigenmode dispersion (wave vs diffusion)';
end


%% ===== RUN =====
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    MinPowerFrac = sProcess.options.minpowerfrac.Value{1};

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        TFmat = in_bst_timefreq(sInput.FileName, 1);   % LoadFull=1

        if iscell(TFmat.Freqs)
            bst_report('Error', sProcess, sInput, 'Frequency bands (cell Freqs) are not supported; use a frequency vector.');
            continue;
        end
        % Convert to power per the file's Measure (do not square already-power data)
        switch lower(TFmat.Measure)
            case 'power'
                Pw = abs(TFmat.TF);
            otherwise   % 'none' (complex) or 'magnitude'
                Pw = abs(TFmat.TF).^2;
        end
        P = reshape(mean(Pw, 2), size(Pw,1), size(Pw,3));   % [K x nFreq], averaged over time
        K = size(P, 1);

        if isempty(TFmat.SurfaceFile)
            bst_report('Error', sProcess, sInput, 'Timefreq has no SurfaceFile; cannot obtain eigenvalues.');
            continue;
        end
        [Em, isComputed] = in_tess_eigenmodes(TFmat.SurfaceFile);
        if ~isComputed
            bst_report('Error', sProcess, sInput, ['No eigenmodes on surface: ' TFmat.SurfaceFile]);
            continue;
        end
        if size(Em.Values, 1) < K
            bst_report('Error', sProcess, sInput, sprintf('Spectrum has %d modes but surface has only %d eigenvalues.', K, size(Em.Values,1)));
            continue;
        end
        lambdas = double(Em.Values(1:K));

        Out = bst_eigenmodes_dispersion(P, lambdas, TFmat.Freqs, 'MinPowerFrac', MinPowerFrac);

        bst_report('Info', sProcess, sInput, sprintf(...
            'Dispersion: regime=%s | c=%.3g m/s, alpha=%.3g | R2wave=%.3f, R2diff=%.3f (margin %.3f)', ...
            Out.Regime, Out.c, Out.alpha, Out.R2wave, Out.R2diff, Out.Margin));

        % Save per-mode features as a matrix file
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));
        M = db_template('matrixmat');
        M.Value       = [Out.PeakFreq'; Out.Bandwidth'];   % [2 x K]
        M.Time        = 1:K;
        M.Description  = {'PeakFreq (Hz)'; 'Bandwidth (Hz)'};
        M.SurfaceFile = TFmat.SurfaceFile;
        M.Comment     = sprintf('EigenDispersion [%s] c=%.3g a=%.3g R2w=%.2f R2d=%.2f | %s', ...
            Out.Regime, Out.c, Out.alpha, Out.R2wave, Out.R2diff, sInput.Comment);
        M.nAvg        = 1;
        M = bst_history('add', M, 'eigenmodes_dispersion', M.Comment);

        OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigendispersion');
        bst_save(OutFile, M, 'v6');
        db_add_data(iStudyOut, OutFile, M);
        OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
    end
end
```

- [ ] **Step 4: Run the options test to verify it passes**

Run: `test_process_eigenmodes_dispersion_options`
Expected: PASS — `ALL TESTS PASSED: test_process_eigenmodes_dispersion_options`.

- [ ] **Step 5: Lint**

`check_matlab_code` on the process file — only standard Brainstorm idioms (`varargout`, stale `%#ok`) acceptable.

- [ ] **Step 6: Live end-to-end validation (self-contained synthetic substrate)**

The options test does not exercise `Run`. Build a tiny throwaway protocol with a sphere surface that has eigenmodes + a synthetic WAVE `(λ,ω)` timefreq referencing it, run the process, and confirm it reports the `wave` regime. In MATLAB (MCP `evaluate_matlab_code`):

```matlab
gui_brainstorm('DeleteProtocol', 'DispSmoke');
gui_brainstorm('CreateProtocol', 'DispSmoke', 0, 0);
[~, iSubject] = db_add_subject('DS', [], 0, 0);
% sphere surface + eigenmodes, saved into the subject's anat folder
sSubject = bst_get('Subject', iSubject);
anatDir = bst_fileparts(file_fullpath(sSubject.FileName));
[V, F] = tess_sphere(642);
Em = tess_eigenmodes(V, F, 'nModes', 40, 'Verbose', 0);
SurfFile = fullfile(anatDir, 'tess_dispsphere.mat');
TessMat = struct('Vertices', V, 'Faces', F, 'Comment', 'dispsphere', 'Eigenmodes', Em);
bst_save(SurfFile, TessMat, 'v7');
iSurf = db_add_surface(iSubject, SurfFile, 'dispsphere'); %#ok<NASGU>
SurfRel = file_short(SurfFile);
% synthetic wave spectrum referencing that surface
K = Em.nModes; nF = 200; Freqs = linspace(1, 50, nF);
lam = double(Em.Values(1:K)); sl = sqrt(max(lam,0));
c0 = 4.0; fpk = c0 * sl / (2*pi);
TF = zeros(K, 1, nF);
for k = 1:K, TF(k,1,:) = exp(-((Freqs - fpk(k)).^2) / (2*2^2)); end
TFmat = db_template('timefreqmat');
TFmat.TF = TF; TFmat.Time = [0 1]; TFmat.Freqs = Freqs;
TFmat.Measure = 'power'; TFmat.Method = 'psd'; TFmat.DataType = 'matrix';
TFmat.SurfaceFile = SurfRel;
TFmat.RowNames = arrayfun(@(k) sprintf('Mode %d', k), (1:K)', 'uni', 0);
TFmat.Comment = 'Synthetic wave (lambda,omega)'; TFmat.nAvg = 1;
sSubject = bst_get('Subject', iSubject);          % refresh after db_add_surface
iStudyDS = sSubject.iStudy(1);                     % the subject's default study
StudyDir = bst_fileparts(file_fullpath(bst_get('Study', iStudyDS).FileName));
fTf = bst_process('GetNewFilename', StudyDir, 'timefreq_synthwave');
bst_save(fTf, TFmat, 'v6'); db_add_data(iStudyDS, fTf, TFmat);
% run dispersion on it
sTf = bst_process('CallProcess', 'process_select_files_timefreq', [], [], 'subjectname', 'DS');
iSel = find(arrayfun(@(s) ~isempty(strfind(s.FileName,'timefreq_synthwave')), sTf), 1);
sD = bst_process('CallProcess', 'process_eigenmodes_dispersion', sTf(iSel), [], 'minpowerfrac', 0);
assert(~isempty(sD) && ~isempty(sD(1).FileName), 'dispersion produced no output');
R = in_bst_matrix(sD(1).FileName);
fprintf('dispersion output: %s | rows=%d\n', R.Comment, size(R.Value,1));
assert(~isempty(strfind(R.Comment, '[wave]')), 'expected wave regime in the output comment');
fprintf('LIVE VALIDATION OK\n');
gui_brainstorm('DeleteProtocol', 'DispSmoke');
```

Expected: prints `LIVE VALIDATION OK` (the synthetic wave is classified `[wave]`). If the surface/timefreq registration glue needs adjustment (e.g. `db_add_surface` path handling, the subject's default study lookup), debug and fix the VALIDATION script (not `Run`) until it passes; if `Run` itself has a bug it surfaces here — fix `Run`, keep the contract. If the substrate registration proves intractable after genuine effort, report DONE_WITH_CONCERNS noting the Run path was not live-exercised (the pure unit test fully covers the dispersion math; the process is thin glue). Always `DeleteProtocol('DispSmoke')` at the end.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_dispersion.m dev/tests/test_process_eigenmodes_dispersion_options.m
git commit -m "$(cat <<'EOF'
Add process_eigenmodes_dispersion: report wave-vs-diffusion regime

Thin process: (lambda,omega) timefreq in -> power per Measure + time-avg,
eigenvalues from SurfaceFile -> bst_eigenmodes_dispersion -> reports
regime/c/alpha/R2 and saves per-mode peak-freq + bandwidth as a matrix.
Validated on a synthetic wave spectrum.

Co-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>
EOF
)"
```

---

## Self-Review

**Spec coverage:**
- Per-mode peak frequency + power-weighted bandwidth → Task 1 (`PeakFreq`, `Bandwidth`). ✓
- Wave fit `f*=a√λ` → `c=2πa`, `R²_wave` → Task 1. ✓
- Diffusion fit `w=bλ` → `α=2πb`, `R²_diff` → Task 1. ✓
- Regime by weighted R² + margin → Task 1 (`Regime`, `Margin`). ✓
- Power-weighting + `MinPowerFrac` + exclude `λ≤0` → Task 1 (`valid`, `w`). ✓
- Process timefreq-in, power per `Measure`, time-average, λ from `SurfaceFile`, report + save features matrix, Index 336.8 → Task 2. ✓
- Discrimination + speed recovery validated → Task 1 test (wave + diffusion). ✓
- Physiological mask / time-varying dispersion / template fitting → NOT built (deferred), as intended. ✓

**Placeholder scan:** No TBD/TODO; every code step has complete code; run steps have exact commands + expected output.

**Type/name consistency:** `bst_eigenmodes_dispersion(P, lambdas, Freqs, 'MinPowerFrac', …) -> Out{.PeakFreq,.Bandwidth,.c,.alpha,.R2wave,.R2diff,.Regime,.Margin}` used identically in Task 1 test and Task 2 Run. The process passes `TFmat.Freqs` (row) and `Em.Values(1:K)`. File prefix `matrix_eigendispersion`; Index 336.8 matches the options test.
