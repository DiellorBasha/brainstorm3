# Eigenmode Wiener Filter Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a self-contained 2-input process (data + empty-room) that estimates a per-`(k,f)` Wiener gain from the noise floor and applies it to the eigenmode-coefficient time series in the frequency domain.

**Architecture:** Reuse M2's already-tested gain math (`bst_eigenmodes_noisefloor.Gain`, extended to honor `Alpha`+`GainFloor`), extract M2's data+noise→coefficients/PSD setup into a shared helper so the gain is computed from the same data it filters, add a pure FFT-domain applier (`bst_bandpass_fft` idiom: mirror → fft → multiply by zero-phase magnitude → ifft → trim), and a thin process to wire them.

**Tech Stack:** MATLAB, Brainstorm process framework. Tests are script-style MATLAB functions run via the MATLAB MCP (`evaluate_matlab_code` to run, `check_matlab_code` for static analysis).

**Spec:** `docs/superpowers/specs/2026-05-28-eigenmode-wiener-design.md`

---

## File Structure

| File | Responsibility | Change |
|------|----------------|--------|
| `toolbox/math/bst_eigenmodes_noisefloor.m` | Joint `(λ,ω)` SNR / clean-PSD / **Wiener gain** | Modify: `Gain` honors `Alpha`+`GainFloor` |
| `toolbox/math/bst_eigenmodes_wiener.m` | Apply a per-`(k,f)` magnitude gain to a coefficient time series via the FFT | **Create** (pure) |
| `toolbox/process/functions/process_eigenmodes_denoise.m` | M2 noise-floor process + **shared coeffs/PSD helper** | Modify: extract `GetCoeffsAndPSD`, refactor `Run` |
| `toolbox/process/functions/process_eigenmodes_wiener.m` | Wiener-filter process (Index 336.95) | **Create** |
| `dev/tests/test_eigenmodes_noisefloor_pure.m` | Guards noise-floor math | Modify: add `Alpha`/`GainFloor` gain cases |
| `dev/tests/test_eigenmodes_wiener_pure.m` | Guards the FFT applier | **Create** |
| `dev/tests/test_process_eigenmodes_wiener_options.m` | Guards the process description/options | **Create** |

**Repo root:** `/Users/diellorbasha/workspace/research/code/brainstorm3`. All test run commands below assume the MATLAB MCP and `cd` into `dev/tests` first (tests `addpath` the repo root and start Brainstorm `nogui` themselves).

**Task order & dependencies:** Task 1 (noisefloor) and Task 2 (pure applier) are independent. Task 3 (helper extraction) is independent of 1 & 2. Task 4 (process) depends on 1, 2, and 3. Task 5 (optional live smoke) depends on 4.

---

## Task 1: Extend the noise-floor Wiener gain (Alpha + GainFloor)

**Files:**
- Modify: `toolbox/math/bst_eigenmodes_noisefloor.m`
- Test: `dev/tests/test_eigenmodes_noisefloor_pure.m`

- [ ] **Step 1: Add the new gain assertions to the existing test (these will fail first)**

In `dev/tests/test_eigenmodes_noisefloor_pure.m`, insert the following block immediately before the final `fprintf('ALL TESTS PASSED: ...')` line:

```matlab
% --- Gain shaping: GainFloor raises the floor, Alpha lowers the gain ---
Gbase = bst_eigenmodes_noisefloor(Pdata, N);                       % Alpha=1, GainFloor=0
Gflr  = bst_eigenmodes_noisefloor(Pdata, N, 'GainFloor', 0.2);
assert(all(Gflr.Gain(:) >= 0.2 - 1e-12), 'GainFloor should raise the gain floor.');
assert(all(Gflr.Gain(:) <= 1 + 1e-12),   'Gain must stay <= 1.');
Gov = bst_eigenmodes_noisefloor(Pdata, N, 'Alpha', 2);
assert(all(Gov.Gain(:) <= Gbase.Gain(:) + 1e-12), 'Over-subtraction (Alpha>1) should not increase the gain.');
% GainFloor outside [0,1] errors
threwGF = false;
try, bst_eigenmodes_noisefloor(Pdata, N, 'GainFloor', 1.5); catch, threwGF = true; end
assert(threwGF, 'GainFloor outside [0,1] should error.');
```

- [ ] **Step 2: Run the test and confirm it fails**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_eigenmodes_noisefloor_pure
```
Expected: FAIL — `GainFloor should raise the gain floor.` (the current `Gain` ignores `GainFloor`; the unknown option is silently dropped, so `Gflr.Gain` equals the unfloored gain whose min is 0).

- [ ] **Step 3: Implement the gain extension**

In `toolbox/math/bst_eigenmodes_noisefloor.m`, change the option-parse block (currently defaults `Alpha = 1; Floor = 0; SnrThresh = 1;` with a `switch`) to also parse `GainFloor`:

```matlab
%% ===== PARSE OPTIONS =====
Alpha = 1; Floor = 0; SnrThresh = 1; GainFloor = 0;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'alpha',     Alpha     = varargin{i+1};
        case 'floor',     Floor     = varargin{i+1};
        case 'snrthresh', SnrThresh = varargin{i+1};
        case 'gainfloor', GainFloor = varargin{i+1};
    end
end
if GainFloor < 0 || GainFloor > 1
    error('bst_eigenmodes_noisefloor: GainFloor must be in [0,1].');
end
if Alpha < 0
    error('bst_eigenmodes_noisefloor: Alpha must be >= 0.');
end
```

Then change the `Out.Gain` line in the `%% ===== COMBINE =====` block from:
```matlab
Out.Gain     = max(Pdata - Nnoise, 0) ./ max(Pdata, eps);
```
to:
```matlab
Out.Gain     = max( max(Pdata - Alpha .* Nnoise, 0) ./ max(Pdata, eps), GainFloor );
```

- [ ] **Step 4: Update the header docstring**

In the header comment of `bst_eigenmodes_noisefloor.m`, update the `Gain(k,f)` line and the OPTIONS list:
- Change the products line to:
  `%       Gain(k,f)     = max(max(Pdata - Alpha*Nnoise, 0)/Pdata, GainFloor)  (Wiener gain in [0,1])`
- Add under OPTIONS, after the `'Floor'` line:
  `%     'GainFloor' : lower bound for the Wiener gain (in [0,1]), default 0.`

- [ ] **Step 5: Run the test and confirm it passes**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_eigenmodes_noisefloor_pure
```
Expected: `ALL TESTS PASSED: test_eigenmodes_noisefloor_pure` (old assertions still pass — defaults `Alpha=1, GainFloor=0` reproduce the prior formula because the inner term is already `>= 0`).

- [ ] **Step 6: Static analysis**

Run (MATLAB MCP `check_matlab_code`) on `toolbox/math/bst_eigenmodes_noisefloor.m`. Expected: no new warnings.

- [ ] **Step 7: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_eigenmodes_noisefloor.m dev/tests/test_eigenmodes_noisefloor_pure.m
git commit -m "$(printf 'Extend eigenmode noise-floor Gain with Alpha + GainFloor\n\nWiener gain now G = max(max(P-Alpha*N,0)/P, GainFloor). Defaults\nAlpha=1, GainFloor=0 reproduce the prior formula byte-for-byte.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: Pure FFT-domain applier `bst_eigenmodes_wiener`

**Files:**
- Create: `toolbox/math/bst_eigenmodes_wiener.m`
- Test: `dev/tests/test_eigenmodes_wiener_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmodes_wiener_pure.m`:

```matlab
function test_eigenmodes_wiener_pure
% Verify FFT-domain application of a per-mode magnitude gain G(k,f).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sfreq = 200; T = 4;
t = 0:(1/sfreq):(T - 1/sfreq);
nTime = numel(t);
K  = 3;
fA = 10; fB = 40;
Coeffs    = repmat(cos(2*pi*fA*t) + cos(2*pi*fB*t), K, 1);   % [K x nTime]
GainFreqs = 0:1:(sfreq/2);                                   % [1 x nGF], ascending
nGF       = numel(GainFreqs);

% Power at frequency f via in-phase/quadrature projection.
powAt = @(x, f) (mean(x .* cos(2*pi*f*t)))^2 + (mean(x .* sin(2*pi*f*t)))^2;

% --- Identity: Gain == 1 returns the input ---
Yid = bst_eigenmodes_wiener(Coeffs, sfreq, ones(K, nGF), GainFreqs);
assert(isequal(size(Yid), size(Coeffs)), 'Output size must match input.');
assert(isreal(Yid), 'Output must be real.');
assert(max(abs(Yid(:) - Coeffs(:))) < 1e-6, 'Identity gain must return the input.');

% --- Null: Gain == 0 returns ~0 ---
Y0 = bst_eigenmodes_wiener(Coeffs, sfreq, zeros(K, nGF), GainFreqs);
assert(max(abs(Y0(:))) < 1e-6, 'Zero gain must return ~0.');

% --- Selectivity: pass fA (gain 1 in 5..15 Hz), kill fB ---
Gsel = zeros(K, nGF);
Gsel(:, GainFreqs >= 5 & GainFreqs <= 15) = 1;
Ysel = bst_eigenmodes_wiener(Coeffs, sfreq, Gsel, GainFreqs);
pA_in  = powAt(Coeffs(1,:), fA); pB_in  = powAt(Coeffs(1,:), fB);
pA_out = powAt(Ysel(1,:),  fA);  pB_out = powAt(Ysel(1,:),  fB);
assert(pA_out > 0.5  * pA_in, 'Passband tone fA should be largely preserved.');
assert(pB_out < 0.05 * pB_in, 'Stopband tone fB should be largely removed.');

% --- Zero-phase: a single passband tone keeps its phase (no quadrature, no inversion) ---
xin   = cos(2*pi*fA*t);
Gpass = double(GainFreqs >= 5 & GainFreqs <= 15);   % [1 x nGF]
yph   = bst_eigenmodes_wiener(xin, sfreq, Gpass, GainFreqs);
c_cos = mean(yph .* cos(2*pi*fA*t));
c_sin = mean(yph .* sin(2*pi*fA*t));
assert(abs(c_sin) < 0.05 * abs(c_cos), 'Zero-phase: quadrature component must be negligible.');
assert(c_cos > 0, 'Zero-phase: in-phase component must stay positive (no inversion).');

% --- Errors ---
threwRow = false;
try, bst_eigenmodes_wiener(Coeffs, sfreq, ones(K+1, nGF), GainFreqs); catch, threwRow = true; end
assert(threwRow, 'Row mismatch between Gain and Coeffs should error.');
threwLen = false;
try, bst_eigenmodes_wiener(Coeffs, sfreq, ones(K, nGF), GainFreqs(1:end-1)); catch, threwLen = true; end
assert(threwLen, 'GainFreqs length not matching Gain columns should error.');

fprintf('ALL TESTS PASSED: test_eigenmodes_wiener_pure\n');
end
```

- [ ] **Step 2: Run the test and confirm it fails**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_eigenmodes_wiener_pure
```
Expected: FAIL — `Undefined function 'bst_eigenmodes_wiener'`.

- [ ] **Step 3: Implement the pure applier**

Create `toolbox/math/bst_eigenmodes_wiener.m`:

```matlab
function Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs, varargin)
% BST_EIGENMODES_WIENER: Apply a per-mode magnitude gain G(k,f) to coefficients via the FFT.
%
% USAGE:  Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs)
%         Filtered = bst_eigenmodes_wiener(Coeffs, sfreq, Gain, GainFreqs, 'Mirror', true)
%
% DESCRIPTION:
%     Applies a frequency-dependent, per-eigenmode magnitude response (e.g. a
%     Wiener gain) to a coefficient time series, in the frequency domain. For
%     each mode k, the gain Gain(k,:) defined on the grid GainFreqs is
%     interpolated onto the (folded) FFT frequency grid and applied as a
%     zero-phase magnitude response: Filtered = real(ifft(fft(Coeffs).*H)).
%     Because the gain is interpolated on the folded (absolute) frequency, H is
%     Hermitian-symmetric, so the output is real and phase is preserved. The
%     signal is mirrored at both ends before the FFT (and trimmed after) to
%     suppress circular-convolution edge artifacts, as in bst_bandpass_fft.
%
% INPUTS:
%     Coeffs    : [K x nTime] real coefficient time series (rows = eigenmodes).
%     sfreq     : sampling frequency (Hz).
%     Gain      : [K x nGainFreq] magnitude gain per mode per frequency, in [0,1].
%     GainFreqs : [1 x nGainFreq] frequencies (Hz) for Gain's columns, ascending.
%
% OPTIONS (name-value):
%     'Mirror'  : reflect the signal at both ends before the FFT, default true.
%
% OUTPUT:
%     Filtered  : [K x nTime] real, same size as Coeffs.
%
% SEE ALSO: bst_eigenmodes_noisefloor, process_eigenmodes_wiener, bst_bandpass_fft

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
Mirror = true;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'mirror', Mirror = logical(varargin{i+1});
    end
end

Coeffs    = double(Coeffs);
Gain      = double(Gain);
GainFreqs = double(GainFreqs(:)');     % row

[K, nTime] = size(Coeffs);
if size(Gain, 1) ~= K
    error('bst_eigenmodes_wiener: Gain must have one row per eigenmode (got %d rows for %d modes).', size(Gain,1), K);
end
if numel(GainFreqs) ~= size(Gain, 2)
    error('bst_eigenmodes_wiener: GainFreqs length (%d) must match the number of Gain columns (%d).', numel(GainFreqs), size(Gain,2));
end

%% ===== MIRROR (edge handling) =====
domirror = Mirror && (nTime >= 2);
if domirror
    x = [fliplr(Coeffs), Coeffs, fliplr(Coeffs)];   % [K x 3*nTime]
else
    x = Coeffs;
end
N = size(x, 2);

%% ===== BUILD ZERO-PHASE MAGNITUDE RESPONSE H(k,f) =====
% Folded (absolute) frequency of each FFT bin, in [0, Nyquist].
fvec = (0:N-1) * (sfreq / N);
hi   = fvec > (sfreq / 2);
fvec(hi) = sfreq - fvec(hi);                          % fold negative-frequency bins
% Interpolate each mode's gain onto the folded grid (columns = modes).
Hcols = interp1(GainFreqs, Gain.', fvec(:), 'linear');   % [N x K], NaN outside GainFreqs
below = fvec(:) < GainFreqs(1);
above = fvec(:) > GainFreqs(end);
if any(below), Hcols(below, :) = repmat(Gain(:,1).',   sum(below), 1); end
if any(above), Hcols(above, :) = repmat(Gain(:,end).', sum(above), 1); end
Hcols = max(0, min(1, Hcols));                        % clamp to [0,1]
H = Hcols.';                                          % [K x N], Hermitian-symmetric

%% ===== APPLY =====
Y = real(ifft(fft(x, [], 2) .* H, [], 2));

%% ===== TRIM MIRROR =====
if domirror
    Filtered = Y(:, nTime + (1:nTime));
else
    Filtered = Y;
end
end
```

- [ ] **Step 4: Run the test and confirm it passes**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_eigenmodes_wiener_pure
```
Expected: `ALL TESTS PASSED: test_eigenmodes_wiener_pure`.

- [ ] **Step 5: Static analysis**

Run `check_matlab_code` on `toolbox/math/bst_eigenmodes_wiener.m`. Expected: no warnings (the `.* H` is element-wise on same-size `[K x N]` matrices — no implicit-expansion lint).

- [ ] **Step 6: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_eigenmodes_wiener.m dev/tests/test_eigenmodes_wiener_pure.m
git commit -m "$(printf 'Add bst_eigenmodes_wiener: FFT-domain per-mode gain applier\n\nZero-phase, Hermitian-symmetric magnitude response interpolated onto\nthe folded FFT grid, with edge mirroring (bst_bandpass_fft idiom).\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: Extract M2's coeffs/PSD setup into a shared helper

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes_denoise.m`
- Test (guard, unchanged): `dev/tests/test_process_eigenmodes_denoise_options.m`

This refactor is behavior-preserving: the denoise process must produce the same output files as before. There is no automated process-level test for M2's `Run`; the structural options test guards `GetDescription`, and the optional live smoke (Task 5) re-checks `Run`.

- [ ] **Step 1: Add the shared subfunction `GetCoeffsAndPSD`**

In `toolbox/process/functions/process_eigenmodes_denoise.m`, add this local function (place it directly after the `Run` function, before the `SaveTF` helper). It is the current `Run` setup moved verbatim, repackaged to return a struct (or `[]` on error):

```matlab
%% ===== SHARED: COMMON-CHANNEL KERNEL, COEFFICIENTS, DATA/NOISE PSDs =====
function S = GetCoeffsAndPSD(sProcess, sInputsA, sInputsB, nModesOpt, WinLen) %#ok<DEFNU>
    S = [];
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
    ChanFileA = bst_get('ChannelFileForStudy', sStudyA.FileName);
    if isempty(ChanFileA)
        bst_report('Error', sProcess, sInputsA, 'No channel file for the data study (Files A).');
        return;
    end
    ChA = in_bst_channel(ChanFileA);
    DA  = in_bst_data(sInputsA(1).FileName);
    if isstruct(DA.F)
        bst_report('Error', sProcess, sInputsA, 'Files A must be imported data (not raw). Import a block first.');
        return;
    end
    iA = good_channel(ChA.Channel, DA.ChannelFlag, 'MEG');
    if isempty(iA), iA = good_channel(ChA.Channel, DA.ChannelFlag, 'EEG'); end
    if isempty(iA)
        bst_report('Error', sProcess, sInputsA, 'No good MEG or EEG channels found in data (Files A).');
        return;
    end

    [sStudyB,~,~,~] = bst_get('Study', sInputsB(1).iStudy);
    ChanFileB = bst_get('ChannelFileForStudy', sStudyB.FileName);
    if isempty(ChanFileB)
        bst_report('Error', sProcess, sInputsB, 'No channel file for the empty-room study (Files B).');
        return;
    end
    ChB = in_bst_channel(ChanFileB);
    DB  = in_bst_data(sInputsB(1).FileName);
    if isstruct(DB.F)
        bst_report('Error', sProcess, sInputsB, 'Files B (empty-room) must be imported data (not raw). Import a block first.');
        return;
    end
    iB = good_channel(ChB.Channel, DB.ChannelFlag, 'MEG');
    if isempty(iB), iB = good_channel(ChB.Channel, DB.ChannelFlag, 'EEG'); end
    if isempty(iB)
        bst_report('Error', sProcess, sInputsB, 'No good MEG or EEG channels found in empty-room (Files B).');
        return;
    end

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

    % ===== PACK =====
    S = struct();
    S.Coeffs      = thD;                  % [K x nTime] data coefficient time series
    S.Pdata       = reshape(TFd, K, []);  % [K x nFreq]
    S.Nnoise      = reshape(TFn, K, []);  % [K x nFreq]
    S.Freqs       = Fv(:)';               % [1 x nFreq]
    S.lambdas     = lambdas;              % [K x 1]
    S.Phi         = Phi;                  % [nVert x K]
    S.SurfaceFile = SurfaceFile;
    S.K           = K;
    S.Time        = DA.Time;              % [1 x nTime]
    S.sfreq       = sfA;
    S.Info        = Info;
end
```

- [ ] **Step 2: Replace the body of `Run` to call the helper**

Replace the entire current `Run` function body (everything from `OutputFiles = {};` through the closing `end` of `Run`, i.e. the option-parse + DATA STUDY + CHANNELS + KERNEL + COEFFICIENTS + COMBINE + SAVE + report blocks) with:

```matlab
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

    S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen);
    if isempty(S), return; end

    % ===== COMBINE =====
    Out = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha', Alpha, 'Floor', FloorFrac, 'SnrThresh', SnrThresh);

    % ===== SAVE =====
    RowNames = cell(S.K,1);
    for k = 1:S.K, RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, S.lambdas(k)); end
    [sStudyOut, iStudyOut] = bst_get('Study', sInputsA(1).iStudy);
    StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

    OutputFiles{end+1} = SaveTF(Out.SNR, S.Freqs, RowNames, S.Time, StudyDir, iStudyOut, ...
        sprintf('EigenSNR (%d modes) | %s', S.K, sInputsA(1).Comment), 'timefreq_eigensnr', S.SurfaceFile); %#ok<AGROW>
    OutputFiles{end+1} = SaveTF(Out.CleanPSD, S.Freqs, RowNames, S.Time, StudyDir, iStudyOut, ...
        sprintf('EigenCleanPSD (%d modes, a=%.1f) | %s', S.K, Alpha, sInputsA(1).Comment), 'timefreq_eigencleanpsd', S.SurfaceFile); %#ok<AGROW>

    bst_report('Info', sProcess, sInputsA, sprintf('Denoise: %d modes (condition %.1f), median reliable-mode cutoff K*=%d at SNR>=%.1f.', ...
        S.K, S.Info.ConditionNumber, round(median(Out.Kstar)), SnrThresh));
end
```

Leave `GetDescription`, `FormatComment`, and `SaveTF` unchanged.

- [ ] **Step 3: Run the options test (structural guard) and confirm it still passes**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_process_eigenmodes_denoise_options
```
Expected: `ALL TESTS PASSED: test_process_eigenmodes_denoise_options`.

- [ ] **Step 4: Confirm the helper resolves via the macro dispatch**

No automated test calls `GetCoeffsAndPSD` through the dispatch without a loaded DB, so verify the dispatch resolves the subfunction name (the failure mode is a typo / it not being reachable via `eval(macro_method)`). The probe treats *any* error other than an "undefined/unknown function" error as success (it means our code was reached):
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3
brainstorm nogui;
resolved = true;
try
    process_eigenmodes_denoise('GetCoeffsAndPSD', struct('Comment','probe'), ...
        struct('iStudy',1,'FileName','x','Comment','y'), struct('iStudy',1,'FileName','z'), 0, 2);
catch ME
    m = lower(ME.message);
    resolved = isempty(strfind(m, 'undefined')) && isempty(strfind(m, 'unknown'));
    fprintf('dispatch error (expected, non-fatal): %s\n', ME.message);
end
fprintf('GetCoeffsAndPSD resolved=%d\n', resolved);
```
Expected: `GetCoeffsAndPSD resolved=1`. If `resolved=0` (an "Undefined function 'GetCoeffsAndPSD'" error), the subfunction is not reachable via the dispatch — fix before proceeding.

- [ ] **Step 5: Static analysis**

Run `check_matlab_code` on `toolbox/process/functions/process_eigenmodes_denoise.m`. Expected: no new warnings (`GetCoeffsAndPSD` carries `%#ok<DEFNU>` since it is reached only via the `eval(macro_method)` dispatch).

- [ ] **Step 6: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/process/functions/process_eigenmodes_denoise.m
git commit -m "$(printf 'Extract shared GetCoeffsAndPSD helper from eigenmode denoise\n\nMoves the data+noise -> common-channel kernel, coefficients and Welch\nPSD setup into a macro-dispatched subfunction so the Wiener filter can\nreuse it. Run is refactored to call it; behavior-preserving.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: New process `process_eigenmodes_wiener`

**Files:**
- Create: `toolbox/process/functions/process_eigenmodes_wiener.m`
- Test: `dev/tests/test_process_eigenmodes_wiener_options.m`

- [ ] **Step 1: Write the failing options test**

Create `dev/tests/test_process_eigenmodes_wiener_options.m`:

```matlab
function test_process_eigenmodes_wiener_options
% Verify the Wiener-filter process is a 2-input data->matrix Sources process.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_eigenmodes_wiener('GetDescription');
assert(strcmp(sProcess.Category, 'Custom'),  'Category must be Custom.');
assert(strcmp(sProcess.SubGroup, 'Sources'), 'SubGroup must be Sources.');
assert(abs(sProcess.Index - 336.95) < 1e-9,  'Index must be 336.95.');
assert(isequal(sProcess.InputTypes,  {'data','raw'}), 'InputTypes must be {data,raw}.');
assert(isequal(sProcess.OutputTypes, {'matrix'}),     'OutputTypes must be {matrix}.');
assert(sProcess.nInputs == 2, 'nInputs must be 2 (data + empty-room).');
for f = {'nmodes','noisewin','alpha','gainfloor','domirror','dorecon','savegain'}
    assert(isfield(sProcess.options, f{1}), 'Missing option: %s', f{1});
end
assert(sProcess.options.alpha.Value{1} == 1,     'Default alpha must be 1.');
assert(sProcess.options.gainfloor.Value{1} == 0, 'Default gainfloor must be 0.');
assert(sProcess.options.domirror.Value == 1,     'Default domirror must be 1 (on).');
assert(sProcess.options.dorecon.Value  == 0,     'Default dorecon must be 0 (off).');
assert(sProcess.options.savegain.Value == 0,     'Default savegain must be 0 (off).');
fprintf('ALL TESTS PASSED: test_process_eigenmodes_wiener_options\n');
end
```

- [ ] **Step 2: Run the test and confirm it fails**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_process_eigenmodes_wiener_options
```
Expected: FAIL — `Undefined function 'process_eigenmodes_wiener'`.

- [ ] **Step 3: Implement the process**

Create `toolbox/process/functions/process_eigenmodes_wiener.m`:

```matlab
function varargout = process_eigenmodes_wiener( varargin )
% PROCESS_EIGENMODES_WIENER: Wiener (noise-floor) spectral filter of eigenmode coefficients.
%
% USAGE:  sProcess = process_eigenmodes_wiener('GetDescription')
%       OutputFiles = process_eigenmodes_wiener('Run', sProcess, sInputsA, sInputsB)
%
% DESCRIPTION:
%     Files A = data recording(s) (surface head model + eigenmodes required).
%     Files B = empty-room recording(s). Builds the data's eigenmode coefficient
%     time series and the data/noise Welch PSDs (shared with the denoise
%     process), forms a per-(k,f) Wiener gain G = max(max(P-Alpha*N,0)/P, Gmin)
%     via bst_eigenmodes_noisefloor, and applies it to the coefficient time
%     series in the frequency domain (bst_eigenmodes_wiener). Both inputs must
%     be imported (not raw).
%
% SEE ALSO: bst_eigenmodes_wiener, bst_eigenmodes_noisefloor, process_eigenmodes_denoise

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
    sProcess.Comment     = 'Eigenmode Wiener filter';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 336.95;   % after the coefficient filter (336.9)
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'matrix'};
    sProcess.nInputs     = 2;        % Files A = data, Files B = empty-room
    sProcess.nMinFiles   = 1;

    sProcess.options.nmodes.Comment   = 'Number of eigenmodes (0 = auto, min of channels and available): ';
    sProcess.options.nmodes.Type      = 'value';
    sProcess.options.nmodes.Value     = {0, '', 0};

    sProcess.options.noisewin.Comment = 'Welch window length: ';
    sProcess.options.noisewin.Type    = 'value';
    sProcess.options.noisewin.Value   = {2, 's', 2};

    sProcess.options.alpha.Comment    = 'Over-subtraction factor alpha (>=1): ';
    sProcess.options.alpha.Type       = 'value';
    sProcess.options.alpha.Value      = {1, '', 2};

    sProcess.options.gainfloor.Comment = 'Gain floor Gmin (0..1): ';
    sProcess.options.gainfloor.Type    = 'value';
    sProcess.options.gainfloor.Value   = {0, '', 2};

    sProcess.options.domirror.Comment = 'Mirror signal edges (reduce FFT edge artifacts)';
    sProcess.options.domirror.Type    = 'checkbox';
    sProcess.options.domirror.Value   = 1;

    sProcess.options.dorecon.Comment  = 'Also reconstruct filtered vertex sources (Q = Phi * theta_filt)';
    sProcess.options.dorecon.Type     = 'checkbox';
    sProcess.options.dorecon.Value    = 0;

    sProcess.options.savegain.Comment = 'Also save the Wiener gain spectrum G(k,f)';
    sProcess.options.savegain.Type    = 'checkbox';
    sProcess.options.savegain.Value   = 0;

    sProcess.options.label_info.Comment = ['<FONT color="#777777">Files A = data, Files B = empty-room. ' ...
        'Estimates a per-mode Wiener gain G(&lambda;,&omega;) from the noise floor and applies it to the ' ...
        'eigenmode coefficient time series in the frequency domain.<BR>Both inputs must be imported; requires ' ...
        'a surface head model + eigenmodes.</FONT>'];
    sProcess.options.label_info.Type    = 'label';
end


%% ===== FORMAT COMMENT =====
function Comment = FormatComment(sProcess) %#ok<DEFNU>
    a = sProcess.options.alpha.Value{1};
    g = sProcess.options.gainfloor.Value{1};
    Comment = sprintf('Eigenmode Wiener filter (a=%.1f, Gmin=%.2f)', a, g);
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
    GainFloor = sProcess.options.gainfloor.Value{1};
    DoMirror  = sProcess.options.domirror.Value;
    DoRecon   = sProcess.options.dorecon.Value;
    SaveGain  = sProcess.options.savegain.Value;

    % Shared setup: coefficients + data/noise PSDs (errors already reported inside)
    S = process_eigenmodes_denoise('GetCoeffsAndPSD', sProcess, sInputsA, sInputsB, nModesOpt, WinLen);
    if isempty(S), return; end

    % Wiener gain from the noise floor, applied in the frequency domain
    NF = bst_eigenmodes_noisefloor(S.Pdata, S.Nnoise, 'Alpha', Alpha, 'GainFloor', GainFloor);
    G  = NF.Gain;   % [K x nFreq] in [0,1]
    Filtered = bst_eigenmodes_wiener(S.Coeffs, S.sfreq, G, S.Freqs, 'Mirror', logical(DoMirror));

    RowNames = cell(S.K,1);
    for k = 1:S.K, RowNames{k} = sprintf('Mode %d (lam=%.1f)', k, S.lambdas(k)); end
    [sStudyOut, iStudyOut] = bst_get('Study', sInputsA(1).iStudy);
    StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

    % --- Filtered coefficients (matrix) ---
    Mout = db_template('matrixmat');
    Mout.Value       = Filtered;
    Mout.Time        = S.Time;
    Mout.Description  = RowNames;
    Mout.SurfaceFile = S.SurfaceFile;
    Mout.Comment     = sprintf('EigenWiener (%d modes, a=%.1f, Gmin=%.2f) | %s', S.K, Alpha, GainFloor, sInputsA(1).Comment);
    Mout.nAvg        = 1;
    Mout = bst_history('add', Mout, 'eigenmodes_wiener', Mout.Comment);
    OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigenwiener');
    bst_save(OutFile, Mout, 'v6');
    db_add_data(iStudyOut, OutFile, Mout);
    OutputFiles{end+1} = file_short(OutFile);

    % --- Optional Wiener gain spectrum (timefreq sidecar; not in OutputFiles) ---
    if SaveGain
        TFmat = db_template('timefreqmat');
        TFmat.TF          = reshape(G, [S.K, 1, numel(S.Freqs)]);
        TFmat.Freqs       = S.Freqs;
        TFmat.Time        = [S.Time(1), S.Time(end)];
        TFmat.RowNames    = RowNames;
        TFmat.Measure     = 'power';
        TFmat.Method      = 'psd';
        TFmat.DataType    = 'matrix';
        TFmat.SurfaceFile = S.SurfaceFile;
        TFmat.Comment     = sprintf('EigenWienerGain (%d modes, a=%.1f, Gmin=%.2f) | %s', S.K, Alpha, GainFloor, sInputsA(1).Comment);
        TFmat.nAvg        = 1;
        TFmat = bst_history('add', TFmat, 'eigenmodes_wiener', TFmat.Comment);
        GainFile = bst_process('GetNewFilename', StudyDir, 'timefreq_eigenwienergain');
        bst_save(GainFile, TFmat, 'v6');
        db_add_data(iStudyOut, GainFile, TFmat);
    end

    % --- Optional filtered vertex reconstruction (results sidecar; not in OutputFiles) ---
    if DoRecon
        Q = S.Phi * Filtered;   % [nVert x nTime]
        ResMat = db_template('resultsmat');
        ResMat.ImageGridAmp  = Q;
        ResMat.ImagingKernel = [];
        ResMat.nComponents   = 1;
        ResMat.Time          = S.Time;
        ResMat.SurfaceFile   = S.SurfaceFile;
        ResMat.HeadModelType = 'surface';
        ResMat.ColormapType  = 'source';
        ResMat.Comment       = sprintf('EigenWiener recon (%d modes) | %s', S.K, sInputsA(1).Comment);
        ResMat.nAvg          = 1;
        ResMat = bst_history('add', ResMat, 'eigenmodes_wiener', ResMat.Comment);
        ReconFile = bst_process('GetNewFilename', StudyDir, 'results_eigenwiener');
        bst_save(ReconFile, ResMat, 'v6');
        db_add_data(iStudyOut, ReconFile, ResMat);
    end

    bst_report('Info', sProcess, sInputsA, sprintf('Wiener-filtered %d eigenmode coefficients (alpha=%.1f, Gmin=%.2f).', S.K, Alpha, GainFloor));
end
```

- [ ] **Step 4: Run the options test and confirm it passes**

Run:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests; test_process_eigenmodes_wiener_options
```
Expected: `ALL TESTS PASSED: test_process_eigenmodes_wiener_options`.

- [ ] **Step 5: Static analysis**

Run `check_matlab_code` on `toolbox/process/functions/process_eigenmodes_wiener.m`. Expected: no warnings (`%#ok<DEFNU>` on `GetDescription`/`FormatComment`/`Run`).

- [ ] **Step 6: Run the full eigenmode test suite (regression)**

Run each in sequence; each must print its `ALL TESTS PASSED` line:
```matlab
cd /Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests
test_eigenmodes_noisefloor_pure
test_eigenmodes_wiener_pure
test_process_eigenmodes_wiener_options
test_process_eigenmodes_denoise_options
test_eigenmodes_filter_gain_pure
test_eigenmodes_filter_pure
```
Expected: six `ALL TESTS PASSED:` lines (the noise-floor, filter, and denoise-options tests confirm no regression from the refactor).

- [ ] **Step 7: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/process/functions/process_eigenmodes_wiener.m dev/tests/test_process_eigenmodes_wiener_options.m
git commit -m "$(printf 'Add process_eigenmodes_wiener (Index 336.95)\n\nTwo-input (data + empty-room) process: reuses the shared coeffs/PSD\nhelper, builds the noise-floor Wiener gain, and applies it via\nbst_eigenmodes_wiener. Optional gain-spectrum and vertex-recon sidecars.\n\nCo-Authored-By: Claude Opus 4.7 (1M context) <noreply@anthropic.com>')"
```

---

## Task 5 (OPTIONAL): Live smoke on a disposable protocol

This validates the full pipeline against the Brainstorm DB. It is optional: if substrate setup proves intractable, report `DONE_WITH_CONCERNS` — the pure + options tests (Tasks 1–4) carry the math and structure. The user has scheduled a dedicated testing/debugging pass after this increment, so this smoke is a sanity check, not the gate.

**Files:** none committed (throwaway protocol + scratch script under `dev/`, deleted after).

- [ ] **Step 1: Build (or reuse) a substrate with a surface head model + eigenmodes**

Reuse the OMEGA sub-0002 substrate used by prior increments (the protocol/subject that already has an imported data block, a surface head model, and precomputed eigenmodes). If present, identify an imported data file (Files A) and an imported empty-room block (Files B). If no empty-room block exists, import a short noise segment, or synthesize one by importing a copy of a baseline segment. Confirm both are imported (not raw) and share channel names.

- [ ] **Step 2: Run the Wiener process via bst_process**

Run (MATLAB MCP `evaluate_matlab_code`, substituting the real file paths/indices found in Step 1):
```matlab
% sFilesA = data block; sFilesB = empty-room block (bst_process('GetInputStruct', ...) or db lookups)
sOut = bst_process('CallProcess', 'process_eigenmodes_wiener', sFilesA, sFilesB, ...
    'nmodes', 0, 'noisewin', 2, 'alpha', 1, 'gainfloor', 0.1, ...
    'domirror', 1, 'dorecon', 1, 'savegain', 1);
M  = in_bst_matrix(sOut(1).FileName);
fprintf('Filtered coeffs: [%d x %d], finite=%d\n', size(M.Value,1), size(M.Value,2), all(isfinite(M.Value(:))));
```
Expected: one matrix output; `Filtered coeffs: [K x nTime], finite=1`. A `timefreq_eigenwienergain` and a `results_eigenwiener` file appear in the tree (sidecars, not in `sOut`).

- [ ] **Step 3: Sanity-check the effect**

Compare the filtered coefficients' high-frequency / low-SNR power against an unfiltered transform of the same data (e.g. via `process_eigenmodes_transform` + `process_fft`), and confirm the saved gain spectrum lies in `[0,1]` with `min >= gainfloor`. This is a qualitative check (no fixed expected numbers for real data).

- [ ] **Step 4: Tear down**

Delete the throwaway protocol/files created in Step 1 (use `gui_brainstorm('DeleteProtocol', ...)` if a dedicated protocol was made; otherwise delete only the files added during the smoke). Do not delete pre-existing OMEGA data. Remove any scratch script. Confirm `git status` shows only the intended source/test changes from Tasks 1–4.

---

## Final review

After Task 4 (and optionally Task 5), dispatch a final code-review subagent over the whole branch diff, then use **superpowers:finishing-a-development-branch**. The user has stated the next phase is a dedicated testing & debugging pass over the eigenmode pipeline, so default to **Merge locally + push to fork** unless the user says otherwise.
