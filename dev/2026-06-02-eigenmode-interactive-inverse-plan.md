# Eigenmode Interactive Inverse Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add "Eigenmode source mapping" as a selectable method in the interactive "Compute sources [2018]" dialog, driven by the same `Compute()` core the process box uses, with Brainstorm's per-modality whitener reused verbatim.

**Architecture:** Eigenmode becomes one more `InverseMethod` (`'eigenmode'`) inside `process_inverse_2018('Compute', …)`. The dialog (`panel_inverse_2018`) gains a dedicated radio + spectral-prior sub-panel + coefficients checkbox, shown only when the active head model has `isEigenmode=1` (standard methods disabled). Two new standalone helpers carry the eigenmode-specific work: `bst_noise_whitener` (verbatim copy of the per-modality whitener) and `bst_eigenmode_reconstruct` (`Φ·M̃`). The process file `process_eigenmodes_inverse` becomes a thin wrapper over `Compute()`. `bst_inverse_linear_2018` is **not modified**.

**Tech Stack:** MATLAB, Brainstorm process-plugin system, Java Swing GUI panels, MATLAB MCP for running tests.

---

## Background references (read before starting)

- Spec: `dev/2026-06-02-eigenmode-interactive-inverse-design.md`
- Solver core: `toolbox/inverse/bst_inverse_eigenmodes.m` — `'SolvePure'` entry at line 32:
  `Kernel = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)`
- Spectral prior: `toolbox/math/bst_eigenmode_prior.m`
- Eigenmodes I/O: `toolbox/io/in_tess_eigenmodes.m` → `[Eigenmodes, isComputed]`, `Eigenmodes.Vectors` `[nVert x nModes]`
- Whitener source (to copy verbatim): `toolbox/inverse/bst_inverse_linear_2018.m`
  - per-modality block: lines **234–446**
  - `CROSS_COVARIANCE_CHANNELTYPES` flag: line **196**
  - subfunction `truncate_and_regularize_covariance`: lines **1098–1234**
  - subfunction `cov1para_local`: lines **1236–1361**
- Interactive dialog: `toolbox/inverse/panel_inverse_2018.m`
- Compute core: `toolbox/process/functions/process_inverse_2018.m`
- Existing batch process: `toolbox/process/functions/process_eigenmodes_inverse.m`
- Test style: `dev/tests/test_inverse_eigenmodes_pure.m` (function script; `addpath` repo root; `brainstorm nogui`; `assert`; final `disp('ALL TESTS PASSED')`)

**Running a test (MATLAB MCP):** use the `run_matlab_file` tool on the absolute path of the test `.m` file, or in the MATLAB console:
`run('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/tests/<test>.m')`
Expected on success: console prints `ALL TESTS PASSED`.

**Branch:** before Task 1, create `feature/eigenmode-interactive-inverse` off `development`.

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git checkout development && git checkout -b feature/eigenmode-interactive-inverse
```

---

## File Structure

| File | Responsibility |
|------|----------------|
| `toolbox/inverse/bst_noise_whitener.m` (new) | Per-modality whitener `iW`, verbatim from `bst_inverse_linear_2018` |
| `toolbox/inverse/bst_eigenmode_reconstruct.m` (new) | `Φ·M̃` cortex kernel from a mode-space kernel |
| `toolbox/inverse/panel_inverse_2018.m` (modify) | Eigenmode radio, prior sub-panel, coefficients checkbox, gating, getters |
| `toolbox/process/functions/process_inverse_2018.m` (modify) | `isEigenmode` detection, `case 'eigenmode'`, coefficients node |
| `toolbox/process/functions/process_eigenmodes_inverse.m` (modify) | Thin wrapper over `Compute()` |
| `dev/tests/test_noise_whitener_pure.m` (new) | Whitening-property + per-modality-structure test |
| `dev/tests/test_eigenmode_reconstruct_pure.m` (new) | `Φ·M̃` shape/value test |
| `dev/tests/test_eigenmodes_inverse_e2e.m` (modify) | Drive the `'eigenmode'` path through `Compute()` |

---

## Task 1: `bst_noise_whitener` — verbatim per-modality whitener

**Files:**
- Create: `toolbox/inverse/bst_noise_whitener.m`
- Test: `dev/tests/test_noise_whitener_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_noise_whitener_pure.m`:

```matlab
function test_noise_whitener_pure
% Verify bst_noise_whitener (verbatim per-modality whitener from bst_inverse_linear_2018):
%   - single modality: iW whitens the covariance (iW*C*iW' ~ I)
%   - two modalities: per-modality whitening, zero cross-modality blocks
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% --- single modality (4 EEG channels): whitening property ---
rng_local = [0.7 0.2 0.1 0.0; 0.2 0.8 0.15 0.05; 0.1 0.15 0.9 0.1; 0.0 0.05 0.1 0.6];
C1 = rng_local * rng_local';                 % SPD
types1 = {'EEG','EEG','EEG','EEG'};
iW1 = bst_noise_whitener(C1, types1, 'reg', 0.1);
assert(isequal(size(iW1), [4 4]), 'iW must be [nCh x nCh].');
W = iW1 * C1 * iW1';
% After whitening the regularized covariance, diagonal dominates and is ~unit scale
assert(all(isfinite(iW1(:))), 'iW must be finite.');
assert(max(abs(W - diag(diag(W))), [], 'all') < 0.5, 'Off-diagonals must be suppressed by whitening.');

% --- two modalities (2 MEG MAG + 2 EEG): cross-modality blocks zeroed ---
C2 = [2.0 0.3 0.4 0.1;
      0.3 1.7 0.2 0.2;
      0.4 0.2 0.9 0.25;
      0.1 0.2 0.25 1.1];
types2 = {'MEG MAG','MEG MAG','EEG','EEG'};
iW2 = bst_noise_whitener(C2, types2, 'reg', 0.1);
% Cross-modality blocks (rows 1:2 vs cols 3:4) must be exactly zero
assert(max(abs(iW2(1:2,3:4)), [], 'all') == 0, 'Cross-modality whitener block must be zero.');
assert(max(abs(iW2(3:4,1:2)), [], 'all') == 0, 'Cross-modality whitener block must be zero.');
assert(all(isfinite(iW2(:))), 'iW must be finite.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_noise_whitener_pure.m` via the MATLAB MCP.
Expected: FAIL — `Undefined function 'bst_noise_whitener'`.

- [ ] **Step 3: Create `bst_noise_whitener.m`**

Create `toolbox/inverse/bst_noise_whitener.m`. The body wraps the **verbatim** per-modality whitener from `bst_inverse_linear_2018.m`. Build a local `OPTIONS`-shaped struct from the args so the pasted block compiles unchanged, paste lines **237–446** verbatim, then return `iWw_noise`. Append the two subfunctions verbatim.

```matlab
function iW = bst_noise_whitener(NoiseCov, ChannelTypes, NoiseMethod, NoiseReg, FourthMoment, nSamples)
% BST_NOISE_WHITENER: Per-modality noise whitener iW = C^(-1/2).
%
% USAGE:  iW = bst_noise_whitener(NoiseCov, ChannelTypes, NoiseMethod, NoiseReg, FourthMoment, nSamples)
%
% DESCRIPTION:
%   Verbatim extraction of the per-modality whitener used by
%   bst_inverse_linear_2018 (lines 234-446 + subfunctions
%   truncate_and_regularize_covariance and cov1para_local). Regularizes and
%   whitens each modality (channel type) separately, with cross-modality
%   covariances zeroed, exactly as the standard 2018 inverse does. Kept as a
%   standalone function so the eigenmode inverse can reuse the same whitener
%   without depending on the standard solver.
%
%   NoiseCov     : [nCh x nCh] noise covariance, good channels only
%   ChannelTypes : 1 x nCh cell of channel type strings (e.g. {'MEG MAG',...})
%   NoiseMethod  : 'reg' | 'diag' | 'none' | 'shrink' | 'median'
%   NoiseReg     : scalar in [0,1] (used by 'reg')
%   FourthMoment : [nCh x nCh] (only required for 'shrink'); default zeros
%   nSamples     : scalar (only required for 'shrink'); default []
%
% NOTE: This duplicates logic in bst_inverse_linear_2018. Do not edit the math
%       here independently; if the standard whitener changes, re-sync verbatim.
%       Future consolidation (make bst_inverse_linear_2018 call this) is tracked
%       as tech debt in dev/2026-06-02-eigenmode-interactive-inverse-design.md.
%
% Authors: Diellor Basha, 2026

if (nargin < 5) || isempty(FourthMoment)
    FourthMoment = zeros(size(NoiseCov));
end
if (nargin < 6)
    nSamples = [];
end
% Assemble the OPTIONS-shaped struct the pasted block expects
OPTIONS.NoiseCovMat.NoiseCov     = NoiseCov;
OPTIONS.NoiseCovMat.FourthMoment = FourthMoment;
OPTIONS.NoiseCovMat.nSamples     = nSamples;
OPTIONS.ChannelTypes             = ChannelTypes;
OPTIONS.NoiseMethod              = NoiseMethod;
OPTIONS.NoiseReg                 = NoiseReg;

% Verbatim from bst_inverse_linear_2018.m line 196
CROSS_COVARIANCE_CHANNELTYPES = false;

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 237-446 =====
% (paste lines 237 through 446 here, unchanged)
% ===== END VERBATIM COPY =====

iW = iWw_noise;
end

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 1098-1234 =====
% (paste truncate_and_regularize_covariance here, unchanged)
% ===== END VERBATIM COPY =====

% ===== BEGIN VERBATIM COPY: bst_inverse_linear_2018.m lines 1236-1361 =====
% (paste cov1para_local here, unchanged)
% ===== END VERBATIM COPY =====
```

Paste instructions (do literally, no edits to the pasted lines):
1. Open `toolbox/inverse/bst_inverse_linear_2018.m`.
2. Copy lines **237–446** into the `BEGIN/END VERBATIM COPY: … 237-446` region. These lines reference `OPTIONS.NoiseCovMat.NoiseCov`, `OPTIONS.ChannelTypes`, `OPTIONS.NoiseMethod`, `OPTIONS.NoiseReg`, `FourthMoment`, `nSamples`, `CROSS_COVARIANCE_CHANNELTYPES` — all defined above — and produce `iWw_noise`.
3. Copy lines **1098–1234** (`truncate_and_regularize_covariance`) into its region.
4. Copy lines **1236–1361** (`cov1para_local`) into its region.

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_noise_whitener_pure.m` via the MATLAB MCP.
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Lint**

Run `/lint-matlab toolbox/inverse/bst_noise_whitener.m` (or `checkcode` via MATLAB MCP). Expected: no errors (Brainstorm idiom warnings are acceptable).

- [ ] **Step 6: Commit**

```bash
git add toolbox/inverse/bst_noise_whitener.m dev/tests/test_noise_whitener_pure.m
git commit -m "feat(inverse): add bst_noise_whitener (verbatim per-modality whitener)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `bst_eigenmode_reconstruct` — cortex kernel `Φ·M̃`

**Files:**
- Create: `toolbox/inverse/bst_eigenmode_reconstruct.m`
- Test: `dev/tests/test_eigenmode_reconstruct_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_eigenmode_reconstruct_pure.m`:

```matlab
function test_eigenmode_reconstruct_pure
% Verify bst_eigenmode_reconstruct: cortex kernel = Phi(:,1:K) * ModeKernel.
% Numeric-Phi mode keeps this DB-free.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nVert = 20; K = 6; nCh = 9;
Phi = zeros(nVert, K+2);                 % surface carries more modes than used
for v = 1:nVert
    for k = 1:(K+2)
        Phi(v,k) = cos(v * k * pi / 17);
    end
end
ModeKernel = reshape(linspace(-1, 1, K*nCh), K, nCh);   % [K x nCh]

KVert = bst_eigenmode_reconstruct(Phi, ModeKernel);
assert(isequal(size(KVert), [nVert nCh]), 'Cortex kernel must be [nVert x nCh].');
assert(max(abs(KVert - Phi(:,1:K) * ModeKernel), [], 'all') < 1e-12, 'Must equal Phi(:,1:K)*ModeKernel.');

% Error when surface has fewer modes than the kernel
threw = false;
try
    bst_eigenmode_reconstruct(Phi(:,1:K-1), ModeKernel);
catch
    threw = true;
end
assert(threw, 'Must error when surface modes < kernel modes.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_eigenmode_reconstruct_pure.m`.
Expected: FAIL — `Undefined function 'bst_eigenmode_reconstruct'`.

- [ ] **Step 3: Create `bst_eigenmode_reconstruct.m`**

```matlab
function ImagingKernel = bst_eigenmode_reconstruct(SurfaceOrPhi, ModeKernel)
% BST_EIGENMODE_RECONSTRUCT: Reconstruct a vertex-space imaging kernel from a
% mode-space kernel: ImagingKernel = Phi(:,1:K) * ModeKernel.
%
% USAGE:  ImagingKernel = bst_eigenmode_reconstruct(SurfaceFile, ModeKernel)
%         ImagingKernel = bst_eigenmode_reconstruct(Phi,         ModeKernel)
%
%   SurfaceOrPhi : cortex SurfaceFile (eigenmodes loaded via in_tess_eigenmodes)
%                  OR a numeric eigenvector matrix Phi [nVert x nModes]
%   ModeKernel   : [K x nGoodCh] mode-space kernel (from bst_inverse_eigenmodes)
%   ImagingKernel: [nVert x nGoodCh] cortex kernel
%
% Authors: Diellor Basha, 2026

K = size(ModeKernel, 1);
if ischar(SurfaceOrPhi)
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceOrPhi);
    if ~isComputed
        error(['No eigenmodes on surface: ' SurfaceOrPhi '. Run "Compute eigenmodes" first.']);
    end
    Phi = double(Eig.Vectors);
else
    Phi = double(SurfaceOrPhi);
end
if size(Phi, 2) < K
    error('Surface has fewer eigenmodes (%d) than kernel modes (%d).', size(Phi,2), K);
end
ImagingKernel = Phi(:, 1:K) * ModeKernel;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_eigenmode_reconstruct_pure.m`.
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_eigenmode_reconstruct.m dev/tests/test_eigenmode_reconstruct_pure.m
git commit -m "feat(inverse): add bst_eigenmode_reconstruct (Phi*M cortex kernel)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Dialog — eigenmode radio, prior sub-panel, coefficients checkbox, gating, getters

**Files:**
- Modify: `toolbox/inverse/panel_inverse_2018.m`

All edits follow the existing Swing patterns in this file. There is no headless unit test for Swing panels; verification is the smoke test in Step 8 plus the e2e in Task 6.

- [ ] **Step 1: Extend `CreatePanel` signature with `isEigenmode`**

At line 32, change:
```matlab
function [bstPanelNew, panelName] = CreatePanel(Modalities, isShared, HeadModelType, nSamplesNoise, nSamplesData) %#ok<DEFNU>
```
to:
```matlab
function [bstPanelNew, panelName] = CreatePanel(Modalities, isShared, HeadModelType, nSamplesNoise, nSamplesData, isEigenmode) %#ok<DEFNU>
```
Immediately after `panelName = 'InverseOptions';` (line 33), add a default so the process-call path (which omits the arg) is safe:
```matlab
    if (nargin < 6) || isempty(isEigenmode)
        isEigenmode = 0;
    end
```
Update the USAGE comment block (lines 4-5) to mention the trailing `isEigenmode` arg on the interactive form.

- [ ] **Step 2: Add the eigenmode radio to the Method panel**

In the `==== PANEL: METHOD ====` block, after the MEM radio creation (after line 147, inside the same block before line 148 `% Default selection`), add:
```matlab
        jRadioMethodEig = gui_component('radio', jPanelMethod, 'br', 'Eigenmode source mapping', jGroupMethod, '', @Method_Callback, []);
```
In the default-selection `switch` (lines 149-154), add a case:
```matlab
            case 'eigenmode', jRadioMethodEig.setSelected(1);
```
After the existing gating (after line 167), add the mutually-exclusive gating:
```matlab
        % Eigenmode head model: only the eigenmode method is valid; disable the others
        if ~isProcess
            if isEigenmode
                jRadioMethodEig.setSelected(1);
                jRadioMethodMn.setEnabled(0);
                jRadioMethodBf.setEnabled(0);
                jRadioMethodDip.setEnabled(0);
                if ~isempty(jRadioMethodMem); jRadioMethodMem.setEnabled(0); end
            else
                jRadioMethodEig.setVisible(0);
                jRadioMethodEig.setEnabled(0);
            end
        else
            jRadioMethodEig.setVisible(0);
            jRadioMethodEig.setEnabled(0);
        end
```

- [ ] **Step 3: Add the eigenmode Measure + Spectral-prior sub-panel**

After the `==== PANEL: MEASURE BEAMFORMER ====` block (after line 197), add a new panel modeled on the MN measure panel:
```matlab
    % ==== PANEL: EIGENMODE (measure + spectral prior) ====
    jPanelEig = gui_river([1,1], [0,6,6,6], 'Eigenmode');
        % Measure
        gui_component('label', jPanelEig, [], 'Measure:', [], '', [], []);
        jGroupEigMeasure = ButtonGroup();
        jRadioEigMne     = gui_component('radio', jPanelEig, [],   'MNE',     jGroupEigMeasure, '', @(h,ev)UpdatePanel(1), []);
        jRadioEigDspm    = gui_component('radio', jPanelEig, [],   'dSPM',    jGroupEigMeasure, '', @(h,ev)UpdatePanel(1), []);
        jRadioEigSloreta = gui_component('radio', jPanelEig, [],   'sLORETA', jGroupEigMeasure, '', @(h,ev)UpdatePanel(1), []);
        jRadioEigDspm.setSelected(1);
        % Spectral prior
        gui_component('label', jPanelEig, 'br', 'Spectral prior:', [], '', [], []);
        jGroupEigPrior = ButtonGroup();
        jRadioEigPriorLog   = gui_component('radio', jPanelEig, [], 'Log (2026)', jGroupEigPrior, '', [], []);
        jRadioEigPriorFlat  = gui_component('radio', jPanelEig, [], 'Flat',       jGroupEigPrior, '', [], []);
        jRadioEigPriorPower = gui_component('radio', jPanelEig, [], 'Power (1/f)',jGroupEigPrior, '', [], []);
        jRadioEigPriorLog.setSelected(1);
        % Apply noise whitening (default ON: replicate the standard inverse procedure)
        jCheckEigWhiten = gui_component('checkbox', jPanelEig, 'br', 'Apply noise whitening (recommended)', [], '', [], []);
        jCheckEigWhiten.setSelected(1);
        % Also save coefficients
        jCheckEigCoeff = gui_component('checkbox', jPanelEig, 'br', 'Also save eigenmode coefficients', [], '', [], []);
    c.gridy = 2;
    jPanelLeft.add(jPanelEig, c);
```
Note: `jPanelMeasureMN`, `jPanelMeasureBf`, and `jPanelEig` all use `c.gridy = 2`; `UpdatePanel` (Step 6) shows exactly one at a time, matching how MN/Bf already share the slot.

- [ ] **Step 4: Register the new controls in the `ctrl` struct**

In the `ctrl = struct(...)` (lines 370-414), after the `'jRadioMethodDip', jRadioMethodDip, ...` line add:
```matlab
            'jRadioMethodEig', jRadioMethodEig, ...
```
After the `'jRadioMethodBfNai', jRadioMethodBfNai, ...` line add:
```matlab
            ... % ==== PANEL: EIGENMODE ====
            'jRadioEigMne',       jRadioEigMne, ...
            'jRadioEigDspm',      jRadioEigDspm, ...
            'jRadioEigSloreta',   jRadioEigSloreta, ...
            'jRadioEigPriorLog',  jRadioEigPriorLog, ...
            'jRadioEigPriorFlat', jRadioEigPriorFlat, ...
            'jRadioEigPriorPower',jRadioEigPriorPower, ...
            'jCheckEigWhiten',    jCheckEigWhiten, ...
            'jCheckEigCoeff',     jCheckEigCoeff, ...
```

- [ ] **Step 5: Update `Method_Callback` so selecting eigenmode is inert for orientation**

No change needed to `Method_Callback` logic itself (it only reacts to MN/Dip). Leave as-is. (Documented here so the implementer does not invent changes.)

- [ ] **Step 6: Update `UpdatePanel` visibility logic**

In `UpdatePanel` (the `if isForced` block, lines 494-506), after `isLinear = ...` add an eigenmode flag and adjust the panel visibility lines. Replace the left-panel visibility lines (498-499) and add eigenmode handling:

Find:
```matlab
            jPanelMeasureMN.setVisible(isLinear && jRadioMethodMn.isSelected());
            jPanelMeasureBf.setVisible(isLinear && jRadioMethodBf.isSelected());
```
Replace with:
```matlab
            isEig = jRadioMethodEig.isSelected();
            jPanelMeasureMN.setVisible(isLinear && jRadioMethodMn.isSelected() && ~isEig);
            jPanelMeasureBf.setVisible(isLinear && jRadioMethodBf.isSelected() && ~isEig);
            jPanelEig.setVisible(isEig);
```
Then adjust the SNR / depth / output visibility so they behave for eigenmode. Find (lines 504-506):
```matlab
            jPanelSnr.setVisible(isLinear && ~jRadioMethodDip.isSelected() && ~jRadioMethodBf.isSelected());
            jPanelDepth.setVisible(isLinear && jRadioMethodMn.isSelected() && ~jRadioMnSloreta.isSelected());
            jPanelOutput.setVisible(isLinear && ~isProcess);
```
Replace with:
```matlab
            jPanelSnr.setVisible(isLinear && ~jRadioMethodDip.isSelected() && ~jRadioMethodBf.isSelected());
            jPanelDepth.setVisible(isLinear && jRadioMethodMn.isSelected() && ~jRadioMnSloreta.isSelected() && ~isEig);
            jPanelOutput.setVisible(isLinear && ~isProcess && ~isEig);
```
(Eigenmode keeps the SNR panel — it uses SNR for regularization — but hides depth weighting and the kernel/full output toggle; output is always a kernel-only cortex node plus the optional coefficients matrix.)

- [ ] **Step 7: Update the getters**

In `GetSelectedMethod` (lines 728-747), add an eigenmode branch before the closing `end`:
```matlab
    elseif isfield(ctrl, 'jRadioMethodEig') && ~isempty(ctrl.jRadioMethodEig) && ctrl.jRadioMethodEig.isSelected()
        Method = 'eigenmode';
        if ctrl.jRadioEigMne.isSelected()
            Measure = 'mne';
        elseif ctrl.jRadioEigDspm.isSelected()
            Measure = 'dspm';
        elseif ctrl.jRadioEigSloreta.isSelected()
            Measure = 'sloreta';
        end
```

In `GetMethodComment` (lines 750-768), add a case in the `switch`:
```matlab
        case 'eigenmode'
            switch (lower(Measure))
                case 'mne',     Comment = 'Eigen-MNE';
                case 'dspm',    Comment = 'Eigen-dSPM';
                case 'sloreta', Comment = 'Eigen-sLORETA';
                otherwise,      Comment = 'Eigenmode';
            end
```

In `GetPanelContents`, inside the `if isLinear` branch, after `[s.InverseMethod, s.InverseMeasure] = GetSelectedMethod(ctrl);` (line 639), add eigenmode-specific fields:
```matlab
        % Eigenmode-specific options
        if strcmpi(s.InverseMethod, 'eigenmode')
            if ctrl.jRadioEigPriorLog.isSelected()
                s.EigenmodePrior = 'log';
            elseif ctrl.jRadioEigPriorFlat.isSelected()
                s.EigenmodePrior = 'flat';
            elseif ctrl.jRadioEigPriorPower.isSelected()
                s.EigenmodePrior = 'power';
            end
            s.SaveCoefficients = ctrl.jCheckEigCoeff.isSelected();
            s.EigenmodeWhiten  = ctrl.jCheckEigWhiten.isSelected();
        end
```

- [ ] **Step 8: Smoke-test the panel builds with `isEigenmode`**

Create a throwaway check (run via MATLAB MCP, do not commit):
```matlab
if ~brainstorm('status'); brainstorm nogui; end
[p, name] = panel_inverse_2018('CreatePanel', {'MEG GRAD'}, 0, 'surface', [], [], 1);
disp(class(p)); disp(name);
```
Expected: prints a panel object class and `InverseOptions` with no error.

- [ ] **Step 9: Lint**

Run `/lint-matlab toolbox/inverse/panel_inverse_2018.m`. Expected: no new errors.

- [ ] **Step 10: Commit**

```bash
git add toolbox/inverse/panel_inverse_2018.m
git commit -m "feat(inverse): add eigenmode method to Compute sources dialog

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `Compute()` — detect eigenmode, dispatch, coefficients node

**Files:**
- Modify: `toolbox/process/functions/process_inverse_2018.m`

- [ ] **Step 1: Add eigenmode-related defaults to `Def_OPTIONS`**

In `Compute` (lines 135-142), add three fields to the `Def_OPTIONS` struct (so the wrapper and field-copy paths are safe):
```matlab
        'ComputeKernel',       1, ...
        'EigenmodePrior',      'log', ...
        'EigenmodeWhiten',     1, ...
        'SaveCoefficients',    0, ...
        'nModes',              0);
```
(Adjust the trailing comma on the previous `'ComputeKernel', 1` line to `, ...`.)
`EigenmodeWhiten` defaults to 1 so every code path (interactive, batch wrapper,
scripted) replicates the standard whitened inverse unless explicitly turned off.

- [ ] **Step 2: Detect `isEigenmode` in the study loop and pass to the dialog**

In the modality/headmodel loop (lines 190-227), after the `HeadModelType` assignment block (after line 206), add detection of an eigenmode head model on the first study:
```matlab
        if (i == 1)
            hmFile_i = sChanStudies(i).HeadModel(sChanStudies(i).iHeadModel).FileName;
            hmFlag_i = in_bst_headmodel(hmFile_i, 0, 'isEigenmode');
            isEigenmode = isfield(hmFlag_i, 'isEigenmode') && ~isempty(hmFlag_i.isEigenmode) && hmFlag_i.isEigenmode;
        end
```
Before the loop (near line 186 where `HeadModelType = 'surface';`), initialize:
```matlab
    isEigenmode = false;
```
At the dialog call (line 245), pass `isEigenmode` as the final arg:
```matlab
        sMethod = gui_show_dialog('Compute sources', @panel_inverse_2018, 1, [], AllMod, isShared, HeadModelType, nSamplesNoise, nSamplesData, isEigenmode);
```

- [ ] **Step 3: Add the `case 'eigenmode'` to the solver switch**

In the `switch( OPTIONS.InverseMethod )` (line 665), add a new case after the `'mem'` case and before `otherwise` (before line 703). At this point `HeadModel.Gain` is the good-channel, SSP-projected, average-referenced eigenmode leadfield `L_tilde` `[nGoodCh x K]`; `OPTIONS.NoiseCovMat.NoiseCov` is the matching good-channel, projected, avg-ref'd covariance; `OPTIONS.ChannelTypes` is set (line 663). SSP/avg-ref are already folded in, so the projector passed to the solver is identity.

```matlab
            case 'eigenmode'
                % Eigenmode leadfield: solve in mode space, reconstruct to cortex.
                % L_tilde is already SSP-projected + avg-ref'd + good-channel only.
                L_tilde = double(HeadModel.Gain);                 % [nGoodCh x K]
                % Eigenvalues + surface come from the composed head model file
                HMeig = in_bst_headmodel(HeadModelFile, 0, 'Eigenvalues', 'SurfaceFile', 'nModes');
                lambdas = double(HMeig.Eigenvalues(:));
                K = size(L_tilde, 2);
                % Optional mode cap
                if isfield(OPTIONS, 'nModes') && ~isempty(OPTIONS.nModes) && OPTIONS.nModes > 0 && OPTIONS.nModes < K
                    K = OPTIONS.nModes;
                    L_tilde = L_tilde(:, 1:K);
                    lambdas = lambdas(1:K);
                end
                % Whitener: ON by default (replicate the standard inverse procedure).
                % Per-modality whitener on the (already-projected) good-channel covariance.
                % When EigenmodeWhiten is off, use identity (pure sensor-space solve);
                % SSP / bad channels / avg-ref are still folded into L_tilde + data.
                doWhiten = ~isfield(OPTIONS, 'EigenmodeWhiten') || isempty(OPTIONS.EigenmodeWhiten) || OPTIONS.EigenmodeWhiten;
                if doWhiten
                    FourthMoment = [];
                    nSmp = [];
                    if isfield(OPTIONS.NoiseCovMat, 'FourthMoment'); FourthMoment = OPTIONS.NoiseCovMat.FourthMoment; end
                    if isfield(OPTIONS.NoiseCovMat, 'nSamples');     nSmp = OPTIONS.NoiseCovMat.nSamples; end
                    iW = bst_noise_whitener(OPTIONS.NoiseCovMat.NoiseCov, OPTIONS.ChannelTypes, ...
                        OPTIONS.NoiseMethod, OPTIONS.NoiseReg, FourthMoment, nSmp);
                else
                    iW = eye(size(L_tilde, 1));
                end
                % SSP already folded into L_tilde and covariance -> identity projector here
                ProjEig = eye(size(L_tilde, 1));
                % Regularization from SNR (fixed SNR field), default 3
                if isfield(OPTIONS, 'SnrFixed') && ~isempty(OPTIONS.SnrFixed)
                    snrVal = OPTIONS.SnrFixed;
                else
                    snrVal = 3;
                end
                ModeKernel = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, ProjEig, ...
                    OPTIONS.InverseMeasure, OPTIONS.EigenmodePrior, 1, snrVal, false);   % [K x nGoodCh]
                % Reconstruct cortex kernel Phi*M~
                Results = struct();
                Results.ImagingKernel = bst_eigenmode_reconstruct(HMeig.SurfaceFile, ModeKernel); % [nVert x nGoodCh]
                Results.ImageGridAmp  = [];
                Results.nComponents   = 1;
                Results.Function      = ['eigenmode_' OPTIONS.InverseMeasure];
                % Stash mode-space data for the optional coefficients node
                Results.EigenModeKernel = ModeKernel;
                Results.Eigenvalues     = lambdas(1:K);
                Results.nModes          = K;
                OPTIONS.FunctionName    = Results.Function;
```

- [ ] **Step 4: Save the optional coefficients matrix node**

The results-node save block runs after the switch (lines 711+). After the cortex results node is registered in the DB (after line 810, at the end of the per-study results handling, before the loop continues), add the coefficients node. Locate the end of the `% ===== REGISTER NEW FILE =====` block (line ~810, after `bst_set('Study', iStudy, sStudy);`) and add:

```matlab
        % ===== EIGENMODE COEFFICIENTS (optional) =====
        if strcmpi(OPTIONS.InverseMethod, 'eigenmode') && isfield(OPTIONS, 'SaveCoefficients') && OPTIONS.SaveCoefficients ...
                && isfield(Results, 'EigenModeKernel') && ~isempty(DataFile)
            DataMatCoeff = in_bst_data(DataFile);
            if ~isstruct(DataMatCoeff.F)   % imported (non-raw) data only
                theta = Results.EigenModeKernel * double(DataMatCoeff.F(GoodChannel, :));  % [K x nTime]
                MatMat = db_template('matrixmat');
                MatMat.Value       = theta;
                MatMat.Time        = DataMatCoeff.Time;
                MatMat.nAvg        = nAvg;  MatMat.Leff = Leff;
                MatMat.SurfaceFile = ResultsMat.SurfaceFile;
                MatMat.Comment     = sprintf('EigenCoeffs %s (%d modes, %s)', upper(OPTIONS.InverseMeasure), Results.nModes, OPTIONS.EigenmodePrior);
                RowNames = cell(Results.nModes, 1);
                for k = 1:Results.nModes
                    RowNames{k} = sprintf('Mode %d (lam=%.3g)', k, Results.Eigenvalues(k));
                end
                MatMat.Description = RowNames;
                MatMat = bst_history('add', MatMat, 'compute', sprintf('Eigenmode coefficients %s, %d modes, prior=%s', OPTIONS.InverseMeasure, Results.nModes, OPTIONS.EigenmodePrior));
                MatFile = bst_process('GetNewFilename', bst_fileparts(file_fullpath(sStudy.FileName)), 'matrix_eigencoeffs');
                bst_save(MatFile, MatMat, 'v6');
                db_add_data(iStudy, MatFile, MatMat);
                OutputFiles{end+1} = file_short(MatFile);
            else
                bst_report('Warning', 'process_inverse_2018', [], 'Eigenmode coefficients require imported data; skipped for raw.');
            end
        end
```
Note: verify `nAvg`, `Leff`, `iStudy`, `sStudy`, `GoodChannel`, `DataFile`, `ResultsMat` are all in scope at that point in `Compute` (they are used by the surrounding results block). If `db_template('matrixmat')` lacks a field used above, mirror the field set in `process_eigenmodes_inverse.m:144-153`.

- [ ] **Step 5: Run the e2e path test (defined in Task 6) — deferred**

This step is validated by Task 6. For now, run a quick scripted sanity check via MATLAB MCP on an existing protocol with an eigenmode head model active (see Task 6, Step 2). Expected: a cortex results node appears.

- [ ] **Step 6: Lint**

Run `/lint-matlab toolbox/process/functions/process_inverse_2018.m`. Expected: no new errors.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_inverse_2018.m
git commit -m "feat(inverse): dispatch eigenmode method in Compute() + coefficients node

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `process_eigenmodes_inverse` — thin wrapper over `Compute()`

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes_inverse.m`

- [ ] **Step 0: Add a whitening option to `GetDescription`**

In `process_eigenmodes_inverse('GetDescription')`, after the `nmodes` option (line 34) and before `outputtype`, add a checkbox so the batch/scripted path exposes the same toggle (default on):
```matlab
    sProcess.options.whiten.Comment = 'Apply noise whitening (recommended)';
    sProcess.options.whiten.Type    = 'checkbox';
    sProcess.options.whiten.Value    = 1;
```

- [ ] **Step 1: Rewrite `Run` to delegate to `Compute()`**

Keep the rest of `GetDescription` (Index 339, SubGroup 'Sources', the method/prior/snr/nmodes/outputtype options) and `FormatComment` as-is. Replace the entire `Run` function body (lines 46-164) with a wrapper that maps the process options into `OPTIONS` and calls the shared core:

```matlab
function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    % Map process options to the shared Compute() OPTIONS
    OPTIONS = struct();
    OPTIONS.InverseMethod    = 'eigenmode';
    OPTIONS.InverseMeasure   = lower(sProcess.options.method.Value);     % mne|dspm|sloreta
    OPTIONS.EigenmodePrior   = lower(sProcess.options.prior.Value);      % log|flat|power
    OPTIONS.SnrFixed         = sProcess.options.snr.Value{1};
    OPTIONS.SnrMethod        = 'fixed';
    OPTIONS.nModes           = sProcess.options.nmodes.Value{1};
    OPTIONS.NoiseMethod      = 'reg';
    OPTIONS.NoiseReg         = 0.1;
    OPTIONS.SourceOrient     = {'fixed'};
    OPTIONS.ComputeKernel    = 1;
    OPTIONS.DisplayMessages  = 0;
    OPTIONS.DataTypes        = [];     % auto-detect MEG else EEG in Compute
    OPTIONS.EigenmodeWhiten  = sProcess.options.whiten.Value;            % default 1
    OPTIONS.SaveCoefficients = ismember(lower(sProcess.options.outputtype.Value), {'coefficients','both'});
    % Resolve study/data indices for the shared core
    iStudies = [sInputs.iStudy];
    iDatas   = [sInputs.iItem];
    [OutputFiles, errMessage] = process_inverse_2018('Compute', iStudies, iDatas, OPTIONS);
    if ~isempty(errMessage)
        bst_report('Error', sProcess, sInputs, errMessage);
    end
end
```
Notes for the implementer:
- The `outputtype` option still offers `coefficients`/`sources`/`both`. With this mapping, `sources` and `both` always produce the cortex node (always created by `Compute`), and `coefficients`/`both` additionally set `SaveCoefficients`. A pure `coefficients` request still yields the cortex node too; if a coefficients-only output is required, add a follow-up to delete the cortex node — out of scope unless requested.
- Confirm `sInputs.iItem` is the correct data index expected by `Compute`'s `iDatas`. Cross-check against how `panel_protocols('TreeInverse', …)` builds `iDatas` (`panel_protocols.m:1355-1396`); if it uses `iItem`, this matches.
- Delete the now-unused helper code (old solve/reconstruct/node-building) entirely.

- [ ] **Step 2: Verify the process still registers and parses**

Run via MATLAB MCP:
```matlab
if ~brainstorm('status'); brainstorm nogui; end
sProc = process_eigenmodes_inverse('GetDescription');
disp(sProc.Comment); disp(sProc.Index);
```
Expected: prints `Eigenmode source mapping` and `339`, no error.

- [ ] **Step 3: Lint**

Run `/lint-matlab toolbox/process/functions/process_eigenmodes_inverse.m`. Expected: no new errors.

- [ ] **Step 4: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_inverse.m
git commit -m "refactor(inverse): make process_eigenmodes_inverse a thin Compute() wrapper

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: e2e test + manual verification

**Files:**
- Modify: `dev/tests/test_eigenmodes_inverse_e2e.m`

- [ ] **Step 1: Read the existing e2e test to learn its protocol setup**

Open `dev/tests/test_eigenmodes_inverse_e2e.m` and identify how it builds/loads a protocol with an eigenmode head model and an imported data file (the existing test drives `process_eigenmodes_inverse` directly). Reuse that setup verbatim.

- [ ] **Step 2: Add a subtest that drives the new `Compute()` path**

After the existing assertions, add a block that calls the shared core directly (DisplayMessages off so no dialog), exercising the interactive code path end-to-end:

```matlab
% ===== Interactive Compute() path: InverseMethod='eigenmode' =====
OPTIONS = struct();
OPTIONS.InverseMethod    = 'eigenmode';
OPTIONS.InverseMeasure   = 'dspm';
OPTIONS.EigenmodePrior   = 'log';
OPTIONS.SnrFixed         = 3;
OPTIONS.SnrMethod        = 'fixed';
OPTIONS.nModes           = 0;
OPTIONS.NoiseMethod      = 'reg';
OPTIONS.NoiseReg         = 0.1;
OPTIONS.SourceOrient     = {'fixed'};
OPTIONS.ComputeKernel    = 1;
OPTIONS.DisplayMessages  = 0;
OPTIONS.DataTypes        = [];
OPTIONS.SaveCoefficients = 1;

[outFiles, errMsg] = process_inverse_2018('Compute', iStudyTest, iDataTest, OPTIONS);
assert(isempty(errMsg), sprintf('Compute eigenmode path failed: %s', errMsg));
assert(~isempty(outFiles), 'Compute eigenmode path produced no output files.');

% A cortex results node must exist and be kernel-only with the right shape
iResNode = find(~cellfun(@isempty, regexp(outFiles, 'results_.*\.mat$', 'once')), 1);
assert(~isempty(iResNode), 'No cortex results node was created.');
ResCheck = in_bst_results(outFiles{iResNode}, 0, 'ImagingKernel', 'nComponents', 'SurfaceFile');
assert(~isempty(ResCheck.ImagingKernel), 'Cortex node must carry an ImagingKernel.');
assert(ResCheck.nComponents == 1, 'Eigenmode cortex node nComponents must be 1.');

% A coefficients matrix node must exist (SaveCoefficients=1, imported data)
iMatNode = find(~cellfun(@isempty, regexp(outFiles, 'matrix_eigencoeffs.*\.mat$', 'once')), 1);
assert(~isempty(iMatNode), 'Coefficients matrix node was not created.');
KWhiten = ResCheck.ImagingKernel;

% ===== Whitening OFF: pure sensor-space solve still produces a valid kernel =====
OPTIONS.EigenmodeWhiten  = 0;
OPTIONS.SaveCoefficients = 0;
[outFiles2, errMsg2] = process_inverse_2018('Compute', iStudyTest, iDataTest, OPTIONS);
assert(isempty(errMsg2), sprintf('Whiten-off eigenmode path failed: %s', errMsg2));
iRes2 = find(~cellfun(@isempty, regexp(outFiles2, 'results_.*\.mat$', 'once')), 1);
assert(~isempty(iRes2), 'Whiten-off path produced no cortex node.');
Res2 = in_bst_results(outFiles2{iRes2}, 0, 'ImagingKernel');
assert(all(isfinite(Res2.ImagingKernel(:))), 'Whiten-off kernel must be finite.');
assert(isequal(size(Res2.ImagingKernel), size(KWhiten)), 'Whiten on/off kernels must share shape.');
assert(max(abs(Res2.ImagingKernel(:) - KWhiten(:))) > 0, 'Whiten on/off should differ (whitener has an effect).');
```
Replace `iStudyTest`/`iDataTest` with the study/data indices the existing test already computes (match its variable names).

- [ ] **Step 3: Run the e2e test**

Run `dev/tests/test_eigenmodes_inverse_e2e.m` via the MATLAB MCP.
Expected: `ALL TESTS PASSED` (or the test's existing success marker).

- [ ] **Step 4: Manual GUI verification**

In the Brainstorm GUI (launch via `/brainstorm-start`):
1. Open a protocol whose active head model is an eigenmode leadfield (`isEigenmode=1`).
2. Right-click a recording → **Compute sources [2018]**.
3. Confirm: only **Eigenmode source mapping** is enabled in the Method panel (MN/LCMV/Dipoles/MEM disabled), the Eigenmode measure + spectral-prior sub-panel shows, and the "Also save eigenmode coefficients" checkbox is present.
4. Run with dSPM + Log, coefficients checked.
5. Confirm a cortex source node appears and opens (3D map + time slider), and an `EigenCoeffs …` matrix node appears with per-mode time series.
6. Right-click a recording in a **non-eigenmode** protocol → Compute sources → confirm the eigenmode option is hidden and standard methods work unchanged.

- [ ] **Step 5: Commit**

```bash
git add dev/tests/test_eigenmodes_inverse_e2e.m
git commit -m "test(inverse): cover interactive eigenmode Compute() path e2e

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] Run all new/changed pure tests in sequence via MATLAB MCP:
  - `dev/tests/test_noise_whitener_pure.m` → `ALL TESTS PASSED`
  - `dev/tests/test_eigenmode_reconstruct_pure.m` → `ALL TESTS PASSED`
  - `dev/tests/test_inverse_eigenmodes_pure.m` (regression, unchanged) → `ALL TESTS PASSED`
  - `dev/tests/test_eigenmodes_inverse_e2e.m` → success marker
- [ ] Confirm `git diff development -- toolbox/inverse/bst_inverse_linear_2018.m` is **empty** (proven solver untouched).
- [ ] Manual GUI checks in Task 6 Step 4 complete.
- [ ] Use `superpowers:finishing-a-development-branch` to decide merge/PR.
```
