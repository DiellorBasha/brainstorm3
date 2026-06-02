# Harmonic Eigenmode Inverse — Phase A (Compute the node)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute an unregularized, whitened eigenmode ("Harmonic") source-imaging kernel and save it as a results node — both interactively (Compute sources dialog) and in batch (process) — carrying the eigenmode-space operator `EigenKernel` so a later phase can drive a consistent time series.

**Architecture:** A pure helper `bst_eigenmodes_harmonic(L,Φ,iW)` returns the eigenmode-space kernel `M̃ = pinv(iW·L·Φ)·iW` (reusing the rank-safe `bst_eigenmodes_transform`). `bst_inverse_eigenmodes` gains a `Method='harmonic'` branch that calls it. Both compute surfaces — the existing `process_eigenmodes_inverse` (batch) and a new `Harmonic` method in the `panel_inverse_2018`/`process_inverse_2018` dialog (interactive) — save a kernel-only results node with `ImagingKernel = Φ·M̃` (vertex-space, displayable), `EigenKernel = M̃`, and `Function = 'eigenmode_harmonic'`.

**Tech Stack:** MATLAB; Brainstorm inverse pipeline (`bst_inverse_eigenmodes`, `bst_inverse_linear_2018`, `process_inverse_2018`, `panel_inverse_2018`), `bst_eigenmodes_transform`, `in_tess_eigenmodes`, `db_template('resultsmat')`.

**Scope note:** Phase B (re-point `view_eigenmodes_timeseries` to read the node + DB-tree launch + remove the figure-popup launch) is a **separate follow-up plan** — it depends on these nodes existing.

---

## Reference facts (already verified in the codebase)

- `[Kernel,Info] = bst_eigenmodes_transform(Gain, Phi)` computes `L̃ = Gain*Phi` then the rank-safe SVD pseudoinverse `pinv(L̃)` `[K×nCh]` (`toolbox/math/bst_eigenmodes_transform.m`). It does **not** whiten.
- `bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, 'Method',...,'nModes',...,'GoodChannel',...)` loads the constrained gain, eigenmodes, builds the whitener `iW`, and (for mne/dspm/sloreta) returns `Results.ImagingKernel` `[K×nGoodCh]` plus `Whitener`, `nModes`, `GoodChannel`, `Eigenvalues`, `SurfaceFile`, `ConditionNumber*`. The method `switch` is the `case 'mne'/'dspm'/'sloreta'` block; methods are validated against `{'mne','dspm','sloreta'}` near the top.
- `process_eigenmodes_inverse.m`: method dropdown at `sProcess.options.method` (radio_linelabel with values `'mne','dspm','sloreta'`); the `'sources'` save path builds `VertexKernel = Phi * InvResults.ImagingKernel`, sets `ResMat.Function = Method`, and saves via `db_add_data`.
- `panel_inverse_2018.CreatePanel`: method radios created at lines ~140–146 (`jRadioMethodMn/Bf/Dip/Mem`); `GetSelectedMethod(ctrl)` maps radios → `(Method,Measure)`; `GetMethodComment(Method,Measure)` returns a short tag; `UpdatePanel(...)` toggles option-group visibility with `setVisible(...)`; `GetPanelContents()` returns `s.InverseMethod` etc.
- `process_inverse_2018.Compute` (line ~666): `switch(OPTIONS.InverseMethod)` with `case {'minnorm','gls','lcmv'}` → `bst_inverse_linear_2018`, `case 'mem'` → `be_main`. Available in scope there: `HeadModelFile`, `HeadModelInit` (`.SurfaceFile`, `.Gain` unconstrained), `GoodChannel` (index vector), `ChannelMat`, `sStudyChannel` (for `NoiseCov(1).FileName`), `nSources`. After the switch: `ResultsMat = db_template('resultsmat'); ResultsMat = struct_copy_fields(ResultsMat, Results, 1);` then fields assembled and `bst_save`.
- `db_template('resultsmat')` has `ImagingKernel, ImageGridAmp, Whitener, nComponents, Function, DataFile, HeadModelFile, HeadModelType, GoodChannel, SurfaceFile, ...` — **no `EigenKernel`** (added explicitly by callers).

---

## File structure

| File | Responsibility | Action |
|------|----------------|--------|
| `toolbox/math/bst_eigenmodes_harmonic.m` | Pure: `M̃ = pinv(iW·L·Φ)·iW` | **Create** |
| `toolbox/inverse/bst_inverse_eigenmodes.m` | `Method='harmonic'` branch | Modify |
| `toolbox/process/functions/process_eigenmodes_inverse.m` | batch `harmonic` option + `EigenKernel`/`Function` on node | Modify |
| `toolbox/inverse/panel_inverse_2018.m` | "Harmonic (eigenmodes)" method radio + gating | Modify |
| `toolbox/process/functions/process_inverse_2018.m` | `case 'harmonic'` → delegate + `EigenKernel` | Modify |
| `dev/tests/test_eigenmodes_harmonic_pure.m` | pure test for the harmonic kernel | **Create** |
| `dev/tests/test_harmonic_inverse_e2e.m` | guarded smoke: compute node, check fields | **Create** |

---

## Task 1: Pure harmonic-kernel helper

**Files:**
- Create: `toolbox/math/bst_eigenmodes_harmonic.m`
- Test: `dev/tests/test_eigenmodes_harmonic_pure.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_harmonic_pure.m`:

```matlab
function test_eigenmodes_harmonic_pure
% Verify the harmonic eigenmode kernel M = pinv(iW*L*Phi)*iW:
%   - reconstructs the whitened compressed system (Phi*M is the min-norm-in-mode-space inverse)
%   - is rank-safe when L*Phi is rank-deficient (no blow-up)
%   - equals bst_eigenmodes_transform(iW*L,Phi)*iW
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Synthetic problem: 8 channels, 12 vertices, 5 modes (deterministic, no rng)
nCh = 8; nV = 12; K = 5;
L   = reshape(1:(nCh*nV), nCh, nV) / 100;        % deterministic [nCh x nV]
Phi = zeros(nV, K);
for k = 1:K, Phi(:,k) = cos((1:nV)' * k / 3); end % deterministic [nV x K]
iW  = diag(linspace(1, 2, nCh));                  % deterministic whitener [nCh x nCh]

M = bst_eigenmodes_harmonic(L, Phi, iW);          % [K x nCh]
assert(isequal(size(M), [K, nCh]), 'M must be [K x nCh].');

% Equals the explicit construction via the transform
[Kt, ~] = bst_eigenmodes_transform(iW*L, Phi);    % pinv(iW*L*Phi) [K x nCh]
Mref = Kt * iW;
assert(max(abs(M(:) - Mref(:))) < 1e-9, 'M must equal pinv(iW*L*Phi)*iW.');

% M maps RAW sensor data to coefficients, so M*(L*Phi) ~ I_K (full column rank):
%   M*(L*Phi) = pinv(iW*L*Phi)*iW*(L*Phi) = pinv(iW*L*Phi)*(iW*L*Phi) = I_K
ID = M * (L * Phi);                                % should be ~ I_K
assert(max(abs(ID(:) - reshape(eye(K),[],1))) < 1e-6, 'M*(L*Phi) must be ~ I_K for full-rank case.');

% Rank-safety: duplicate two mode columns so L*Phi is rank-deficient -> no Inf/NaN
Phi2 = Phi; Phi2(:,5) = Phi2(:,1);
M2 = bst_eigenmodes_harmonic(L, Phi2, iW);
assert(all(isfinite(M2(:))), 'Rank-deficient input must not produce Inf/NaN (rank-safe pinv).');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run the test, confirm it FAILS** with `Undefined function 'bst_eigenmodes_harmonic'`.

Run (MATLAB MCP `run_matlab_file` on the test path), or:
```matlab
run('dev/tests/test_eigenmodes_harmonic_pure.m')
```

- [ ] **Step 3: Create `toolbox/math/bst_eigenmodes_harmonic.m`:**

```matlab
function M = bst_eigenmodes_harmonic(L, Phi, iW)
% BST_EIGENMODES_HARMONIC: Unregularized whitened eigenmode imaging kernel.
%
% USAGE:  M = bst_eigenmodes_harmonic(L, Phi, iW)
%
% DESCRIPTION:
%     Returns the eigenmode-space "Harmonic" kernel
%         M = pinv(iW * L * Phi) * iW            [K x nCh]
%     i.e. the rank-safe (SVD) pseudoinverse of the whitened compressed lead
%     field, applied to the whitener so it maps RAW recordings to eigenmode
%     coefficients. There is no source/inverse regularization (no Tikhonov, no
%     prior) -- only the rank-safe singular-value floor of bst_eigenmodes_transform.
%     The vertex-space source map is Phi*M; the eigenmode time series is M*Data.
%
% INPUTS:
%     L   : [nCh x nVert] constrained (fixed-orientation) lead field, good channels.
%     Phi : [nVert x K]   eigenmode matrix (caller truncates to K modes).
%     iW  : [nCh x nCh]   noise whitener (pass eye(nCh) for no whitening).
%
% OUTPUT:
%     M   : [K x nCh] harmonic eigenmode kernel.
%
% SEE ALSO: bst_eigenmodes_transform, bst_inverse_eigenmodes

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

% Rank-safe pseudoinverse of the whitened compressed lead field, then de-whiten
% the data side so the kernel applies to raw recordings.
[Kt, ~] = bst_eigenmodes_transform(iW * L, Phi);   % pinv(iW*L*Phi)  [K x nCh]
M = Kt * iW;                                        % [K x nCh]
end
```

- [ ] **Step 4: Run the test, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_harmonic.m dev/tests/test_eigenmodes_harmonic_pure.m
git commit -m "Harmonic inverse: pure whitened eigenmode kernel helper"
```

---

## Task 2: `harmonic` method in `bst_inverse_eigenmodes`

**Files:**
- Modify: `toolbox/inverse/bst_inverse_eigenmodes.m`

- [ ] **Step 1: Allow the method.** Find the validation line near the top:
```matlab
if ~ismember(Method, {'mne', 'dspm', 'sloreta'})
    errMsg = ['Unknown method: ' Method '. Use mne, dspm, or sloreta.'];
    return;
end
```
Replace with:
```matlab
if ~ismember(Method, {'mne', 'dspm', 'sloreta', 'harmonic'})
    errMsg = ['Unknown method: ' Method '. Use mne, dspm, sloreta, or harmonic.'];
    return;
end
```

- [ ] **Step 2: Add the `harmonic` branch to the kernel `switch`.** Find the method switch (the `switch Method` with `case 'mne' / 'dspm' / 'sloreta'`). Add a new case at the top of it:
```matlab
    case 'harmonic'
        % Unregularized, whitened eigenmode reconstruction (no SNR/prior).
        % iW is the whitener built above; L_tilde = L*Phi is already computed.
        Kernel = bst_eigenmodes_harmonic(L, Phi, iW);   % [K x nChannels]
        EigenGains = ones(K, 1);
```
This uses the in-scope `L` (good-channel constrained gain), `Phi`, `iW`, and `K` already defined earlier in the function. (For `harmonic`, the earlier `SourcePrior`/`L_ws`/SVD/`Lambda` computations are still harmless to run; the branch ignores them.)

- [ ] **Step 3: Verify it parses and runs the existing transform path.** (MATLAB MCP `check_matlab_code` on the file → no new issues; `clear functions; which('bst_inverse_eigenmodes')` → path.) Re-run Task 1's pure test (helper unaffected) → `ALL TESTS PASSED`.

- [ ] **Step 4: Commit**
```bash
git add toolbox/inverse/bst_inverse_eigenmodes.m
git commit -m "Harmonic inverse: add Method='harmonic' to bst_inverse_eigenmodes"
```

---

## Task 3: Batch — `harmonic` option in `process_eigenmodes_inverse`

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes_inverse.m`

- [ ] **Step 1: Add `harmonic` to the method dropdown.** Replace the `sProcess.options.method` block:
```matlab
    sProcess.options.method.Comment = {'MNE (minimum norm)', 'dSPM (noise-normalized)', ...
                                       'sLORETA (standardized)', ...
                                       'Inverse method:'; ...
                                       'mne', 'dspm', 'sloreta', ''};
    sProcess.options.method.Type    = 'radio_linelabel';
    sProcess.options.method.Value   = 'dspm';
```
with:
```matlab
    sProcess.options.method.Comment = {'MNE (minimum norm)', 'dSPM (noise-normalized)', ...
                                       'sLORETA (standardized)', 'Harmonic (unregularized)', ...
                                       'Inverse method:'; ...
                                       'mne', 'dspm', 'sloreta', 'harmonic', ''};
    sProcess.options.method.Type    = 'radio_linelabel';
    sProcess.options.method.Value   = 'dspm';
```

- [ ] **Step 2: Tag the saved sources node as Harmonic and carry `EigenKernel`.** In `Run`, the `'sources'` save block builds `VertexKernel = Phi * InvResults.ImagingKernel` and a `ResMat`. Locate the assignment:
```matlab
            ResMat.Function      = Method;
```
Immediately after it, add:
```matlab
            if strcmpi(Method, 'harmonic')
                ResMat.Function   = 'eigenmode_harmonic';
                ResMat.EigenKernel = InvResults.ImagingKernel;   % M̃ [K x nGoodChannels], for the time series
            end
```
(There are two `ResMat.Function = Method;` style assignments if both raw and imported branches save a sources node — apply the same three lines after **each** `ResMat.Function`/`ResMat.Function      =` assignment that precedes a sources-node `bst_save`. Search the file for `.Function` to find them.)

- [ ] **Step 3: Fix the comment for harmonic** so it doesn't print "alpha". In `FormatComment` and in the `ResMat.Comment = sprintf('EigenInv %s ...')` lines, harmonic has no alpha/SNR. Minimal change: after computing `alphaStr`, force it empty for harmonic. Find:
```matlab
        alphaStr = '';
        if PriorAlpha > 0
            alphaStr = sprintf(', a=%.1f', PriorAlpha);
        end
```
Replace with:
```matlab
        alphaStr = '';
        if PriorAlpha > 0 && ~strcmpi(Method, 'harmonic')
            alphaStr = sprintf(', a=%.1f', PriorAlpha);
        end
```

- [ ] **Step 4: Verify.** (MATLAB MCP `check_matlab_code` on the file → no new issues; `clear functions; which('process_eigenmodes_inverse')` → path.)

- [ ] **Step 5: Commit**
```bash
git add toolbox/process/functions/process_eigenmodes_inverse.m
git commit -m "Harmonic inverse: batch process option saves eigenmode_harmonic node with EigenKernel"
```

---

## Task 4: Interactive — "Harmonic (eigenmodes)" method radio

**Files:**
- Modify: `toolbox/inverse/panel_inverse_2018.m`

- [ ] **Step 1: Add the radio.** In `CreatePanel`, after the `jRadioMethodDip` line and before the MEM block (the `if ~isProcess ... jRadioMethodMem ...`), insert:
```matlab
        jRadioMethodHarm = gui_component('radio', jPanelMethod, 'br', 'Harmonic (eigenmodes)', jGroupMethod, '', @Method_Callback, []);
```
Then, in the same scope where MEM is disabled for non-surface (`if ~isempty(jRadioMethodMem) && ~isProcess && (~strcmpi(HeadModelType, 'surface') || isShared)`), add an equivalent guard so Harmonic is **surface-only**:
```matlab
        % Harmonic (eigenmode) inverse requires a surface head model
        if ~isProcess && ~strcmpi(HeadModelType, 'surface')
            jRadioMethodHarm.setEnabled(0);
        end
```
And in the default-selection `switch lower(OPTIONS.InverseMethod)`, add:
```matlab
            case 'harmonic', jRadioMethodHarm.setSelected(1);
```

- [ ] **Step 2: Register the handle in the `ctrl` struct.** In the struct returned to `gui_river`/`bst_get('PanelControls')` (the big `ctrl = struct('jRadioMethodMn', jRadioMethodMn, ...)` near line ~374), add:
```matlab
            'jRadioMethodHarm', jRadioMethodHarm, ...
```

- [ ] **Step 3: Map the radio in `GetSelectedMethod`.** Add a branch:
```matlab
    elseif isfield(ctrl, 'jRadioMethodHarm') && ~isempty(ctrl.jRadioMethodHarm) && ctrl.jRadioMethodHarm.isSelected()
        Method  = 'harmonic';
        Measure = 'amplitude';
```
(Place it among the existing `elseif ctrl.jRadioMethodDip... / jRadioMethodBf...` branches.)

- [ ] **Step 4: Give it a tag in `GetMethodComment`.** In the `switch (lower(Method))`, add:
```matlab
        case 'harmonic'
            Comment = 'Harmonic';
```

- [ ] **Step 5: Gate the option groups in `UpdatePanel`.** The visibility block uses `isLinear` + which radio is selected. Add a local `isHarm = ctrl.jRadioMethodHarm.isSelected();` near the top of `UpdatePanel` (after `ctrl` is fetched), and AND `~isHarm` into the `setVisible(...)` calls for the groups Harmonic must not show: Min-norm measure (`jPanelMeasureMN`), Beamformer measure (`jPanelMeasureBf`), SNR (`jPanelSnr`), and Depth weighting (`jPanelDepth`). For example:
```matlab
        jPanelMeasureMN.setVisible(isLinear && ctrl.jRadioMethodMn.isSelected() && ~isHarm);
        jPanelSnr.setVisible(isLinear && ~ctrl.jRadioMethodDip.isSelected() && ~ctrl.jRadioMethodBf.isSelected() && ~isHarm);
        jPanelDepth.setVisible(isLinear && ctrl.jRadioMethodMn.isSelected() && ~ctrl.jRadioMnSloreta.isSelected() && ~isHarm);
```
The Output (kernel/full), Sensors, Noise-cov-regularization, and Source-model panels stay visible (Harmonic still needs Output=kernel and the channel selection; the orientation is forced constrained — see Task 5). Keep the existing logic for all other methods intact; only add the `~isHarm` conjuncts.

- [ ] **Step 6: Verify it parses.** (MATLAB MCP `check_matlab_code` on `panel_inverse_2018.m` → no NEW issues from the added lines; `clear functions; which('panel_inverse_2018')` → path.) Full GUI behavior is validated manually + by Task 6's e2e (which calls `process_inverse_2018.Compute` directly, bypassing the dialog).

- [ ] **Step 7: Commit**
```bash
git add toolbox/inverse/panel_inverse_2018.m
git commit -m "Harmonic inverse: add Harmonic method to the Compute sources dialog"
```

---

## Task 5: Interactive — `case 'harmonic'` in `process_inverse_2018.Compute`

**Files:**
- Modify: `toolbox/process/functions/process_inverse_2018.m`

- [ ] **Step 1: Add the method branch.** In `Compute`, the `switch( OPTIONS.InverseMethod )` (line ~666) has `case {'minnorm','gls','lcmv'}` and `case 'mem'`. Add before `otherwise`:
```matlab
            case 'harmonic'
                % Unregularized whitened eigenmode (Harmonic) inverse.
                SurfaceFile  = HeadModelInit.SurfaceFile;
                NoiseCovFile = '';
                if ~isempty(sStudyChannel.NoiseCov) && ~isempty(sStudyChannel.NoiseCov(1).FileName)
                    NoiseCovFile = sStudyChannel.NoiseCov(1).FileName;
                end
                % Logical good-channel mask over ALL channels (bst_inverse_eigenmodes expects this)
                nAllChan = length(ChannelMat.Channel);
                GoodMask = false(nAllChan, 1);
                GoodMask(GoodChannel) = true;
                [InvE, errE] = bst_inverse_eigenmodes(HeadModelFile, SurfaceFile, NoiseCovFile, ...
                    'Method', 'harmonic', 'nModes', OPTIONS.nModes, 'GoodChannel', GoodMask);
                if ~isempty(errE)
                    errMessage = [errMessage errE 10];
                    break;
                end
                % Vertex-space kernel for cortex display
                [Eig, ~] = in_tess_eigenmodes(SurfaceFile);
                Kuse = InvE.nModes;
                PhiUse = double(Eig.Vectors(:, 1:Kuse));
                Results = struct();
                Results.ImagingKernel = PhiUse * InvE.ImagingKernel;   % [nVert x nGoodCh]
                Results.ImageGridAmp  = [];
                Results.nComponents   = 1;
                Results.Whitener      = InvE.Whitener;
                Results.EigenKernel   = InvE.ImagingKernel;            % M̃ [K x nGoodCh]
                OPTIONS.FunctionName  = 'eigenmode_harmonic';
```

- [ ] **Step 2: Carry `EigenKernel` past `struct_copy_fields`.** `db_template('resultsmat')` has no `EigenKernel`, and `struct_copy_fields` may drop unknown fields. After the existing line:
```matlab
        ResultsMat = struct_copy_fields(ResultsMat, Results, 1);
```
add:
```matlab
        if strcmpi(OPTIONS.InverseMethod, 'harmonic') && isfield(Results, 'EigenKernel')
            ResultsMat.EigenKernel = Results.EigenKernel;
        end
```

- [ ] **Step 3: Ensure `OPTIONS.nModes` exists.** The dialog/process options may not define `nModes`. Near the top of `Compute`, where `OPTIONS` defaults are set (the `Def_OPTIONS`/merge area ~line 135), add a default after the struct is built:
```matlab
    if ~isfield(OPTIONS, 'nModes') || isempty(OPTIONS.nModes)
        OPTIONS.nModes = 0;   % 0 = auto (min of channels and available modes)
    end
```
(And in Task 4, `GetPanelContents` should set `s.nModes` from the Harmonic `nModes` field if you add one; `0`/auto is an acceptable default for this phase, so a dedicated field is optional.)

- [ ] **Step 4: Verify.** (MATLAB MCP `check_matlab_code` on `process_inverse_2018.m` → no new issues; `clear functions; which('process_inverse_2018')` → path.)

- [ ] **Step 5: Commit**
```bash
git add toolbox/process/functions/process_inverse_2018.m
git commit -m "Harmonic inverse: compute branch delegates to bst_inverse_eigenmodes + saves EigenKernel"
```

---

## Task 6: End-to-end smoke (guarded)

**Files:**
- Create: `dev/tests/test_harmonic_inverse_e2e.m`

- [ ] **Step 1: Write the guarded e2e harness** `dev/tests/test_harmonic_inverse_e2e.m`:

```matlab
function test_harmonic_inverse_e2e
% Smoke: compute a Harmonic inverse node and verify its consistency contract:
%   - node.Function == 'eigenmode_harmonic'
%   - node.EigenKernel is [K x nGoodCh]
%   - node.ImagingKernel == Phi * node.EigenKernel  (vertex map = Phi * coefficients)
% Requires a protocol with a surface head model + eigenmodes + noise cov + a data file.
% Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol, 'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.');
    return;
end
% Find a study with head model (surface) + eigenmodes + noise cov + data
target = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && (s.iHeadModel >= 1) ...
            && (length(s.HeadModel) >= s.iHeadModel) && isfield(s,'Data') && ~isempty(s.Data) ...
            && isfield(s,'NoiseCov') && ~isempty(s.NoiseCov) && ~isempty(s.NoiseCov(1).FileName)
        try
            hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType', 'SurfaceFile');
            if strcmpi(hm.HeadModelType, 'surface')
                [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
                if isEig
                    target = struct('study', s, 'surf', hm.SurfaceFile, ...
                        'hmFile', s.HeadModel(s.iHeadModel).FileName, ...
                        'ncFile', s.NoiseCov(1).FileName);
                    break;
                end
            end
        catch
        end
    end
end
if isempty(target)
    disp('SKIP: no study with surface head model + eigenmodes + noise cov.');
    return;
end

% Compute the harmonic kernel directly via bst_inverse_eigenmodes (engine under both surfaces)
[InvE, errE] = bst_inverse_eigenmodes(target.hmFile, target.surf, target.ncFile, ...
    'Method', 'harmonic', 'nModes', 0);
assert(isempty(errE), ['bst_inverse_eigenmodes harmonic failed: ' errE]);
assert(~isempty(InvE.ImagingKernel), 'Harmonic kernel must be non-empty.');
K = InvE.nModes;
assert(size(InvE.ImagingKernel, 1) == K, 'EigenKernel must have K rows.');

% Vertex map = Phi * M̃
[Eig, ~] = in_tess_eigenmodes(target.surf);
Phi = double(Eig.Vectors(:, 1:K));
Vmap = Phi * InvE.ImagingKernel;
assert(size(Vmap, 1) == size(Eig.Vectors, 1), 'Vertex kernel must have nVert rows.');
assert(all(isfinite(InvE.ImagingKernel(:))), 'Harmonic kernel must be finite (rank-safe).');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it** → expect `ALL TESTS PASSED` (if suitable data) or a `SKIP:` line. Must not error.

- [ ] **Step 3: Commit**
```bash
git add dev/tests/test_harmonic_inverse_e2e.m
git commit -m "Harmonic inverse: end-to-end smoke (kernel + Phi*M̃ consistency)"
```

---

## Self-Review Notes

- **Spec coverage:** unit 1 (`harmonic` method) → Tasks 1–2; unit 2 (node fields `ImagingKernel=Φ·M̃`, `EigenKernel=M̃`, `Function`) → Tasks 3 & 5; unit 3 (interactive dialog) → Tasks 4–5; unit 4 (batch) → Task 3. Units 5–6 (viewer re-point, tree launch) are **Phase B** (separate plan), as flagged.
- **`EigenKernel` field:** added explicitly by each caller (batch Task 3, interactive Task 5 Step 2) because `db_template('resultsmat')` lacks it — no core template change.
- **Whitener:** Harmonic uses `bst_inverse_eigenmodes`' own whitener from the noise-cov file; the dialog's noise/SNR controls are inert for it (gated off in Task 4). If no noise cov, `iW=I` (unwhitened) — `bst_inverse_eigenmodes` already handles this.
- **GUI caveat:** Task 4's Swing edits follow the file's existing radio/`UpdatePanel` patterns; the implementer must match the exact surrounding text (the spec gives anchors, not full-file context). Task 6's e2e validates the math via the engine directly, independent of the dialog.
