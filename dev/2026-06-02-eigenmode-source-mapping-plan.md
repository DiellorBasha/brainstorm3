# Eigenmode Source Mapping Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replicate the GBF eigenmode source-mapping method in Brainstorm by splitting it into a standalone forward composer (`L̃ = L·Φ`) and a dedicated mode-space inverse whose spatial prior is the eigenvalue spectrum used as the source covariance `R`.

**Architecture:** Three clean stages along Brainstorm's existing forward/inverse seam. STAGE 1 (`toolbox/forward/`) composes any base head model with the surface eigenmodes into a `headmodel_eigenmode_*.mat` node (`Gain = L·Φ`). STAGE 2 (`toolbox/inverse/`) consumes that composed leadfield, applies Brainstorm's standard data cleaning (bad channels, SSP projectors, noise whitening), builds the spectral prior `R` from the eigenvalues, and solves the regularized MAP estimate. STAGE 3 emits a standard coefficient matrix node and a kernel-only cortex results node. The old fused `bst_inverse_eigenmodes`, the harmonic/transform paths, and the bespoke time-series viewer are retired; the eigenspectrum tools are untouched.

**Tech Stack:** MATLAB / Brainstorm. Tests are MATLAB functions under `dev/tests/` following the existing `test_*_pure.m` (synthetic, deterministic, no DB) and `test_*_e2e.m` (skips cleanly without a protocol) convention, run via the MATLAB MCP (`run_matlab_test_file`).

**Spec:** `dev/2026-06-02-eigenmode-source-mapping-design.md`

---

## File structure

| File | Action | Responsibility |
|---|---|---|
| `toolbox/forward/bst_eigenmode_leadfield.m` | Create | Engine: base `HeadModel` struct + `Eigenmodes` → composed head model struct (`Gain = L·Φ`, constrained). |
| `toolbox/process/functions/process_eigenmode_leadfield.m` | Create | Process: pick base head model + surface eigenmodes, write `headmodel_eigenmode_*.mat` node. |
| `toolbox/math/bst_eigenmode_prior.m` | Create | Build the diagonal source prior `R` from eigenvalues (log / flat / power). |
| `toolbox/inverse/bst_inverse_eigenmodes.m` | Rewrite | Consume composed leadfield + cleaning + prior + regularized solve → `M̃`. |
| `toolbox/process/functions/process_eigenmodes_inverse.m` | Rewrite | Run the inverse; emit coefficient matrix + cortex results nodes. |
| `toolbox/inverse/bst_resolution_metrics.m` | Create | Resolution-matrix point-spread metrics (validation Level 1). |
| `toolbox/script/tutorial_eigenmodes_validation.m` | Create | Validation harness: Levels 1–3, writes a results report. |
| `toolbox/math/bst_eigenmodes_harmonic.m` | Delete | Retired (flat-`R` + unreg switch replaces it). |
| `toolbox/math/bst_eigenmodes_transform.m` | **Keep** | Shared eigenspectrum infra: `process_eigenmodes_transform`/`_denoise` depend on it for PSD coefficient building. (Plan originally listed it for deletion — corrected during execution.) |
| `toolbox/gui/view_eigenmodes_timeseries.m` | Delete | Retired (coefficients are a standard matrix node). |
| `toolbox/inverse/panel_inverse_2018.m`, `toolbox/tree/tree_callbacks.m`, `toolbox/core/bst_figures.m` | Edit | Remove the "Harmonic (eigenmodes)" option/menu/figure hooks wired into the standard inverse GUI. |
| `dev/tests/test_eigenmodes_harmonic_pure.m`, `test_harmonic_inverse_e2e.m`, `test_eigenmodes_transform_pure.m`, `test_view_eigenmodes_timeseries_pure.m`, `test_eigenmode_timeseries_e2e.m`, `test_eigenmode_viewer_*.m` | Delete | Tests for retired units. |

**Reused helpers (do not modify):** `bst_gain_orient(Gain, GridOrient, GridAtlas)`, `in_tess_eigenmodes(SurfaceFile)` → `(Eigenmodes{.Vectors[nV×nModes], .Values[nModes×1], .nModes}, isComputed)`, `bst_whitener(NoiseCov, ChannelFile, DataTypes, ChannelFlag)`, `tess_laplacian(Vertices, Faces)` → `[L, M]`, `in_bst_headmodel`, `in_bst_channel`, `db_template`, `db_add_data`, `good_channel`.

---

## Phase 1 — Forward composer

### Task 1: `bst_eigenmode_leadfield` engine

**Files:**
- Create: `toolbox/forward/bst_eigenmode_leadfield.m`
- Test: `dev/tests/test_eigenmode_leadfield_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_eigenmode_leadfield_pure
% Verify the forward composer:
%   - composed Gain == constrained(L) * Phi(:,1:K), exactly
%   - constrained extraction matches bst_gain_orient on unconstrained gain
%   - K is clamped to available modes; eigenvalues carried through
%   - output is a valid composed head-model struct (isEigenmode, SurfaceFile, nModes)
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Synthetic: 6 channels, 5 vertices, unconstrained gain (3 cols/vertex), 4 modes.
nCh = 6; nV = 5; K = 4;
GainU = zeros(nCh, 3*nV);
for i = 1:nCh
    for j = 1:3*nV
        GainU(i,j) = sin(i * j * pi / 17);     % deterministic, full rank
    end
end
% Surface normals (unit) per vertex -> constrained leadfield via bst_gain_orient
GridOrient = zeros(nV,3);
for v = 1:nV
    n = [cos(v), sin(v), 0.5];
    GridOrient(v,:) = n / norm(n);
end
Lc = bst_gain_orient(GainU, GridOrient);        % [nCh x nV] constrained reference

% Eigenmodes struct (deterministic DCT-II basis, M-orthonormal not required for shape test)
Vectors = zeros(nV, K+1);                        % one extra mode to test clamping
for k = 1:(K+1)
    for n = 1:nV
        Vectors(n,k) = cos(pi/nV * (n-0.5) * (k-1));
    end
end
Values = (0:K)';                                 % ascending eigenvalues, includes DC=0
Eig = struct('Vectors', Vectors, 'Values', Values, 'nModes', K+1);

% Base head-model struct mimicking in_bst_headmodel output
HeadModel = struct('Gain', GainU, 'GridOrient', GridOrient, ...
    'GridLoc', zeros(nV,3), 'GridAtlas', [], ...
    'HeadModelType', 'surface', 'SurfaceFile', 'tess_cortex_test.mat', ...
    'Comment', 'OS-MEG test');

CompHM = bst_eigenmode_leadfield(HeadModel, Eig, 'nModes', K);

% --- assertions ---
assert(isequal(size(CompHM.Gain), [nCh, K]), 'Composed Gain must be [nCh x K].');
Lref = Lc * Vectors(:,1:K);
assert(max(abs(CompHM.Gain(:) - Lref(:))) < 1e-9, 'Composed Gain must equal constrained(L)*Phi.');
assert(CompHM.nModes == K, 'nModes must be clamped to K.');
assert(isequal(CompHM.Eigenvalues(:), Values(1:K)), 'Eigenvalues must be carried through (first K).');
assert(isfield(CompHM,'isEigenmode') && CompHM.isEigenmode == 1, 'isEigenmode flag must be set.');
assert(strcmp(CompHM.SurfaceFile, 'tess_cortex_test.mat'), 'SurfaceFile must be carried through.');
assert(isempty(CompHM.GridLoc) && isempty(CompHM.GridOrient), 'GridLoc/GridOrient must be empty on composed HM.');

% Clamp test: asking for more modes than available returns all available
CompHM2 = bst_eigenmode_leadfield(HeadModel, Eig, 'nModes', 999);
assert(CompHM2.nModes == (K+1), 'nModes must clamp to available count.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `run_matlab_test_file`: `dev/tests/test_eigenmode_leadfield_pure.m`
Expected: FAIL — `Undefined function 'bst_eigenmode_leadfield'`.

- [ ] **Step 3: Write minimal implementation**

```matlab
function CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, varargin)
% BST_EIGENMODE_LEADFIELD: Compose a leadfield into the LBO eigenmode basis.
%
% USAGE:  CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, 'nModes', K)
%
% DESCRIPTION:
%     Forward solution for eigenmode-space source mapping (GBF). Takes a base
%     surface head model and the precomputed surface eigenmodes and returns a
%     composed head-model struct whose Gain is the eigenmode leadfield
%         L̃ = L · Φ            [nChannels × K]
%     where L is the constrained (surface-normal) leadfield and Φ = Eigenmodes
%     truncated to K modes. Each column of L̃ is the sensor topography of one
%     eigenmode. This is strictly a forward operation; the inverse is separate.
%
% INPUTS:
%     HeadModel  : base head-model struct (from in_bst_headmodel, ApplyOrient=0):
%                  .Gain [nCh × 3*nVert] unconstrained, .GridOrient [nVert×3],
%                  .GridAtlas, .HeadModelType, .SurfaceFile, .Comment
%     Eigenmodes : struct from in_tess_eigenmodes: .Vectors [nVert×nModes],
%                  .Values [nModes×1], .nModes
% OPTIONS:
%     'nModes' : number of leading modes to keep (default: all; clamped to available)
%
% OUTPUT:
%     CompHM : composed head-model struct ready to save as headmodel_eigenmode_*.mat
%
% Authors: Diellor Basha, 2026

nModes = [];
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nmodes', nModes = varargin{i+1};
    end
end

% Constrained leadfield: [nCh × nVert]
Lc = bst_gain_orient(double(HeadModel.Gain), HeadModel.GridOrient, ...
    getfield_default(HeadModel, 'GridAtlas', []));

Phi    = double(Eigenmodes.Vectors);     % [nVert × nModesAll]
Values = double(Eigenmodes.Values(:));   % [nModesAll × 1]
nVert  = size(Phi, 1);
nModesAll = size(Phi, 2);

if size(Lc,2) ~= nVert
    error('bst_eigenmode_leadfield:VertexMismatch', ...
        'Leadfield has %d sources but eigenmodes have %d vertices.', size(Lc,2), nVert);
end

% Clamp K
if isempty(nModes) || nModes <= 0
    K = nModesAll;
else
    K = min(nModes, nModesAll);
end
Phi = Phi(:, 1:K);

% Compose: L̃ = L · Φ   [nCh × K]
L_tilde = Lc * Phi;

% Build composed head-model struct
CompHM = HeadModel;
CompHM.Gain        = L_tilde;
CompHM.GridLoc     = [];
CompHM.GridOrient  = [];
CompHM.GridAtlas   = [];
CompHM.isEigenmode = 1;
CompHM.nModes      = K;
CompHM.Eigenvalues = Values(1:K);
CompHM.HeadModelType = 'surface';
CompHM.Comment     = sprintf('Eigenmode leadfield (%d modes) | %s', K, ...
    getfield_default(HeadModel, 'Comment', ''));
end

function v = getfield_default(s, f, d)
if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_eigenmode_leadfield_pure.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/forward/bst_eigenmode_leadfield.m dev/tests/test_eigenmode_leadfield_pure.m
git commit -m "$(printf 'Eigenmode forward: bst_eigenmode_leadfield composes L*Phi\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

### Task 2: `process_eigenmode_leadfield` process node

**Files:**
- Create: `toolbox/process/functions/process_eigenmode_leadfield.m`
- Test: `dev/tests/test_eigenmode_leadfield_e2e.m`

- [ ] **Step 1: Write the failing e2e test**

```matlab
function test_eigenmode_leadfield_e2e
% Smoke: build a composed eigenmode head model from a real study and verify it
% is a valid headmodel node with Gain=[nCh x K] and carried metadata.
% Skips cleanly without a suitable protocol.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sProtocol = bst_get('ProtocolStudies');
if isempty(sProtocol) || ~isfield(sProtocol,'Study') || isempty(sProtocol.Study)
    disp('SKIP: no protocol loaded.'); return;
end
iStudyTarget = [];
for iS = 1:numel(sProtocol.Study)
    s = sProtocol.Study(iS);
    if isfield(s,'iHeadModel') && ~isempty(s.iHeadModel) && s.iHeadModel >= 1 ...
            && numel(s.HeadModel) >= s.iHeadModel
        try
            hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0, 'HeadModelType','SurfaceFile');
            if strcmpi(hm.HeadModelType,'surface')
                [~, isEig] = in_tess_eigenmodes(hm.SurfaceFile);
                if isEig; iStudyTarget = iS; break; end
            end
        catch; end
    end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with surface head model + eigenmodes.'); return;
end

OutFiles = process_eigenmode_leadfield('Run', ...
    struct('options', struct('nmodes', struct('Value',{{0,'',0}}))), ...
    struct('iStudy', iStudyTarget, 'FileName', '', 'Comment', 'test'));
assert(~isempty(OutFiles), 'Process must produce a head model file.');
HM = in_bst_headmodel(OutFiles{1}, 0);
assert(isfield(HM,'isEigenmode') && HM.isEigenmode==1, 'Output must be flagged isEigenmode.');
assert(size(HM.Gain,2) == HM.nModes, 'Gain columns must equal nModes.');
assert(~isempty(HM.Eigenvalues), 'Eigenvalues must be stored.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_eigenmode_leadfield_e2e.m`
Expected: FAIL — `Undefined function 'process_eigenmode_leadfield'` (or SKIP if no protocol; develop against the OMEGA protocol where it will run).

- [ ] **Step 3: Write the implementation**

```matlab
function varargout = process_eigenmode_leadfield( varargin )
% PROCESS_EIGENMODE_LEADFIELD: Compose a base head model into the eigenmode basis.
%
% Produces a headmodel_eigenmode_*.mat node whose Gain = L*Phi (the eigenmode
% leadfield). Requires a surface head model and precomputed surface eigenmodes.
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Compute eigenmode leadfield';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = all available): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Composes the standard leadfield with the ' ...
        'surface eigenmodes (L*Phi).<BR>Requires a surface head model and precomputed eigenmodes.</FONT>'];
    sProcess.options.label_info.Type = 'label';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    n = sProcess.options.nmodes.Value{1};
    if n > 0; Comment = sprintf('Eigenmode leadfield (%d modes)', n);
    else;     Comment = 'Eigenmode leadfield (all modes)'; end
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    nModes = sProcess.options.nmodes.Value{1};

    [sStudy, iStudy] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.'); return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;

    % Load base head model (unconstrained gain + orientations)
    HeadModel = in_bst_headmodel(HeadModelFile, 0);
    if ~strcmpi(HeadModel.HeadModelType, 'surface')
        bst_report('Error', sProcess, sInputs, 'Eigenmode leadfield requires a surface head model.'); return;
    end

    % Load eigenmodes from the head model's surface
    [Eigenmodes, isComputed] = in_tess_eigenmodes(HeadModel.SurfaceFile);
    if ~isComputed
        bst_report('Error', sProcess, sInputs, ...
            ['No eigenmodes on surface: ' HeadModel.SurfaceFile '. Run "Compute eigenmodes" first.']); return;
    end

    % Vertex-count consistency (head model is unconstrained: 3 cols/vertex)
    nVertHM = size(HeadModel.Gain, 2) / 3;
    if nVertHM ~= size(Eigenmodes.Vectors, 1)
        bst_report('Error', sProcess, sInputs, sprintf( ...
            ['Head model has %d vertices but eigenmodes have %d.\nRecompute the head model ' ...
             '(computing eigenmodes may have repaired the mesh).'], nVertHM, size(Eigenmodes.Vectors,1))); return;
    end

    % Compose
    CompHM = bst_eigenmode_leadfield(HeadModel, Eigenmodes, 'nModes', nModes);
    CompHM = bst_history('add', CompHM, 'eigenmode_leadfield', ...
        sprintf('Composed eigenmode leadfield: %d modes from %s', CompHM.nModes, HeadModelFile));

    % Save as a new head model node
    StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_eigenmode');
    bst_save(OutputFile, CompHM, 'v7');

    % Register in the DB. Build the study-level descriptor by cloning the shape of
    % the existing head-model descriptor (avoids depending on a db_template case
    % name), then overwrite the fields we control.
    sHeadModel = sStudy.HeadModel(sStudy.iHeadModel);
    sHeadModel.FileName      = file_short(OutputFile);
    sHeadModel.Comment       = CompHM.Comment;
    sHeadModel.HeadModelType = 'surface';
    iHM = length(sStudy.HeadModel) + 1;
    sStudy.HeadModel(iHM) = sHeadModel;
    sStudy.iHeadModel     = iHM;
    bst_set('Study', iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', iStudy);

    OutputFiles{end+1} = file_short(OutputFile);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_eigenmode_leadfield_e2e.m` against the OMEGA protocol.
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmode_leadfield.m dev/tests/test_eigenmode_leadfield_e2e.m
git commit -m "$(printf 'Eigenmode forward: process node writes headmodel_eigenmode_*.mat\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Phase 2 — Mode-space inverse

### Task 3: `bst_eigenmode_prior` (spectral prior R)

**Files:**
- Create: `toolbox/math/bst_eigenmode_prior.m`
- Test: `dev/tests/test_eigenmode_prior_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_eigenmode_prior_pure
% Verify the spectral prior R from eigenvalues:
%   - 'flat'  -> all ones
%   - 'power' -> normalized lambda^(-alpha), DC handled
%   - 'log'   -> positive, decreasing in lambda, ratio-preserving, max==1
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

lambdas = [0; 1; 4; 9; 16; 25];     % includes DC=0; ascending
K = 5;

% flat
Rf = bst_eigenmode_prior(lambdas, K, 'flat', 0);
assert(isequal(size(Rf), [K 1]), 'R must be [K x 1].');
assert(all(Rf == 1), 'flat prior must be all ones.');

% power (alpha=1): proportional to lambda^-1 after DC swap, normalized to max 1
Rp = bst_eigenmode_prior(lambdas, K, 'power', 1);
assert(all(Rp > 0), 'power prior must be positive.');
assert(abs(max(Rp) - 1) < 1e-12, 'prior must be normalized to max 1.');
assert(all(diff(Rp) <= 1e-12), 'power prior must be non-increasing in lambda.');

% log (2026): positive, decreasing, normalized; smoother modes favored
Rl = bst_eigenmode_prior(lambdas, K, 'log', 0);
assert(all(Rl > 0), 'log prior must be positive (lambda normalized into (0,1)).');
assert(abs(max(Rl) - 1) < 1e-12, 'log prior must be normalized to max 1.');
assert(all(diff(Rl) <= 1e-9), 'log prior must be non-increasing in lambda.');

% ratio-preservation invariant: scaling all eigenvalues by c only shifts log-space,
% so the resulting log prior is unchanged after max-normalization.
Rl_scaled = bst_eigenmode_prior(lambdas * 7.3, K, 'log', 0);
assert(max(abs(Rl - Rl_scaled)) < 1e-9, 'log prior must be invariant to uniform eigenvalue scaling.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_eigenmode_prior_pure.m`
Expected: FAIL — `Undefined function 'bst_eigenmode_prior'`.

- [ ] **Step 3: Write the implementation**

```matlab
function R = bst_eigenmode_prior(lambdas, K, priorType, alpha)
% BST_EIGENMODE_PRIOR: Diagonal source-covariance prior R from LBO eigenvalues.
%
% USAGE:  R = bst_eigenmode_prior(lambdas, K, priorType, alpha)
%
% DESCRIPTION:
%     Returns the diagonal of the source-covariance prior R [K x 1] in eigenmode
%     space. R plays the role of the source covariance in the standard inverse
%     J = R*L'*(L*R*L' + lambda*C)^-1*d, replacing depth weighting. Larger R_k =
%     more prior variance for mode k.
%
%     priorType:
%       'flat'  : R = ones(K,1)                         (no spectral prior)
%       'power' : R ∝ lambda_k^(-alpha)                 (legacy 1/f-like)
%       'log'   : R ∝ -log(lambda_k / lambda_ref)   (2026 GBF log prior; covariance)
%
%     For 'log', eigenvalues are normalized into (0,1) by lambda_ref (the first
%     discarded eigenvalue, i.e. lambda(K+1), else lambda(K)*(1+eps)). This is a
%     pure shift in log-space: it preserves eigenvalue ratios/ordering and only
%     guarantees -1/log(.) > 0. The DC mode (lambda~0) is swapped to lambda(2).
%     R is normalized so max(R) = 1; absolute scale is absorbed by the global
%     regularizer in the inverse.
%
% Authors: Diellor Basha, 2026

lambdas = double(lambdas(:));
nAvail  = numel(lambdas);
K = max(1, min(K, nAvail));

% DC handling: replace a (near-)zero leading eigenvalue with the next one
lam = lambdas;
if lam(1) <= max(lam) * 1e-12 && nAvail >= 2
    lam(1) = lam(2);
end

switch lower(priorType)
    case 'flat'
        R = ones(K, 1);
        return;

    case 'power'
        lamK = lam(1:K);
        lamK = max(lamK, max(lamK) * 1e-12);
        R = lamK .^ (-alpha);

    case 'log'
        % Reference scale = first discarded eigenvalue (else just above max used)
        if nAvail >= K+1
            lamRef = lambdas(K+1);
        else
            lamRef = lam(K) * (1 + 1e-6);
        end
        lamRef = max(lamRef, lam(K) * (1 + 1e-12));   % ensure strictly > lam(K)
        lamTilde = lam(1:K) / lamRef;                  % in (0,1)
        lamTilde = min(lamTilde, 1 - 1e-12);
        R = -log(lamTilde);                            % = log(lamRef/lam) > 0, decreasing in lambda
                                                       % (R is the covariance; GBF's -1/log lambda is the precision)

    otherwise
        error('bst_eigenmode_prior:UnknownPrior', 'Unknown priorType: %s', priorType);
end

% Normalize to max 1
R = R / max(R);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_eigenmode_prior_pure.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmode_prior.m dev/tests/test_eigenmode_prior_pure.m
git commit -m "$(printf 'Eigenmode inverse: spectral prior R from eigenvalues (log/flat/power)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

### Task 4: Rewrite `bst_inverse_eigenmodes` (mode-space solver)

**Files:**
- Rewrite: `toolbox/inverse/bst_inverse_eigenmodes.m`
- Test: `dev/tests/test_inverse_eigenmodes_pure.m`

The new engine takes the **composed** head model file (Gain already `= L·Φ`, with `.Eigenvalues`, `.nModes`), the noise covariance file, the channel file, and a good-channel mask. It cleans (projector + whitener), builds `R` via `bst_eigenmode_prior`, solves the regularized MAP estimate, and returns `M̃ [K × nGoodChannels]`.

New signature:
```matlab
[Results, errMsg] = bst_inverse_eigenmodes(CompHeadModelFile, NoiseCovFile, ChannelFile, GoodChannel, varargin)
%   OPTIONS: 'Method' {'mne','dspm','sloreta'}, 'Prior' {'log','flat','power'},
%            'Alpha' (power exponent), 'SNR', 'Unreg' (logical: pinv, ignore SNR),
%            'DataTypes' (cell, default {'MEG','MEG MAG','MEG GRAD'})
```

- [ ] **Step 1: Write the failing pure test** (math core: no DB, exercise the solver via an injected-operator path)

> The DB-coupled cleaning (channel file, noise cov file) is covered by the e2e test in Task 5. This pure test exercises the math kernel directly through a thin internal entry `bst_inverse_eigenmodes('SolvePure', ...)` so the solver is unit-tested without a protocol.

```matlab
function test_inverse_eigenmodes_pure
% Verify the mode-space solver math (no DB), via the 'SolvePure' entry:
%   Kernel = SolvePure(L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)
%   - shapes [K x nCh]; finite
%   - projector folding: Kernel annihilates the projected-out direction
%   - harmonic limit: flat prior + Unreg reproduces pinv(iW*L_tilde)*iW
%   - dSPM rows are unit noise-normalized in whitened space
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 8; K = 5;
L_tilde = zeros(nCh, K);
for i = 1:nCh
    for j = 1:K
        L_tilde(i,j) = sin(i * j * pi / 13);
    end
end
lambdas = (0:K-1)';                       % ascending, DC=0
iW   = diag(linspace(1, 2, nCh));         % whitener
Proj = eye(nCh);                          % no projector for base cases

% --- harmonic limit: flat + Unreg == pinv(iW*L)*iW ---
Kh = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'mne', 'flat', 0, 3, true);
Mref = pinv(iW * L_tilde) * iW;
assert(isequal(size(Kh), [K nCh]), 'Kernel must be [K x nCh].');
assert(max(abs(Kh(:) - Mref(:))) < 1e-8, 'flat+Unreg must equal pinv(iW*L)*iW.');

% --- regularized MNE with log prior: finite, correct shape ---
Km = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'mne', 'log', 0, 3, false);
assert(all(isfinite(Km(:))), 'Regularized MNE kernel must be finite.');

% --- projector folding: a rank-1 projector removes channel direction p ---
p = zeros(nCh,1); p(3) = 1;               % project out channel 3 subspace
Proj1 = eye(nCh) - (p*p')/(p'*p);
Kp = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj1, 'mne', 'flat', 0, 3, false);
assert(max(abs(Kp * p)) < 1e-9, 'Kernel must annihilate the projected-out direction.');

% --- dSPM: rows unit-normalized by noise std in whitened space ---
Kd = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, 'dspm', 'flat', 0, 3, false);
% noise std per row in whitened space = sqrt(sum((row*inv(iW)).^2))? dSPM normalizes
% by sqrt(sum(K_white.^2,2)); after that the whitened-noise row norm is ~1.
Kd_white = Kd / iW;                        % undo final de-whitening
rn = sqrt(sum(Kd_white.^2, 2));
assert(max(abs(rn - 1)) < 1e-6, 'dSPM rows must have unit noise norm in whitened space.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_inverse_eigenmodes_pure.m`
Expected: FAIL — old `bst_inverse_eigenmodes` has no `'SolvePure'` entry (errors on the string head-model path).

- [ ] **Step 3: Write the rewritten implementation**

```matlab
function [Results, errMsg] = bst_inverse_eigenmodes(varargin)
% BST_INVERSE_EIGENMODES: Mode-space MEG/EEG inverse on a composed eigenmode leadfield.
%
% USAGE:
%   [Results, errMsg] = bst_inverse_eigenmodes(CompHeadModelFile, NoiseCovFile, ChannelFile, GoodChannel, ...)
%   Kernel            = bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)
%
% DESCRIPTION:
%   Consumes a composed eigenmode head model (Gain = L*Phi, with .Eigenvalues and
%   .nModes; produced by bst_eigenmode_leadfield) and solves the regularized MAP
%   estimate in mode space:
%       M̃ = R L̃' (L̃ R L̃' + lambda C)^-1        [K x nGoodChannels]
%   where R = bst_eigenmode_prior(lambdas, K, Prior, Alpha) is the spectral source
%   covariance (replaces depth weighting) and C is the noise covariance (applied as
%   the whitener iW). SSP projectors and bad channels are folded into the clean
%   operator iW*Proj exactly as the standard inverse does.
%
%   The math core is exposed as ('SolvePure', ...) for unit testing.
%
% OPTIONS (name-value):
%   'Method'    : 'mne' (default) | 'dspm' | 'sloreta'
%   'Prior'     : 'log' (default) | 'flat' | 'power'
%   'Alpha'     : power-prior exponent (default 1)
%   'SNR'       : signal-to-noise ratio for regularization (default 3)
%   'Unreg'     : logical; if true, ignore SNR and use rank-safe pinv (default false)
%   'nModes'    : cap on modes used (default: all in the composed model)
%   'DataTypes' : channel types for the whitener (default {'MEG','MEG MAG','MEG GRAD'})
%
% Authors: Diellor Basha, 2026

% Dispatch the pure math entry
if ischar(varargin{1}) && strcmpi(varargin{1}, 'SolvePure')
    [L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg] = varargin{2:10};
    K = size(L_tilde, 2);
    R = bst_eigenmode_prior(lambdas, K, Prior, Alpha);
    Results = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg);
    errMsg = '';
    return;
end

Results = []; errMsg = '';
[CompHeadModelFile, NoiseCovFile, ChannelFile, GoodChannel] = varargin{1:4};

% Options
Method = 'mne'; Prior = 'log'; Alpha = 1; SNR = 3; Unreg = false; nModes = [];
DataTypes = {'MEG','MEG MAG','MEG GRAD'};
for i = 5:2:numel(varargin)
    switch lower(varargin{i})
        case 'method',    Method    = lower(varargin{i+1});
        case 'prior',     Prior     = lower(varargin{i+1});
        case 'alpha',     Alpha     = varargin{i+1};
        case 'snr',       SNR       = varargin{i+1};
        case 'unreg',     Unreg     = logical(varargin{i+1});
        case 'nmodes',    nModes    = varargin{i+1};
        case 'datatypes', DataTypes = varargin{i+1};
    end
end
if ~ismember(Method, {'mne','dspm','sloreta'})
    errMsg = ['Unknown method: ' Method '. Use mne, dspm, or sloreta.']; return;
end

% Load composed head model
HM = in_bst_headmodel(CompHeadModelFile, 0);
if ~isfield(HM, 'isEigenmode') || ~HM.isEigenmode
    errMsg = 'Head model is not an eigenmode leadfield. Run "Compute eigenmode leadfield" first.'; return;
end
L_all = double(HM.Gain);                          % [nAllCh x K]
lambdas = double(HM.Eigenvalues(:));
nAllCh = size(L_all, 1);

% Good channels
if isempty(GoodChannel) || numel(GoodChannel) ~= nAllCh
    iGood = (1:nAllCh)';
else
    iGood = find(GoodChannel(:));
end
L_tilde = L_all(iGood, :);
K = size(L_tilde, 2);
if ~isempty(nModes) && nModes > 0 && nModes < K
    K = nModes; L_tilde = L_tilde(:, 1:K); lambdas = lambdas(1:K);
end
nCh = numel(iGood);

% Whitener iW from noise covariance (Brainstorm convention)
if ~isempty(NoiseCovFile)
    NC = load(file_fullpath(NoiseCovFile));
    NoiseCov = NC.NoiseCov(iGood, iGood);
    iW = bst_whitener(NoiseCov, ChannelFile, DataTypes, []);   % [nCh x nCh]
    if isempty(iW) || ~isequal(size(iW), [nCh nCh])
        iW = eye(nCh);
    end
else
    iW = eye(nCh);
end

% Projector from the channel file (SSP), restricted to good channels
Proj = eigenmode_projector(ChannelFile, iGood, nCh);

% Build prior and solve
R = bst_eigenmode_prior(lambdas, K, Prior, Alpha);
Kernel = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg);

Results = struct();
Results.ImagingKernel = Kernel;        % [K x nGoodChannels]
Results.nModes        = K;
Results.Method        = Method;
Results.Prior         = Prior;
Results.SNR           = SNR;
Results.Unreg         = Unreg;
Results.GoodChannel   = iGood;
Results.Eigenvalues   = lambdas(1:K);
Results.SourcePrior   = R;
Results.SurfaceFile   = HM.SurfaceFile;
end

% ---- core solver (whitened, projected, mode-space MAP) ----
function Kernel = solve_modespace(L_tilde, R, iW, Proj, Method, SNR, Unreg)
% Clean operator: project then whiten. Folded into the final kernel so it maps RAW data.
Cop  = iW * Proj;                         % [nCh x nCh]
Lc   = Cop * L_tilde;                     % cleaned compressed leadfield [nCh x K]
nCh  = size(Lc, 1);
K    = size(Lc, 2);
sP   = sqrt(R(:));                        % column scaling by sqrt(prior)
Lws  = Lc .* sP.';                        % [nCh x K]

if Unreg
    % Rank-safe pseudoinverse (harmonic limit). Prior column scaling cancels.
    Kernel = pinv(Lc) * Cop;
    return;
end

[U, S, V] = svd(Lws, 'econ');
s = diag(S);
Lambda = sum(s.^2) / (nCh * SNR^2);
alpha  = s ./ (s.^2 + Lambda);            % MNE filter factors
Kmne_white = diag(sP) * V * diag(alpha) * U';   % [K x nCh], whitened-data side

switch Method
    case 'mne'
        Kernel = Kmne_white * Cop;
    case 'dspm'
        nn = sqrt(sum(Kmne_white.^2, 2));
        nn(nn < max(nn)*1e-10) = max(nn)*1e-10;
        Kernel = (Kmne_white ./ nn) * Cop;
    case 'sloreta'
        sa = s .* alpha;
        ResDiag = sum((V.^2) .* (sa'), 2);       % diag of resolution matrix
        ResDiag(ResDiag < max(ResDiag)*1e-10) = max(ResDiag)*1e-10;
        Kernel = (Kmne_white ./ sqrt(ResDiag)) * Cop;
end
end

% ---- SSP projector assembly from a channel file ----
function Proj = eigenmode_projector(ChannelFile, iGood, nCh)
Proj = eye(nCh);
if isempty(ChannelFile); return; end
try
    ChannelMat = in_bst_channel(ChannelFile);
catch
    return;
end
if ~isfield(ChannelMat, 'Projector') || isempty(ChannelMat.Projector); return; end
P = ChannelMat.Projector;
if isstruct(P)
    U = [];
    for k = 1:numel(P)
        if isfield(P(k),'Status') && P(k).Status == 1 && ~isempty(P(k).Components)
            comps = P(k).Components(iGood, :);
            if isfield(P(k),'CompMask') && ~isempty(P(k).CompMask)
                comps = comps(:, logical(P(k).CompMask));
            end
            U = [U, comps]; %#ok<AGROW>
        end
    end
    if ~isempty(U)
        [Uo, ~] = qr(U, 0);
        Proj = eye(nCh) - Uo*Uo';
    end
elseif isnumeric(P) && isequal(size(P), [numel(iGood) numel(iGood)]+[0 0])
    Proj = P;     % already a matrix over all channels
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_inverse_eigenmodes_pure.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_inverse_eigenmodes.m dev/tests/test_inverse_eigenmodes_pure.m
git commit -m "$(printf 'Eigenmode inverse: rewrite to consume composed leadfield + spectral-prior MAP\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

### Task 5: Rewrite `process_eigenmodes_inverse` (emit nodes)

**Files:**
- Rewrite: `toolbox/process/functions/process_eigenmodes_inverse.m`
- Test: `dev/tests/test_eigenmodes_inverse_e2e.m`

- [ ] **Step 1: Write the failing e2e test**

```matlab
function test_eigenmodes_inverse_e2e
% Smoke: run the eigenmode inverse on a study that has a composed eigenmode head
% model + noise cov + imported data; verify a coefficient matrix node and a
% kernel-only cortex results node are produced with correct shapes. Skips cleanly.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study with an eigenmode head model + noise cov + a data file
iStudyTarget = []; iDataFile = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || s.iHeadModel < 1 || numel(s.HeadModel) < s.iHeadModel; continue; end
    try
        hm = in_bst_headmodel(s.HeadModel(s.iHeadModel).FileName, 0);
    catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode==1 ...
            && isfield(s,'NoiseCov') && ~isempty(s.NoiseCov) && ~isempty(s.NoiseCov(1).FileName) ...
            && ~isempty(s.Data)
        iStudyTarget = iS; iDataFile = 1; break;
    end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with composed eigenmode head model + noise cov + data.'); return;
end

sInput = struct('iStudy', iStudyTarget, 'FileName', sStudies.Study(iStudyTarget).Data(iDataFile).FileName, ...
                'Comment', 'eigtest');
sProcess = struct('options', struct( ...
    'method',     struct('Value','dspm'), ...
    'prior',      struct('Value','log'), ...
    'snr',        struct('Value',{{3,'',1}}), ...
    'nmodes',     struct('Value',{{0,'',0}}), ...
    'outputtype', struct('Value','both')));
OutFiles = process_eigenmodes_inverse('Run', sProcess, sInput);
assert(~isempty(OutFiles), 'Inverse must produce output files.');

% At least one results (cortex) and one matrix (coefficients)
hasRes = any(~cellfun(@isempty, regexp(OutFiles, 'results_', 'once')));
hasMat = any(~cellfun(@isempty, regexp(OutFiles, 'matrix_',  'once')));
assert(hasRes, 'Expected a cortex results node.');
assert(hasMat, 'Expected a coefficient matrix node.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_eigenmodes_inverse_e2e.m`
Expected: FAIL (old process expects a raw head model and lacks the new options/output contract) or SKIP without a composed head model. Develop against OMEGA after running Task 2's process there.

- [ ] **Step 3: Write the rewritten process** (replace the whole file)

```matlab
function varargout = process_eigenmodes_inverse( varargin )
% PROCESS_EIGENMODES_INVERSE: Mode-space source mapping on a composed eigenmode leadfield.
%
% Consumes a composed eigenmode head model (headmodel_eigenmode_*.mat, from
% process_eigenmode_leadfield) and produces eigenmode coefficients (matrix) and/or
% reconstructed cortex sources (kernel-only results: ImagingKernel = Phi*M̃).
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Eigenmode source mapping';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.method.Comment = {'MNE (minimum norm)', 'dSPM (noise-normalized)', ...
        'sLORETA (standardized)', 'Inverse method:'; 'mne', 'dspm', 'sloreta', ''};
    sProcess.options.method.Type  = 'radio_linelabel';
    sProcess.options.method.Value = 'dspm';
    sProcess.options.prior.Comment = {'Log (2026)', 'Flat (none)', 'Power (1/f)', 'Spectral prior:'; ...
        'log', 'flat', 'power', ''};
    sProcess.options.prior.Type  = 'radio_linelabel';
    sProcess.options.prior.Value = 'log';
    sProcess.options.snr.Comment = 'Signal-to-noise ratio: ';
    sProcess.options.snr.Type    = 'value';
    sProcess.options.snr.Value   = {3, '', 1};
    sProcess.options.nmodes.Comment = 'Number of eigenmodes (0 = all in head model): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.outputtype.Comment = {'Coefficients (matrix)', 'Sources (results)', 'Both', 'Output:'; ...
        'coefficients', 'sources', 'both', ''};
    sProcess.options.outputtype.Type  = 'radio_linelabel';
    sProcess.options.outputtype.Value = 'both';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sprintf('Eigenmode %s (%s prior)', upper(sProcess.options.method.Value), ...
        sProcess.options.prior.Value);
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    Method     = lower(sProcess.options.method.Value);
    Prior      = lower(sProcess.options.prior.Value);
    SNR        = sProcess.options.snr.Value{1};
    nModes     = sProcess.options.nmodes.Value{1};
    OutputType = lower(sProcess.options.outputtype.Value);

    [sStudy, ~] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.iHeadModel) || sStudy.iHeadModel < 1
        bst_report('Error', sProcess, sInputs, 'No head model for this study.'); return;
    end
    HeadModelFile = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    HM = in_bst_headmodel(HeadModelFile, 0);
    if ~isfield(HM,'isEigenmode') || ~HM.isEigenmode
        bst_report('Error', sProcess, sInputs, ...
            'Active head model is not an eigenmode leadfield. Run "Compute eigenmode leadfield" first.'); return;
    end

    NoiseCovFile = '';
    if ~isempty(sStudy.NoiseCov) && ~isempty(sStudy.NoiseCov(1).FileName)
        NoiseCovFile = sStudy.NoiseCov(1).FileName;
    else
        bst_report('Warning', sProcess, sInputs, 'No noise covariance: using identity whitening.');
    end
    ChannelFile = bst_get('ChannelFileForStudy', sStudy.FileName);
    ChannelMat  = in_bst_channel(ChannelFile);

    % Good channels from the first input's flags (MEG, else EEG)
    DataFlag = in_bst_data(sInputs(1).FileName, 'ChannelFlag');
    ChannelFlag = ones(numel(ChannelMat.Channel), 1);
    if isfield(DataFlag,'ChannelFlag') && ~isempty(DataFlag.ChannelFlag)
        ChannelFlag = DataFlag.ChannelFlag;
    end
    iSel = good_channel(ChannelMat.Channel, ChannelFlag, 'MEG');
    if isempty(iSel); iSel = good_channel(ChannelMat.Channel, ChannelFlag, 'EEG'); end
    if isempty(iSel); bst_report('Error', sProcess, sInputs, 'No good MEG/EEG channels.'); return; end
    GoodChannel = false(numel(ChannelMat.Channel), 1); GoodChannel(iSel) = true;

    % Solve
    [Inv, errMsg] = bst_inverse_eigenmodes(HeadModelFile, NoiseCovFile, ChannelFile, GoodChannel, ...
        'Method', Method, 'Prior', Prior, 'SNR', SNR, 'nModes', nModes);
    if ~isempty(errMsg); bst_report('Error', sProcess, sInputs, errMsg); return; end
    K = Inv.nModes;

    % Eigenmodes (Phi) for reconstruction
    [Eig, ~] = in_tess_eigenmodes(HM.SurfaceFile);
    Phi = double(Eig.Vectors(:, 1:K));               % [nVert x K]
    lambdas = Inv.Eigenvalues;

    for iInput = 1:numel(sInputs)
        sInput = sInputs(iInput);
        DataMat = in_bst_data(sInput.FileName);
        isRaw = isstruct(DataMat.F);
        [sStudyOut, iStudyOut] = bst_get('Study', sInput.iStudy);
        StudyDir = bst_fileparts(file_fullpath(sStudyOut.FileName));

        % Cortex results node (kernel-only): ImagingKernel = Phi * M̃
        if ismember(OutputType, {'sources','both'})
            ResMat = db_template('resultsmat');
            ResMat.ImagingKernel = Phi * Inv.ImagingKernel;     % [nVert x nGoodCh]
            ResMat.ImageGridAmp  = [];
            ResMat.nComponents   = 1;
            ResMat.Comment       = sprintf('Eigenmode %s (%d modes, %s) | %s', ...
                upper(Method), K, Prior, sInput.Comment);
            ResMat.Function      = ['eigenmode_' Method];
            ResMat.Time          = DataMat.Time;
            ResMat.DataFile      = sInput.FileName;
            ResMat.HeadModelFile = HeadModelFile;
            ResMat.HeadModelType = 'surface';
            ResMat.SurfaceFile   = HM.SurfaceFile;
            ResMat.GoodChannel   = iSel;
            ResMat.ChannelFlag   = ChannelFlag;
            ResMat.nAvg          = DataMat.nAvg; ResMat.Leff = DataMat.Leff;
            ResMat = bst_history('add', ResMat, 'eigenmodes_inverse', ...
                sprintf('Eigenmode %s, %d modes, prior=%s, SNR=%.1f', Method, K, Prior, SNR));
            OutFile = bst_process('GetNewFilename', StudyDir, 'results_eigeninverse');
            bst_save(OutFile, ResMat, 'v6');
            db_add_data(iStudyOut, OutFile, ResMat);
            OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
        end

        % Coefficients matrix node: theta = M̃ * d (imported data only)
        if ismember(OutputType, {'coefficients','both'}) && ~isRaw
            theta = Inv.ImagingKernel * double(DataMat.F(iSel, :));   % [K x nTime]
            MatMat = db_template('matrixmat');
            MatMat.Value       = theta;
            MatMat.Time        = DataMat.Time;
            MatMat.nAvg        = DataMat.nAvg; MatMat.Leff = DataMat.Leff;
            MatMat.SurfaceFile = HM.SurfaceFile;
            MatMat.Comment     = sprintf('EigenCoeffs %s (%d modes, %s) | %s', ...
                upper(Method), K, Prior, sInput.Comment);
            RowNames = cell(K,1);
            for k = 1:K; RowNames{k} = sprintf('Mode %d (lam=%.3g)', k, lambdas(k)); end
            MatMat.Description = RowNames;
            MatMat = bst_history('add', MatMat, 'eigenmodes_inverse', ...
                sprintf('Eigenmode coefficients %s, %d modes, prior=%s', Method, K, Prior));
            OutFile = bst_process('GetNewFilename', StudyDir, 'matrix_eigencoeffs');
            bst_save(OutFile, MatMat, 'v6');
            db_add_data(iStudyOut, OutFile, MatMat);
            OutputFiles{end+1} = file_short(OutFile); %#ok<AGROW>
        elseif ismember(OutputType, {'coefficients','both'}) && isRaw
            bst_report('Warning', sProcess, sInput, 'Coefficients require imported data; skipped for raw.');
        end
    end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_eigenmodes_inverse_e2e.m` against OMEGA (after composing a head model via Task 2).
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_inverse.m dev/tests/test_eigenmodes_inverse_e2e.m
git commit -m "$(printf 'Eigenmode inverse: process emits coefficient + cortex nodes from composed HM\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Phase 3 — Retire the old paths

### Task 6: Delete retired source-mapping units and their tests

**Files:**
- Delete: `toolbox/math/bst_eigenmodes_harmonic.m`, `toolbox/math/bst_eigenmodes_transform.m`, `toolbox/gui/view_eigenmodes_timeseries.m`
- Delete tests: `dev/tests/test_eigenmodes_harmonic_pure.m`, `dev/tests/test_harmonic_inverse_e2e.m`, `dev/tests/test_eigenmodes_transform_pure.m`, `dev/tests/test_view_eigenmodes_timeseries_pure.m`, `dev/tests/test_eigenmode_timeseries_e2e.m`, `dev/tests/test_eigenmode_viewer_e2e.m`, `dev/tests/test_eigenmode_viewer_synth.m`, `dev/tests/test_view_eigenmodes_pure.m`
- Modify: any caller referencing the deleted symbols (tree callbacks, panels, `EigenKernel` node handling).

- [ ] **Step 1: Find all references to retired symbols**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
grep -rn -E 'bst_eigenmodes_harmonic|bst_eigenmodes_transform|view_eigenmodes_timeseries|EigenKernel|eigenmode_harmonic' toolbox dev | grep -v -E '2026-06-02-eigenmode'
```
Expected: a list of call sites (tree_callbacks, node menus, the retired tests). Record them.

- [ ] **Step 2: Remove the call sites**

For each non-test hit (e.g. a context-menu entry in `toolbox/tree/tree_callbacks.m` that opens `view_eigenmodes_timeseries`, or any `ResMat.EigenKernel`/`'eigenmode_harmonic'` branch), delete the menu/handler block. The coefficients are now a standard matrix node; no bespoke viewer entry is needed. Make each edit minimal and leave a comment only where a non-obvious branch was removed.

- [ ] **Step 3: Delete the retired files**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git rm toolbox/math/bst_eigenmodes_harmonic.m toolbox/math/bst_eigenmodes_transform.m toolbox/gui/view_eigenmodes_timeseries.m
git rm dev/tests/test_eigenmodes_harmonic_pure.m dev/tests/test_harmonic_inverse_e2e.m \
       dev/tests/test_eigenmodes_transform_pure.m dev/tests/test_view_eigenmodes_timeseries_pure.m \
       dev/tests/test_eigenmode_timeseries_e2e.m dev/tests/test_eigenmode_viewer_e2e.m \
       dev/tests/test_eigenmode_viewer_synth.m dev/tests/test_view_eigenmodes_pure.m
```

- [ ] **Step 4: Verify no dangling references remain**

Run:
```bash
grep -rn -E 'bst_eigenmodes_harmonic|bst_eigenmodes_transform|view_eigenmodes_timeseries|EigenKernel|eigenmode_harmonic' toolbox dev | grep -v -E '2026-06-02-eigenmode'
```
Expected: no output.

- [ ] **Step 5: Sanity-run the eigenspectrum tests (must be untouched)**

Run via MCP: `dev/tests/test_eigenmodes_project_pure.m`, `dev/tests/test_eigenmodes_filter_pure.m`, `dev/tests/test_bst_eigenmodes_modekernel_pure.m`.
Expected: each prints `ALL TESTS PASSED` (these analyze surface fields and must not be affected).

- [ ] **Step 6: Commit**

```bash
git commit -m "$(printf 'Eigenmode: retire harmonic/transform paths and bespoke viewer\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Phase 4 — Validation harness

### Task 7: `bst_resolution_metrics` (Level 1)

**Files:**
- Create: `toolbox/inverse/bst_resolution_metrics.m`
- Test: `dev/tests/test_resolution_metrics_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_resolution_metrics_pure
% Verify resolution-matrix point-spread metrics on a synthetic kernel/leadfield:
%   - identity resolution -> zero localization error, minimal dispersion
%   - shapes; finite; localization error in mm
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% 4 vertices on a line 10mm apart; perfect kernel (Res = I) -> LE = 0
GridLoc = [0 0 0; 10 0 0; 20 0 0; 30 0 0] / 1000;   % meters
L = magic(4); L = L + 4*eye(4);                      % invertible leadfield [nCh=4 x nV=4]
Kern = inv(L);                                        % perfect inverse -> Res = I

M = bst_resolution_metrics(Kern, L, GridLoc);
assert(isfield(M,'LocError') && isfield(M,'SpatialDispersion'), 'Must return LE and SD.');
assert(numel(M.LocError) == 4, 'LE must be per-vertex.');
assert(max(M.LocError) < 1e-9, 'Perfect inverse must give zero localization error.');
assert(all(isfinite(M.SpatialDispersion)), 'SD must be finite.');

% A blurred kernel (average of neighbors) increases dispersion and LE
Blur = [0.5 0.5 0 0; 0.3 0.4 0.3 0; 0 0.3 0.4 0.3; 0 0 0.5 0.5];
Mb = bst_resolution_metrics(Blur * inv(L) * L, L, GridLoc);  %#ok<MINV>
assert(mean(Mb.SpatialDispersion) >= mean(M.SpatialDispersion), 'Blur must not reduce dispersion.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run: `dev/tests/test_resolution_metrics_pure.m`
Expected: FAIL — `Undefined function 'bst_resolution_metrics'`.

- [ ] **Step 3: Write the implementation**

```matlab
function Metrics = bst_resolution_metrics(Kernel, Leadfield, GridLoc)
% BST_RESOLUTION_METRICS: Point-spread metrics from the resolution matrix.
%
% USAGE:  Metrics = bst_resolution_metrics(Kernel, Leadfield, GridLoc)
%
% INPUTS:
%   Kernel    : [nSrc x nCh]   imaging kernel (vertex space)
%   Leadfield : [nCh  x nSrc]  constrained leadfield (vertex space)
%   GridLoc   : [nSrc x 3]     source positions in meters
%
% OUTPUT (per-vertex, [nSrc x 1] unless noted):
%   .Res                resolution matrix [nSrc x nSrc] = Kernel*Leadfield
%   .LocError           localization error (mm): distance from true source to PSF peak
%   .SpatialDispersion  PSF spread (mm): sqrt(sum(d^2 * psf^2)/sum(psf^2))
%   .OverallAmplitude   peak PSF amplitude per source (depth-bias indicator)
%
% Each column j of Res is the point-spread function (PSF) for source j.
%
% Authors: Diellor Basha, 2026

Res = Kernel * Leadfield;                 % [nSrc x nSrc]
nSrc = size(Res, 2);
LocError = zeros(nSrc, 1);
SpatialDispersion = zeros(nSrc, 1);
OverallAmplitude  = zeros(nSrc, 1);
mm = 1000;                                % meters -> mm

for j = 1:nSrc
    psf = abs(Res(:, j));
    [pk, iPk] = max(psf);
    OverallAmplitude(j) = pk;
    % Localization error: true source location j vs PSF peak location
    LocError(j) = norm(GridLoc(iPk, :) - GridLoc(j, :)) * mm;
    % Spatial dispersion about the true source location
    d = sqrt(sum((GridLoc - GridLoc(j, :)).^2, 2)) * mm;   % [nSrc x 1]
    w = psf.^2;
    if sum(w) > 0
        SpatialDispersion(j) = sqrt(sum(d.^2 .* w) / sum(w));
    end
end

Metrics = struct('Res', Res, 'LocError', LocError, ...
    'SpatialDispersion', SpatialDispersion, 'OverallAmplitude', OverallAmplitude);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run: `dev/tests/test_resolution_metrics_pure.m`
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_resolution_metrics.m dev/tests/test_resolution_metrics_pure.m
git commit -m "$(printf 'Eigenmode validation: resolution-matrix point-spread metrics\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

### Task 8: Validation harness script (Levels 1–3)

**Files:**
- Create: `toolbox/script/tutorial_eigenmodes_validation.m`

This is an integration script (run via MCP / MATLAB, not unit-TDD). It assumes the OMEGA tutorial protocol is loaded with its two subjects, an Elekta phantom protocol available for Level 3 ground truth, and standard inverses already runnable. It writes a markdown report + figures next to the existing `dev/tests/omega-icosphere-sourcemap-results.md`.

- [ ] **Step 1: Write the harness script**

```matlab
function tutorial_eigenmodes_validation(ReportDir)
% TUTORIAL_EIGENMODES_VALIDATION: Benchmark eigenmode source mapping vs Brainstorm defaults.
%
% Levels:
%   1. Resolution-matrix metrics (analytic): eigenmode-MAP vs wMNE/dSPM/sLORETA.
%   2. Ground-truth simulation sweep (depth x SNR) on real OMEGA geometry.
%   3. Real data: Elekta phantom localization error (ground truth) + OMEGA GBF-vs-dSPM.
%
% Writes <ReportDir>/eigenmode-validation-results.md and PNG figures.
%
% USAGE:  tutorial_eigenmodes_validation                 % default report dir = dev/tests
%         tutorial_eigenmodes_validation('/path/out')
%
% Authors: Diellor Basha, 2026

if nargin < 1 || isempty(ReportDir)
    ReportDir = bst_fullfile(fileparts(fileparts(mfilename('fullpath'))), '..', 'dev', 'tests');
end
if ~brainstorm('status'); brainstorm nogui; end
lines = {sprintf('# Eigenmode source mapping — validation (%s)', datestr(now))}; %#ok<TNOW1,DATST>

% ---------- LEVEL 1: resolution metrics ----------
lines{end+1} = sprintf('\n## Level 1 — Resolution matrix');
L1 = level1_resolution();   % struct with per-method LocError/SD/Amplitude summaries
for m = 1:numel(L1)
    lines{end+1} = sprintf('- %s: LE median %.1f mm, SD median %.1f mm, depth-bias slope %.3g', ...
        L1(m).name, median(L1(m).LocError), median(L1(m).SpatialDispersion), L1(m).depthSlope);
end

% ---------- LEVEL 2: simulation sweep ----------
lines{end+1} = sprintf('\n## Level 2 — Ground-truth simulation (depth x SNR)');
L2 = level2_simulation();   % table: method x depth x SNR -> DLE, AUC, corr
lines = [lines, format_table(L2)];

% ---------- LEVEL 3: real data ----------
lines{end+1} = sprintf('\n## Level 3 — Real data');
L3p = level3_phantom();     % phantom localization error vs true dipole positions
if ~isempty(L3p); lines = [lines, L3p]; else; lines{end+1} = '- Phantom: SKIP (protocol not found).'; end
L3o = level3_omega();       % OMEGA GBF-vs-dSPM qualitative agreement (2 subjects)
if ~isempty(L3o); lines = [lines, L3o]; else; lines{end+1} = '- OMEGA: SKIP (protocol not found).'; end

% ---------- write report ----------
outFile = bst_fullfile(ReportDir, 'eigenmode-validation-results.md');
fid = fopen(outFile, 'w'); fprintf(fid, '%s\n', lines{:}); fclose(fid);
fprintf('Validation report written to %s\n', outFile);
end

% ===== Level 1: build all four kernels on the same OMEGA head model =====
function out = level1_resolution()
out = struct('name',{},'LocError',{},'SpatialDispersion',{},'depthSlope',{});
T = find_omega_study_with_eigenmode();      % helper below; returns files or []
if isempty(T); return; end
% Constrained vertex leadfield from the BASE (non-eigenmode) head model
baseHM = in_bst_headmodel(T.baseHmFile, 1);             % ApplyOrient=1 -> [nCh x nVert]
Lc = double(baseHM.Gain(T.iGood, :));
GridLoc = baseHM.GridLoc;
% Eigenmode-MAP kernel (vertex space = Phi * M̃)
[Inv, ~] = bst_inverse_eigenmodes(T.compHmFile, T.ncFile, T.chFile, T.goodMask, ...
    'Method','dspm','Prior','log','SNR',3);
[Eig, ~] = in_tess_eigenmodes(T.surf); Phi = double(Eig.Vectors(:,1:Inv.nModes));
Keig = Phi * Inv.ImagingKernel;
% Standard kernels via bst_inverse_linear_2018 (wMNE/dSPM/sLORETA) on the base HM
Kstd = standard_kernels(T);                 % struct .mne .dspm .sloreta [nVert x nCh]
methods = {'eigenmode-dSPM', Keig; 'wMNE', Kstd.mne; 'dSPM', Kstd.dspm; 'sLORETA', Kstd.sloreta};
depth = vertex_depth(GridLoc);              % mm from inner-skull/centroid proxy
for m = 1:size(methods,1)
    Mx = bst_resolution_metrics(methods{m,2}, Lc, GridLoc);
    p = polyfit(depth, Mx.OverallAmplitude, 1);
    out(end+1) = struct('name',methods{m,1}, 'LocError',Mx.LocError, ...
        'SpatialDispersion',Mx.SpatialDispersion, 'depthSlope',p(1)); %#ok<AGROW>
end
end

% ===== Level 2: simulate known sources on OMEGA geometry, reconstruct, score =====
function T2 = level2_simulation()
T2 = struct('rows', {{}});
S = find_omega_study_with_eigenmode(); if isempty(S); return; end
depths = {'superficial','deep'}; snrs = [1 3 6 10];
rows = {};
for d = 1:numel(depths)
    for q = 1:numel(snrs)
        sim = simulate_sources_on_study(S, depths{d}, snrs(q));   % uses process_simulate_*
        recon = reconstruct_all(S, sim);                          % eigenmode + std methods
        for m = 1:numel(recon)
            rows{end+1} = struct('method',recon(m).name, 'depth',depths{d}, 'snr',snrs(q), ...
                'DLE',recon(m).DLE, 'AUC',recon(m).AUC, 'corr',recon(m).corr); %#ok<AGROW>
        end
    end
end
T2.rows = rows;
end

% ===== Level 3a: Elekta phantom localization error vs known dipoles =====
function lines = level3_phantom()
lines = {};
P = find_phantom_study(); if isempty(P); return; end
truth = P.dipoleXYZ;                          % [nDip x 3] known positions (m)
recon = reconstruct_all_phantom(P);           % eigenmode + dSPM + LCMV
for m = 1:numel(recon)
    lines{end+1} = sprintf('- Phantom %s: median LE %.1f mm (n=%d dipoles)', ...
        recon(m).name, median(recon(m).LE_mm), size(truth,1)); %#ok<AGROW>
end
end

% ===== Level 3b: OMEGA GBF-vs-dSPM agreement across 2 subjects =====
function lines = level3_omega()
lines = {};
subs = find_omega_subjects();   % up to 2; [] if none
if isempty(subs); return; end
for s = 1:numel(subs)
    a = subject_map(subs(s), 'eigenmode-dspm');
    b = subject_map(subs(s), 'dspm');
    r = corr(a(:), b(:));
    lines{end+1} = sprintf('- OMEGA subject %s: spatial corr(eigenmode-dSPM, dSPM) = %.3f', ...
        subs(s).Name, r); %#ok<AGROW>
end
end
```

> **Implementation note for the engineer:** the private helpers referenced above
> (`find_omega_study_with_eigenmode`, `standard_kernels`, `vertex_depth`,
> `simulate_sources_on_study`, `reconstruct_all`, `find_phantom_study`,
> `reconstruct_all_phantom`, `find_omega_subjects`, `subject_map`, `format_table`)
> are thin wrappers over existing Brainstorm calls: `bst_inverse_linear_2018`
> (with `InverseMeasure` ∈ {amplitude,dspm2018,sloreta}) for the standard kernels;
> `process_simulate_sources`/`process_simulate_recordings` for Level 2 ground truth;
> `bst_get`/`in_bst_*` for file discovery; `process_dipole_scanning` for phantom
> peak extraction. Build each helper test-first as a small local function and run
> the script via the MATLAB MCP against the loaded OMEGA + phantom protocols.
> Each helper returns `[]`/empty so the script SKIPs cleanly when a protocol is
> absent. Define every helper before first use; do not leave any unimplemented.

- [ ] **Step 2: Run the harness against OMEGA (+ phantom if loaded)**

Run via MCP: `tutorial_eigenmodes_validation`
Expected: prints "Validation report written to …"; report contains Level 1 numbers and Level 2 table; Levels 3 either populated or cleanly SKIPped.

- [ ] **Step 3: Review the report against acceptance criteria**

Open `dev/tests/eigenmode-validation-results.md`. Confirm: Level 1 eigenmode-dSPM spatial dispersion ≤ wMNE and a smaller depth-bias slope; Level 2 eigenmode DLE/AUC competitive-or-better than dSPM/sLORETA across SNR with less depth bias than wMNE; Level 3 phantom LE not worse than dSPM. Record any deviations as follow-up issues.

- [ ] **Step 4: Commit**

```bash
git add toolbox/script/tutorial_eigenmodes_validation.m dev/tests/eigenmode-validation-results.md
git commit -m "$(printf 'Eigenmode validation: harness (resolution + simulation + phantom/OMEGA)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Final verification

- [ ] Run the full new pure-test suite via MCP and confirm each prints `ALL TESTS PASSED`:
  `test_eigenmode_leadfield_pure`, `test_eigenmode_prior_pure`, `test_inverse_eigenmodes_pure`, `test_resolution_metrics_pure`.
- [ ] Run the e2e tests against OMEGA: `test_eigenmode_leadfield_e2e`, `test_eigenmodes_inverse_e2e` → `ALL TESTS PASSED`.
- [ ] Run the retained eigenspectrum tests → still `ALL TESTS PASSED` (no regressions).
- [ ] Confirm `grep` for retired symbols returns nothing (Task 6 Step 4).
- [ ] Update `doc/updates.txt` with a dated line summarizing the eigenmode forward/inverse split, and commit.

---

## Notes & decisions carried from the spec

- **Default `K`:** the composer keeps all available modes; the inverse defaults to all modes in the composed model. (Spec §6 left open whether to clamp to `nChannels`. The solver is well-posed for `K > nChannels` because the spectral prior `R` regularizes — GBF's regime — so we do **not** clamp by default. The harmonic/`Unreg` path *does* require `K ≤ nChannels`; its test uses `K=5 < nCh=8`.)
- **`λ_ref` for the log prior:** first discarded eigenvalue `λ_{K+1}`, else `λ_K·(1+ε)`. Implemented in `bst_eigenmode_prior`.
- **Whitening/projector reuse:** `bst_whitener` for the whitener; SSP assembled from `ChannelMat.Projector` in `eigenmode_projector`. Order: project → whiten, folded into the kernel so it maps RAW data.
