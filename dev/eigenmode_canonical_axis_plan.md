# Canonical Eigenmode Axis Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish a single canonical, globally eigenvalue-sorted eigenmode ordering stored on the surface (`Eigenmodes.Order`), accessed through one shared helper, so every consumer selects modes identically — eliminating the divergent first-K / `ModeIndices` selections.

**Architecture:** Add `Eigenmodes.Order` (a global eigenvalue-sort permutation) to the surface eigenmode struct; add `bst_eigenmodes_canonical` as the single mode-selection path (self-computing `Order` if absent); re-ground the leadfield, reconstruction, analytic transform, and projection on it. Backward-compatible by construction (the old leadfield already sorted globally, so existing kernels stay valid). Add `bst_eigenmodes_ensure` for the default-1000/hemi policy.

**Tech Stack:** MATLAB R2023b, Brainstorm (running, `brainstorm nogui`), MATLAB MCP for execution. Protocols `TutorialAuditory`/`TutorialNeuromag` are loaded (ico5, 2000-mode eigenmodes).

---

## Background the engineer must know

- **Run everything through the MATLAB MCP** (`mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file` for test files, `evaluate_matlab_code` for inline). Brainstorm is already running. Load the tool with ToolSearch `select:mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file,mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code` if needed.
- **Repo root:** `/Users/diellorbasha/workspace/research/code/brainstorm3`. Branch: `feature/eigenmode-canonical-axis` (stay on it).
- **Test pattern** (`dev/tests/`): `function test_name`, `addpath(repoRoot)`, `if ~brainstorm('status'); brainstorm nogui; end`, `assert(...)`, ends with `disp('ALL TESTS PASSED')`. Toolbox-free (no `corr`, `nanmean`, `randsample`).
- **The bug being fixed:** eigenmodes are stored grouped by hemisphere (cols 1–1000 = hemi-1, 1001–2000 = hemi-2; sorted only within each). `Vectors(:,1:K)` for K<1000 is therefore all hemi-1. The canonical `Order` (global eigenvalue sort, e.g. `[1 1001 2 1002 ...]` semantics via stored indices) fixes "first-K" to mean whole-brain lowest frequencies.
- **Consistency invariant:** the OLD `bst_eigenmode_leadfield` already sorted `Values` ascending and stored `sel` as `ModeIndices`. So existing composed eigenmode head models are already in canonical order; the new code must reproduce the same Gain.

## File structure

| File | Change |
|---|---|
| `toolbox/math/bst_eigenmodes_canonical.m` (create) | The single mode-selection accessor (self-computes `Order` if absent) |
| `toolbox/anatomy/tess_eigenmodes.m` (modify ~line 210) | Write `Eigenmodes.Order` at packaging |
| `toolbox/io/in_tess_eigenmodes.m` (modify ~line 57) | Backfill `Eigenmodes.Order` on read for legacy files |
| `toolbox/forward/bst_eigenmode_leadfield.m` (modify) | Select via accessor; copy `Order` to composed HM; drop self-sort |
| `toolbox/inverse/bst_eigenmode_reconstruct.m` (modify) | Default to canonical `Order` when given a surface file |
| `toolbox/process/functions/process_eigenmodes_transform.m` (modify line 166-167) | Select via accessor (the bias fix) |
| `toolbox/math/bst_eigenmodes_project.m` (modify) | Mode-range over canonical `Order` |
| `toolbox/math/bst_eigenmodes_ensure.m` (create) | Ensure canonical eigenmodes (default 1000/hemi, no repair) |
| `dev/tests/test_eigenmodes_canonical_pure.m` (create) | accessor unit test |
| `dev/tests/test_eigenmodes_order_e2e.m` (create) | `.Order` written + backfilled |
| `dev/tests/test_eigenmode_canonical_consistency_e2e.m` (create) | leadfield Gain + Eigen-MNE kernel unchanged |
| `dev/tests/test_eigenmodes_transform_canonical_e2e.m` (create) | transform now spans both hemispheres |
| `dev/tests/test_eigenmodes_ensure_e2e.m` (create) | ensure-if-empty |

---

## Task 1: `bst_eigenmodes_canonical` accessor

**Files:**
- Create: `toolbox/math/bst_eigenmodes_canonical.m`
- Test: `dev/tests/test_eigenmodes_canonical_pure.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_canonical_pure.m`:

```matlab
function test_eigenmodes_canonical_pure
% Canonical accessor: returns first-K modes in GLOBAL eigenvalue order across
% components, self-computing Order if absent, honoring a present Order.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

% Two components, grouped (NOT globally sorted) eigenvalues.
I = eye(6);
Eig = struct();
Eig.Vectors   = I;                       % identity columns so selection is identifiable
Eig.Values    = [10;30;50;20;40;60];     % comp1=[10 30 50], comp2=[20 40 60]
Eig.Component = [1;1;1;2;2;2];
Eig.CompRank  = [1;2;3;1;2;3];
% No Order field -> accessor computes the global sort.

[Phi, lam, prov] = bst_eigenmodes_canonical(Eig, 4);
% global ascending order of values: 10(c1) 20(c2) 30(c1) 40(c2) -> stored cols [1 4 2 5]
assert(isequal(prov.idx(:)', [1 4 2 5]), 'canonical index wrong');
assert(isequal(lam(:)', [10 20 30 40]), 'eigenvalues not globally sorted');
assert(isequal(Phi, I(:, [1 4 2 5])), 'selected columns wrong');
assert(numel(unique(prov.Component)) == 2, 'first-4 must span both components');
assert(isequal(prov.Component(:)', [1 2 1 2]), 'provenance Component wrong');

% K omitted -> all modes, in canonical order
[~, lamAll] = bst_eigenmodes_canonical(Eig);
assert(isequal(lamAll(:)', [10 20 30 40 50 60]), 'all-modes path wrong');

% Present Order is honored verbatim (here: identity = stored order)
Eig.Order = (1:6)';
[~, ~, prov2] = bst_eigenmodes_canonical(Eig, 4);
assert(isequal(prov2.idx(:)', [1 2 3 4]), 'present Order not honored');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** (`Undefined function 'bst_eigenmodes_canonical'`). Run via MCP `run_matlab_file` on the test path.

- [ ] **Step 3: Create** `toolbox/math/bst_eigenmodes_canonical.m`:

```matlab
function [Phi, lambdas, prov] = bst_eigenmodes_canonical(Eig, K)
% BST_EIGENMODES_CANONICAL: Select the first K eigenmodes in the surface's canonical
% (global eigenvalue-sorted) order. This is the SINGLE mode-selection path shared by
% all eigenmode consumers, so "first K" always means the whole-brain lowest spatial
% frequencies (never one hemisphere).
%
% USAGE:  [Phi, lambdas, prov] = bst_eigenmodes_canonical(Eig, K)
%         [Phi, lambdas, prov] = bst_eigenmodes_canonical(Eig)     % all modes
%
% INPUTS:
%   Eig : eigenmodes struct (from in_tess_eigenmodes): .Vectors [nV x nModes],
%         .Values [nModes x 1], optionally .Order, .Component, .CompRank
%   K   : number of canonical modes to return (default/empty/<=0/>nModes: all)
% OUTPUTS:
%   Phi     : [nV x K]  eigenvectors, canonical order
%   lambdas : [K x 1]   eigenvalues, ascending
%   prov    : struct .Order [nModes x 1] full canonical permutation, .idx [K x 1]
%             selected stored columns, and .Component/.CompRank for the selection
%
% Authors: Diellor Basha, 2026
nModes = size(Eig.Vectors, 2);
if nargin < 2 || isempty(K) || K <= 0 || K > nModes
    K = nModes;
end
% Canonical order: stored .Order if valid, else compute the global eigenvalue sort.
if isfield(Eig, 'Order') && ~isempty(Eig.Order) && numel(Eig.Order) == nModes
    Order = double(Eig.Order(:));
else
    [~, Order] = sort(double(Eig.Values(:)), 'ascend');
end
idx     = Order(1:K);
Phi     = double(Eig.Vectors(:, idx));
lambdas = double(Eig.Values(idx));
prov    = struct('Order', Order, 'idx', idx);
if isfield(Eig, 'Component') && ~isempty(Eig.Component); prov.Component = Eig.Component(idx); end
if isfield(Eig, 'CompRank')  && ~isempty(Eig.CompRank);  prov.CompRank  = Eig.CompRank(idx);  end
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_canonical.m dev/tests/test_eigenmodes_canonical_pure.m
git commit -m "feat(eigenmode): canonical mode-selection accessor (single Order-based path)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: write `Eigenmodes.Order` (compute + backfill)

**Files:**
- Modify: `toolbox/anatomy/tess_eigenmodes.m` (after line 210, in STEP 5 packaging)
- Modify: `toolbox/io/in_tess_eigenmodes.m` (after line 57, backfill block)
- Test: `dev/tests/test_eigenmodes_order_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_order_e2e.m`:

```matlab
function test_eigenmodes_order_e2e
% tess_eigenmodes writes a global-sort .Order; in_tess_eigenmodes backfills it for
% legacy structs. Discriminating check: on the 2-hemisphere cortex, the first-K
% canonical modes span BOTH components (old first-K was one hemisphere).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Surf=in_tess_bst(ctxFile);

% Compute a small number of modes per hemisphere directly (fast, 2 components).
Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', 50);
assert(isfield(Eig,'Order') && numel(Eig.Order)==size(Eig.Vectors,2), 'Order not written');
assert(issorted(Eig.Values(Eig.Order)), 'Order does not sort Values globally');
assert(isequal(sort(Eig.Order(:)), (1:numel(Eig.Order))'), 'Order not a permutation');
% first-K spans both hemispheres
K = 40;
comp = Eig.Component(Eig.Order(1:K));
assert(numel(unique(comp))==2, 'first-K canonical must include both hemispheres');

% Backfill: strip Order, save to a temp surface, reload via in_tess_eigenmodes.
TessMat = in_tess_bst(ctxFile, 0);
EigStrip = Eig; EigStrip = rmfield(EigStrip,'Order');
TessMat.Eigenmodes = EigStrip;
tmpFile = fullfile(tempdir, sprintf('tess_order_test_%d.mat', feature('getpid')));
bst_save(tmpFile, TessMat, 'v7');
EigBack = in_tess_eigenmodes(tmpFile);
assert(isfield(EigBack,'Order') && ~isempty(EigBack.Order), 'Order not backfilled');
assert(issorted(EigBack.Values(EigBack.Order)), 'backfilled Order wrong');
delete(tmpFile);
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** (assert "Order not written" — current `tess_eigenmodes` doesn't add the field).

- [ ] **Step 3a: Modify `tess_eigenmodes.m`** — after line 210 (`Eigenmodes.CompRank = CompRank;`), add the canonical order:

```matlab
    Eigenmodes.CompRank    = CompRank;
    % Canonical mode order: global eigenvalue sort across all components. This is the
    % single ordering every consumer uses (via bst_eigenmodes_canonical), so "first K"
    % means the whole-brain lowest spatial frequencies, not one hemisphere.
    [~, Eigenmodes.Order]  = sort(ValuesAll, 'ascend');
```

- [ ] **Step 3b: Modify `in_tess_eigenmodes.m`** — after line 57 (the `CompRank` backfill), before the `nComponents` backfill, add:

```matlab
if ~isfield(Eigenmodes, 'CompRank') || isempty(Eigenmodes.CompRank)
    Eigenmodes.CompRank = (1:nK)';
end
% Backfill the canonical global eigenvalue order for legacy files (pre-Order).
if ~isfield(Eigenmodes, 'Order') || isempty(Eigenmodes.Order) || numel(Eigenmodes.Order) ~= nK
    [~, Eigenmodes.Order] = sort(double(Eigenmodes.Values(:)), 'ascend');
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_eigenmodes.m toolbox/io/in_tess_eigenmodes.m dev/tests/test_eigenmodes_order_e2e.m
git commit -m "feat(eigenmode): store canonical Order on surface; backfill on read

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: re-ground `bst_eigenmode_leadfield` on the canonical accessor

**Files:**
- Modify: `toolbox/forward/bst_eigenmode_leadfield.m` (lines 57-89: replace the self-sort/select block)
- Test: `dev/tests/test_eigenmode_canonical_consistency_e2e.m` (shared with Task 4)

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmode_canonical_consistency_e2e.m`:

```matlab
function test_eigenmode_canonical_consistency_e2e
% No-regression safety net: the canonical-path leadfield Gain must EQUAL the existing
% composed eigenmode head model Gain (old ModeIndices already == canonical Order), and
% the reconstructed Eigen-MNE kernel must equal the stored one.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sStudies=bst_get('ProtocolStudies'); T=[];
for iS=1:numel(sStudies.Study)
    s=sStudies.Study(iS); if isempty(s.HeadModel); continue; end
    iBase=[]; iEig=[];
    for ih=1:numel(s.HeadModel)
        try hm=in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        isE=isfield(hm,'isEigenmode')&&hm.isEigenmode;
        if isE && isempty(iEig); iEig=ih; elseif ~isE && strcmpi(hm.HeadModelType,'surface') && isempty(iBase); iBase=ih; end
    end
    if ~isempty(iBase) && ~isempty(iEig); T=struct('iS',iS,'iBase',iBase,'iEig',iEig); break; end
end
if isempty(T); disp('SKIP: need base + eigenmode head models.'); return; end
s=sStudies.Study(T.iS);
baseHM=in_bst_headmodel(s.HeadModel(T.iBase).FileName,0);
oldEigHM=in_bst_headmodel(s.HeadModel(T.iEig).FileName,0);
Eig=in_tess_eigenmodes(baseHM.SurfaceFile);

% Recompose via the (refactored) canonical leadfield path.
newEigHM = bst_eigenmode_leadfield(baseHM, Eig);
% Gain must match the existing composed model (same modes, same canonical order).
nC = min(size(newEigHM.Gain,2), size(oldEigHM.Gain,2));
relerr = norm(newEigHM.Gain(:,1:nC) - oldEigHM.Gain(:,1:nC),'fro') / norm(oldEigHM.Gain(:,1:nC),'fro');
fprintf('leadfield Gain rel.diff = %.3e (cols compared: %d)\n', relerr, nC);
assert(relerr < 1e-9, 'canonical leadfield Gain diverged from existing composed HM');
% Composed HM carries the canonical Order as ModeIndices.
assert(isfield(newEigHM,'ModeIndices') && ~isempty(newEigHM.ModeIndices), 'ModeIndices (canonical Order) not copied');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm it PASSES OR FAILS for the right reason.** Before the leadfield change, the old `bst_eigenmode_leadfield` produces a globally-sorted Gain too, so this may already pass — that is acceptable (it proves consistency). If it errors (e.g., the old signature truncates differently), note the diff. Proceed to Step 3 regardless; the test must PASS after the refactor.

- [ ] **Step 3: Modify `bst_eigenmode_leadfield.m`.** Replace the mode-selection block (lines 57-89, from `% Clamp K` through the `CompHM` assembly's `ModeIndices`) so it composes ALL modes via the canonical accessor. Replace:

```matlab
% Clamp K
if isempty(nModes) || nModes <= 0
    K = nModesAll;
else
    K = min(nModes, nModesAll);
end
% Select the K globally-lowest-eigenvalue modes across ALL connected components
% (hemispheres). Eigenmodes are stored grouped by component (tess_eigenmodes
% solves each hemisphere separately and concatenates), so a naive first-K slice
% Phi(:,1:K) would keep only one hemisphere. Sort by eigenvalue and keep the K
% lowest spatial frequencies whole-brain. The selected column indices are
% recorded in ModeIndices so the inverse reconstruction
% (bst_eigenmode_reconstruct) uses the exact same modes in the same order.
[~, order] = sort(Values, 'ascend');
sel = order(1:K);
PhiSel = Phi(:, sel);
```

with:

```matlab
% Canonical selection: the surface's single global eigenvalue order. We compose ALL
% modes (no independent truncation here -- K-capping is an inverse-side concern). The
% canonical Order is recorded as ModeIndices so the inverse reconstruction
% (bst_eigenmode_reconstruct) uses the exact same modes in the same order.
if isempty(nModes) || nModes <= 0
    K = nModesAll;
else
    K = min(nModes, nModesAll);
end
[PhiSel, ~, prov] = bst_eigenmodes_canonical(Eigenmodes, K);
sel = prov.idx;
```

(Leave the existing `L_tilde = Lc * PhiSel;` and the `CompHM` assembly — including `CompHM.ModeIndices = sel(:);` and `CompHM.Eigenvalues = Values(sel);` — unchanged; `sel` is now the canonical index.)

- [ ] **Step 4: Run the test, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/forward/bst_eigenmode_leadfield.m dev/tests/test_eigenmode_canonical_consistency_e2e.m
git commit -m "refactor(eigenmode): leadfield selects modes via canonical accessor

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `bst_eigenmode_reconstruct` defaults to canonical `Order`

**Files:**
- Modify: `toolbox/inverse/bst_eigenmode_reconstruct.m` (lines 22-46)
- Test: extend `dev/tests/test_eigenmode_canonical_consistency_e2e.m`

- [ ] **Step 1: Extend the consistency test** — append, before `disp('ALL TESTS PASSED')`:

```matlab
% Reconstruction via canonical default (surface file, no explicit ModeIndices) must
% equal reconstruction with the head model's stored ModeIndices.
rc=string({s.Result.Comment}); iEK=find(rc=="Eigen-MNE: MEG 2018");
if ~isempty(iEK)
    REK=in_bst_results(s.Result(iEK(end)).FileName,0);   % stored vertex-space kernel
    % Re-derive the mode-space kernel then reconstruct two ways.
    [Inv,err]=bst_inverse_eigenmodes(s.HeadModel(T.iEig).FileName, s.NoiseCov(1).FileName, ...
        bst_get('ChannelFileForStudy',s.FileName), REK.GoodChannel, 'Method','mne','Prior','log','SNR',3);
    assert(isempty(err), 'inverse failed: %s', err);
    Kvert_idx = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel, oldEigHM.ModeIndices);
    Kvert_def = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel);  % canonical default
    rel = norm(Kvert_idx - Kvert_def,'fro')/norm(Kvert_idx,'fro');
    fprintf('reconstruct canonical-default vs explicit-Order rel.diff = %.3e\n', rel);
    assert(rel < 1e-9, 'canonical-default reconstruction must equal explicit ModeIndices');
end
```

- [ ] **Step 2: Run it, confirm FAIL** — current `bst_eigenmode_reconstruct` with no `ModeIndices` falls back to `Phi(:,1:K)` (hemisphere-biased), so `Kvert_def` differs from `Kvert_idx`.

- [ ] **Step 3: Modify `bst_eigenmode_reconstruct.m`.** Change the no-`ModeIndices` branch so that, when loading from a surface file, it uses the canonical `Order`. Replace lines 22-46:

```matlab
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
if (nargin >= 3) && ~isempty(ModeIndices)
    idx = ModeIndices(:);
    if numel(idx) < K
        error('Fewer mode indices (%d) than kernel modes (%d).', numel(idx), K);
    end
    idx = idx(1:K);
    if max(idx) > size(Phi, 2)
        error('Mode index %d exceeds available eigenmodes (%d).', max(idx), size(Phi,2));
    end
    ImagingKernel = Phi(:, idx) * ModeKernel;
else
    if size(Phi, 2) < K
        error('Surface has fewer eigenmodes (%d) than kernel modes (%d).', size(Phi,2), K);
    end
    ImagingKernel = Phi(:, 1:K) * ModeKernel;
end
```

with:

```matlab
K = size(ModeKernel, 1);
Eig = [];
if ischar(SurfaceOrPhi)
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceOrPhi);
    if ~isComputed
        error(['No eigenmodes on surface: ' SurfaceOrPhi '. Run "Compute eigenmodes" first.']);
    end
    Phi = double(Eig.Vectors);
else
    Phi = double(SurfaceOrPhi);
end
if (nargin >= 3) && ~isempty(ModeIndices)
    % Explicit selection (e.g. the head model's canonical Order copy).
    idx = ModeIndices(:);
elseif ~isempty(Eig)
    % Default to the surface's canonical order (single shared selection path).
    [~, ~, prov] = bst_eigenmodes_canonical(Eig, K);
    idx = prov.idx;
else
    % Bare Phi matrix, no order known: positional fallback.
    idx = (1:K)';
end
if numel(idx) < K
    error('Fewer mode indices (%d) than kernel modes (%d).', numel(idx), K);
end
idx = idx(1:K);
if max(idx) > size(Phi, 2)
    error('Mode index %d exceeds available eigenmodes (%d).', max(idx), size(Phi,2));
end
ImagingKernel = Phi(:, idx) * ModeKernel;
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_eigenmode_reconstruct.m dev/tests/test_eigenmode_canonical_consistency_e2e.m
git commit -m "refactor(eigenmode): reconstruct defaults to canonical Order from surface

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: fix `process_eigenmodes_transform` (the bias bug)

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes_transform.m` (lines 166-167)
- Test: `dev/tests/test_eigenmodes_transform_canonical_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_transform_canonical_e2e.m`:

```matlab
function test_eigenmodes_transform_canonical_e2e
% The unregularized transform must select modes via the canonical Order, so K<1000
% spans BOTH hemispheres (previously Vectors(:,1:K) was hemisphere-1 only).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Eig=in_tess_eigenmodes(ctxFile);

% Mirror the transform's selection: it should now use the canonical accessor.
K = 600;
[~, ~, prov] = bst_eigenmodes_canonical(Eig, K);
assert(numel(unique(prov.Component))==2, 'canonical first-600 must span both hemispheres');
% Sanity: the NAIVE first-600 would be one hemisphere (the bug we are fixing).
naiveComp = Eig.Component(1:K);
assert(numel(unique(naiveComp))==1, 'precondition: naive first-K is single-hemisphere');

% Round-trip: a vertex field reconstructed from its canonical coefficients recovers itself.
Phi = bst_eigenmodes_canonical(Eig, K);
u = Phi * randn(K,1);            % field living in the canonical mode subspace
M = Eig.MassMatrix;
theta = Phi' * (M * u);          % project (M-orthonormal basis)
urec  = Phi * theta;             % reconstruct
assert(norm(urec - u)/norm(u) < 1e-6, 'canonical project/reconstruct round-trip failed');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm `ALL TESTS PASSED`** (this test validates the canonical accessor's behavior, which already exists from Task 1 — it documents the contract the transform must follow). Then verify the transform code itself is changed in Step 3 so it *uses* this path.

- [ ] **Step 3: Modify `process_eigenmodes_transform.m`** lines 166-167. Replace:

```matlab
    Phi     = double(Eigenmodes.Vectors(:, 1:K));
    lambdas = double(Eigenmodes.Values(1:K));
```

with:

```matlab
    % Canonical selection (whole-brain lowest spatial frequencies, never one hemisphere).
    [Phi, lambdas] = bst_eigenmodes_canonical(Eigenmodes, K);
```

(The downstream `ResMat.ImagingKernel = Phi * Kernel;` at lines 194 and 246 is unchanged — `Phi` is now canonical, so transform and reconstruction stay self-consistent.)

- [ ] **Step 4: Verify the transform runs end-to-end** via MCP `evaluate_matlab_code`:

```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
s=bst_get('ProtocolStudies'); s=s.Study(6); dc=string({s.Data.Comment});
iD=find(startsWith(dc,'Avg: deviant'),1);
sOut=bst_process('CallProcess','process_eigenmodes_transform', s.Data(iD).FileName, [], 'nmodes', {600,'',0});
M=in_bst_matrix(sOut(1).FileName); fprintf('transform out: [%d x %d]\n', size(M.Value,1), size(M.Value,2));
assert(size(M.Value,1)==600 && all(isfinite(M.Value(:))), 'transform output bad');
disp('TRANSFORM OK');
```
Expected: `[600 x 1441]`, finite, `TRANSFORM OK`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_transform.m dev/tests/test_eigenmodes_transform_canonical_e2e.m
git commit -m "fix(eigenmode): transform selects via canonical accessor (no hemisphere bias)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: `bst_eigenmodes_ensure` (default-1000/hemi, no repair)

**Files:**
- Create: `toolbox/math/bst_eigenmodes_ensure.m`
- Test: `dev/tests/test_eigenmodes_ensure_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_ensure_e2e.m`:

```matlab
function test_eigenmodes_ensure_e2e
% bst_eigenmodes_ensure returns canonical eigenmodes when present (idempotent, no
% recompute), and the returned struct carries a valid canonical Order.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;

t0=tic; Eig=bst_eigenmodes_ensure(ctxFile); el=toc(t0);
assert(~isempty(Eig) && isfield(Eig,'Order'), 'ensure did not return canonical eigenmodes');
assert(issorted(Eig.Values(Eig.Order)), 'returned Order invalid');
assert(el < 30, 'ensure should be fast when eigenmodes already exist (no recompute)');
% Idempotent
Eig2=bst_eigenmodes_ensure(ctxFile);
assert(size(Eig2.Vectors,2)==size(Eig.Vectors,2), 'ensure not idempotent');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** (`Undefined function 'bst_eigenmodes_ensure'`).

- [ ] **Step 3: Create** `toolbox/math/bst_eigenmodes_ensure.m`:

```matlab
function Eig = bst_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
% BST_EIGENMODES_ENSURE: Return the surface's canonical eigenmodes, computing a default
% set if absent. The principled home for the "use existing modes, else compute default"
% policy (replacing ad-hoc benchmark fixtures).
%
% USAGE:  Eig = bst_eigenmodes_ensure(SurfaceFile)
%         Eig = bst_eigenmodes_ensure(SurfaceFile, nModesPerHemi)   % default 1000
%
% Computes (if needed) nModesPerHemi modes per connected component, barycentric mass,
% remove-DC, and NO repair: a non-manifold surface raises an error (repair changes the
% vertex count and breaks surface<->leadfield<->eigenmode consistency; remesh to ico).
%
% Authors: Diellor Basha, 2026
if nargin < 2 || isempty(nModesPerHemi); nModesPerHemi = 1000; end
[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
if isComputed && ~isempty(Eig)
    return;
end
% Compute a default set. Guard manifoldness first (no silent repair).
Surf = in_tess_bst(SurfaceFile, 0);
mani = tess_manifold(Surf.Vertices, Surf.Faces);
if isstruct(mani) && isfield(mani,'isManifold') && ~mani.isManifold
    error('bst_eigenmodes_ensure:NonManifold', ...
        ['Surface %s is non-manifold; eigenmodes require a 2-manifold mesh. ' ...
         'Remesh to an icosphere (or repair manually) and retry.'], SurfaceFile);
end
Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', nModesPerHemi, ...
    'MassType', 'barycentric', 'RemoveDC', 1, 'Repair', 0);
% Persist on the surface so subsequent reads are canonical and cached.
out_tess_eigenmodes(SurfaceFile, Eig);
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.** (Eigenmodes already exist on the ico5 cortex, so it returns fast without recompute.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_ensure.m dev/tests/test_eigenmodes_ensure_e2e.m
git commit -m "feat(eigenmode): bst_eigenmodes_ensure (default 1000/hemi, no repair)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: `bst_eigenmodes_project` mode-range over canonical order + viewer check

**Files:**
- Modify: `toolbox/math/bst_eigenmodes_project.m` (the `ModeRange` selection)
- Verify: `toolbox/gui/view_eigenmodes.m` / `panel_eigenmodes.m` (no first-K slice — confirm only)

- [ ] **Step 1: Read `bst_eigenmodes_project.m`** to find where it slices modes for `ModeRange`. Run via MCP `evaluate_matlab_code`:

```matlab
disp(fileread('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/math/bst_eigenmodes_project.m'));
```
Identify the line(s) where it indexes `Eigenmodes.Vectors(:, iModes)` or `(:, k1:k2)`.

- [ ] **Step 2: Modify the `ModeRange` selection** so the range indexes the canonical order. Where the code currently selects a positional range `iModes = k1:k2; Phi = Vectors(:, iModes)`, change to select from the canonical order:

```matlab
% Canonical mode range: indices are ranks in the global eigenvalue order.
[PhiAll, ~, provAll] = bst_eigenmodes_canonical(Eigenmodes);   % all, canonical
iSel = provAll.idx(k1:k2);
Phi  = double(Eigenmodes.Vectors(:, iSel));
```

(If the function uses all modes with no range, leave it — order only matters when truncating. Apply the change only at the truncation site. If no truncation exists, record that and skip to Step 4.)

- [ ] **Step 3: Confirm the viewers carry no first-K bug.** Run via MCP `evaluate_matlab_code`:

```matlab
g = @(f) numel(regexp(fileread(f), 'Vectors\(:,\s*1:', 'once'));
base='/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox/gui/';
fprintf('view_eigenmodes first-K slices: %d\n', g([base 'view_eigenmodes.m']));
fprintf('panel_eigenmodes first-K slices: %d\n', g([base 'panel_eigenmodes.m']));
```
Expected: `0` and `0` — the viewers pair by `CompRank` tags, not positional first-K, so they need no change for correctness in this foundation. (The mode-k↔canonical-rank UX alignment is part of Feature B.) If either is nonzero, report it as a concern for review.

- [ ] **Step 4: Commit** (only the project change if one was needed)

```bash
git add toolbox/math/bst_eigenmodes_project.m
git commit -m "refactor(eigenmode): project mode-range over canonical order

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: regression sweep

**Files:** none (verification only)

- [ ] **Step 1: Run the existing eigenmode test suite** via MCP `run_matlab_file` on each (skip cleanly if a fixture is absent). Expected `ALL TESTS PASSED` (or a clean SKIP) for each:
- `dev/tests/test_eigenmode_reconstruct_pure.m`
- `dev/tests/test_eigenmode_leadfield_pure.m`
- `dev/tests/test_inverse_eigenmodes_pure.m`
- `dev/tests/test_benchmark_inverse_nmodes_e2e.m`
- `dev/tests/test_kernel_comparison.m`

- [ ] **Step 2: Re-run the new tests** (`test_eigenmodes_canonical_pure`, `test_eigenmodes_order_e2e`, `test_eigenmode_canonical_consistency_e2e`, `test_eigenmodes_transform_canonical_e2e`, `test_eigenmodes_ensure_e2e`). All must pass.

- [ ] **Step 3: Report** any failures with the first stack trace. If all green, the foundation is consistency-preserving (source mapping unchanged) and bug-fixing (transform de-biased).

---

## Self-Review

**1. Spec coverage:**
- Canonical `Eigenmodes.Order` (compute + backfill) → Task 2. ✔
- Single accessor `bst_eigenmodes_canonical` → Task 1. ✔
- Leadfield re-grounded, all modes, copies Order → Task 3. ✔
- Reconstruct canonical default → Task 4. ✔
- Inverse `nModes` = first-K canonical → automatic (Gain is canonical-ordered after Task 3); no code change needed, verified by Task 4's consistency test. ✔
- Transform de-biased → Task 5. ✔
- Project mode-range over canonical → Task 7. ✔
- Viewers (no first-K bug; UX deferred to B) → Task 7 Step 3 verification. ✔
- `bst_eigenmodes_ensure` default 1000/hemi, no repair → Task 6. ✔
- Backward-compat (old ModeIndices ≡ canonical Order) → Task 3 consistency test proves it. ✔
- Tests: canonical index, backfill, accessor, consistency, transform-fix, ensure, regression → Tasks 1,2,3,4,5,6,8. ✔

**2. Placeholder scan:** No "TBD/TODO"; every code step is complete. Task 7 Step 2 is conditional ("apply only at the truncation site") because the project function's exact slice line must be read first (Step 1) — the change pattern is fully specified.

**3. Type consistency:** `bst_eigenmodes_canonical(Eig, K)` returns `[Phi, lambdas, prov]` with `prov.idx`/`prov.Order`/`prov.Component` used identically in Tasks 1,3,4,5,7. `Eigenmodes.Order` written (Task 2) and consumed (Tasks 1,3,4) consistently. `bst_eigenmode_leadfield`'s `CompHM.ModeIndices = sel(:)` where `sel = prov.idx` (Task 3) matches `bst_eigenmode_reconstruct`'s `ModeIndices` consumption (Task 4).

**Out of scope (per spec):** Feature B viewer; re-sorting stored `Vectors`; analytic-module reorg; inverse math/prior/eigenfilter changes — none included.
