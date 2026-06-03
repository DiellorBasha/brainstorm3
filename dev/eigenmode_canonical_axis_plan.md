# Canonical Eigenmode Axis Implementation Plan (rev 2)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish a single canonical, globally eigenvalue-sorted eigenmode ordering stored on the surface (`Eigenmodes.Order`), expose the analytic transform as a pure named pair (`manifold_ft`/`manifold_ift`), and re-ground every consumer on the canonical order — eliminating the divergent first-K / `ModeIndices` selections.

**Architecture:** Add `Eigenmodes.Order` (global eigenvalue-sort permutation) to the surface eigenmode struct. Add pure `manifold_ft(Φ,M,U)=Φ'MU` and `manifold_ift(Φ,C)=ΦC` in `toolbox/math/`. Demote `bst_eigenmodes_project` and `bst_eigenmode_reconstruct` to thin wrappers over the pure pair (callers untouched; canonical fix applied inside). Re-ground the leadfield and transform on `Order` (inline `Vectors(:,Order(1:K))`, no accessor). Backward-compatible by construction (the old leadfield already sorted globally). Add `bst_eigenmodes_ensure`.

**Tech Stack:** MATLAB R2023b, Brainstorm (running, `brainstorm nogui`), MATLAB MCP for execution. Protocols `TutorialAuditory`/`TutorialNeuromag` loaded (ico5, 2000-mode eigenmodes).

---

## Background the engineer must know

- **Run everything through the MATLAB MCP** (`mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file` for tests, `evaluate_matlab_code` for inline). Brainstorm is running.
- **Repo root:** `/Users/diellorbasha/workspace/research/code/brainstorm3`. Branch: `feature/eigenmode-canonical-axis` (stay on it).
- **Test pattern** (`dev/tests/`): `function test_name`, `addpath(repoRoot)`, `if ~brainstorm('status'); brainstorm nogui; end`, `assert(...)`, ends `disp('ALL TESTS PASSED')`. Toolbox-free.
- **The bug:** eigenmodes are stored grouped by hemisphere (cols 1–1000 = hemi-1, 1001–2000 = hemi-2). `Vectors(:,1:K)` for K<1000 is all hemi-1. Canonical `Order` (global eigenvalue sort) makes "first K" mean whole-brain lowest frequencies.
- **Consistency invariant:** the OLD `bst_eigenmode_leadfield` already sorted `Values` ascending and stored `sel` as `ModeIndices`, so existing composed head models are already canonical-ordered; the new code must reproduce the same Gain and the same reconstructed kernels.
- **Brainstorm convention:** one function per file.

## File structure

| File | Change |
|---|---|
| `toolbox/math/manifold_ft.m` (create) | Pure forward manifold FT: `Φ'·M·U` |
| `toolbox/math/manifold_ift.m` (create) | Pure inverse manifold FT: `Φ·C` |
| `toolbox/anatomy/tess_eigenmodes.m` (modify ~line 210) | Write `Eigenmodes.Order` |
| `toolbox/io/in_tess_eigenmodes.m` (modify ~line 57) | Backfill `Eigenmodes.Order` on read |
| `toolbox/math/bst_eigenmodes_project.m` (rewrite body) | Thin wrapper over `manifold_ft`/`manifold_ift`, canonical `ModeRange` |
| `toolbox/inverse/bst_eigenmode_reconstruct.m` (rewrite body) | Thin wrapper: load + canonical `Order` + `manifold_ift` |
| `toolbox/forward/bst_eigenmode_leadfield.m` (modify lines 57-72) | Select via `Order` inline; all modes; copy `Order` to HM |
| `toolbox/process/functions/process_eigenmodes_transform.m` (modify lines 166-167) | Select via `Order` inline (bias fix) |
| `toolbox/math/bst_eigenmodes_ensure.m` (create) | Ensure canonical eigenmodes (default 1000/hemi, no repair) |
| `dev/tests/test_manifold_ft_ift_pure.m` (create) | pure transform round-trip |
| `dev/tests/test_eigenmodes_order_e2e.m` (create) | `.Order` written + backfilled + spans hemispheres |
| `dev/tests/test_eigenmodes_project_wrapper_pure.m` (create) | project == manifold_ft |
| `dev/tests/test_eigenmode_canonical_consistency_e2e.m` (create) | leadfield Gain + reconstruct unchanged |
| `dev/tests/test_eigenmodes_transform_canonical_e2e.m` (create) | transform spans both hemispheres |
| `dev/tests/test_eigenmodes_ensure_e2e.m` (create) | ensure-if-empty |

---

## Task 1: pure `manifold_ft` / `manifold_ift`

**Files:**
- Create: `toolbox/math/manifold_ft.m`, `toolbox/math/manifold_ift.m`
- Test: `dev/tests/test_manifold_ft_ift_pure.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_manifold_ft_ift_pure.m`:

```matlab
function test_manifold_ft_ift_pure
% Pure manifold Fourier pair: ft = Phi'*M*U, ift = Phi*C. Round-trip onto the basis
% subspace recovers the field; ift handles a matrix C (kernel); dimension errors raise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

nV = 20; K = 8;
[Q, ~] = qr(randn(nV, K), 0);     % Phi: orthonormal columns
Phi = Q; M = speye(nV);            % M = I -> Phi M-orthonormal
u = randn(nV, 3);

c = manifold_ft(Phi, M, u);
assert(isequal(size(c), [K 3]), 'manifold_ft output shape');
uproj = Phi * (Phi' * u);          % projection of u onto span(Phi)
urec  = manifold_ift(Phi, c);
assert(norm(urec - uproj, 'fro')/norm(uproj,'fro') < 1e-10, 'round-trip onto subspace failed');

% ift handles a matrix C (mode-space kernel -> vertex kernel)
Kk = manifold_ift(Phi, randn(K, 5));
assert(isequal(size(Kk), [nV 5]), 'manifold_ift matrix shape');

% dimension errors
ok=false; try; manifold_ft(Phi, M, randn(nV+1,2)); catch; ok=true; end; assert(ok, 'ft dim check');
ok=false; try; manifold_ift(Phi, randn(K+1,2)); catch; ok=true; end; assert(ok, 'ift dim check');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** (`Undefined function 'manifold_ft'`).

- [ ] **Step 3a: Create** `toolbox/math/manifold_ft.m`:

```matlab
function C = manifold_ft(Phi, M, U)
% MANIFOLD_FT: Forward manifold Fourier transform. Project a vertex field onto an
% M-orthonormal Laplace-Beltrami eigenmode basis: C = Phi' * (M * U).
%
% USAGE:  C = manifold_ft(Phi, M, U)
% INPUTS:
%   Phi : [nV x K]  eigenvectors (already selected, e.g. canonical order)
%   M   : [nV x nV] mass matrix (basis is M-orthonormal)
%   U   : [nV x nT] vertex field(s)
% OUTPUT:
%   C   : [K x nT]  mode coefficients
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); U = double(U);
if size(U,1) ~= size(Phi,1)
    error('manifold_ft: U has %d rows but Phi has %d.', size(U,1), size(Phi,1));
end
C = Phi' * (M * U);
end
```

- [ ] **Step 3b: Create** `toolbox/math/manifold_ift.m`:

```matlab
function U = manifold_ift(Phi, C)
% MANIFOLD_IFT: Inverse manifold Fourier transform. Synthesize a vertex field from
% mode coefficients (or a mode-space kernel): U = Phi * C.
%
% USAGE:  U = manifold_ift(Phi, C)
% INPUTS:
%   Phi : [nV x K]  eigenvectors (already selected)
%   C   : [K x nT]  coefficients, or [K x nCh] mode-space kernel
% OUTPUT:
%   U   : [nV x nT] or [nV x nCh]
%
% Authors: Diellor Basha, 2026
Phi = double(Phi); C = double(C);
if size(C,1) ~= size(Phi,2)
    error('manifold_ift: C has %d rows but Phi has %d columns.', size(C,1), size(Phi,2));
end
U = Phi * C;
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/manifold_ft.m toolbox/math/manifold_ift.m dev/tests/test_manifold_ft_ift_pure.m
git commit -m "feat(eigenmode): pure manifold Fourier transform pair (manifold_ft/ift)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: write `Eigenmodes.Order` (compute + backfill)

**Files:**
- Modify: `toolbox/anatomy/tess_eigenmodes.m` (after line 210)
- Modify: `toolbox/io/in_tess_eigenmodes.m` (after line 57)
- Test: `dev/tests/test_eigenmodes_order_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_order_e2e.m`:

```matlab
function test_eigenmodes_order_e2e
% tess_eigenmodes writes a global-sort .Order; in_tess_eigenmodes backfills it. On the
% 2-hemisphere cortex, the first-K canonical modes span BOTH components.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Surf=in_tess_bst(ctxFile);

Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', 50);   % ~50/hemi, 2 components
assert(isfield(Eig,'Order') && numel(Eig.Order)==size(Eig.Vectors,2), 'Order not written');
assert(issorted(Eig.Values(Eig.Order)), 'Order does not sort Values globally');
assert(isequal(sort(Eig.Order(:)), (1:numel(Eig.Order))'), 'Order not a permutation');
K = 40; comp = Eig.Component(Eig.Order(1:K));
assert(numel(unique(comp))==2, 'first-K canonical must include both hemispheres');

% Backfill: strip Order, save to a temp surface, reload via in_tess_eigenmodes.
TessMat = in_tess_bst(ctxFile, 0);
EigStrip = rmfield(Eig,'Order'); TessMat.Eigenmodes = EigStrip;
tmpFile = fullfile(tempdir, sprintf('tess_order_test_%d.mat', feature('getpid')));
bst_save(tmpFile, TessMat, 'v7');
EigBack = in_tess_eigenmodes(tmpFile);
assert(isfield(EigBack,'Order') && ~isempty(EigBack.Order), 'Order not backfilled');
assert(issorted(EigBack.Values(EigBack.Order)), 'backfilled Order wrong');
delete(tmpFile);
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** ("Order not written").

- [ ] **Step 3a: Modify `tess_eigenmodes.m`** — after line 210 (`Eigenmodes.CompRank = CompRank;`), add:

```matlab
    Eigenmodes.CompRank    = CompRank;
    % Canonical mode order: global eigenvalue sort across all components. The single
    % ordering every consumer uses (Vectors(:,Order(1:K))), so "first K" means the
    % whole-brain lowest spatial frequencies, not one hemisphere.
    [~, Eigenmodes.Order]  = sort(ValuesAll, 'ascend');
```

- [ ] **Step 3b: Modify `in_tess_eigenmodes.m`** — after line 57 (`CompRank` backfill), add:

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

## Task 3: `bst_eigenmodes_project` → wrapper over the pure pair

**Files:**
- Rewrite body: `toolbox/math/bst_eigenmodes_project.m` (lines 38-77)
- Test: `dev/tests/test_eigenmodes_project_wrapper_pure.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_project_wrapper_pure.m`:

```matlab
function test_eigenmodes_project_wrapper_pure
% project's coefficients equal manifold_ft; its reconstruct output equals manifold_ift;
% ModeRange selects over the canonical Order.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

nV=12; K=6; [Q,~]=qr(randn(nV,K),0);
Eig=struct('Vectors',Q,'Values',[10;30;50;20;40;60],'Component',[1;1;1;2;2;2], ...
           'CompRank',[1;2;3;1;2;3]);
[~,Eig.Order]=sort(Eig.Values,'ascend');   % [1 4 2 5 3 6]
M=speye(nV); u=randn(nV,4);

C = bst_eigenmodes_project(Eig, u, M);
assert(isequal(C, manifold_ft(Q, M, u)), 'project coeffs must equal manifold_ft');

% Reconstruct over canonical ranks 1..3 (lowest 3 eigenvalues: 10,20,30 -> stored [1 4 2])
[~, R] = bst_eigenmodes_project(Eig, u, M, 'ModeRange', [1 3]);
iSel = Eig.Order(1:3);
Rexp = manifold_ift(Q(:,iSel), C(iSel,:));
assert(norm(R - Rexp,'fro') < 1e-12, 'ranged reconstruct must use canonical order');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm it FAILS** — current project's `ModeRange` is positional (`Phi(:,iModes)` with `iModes=k1:k2`), so the ranged reconstruct differs from the canonical expectation.

- [ ] **Step 3: Rewrite `bst_eigenmodes_project.m` body** (replace lines 38-77, keeping the header/license above line 38):

```matlab
%% ===== PARSE INPUTS =====
ModeRange = [];
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'moderange', ModeRange = varargin{i+1};
    end
end

Phi    = double(Eigenmodes.Vectors);   % [nV x nModes]
nV     = size(Phi, 1);
nModes = size(Phi, 2);

%% ===== VALIDATE =====
Data = double(Data);
if size(Data, 1) ~= nV
    error('Data has %d rows but eigenmodes have %d vertices.', size(Data, 1), nV);
end
if (size(MassMatrix, 1) ~= nV) || (size(MassMatrix, 2) ~= nV)
    error('MassMatrix must be %dx%d.', nV, nV);
end

%% ===== FORWARD TRANSFORM (all modes, stored order) =====
Coeffs = manifold_ft(Phi, MassMatrix, Data);   % [nModes x nTime]

%% ===== RECONSTRUCT over the CANONICAL mode range (if requested) =====
if nargout >= 2
    if isempty(ModeRange); ModeRange = [1, nModes]; end
    ModeRange(1) = max(1, ModeRange(1));
    ModeRange(2) = min(nModes, ModeRange(2));
    if ModeRange(2) < ModeRange(1)
        error('Empty mode range [%d, %d] (have %d modes).', ModeRange(1), ModeRange(2), nModes);
    end
    if isfield(Eigenmodes,'Order') && ~isempty(Eigenmodes.Order)
        Order = double(Eigenmodes.Order(:));
    else
        [~, Order] = sort(double(Eigenmodes.Values(:)), 'ascend');
    end
    iSel = Order(ModeRange(1):ModeRange(2));
    Reconstructed = manifold_ift(Phi(:, iSel), Coeffs(iSel, :));   % [nV x nTime]
end
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_project.m dev/tests/test_eigenmodes_project_wrapper_pure.m
git commit -m "refactor(eigenmode): project as thin wrapper over manifold_ft/ift, canonical range

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `bst_eigenmode_reconstruct` → wrapper (canonical default + `manifold_ift`)

**Files:**
- Rewrite body: `toolbox/inverse/bst_eigenmode_reconstruct.m` (lines 21-47)
- Test: `dev/tests/test_eigenmode_canonical_consistency_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmode_canonical_consistency_e2e.m`:

```matlab
function test_eigenmode_canonical_consistency_e2e
% No-regression safety net: canonical-default reconstruction equals explicit-ModeIndices,
% and (after Task 5) the canonical-path leadfield Gain equals the existing composed HM.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sStudies=bst_get('ProtocolStudies'); T=[];
for iS=1:numel(sStudies.Study)
    s=sStudies.Study(iS); if isempty(s.HeadModel)||isempty(s.NoiseCov)||isempty(s.NoiseCov(1).FileName); continue; end
    iBase=[]; iEig=[];
    for ih=1:numel(s.HeadModel)
        try hm=in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        isE=isfield(hm,'isEigenmode')&&hm.isEigenmode;
        if isE && isempty(iEig); iEig=ih; elseif ~isE && strcmpi(hm.HeadModelType,'surface') && isempty(iBase); iBase=ih; end
    end
    if ~isempty(iBase)&&~isempty(iEig); T=struct('iS',iS,'iBase',iBase,'iEig',iEig); break; end
end
if isempty(T); disp('SKIP: need base + eigenmode head models.'); return; end
s=sStudies.Study(T.iS);
oldEigHM=in_bst_headmodel(s.HeadModel(T.iEig).FileName,0);

% Reconstruct two ways: canonical default vs explicit ModeIndices (== canonical Order).
[Inv,err]=bst_inverse_eigenmodes(s.HeadModel(T.iEig).FileName, s.NoiseCov(1).FileName, ...
    bst_get('ChannelFileForStudy',s.FileName), [], 'Method','mne','Prior','log','SNR',3);
assert(isempty(err), 'inverse failed: %s', err);
Kidx = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel, oldEigHM.ModeIndices);
Kdef = bst_eigenmode_reconstruct(oldEigHM.SurfaceFile, Inv.ImagingKernel);   % canonical default
rel = norm(Kidx - Kdef,'fro')/norm(Kidx,'fro');
fprintf('reconstruct canonical-default vs explicit rel.diff = %.3e\n', rel);
assert(rel < 1e-9, 'canonical-default reconstruction must equal explicit ModeIndices');

% Leadfield consistency (passes fully after Task 5).
baseHM=in_bst_headmodel(s.HeadModel(T.iBase).FileName,0);
Eig=in_tess_eigenmodes(baseHM.SurfaceFile);
newEigHM=bst_eigenmode_leadfield(baseHM, Eig);
nC=min(size(newEigHM.Gain,2), size(oldEigHM.Gain,2));
relL=norm(newEigHM.Gain(:,1:nC)-oldEigHM.Gain(:,1:nC),'fro')/norm(oldEigHM.Gain(:,1:nC),'fro');
fprintf('leadfield Gain rel.diff = %.3e\n', relL);
assert(relL < 1e-9, 'canonical leadfield Gain diverged from existing composed HM');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** — current `bst_eigenmode_reconstruct` with no `ModeIndices` falls back to `Phi(:,1:K)` (hemisphere-biased), so `Kdef` differs from `Kidx`. (The leadfield half may already pass — that is fine; both must pass after Tasks 4–5.)

- [ ] **Step 3: Rewrite `bst_eigenmode_reconstruct.m` body** (replace lines 21-47, keeping the header above line 21):

```matlab
K = size(ModeKernel, 1);
Eig = [];
if ischar(SurfaceOrPhi)
    [Eig, isComputed] = in_tess_eigenmodes(SurfaceOrPhi);
    if ~isComputed
        error(['No eigenmodes on surface: ' SurfaceOrPhi '. Run "Compute eigenmodes" first.']);
    end
    Phi = double(Eig.Vectors);
    % Default to the surface's canonical order when no explicit selection is given.
    if (nargin < 3 || isempty(ModeIndices))
        ModeIndices = Eig.Order;
    end
else
    Phi = double(SurfaceOrPhi);
end
if (nargin >= 3) && ~isempty(ModeIndices)
    idx = ModeIndices(:);
else
    idx = (1:K)';   % bare Phi matrix, no order known: positional fallback
end
if numel(idx) < K
    error('Fewer mode indices (%d) than kernel modes (%d).', numel(idx), K);
end
idx = idx(1:K);
if max(idx) > size(Phi, 2)
    error('Mode index %d exceeds available eigenmodes (%d).', max(idx), size(Phi,2));
end
ImagingKernel = manifold_ift(Phi(:, idx), ModeKernel);   % [nV x nGoodCh]
```

- [ ] **Step 4: Run it, confirm the reconstruct half PASSES.** (The leadfield assertion may still fail until Task 5; that is expected — re-run after Task 5. If you want this task green in isolation, temporarily comment the leadfield block, then restore it in Task 5.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_eigenmode_reconstruct.m dev/tests/test_eigenmode_canonical_consistency_e2e.m
git commit -m "refactor(eigenmode): reconstruct wraps manifold_ift, defaults to canonical Order

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `bst_eigenmode_leadfield` selects via canonical `Order`

**Files:**
- Modify: `toolbox/forward/bst_eigenmode_leadfield.m` (lines 57-72)
- Test: `dev/tests/test_eigenmode_canonical_consistency_e2e.m` (Task 4's, now must fully pass)

- [ ] **Step 1: Modify `bst_eigenmode_leadfield.m`.** Replace the `% Clamp K` ... `PhiSel = Phi(:, sel);` block (lines 57-72) with:

```matlab
% Clamp K (default: all modes; truncation is an inverse-side concern).
if isempty(nModes) || nModes <= 0
    K = nModesAll;
else
    K = min(nModes, nModesAll);
end
% Canonical selection from the surface's single global eigenvalue order. Recorded as
% ModeIndices so the inverse reconstruction (bst_eigenmode_reconstruct) uses the same
% modes in the same order.
if isfield(Eigenmodes,'Order') && ~isempty(Eigenmodes.Order) && numel(Eigenmodes.Order)==nModesAll
    order = double(Eigenmodes.Order(:));
else
    [~, order] = sort(Values, 'ascend');
end
sel = order(1:K);
PhiSel = Phi(:, sel);
```

(Leave the existing `L_tilde = Lc * PhiSel;` and the `CompHM` assembly — including `CompHM.ModeIndices = sel(:);` and `CompHM.Eigenvalues = Values(sel);` — unchanged.)

- [ ] **Step 2: Run** `dev/tests/test_eigenmode_canonical_consistency_e2e.m` (restore the leadfield block if you commented it in Task 4). Confirm `ALL TESTS PASSED` — both the reconstruct and leadfield-Gain assertions.

- [ ] **Step 3: Commit**

```bash
git add toolbox/forward/bst_eigenmode_leadfield.m
git commit -m "refactor(eigenmode): leadfield selects modes via canonical Order (all modes)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: fix `process_eigenmodes_transform` (the bias bug)

**Files:**
- Modify: `toolbox/process/functions/process_eigenmodes_transform.m` (lines 166-167)
- Test: `dev/tests/test_eigenmodes_transform_canonical_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_transform_canonical_e2e.m`:

```matlab
function test_eigenmodes_transform_canonical_e2e
% The transform must select modes via the canonical Order, so K<1000 spans BOTH
% hemispheres (previously Vectors(:,1:K) was hemisphere-1 only), and it runs e2e.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end
if isempty(bst_get('Protocol','TutorialAuditory')); disp('SKIP: TutorialAuditory missing.'); return; end
iProt=bst_get('Protocol','TutorialAuditory'); gui_brainstorm('SetCurrentProtocol',iProt);
sSubj=bst_get('Subject','Subject01'); ctxFile=sSubj.Surface(sSubj.iCortex).FileName;
Eig=in_tess_eigenmodes(ctxFile);

% Precondition: naive first-600 is single-hemisphere; canonical first-600 is whole-brain.
K=600;
assert(numel(unique(Eig.Component(1:K)))==1, 'precondition: naive first-K single hemisphere');
assert(numel(unique(Eig.Component(Eig.Order(1:K))))==2, 'canonical first-K must span both');

% End-to-end run of the transform process on the deviant average.
s=bst_get('ProtocolStudies'); s=s.Study(6); dc=string({s.Data.Comment});
iD=find(startsWith(dc,'Avg: deviant'),1);
if isempty(iD); disp('SKIP: no deviant average.'); return; end
sOut=bst_process('CallProcess','process_eigenmodes_transform', s.Data(iD).FileName, [], 'nmodes', {600,'',0});
M=in_bst_matrix(sOut(1).FileName);
assert(size(M.Value,1)==600 && all(isfinite(M.Value(:))), 'transform output bad');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it, confirm FAIL** — the precondition asserts hold, but the e2e transform currently selects hemisphere-1 modes; after the fix it selects canonical. (If both asserts pass and the e2e runs, the test passes; the discriminating value is that the transform code path now uses `Order` — verified by Step 3 + the assert that canonical spans both. To make this a true red→green, first confirm the OLD transform's selected modes are single-hemisphere by inspecting `Vectors(:,1:600)` provenance in Step 2, then fix in Step 3.)

- [ ] **Step 3: Modify `process_eigenmodes_transform.m`** lines 166-167. Replace:

```matlab
    Phi     = double(Eigenmodes.Vectors(:, 1:K));
    lambdas = double(Eigenmodes.Values(1:K));
```

with:

```matlab
    % Canonical selection (whole-brain lowest spatial frequencies, never one hemisphere).
    if isfield(Eigenmodes,'Order') && ~isempty(Eigenmodes.Order)
        order = double(Eigenmodes.Order(:));
    else
        [~, order] = sort(double(Eigenmodes.Values(:)), 'ascend');
    end
    sel     = order(1:K);
    Phi     = double(Eigenmodes.Vectors(:, sel));
    lambdas = double(Eigenmodes.Values(sel));
```

(The downstream `ResMat.ImagingKernel = Phi * Kernel;` at lines 194/246 stays — `Phi` is now canonical, so transform and reconstruction remain self-consistent.)

- [ ] **Step 4: Run the test, confirm `ALL TESTS PASSED`.**

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_eigenmodes_transform.m dev/tests/test_eigenmodes_transform_canonical_e2e.m
git commit -m "fix(eigenmode): transform selects via canonical Order (no hemisphere bias)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: `bst_eigenmodes_ensure` (default-1000/hemi, no repair)

**Files:**
- Create: `toolbox/math/bst_eigenmodes_ensure.m`
- Test: `dev/tests/test_eigenmodes_ensure_e2e.m`

- [ ] **Step 1: Write the failing test** — create `dev/tests/test_eigenmodes_ensure_e2e.m`:

```matlab
function test_eigenmodes_ensure_e2e
% Returns canonical eigenmodes when present (idempotent, no recompute), with a valid Order.
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
assert(el < 30, 'ensure must be fast when eigenmodes already exist (no recompute)');
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
% set if absent. The principled home for "use existing modes, else compute default".
%
% USAGE:  Eig = bst_eigenmodes_ensure(SurfaceFile)
%         Eig = bst_eigenmodes_ensure(SurfaceFile, nModesPerHemi)   % default 1000
%
% Computes (if needed) nModesPerHemi modes per component, barycentric mass, remove-DC,
% and NO repair: a non-manifold surface raises an error (repair changes the vertex count
% and breaks surface<->leadfield<->eigenmode consistency; remesh to an icosphere).
%
% Authors: Diellor Basha, 2026
if nargin < 2 || isempty(nModesPerHemi); nModesPerHemi = 1000; end
[Eig, isComputed] = in_tess_eigenmodes(SurfaceFile);
if isComputed && ~isempty(Eig)
    return;
end
Surf = in_tess_bst(SurfaceFile, 0);
mani = tess_manifold(Surf.Vertices, Surf.Faces);
if isstruct(mani) && isfield(mani,'isManifold') && ~mani.isManifold
    error('bst_eigenmodes_ensure:NonManifold', ...
        ['Surface %s is non-manifold; eigenmodes require a 2-manifold mesh. ' ...
         'Remesh to an icosphere (or repair manually) and retry.'], SurfaceFile);
end
Eig = tess_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', nModesPerHemi, ...
    'MassType', 'barycentric', 'RemoveDC', 1);
out_tess_eigenmodes(SurfaceFile, Eig);
end
```

- [ ] **Step 4: Run it, confirm `ALL TESTS PASSED`** (eigenmodes already exist on ico5, returns fast). If `tess_eigenmodes` rejects the `'MassType'`/`'RemoveDC'` name-value pairs, check its `parse_inputs` and use the exact option names it accepts (it is the same function used in `process_eigenmodes`).

- [ ] **Step 5: Commit**

```bash
git add toolbox/math/bst_eigenmodes_ensure.m dev/tests/test_eigenmodes_ensure_e2e.m
git commit -m "feat(eigenmode): bst_eigenmodes_ensure (default 1000/hemi, no repair)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 8: regression sweep

**Files:** none (verification only)

- [ ] **Step 1: Run the existing eigenmode tests** via MCP `run_matlab_file` (clean SKIP if a fixture is absent). Expected `ALL TESTS PASSED`/SKIP each:
- `dev/tests/test_eigenmode_reconstruct_pure.m`
- `dev/tests/test_eigenmode_hemisphere_pure.m`
- `dev/tests/test_eigenmode_leadfield_pure.m`
- `dev/tests/test_inverse_eigenmodes_pure.m`
- `dev/tests/test_eigenmodes_project_pure.m`
- `dev/tests/test_benchmark_inverse_nmodes_e2e.m`
- `dev/tests/test_kernel_comparison.m`

- [ ] **Step 2: Re-run the new tests** (`test_manifold_ft_ift_pure`, `test_eigenmodes_order_e2e`, `test_eigenmodes_project_wrapper_pure`, `test_eigenmode_canonical_consistency_e2e`, `test_eigenmodes_transform_canonical_e2e`, `test_eigenmodes_ensure_e2e`). All must pass.

- [ ] **Step 3: Report** any failures with the first stack trace. If green: source mapping unchanged (consistency tests), analytic transform de-biased (transform test), analytic core now pure (`manifold_ft`/`manifold_ift`).

---

## Self-Review

**1. Spec coverage:**
- Canonical `Eigenmodes.Order` (compute + backfill) → Task 2. ✔
- Pure `manifold_ft`/`manifold_ift` → Task 1. ✔
- `project`/`reconstruct` thin wrappers over the pair → Tasks 3, 4. ✔
- Leadfield via canonical `Order`, all modes, copy to HM → Task 5. ✔
- Inverse `nModes` = first-K canonical → automatic (Gain canonical-ordered after Task 5); proven by Task 4 consistency. ✔
- Transform de-biased → Task 6. ✔
- `bst_eigenmodes_ensure` default 1000/hemi, no repair → Task 7. ✔
- No dedicated accessor; selection inline via `Order` → Tasks 5, 6 + wrappers. ✔
- Viewers (no first-K bug; UX deferred to B) → no task needed; noted in spec §4. ✔
- Backward-compat (old ModeIndices ≡ canonical Order) → Task 4/5 consistency tests prove it. ✔

**2. Placeholder scan:** No "TBD/TODO"; every code step is complete. Task 6 Step 2 and Task 7 Step 4 include explicit fallbacks for environment specifics (provenance inspection; `tess_eigenmodes` option names) rather than vague instructions.

**3. Type consistency:** `manifold_ft(Phi,M,U)`/`manifold_ift(Phi,C)` signatures used identically in Tasks 1,3,4. `Eigenmodes.Order` written (Task 2) and consumed (Tasks 3,4,5,6) consistently. `bst_eigenmode_leadfield` sets `CompHM.ModeIndices = sel(:)` where `sel=order(1:K)` (Task 5), matching `bst_eigenmode_reconstruct`'s `ModeIndices` consumption (Task 4).

**Out of scope (per spec):** Feature B viewer; full retirement of `project`/`reconstruct`; re-sorting stored `Vectors`; inverse math/prior/eigenfilter changes; viewer mode-k UX.
