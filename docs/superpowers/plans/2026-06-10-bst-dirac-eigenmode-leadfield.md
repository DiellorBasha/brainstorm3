# Dirac Eigenmode Leadfield (Phase B — `bst_dirac_eigenmode_leadfield`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compose the **unconstrained** leadfield into the Dirac eigenbasis: embed each per-vertex 3-vector gain as a pure-imaginary quaternion and project onto the stored `DiracEigen`, producing a composed eigenmode head model that the existing `bst_inverse_eigenmodes` consumes unchanged.

**Architecture:** New forward composer `bst_dirac_eigenmode_leadfield` (sibling to `bst_eigenmode_leadfield`). Pure function: takes the base unconstrained head model `Gain [nCh×3·nVert]` and the per-hemisphere `DiracEigen` (from Phase A), and per hemisphere builds `Ψ [4Vₕ×nCh]` (w=0, ambient x/y/z in the imaginary part), projects `L̃ₕ = Ψₕᵀ·Bₕ·Φ_D,ₕ [nCh×K]` (`Bₕ=kron(Mass,I₄)`), and stacks → `CompHM.Gain [nCh×2K]`, `.Eigenvalues [2K×1]`, `.nModes`, `.ModeHemisphere`, `.HemiGlobalVertices`. No nxr call (uses the stored `Mass`).

**Tech Stack:** MATLAB (Brainstorm toolbox), `matlab.unittest`-free script tests (matching `test_eigenmode_leadfield_pure`), run via the MATLAB MCP. Phase A (`tess_dirac_eigenmodes` → `TessMat.DiracEigen`) is merged.

**Reference facts (grounded):**
- Scalar sibling `toolbox/forward/bst_eigenmode_leadfield.m`: constrains `Lc=bst_gain_orient(Gain,GridOrient)`, composes `L̃=Lc·Φ [nCh×K]`, builds `CompHM=HeadModel` with `Gain=L̃`, `GridLoc/Orient/Atlas=[]`, `isEigenmode=1`, `nModes`, `Eigenvalues`, `ModeIndices`, `HeadModelType='surface'`, `Comment`. The Dirac version does NOT constrain — it consumes the unconstrained 3-vector gain directly.
- `DiracEigen(hh)` (Phase A) = `{Vectors[4Vₕ×K], Values[K×1], Mass[Vₕ×Vₕ], nModes, Order, Tau, GlobalVertices, Hemisphere, Provenance}`; `Vectors` is B-orthonormal (`B=kron(Mass,I₄)`), vertex-interleaved `4v+c`, order `[w,x,y,z]`.
- Unconstrained `Gain [nCh×3·nVert]`: source `s` occupies columns `3*(s-1)+(1:3)` (x,y,z, ambient world coords).
- `bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg)` → `Kernel [K×nCh]` (here `K=2K`). It reads `HM.Gain [nAllCh×K]`, `HM.Eigenvalues`.
- Test style: synthetic script function ending `disp('ALL TESTS PASSED')`, run by calling the function name (NOT `runtests`).

**MATLAB-session discipline:** before each run, `rehash; clear <names>;`. **Never** a bare `clear` (wipes Brainstorm `GlobalData`). These tests are fully synthetic (no eigs, no real cortex) → fast.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/forward/bst_dirac_eigenmode_leadfield.m` | **Create.** Unconstrained-leadfield → Dirac-eigenbasis composer. |
| `dev/tests/test_dirac_eigenmode_leadfield_pure.m` | **Create.** Synthetic composition + projection-identity + error tests. |
| `dev/tests/test_dirac_eigenmode_leadfield_inverse.m` | **Create.** Integration: composed model through `bst_inverse_eigenmodes('SolvePure',…)`. |

---

## Task 1: `bst_dirac_eigenmode_leadfield` + pure test

**Files:**
- Create: `toolbox/forward/bst_dirac_eigenmode_leadfield.m`
- Test: `dev/tests/test_dirac_eigenmode_leadfield_pure.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_dirac_eigenmode_leadfield_pure.m`:

```matlab
function test_dirac_eigenmode_leadfield_pure
% Verify the Dirac forward composer (synthetic, no eigs/cortex):
%   - composed Gain is [nCh x 2K] and equals the per-hemisphere quaternion projection
%   - embedding places the ambient x/y/z gain in the imaginary quaternion slots (w=0)
%   - Eigenvalues / ModeHemisphere / nModes carried through; basis stacked L then R
%   - non-surface and non-unconstrained head models are rejected
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 6; vH1 = [1 2 3 4]; vH2 = [5 6 7]; nVert = 7; K = 3;

% deterministic unconstrained gain [nCh x 3*nVert]
GainU = zeros(nCh, 3*nVert);
for i = 1:nCh
    for j = 1:3*nVert
        GainU(i,j) = sin(i*j*pi/23);
    end
end

% synthetic per-hemisphere DiracEigen with Mass = I (so B = I; Vectors orthonormal)
hemis = {vH1(:), vH2(:)}; tags = {'L','R'};
DE = struct('Vectors',{[]},'Values',{[]},'Mass',{[]},'nModes',{[]},'Order',{[]}, ...
            'Tau',{[]},'GlobalVertices',{[]},'Hemisphere',{[]},'Provenance',{[]});
for hh = 1:2
    nVh = numel(hemis{hh});
    raw = reshape(cos((1:(4*nVh*(K+2)))*0.7), 4*nVh, K+2);   % deterministic, full rank
    Phi = orth(raw); Phi = Phi(:, 1:K);                      % orthonormal columns (B=I)
    DE(hh).Vectors        = Phi;
    DE(hh).Values         = (1:K)' + 0.1*hh;
    DE(hh).Mass           = speye(nVh);
    DE(hh).nModes         = K;
    DE(hh).Order          = (1:K)';
    DE(hh).Tau            = 0.5;
    DE(hh).GlobalVertices = hemis{hh};
    DE(hh).Hemisphere     = tags{hh};
    DE(hh).Provenance     = struct('Backend','nxr');
end

HeadModel = struct('Gain', GainU, 'GridLoc', zeros(nVert,3), 'GridOrient', zeros(nVert,3), ...
    'GridAtlas', [], 'HeadModelType', 'surface', 'SurfaceFile', 'tess_cortex_test.mat', ...
    'Comment', 'OS-MEG test');

CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DE, 'nModes', K);

% --- shape + carried metadata ---
assert(isequal(size(CompHM.Gain), [nCh, 2*K]), 'Composed Gain must be [nCh x 2K].');
assert(CompHM.nModes == 2*K, 'nModes must be 2K.');
assert(isequal(CompHM.Eigenvalues(:), [DE(1).Values; DE(2).Values]), 'Eigenvalues stacked L then R.');
assert(isequal(CompHM.ModeHemisphere(:), [ones(K,1); 2*ones(K,1)]), 'ModeHemisphere L then R.');
assert(CompHM.isDiracEigenmode == 1, 'isDiracEigenmode flag set.');
assert(isempty(CompHM.GridLoc) && isempty(CompHM.GridOrient), 'GridLoc/Orient cleared.');
assert(strcmp(CompHM.SurfaceFile, 'tess_cortex_test.mat'), 'SurfaceFile carried through.');

% --- projection identity with an INDEPENDENT explicit embedding (oracle) ---
for hh = 1:2
    vH = DE(hh).GlobalVertices; nVh = numel(vH);
    Psi = zeros(4*nVh, nCh);
    for vloc = 1:nVh
        s = vH(vloc);
        g = GainU(:, 3*(s-1)+(1:3));     % [nCh x 3] ambient
        Psi(4*(vloc-1)+1, :) = 0;        % w
        Psi(4*(vloc-1)+2, :) = g(:,1)';  % x
        Psi(4*(vloc-1)+3, :) = g(:,2)';  % y
        Psi(4*(vloc-1)+4, :) = g(:,3)';  % z
    end
    B = kron(DE(hh).Mass, speye(4));
    Lref = Psi' * (B * DE(hh).Vectors);  % [nCh x K]
    block = CompHM.Gain(:, (hh-1)*K + (1:K));
    assert(max(abs(block(:) - Lref(:))) < 1e-9, sprintf('Hemisphere %d projection mismatch.', hh));
end

% --- errors ---
ok = false;
HMvol = HeadModel; HMvol.HeadModelType = 'volume';
try, bst_dirac_eigenmode_leadfield(HMvol, DE); catch ME, ok = strcmp(ME.identifier,'bst_dirac_eigenmode_leadfield:NotSurface'); end
assert(ok, 'Volume head model must raise NotSurface.');

ok = false;
HMbad = HeadModel; HMbad.Gain = GainU(:, 1:end-1);   % cols not a multiple of 3
try, bst_dirac_eigenmode_leadfield(HMbad, DE); catch ME, ok = strcmp(ME.identifier,'bst_dirac_eigenmode_leadfield:NotUnconstrained'); end
assert(ok, 'Non-unconstrained gain must raise NotUnconstrained.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_dirac_eigenmode_leadfield_pure bst_dirac_eigenmode_leadfield; disp('rehashed');
test_dirac_eigenmode_leadfield_pure
```
Expected: ERROR `Undefined function 'bst_dirac_eigenmode_leadfield'`.

- [ ] **Step 3: Write `toolbox/forward/bst_dirac_eigenmode_leadfield.m`**

```matlab
function CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen, varargin)
% BST_DIRAC_EIGENMODE_LEADFIELD: Compose the UNCONSTRAINED leadfield into the Dirac eigenbasis.
%
% USAGE:  CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen, 'nModes', K)
%
% DESCRIPTION:
%     Forward step for Dirac eigenmode source mapping. Unlike the scalar LBO
%     composer (which constrains the leadfield to the surface normal), this keeps
%     the full unconstrained 3-vector gain and expands it in the curvature-aware
%     Dirac eigenbasis. Each per-vertex gain 3-vector (ambient/world coords) is
%     embedded as a pure-imaginary quaternion psi = [0, gx, gy, gz] and projected
%     onto the per-hemisphere Dirac eigenvectors:
%         L~_h = Psi_h' * B_h * Phi_h      [nCh x K],   B_h = kron(Mass_h, I4)
%     The two hemispheres are stacked into a composed head model
%         Gain = [L~_L, L~_R]              [nCh x 2K]
%     consumed unchanged by bst_inverse_eigenmodes.
%
% INPUTS:
%     HeadModel  : base UNCONSTRAINED surface head model (in_bst_headmodel, ApplyOrient=0):
%                  .Gain [nCh x 3*nVert], .HeadModelType, .SurfaceFile, .Comment
%     DiracEigen : 1x2 per-hemisphere struct (TessMat.DiracEigen, from tess_dirac_eigenmodes):
%                  .Vectors [4Vh x K], .Values [K x 1], .Mass [Vh x Vh],
%                  .nModes, .Tau, .GlobalVertices
% OPTIONS:
%     'nModes' : modes per hemisphere to keep (default all; clamped to available)
%
% OUTPUT:
%     CompHM : composed head-model struct (Gain [nCh x 2K], Eigenvalues [2K x 1],
%              nModes=2K, ModeHemisphere [2K x 1], HemiGlobalVertices {L,R},
%              isEigenmode/isDiracEigenmode flags) ready to save / feed the inverse.
%
% Authors: Diellor Basha, 2026

    nModes = [];
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'nmodes', nModes = varargin{i+1};
        end
    end

    % surface head model only
    if isfield(HeadModel,'HeadModelType') && ~isempty(HeadModel.HeadModelType) ...
            && ~strcmpi(HeadModel.HeadModelType,'surface')
        error('bst_dirac_eigenmode_leadfield:NotSurface', ...
            'Dirac eigenmode leadfield requires a surface head model (got ''%s'').', HeadModel.HeadModelType);
    end

    G = double(HeadModel.Gain);          % [nCh x 3*nVert] unconstrained
    nCh = size(G,1);
    if mod(size(G,2), 3) ~= 0
        error('bst_dirac_eigenmode_leadfield:NotUnconstrained', ...
            'Dirac leadfield requires an unconstrained head model: Gain must be [nCh x 3*nVert].');
    end

    if numel(DiracEigen) ~= 2
        error('bst_dirac_eigenmode_leadfield:badBasis', ...
            'DiracEigen must be a 1x2 per-hemisphere struct array (from tess_dirac_eigenmodes).');
    end

    Kfull = min([DiracEigen.nModes]);
    if isempty(nModes) || nModes <= 0
        K = Kfull;
    else
        K = min(nModes, Kfull);
    end

    Lblk = cell(1,2); vblk = cell(1,2); hblk = cell(1,2);
    for hh = 1:2
        vH  = DiracEigen(hh).GlobalVertices(:);
        nVh = numel(vH);
        Phi = double(DiracEigen(hh).Vectors);
        if size(Phi,1) ~= 4*nVh
            error('bst_dirac_eigenmode_leadfield:shapeMismatch', ...
                'Hemisphere %d: Vectors has %d rows, expected 4*nV=%d.', hh, size(Phi,1), 4*nVh);
        end
        Phi  = Phi(:, 1:K);
        Vals = double(DiracEigen(hh).Values(:)); Vals = Vals(1:K);

        % embed unconstrained gain as a pure-imaginary quaternion field (w=0)
        Psi = zeros(4*nVh, nCh);
        Psi(2:4:end, :) = G(:, (vH-1)*3 + 1).';   % i (x)
        Psi(3:4:end, :) = G(:, (vH-1)*3 + 2).';   % j (y)
        Psi(4:4:end, :) = G(:, (vH-1)*3 + 3).';   % k (z)   (w rows 1:4:end stay 0)

        B = kron(DiracEigen(hh).Mass, speye(4));   % [4Vh x 4Vh]
        Lblk{hh} = Psi' * (B * Phi);               % [nCh x K]
        vblk{hh} = Vals;
        hblk{hh} = hh * ones(K,1);
    end

    CompHM = HeadModel;
    CompHM.Gain               = [Lblk{1}, Lblk{2}];     % [nCh x 2K]
    CompHM.GridLoc            = [];
    CompHM.GridOrient         = [];
    CompHM.GridAtlas          = [];
    CompHM.isEigenmode        = 1;
    CompHM.isDiracEigenmode   = 1;
    CompHM.nModes             = 2*K;
    CompHM.Eigenvalues        = [vblk{1}; vblk{2}];     % [2K x 1]
    CompHM.ModeHemisphere     = [hblk{1}; hblk{2}];     % [2K x 1] hemisphere index per mode
    CompHM.HemiGlobalVertices = {DiracEigen(1).GlobalVertices(:), DiracEigen(2).GlobalVertices(:)};
    CompHM.HeadModelType      = 'surface';
    CompHM.Comment            = sprintf('Dirac eigenmode leadfield (%d modes, tau=%.3g) | %s', ...
        2*K, DiracEigen(1).Tau, local_default(HeadModel,'Comment',''));
end

function v = local_default(s, f, d)
    if isfield(s, f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_dirac_eigenmode_leadfield_pure bst_dirac_eigenmode_leadfield; disp('rehashed');
test_dirac_eigenmode_leadfield_pure
```
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/forward/bst_dirac_eigenmode_leadfield.m dev/tests/test_dirac_eigenmode_leadfield_pure.m
git commit -m "feat(bst-dirac-eigenmode-leadfield): unconstrained leadfield -> Dirac eigenbasis composer"
```

---

## Task 2: Integration through `bst_inverse_eigenmodes`

**Files:**
- Create: `dev/tests/test_dirac_eigenmode_leadfield_inverse.m`

- [ ] **Step 1: Write the integration test**

Create `dev/tests/test_dirac_eigenmode_leadfield_inverse.m`:

```matlab
function test_dirac_eigenmode_leadfield_inverse
% The composed Dirac head model must be consumed unchanged by the existing
% mode-space inverse: bst_inverse_eigenmodes('SolvePure', ...) -> Kernel [2K x nCh].
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 8; vH1 = (1:5)'; vH2 = (6:9)'; nVert = 9; K = 3;
GainU = zeros(nCh, 3*nVert);
for i = 1:nCh
    for j = 1:3*nVert
        GainU(i,j) = cos(i*j*pi/29);
    end
end
hemis = {vH1, vH2}; tags = {'L','R'};
DE = struct('Vectors',{[]},'Values',{[]},'Mass',{[]},'nModes',{[]},'Order',{[]}, ...
            'Tau',{[]},'GlobalVertices',{[]},'Hemisphere',{[]},'Provenance',{[]});
for hh = 1:2
    nVh = numel(hemis{hh});
    Phi = orth(reshape(sin((1:(4*nVh*(K+2)))*0.9), 4*nVh, K+2)); Phi = Phi(:,1:K);
    DE(hh).Vectors=Phi; DE(hh).Values=(1:K)'+0.1*hh; DE(hh).Mass=speye(nVh);
    DE(hh).nModes=K; DE(hh).Order=(1:K)'; DE(hh).Tau=0.5;
    DE(hh).GlobalVertices=hemis{hh}; DE(hh).Hemisphere=tags{hh}; DE(hh).Provenance=struct('Backend','nxr');
end
HeadModel = struct('Gain',GainU,'GridLoc',zeros(nVert,3),'GridOrient',zeros(nVert,3), ...
    'GridAtlas',[],'HeadModelType','surface','SurfaceFile','tess_cortex_test.mat','Comment','test');

CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DE, 'nModes', K);

% Feed the composed leadfield into the mode-space inverse (pure-math entry).
iW   = eye(nCh);                 % trivial whitener
Proj = eye(nCh);                 % no SSP
Kernel = bst_inverse_eigenmodes('SolvePure', CompHM.Gain, CompHM.Eigenvalues, ...
    iW, Proj, 'mne', 'log', 1, 3, false);

assert(isequal(size(Kernel), [CompHM.nModes, nCh]), 'Kernel must be [2K x nChannels].');
assert(all(isfinite(Kernel(:))), 'Kernel must be finite.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run to verify it passes**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_dirac_eigenmode_leadfield_inverse bst_dirac_eigenmode_leadfield; disp('rehashed');
test_dirac_eigenmode_leadfield_inverse
```
Expected: `ALL TESTS PASSED`. (If the `SolvePure` argument order differs in the installed `bst_inverse_eigenmodes`, fix the TEST call to match the real signature — `[L_tilde, lambdas, iW, Proj, Method, Prior, Alpha, SNR, Unreg]` — and note it.)

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_dirac_eigenmode_leadfield_inverse.m
git commit -m "test(bst-dirac-eigenmode-leadfield): integration through bst_inverse_eigenmodes"
```

---

## Final verification

- [ ] **Run both tests; expect both `ALL TESTS PASSED`.**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_dirac_eigenmode_leadfield_pure test_dirac_eigenmode_leadfield_inverse bst_dirac_eigenmode_leadfield; disp('rehashed');
test_dirac_eigenmode_leadfield_pure
test_dirac_eigenmode_leadfield_inverse
```

- [ ] **Then complete via superpowers:finishing-a-development-branch.**

---

## Notes for the implementer

- **Unconstrained, not constrained:** unlike `bst_eigenmode_leadfield`, do NOT call `bst_gain_orient`. The whole point is to expand the full 3-vector gain in the vector (Dirac) basis.
- **Embedding convention:** `Ψ` rows `4(v-1)+{2,3,4}` carry ambient `{gx,gy,gz}`; the `w` row `4(v-1)+1` stays 0. Component order `[w,x,y,z]` must match the Phase A `DiracEigen` layout.
- **No nxr call:** `B` is rebuilt from the stored `Mass` (`kron(Mass,I₄)`), so Phase B is pure linear algebra.
- **Reconstruction is downstream:** `ModeHemisphere` + `HemiGlobalVertices` are stored so a later step can map mode coefficients back via `imag(Φ_D·coeffs)` to per-vertex 3-vector currents (`nComponents=3`). That results/visualization wiring is out of scope here.
- **Caller wiring** (not in this plan): a process loads the base unconstrained head model + `TessMat.DiracEigen` (computing it via `tess_dirac_eigenmodes` if absent) and saves `CompHM` as a `headmodel_*.mat`, exactly as the scalar `process_eigenmode_leadfield` does.
