# Connection Laplacian (M1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `tess_connection_laplacian`, a pure `(Vertices, Faces)` function returning the complex-Hermitian vertex connection Laplacian (Levi-Civita, n-RoSy bundle) and its companion lumped mass, via the nxr-compute plugin.

**Architecture:** Thin MATLAB wrapper over `nxr.manifold.operator.connectionLaplacian` (requested in `'complex'` format) + `nxr.manifold.context`. nxr-only (no MATLAB fallback), mirroring `tess_tangents`. The operator is intrinsic to the mesh; phase readout, registration alignment, eigenmodes, and storage are out of scope (later milestones). Spectral sanity is checked with MATLAB's native `eigs` because nxr's `solve` rejects complex matrices.

**Tech Stack:** MATLAB, Brainstorm (`in_tess_bst`, `bst_plugin`, `bst_get`, `tess_manifold`, `tess_sphere`), nxr-compute MEX plugin (`+nxr` package). Tests are function-style scripts under `dev/tests/`, run via the MATLAB MCP.

**Reference spec:** `dev/connection_laplacian_integration.md`

---

## File Structure

| Path | Responsibility |
|---|---|
| `toolbox/anatomy/tess_connection_laplacian.m` | The operator function + a local `nxr_is_loaded` subfunction (copy of `tess_laplacian`'s; no shared-helper refactor). |
| `dev/tests/test_connection_laplacian_smoke.m` | Assemble on real `cortex_20484V`; assert shape/complex/Hermitian/mass/Info. |
| `dev/tests/test_connection_laplacian_spectrum.m` | Smallest-`k` modes via MATLAB `eigs`; eigenvalues real & ≥ 0. |
| `dev/tests/test_connection_laplacian_nsym.m` | `nSym = 1,2,4` each assemble and stay Hermitian. |
| `dev/tests/test_connection_laplacian_guard.m` | nxr unloaded → `nxrNotLoaded` error (reloads in cleanup). |

`tess_laplacian.m` is **not** modified.

---

## Task 1: The operator function + smoke test

**Files:**
- Create: `toolbox/anatomy/tess_connection_laplacian.m`
- Test: `dev/tests/test_connection_laplacian_smoke.m`

- [ ] **Step 1: Write the failing smoke test**

Create `dev/tests/test_connection_laplacian_smoke.m`:

```matlab
function test_connection_laplacian_smoke
% Smoke test for tess_connection_laplacian on the real TutorialAnatomy cortex
% (cortex_20484V): the complex-Hermitian vertex connection Laplacian, its lumped
% mass, and the Info struct. Pure (V,F) operator — no storage, no readout.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Require + load nxr-compute (operator is nxr-only) ---
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required for this test: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

% --- Resolve the real cortex_20484V from the current protocol ---
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: cortex_20484V not found in the current protocol (load TutorialAnatomy).\n');
    return;
end
fprintf('Source cortex: %s\n', SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

% --- Assemble (defaults: vertex domain, nSym=1, complex) ---
[K, M, Info] = tess_connection_laplacian(V, F);

% --- Operator shape / type / symmetry ---
assert(isequal(size(K), [nV nV]), 'K must be nV x nV.');
assert(issparse(K), 'K must be sparse.');
assert(~isreal(K), 'K must be complex (nonzero imaginary part on a curved cortex).');
assert(max(max(abs(K - K'))) < 1e-9, 'K must be Hermitian (conjugate-symmetric).');

% --- Mass ---
assert(isequal(size(M), [nV nV]), 'M must be nV x nV.');
assert(isdiag(M), 'M must be diagonal (lumped vertex mass).');
assert(all(diag(M) > 0), 'M diagonal entries must be positive.');

% --- Info ---
assert(Info.nSym == 1, 'Info.nSym should default to 1.');
assert(strcmp(Info.Domain, 'vertex'), 'Info.Domain should be ''vertex''.');
assert(strcmp(Info.Format, 'complex'), 'Info.Format should be ''complex''.');
assert(strcmp(Info.Backend, 'nxr'), 'Info.Backend should be ''nxr''.');
assert(Info.baseDim == nV, 'Info.baseDim should equal nV.');

fprintf('PASSED: connection Laplacian assembled (nV=%d): sparse, complex, Hermitian; lumped mass; Info OK.\n', nV);
fprintf('ALL TESTS PASSED: test_connection_laplacian_smoke\n');
end


function SurfaceFile = find_cortex_20484V()
% Return the FileName of a cortex surface whose name contains 'cortex_20484V'
% in the current protocol, or '' if none (e.g. wrong protocol loaded).
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if strcmpi(surf(iF).SurfaceType, 'Cortex') && ...
           ~isempty(strfind(lower(surf(iF).FileName), 'cortex_20484v'))
            SurfaceFile = surf(iF).FileName;
            return;
        end
    end
end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_connection_laplacian_smoke.m`
Expected: FAIL — `Undefined function or variable 'tess_connection_laplacian'`.

- [ ] **Step 3: Write the function**

Create `toolbox/anatomy/tess_connection_laplacian.m`:

```matlab
function [K, M, Info] = tess_connection_laplacian(Vertices, Faces, varargin)
% TESS_CONNECTION_LAPLACIAN: Complex-Hermitian vertex connection Laplacian (nxr).
%
% USAGE:  [K, M, Info] = tess_connection_laplacian(Vertices, Faces)
%         [K, M, Info] = tess_connection_laplacian(Vertices, Faces, 'nSym', 1)
%
% DESCRIPTION:
%     Assembles the connection Laplacian on the n-RoSy tangent bundle of a
%     triangle mesh, using nxr-compute (geometry-central). The connection is the
%     intrinsic discrete Levi-Civita connection: transport rotations along each
%     halfedge are raised to the nSym power. The operator is complex Hermitian and
%     positive (semi-)definite; its eigenvalues are intrinsic to the mesh and its
%     eigenvector phase is gauge-dependent (a globally consistent reference frame
%     is required to read out comparable phase — a later milestone).
%
%     This is a PURE OPERATOR: it depends on (Vertices, Faces) alone and does not
%     consume or store any tangent frame. nxr is REQUIRED (no MATLAB fallback);
%     the operator needs geometry-central's Levi-Civita transport vectors.
%
% INPUTS:
%     Vertices : [nV x 3] vertex positions.
%     Faces    : [nF x 3] triangle vertex indices (1-based, Brainstorm convention).
%
% OPTIONS:
%     nSym           : Bundle symmetry (default 1). 1 = true vector field (carries
%                      phase); 2 = line field; 4 = cross field.
%     Domain         : 'vertex' (default) | 'face' | 'edge'. Vertex is the target;
%                      M is returned only for the vertex domain.
%     Regularization : Diagonal epsilon added by nxr for strict positive-
%                      definiteness (default 1e-8).
%     CheckManifold  : (logical) Pre-check the mesh with tess_manifold for a
%                      friendlier error (default false).
%
% OUTPUTS:
%     K    : [N x N] complex Hermitian sparse connection Laplacian (N = nV for the
%            vertex domain), assembled as K_real + 1i*K_imag.
%     M    : [N x N] real diagonal lumped vertex mass (area/3 per vertex) for the
%            vertex domain; [] for 'face'/'edge' (masses for those domains are a
%            later-milestone concern).
%     Info : struct with fields nSym, Domain, Regularization, baseDim, Format
%            ('complex'), Backend ('nxr').
%
% NOTE: nxr's eigensolver ('solve') is real-only and cannot consume this complex
%     K. Solve the generalized problem K*phi = lambda*M*phi with MATLAB eigs.
%
% SEE ALSO: tess_laplacian, tess_tangents, tess_manifold

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

%% ===== PARSE INPUTS =====
nSym           = 1;
Domain         = 'vertex';
Regularization = 1e-8;
CheckManifold  = false;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'nsym',           nSym = varargin{i+1};
        case 'domain',         Domain = lower(varargin{i+1});
        case 'regularization', Regularization = varargin{i+1};
        case 'checkmanifold',  CheckManifold = varargin{i+1};
    end
end

% Validate
if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Vertices and Faces must have 3 columns.');
end
if ~ismember(Domain, {'vertex', 'face', 'edge'})
    error('tess_connection_laplacian:badDomain', ...
        'Domain must be ''vertex'', ''face'', or ''edge'' (got ''%s'').', Domain);
end

%% ===== REQUIRE NXR (no MATLAB fallback) =====
if ~nxr_is_loaded()
    error('tess_connection_laplacian:nxrNotLoaded', ...
        ['The connection Laplacian requires the nxr-compute plugin (no MATLAB ' ...
         'fallback). Install/load it via bst_plugin(''Install'',''nxr-compute'').']);
end

%% ===== OPTIONAL MANIFOLD CHECK =====
if CheckManifold
    [~, ~, isManifold] = tess_manifold(Vertices, double(Faces), 'Repair', 0, 'Verbose', 0);
    if ~isManifold
        error('tess_connection_laplacian:NonManifold', ...
            ['Input mesh is not a clean 2-manifold; the connection Laplacian ' ...
             'requires one. Validate/repair with tess_manifold.']);
    end
end

%% ===== ASSEMBLE (complex Hermitian) =====
% Faces are passed 1-based; the nxr MEX marshalling subtracts 1 (proven to
% machine precision by the tess_laplacian parity test).
mctx = nxr.manifold.context(Vertices, Faces);
opts = struct('domain', Domain, 'nSym', nSym, ...
              'regularization', Regularization, 'format', 'complex');
CL = nxr.manifold.operator.connectionLaplacian(mctx, opts);

% nxr returns the complex Hermitian operator as parallel real sparse parts.
K = CL.K_real + 1i * CL.K_imag;
% Clear floating-point asymmetry so downstream eigs() detects a Hermitian
% operator and returns real eigenvalues.
K = (K + K') / 2;

%% ===== MASS (vertex domain only) =====
% Lumped vertex mass (area/3 per vertex) = ctx.M. For the complex bundle the
% inner product is <u,v> = sum_i area_i * conj(u_i) * v_i, so the mass is the
% ordinary real diagonal lumped mass (no block duplication for complex format).
% Face/edge-domain masses are a later-milestone concern.
if strcmp(Domain, 'vertex')
    M = mctx.M;
else
    M = [];
end

%% ===== INFO =====
Info = struct();
Info.nSym           = nSym;
Info.Domain         = Domain;
Info.Regularization = Regularization;
Info.baseDim        = CL.baseDim;
Info.Format         = 'complex';
Info.Backend        = 'nxr';

end


function tf = nxr_is_loaded()
% True only if the nxr-compute plugin is currently loaded (cheap, in-memory).
% Local copy of the same guard in tess_laplacian (intentionally not shared).
tf = false;
try
    PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
    tf = ~isempty(PlugDesc) && isfield(PlugDesc, 'isLoaded') && PlugDesc.isLoaded;
catch
    tf = false;   % bst_plugin unavailable / any error -> treat as not loaded
end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_connection_laplacian_smoke.m`
Expected: PASS — prints `ALL TESTS PASSED: test_connection_laplacian_smoke` (or `SKIP:` if the cortex isn't in the loaded protocol — in that case load TutorialAnatomy and re-run).

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_connection_laplacian.m dev/tests/test_connection_laplacian_smoke.m
git commit -m "feat(connection-laplacian): tess_connection_laplacian operator + smoke test"
```

---

## Task 2: Spectral sanity test

**Files:**
- Test: `dev/tests/test_connection_laplacian_spectrum.m`

- [ ] **Step 1: Write the test**

Create `dev/tests/test_connection_laplacian_spectrum.m`:

```matlab
function test_connection_laplacian_spectrum
% Spectral sanity for tess_connection_laplacian on the real cortex_20484V.
% nxr's 'solve' is real-only and cannot consume the complex K, so the
% generalized problem K*phi = lambda*M*phi is solved with MATLAB eigs (which
% handles a Hermitian K and a Hermitian-PD M). This is a light operator-level
% sanity check; full spectral validation belongs to the eigenmode milestone.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required for this test: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: cortex_20484V not found in the current protocol (load TutorialAnatomy).\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);

[K, M] = tess_connection_laplacian(V, F);

% Smallest k modes of the complex Hermitian generalized eigenproblem.
% (The whole-mesh operator is block-diagonal across the two hemispheres.)
k = 6;
[Phi, Lam] = eigs(K, M, k, 'smallestabs');
lam = diag(Lam);

assert(isequal(size(Phi), [size(K,1) k]), 'Phi must be N x k.');
assert(max(abs(imag(lam))) < 1e-6 * max(abs(real(lam)) + 1), ...
    'Eigenvalues of a Hermitian operator must be real.');
lam = sort(real(lam), 'ascend');
assert(all(lam > -1e-6), 'Connection-Laplacian eigenvalues must be >= 0.');

fprintf('PASSED: %d smallest eigenvalues real and >= 0; range [%.3g, %.3g].\n', ...
    k, lam(1), lam(end));
fprintf('ALL TESTS PASSED: test_connection_laplacian_spectrum\n');
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if strcmpi(surf(iF).SurfaceType, 'Cortex') && ...
           ~isempty(strfind(lower(surf(iF).FileName), 'cortex_20484v'))
            SurfaceFile = surf(iF).FileName;
            return;
        end
    end
end
end
```

- [ ] **Step 2: Run the test to verify it passes**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_connection_laplacian_spectrum.m`
Expected: PASS — prints eigenvalue range and `ALL TESTS PASSED` (the function from Task 1 already provides `K`/`M`; this test validates its spectrum). `SKIP:` if the cortex isn't present.

- [ ] **Step 3: Commit**

```bash
git add dev/tests/test_connection_laplacian_spectrum.m
git commit -m "test(connection-laplacian): spectral sanity via MATLAB eigs"
```

---

## Task 3: nSym-variants test

**Files:**
- Test: `dev/tests/test_connection_laplacian_nsym.m`

- [ ] **Step 1: Write the test**

Create `dev/tests/test_connection_laplacian_nsym.m`:

```matlab
function test_connection_laplacian_nsym
% tess_connection_laplacian assembles for each n-RoSy symmetry (vector/line/cross)
% and the operator stays complex Hermitian. Real cortex_20484V.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required for this test: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: cortex_20484V not found in the current protocol (load TutorialAnatomy).\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

for nSym = [1 2 4]
    [K, ~, Info] = tess_connection_laplacian(V, F, 'nSym', nSym);
    assert(Info.nSym == nSym, 'Info.nSym should be %d.', nSym);
    assert(isequal(size(K), [nV nV]), 'K must be nV x nV for nSym=%d.', nSym);
    assert(max(max(abs(K - K'))) < 1e-9, 'K must be Hermitian for nSym=%d.', nSym);
    fprintf('PASSED: nSym=%d assembled, Hermitian.\n', nSym);
end

fprintf('ALL TESTS PASSED: test_connection_laplacian_nsym\n');
end


function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if strcmpi(surf(iF).SurfaceType, 'Cortex') && ...
           ~isempty(strfind(lower(surf(iF).FileName), 'cortex_20484v'))
            SurfaceFile = surf(iF).FileName;
            return;
        end
    end
end
end
```

- [ ] **Step 2: Run the test to verify it passes**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_connection_laplacian_nsym.m`
Expected: PASS — three `PASSED: nSym=...` lines then `ALL TESTS PASSED`. `SKIP:` if the cortex isn't present.

- [ ] **Step 3: Commit**

```bash
git add dev/tests/test_connection_laplacian_nsym.m
git commit -m "test(connection-laplacian): nSym=1/2/4 assemble and stay Hermitian"
```

---

## Task 4: Backend-guard test

**Files:**
- Test: `dev/tests/test_connection_laplacian_guard.m`

- [ ] **Step 1: Write the test**

Create `dev/tests/test_connection_laplacian_guard.m`:

```matlab
function test_connection_laplacian_guard
% tess_connection_laplacian must error cleanly (nxrNotLoaded) when the nxr-compute
% plugin is not loaded, since there is no MATLAB fallback. Uses a synthetic sphere
% (the guard fires before any assembly). Restores the plugin's loaded state after.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Record whether nxr is installed/loaded so we can restore it afterward.
PlugDesc = bst_plugin('GetInstalled', 'nxr-compute');
wasLoaded = ~isempty(PlugDesc) && isfield(PlugDesc, 'isLoaded') && PlugDesc.isLoaded;
if isempty(PlugDesc)
    fprintf('SKIP: nxr-compute not installed; cannot exercise the loaded/unloaded guard.\n');
    return;
end
% Ensure it is reloaded at the end regardless of assertion outcome.
restore = onCleanup(@() bst_plugin('Load', 'nxr-compute'));

% Unload so the guard precondition (not loaded) holds.
if wasLoaded
    bst_plugin('Unload', 'nxr-compute');
end

[V, F] = tess_sphere(162);
threw = false;
try
    tess_connection_laplacian(V, F);
catch ME
    threw = strcmp(ME.identifier, 'tess_connection_laplacian:nxrNotLoaded');
    if ~threw
        rethrow(ME);
    end
end
assert(threw, 'tess_connection_laplacian did not raise nxrNotLoaded when the plugin was unloaded.');

fprintf('PASSED: errors with nxrNotLoaded when nxr-compute is not loaded.\n');
fprintf('ALL TESTS PASSED: test_connection_laplacian_guard\n');
end
```

- [ ] **Step 2: Run the test to verify it passes**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_connection_laplacian_guard.m`
Expected: PASS — prints `PASSED: errors with nxrNotLoaded ...` then `ALL TESTS PASSED`. The `onCleanup` reloads nxr-compute so the session is left as found. `SKIP:` if nxr-compute is not installed at all.

- [ ] **Step 3: Commit**

```bash
git add dev/tests/test_connection_laplacian_guard.m
git commit -m "test(connection-laplacian): nxrNotLoaded guard when plugin unloaded"
```

---

## Notes for the implementer

- **Run tests via the MATLAB MCP** (`run_matlab_file` on the test path), not `runtests` — these are plain function-style scripts, matching the repo idiom (e.g. `dev/tests/test_tess_tangents.m`).
- **Live MATLAB session:** never use `clear` (it wipes `GlobalData` and hangs the session); edited `.m` files auto-reload. Unloading/reloading a plugin (Task 4) is fine.
- **Protocol prerequisite:** the cortex tests need the `TutorialAnatomy` protocol loaded with `cortex_20484V`. They `SKIP` cleanly otherwise — load that protocol and re-run if you see a SKIP.
- **`strfind` vs `contains`:** `strfind` is used for broad MATLAB-version compatibility, matching the existing test idiom.
```
