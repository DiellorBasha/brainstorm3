# Connection Eigenmodes (M2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `ConnEigenmodes`, a canonical intrinsic tess-file axis holding the complex eigenmodes of the connection Laplacian (vector-field basis), as a parallel quartet to the scalar `Eigenmodes` axis.

**Architecture:** `tess_conn_eigenmodes` assembles `[K, M]` via `tess_connection_laplacian(V,F,'nSym',1)` and solves per connected component with MATLAB `eigs` on the block `K(idx,idx)` (the operator is block-diagonal across hemispheres). `out_/in_tess_conn_eigenmodes` persist/load `TessMat.ConnEigenmodes` (complex `single` vectors, `double` sparse operators). `bst_conn_eigenmodes_ensure` reuses-or-computes, matching the scalar axis's per-component count. The scalar functions and `tess_connection_laplacian` are not modified.

**Tech Stack:** MATLAB, Brainstorm (`in_tess_bst`, `bst_save`, `bst_history`, `file_fullpath`, `tess_vertconn`, `tess_manifold`, `tess_eigenmodes`, `bst_eigenmodes_ensure`, `bst_plugin`), nxr-compute (via `tess_connection_laplacian`), MATLAB `eigs`.

**Reference spec:** `dev/connection_eigenmodes_integration.md`

---

## File Structure

| Path | Responsibility |
|---|---|
| `toolbox/anatomy/tess_conn_eigenmodes.m` | Compute the `ConnEigenmodes` struct (per-component complex eigs). |
| `toolbox/io/out_tess_conn_eigenmodes.m` | Write `TessMat.ConnEigenmodes` (complex single vectors, double operators). |
| `toolbox/io/in_tess_conn_eigenmodes.m` | Read it back, cast to double complex, backfill metadata. |
| `toolbox/anatomy/bst_conn_eigenmodes_ensure.m` | Reuse-or-compute; match scalar per-component count. |
| `dev/tests/test_conn_eigenmodes_compute.m` | Compute correctness (shape/type/eigenvalues/block/orthonormality). |
| `dev/tests/test_conn_eigenmodes_roundtrip.m` | `out`→`in` preserves complex vectors + operator. |
| `dev/tests/test_conn_eigenmodes_ensure.m` | Match-scalar derivation + idempotent reuse. |

No changes to `tess_eigenmodes`, `in_/out_tess_eigenmodes`, `bst_eigenmodes_ensure`, or `tess_connection_laplacian`.

**Test environment notes (apply to every test):** run via the MATLAB MCP (`run_matlab_file`) against the live session (Brainstorm running, TutorialAuditory loaded, Subject01 has a 20484-vertex cortex). Never run MATLAB `clear`. Each test resolves the cortex by vertex count (resolver duplicated per file, matching the repo idiom) and SKIPs cleanly if absent. Tests use a small `nModes` for speed.

---

## Task 1: `tess_conn_eigenmodes` + compute test

**Files:**
- Create: `toolbox/anatomy/tess_conn_eigenmodes.m`
- Test: `dev/tests/test_conn_eigenmodes_compute.m`

- [ ] **Step 1: Write the failing compute test**

Create `dev/tests/test_conn_eigenmodes_compute.m`:

```matlab
function test_conn_eigenmodes_compute
% Compute correctness for tess_conn_eigenmodes on a real 20484-vertex cortex:
% complex eigenvectors, real positive eigenvalues, per-component block structure,
% canonical Order, and M-orthonormality within a component.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
nV = size(V, 1);

k = 20;
[ConnEig, K, M] = tess_conn_eigenmodes(V, F, 'nModes', k);

% --- Type / shape ---
assert(~isreal(ConnEig.Vectors), 'Vectors must be complex.');
assert(size(ConnEig.Vectors, 1) == nV, 'Vectors must have nV rows.');
assert(ConnEig.nModes == size(ConnEig.Vectors, 2), 'nModes must equal column count.');
assert(isreal(ConnEig.Values) && all(ConnEig.Values > 0), 'Values must be real and > 0.');

% --- Metadata ---
assert(strcmp(ConnEig.OperatorType, 'Connection-LeviCivita'), 'OperatorType mismatch.');
assert(ConnEig.nSym == 1, 'nSym should be 1.');
assert(ConnEig.nRemoved == 0, 'No DC mode is removed for the connection bundle.');
assert(issparse(ConnEig.ConnLaplacian) && ~isreal(ConnEig.ConnLaplacian), 'ConnLaplacian must be complex sparse.');
assert(isequal(size(K), [nV nV]) && ~isreal(K), 'Returned K must be nV x nV complex.');

% --- Canonical order sorts Values ascending ---
assert(issorted(ConnEig.Values(ConnEig.Order)), 'Order must sort Values ascending.');

% --- Per-component checks (first component) ---
compId = conncomp(graph(tess_vertconn(V, F)));
assert(ConnEig.nComponents == max(compId), 'nComponents mismatch.');
c    = 1;
idx  = find(compId == c);
cols = find(ConnEig.Component == c);
Uc   = ConnEig.Vectors(idx, cols);
Mc   = M(idx, idx);
% Block structure: a component's modes vanish off that component.
other = setdiff((1:nV)', idx);
assert(max(max(abs(ConnEig.Vectors(other, cols)))) < 1e-10, 'Component modes must vanish off-component.');
% M-orthonormality (Hermitian inner product; '' is conjugate transpose).
G = Uc' * Mc * Uc;
assert(max(max(abs(G - eye(size(G))))) < 1e-6, 'Per-component modes must be M-orthonormal.');

fprintf('PASSED: %d-mode connection eigenbasis (nV=%d, %d components): complex, real>0 eigenvalues, block-structured, M-orthonormal.\n', ...
    ConnEig.nModes, nV, ConnEig.nComponents);
fprintf('ALL TESTS PASSED: test_conn_eigenmodes_compute\n');
end


function SurfaceFile = find_cortex_20484V()
% Return a Cortex surface with exactly 20484 vertices in the current protocol,
% preferring one with a FreeSurfer registration sphere; '' if none.
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
```

- [ ] **Step 2: Run the test to verify it FAILS**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_compute.m`
Expected: FAIL — `Unrecognized function or variable 'tess_conn_eigenmodes'`.

- [ ] **Step 3: Write the function**

Create `toolbox/anatomy/tess_conn_eigenmodes.m`:

```matlab
function [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces, varargin)
% TESS_CONN_EIGENMODES: Eigenmodes of the connection Laplacian (vector-field basis).
%
% USAGE:  [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces)
%         [ConnEig, K, M] = tess_conn_eigenmodes(Vertices, Faces, 'nModes', 300)
%
% DESCRIPTION:
%     Computes the smallest eigenmodes of the complex-Hermitian connection
%     Laplacian (n-RoSy bundle, nSym=1) of a triangle mesh, solving the
%     generalized problem K*phi = lambda*M*phi independently per connected
%     component (the operator is block-diagonal across components, e.g. the two
%     cortical hemispheres). The result is a vector-field spectral basis, the
%     tangent-field sibling of the scalar Laplace-Beltrami eigenmodes
%     (tess_eigenmodes).
%
%     The operator and mass come from tess_connection_laplacian (nxr-compute,
%     no MATLAB fallback). nxr's own eigensolver is real-only, so the complex
%     Hermitian eigenproblem is solved with MATLAB eigs. There is NO DC/zero
%     mode (no globally consistent parallel vector field on a curved closed
%     surface), so none is removed. Eigenvectors are stored raw, in nxr's
%     intrinsic per-vertex frames; phase readout/gauge alignment is a later step.
%
% INPUTS:
%     Vertices : [nV x 3] vertex positions.
%     Faces    : [nF x 3] triangle vertex indices (1-based).
%
% OPTIONS:
%     nModes         : Modes per connected component (default 300).
%     nSym           : n-RoSy symmetry (default 1; true vector field).
%     Regularization : Operator diagonal epsilon (default 1e-8).
%     Tolerance      : eigs convergence tolerance (default 1e-10).
%     Verbose        : Print progress (default 1).
%
% OUTPUTS:
%     ConnEig : struct with fields Vectors [nV x nModes] complex (M-orthonormal,
%               block-structured per component), Values [nModes x 1] real, nModes,
%               Component, CompRank, Order, nComponents, MassMatrix (lumped, sparse),
%               ConnLaplacian (K, complex sparse), OperatorType, nSym,
%               Regularization, Sigma, Tolerance, nRemoved (=0), ComputeTime.
%     K       : [nV x nV] complex Hermitian connection Laplacian (whole mesh).
%     M       : [nV x nV] real diagonal lumped vertex mass (whole mesh).
%
% SEE ALSO: tess_connection_laplacian, tess_eigenmodes, out_tess_conn_eigenmodes

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
nModes         = 300;
nSym           = 1;
Regularization = 1e-8;
Tolerance      = 1e-10;
Verbose        = 1;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'nmodes',         nModes = varargin{i+1};
        case 'nsym',           nSym = varargin{i+1};
        case 'regularization', Regularization = varargin{i+1};
        case 'tolerance',      Tolerance = varargin{i+1};
        case 'verbose',        Verbose = varargin{i+1};
    end
end
tStart = tic;

if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Vertices and Faces must have 3 columns.');
end
Faces = double(Faces);
nV    = size(Vertices, 1);

%% ===== WHOLE-MESH OPERATOR + MASS (complex Hermitian) =====
[K, M] = tess_connection_laplacian(Vertices, Faces, 'nSym', nSym, 'Regularization', Regularization);

%% ===== CONNECTED COMPONENTS =====
% Solve each disconnected component (e.g. a hemisphere) independently. The
% connection Laplacian has no cross-component entries, so K(idx,idx) is the exact
% component operator; restricting beats rebuilding sub-meshes.
compId = conncomp(graph(tess_vertconn(Vertices, Faces)));   % [1 x nV]
nComp  = max(compId);
if Verbose
    fprintf('BST> tess_conn_eigenmodes: %d connected component(s).\n', nComp);
end

%% ===== SOLVE PER COMPONENT =====
VectorsAll = complex(zeros(nV, 0));
ValuesAll  = zeros(0, 1);
Component  = zeros(0, 1);
CompRank   = zeros(0, 1);
eigsOpts   = struct('tol', Tolerance);
for c = 1:nComp
    vIdx = find(compId == c);
    nvc  = numel(vIdx);
    kC   = min(nModes, nvc - 2);            % leave an eigs margin
    Kc   = K(vIdx, vIdx);
    Mc   = M(vIdx, vIdx);
    [Uc, Dc] = eigs(Kc, Mc, kC, 'smallestabs', eigsOpts);
    lam = real(diag(Dc));
    [lam, ord] = sort(lam, 'ascend');
    Uc  = Uc(:, ord);
    % Guarantee M-orthonormal magnitude (eigs returns B-normalized; this guards it).
    nrm = sqrt(real(sum(conj(Uc) .* (Mc * Uc), 1)));
    Uc  = Uc ./ nrm;
    Ufull = complex(zeros(nV, kC));
    Ufull(vIdx, :) = Uc;
    VectorsAll = [VectorsAll, Ufull];                       %#ok<AGROW>
    ValuesAll  = [ValuesAll;  lam(:)];                      %#ok<AGROW>
    Component  = [Component;  c * ones(kC, 1)];             %#ok<AGROW>
    CompRank   = [CompRank;   (1:kC)'];                     %#ok<AGROW>
    if Verbose
        fprintf('BST> tess_conn_eigenmodes: component %d: %d modes, range [%.3g, %.3g].\n', ...
            c, kC, lam(1), lam(end));
    end
end

%% ===== PACKAGE =====
ConnEig = struct();
ConnEig.Vectors       = VectorsAll;
ConnEig.Values        = ValuesAll;
ConnEig.nModes        = numel(ValuesAll);
ConnEig.Component     = Component;
ConnEig.CompRank      = CompRank;
% Canonical global order: ascending eigenvalue across all components, so
% Vectors(:,Order(1:K)) are the whole-brain lowest-frequency modes.
[~, ConnEig.Order]    = sort(ValuesAll, 'ascend');
ConnEig.nComponents   = nComp;
ConnEig.MassMatrix    = M;                       % lumped vertex mass (basis is M-orthonormal)
ConnEig.ConnLaplacian = K;                       % complex Hermitian operator (reuse downstream)
ConnEig.OperatorType  = 'Connection-LeviCivita';
ConnEig.nSym          = nSym;
ConnEig.Regularization = Regularization;
ConnEig.Sigma         = 'smallestabs';
ConnEig.Tolerance     = Tolerance;
ConnEig.nRemoved      = 0;                        % no DC mode in the connection bundle
ConnEig.ComputeTime   = toc(tStart);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_compute.m`
Expected: PASS — prints the per-component ranges and `ALL TESTS PASSED: test_conn_eigenmodes_compute`. A `SKIP:` is NOT acceptable (a 20484V cortex is loaded). If `eigs` errors on the `'smallestabs', eigsOpts` form, retry with `eigs(Kc, Mc, kC, 'smallestabs')` (no opts) and set `ConnEig.Tolerance = 1e-14`; report the change. Do not otherwise alter the operator.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_conn_eigenmodes.m dev/tests/test_conn_eigenmodes_compute.m
git commit -m "feat(conn-eigenmodes): tess_conn_eigenmodes (per-component complex eigs) + compute test"
```

---

## Task 2: `out_`/`in_tess_conn_eigenmodes` + round-trip test

**Files:**
- Create: `toolbox/io/out_tess_conn_eigenmodes.m`
- Create: `toolbox/io/in_tess_conn_eigenmodes.m`
- Test: `dev/tests/test_conn_eigenmodes_roundtrip.m`

- [ ] **Step 1: Write the failing round-trip test**

Create `dev/tests/test_conn_eigenmodes_roundtrip.m`:

```matlab
function test_conn_eigenmodes_roundtrip
% out_tess_conn_eigenmodes then in_tess_conn_eigenmodes preserves the complex
% eigenvectors (single round-trip) and the complex operator. Works on a temp COPY
% so the DB surface is not mutated.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

srcFile = find_cortex_20484V();
if isempty(srcFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
% Temp copy named tess_cortex_*.mat so file_fullpath accepts the absolute path.
tmpFile = fullfile(tempdir, 'tess_cortex_conneig_roundtrip.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);

% Not computed yet.
[~, isComputed0] = in_tess_conn_eigenmodes(tmpFile);
assert(~isComputed0, 'ConnEigenmodes should be absent before computing.');

% Compute + store.
ConnEig = tess_conn_eigenmodes(V, F, 'nModes', 20);
out_tess_conn_eigenmodes(tmpFile, ConnEig, V, F);

% Load back.
[ConnEig2, isComputed] = in_tess_conn_eigenmodes(tmpFile);
assert(isComputed, 'ConnEigenmodes should be present after store.');
assert(~isreal(ConnEig2.Vectors), 'Loaded Vectors must be complex.');
assert(isa(ConnEig2.Vectors, 'double'), 'Loaded Vectors must be cast to double.');
assert(isequal(size(ConnEig2.Vectors), size(ConnEig.Vectors)), 'Vector shape must be preserved.');
% Stored as single, loaded as double(single): compare against the single-cast original.
d = max(max(abs(double(single(ConnEig.Vectors)) - ConnEig2.Vectors)));
assert(d < 1e-5, 'Vectors must survive the single round-trip (max diff = %g).', d);
assert(isfield(ConnEig2, 'ConnLaplacian') && issparse(ConnEig2.ConnLaplacian) ...
       && ~isreal(ConnEig2.ConnLaplacian), 'ConnLaplacian must round-trip as complex sparse.');
assert(strcmp(ConnEig2.OperatorType, 'Connection-LeviCivita'), 'OperatorType must round-trip.');
assert(ConnEig2.nModes == ConnEig.nModes, 'nModes must round-trip.');

fprintf('PASSED: ConnEigenmodes round-trip (%d modes; vector max diff %g; complex operator preserved).\n', ...
    ConnEig2.nModes, d);
fprintf('ALL TESTS PASSED: test_conn_eigenmodes_roundtrip\n');
end


function SurfaceFile = find_cortex_20484V()
% Return a Cortex surface with exactly 20484 vertices in the current protocol,
% preferring one with a FreeSurfer registration sphere; '' if none.
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
```

- [ ] **Step 2: Run the test to verify it FAILS**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_roundtrip.m`
Expected: FAIL — `Unrecognized function or variable 'in_tess_conn_eigenmodes'` (or `out_...`).

- [ ] **Step 3a: Write `out_tess_conn_eigenmodes`**

Create `toolbox/io/out_tess_conn_eigenmodes.m`:

```matlab
function TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces, isInteractive)
% OUT_TESS_CONN_EIGENMODES: Save connection-Laplacian eigenmodes into a surface file.
%
% USAGE:  TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces)
%         TessMat = out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Vertices, Faces, isInteractive)
%
% SEE ALSO: in_tess_conn_eigenmodes, tess_conn_eigenmodes

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

if nargin < 5
    isInteractive = true;
end

%% ===== VALIDATE INPUTS =====
if ~isstruct(ConnEig) || ~isfield(ConnEig, 'Vectors') || ~isfield(ConnEig, 'Values')
    error('ConnEig must be a structure with Vectors and Values fields.');
end
if size(ConnEig.Vectors, 1) ~= size(Vertices, 1)
    error('ConnEig.Vectors (%d rows) must match Vertices (%d rows).', ...
        size(ConnEig.Vectors, 1), size(Vertices, 1));
end

%% ===== LOAD SURFACE FILE =====
SurfaceFileFull = file_fullpath(SurfaceFile);
if ~file_exist(SurfaceFileFull)
    error('Surface file not found: %s', SurfaceFileFull);
end
TessMat = load(SurfaceFileFull);

%% ===== STORE (complex single vectors; double sparse operators) =====
Store = struct();
Store.Vectors     = single(ConnEig.Vectors);     % complex single
Store.Values      = ConnEig.Values(:);
Store.nModes      = ConnEig.nModes;
if isfield(ConnEig, 'Component'),   Store.Component   = ConnEig.Component(:);   end
if isfield(ConnEig, 'CompRank'),    Store.CompRank    = ConnEig.CompRank(:);    end
if isfield(ConnEig, 'Order'),       Store.Order       = ConnEig.Order(:);       end
if isfield(ConnEig, 'nComponents'), Store.nComponents = ConnEig.nComponents;    end
% Sparse stays double (MATLAB has no single sparse); negligible vs Vectors.
if isfield(ConnEig, 'MassMatrix')    && ~isempty(ConnEig.MassMatrix),    Store.MassMatrix    = ConnEig.MassMatrix;    end
if isfield(ConnEig, 'ConnLaplacian') && ~isempty(ConnEig.ConnLaplacian), Store.ConnLaplacian = ConnEig.ConnLaplacian; end
Store.OperatorType   = ConnEig.OperatorType;
Store.nSym           = ConnEig.nSym;
Store.Regularization = ConnEig.Regularization;
Store.Sigma          = ConnEig.Sigma;
Store.Tolerance      = ConnEig.Tolerance;
Store.nRemoved       = ConnEig.nRemoved;
Store.ComputeTime    = ConnEig.ComputeTime;
Store.ComputeDate    = datestr(now, 'yyyy-mm-dd HH:MM:SS');
TessMat.ConnEigenmodes = Store;

%% ===== HISTORY + SINGLE SAVE =====
TessMat = bst_history('add', TessMat, 'conn_eigenmodes', ...
    sprintf('Computed %d connection-Laplacian eigenmodes (nSym=%d, reg=%.1e)', ...
        ConnEig.nModes, ConnEig.nSym, ConnEig.Regularization));
bst_save(SurfaceFileFull, TessMat, 'v7');

if isInteractive
    fprintf('BST> Saved %d connection eigenmodes to: %s\n', ConnEig.nModes, SurfaceFile);
end
end
```

- [ ] **Step 3b: Write `in_tess_conn_eigenmodes`**

Create `toolbox/io/in_tess_conn_eigenmodes.m`:

```matlab
function [ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile)
% IN_TESS_CONN_EIGENMODES: Load connection-Laplacian eigenmodes from a surface file.
%
% USAGE:  [ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile)
%
% DESCRIPTION:
%     Loads the embedded ConnEigenmodes field via in_tess_bst. Returns [] and
%     false if connection eigenmodes have not been computed for this surface.
%
% SEE ALSO: out_tess_conn_eigenmodes, tess_conn_eigenmodes, in_tess_bst

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

ConnEig    = [];
isComputed = false;

% isComputeMissing=0: a frequent read must not trigger curvature/normal recompute.
TessMat = in_tess_bst(SurfaceFile, 0);

if ~isfield(TessMat, 'ConnEigenmodes') || isempty(TessMat.ConnEigenmodes)
    return;
end

ConnEig = TessMat.ConnEigenmodes;
if isfield(ConnEig, 'Vectors') && isa(ConnEig.Vectors, 'single')
    ConnEig.Vectors = double(ConnEig.Vectors);   % preserves complex
end
nK = size(ConnEig.Vectors, 2);
if ~isfield(ConnEig, 'Component') || isempty(ConnEig.Component)
    ConnEig.Component = ones(nK, 1);
end
if ~isfield(ConnEig, 'CompRank') || isempty(ConnEig.CompRank)
    ConnEig.CompRank = (1:nK)';
end
if ~isfield(ConnEig, 'Order') || isempty(ConnEig.Order) || numel(ConnEig.Order) ~= nK
    [~, ConnEig.Order] = sort(double(ConnEig.Values(:)), 'ascend');
end
if ~isfield(ConnEig, 'nComponents') || isempty(ConnEig.nComponents)
    ConnEig.nComponents = max(ConnEig.Component(:));
end
isComputed = true;
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_roundtrip.m`
Expected: PASS — prints the round-trip diff and `ALL TESTS PASSED: test_conn_eigenmodes_roundtrip`. A `SKIP:` is NOT acceptable.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/io/out_tess_conn_eigenmodes.m toolbox/io/in_tess_conn_eigenmodes.m dev/tests/test_conn_eigenmodes_roundtrip.m
git commit -m "feat(conn-eigenmodes): in_/out_tess_conn_eigenmodes persistence + round-trip test"
```

---

## Task 3: `bst_conn_eigenmodes_ensure` + ensure test

**Files:**
- Create: `toolbox/anatomy/bst_conn_eigenmodes_ensure.m`
- Test: `dev/tests/test_conn_eigenmodes_ensure.m`

- [ ] **Step 1: Write the failing ensure test**

Create `dev/tests/test_conn_eigenmodes_ensure.m`:

```matlab
function test_conn_eigenmodes_ensure
% bst_conn_eigenmodes_ensure: (1) with no count, derives the per-component count
% from the surface's scalar Eigenmodes axis (match-scalar); (2) reuses an existing
% ConnEigenmodes axis idempotently. Works on a temp COPY; uses small mode counts.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute plugin required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

srcFile = find_cortex_20484V();
if isempty(srcFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol (e.g. load TutorialAuditory / Subject01).\n');
    return;
end
tmpFile = fullfile(tempdir, 'tess_cortex_conneig_ensure.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V = TessMat.Vertices;
F = double(TessMat.Faces);

% --- Seed a SMALL scalar Eigenmodes axis so the match-scalar derivation is fast ---
sEig = tess_eigenmodes(V, F, 'nModes', 15, 'MassType', 'barycentric', 'RemoveDC', 1);
out_tess_eigenmodes(tmpFile, sEig, V, F);
sEigStored = in_tess_eigenmodes(tmpFile);
expectedPerHemi = max(1, round(sEigStored.nModes / sEigStored.nComponents));

% --- Ensure with no count: derives + computes + stores ---
ConnEig = bst_conn_eigenmodes_ensure(tmpFile);
assert(ConnEig.nComponents == sEigStored.nComponents, 'Component count should match the mesh.');
for c = 1:ConnEig.nComponents
    nc = sum(ConnEig.Component == c);
    assert(nc == expectedPerHemi, ...
        'Component %d: expected %d connection modes (match scalar), got %d.', c, expectedPerHemi, nc);
end
fprintf('PASSED: match-scalar count = %d modes/component.\n', expectedPerHemi);

% --- Second call reuses (idempotent, fast: load not recompute) ---
tReuse = tic;
ConnEig2 = bst_conn_eigenmodes_ensure(tmpFile);
elapsed = toc(tReuse);
assert(ConnEig2.nModes == ConnEig.nModes, 'Reuse must return the same nModes.');
assert(elapsed < 5, 'Reuse should be fast (no recompute); took %.1fs.', elapsed);
fprintf('PASSED: reuse is idempotent (%.2fs, %d modes).\n', elapsed, ConnEig2.nModes);

fprintf('ALL TESTS PASSED: test_conn_eigenmodes_ensure\n');
end


function SurfaceFile = find_cortex_20484V()
% Return a Cortex surface with exactly 20484 vertices in the current protocol,
% preferring one with a FreeSurfer registration sphere; '' if none.
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects)
    return;
end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex')
            continue;
        end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Vertices', 'Reg');
        catch
            continue;
        end
        if size(T.Vertices, 1) ~= 20484
            continue;
        end
        hasReg = isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
                 && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices);
        if hasReg
            SurfaceFile = surf(iF).FileName;
            return;
        elseif isempty(fallback)
            fallback = surf(iF).FileName;
        end
    end
end
if isempty(SurfaceFile)
    SurfaceFile = fallback;
end
end
```

- [ ] **Step 2: Run the test to verify it FAILS**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_ensure.m`
Expected: FAIL — `Unrecognized function or variable 'bst_conn_eigenmodes_ensure'`.

- [ ] **Step 3: Write `bst_conn_eigenmodes_ensure`**

Create `toolbox/anatomy/bst_conn_eigenmodes_ensure.m`:

```matlab
function ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
% BST_CONN_EIGENMODES_ENSURE: Return a surface's canonical connection eigenmodes,
% computing a default set if absent. Sibling of bst_eigenmodes_ensure for the
% vector-field (connection-Laplacian) axis.
%
% USAGE:  ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile)
%         ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile, nModesPerHemi)
%
% DESCRIPTION:
%     If the surface already carries ConnEigenmodes, they are returned as-is.
%     Otherwise a default set is computed and stored. When nModesPerHemi is not
%     given, the per-component count is matched to the surface's scalar Eigenmodes
%     axis (round(nModes / nComponents)), ensuring that axis exists first via
%     bst_eigenmodes_ensure. NO repair is attempted: a non-manifold surface raises
%     an error (repair changes the vertex count and breaks surface<->eigenmode
%     consistency). Remesh to an icosphere instead.
%
% SEE ALSO: tess_conn_eigenmodes, in_tess_conn_eigenmodes, out_tess_conn_eigenmodes,
%           bst_eigenmodes_ensure, tess_manifold

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

% Reuse existing connection eigenmodes when present.
[ConnEig, isComputed] = in_tess_conn_eigenmodes(SurfaceFile);
if isComputed && ~isempty(ConnEig)
    return;
end

% Derive the per-component count from the scalar axis if not given (ensuring the
% scalar axis exists first; bst_eigenmodes_ensure reuses it when already present).
if nargin < 2 || isempty(nModesPerHemi)
    sEig = bst_eigenmodes_ensure(SurfaceFile);
    nModesPerHemi = max(1, round(sEig.nModes / max(1, sEig.nComponents)));
end

% Compute a default set. Guard manifoldness first -- no silent repair.
Surf = in_tess_bst(SurfaceFile, 0);
mani = tess_manifold(Surf.Vertices, Surf.Faces);
if isstruct(mani) && isfield(mani, 'isManifold') && ~mani.isManifold
    error('bst_conn_eigenmodes_ensure:NonManifold', ...
        ['Surface %s is non-manifold; connection eigenmodes require a 2-manifold mesh. ' ...
         'Remesh to an icosphere (or repair manually) and retry.'], SurfaceFile);
end

ConnEig = tess_conn_eigenmodes(Surf.Vertices, Surf.Faces, 'nModes', nModesPerHemi);
out_tess_conn_eigenmodes(SurfaceFile, ConnEig, Surf.Vertices, Surf.Faces);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_conn_eigenmodes_ensure.m`
Expected: PASS — prints the match-scalar count and reuse timing, then `ALL TESTS PASSED: test_conn_eigenmodes_ensure`. A `SKIP:` is NOT acceptable.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/bst_conn_eigenmodes_ensure.m dev/tests/test_conn_eigenmodes_ensure.m
git commit -m "feat(conn-eigenmodes): bst_conn_eigenmodes_ensure (match scalar count) + ensure test"
```

---

## Notes for the implementer

- **Run tests via the MATLAB MCP** (`run_matlab_file` on the absolute test path), not `runtests`.
- **Never `clear`** in the live MATLAB session (wipes GlobalData). Edited `.m` files auto-reload.
- **Protocol prerequisite:** a 20484-vertex cortex must be loaded (on this machine: TutorialAuditory / Subject01). Tests SKIP cleanly otherwise.
- **Do not modify** `tess_eigenmodes`, `in_/out_tess_eigenmodes`, `bst_eigenmodes_ensure`, or `tess_connection_laplacian`.
- **`eigs` complex Hermitian:** the operator is Hermitian (`tess_connection_laplacian` symmetrizes) and `M` is real-diagonal positive, so `eigs(Kc, Mc, k, 'smallestabs')` returns real eigenvalues. If the `eigsOpts` form errors, use the no-opts form (Task 1 Step 4 note).
```
