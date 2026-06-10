# Dirac Eigenbasis (Phase A — `tess_dirac_eigenmodes`) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Compute the per-hemisphere **Dirac eigenbasis** of the cortex and store it on the surface file as `TessMat.DiracEigen` (a `1×2` per-hemisphere struct array), for later expansion of the unconstrained leadfield (Phase B).

**Architecture:** New self-contained writer `tess_dirac_eigenmodes` (sibling to `tess_eigenmodes`): split with `tess_hemisplit` (atlas L/R, never `conncomp`); per hemisphere pull `L = operators(h,'dirac',τ)` and `Mg = operators(h,'mass','galerkin')` from nxr, form `B = kron(Mg,I₄)`, `eigs(L,B,K,'smallestabs')`, B-orthonormalize (Cholesky of `ΦᵀBΦ`), and store `Vectors`/`Values`/`Mass`/provenance per hemisphere via `bst_save('v7')`. MATLAB eigs (no new nxr command), matching `tess_eigenmodes`.

**Tech Stack:** MATLAB (Brainstorm toolbox), `nxr-compute` MEX (`operators … dirac τ`, `operators … mass galerkin`), `matlab.unittest`, run via the MATLAB MCP.

**Reference facts (grounded):**
- nxr eigenproblem: `L(τ)φ = λ·(Mg⊗I₄)φ`; `L = operators(h,'dirac',τ)` is `[4Vₕ×4Vₕ]` real symmetric (vertex-interleaved `4v+c`, order `[w,x,y,z]`); eigenvalues real, in **4-fold quaternionic multiplets**.
- `tess_eigenmodes` pattern: `[U,D]=eigs(L,M,nReq,'smallestabs',opts)` → sort → trim → `mOrthonormalize(U,M)` (Cholesky of `G=UᵀMU`, `U=U/R`); stores `Vectors/Values/nModes/Order/nComponents`.
- `tess_frame` pattern (per-hemi): `in_tess_bst` → atlas-label guard → `tess_hemisplit` → local submesh (`map(Fcs(fMask,:))`, `find(fMask)`) → `nxr_compute('create'/'destroy')` → `bst_save(file_fullpath,...,'v7')` + `bst_history`.
- Defaults: `Tau=0.5`, `K=400`. No DC-mode removal (the smoothest vector modes carry leadfield content — keep all K).

**MATLAB-session discipline:** before each test run, `rehash; clear test_tess_dirac_eigenmodes tess_dirac_eigenmodes;`. **Never** a bare `clear` (wipes Brainstorm `GlobalData`). Run via the MATLAB MCP. Tests that write to disk back up the `.mat` and restore via `onCleanup`. Tests use a small `K` (e.g. 40) for speed; the default stays 400.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/anatomy/tess_dirac_eigenmodes.m` | **Create.** Per-hemisphere Dirac eigenbasis → `TessMat.DiracEigen` (1×2). |
| `dev/tests/test_tess_dirac_eigenmodes.m` | **Create.** Property tests on the real 20484-vertex cortex. |

---

## Task 0: Prerequisite — ensure the installed MEX exposes `operators … dirac`

**Files:** none (verify + build/install only).

- [ ] **Step 1: Check whether the installed plugin binary has the dirac operator**

In MATLAB (MCP), load the **installed** plugin and probe:
```matlab
clear nxr_compute
% drop any repo-build copy so the installed plugin binary loads
warning('off','MATLAB:rmpath:DirNotFound');
rmpath('/Users/diellorbasha/workspace/research/code/nxr-compute/build/Release');
clear nxr_compute
addpath('/Users/diellorbasha/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b');
disp(which('nxr_compute'));
V=[1 0 0;-1 0 0;0 1 0;0 -1 0;0 0 1;0 0 -1]; F=[1 3 5;3 2 5;2 4 5;4 1 5;3 1 6;2 3 6;4 2 6;1 4 6];
h=nxr_compute('create',V,F); ok=true;
try, L=nxr_compute('operators',h,'dirac',0.5); catch ME, ok=false; disp(ME.message); end
nxr_compute('destroy',h);
fprintf('installed dirac present: %d\n', ok);
```
Only ever `clear nxr_compute` — never bare `clear`.

- [ ] **Step 2: If absent, build + install from nxr `main`**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
git status --short && git rev-parse --abbrev-ref HEAD   # expect main, clean-ish
bash scripts/build.sh Release
PLUGDIR="$HOME/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b"
cp "$PLUGDIR/nxr_compute.mexmaca64" "$PLUGDIR/nxr_compute.mexmaca64.bak.predirac-20260610"
cp build/Release/nxr_compute.mexmaca64 "$PLUGDIR/nxr_compute.mexmaca64"
ls -la "$PLUGDIR/"nxr_compute.mexmaca64*
```
Then re-run Step 1's probe and confirm `installed dirac present: 1`. (If the build sandbox blocks, retry the Bash with `dangerouslyDisableSandbox`.)

- [ ] **Step 3: No commit** (the installed binary is outside the repo). Proceed to Task 1.

---

## Task 1: `tess_dirac_eigenmodes` — compute + store the Dirac eigenbasis

**Files:**
- Create: `toolbox/anatomy/tess_dirac_eigenmodes.m`
- Test: `dev/tests/test_tess_dirac_eigenmodes.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_dirac_eigenmodes.m`:

```matlab
function tests = test_tess_dirac_eigenmodes
% Property tests for tess_dirac_eigenmodes (per-hemisphere Dirac eigenbasis).
tests = functiontests(localfunctions);
end

function SurfaceFile = local_cortex()
    if ~brainstorm('status'); brainstorm nogui; end
    SurfaceFile = local_find_cortex(20484);
end

function s = local_find_cortex(nVertTarget)
    s = '';
    P = bst_get('ProtocolSubjects');
    subj = P.Subject;
    if isfield(P,'DefaultSubject') && ~isempty(P.DefaultSubject)
        subj = [P.DefaultSubject, subj];
    end
    for k = 1:numel(subj)
        surfs = subj(k).Surface;
        for i = 1:numel(surfs)
            if strcmpi(surfs(i).SurfaceType, 'Cortex')
                m = load(file_fullpath(surfs(i).FileName), 'Vertices');
                if size(m.Vertices,1) == nVertTarget
                    s = surfs(i).FileName; return;
                end
            end
        end
    end
    error('No %d-vertex cortex found in the loaded protocol.', nVertTarget);
end

function test_store_1x2_shapes(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    DE = tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    verifyEqual(tc, size(DE), [1 2]);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'DiracEigen') && isequal(size(T.DiracEigen),[1 2]));
    for hh = 1:2
        nVh = numel(T.DiracEigen(hh).GlobalVertices);
        verifyEqual(tc, size(T.DiracEigen(hh).Vectors), [4*nVh, K]);
        verifyEqual(tc, numel(T.DiracEigen(hh).Values), K);
        verifyEqual(tc, T.DiracEigen(hh).Tau, 0.5);
        verifyEqual(tc, T.DiracEigen(hh).Provenance.Backend, 'nxr');
    end
    verifyEqual(tc, T.DiracEigen(1).Hemisphere, 'L');
    verifyEqual(tc, T.DiracEigen(2).Hemisphere, 'R');
end

function test_B_orthonormal(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        Phi = T.DiracEigen(hh).Vectors;
        B   = kron(T.DiracEigen(hh).Mass, speye(4));
        G   = Phi' * (B * Phi);
        verifyLessThan(tc, full(max(abs(G - speye(K)), [], 'all')), 1e-7);
    end
end

function test_values_ascending_nonneg_multiplets(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        lam = T.DiracEigen(hh).Values;
        verifyGreaterThanOrEqual(tc, min(lam), -1e-9);
        verifyTrue(tc, issorted(lam));
        % 4-fold quaternionic multiplets: first cluster of 4 is tight
        verifyLessThan(tc, lam(4) - lam(1), 1e-3 * (1 + abs(lam(4))));
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_tess_dirac_eigenmodes tess_dirac_eigenmodes; disp('rehashed');
runtests('dev/tests/test_tess_dirac_eigenmodes.m')
```
Expected: ERROR `Undefined function 'tess_dirac_eigenmodes'`.

- [ ] **Step 3: Write `toolbox/anatomy/tess_dirac_eigenmodes.m`**

```matlab
function DiracEigen = tess_dirac_eigenmodes(SurfaceFile, varargin)
% TESS_DIRAC_EIGENMODES: Per-hemisphere Dirac eigenbasis stored on the surface.
%
% USAGE:  DiracEigen = tess_dirac_eigenmodes(SurfaceFile)
%         DiracEigen = tess_dirac_eigenmodes(SurfaceFile, 'Tau',0.5, 'K',400, ...
%                                            'NoSave',1, 'ForceRecompute',1)
%
% DESCRIPTION:
%     Splits the cortex by hemisphere (tess_hemisplit, atlas L/R; never conncomp)
%     and, per hemisphere, solves the generalized Dirac eigenproblem
%         L(Tau) * phi = lambda * B * phi,   B = kron(Mg, I4)
%     where L(Tau)=operators(h,'dirac',Tau) is the [4Vh x 4Vh] relative-Dirac
%     family and Mg=operators(h,'mass','galerkin') is the vertex Galerkin mass.
%     The K smallest-magnitude eigenpairs are B-orthonormalized and stored as
%     TessMat.DiracEigen, a 1x2 per-hemisphere struct array ((1)=L,(2)=R), each:
%       .Vectors        [4Vh x K]  B-orthonormal (vertex-interleaved 4v+c, [w,x,y,z])
%       .Values         [K x 1]    ascending, >= 0
%       .Mass           [Vh x Vh]  galerkin vertex mass (so B = kron(Mass,I4))
%       .nModes, .Order, .Tau, .GlobalVertices, .Hemisphere, .Provenance
%
%     Eigenvalues come in 4-fold quaternionic multiplets; only the B-orthonormal
%     subspace matters (no multiplet canonicalization).
%
% OPTIONS:  'Tau' (0.5) | 'K' (400) | 'NoSave' (false) | 'ForceRecompute' (false)
%
% Requires the nxr-compute plugin (operators 'dirac'/'mass').
%
% SEE ALSO: tess_eigenmodes, tess_frame, tess_hemisplit

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

    % --- options ---
    Tau=0.5; K=400; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'tau',            Tau=varargin{i+1};
            case 'k',              K=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end
    if ~(isscalar(Tau) && Tau>=0 && Tau<=1)
        error('tess_dirac_eigenmodes:badTau', 'Tau must be a scalar in [0,1].');
    end

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return: matching Tau and K already stored ---
    if ~ForceRecompute && isfield(TessMat,'DiracEigen') && ~isempty(TessMat.DiracEigen) ...
       && isequal(size(TessMat.DiracEigen),[1 2]) ...
       && isequal(TessMat.DiracEigen(1).Tau, Tau) ...
       && (TessMat.DiracEigen(1).nModes == K)
        DiracEigen = TessMat.DiracEigen;
        return;
    end

    % --- require nxr-compute ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_dirac_eigenmodes:nxrUnavailable', 'requires nxr-compute: %s', errMsg);
    end

    % --- require Structures atlas with L/R labels (atlas split, never conncomp) ---
    hasLabels = false;
    if isfield(TessMat,'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts = TessMat.Atlas(iStruct).Scouts;
            labels = {scouts.Label}; regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasLabels = (any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'))) && ...
                        (any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R')));
        end
    end
    if ~hasLabels
        error('tess_dirac_eigenmodes:noHemisphereLabels', ...
            'Surface has no Structures atlas with left/right hemisphere labels.');
    end

    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_dirac_eigenmodes:connectedHemispheres', ...
            'Hemispheres are connected; each is eigensolved as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','dirac-eigenmodes', 'NxrVersion',nxrVer, ...
                  'Tau',Tau, 'K',K, 'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    clear DE
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_dirac_eigenmodes:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));
        nVh  = numel(vH);

        h  = nxr_compute('create', Vloc, Floc);
        L  = nxr_compute('operators', h, 'dirac', Tau);     % [4nVh x 4nVh]
        Mg = nxr_compute('operators', h, 'mass', 'galerkin'); % [nVh x nVh]
        nxr_compute('destroy', h);

        B = kron(Mg, speye(4));                              % [4nVh x 4nVh]
        if K > 4*nVh - 2
            error('tess_dirac_eigenmodes:tooManyModes', ...
                'K=%d exceeds 4*nV-2=%d on hemisphere %s.', K, 4*nVh-2, tags{hh});
        end
        nRequest = min(K + 8, 4*nVh - 2);
        opts = struct('tol', 1e-6, 'maxit', 1000, 'disp', 0);
        [U, D] = eigs(L, B, nRequest, 'smallestabs', opts);
        lam = real(diag(D));
        [lam, idx] = sort(lam, 'ascend');
        U = real(U(:, idx));
        % trim to K, B-orthonormalize, clamp tiny negatives
        U = U(:, 1:K); lam = lam(1:K);
        U = local_b_orthonormalize(U, B);
        lam(lam < 0) = 0;

        s = struct();
        s.Vectors        = U;
        s.Values         = lam;
        s.Mass           = Mg;
        s.nModes         = K;
        s.Order          = (1:K)';
        s.Tau            = Tau;
        s.GlobalVertices = vH;
        s.Hemisphere     = tags{hh};
        s.Provenance     = prov;
        DE(hh) = s;   %#ok<AGROW>
    end

    DiracEigen = DE;

    if ~NoSave
        TessMat_full = load(TessFile);
        TessMat_full.DiracEigen = DiracEigen;
        TessMat_full = bst_history('add', TessMat_full, 'dirac-eigenmodes', ...
            sprintf('Stored Dirac eigenbasis (tau=%.3g, K=%d) as per-hemisphere DiracEigen.', Tau, K));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end

% ----------------------------------------------------------------------------
function U = local_b_orthonormalize(U, B)
% Enforce U' * B * U = I (Cholesky of the Gram matrix; eig fallback).
    G = U' * (B * U); G = (G + G')/2;
    [R, flag] = chol(G);
    if flag == 0
        U = U / R;
    else
        [V, D] = eig(full(G)); d = diag(D);
        d = max(d, max(d) * 1e-12);
        U = U * (V * diag(1 ./ sqrt(d)) * V');
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_tess_dirac_eigenmodes tess_dirac_eigenmodes; disp('rehashed');
runtests('dev/tests/test_tess_dirac_eigenmodes.m')
```
Expected: 3 tests PASS. (eigs on `4Vₕ≈40k` for `K=40` runs in seconds–tens of seconds per hemisphere.)

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_dirac_eigenmodes.m dev/tests/test_tess_dirac_eigenmodes.m
git commit -m "feat(tess-dirac-eigenmodes): per-hemisphere Dirac eigenbasis stored on surface"
```

---

## Task 2: Cache-return + NoSave + Tau/K-mismatch recompute

**Files:**
- Modify: `dev/tests/test_tess_dirac_eigenmodes.m` (add tests; implementation already supports these)

- [ ] **Step 1: Add the tests**

Append to `dev/tests/test_tess_dirac_eigenmodes.m`:

```matlab
function test_cache_return_no_recompute(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    TF = load(TessFile);
    TF.DiracEigen(1).Provenance.ComputeDate = 'SENTINEL';
    bst_save(TessFile, TF, 'v7');

    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5);   % must cache-return
    T = in_tess_bst(SurfaceFile, 0);
    verifyEqual(tc, T.DiracEigen(1).Provenance.ComputeDate, 'SENTINEL');
end

function test_tau_mismatch_recomputes(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    K = 40;
    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.5, 'ForceRecompute', 1);
    TF = load(TessFile);
    TF.DiracEigen(1).Provenance.ComputeDate = 'SENTINEL';
    bst_save(TessFile, TF, 'v7');

    tess_dirac_eigenmodes(SurfaceFile, 'K', K, 'Tau', 0.75);  % different Tau -> recompute
    T = in_tess_bst(SurfaceFile, 0);
    verifyNotEqual(tc, T.DiracEigen(1).Provenance.ComputeDate, 'SENTINEL');
    verifyEqual(tc, T.DiracEigen(1).Tau, 0.75);
end

function test_nosave_returns_without_writing(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    TF = load(TessFile);
    if isfield(TF,'DiracEigen'), TF = rmfield(TF,'DiracEigen'); end
    bst_save(TessFile, TF, 'v7');

    DE = tess_dirac_eigenmodes(SurfaceFile, 'K', 40, 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(DE), [1 2]);
    T = in_tess_bst(SurfaceFile, 0);
    verifyFalse(tc, isfield(T,'DiracEigen'));
end
```

- [ ] **Step 2: Run the full test file**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_tess_dirac_eigenmodes tess_dirac_eigenmodes; disp('rehashed');
runtests('dev/tests/test_tess_dirac_eigenmodes.m')
```
Expected: 6 tests PASS.

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_tess_dirac_eigenmodes.m
git commit -m "test(tess-dirac-eigenmodes): cache-return, tau-mismatch recompute, NoSave"
```

---

## Final verification

- [ ] **Run the test file; expect 6 passing.**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_tess_dirac_eigenmodes tess_dirac_eigenmodes; disp('rehashed');
r = runtests('dev/tests/test_tess_dirac_eigenmodes.m');
fprintf('\nTOTAL: %d passed, %d failed\n', sum([r.Passed]), sum([r.Failed]));
assert(all([r.Passed]) && ~any([r.Failed]), 'Some tests failed.');
```

- [ ] **Then complete via superpowers:finishing-a-development-branch**, and proceed to Plan 2 (Phase B — `bst_dirac_eigenmode_leadfield`).

---

## Notes for the implementer

- **Atlas split, never `conncomp`:** unlike `tess_eigenmodes` (which uses `conncomp`), this writer must use `tess_hemisplit` (project hard rule). The guard errors clearly if the Structures L/R atlas is missing.
- **`Mass` is stored** (the small `[Vₕ×Vₕ]` galerkin mass), not `B` — Phase B rebuilds `B=kron(Mass,I₄)` with no nxr call.
- **No DC removal:** keep all K smallest modes (the smoothest vector modes carry leadfield content, unlike scalar source patterns).
- **eigs cost:** `4Vₕ≈40k`; `'smallestabs'` does shift-invert (a sparse factorization). For the default `K=400` this is the heavy step (≈minutes/hemisphere) — acceptable, one-time, cached on the surface. Tests use `K=40`.
- **4-fold multiplets:** do not canonicalize; `local_b_orthonormalize` makes the basis B-orthonormal, which is all Phase B needs.
