# Phase 2 — tess_eigen Implementation Plan (Milestone 1 gate)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `tess_eigen` to the clean branch — per-hemisphere shift-invert LBO eigensolve with the eigenbasis embedded in the surface file (`SurfaceMat.Eigen.<Variant>`, spec §3.3–3.5) — completing Milestone 1 (Gate 2).

**Architecture:** One new function on the clean branch consuming Phase 1's
`tess_operators` pencil. Solver method is USER-DECIDED: shift-invert
(`eigs(A + tau*B, B, K, 'smallestabs')`, eigenvalues recovered as mu - tau);
Task 1 pins the shift scale and reproducibility settings by experiment.
Storage follows the approved data model exactly: named-variant slots, machine
readable `.Operator` recipe + `.Solver` record, mrimask append-save pattern,
`bst_history` rows, no DB calls. Dirac remains OUT (Open Question 2).

**Tech Stack:** MATLAB R2023b (`matlab -batch`), git worktree `~/workspace/research/code/brainstorm3-clean`, harness `dev/verify/phase2/`.

**Spec:** `docs/superpowers/specs/2026-08-07-cortical-flow-distillation-design.md` (§3.3, §3.4, §3.5, §4 Phase 2, §6 OQ1)

## Global Constraints

- Clean-branch (worktree) commits: upstream style subjects, NO Claude
  trailers, explicit path staging only, canonical Brainstorm headers
  (copy boilerplate from `toolbox/anatomy/tess_hemisplit.m`, then
  `% Authors: Diellor Basha, 2026`).
- Worktree `~/workspace/research/code/brainstorm3-clean` on
  `feature/cortical-flow-core` (7 commits ahead of master; HEAD b08b2a42).
  Never switch branches in the main checkout.
- Harness: `dev/verify/phase2/` in the main checkout (`git add -f` at commit
  time, Task 5 only). Runner: `dev/verify/phase0/run_matlab.sh` — Task 1
  extends its addpath to phase2 (same one-line pattern as phase1).
- Real-cortex surface (clean-branch-readable, schema 5.03):
  `/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat`
  Tests that WRITE the Eigen field operate on a COPY of this file in the
  harness dir — never mutate the Gate-0 original (the Phase 1 oracle's
  `meta.SurfaceFile` points at it).
- The NEVER-conncomp rule and the {L,R}/sorted-vH conventions from Phase 1
  bind here; `tess_eigen` reaches geometry ONLY through
  `tess_operators(Surface, recipe)`.
- Solver method is fixed (user decision): shift-invert. Plain
  `'smallestabs'` at sigma=0 appears in the Task 1 experiment as a
  comparison arm only, never in shipped code.
- MATLAB runs: Bash timeout 600000 ms. Scratch under
  `/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad/`.
- Success metrics used throughout (K modes, per hemisphere):
  relative pencil residual `max_k ||A*phi_k - lambda_k*B*phi_k|| / ||(A + B)*phi_k||`
  <= 1e-10; B-orthonormality `norm(Phi'*B*Phi - eye(K), 'fro')` <= 1e-10;
  lambda_1 in [-1e-8, 1e-8]*sigmaScale and ascending order.

---

### Task 1: Shift experiment + phase2 harness scaffold

**Files:**
- Create (main checkout): `dev/verify/phase2/experiment_shift_invert.m`
- Create (by running it): `dev/verify/phase2/shift_experiment_results.md`
- Modify (main checkout): `dev/verify/phase0/run_matlab.sh` (addpath phase2)

**Interfaces:**
- Produces: the pinned solver parameters consumed verbatim by Task 2's
  implementation: `TauRel` (relative shift; tau = TauRel * sigmaScale where
  sigmaScale = trace(A)/trace(B)), StartVector convention, and the measured
  evidence justifying them. Written at the top of shift_experiment_results.md
  as `PINNED: TauRel=<value>, StartVector=deterministic-ones`.

- [ ] **Step 1: Extend the runner** — in `dev/verify/phase0/run_matlab.sh`
  change `addpath('$HARNESS'); addpath('$HARNESS/../phase1')` to also append
  `addpath('$HARNESS/../phase2')`.

- [ ] **Step 2: Write the experiment script**

Create `dev/verify/phase2/experiment_shift_invert.m`. It uses the Phase 1
oracle pencil (real cortex, both hemispheres) — no new geometry code:

```matlab
% EXPERIMENT_SHIFT_INVERT: pin the LBO shift-invert parameters (Open Question 1).
% Arms: sigma=0 'smallestabs' (comparison only) vs shift-invert at
% TauRel in {1e-6, 1e-4, 1e-2}. K=400. Metrics: residual, B-orthonormality,
% lambda_1 zero-mode recovery, run-to-run reproducibility (2 runs / arm,
% principal-angle subspace correlation), cross-arm subspace agreement, runtime.
S = load(fullfile(fileparts(mfilename('fullpath')), '..', 'phase1', 'oracle_lbo_sub0002.mat'));
K = 400;
TauRelList = [1e-6, 1e-4, 1e-2];
out = fopen(fullfile(fileparts(mfilename('fullpath')), 'shift_experiment_results.md'), 'w');
fprintf(out, '# Shift-invert experiment (K=%d, real cortex pencil)\n\n', K);
for hh = 1:2
    A = S.A{hh}; B = S.B{hh}; n = size(A,1);
    sigmaScale = full(sum(diag(A))) / full(sum(diag(B)));
    fprintf(out, '## Hemisphere %d (n=%d, sigmaScale=%.6g)\n\n', hh, n, sigmaScale);
    % --- comparison arm: plain smallestabs, sigma=0 (known-bug arm) ---
    armPhi = {};
    for arm = 0:numel(TauRelList)
        if arm == 0
            label = 'sigma0-smallestabs';
            solver = @(sv) eigsArm(A, B, K, 0, sv);
        else
            label = sprintf('shift-invert TauRel=%g', TauRelList(arm));
            solver = @(sv) eigsArm(A, B, K, TauRelList(arm)*sigmaScale, sv);
        end
        try
            t0 = tic;
            [Phi1, L1] = solver(ones(n,1));
            t1 = toc(t0);
            [Phi2, L2] = solver(ones(n,1));     % identical start: determinism probe
            [Phi3, L3] = solver(rand(n,1));     % different start: robustness probe
            res  = pencilResidual(A, B, Phi1, L1);
            orth = norm(Phi1' * B * Phi1 - eye(K), 'fro');
            rep12 = subspaceCorr(Phi1, Phi2, B);
            rep13 = subspaceCorr(Phi1, Phi3, B);
            fprintf(out, '- %s: t=%.1fs, maxres=%.3g, orth=%.3g, lambda1=%.3g, det-rep=%.3g, rand-rep=%.3g\n', ...
                label, t1, res, orth, L1(1), 1-rep12, 1-rep13);
            armPhi{end+1} = Phi1; %#ok<AGROW>
        catch err
            fprintf(out, '- %s: FAILED (%s)\n', label, err.message);
            armPhi{end+1} = []; %#ok<AGROW>
        end
    end
    % --- cross-arm agreement among the shift-invert arms ---
    for a = 2:numel(armPhi)
        for b = a+1:numel(armPhi)
            if ~isempty(armPhi{a}) && ~isempty(armPhi{b})
                fprintf(out, '- cross-arm subspace disagreement arms %d vs %d: %.3g\n', ...
                    a-1, b-1, 1 - subspaceCorr(armPhi{a}, armPhi{b}, B));
            end
        end
    end
    fprintf(out, '\n');
end
fprintf(out, 'PINNED: <filled in by the implementer from the numbers above>\n');
fclose(out);
disp('EXPERIMENT DONE');

function [Phi, Lambda] = eigsArm(A, B, K, tau, startVec)
    opts.v0 = startVec; opts.issym = true;
    [Phi, Mu] = eigs(A + tau*B, B, K, 'smallestabs', opts);
    [Lambda, iSort] = sort(diag(Mu) - tau);
    Phi = Phi(:, iSort);
end
function r = pencilResidual(A, B, Phi, Lambda)
    R = A*Phi - B*Phi*diag(Lambda);
    N = (A + B)*Phi;
    r = max(sqrt(sum(R.^2,1)) ./ sqrt(sum(N.^2,1)));
end
function c = subspaceCorr(P1, P2, B)
    % mean squared cosine of principal angles between B-orthonormal bases
    M = P1' * B * P2;
    c = mean(svd(M).^2);
end
```

- [ ] **Step 3: Run it** (`dev/verify/phase0/run_matlab.sh` with the
  script's absolute path; runtime ~5–20 min). Expected: `EXPERIMENT DONE`,
  results file populated with all arms.

- [ ] **Step 4: Pin the parameters.** Edit the `PINNED:` line based on the
  evidence. Decision rule: choose the smallest TauRel whose arms show
  (a) residual and orthonormality <= 1e-10, (b) deterministic-start
  reproducibility disagreement <= 1e-12, (c) cross-arm disagreement with the
  other shift-invert arms <= 1e-10 (i.e., the choice of tau is immaterial —
  then prefer 1e-4 as the middle default). If the sigma=0 arm FAILS or
  disagrees, record it as confirmation of the known bug; if it succeeds and
  agrees on this well-conditioned Galerkin pencil, note that too (the shift
  still ships — user decision — but the record matters for the spec's OQ1).
  If NO arm meets the bar, STOP and report BLOCKED with the table.

(No commits this task.)

---

### Task 2: `tess_eigen` on the clean branch

**Files:**
- Create (worktree): `toolbox/anatomy/tess_eigen.m`
- Create (main checkout): `dev/verify/phase2/test_tess_eigen_sphere.m`

**Interfaces:**
- Consumes: `tess_operators` (Phase 1), pinned `TauRel` from Task 1.
- Produces (later tasks + Phase 3 consume):
  `EigenEntry = tess_eigen(SurfaceFile, OperatorName, 'nModes', K, 'ForceRecompute', 0/1, 'NoSave', 0/1)`
  where `SurfaceFile` is a surface .mat path (absolute or protocol-relative)
  and `OperatorName` is `'Laplace-Beltrami'` or a recipe struct. Returns the
  stored-format entry:
  `.Phi {1x2}` (B-orthonormal, {L,R}, hemi-local sorted rows), `.Lambda {1x2}`
  (ascending), `.nModes`, `.Operator` (recipe: .Name, .Tau, .Assembly),
  `.Solver` (.Method='shift-invert', .TauRel, .SigmaScale [1x2], .StartVector,
  .MatlabVersion, .Date). Cache field on disk: `SurfaceMat.Eigen.LaplaceBeltrami`
  (variant name mapped by local `VariantToField`: strip non-alphanumerics).

- [ ] **Step 1: Write the failing sphere test (main checkout)**

Create `dev/verify/phase2/test_tess_eigen_sphere.m` — a synthetic TWO-sphere
surface (one per "hemisphere") makes the analytic spectrum checkable AND
exercises the atlas split:

```matlab
% TEST_TESS_EIGEN_SPHERE: analytic sphere spectrum + storage round-trip.
scratch = '/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad';
TestFile = fullfile(scratch, 'tess_twosphere_unit.mat');
% --- build: two unit spheres, offset, labeled lh/rh via Structures atlas ---
[Vs, Fs] = tess_sphere(2562);
n1 = size(Vs,1);
TessMat = db_template('surfacemat');
TessMat.Comment  = 'twosphere_unit';
TessMat.Vertices = [Vs; Vs + repmat([3 0 0], n1, 1)];
TessMat.Faces    = [Fs; Fs + n1];
TessMat.Atlas(2).Name = 'Structures';
TessMat.Atlas(2).Scouts(1) = db_template('scout');
TessMat.Atlas(2).Scouts(1).Label = 'lh';  TessMat.Atlas(2).Scouts(1).Region = 'LU';
TessMat.Atlas(2).Scouts(1).Vertices = 1:n1;  TessMat.Atlas(2).Scouts(1).Seed = 1;
TessMat.Atlas(2).Scouts(2) = db_template('scout');
TessMat.Atlas(2).Scouts(2).Label = 'rh';  TessMat.Atlas(2).Scouts(2).Region = 'RU';
TessMat.Atlas(2).Scouts(2).Vertices = n1+1:2*n1;  TessMat.Atlas(2).Scouts(2).Seed = n1+1;
bst_save(TestFile, TessMat, 'v7');
% --- solve K=50 ---
E = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 50);
for hh = 1:2
    Lam = E.Lambda{hh};
    assert(abs(Lam(1)) < 1e-8, 'lambda_1 must be ~0 (constant mode)');
    assert(issorted(Lam), 'eigenvalues must ascend');
    % analytic: lambda = l(l+1) on the unit sphere, multiplicity 2l+1
    lExp = []; for l = 0:6, lExp = [lExp, repmat(l*(l+1), 1, 2*l+1)]; end %#ok<AGROW>
    lExp = lExp(1:50)';
    relErr = abs(Lam - lExp) ./ max(lExp, 1);
    assert(max(relErr) < 0.05, 'sphere spectrum deviates >5%% (max rel err %g)', max(relErr));
end
% --- storage round-trip: field embedded, reusable, truncatable ---
S2 = load(TestFile);
assert(isfield(S2, 'Eigen') && isfield(S2.Eigen, 'LaplaceBeltrami'), 'Eigen field not embedded');
assert(S2.Eigen.LaplaceBeltrami.nModes == 50);
assert(~isempty(S2.History) && any(strcmpi(S2.History(:,2), 'eigen')), 'History row missing');
E2 = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 30);   % reuse + truncate
assert(size(E2.Phi{1}, 2) == 30 && isequal(E2.Phi{1}, E.Phi{1}(:,1:30)), 'truncated reuse failed');
S3 = load(TestFile);
assert(S3.Eigen.LaplaceBeltrami.nModes == 50, 'reuse must not shrink the stored basis');
E3 = tess_eigen(TestFile, 'Laplace-Beltrami', 'nModes', 50, 'ForceRecompute', 1);  % replace slot
S4 = load(TestFile);
assert(numel(fieldnames(S4.Eigen)) == 1, 'one slot per variant');
% --- residual + orthonormality on the recomputed basis ---
[Op, Ms] = tess_operators(TestFile, 'Laplace-Beltrami');
for hh = 1:2
    R = Op{hh}*E3.Phi{hh} - Ms{hh}*E3.Phi{hh}*diag(E3.Lambda{hh});
    N = (Op{hh} + Ms{hh})*E3.Phi{hh};
    assert(max(sqrt(sum(R.^2,1))./sqrt(sum(N.^2,1))) < 1e-10, 'pencil residual too large');
    assert(norm(E3.Phi{hh}'*Ms{hh}*E3.Phi{hh} - eye(50), 'fro') < 1e-10, 'not B-orthonormal');
end
delete(TestFile);
disp('test_tess_eigen_sphere PASSED');
```

- [ ] **Step 2: Run — must FAIL** (`Unrecognized function ... 'tess_eigen'`;
  the clean branch has none).

- [ ] **Step 3: Implement `tess_eigen` (worktree)**

Create `toolbox/anatomy/tess_eigen.m` (canonical header; help documents the
embedded-cache contract, one-slot-per-variant replace semantics, the row
convention pointer to tess_operators, and the mrimask save pattern):

```matlab
function EigenEntry = tess_eigen(SurfaceFile, OperatorName, varargin)
    % --- options ---
    nModes = 400; ForceRecompute = false; NoSave = false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case {'nmodes','k'},    nModes = varargin{i+1};
            case 'forcerecompute',  ForceRecompute = logical(varargin{i+1});
            case 'nosave',          NoSave = logical(varargin{i+1});
            otherwise, error('tess_eigen:badOption', 'Unknown option "%s".', varargin{i});
        end
    end
    % --- recipe normalization ---
    if isstruct(OperatorName)
        Recipe = OperatorName;
    else
        Recipe = struct('Name', OperatorName, 'Tau', []);
    end
    if ~isfield(Recipe, 'Name')
        error('tess_eigen:unknownVariant', 'Operator recipe must carry a .Name field.');
    end
    FieldName = VariantToField(Recipe.Name);
    FullFile  = file_fullpath(SurfaceFile);
    TessMat   = load(FullFile);
    nVertices = size(TessMat.Vertices, 1);
    % --- find-or-create against the embedded cache ---
    if ~ForceRecompute && isfield(TessMat, 'Eigen') && isfield(TessMat.Eigen, FieldName)
        Cached = TessMat.Eigen.(FieldName);
        if RecipeMatches(Cached.Operator, Recipe) && (Cached.nModes >= nModes) ...
                && CacheConsistent(Cached, nVertices)
            EigenEntry = TruncateEntry(Cached, nModes);
            return;
        end
    end
    % --- assemble pencil + solve per hemisphere ---
    [Operator, Mass] = tess_operators(TessMat, Recipe);
    TauRel = 1e-4;   % pinned by the Phase 2 shift experiment
    Phi = cell(1,2); Lambda = cell(1,2); SigmaScale = [0 0];
    for hh = 1:2
        A = Operator{hh}; B = Mass{hh}; n = size(A,1);
        if nModes >= n
            error('tess_eigen:tooManyModes', ...
                'Requested %d modes but hemisphere has only %d vertices.', nModes, n);
        end
        SigmaScale(hh) = full(sum(diag(A))) / full(sum(diag(B)));
        tau = TauRel * SigmaScale(hh);
        opts.v0 = ones(n, 1);        % deterministic start (reproducible runs)
        opts.issym = true;
        [V, Mu, eigsFlag] = eigs(A + tau*B, B, nModes, 'smallestabs', opts);
        if eigsFlag ~= 0
            error('tess_eigen:noConvergence', ...
                'eigs did not converge for hemisphere %d (flag %d).', hh, eigsFlag);
        end
        [lam, iSort] = sort(diag(Mu) - tau);
        V = V(:, iSort);
        % sign convention: largest-magnitude component positive
        for k = 1:nModes
            [~, iMax] = max(abs(V(:,k)));
            if V(iMax,k) < 0, V(:,k) = -V(:,k); end
        end
        Phi{hh} = V; Lambda{hh} = lam;
    end
    % --- build the entry ---
    EigenEntry = struct();
    EigenEntry.Phi     = Phi;
    EigenEntry.Lambda  = Lambda;
    EigenEntry.nModes  = nModes;
    EigenEntry.Operator = struct('Name', Recipe.Name, 'Tau', [], ...
        'Assembly', 'tess_laplacian/tess_massmatrix (galerkin) v1.0');
    EigenEntry.Solver  = struct('Method', 'shift-invert', 'TauRel', TauRel, ...
        'SigmaScale', SigmaScale, 'StartVector', 'ones', ...
        'MatlabVersion', version, 'Date', datestr(now));
    % --- persist (mrimask pattern): append-save Eigen only + History row ---
    if ~NoSave
        if isfield(TessMat, 'Eigen'), s.Eigen = TessMat.Eigen; end
        s.Eigen.(FieldName) = EigenEntry;
        bst_save(FullFile, s, 'v7', 1);
        bst_history('add', FullFile, 'eigen', sprintf( ...
            '%s: %d modes/hemisphere, shift-invert TauRel=%g, tess_eigen v1.0', ...
            Recipe.Name, nModes, TauRel));
        PatchLoadedSurface(FullFile, s.Eigen);
    end
end

function f = VariantToField(name)
    f = regexprep(char(name), '[^A-Za-z0-9]', '');
    if isempty(f) || ~isletter(f(1))
        error('tess_eigen:unknownVariant', 'Cannot map variant "%s" to a field name.', char(name));
    end
end
function ok = RecipeMatches(stored, requested)
    ok = strcmpi(stored.Name, requested.Name);
    if ok && isfield(requested, 'Tau') && ~isempty(requested.Tau)
        ok = isfield(stored, 'Tau') && ~isempty(stored.Tau) && (stored.Tau == requested.Tau);
    end
end
function ok = CacheConsistent(entry, nVertices)
    ok = (size(entry.Phi{1},1) + size(entry.Phi{2},1)) == nVertices;
end
function e = TruncateEntry(e, nModes)
    for hh = 1:2
        e.Phi{hh}    = e.Phi{hh}(:, 1:nModes);
        e.Lambda{hh} = e.Lambda{hh}(1:nModes);
    end
    e.nModes = nModes;
end
function PatchLoadedSurface(FullFile, EigenField)
    % mrimask pattern: update the in-memory copy if this surface is loaded
    global GlobalData;
    if isempty(GlobalData) || ~isfield(GlobalData, 'Surface'), return; end
    iSurf = find(file_compare({GlobalData.Surface.FileName}, file_short(FullFile)));
    if ~isempty(iSurf)
        GlobalData.Surface(iSurf).Eigen = EigenField;
    end
end
```

Adjust `TauRel` to the value Task 1 pinned if it differs from 1e-4. Note
`CacheConsistent` uses total row count = nVertices (LBO scalar; the 4x Dirac
factor arrives with the Dirac variant later).

- [ ] **Step 4: Run the sphere test — must PASS.**

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_eigen.m
git commit -m "Anatomy: Add tess_eigen, embedded per-hemisphere Laplace-Beltrami eigenbases"
```

---

### Task 3: Real-cortex verification (residuals, storage semantics, size)

**Files:**
- Create (main checkout): `dev/verify/phase2/test_tess_eigen_cortex.m`

**Interfaces:**
- Consumes: `tess_eigen` (Task 2), Phase 1 oracle (pencil for residuals),
  the Gate-0 cortex (COPIED — never in place).

- [ ] **Step 1: Write the test**

Create `dev/verify/phase2/test_tess_eigen_cortex.m`:

```matlab
% TEST_TESS_EIGEN_CORTEX: K=400/hemisphere on the real ico5 cortex (a COPY).
scratch = '/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad';
Orig = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
Work = fullfile(scratch, 'tess_cortex_eigen_work.mat');
copyfile(Orig, Work);  fileattrib(Work, '+w');
info0 = dir(Work);
t0 = tic;
E = tess_eigen(Work, 'Laplace-Beltrami', 'nModes', 400);
tSolve = toc(t0);
fprintf('solve time: %.1f s\n', tSolve);
% --- residual + orthonormality vs the INDEPENDENT nxr oracle pencil ---
S = load('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase1/oracle_lbo_sub0002.mat');
for hh = 1:2
    A = S.A{hh}; B = S.B{hh};
    R = A*E.Phi{hh} - B*E.Phi{hh}*diag(E.Lambda{hh});
    N = (A + B)*E.Phi{hh};
    res = max(sqrt(sum(R.^2,1)) ./ sqrt(sum(N.^2,1)));
    orth = norm(E.Phi{hh}'*B*E.Phi{hh} - eye(400), 'fro');
    fprintf('hemi %d: residual(vs nxr pencil)=%.3g, orth=%.3g, lambda1=%.3g, lambda400=%.6g\n', ...
        hh, res, orth, E.Lambda{hh}(1), E.Lambda{hh}(400));
    assert(res < 1e-8, 'residual vs independent nxr pencil too large');   % cross-implementation
    assert(orth < 1e-10, 'not B-orthonormal');
    assert(abs(E.Lambda{hh}(1)) < 1e-6 * E.Lambda{hh}(400), 'zero mode not recovered');
end
% --- storage: file growth, reuse speed, replace semantics ---
info1 = dir(Work);
fprintf('file: %.1f MB -> %.1f MB\n', info0.bytes/1e6, info1.bytes/1e6);
% 2x10242x400 doubles = 65.5 MB raw; -v7 gzip on smooth eigenvectors keeps most of it
assert(info1.bytes > info0.bytes + 40e6, 'Eigen field looks too small for 2x10242x400 doubles');
t1 = tic; E2 = tess_eigen(Work, 'Laplace-Beltrami', 'nModes', 400); tReuse = toc(t1);
fprintf('reuse time: %.2f s\n', tReuse);
assert(tReuse < tSolve / 10, 'cache reuse should be much faster than solving');
assert(isequal(E2.Phi{1}, E.Phi{1}), 'reuse must return the stored basis');
% --- in_tess_bst passes the field through untouched ---
T = in_tess_bst(Work, 0);
assert(isfield(T, 'Eigen') && isfield(T.Eigen, 'LaplaceBeltrami'), 'in_tess_bst dropped Eigen');
delete(Work);
disp('test_tess_eigen_cortex PASSED');
```

- [ ] **Step 2: Run — must PASS** (solve expected ~1–10 min for 2x400
  modes; note the time). The residual is checked against the INDEPENDENT
  nxr-built pencil, so this single number certifies the whole chain
  (primitives + solver) end to end.

- [ ] **Step 3: If the residual bar fails only marginally** (1e-8..1e-6),
  investigate eigs tolerance (`opts.tol`) before touching anything else;
  a genuine convention/structure failure is BLOCKED per the standing rule.

(No commits — harness only; the clean branch is untouched by this task.)

---

### Task 4: Eigenmode viewer demo (Gate 2 live artifact)

**Files:**
- Create (main checkout): `dev/verify/phase2/make_eigenmode_movie.m`

**Interfaces:**
- Consumes: the isolated protocol (clean-branch-readable), `tess_eigen`.
- Produces: a standard Brainstorm results file "LBO eigenmodes (mode index as
  time)" under the isolated protocol's sub-0002 @intra study — displayable
  with stock Brainstorm 3D viewers, zero new GUI code — plus snapshot PNGs
  for the Gate 2 report.

- [ ] **Step 1: Write the script**

Create `dev/verify/phase2/make_eigenmode_movie.m`:

```matlab
% MAKE_EIGENMODE_MOVIE: writes the first 100 LBO eigenmodes as a Brainstorm
% results file (Time axis = mode index) in the isolated protocol, and saves
% snapshot PNGs of selected modes for the Gate 2 report.
nShow = 100;
SurfaceFile = 'sub-0002/tess_cortex_pial_low.mat';   % protocol-relative
E = tess_eigen(SurfaceFile, 'Laplace-Beltrami', 'nModes', max(400, nShow));
Full = file_fullpath(SurfaceFile);
T = load(Full, 'Vertices');
nV = size(T.Vertices, 1);
[rH, lH] = tess_hemisplit(in_tess_bst(SurfaceFile));
Maps = zeros(nV, nShow);
Maps(sort(lH(:)), :) = E.Phi{1}(:, 1:nShow);
Maps(sort(rH(:)), :) = E.Phi{2}(:, 1:nShow);
% results structure (source map on cortex, mode index as "time")
ResultsMat = db_template('resultsmat');
ResultsMat.Comment       = sprintf('LBO eigenmodes 1-%d (mode=frame)', nShow);
ResultsMat.ImageGridAmp  = Maps;
ResultsMat.Time          = 1:nShow;
ResultsMat.SurfaceFile   = SurfaceFile;
ResultsMat.HeadModelType = 'surface';
ResultsMat.Function      = 'eigenmodes';
ResultsMat.nComponents   = 1;
% save into sub-0002 @intra
[sSubject, iSubject] = bst_get('Subject', 'sub-0002');
sStudy = bst_get('AnalysisIntraStudy', iSubject);
[~, iStudy] = bst_get('Study', sStudy.FileName);
OutFile = bst_process('GetNewFilename', bst_fileparts(sStudy.FileName), 'results_lbo_eigenmodes');
bst_save(OutFile, ResultsMat, 'v7');
db_add_data(iStudy, OutFile, ResultsMat);
fprintf('EIGENMODE MOVIE: %s\n', file_short(OutFile));
% snapshots of modes 2, 10, 50 for the report
outDir = '/Users/diellorbasha/workspace/research/code/brainstorm3/dev/verify/phase2';
hFig = view_surface_data(SurfaceFile, file_short(OutFile));
for m = [2 10 50]
    panel_time('SetCurrentTime', m);
    frame = bst_call(@out_figure_image, hFig);
    imwrite(frame, fullfile(outDir, sprintf('eigenmode_%03d.png', m)));
end
close(hFig);
disp('make_eigenmode_movie DONE');
```

- [ ] **Step 2: Run via the isolated runner.** The tess_eigen call here also
  permanently embeds the K=400 basis in the ISOLATED protocol's cortex (this
  is intentional — it becomes the standing demo protocol; the Phase 1 oracle
  is unaffected since its A/B/vH are already extracted into the oracle .mat;
  note this in the report). If `view_surface_data`/snapshot steps fail
  headless, save the results file anyway, skip PNGs, and record that the
  visual check happens in the user's Gate-2 GUI session instead — the movie
  file is the deliverable, the PNGs are a bonus.

- [ ] **Step 3: Verify**: re-run `dev/verify/phase1/test_tess_operators.m`
  (parity suite) to prove the embedded Eigen field on the isolated cortex
  breaks nothing downstream (in particular the oracle-based tests still
  load the surface fine with the extra field present).

(No commits — harness only.)

---

### Task 5: Gate 2 (= Milestone 1) report + harness commit + final review

**Files:**
- Create (main checkout): `dev/verify/phase2/gate2_report.md`
- Commit (main checkout, lab branch): `dev/verify/phase2/*` + runner edit

- [ ] **Step 1: Re-run the two clean-branch tests, capturing logs**

```bash
cd ~/workspace/research/code/brainstorm3
for t in test_tess_eigen_sphere test_tess_eigen_cortex; do
  dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase2/$t.m" 2>&1 | tee "dev/verify/phase2/$t.log"
done
```

- [ ] **Step 2: Write `dev/verify/phase2/gate2_report.md`** — concrete:
  clean-branch commit SHA + subject; the pinned solver parameters and the
  experiment table summary (incl. what happened to the sigma=0 arm — this
  documents Open Question 1's resolution); sphere-spectrum max rel err;
  cortex residual/orthonormality numbers vs the independent nxr pencil;
  solve vs reuse timings and file-size growth; the eigenmode-movie file path
  + PNGs (or the headless-skip note); deviations.

- [ ] **Step 3: Commit harness on the lab branch**

```bash
cd ~/workspace/research/code/brainstorm3
git add -f dev/verify/phase2/
git add dev/verify/phase0/run_matlab.sh
git commit -m "verify(phase2): Gate 2 harness — shift experiment, tess_eigen tests, eigenmode demo

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

- [ ] **Step 4: Present Gate 2 = Milestone 1 to the user.** Show:
  gate2_report.md, `git -C ~/workspace/research/code/brainstorm3-clean log --oneline master..`
  (8 commits), the eigenmode snapshots, and how to view the eigenmode movie
  in the isolated GUI (`dev/verify/phase0/launch_gui_clean.sh`, open
  sub-0002 > @intra > "LBO eigenmodes"). Note the two follow-on decisions:
  Open Question 2 (Dirac variant) and Phase 3 scoping (filtering framework
  vs Dirac kernel projection first). Do NOT push.
