# Phase 1 — LBO Operator Primitives Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add pure-MATLAB `tess_massmatrix` (Galerkin mass), `tess_laplacian` (cotangent stiffness), and compute-only `tess_operators` (per-hemisphere Laplace-Beltrami pencil) to the clean branch, verified to near machine precision against the nxr-built operators (Gate 1).

**Architecture:** Two I/O-free primitives + one pencil-layer dispatcher, per the
approved spec (§3.1–3.2). Parity oracle generated ONCE on the dev side (nxr,
main checkout) from the real sub-0002 ico5 cortex; clean-branch tests compare
against the saved oracle file. Dirac is OUT of this plan — gated on the
pinned Open Question 2 (exact Dirac variant), to be planned separately after
that decision.

**Tech Stack:** MATLAB R2023b (`matlab -batch`), Brainstorm, git worktree `~/workspace/research/code/brainstorm3-clean`.

**Spec:** `docs/superpowers/specs/2026-08-07-cortical-flow-distillation-design.md` (§3.1, §3.2, §4 Phase 1, §6)

## Global Constraints

- Clean-branch (worktree) commits: upstream Brainstorm style subjects, NO
  Claude trailers, explicit path staging only (never `git add -A`), one
  logical unit per commit.
- Worktree = `~/workspace/research/code/brainstorm3-clean`, branch
  `feature/cortical-flow-core` (currently 4 commits ahead of master). Never
  switch branches in the main checkout (`~/workspace/research/code/brainstorm3`).
- Harness for THIS phase: `dev/verify/phase1/` in the main checkout
  (gitignored via `.git/info/exclude` `dev/` — use `git add -f` when
  committing). Reuse `dev/verify/phase0/run_matlab.sh` for clean-branch runs;
  Task 1 extends it to addpath phase1 too.
- Oracle runs (nxr) execute in the MAIN checkout's MATLAB (dev code, real
  user config) — the user's Brainstorm GUI must be closed during those runs.
  Clean-branch runs stay fully isolated (phase0 runner).
- The parity surface is the user's protocol cortex:
  `/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat`
  — loaded DIRECTLY from file (never through the DB), read-only.
- Numerical convention authority: the ORACLE. Task 1 documents nxr's exact
  sign/scale conventions; Tasks 2–3 match them exactly (units are meters —
  Brainstorm vertex coordinates — so mass entries are m²).
- Parity tolerance: `norm(ours - oracle, 'fro') / norm(oracle, 'fro') <= 1e-12`
  per hemisphere per matrix. A convention mismatch (sign, lumping, scaling)
  is NOT a tolerance failure to paper over — report BLOCKED with the observed
  relationship (e.g. "ours == -oracle", "ours == diag(sum(oracle,2))").
- Hemisphere order everywhere: index 1 = LEFT (`lH`), index 2 = RIGHT (`rH`),
  local vertex order = `tess_hemisplit` output order (sorted global indices)
  — matches the PoC.
- MATLAB runs: Bash timeout 600000 ms; scratch files under
  `/private/tmp/claude-501/-Users-diellorbasha-workspace-research-code-brainstorm3/8da2056d-09c9-4f92-8237-a34a82e5f5e8/scratchpad/`.

---

### Task 1: Parity oracle (dev side) + phase1 harness scaffold

**Files:**
- Create (main checkout): `dev/verify/phase1/make_oracle_lbo.m`
- Create (by running it): `dev/verify/phase1/oracle_lbo_sub0002.mat`
- Modify (main checkout): `dev/verify/phase0/run_matlab.sh` (addpath phase1)

**Interfaces:**
- Produces: `oracle_lbo_sub0002.mat` with variables per hemisphere h∈{1=L,2=R}:
  `A{h}` (nxr laplacian/cotan, sparse [nVh x nVh]), `B{h}` (nxr mass/galerkin),
  `vH{h}` (global vertex indices, tess_hemisplit order), plus `meta` struct
  (SurfaceFile path, nVertices, nxr version, date, and a `conventions` field
  filled in Step 4). All later tasks consume this file read-only.

- [ ] **Step 1: Write the oracle generator (main checkout)**

Create `dev/verify/phase1/make_oracle_lbo.m`:

```matlab
% MAKE_ORACLE_LBO: nxr-built LBO pencil for the parity surface (dev side, run ONCE).
SurfaceFile = '/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db/omega-tutorial-cortical-flow/anat/sub-0002/tess_cortex_pial_low.mat';
OutFile = fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat');
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute unavailable: %s', errMsg);
TessMat = load(SurfaceFile);
[rH, lH, isConn] = tess_hemisplit(TessMat);
assert(~isConn, 'hemispheres connected');
hemis = {lH(:), rH(:)};
Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces);
nVtot = size(Vtx, 1);
A = cell(1,2); B = cell(1,2); vH = cell(1,2);
for hh = 1:2
    v = hemis{hh};
    isV = false(nVtot,1); isV(v) = true;
    fMask = all(isV(Fcs), 2);
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    h = nxr_compute('create', Vtx(v,:), mapV(Fcs(fMask,:)));
    A{hh} = nxr_compute('operators', h, 'laplacian', 'cotan');
    B{hh} = nxr_compute('operators', h, 'mass', 'galerkin');
    nxr_compute('delete', h);
    vH{hh} = v;
end
nxrVer = ''; try, nxrVer = nxr_compute('version'); catch, end %#ok<CTCH>
meta = struct('SurfaceFile', SurfaceFile, 'nVertices', nVtot, ...
              'NxrVersion', nxrVer, 'Date', datestr(now), 'conventions', '');
save(OutFile, 'A', 'B', 'vH', 'meta', '-v7');
fprintf('ORACLE WRITTEN: %s  (nV=%d, L=%d, R=%d)\n', OutFile, nVtot, numel(vH{1}), numel(vH{2}));
```

Note: if `nxr_compute('delete', h)` is not a valid nxr verb, drop that line
(check `help nxr_compute` or the dev tess_operators, which does not delete
handles either).

- [ ] **Step 2: Run it with the MAIN checkout's MATLAB (dev code, nxr available)**

User's Brainstorm GUI must be closed. From the main checkout:

```bash
cd ~/workspace/research/code/brainstorm3
/Applications/MATLAB_R2023b.app/bin/matlab -batch "cd('$PWD'); brainstorm server; run('$PWD/dev/verify/phase1/make_oracle_lbo.m'); brainstorm stop; exit(0);"
```

Expected: `ORACLE WRITTEN: ... (nV=20484, L=10242, R=10242)`.

- [ ] **Step 3: Characterize the oracle's numerical conventions**

Run a short probe (same MATLAB invocation pattern, or extend Step 1's script)
and record results:

```matlab
S = load('dev/verify/phase1/oracle_lbo_sub0002.mat');
for hh = 1:2
    Ah = S.A{hh}; Bh = S.B{hh};
    fprintf('hemi %d: A sym err %g | A*1 max %g | diagA sign %d | min eig sign(probe) %g\n', ...
        hh, norm(Ah-Ah','fro'), max(abs(Ah*ones(size(Ah,1),1))), sign(full(Ah(1,1))), ...
        min(eig(full(Ah(1:200,1:200)))));
    fprintf('hemi %d: B sym err %g | sum(B(:)) = %g (total area m^2) | B offdiag sign %d\n', ...
        hh, norm(Bh-Bh','fro'), full(sum(Bh(:))), sign(max(full(Bh(1,2:end)))));
end
```

Record in `meta.conventions` (re-save the .mat) and in your report: A symmetric
yes/no, A*1≈0 yes/no, diagonal sign (+ = PSD stiffness convention, − =
negative-semidefinite Laplacian), B total = plausible cortex hemisphere area
in m² (~0.01 m² order), B consistent-Galerkin (off-diagonals positive) vs
lumped (diagonal). These EXACT conventions are what Tasks 2–3 implement.

- [ ] **Step 4: Extend the phase0 runner to addpath phase1**

In `dev/verify/phase0/run_matlab.sh`, change the final `-batch` string:
`addpath('$HARNESS')` → `addpath('$HARNESS'); addpath('$HARNESS/../phase1')`.

- [ ] **Step 5: Sanity-run one phase0 test through the modified runner**

```bash
cd ~/workspace/research/code/brainstorm3
dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase0/test_tess_repair_unit.m"
```

Expected: still `test_tess_repair_unit PASSED`. (No commits this task; the
phase1 harness is committed in Task 5.)

---

### Task 2: `tess_massmatrix` (clean branch)

**Files:**
- Create (worktree): `toolbox/anatomy/tess_massmatrix.m`
- Create (main checkout): `dev/verify/phase1/test_tess_massmatrix.m`

**Interfaces:**
- Consumes: `oracle_lbo_sub0002.mat` (Task 1).
- Produces: `B = tess_massmatrix(Vertices, Faces)` — sparse [nV x nV]
  Galerkin mass, no I/O, no validation side effects. Task 4 calls it
  per hemisphere.

- [ ] **Step 1: Write the failing test (main checkout)**

Create `dev/verify/phase1/test_tess_massmatrix.m`:

```matlab
% TEST_TESS_MASSMATRIX: analytic single-triangle + sphere-area + oracle parity.
% --- analytic: unit right triangle, Area = 0.5 ---
V = [0 0 0; 1 0 0; 0 1 0]; F = [1 2 3];
B = tess_massmatrix(V, F);
Bexp = 0.5 * [1/6 1/12 1/12; 1/12 1/6 1/12; 1/12 1/12 1/6];
assert(norm(full(B) - Bexp, 'fro') < 1e-15, 'single-triangle Galerkin entries wrong');
% --- sphere: total mass = surface area ---
[Vs, Fs] = tess_sphere(2562);          % unit sphere approximation
Bs = tess_massmatrix(Vs, Fs);
v1 = Vs(Fs(:,1),:); v2 = Vs(Fs(:,2),:); v3 = Vs(Fs(:,3),:);
totalArea = sum(0.5 * sqrt(sum(cross(v2-v1, v3-v1, 2).^2, 2)));
assert(abs(full(sum(Bs(:))) - totalArea) < 1e-12 * totalArea, 'total mass ~= mesh area');
assert(norm(Bs - Bs', 'fro') == 0, 'mass must be exactly symmetric');
% --- parity vs nxr oracle on the real cortex ---
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
T = load(S.meta.SurfaceFile);
Vtx = double(T.Vertices); Fcs = double(T.Faces); nVtot = size(Vtx,1);
for hh = 1:2
    v = S.vH{hh};
    isV = false(nVtot,1); isV(v) = true;
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    Bh = tess_massmatrix(Vtx(v,:), mapV(Fcs(all(isV(Fcs),2),:)));
    relErr = norm(Bh - S.B{hh}, 'fro') / norm(S.B{hh}, 'fro');
    fprintf('hemi %d mass parity rel err: %g\n', hh, relErr);
    assert(relErr <= 1e-12, 'mass parity failed (hemi %d): %g', hh, relErr);
end
disp('test_tess_massmatrix PASSED');
```

If Task 1's convention probe showed nxr's galerkin is NOT the consistent-FEM
(Area/6, Area/12) form, adjust the ANALYTIC expectation to nxr's documented
convention (from meta.conventions) BEFORE implementing — the oracle is the
authority; the analytic block must encode the same convention.

- [ ] **Step 2: Run — must FAIL (function absent)**

```bash
cd ~/workspace/research/code/brainstorm3
dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase1/test_tess_massmatrix.m"
```

Expected: `Unrecognized function ... 'tess_massmatrix'`.

- [ ] **Step 3: Implement (worktree)**

Create `toolbox/anatomy/tess_massmatrix.m` with a standard Brainstorm header
(copyright block copied from a neighboring toolbox/anatomy function, author
Diellor Basha 2026, concise help text stating: Galerkin/consistent FEM mass,
units follow input coordinates squared, mesh assumed triangular):

```matlab
function B = tess_massmatrix(Vertices, Faces)
    nVertices = size(Vertices, 1);
    v1 = Vertices(Faces(:,1),:); v2 = Vertices(Faces(:,2),:); v3 = Vertices(Faces(:,3),:);
    faceArea = 0.5 * sqrt(sum(cross(v2 - v1, v3 - v1, 2).^2, 2));
    iRow = [Faces(:,1); Faces(:,2); Faces(:,3); ...
            Faces(:,1); Faces(:,2); Faces(:,3); ...
            Faces(:,2); Faces(:,3); Faces(:,1)];
    iCol = [Faces(:,1); Faces(:,2); Faces(:,3); ...
            Faces(:,2); Faces(:,3); Faces(:,1); ...
            Faces(:,1); Faces(:,2); Faces(:,3)];
    sVal = [faceArea/6; faceArea/6; faceArea/6; ...
            faceArea/12; faceArea/12; faceArea/12; ...
            faceArea/12; faceArea/12; faceArea/12];
    B = sparse(iRow, iCol, sVal, nVertices, nVertices);
end
```

(Adapt entries if the oracle convention differs — see Step 1 note.)

- [ ] **Step 4: Run the test — must PASS**

Same command as Step 2. Expected `test_tess_massmatrix PASSED` with parity
rel err printed (~1e-16). If parity fails with a STRUCTURED discrepancy
(pure sign flip, exact lumping, constant scale), report BLOCKED with the
observed relationship — do not fudge tolerances.

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_massmatrix.m
git commit -m "Anatomy: Add tess_massmatrix, Galerkin mass matrix for triangular surfaces"
```

---

### Task 3: `tess_laplacian` (clean branch)

**Files:**
- Create (worktree): `toolbox/anatomy/tess_laplacian.m`
- Create (main checkout): `dev/verify/phase1/test_tess_laplacian.m`

**Interfaces:**
- Consumes: `oracle_lbo_sub0002.mat`.
- Produces: `A = tess_laplacian(Vertices, Faces)` — sparse [nV x nV]
  cotangent stiffness in the ORACLE's sign convention. Task 4 calls it per
  hemisphere.

- [ ] **Step 1: Write the failing test (main checkout)**

Create `dev/verify/phase1/test_tess_laplacian.m` (written for the PSD
convention `A(i,i) = +sum(w)`, `A(i,j) = -w_ij`; flip expectations if Task
1's probe recorded the negative convention):

```matlab
% TEST_TESS_LAPLACIAN: analytic right triangle + null space + oracle parity.
V = [0 0 0; 1 0 0; 0 1 0]; F = [1 2 3];
A = tess_laplacian(V, F);
% right angle at v1 (cot 0), 45 deg at v2,v3 (cot 1):
% w(edge v2v3) = cot(v1)/2 = 0 ; w(v1v3) = cot(v2)/2 = 0.5 ; w(v1v2) = cot(v3)/2 = 0.5
Aexp = [1 -0.5 -0.5; -0.5 0.5 0; -0.5 0 0.5];
assert(norm(full(A) - Aexp, 'fro') < 1e-14, 'right-triangle cotan weights wrong');
% --- sphere checks ---
[Vs, Fs] = tess_sphere(2562);
As = tess_laplacian(Vs, Fs);
assert(norm(As - As', 'fro') < 1e-14, 'stiffness must be symmetric');
assert(max(abs(As * ones(size(As,1),1))) < 1e-12, 'constants must be in the null space');
e = eigs(As, 3, 'smallestreal');
assert(min(e) > -1e-10, 'stiffness must be positive semidefinite');
% --- parity vs nxr oracle on the real cortex ---
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
T = load(S.meta.SurfaceFile);
Vtx = double(T.Vertices); Fcs = double(T.Faces); nVtot = size(Vtx,1);
for hh = 1:2
    v = S.vH{hh};
    isV = false(nVtot,1); isV(v) = true;
    mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
    Ah = tess_laplacian(Vtx(v,:), mapV(Fcs(all(isV(Fcs),2),:)));
    relErr = norm(Ah - S.A{hh}, 'fro') / norm(S.A{hh}, 'fro');
    fprintf('hemi %d stiffness parity rel err: %g\n', hh, relErr);
    assert(relErr <= 1e-12, 'stiffness parity failed (hemi %d): %g', hh, relErr);
end
disp('test_tess_laplacian PASSED');
```

- [ ] **Step 2: Run — must FAIL (function absent)**

```bash
cd ~/workspace/research/code/brainstorm3
dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase1/test_tess_laplacian.m"
```

- [ ] **Step 3: Implement (worktree)**

Create `toolbox/anatomy/tess_laplacian.m` (standard header as in Task 2;
help text: cotangent stiffness, natural/Neumann boundary, negative weights
from obtuse triangles are kept by design):

```matlab
function A = tess_laplacian(Vertices, Faces)
    nVertices = size(Vertices, 1);
    % Edge vectors opposite each corner: eK is the edge NOT touching corner K
    e1 = Vertices(Faces(:,3),:) - Vertices(Faces(:,2),:);
    e2 = Vertices(Faces(:,1),:) - Vertices(Faces(:,3),:);
    e3 = Vertices(Faces(:,2),:) - Vertices(Faces(:,1),:);
    doubleArea = sqrt(sum(cross(e1, e2, 2).^2, 2));   % |e1 x e2| = 2*faceArea
    % cot at corner K = -(eI . eJ) / (2*Area), {I,J} the edges adjacent to K
    cot1 = -sum(e2 .* e3, 2) ./ doubleArea;
    cot2 = -sum(e3 .* e1, 2) ./ doubleArea;
    cot3 = -sum(e1 .* e2, 2) ./ doubleArea;
    % Half-cotangent weight accumulates on the edge OPPOSITE the corner
    iRow = [Faces(:,2); Faces(:,3); Faces(:,3); Faces(:,1); Faces(:,1); Faces(:,2)];
    iCol = [Faces(:,3); Faces(:,2); Faces(:,1); Faces(:,3); Faces(:,2); Faces(:,1)];
    sVal = -0.5 * [cot1; cot1; cot2; cot2; cot3; cot3];
    A = sparse(iRow, iCol, sVal, nVertices, nVertices);
    A = A - spdiags(sum(A, 2), 0, nVertices, nVertices);   % rows sum to zero
end
```

(If the oracle uses the negative convention, negate `A` at the end instead —
one line, matching Step 1's adjusted expectations.)

- [ ] **Step 4: Run the test — must PASS**

Expected `test_tess_laplacian PASSED`, parity ~1e-16. Same BLOCKED rule as
Task 2 for structured mismatches.

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_laplacian.m
git commit -m "Anatomy: Add tess_laplacian, cotangent stiffness matrix for triangular surfaces"
```

---

### Task 4: `tess_operators` — compute-only pencil layer (clean branch)

**Files:**
- Create (worktree): `toolbox/anatomy/tess_operators.m`
- Create (main checkout): `dev/verify/phase1/test_tess_operators.m`

**Interfaces:**
- Consumes: `tess_massmatrix`, `tess_laplacian`, upstream `tess_hemisplit`.
- Produces (Phase 2 `tess_eigen` will consume this exact signature):
  `[Operator, Mass, vH] = tess_operators(Surface, OperatorName)` where
  `Surface` is a surface-file path OR a loaded surface struct, `OperatorName`
  is `'Laplace-Beltrami'` (string) OR a recipe struct with field `.Name`
  (spec §3.4 reproducibility invariant: `tess_operators(SurfaceFile,
  S.Eigen.<v>.Operator)` rebuilds the pencil). Returns 1x2 cells
  `{L, R}`: `Operator{h}` stiffness, `Mass{h}` Galerkin mass, `vH{h}` global
  vertex indices in hemisphere-local row order. Unknown variant names error
  with `tess_operators:unknownVariant` (Dirac deliberately NOT accepted yet —
  pending the spec's Open Question 2).

- [ ] **Step 1: Write the failing test (main checkout)**

Create `dev/verify/phase1/test_tess_operators.m`:

```matlab
% TEST_TESS_OPERATORS: pencil assembly on the real cortex + guards.
S = load(fullfile(fileparts(mfilename('fullpath')), 'oracle_lbo_sub0002.mat'));
SurfaceFile = S.meta.SurfaceFile;
% --- happy path: parity of the assembled pencil, both hemispheres ---
[Op, Ms, vH] = tess_operators(SurfaceFile, 'Laplace-Beltrami');
for hh = 1:2
    assert(isequal(vH{hh}(:), S.vH{hh}(:)), 'hemisphere vertex order mismatch');
    assert(norm(Op{hh} - S.A{hh}, 'fro') / norm(S.A{hh}, 'fro') <= 1e-12, 'stiffness parity');
    assert(norm(Ms{hh} - S.B{hh}, 'fro') / norm(S.B{hh}, 'fro') <= 1e-12, 'mass parity');
end
% --- recipe-struct call (reproducibility invariant) ---
recipe = struct('Name', 'Laplace-Beltrami', 'Tau', []);
[Op2, Ms2] = tess_operators(SurfaceFile, recipe);
assert(isequal(Op2{1}, Op{1}) && isequal(Ms2{2}, Ms{2}), 'recipe call must reproduce pencil');
% --- guard: unknown variant ---
ok = false;
try, tess_operators(SurfaceFile, 'Dirac'); catch err, ok = strcmp(err.identifier, 'tess_operators:unknownVariant'); end
assert(ok, 'Dirac must error with unknownVariant until the variant decision');
% --- guard: missing Structures atlas ---
T = load(SurfaceFile); T2 = rmfield(T, 'Atlas');
ok = false;
try, tess_operators(T2, 'Laplace-Beltrami'); catch err, ok = strcmp(err.identifier, 'tess_operators:noHemisphereLabels'); end
assert(ok, 'missing atlas must raise noHemisphereLabels');
% --- guard: non-manifold hemisphere (duplicate a face) ---
T3 = load(SurfaceFile); T3.Faces = [T3.Faces; T3.Faces(1, [2 1 3])];
ok = false;
try, tess_operators(T3, 'Laplace-Beltrami'); catch err, ok = strcmp(err.identifier, 'tess_operators:nonManifold'); end
assert(ok, 'non-manifold input must raise nonManifold');
disp('test_tess_operators PASSED');
```

- [ ] **Step 2: Run — must FAIL** (same runner command pattern; expected
`Unrecognized function ... 'tess_operators'` — note the clean branch has no
tess_operators; only dev does).

- [ ] **Step 3: Implement (worktree)**

Create `toolbox/anatomy/tess_operators.m` (standard header; help documents
the compute-only contract, the {L,R} order, the hemi-local row convention =
tess_hemisplit sorted order, and that operators are cheap and never stored):

```matlab
function [Operator, Mass, vH] = tess_operators(Surface, OperatorName)
    % --- load surface (file path or struct) ---
    if ischar(Surface)
        TessMat = load(file_fullpath(Surface));
    else
        TessMat = Surface;
    end
    % --- variant (string or recipe struct with .Name) ---
    if isstruct(OperatorName)
        Variant = OperatorName.Name;
    else
        Variant = OperatorName;
    end
    if ~strcmpi(Variant, 'Laplace-Beltrami')
        error('tess_operators:unknownVariant', ...
            'Unknown operator variant "%s" (supported: Laplace-Beltrami).', char(Variant));
    end
    % --- guard: Structures atlas with L/R labels (message pattern from the PoC) ---
    hasLabels = false;
    if isfield(TessMat, 'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts  = TessMat.Atlas(iStruct).Scouts;
            labels  = {scouts.Label};
            regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasLabels = (any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'))) && ...
                        (any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R')));
        end
    end
    if ~hasLabels
        error('tess_operators:noHemisphereLabels', ...
            ['Surface has no Structures atlas with left/right hemisphere labels ' ...
             '(required for the atlas-based hemisphere split).']);
    end
    % --- hemisphere split (atlas-based, never conncomp) ---
    [rH, lH, isConnected] = tess_hemisplit(TessMat);
    if isConnected
        error('tess_operators:connectedHemispheres', ...
            'Hemispheres are connected; operators require independent hemisphere meshes.');
    end
    hemis = {lH(:), rH(:)};
    Vertices = double(TessMat.Vertices);
    Faces    = double(TessMat.Faces);
    nVtot    = size(Vertices, 1);
    Operator = cell(1,2); Mass = cell(1,2); vH = cell(1,2);
    for hh = 1:2
        v = hemis{hh};
        isV = false(nVtot,1); isV(v) = true;
        mapV = zeros(nVtot,1); mapV(v) = 1:numel(v);
        FacesLocal = mapV(Faces(all(isV(Faces),2), :));
        CheckManifold(FacesLocal, numel(v));
        Operator{hh} = tess_laplacian(Vertices(v,:), FacesLocal);
        Mass{hh}     = tess_massmatrix(Vertices(v,:), FacesLocal);
        vH{hh}       = v;
    end
end

function CheckManifold(Faces, nVertices)
    % Closed-2-manifold guard: every undirected edge in exactly 2 faces,
    % every directed edge unique (consistent orientation), no unused vertices.
    directedEdges = [Faces(:,[1 2]); Faces(:,[2 3]); Faces(:,[3 1])];
    if size(unique(directedEdges, 'rows'), 1) ~= size(directedEdges, 1)
        error('tess_operators:nonManifold', ...
            'Hemisphere mesh has inconsistent orientation or duplicate faces.');
    end
    [~, ~, iEdge] = unique(sort(directedEdges, 2), 'rows');
    edgeCount = accumarray(iEdge, 1);
    if any(edgeCount ~= 2)
        error('tess_operators:nonManifold', ...
            'Hemisphere mesh is not a closed 2-manifold (%d boundary/non-manifold edges).', ...
            sum(edgeCount ~= 2));
    end
    if any(~ismember(1:nVertices, Faces(:)))
        error('tess_operators:nonManifold', 'Hemisphere mesh has unreferenced vertices.');
    end
end
```

Note: the duplicate-face corruption in the test creates a repeated directed
edge → caught by the orientation check. If the icosphere cortex legitimately
fails `edgeCount ~= 2` (it should not — icosphere is closed), STOP and
report; do not weaken the guard.

- [ ] **Step 4: Run the test — must PASS**

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_operators.m
git commit -m "Anatomy: Add tess_operators, per-hemisphere Laplace-Beltrami operator pencil"
```

---

### Task 5: Gate 1 report + harness commit + final review

**Files:**
- Create (main checkout): `dev/verify/phase1/gate1_report.md`
- Commit (main checkout, lab branch): `dev/verify/phase1/*` + the run_matlab.sh edit

**Interfaces:**
- Consumes: all Task 2–4 tests.
- Produces: Gate 1 evidence for the user.

- [ ] **Step 1: Re-run the three clean-branch tests, capturing logs**

```bash
cd ~/workspace/research/code/brainstorm3
for t in test_tess_massmatrix test_tess_laplacian test_tess_operators; do
  dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase1/$t.m" 2>&1 | tee "dev/verify/phase1/$t.log"
done
grep -c "PASSED" dev/verify/phase1/*.log
```

Expected: all three logs PASSED.

- [ ] **Step 2: Write `dev/verify/phase1/gate1_report.md`**

Concrete content from logs and the oracle meta: clean-branch commit SHAs +
subjects (3 new commits), per-test result + parity rel-err numbers, the
oracle's recorded conventions (sign, Galerkin form, total areas in m²), the
guard behaviors verified, and any deviations.

- [ ] **Step 3: Commit harness on the lab branch (main checkout)**

```bash
cd ~/workspace/research/code/brainstorm3
git add -f dev/verify/phase1/ 
git add dev/verify/phase0/run_matlab.sh
git commit -m "verify(phase1): Gate 1 harness — nxr parity oracle + LBO primitive tests

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

(Exclude nothing here — the oracle .mat is a few MB of sparse matrices and
belongs with the evidence. If it exceeds ~50 MB, gitignore it and record its
SHA-256 in the report instead.)

- [ ] **Step 4: Present Gate 1 to the user**

Show: gate1_report.md, `git -C ~/workspace/research/code/brainstorm3-clean log --oneline master..` (7 commits total), the parity numbers, and the two follow-on decisions now due: Open Question 2 (Dirac variant — blocks the Dirac operator task) and Phase 2 planning (tess_eigen + Open Question 1 solver experiments). Do NOT push anything.
