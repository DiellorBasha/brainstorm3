# Sign-Ambiguity Continuity Test Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a two-component offline analysis under `dev/benchmarks/sign_ambiguity/` that demonstrates the connection-Laplacian (Fiedler) frame yields a continuous, gauge-consistent current phase across cortical folds where the constrained source map's sign flips.

**Architecture:** Five focused helper functions plus one driver, all in `dev/benchmarks/sign_ambiguity/`. Component 1 (`sa_smoothness`) is data-free and operates on the eigenvector phases and the connection Laplacian's off-diagonal transport. Component 2 (`sa_continuity`) applies the unconstrained MEG kernel and compares phase continuity across facing sulcal-wall vertex pairs under three frames. Shared frame construction lives in `sa_frames`; sulcal-wall pair detection in `sa_sulcal_walls`; rendering in `sa_figures`. Each helper has a focused smoke test under `dev/tests/` that SKIPs when prerequisites (real cortex, kernel) are absent.

**Tech Stack:** MATLAB, Brainstorm toolbox (`in_tess_bst`, `in_bst_results`, `tess_*`), nxr-compute plugin (`nxr.manifold.context`, `vertexFrame`), existing M2/M3 functions (`bst_conn_eigenmodes_ensure`, `bst_conn_phase`, `bst_tangent_face2vertex`). Tests run via the MATLAB MCP (`run_matlab_test_file` / `run_matlab_file`).

**Spec:** `dev/sign_ambiguity_continuity_test.md`

**Reference facts the engineer needs (verified against the codebase):**
- `bst_conn_eigenmodes_ensure(SurfaceFile)` returns a `ConnEig` struct with fields: `Vectors` `[nV×nK]` complex, `Values`, `Component` `[nK×1]` (which connected component / hemisphere each mode belongs to), `CompRank` `[nK×1]` (within-component rank; `CompRank==1` is the Fiedler mode), `Order`, `nComponents`, `MassMatrix`, `ConnLaplacian` `[nV×nV]` complex sparse (the operator `K`).
- A connection eigenvector column `z = ConnEig.Vectors(:,col)` is **nonzero only on its own component** (block structure is exact-zero off it). Find a component's Fiedler column with `find(ConnEig.Component==c & ConnEig.CompRank==1, 1)`.
- The connection Laplacian `K` (= `ConnEig.ConnLaplacian`) has off-diagonal entries whose argument encodes the Levi-Civita parallel transport: for an edge `(i,j)`, the transport is `R_ij = -K(i,j)/|K(i,j)|` (the off-diagonal is `-w_ij·R_ij`). A field parallel-transported along the edge satisfies `z(i) = R_ij·z(j)`, so the **covariant phase increment** `wrap(arg(z_i) − arg(z_j) − arg(R_ij))` is ~0 for a smooth field. This is the gauge-correct `δf`.
- `R = bst_conn_phase(ConnEig, vFrame, 'Rank',1, 'FsFrame',FsFrame, 'nSing',2)` returns `R.Field` `[nV×3]` (gauge-independent 3D Fiedler tangent field), `R.Magnitude`, `R.Phase`, `R.Singularities` (singularity vertex indices, `nSing` per component).
- `vFrame = nxr.manifold.measure.vertexFrame(mctx)` where `mctx = nxr.manifold.context(Vtx, double(Fcs))`; `vFrame` has `e1`, `e2`, `normals` `[nV×3]`.
- tess_tangents frame at vertices: `[Uf,~] = tess_tangents(SurfaceFile,'NoSave',1)` (per-face), then `[Uv,Vv] = bst_tangent_face2vertex(double(Faces), Uf, VertNormals)` (per-vertex).
- `Results = in_bst_results(ResultsFile, 1)` applies `Kernel*Data` and returns full `Results.ImageGridAmp` `[3·nV×nTime]`, `Results.nComponents==3`, `Results.Time` `[1×nTime]`. Unconstrained ordering is component-interleaved per vertex: `reshape(ImageGridAmp(:,ti), 3, []).'` gives `J` `[nV×3]` at time index `ti`.
- `SulciMap`: `TessMat.SulciMap` if present (binary `[nV×1]`, 1=sulcal), else `tess_sulcimap(TessMat)`.
- The real test cortex is the 20484-vertex cortex in the current protocol. Tests must locate it and SKIP if absent (idiom in `dev/tests/test_view_connection_phase.m::find_cortex_20484V`).
- The unconstrained kernel for Component 2 lives in the TutorialAuditory protocol, Subject01 deviant study: a results file with `Comment` starting `MN: MEG(Unconstr)` and `nComponents==3`.

**Branch:** `feature/sign-ambiguity-continuity` (already created; spec already committed).

**Toolbox:** Statistics and Machine Learning Toolbox is available — `rangesearch`/`knnsearch` may be used freely.

**Output directory produced at runtime:** `dev/benchmarks/sign_ambiguity/sign_ambiguity_run/` (figures + stats). Do not create it by hand; the driver creates it.

---

## Task 1: Shared frame construction — `sa_frames`

Builds the three vertex frames (Fiedler, tess_tangents, global-xyz) and the Fiedler decode once, for reuse by both components.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_frames.m`
- Test: `dev/tests/test_sa_frames.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_frames.m`:

```matlab
function test_sa_frames
% Smoke test: sa_frames returns three per-vertex tangent frames + the Fiedler
% decode on the real DB cortex. SKIP if no 20484V cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end

F = sa_frames(SurfaceFile);
nV = size(F.Vtx, 1);
for name = {'fiedler','tang','xyz'}
    fr = F.(name{1});
    assert(isequal(size(fr.e1), [nV 3]) && isequal(size(fr.e2), [nV 3]), ...
        'Frame %s must be [nV x 3].', name{1});
    % Where defined, e1 and e2 are unit and orthogonal.
    nrm = sqrt(sum(fr.e1.^2, 2));
    def = nrm > 1e-6;
    assert(all(abs(nrm(def) - 1) < 1e-3), 'Frame %s e1 not unit.', name{1});
    dt  = abs(sum(fr.e1(def,:) .* fr.e2(def,:), 2));
    assert(all(dt < 1e-3), 'Frame %s e1/e2 not orthogonal.', name{1});
end
assert(isequal(size(F.R.Field), [nV 3]), 'Fiedler field must be [nV x 3].');
assert(~isempty(F.ConnEig.ConnLaplacian), 'sa_frames must expose ConnLaplacian K.');
fprintf('PASSED: test_sa_frames (%d vertices).\n', nV);
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
fallback = '';
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices, 1) ~= 20484, continue; end
        if isempty(fallback), fallback = surf(iF).FileName; end
        SurfaceFile = surf(iF).FileName; return;
    end
end
if isempty(SurfaceFile), SurfaceFile = fallback; end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MATLAB MCP `run_matlab_test_file` on `dev/tests/test_sa_frames.m`.
Expected: FAIL — `Undefined function 'sa_frames'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_frames.m`:

```matlab
function F = sa_frames(SurfaceFile)
% SA_FRAMES: Build the three per-vertex tangent frames used by the sign-
% ambiguity continuity test, plus the Fiedler decode and the connection
% Laplacian operator.
%
% USAGE:  F = sa_frames(SurfaceFile)
%
% OUTPUT F (struct):
%   .Vtx, .Fcs, .Nv      : geometry ([nV x 3], [nF x 3] double, [nV x 3])
%   .VertConn            : sparse adjacency [nV x nV]
%   .ConnEig             : connection eigenmodes (incl. .ConnLaplacian = K)
%   .R                   : bst_conn_phase decode (Rank 1 Fiedler), with FsFrame
%   .fiedler             : struct e1,e2 ([nV x 3])  -- Fiedler frame (e1 = field dir)
%   .tang                : struct e1,e2 ([nV x 3])  -- tess_tangents (FS) frame
%   .xyz                 : struct e1,e2 ([nV x 3])  -- global-xyz projected frame
%   .vFrame              : nxr per-vertex frame (e1,e2,normals)
%
% Frames are zero-filled where undefined (e.g. Fiedler near singularities).

TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
Nv  = TessMat.VertNormals;
nV  = size(Vtx, 1);
if isfield(TessMat, 'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(Vtx, Fcs);
end

% Connection eigenmodes + nxr gauge frame.
ConnEig = bst_conn_eigenmodes_ensure(SurfaceFile);
mctx    = nxr.manifold.context(Vtx, Fcs);
vFrame  = nxr.manifold.measure.vertexFrame(mctx);

% tess_tangents (FreeSurfer trivial-connection) frame at vertices.
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
FsFrame  = struct('e1', Uv, 'e2', Vv);

% Fiedler decode (gauge-independent 3D field) + FS-gauge phase.
R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

% Fiedler frame: e1 = normalized field direction, e2 = n x e1 (zero where |field|~0).
fld = R.Field;
mag = sqrt(sum(fld.^2, 2));
def = mag > 1e-9;
fe1 = zeros(nV, 3);
fe1(def, :) = fld(def, :) ./ mag(def);
fe2 = cross(Nv, fe1, 2);
fiedler = struct('e1', fe1, 'e2', fe2);

% Global-xyz frame: world [1 0 0] projected to each tangent plane (fallback
% [0 1 0] where n is ~parallel to x), e2 = n x e1.
ref = repmat([1 0 0], nV, 1);
bad = abs(sum(Nv .* ref, 2)) > 0.95;        % n nearly parallel to x
ref(bad, :) = repmat([0 1 0], sum(bad), 1);
xe1 = ref - (sum(ref .* Nv, 2)) .* Nv;       % project to tangent plane
xn  = sqrt(sum(xe1.^2, 2));
xok = xn > 1e-9;
xe1(xok, :) = xe1(xok, :) ./ xn(xok);
xe2 = cross(Nv, xe1, 2);
xyz = struct('e1', xe1, 'e2', xe2);

F = struct('Vtx', Vtx, 'Fcs', Fcs, 'Nv', Nv, 'VertConn', VertConn, ...
           'ConnEig', ConnEig, 'R', R, 'fiedler', fiedler, ...
           'tang', FsFrame, 'xyz', xyz, 'vFrame', vFrame);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_frames.m` via MCP.
Expected: PASS `PASSED: test_sa_frames (20484 vertices).` (or `SKIP:` if no cortex — acceptable, note it).

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_frames.m dev/tests/test_sa_frames.m
git commit -m "feat(sign-ambiguity): sa_frames shared frame construction

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: Sulcal-wall pair detection — `sa_sulcal_walls`

Finds facing-wall vertex pairs restricted to sulcal vertices.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_sulcal_walls.m`
- Test: `dev/tests/test_sa_sulcal_walls.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_sulcal_walls.m`:

```matlab
function test_sa_sulcal_walls
% sa_sulcal_walls returns facing sulcal-wall pairs; every pair satisfies the
% sulcal / 3mm / anti-normal / non-adjacent criteria. SKIP if no cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end

TessMat = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices; Fcs = double(TessMat.Faces); Nv = TessMat.VertNormals;
if isfield(TessMat,'VertConn') && ~isempty(TessMat.VertConn)
    VertConn = TessMat.VertConn;
else
    VertConn = tess_vertconn(Vtx, Fcs);
end
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

opts = struct('MaxDist', 0.003, 'NormalDot', -0.7, 'Nring', 3);
pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts);
assert(~isempty(pairs), 'Expected a non-empty set of sulcal-wall pairs.');
assert(size(pairs,2) == 2, 'pairs must be [nPairs x 2].');

i = pairs(:,1); j = pairs(:,2);
% Both endpoints sulcal.
assert(all(SulciMap(i) > 0 & SulciMap(j) > 0), 'Both endpoints must be sulcal.');
% Within MaxDist.
d = sqrt(sum((Vtx(i,:) - Vtx(j,:)).^2, 2));
assert(all(d <= opts.MaxDist + 1e-9), 'Pairs must be within MaxDist.');
% Anti-aligned normals.
nd = sum(Nv(i,:) .* Nv(j,:), 2);
assert(all(nd < opts.NormalDot + 1e-9), 'Pairs must have anti-aligned normals.');
% Not mesh-adjacent within Nring.
A = double(VertConn > 0);
Aring = A; P = A;
for r = 2:opts.Nring, P = double((P*A) > 0); Aring = double((Aring + P) > 0); end
lin = sub2ind(size(Aring), i, j);
assert(all(Aring(lin) == 0), 'Pairs must not be within the N-ring.');
fprintf('PASSED: test_sa_sulcal_walls (%d pairs).\n', size(pairs,1));
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sa_sulcal_walls.m`.
Expected: FAIL — `Undefined function 'sa_sulcal_walls'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_sulcal_walls.m`:

```matlab
function pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts)
% SA_SULCAL_WALLS: Detect facing sulcal-wall vertex pairs.
%
% USAGE:  pairs = sa_sulcal_walls(Vtx, Nv, SulciMap, VertConn, opts)
%
% A pair (i,j) qualifies when:
%   (1) both i and j are sulcal vertices (SulciMap>0),
%   (2) ||Vtx_i - Vtx_j|| <= opts.MaxDist,
%   (3) Nv_i . Nv_j < opts.NormalDot  (anti-aligned -> facing walls),
%   (4) i and j are NOT within each other's opts.Nring mesh neighborhood
%       (so they are across a fold, not neighbors on the same wall).
%
% opts fields: MaxDist (m, e.g. 0.003), NormalDot (e.g. -0.7), Nring (e.g. 3).
% OUTPUT: pairs [nPairs x 2], unique unordered (i<j), sorted.

if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts,'MaxDist'),   opts.MaxDist   = 0.003; end
if ~isfield(opts,'NormalDot'), opts.NormalDot = -0.7;  end
if ~isfield(opts,'Nring'),     opts.Nring     = 3;     end

% (1) Restrict to sulcal vertices.
iSulci = find(SulciMap > 0);
if isempty(iSulci), pairs = zeros(0,2); return; end
Vs = Vtx(iSulci, :);

% (2) Range search among sulcal vertices (3D proximity).
nb = rangesearch(Vs, Vs, opts.MaxDist);

% N-ring adjacency (global indices) for criterion (4).
A = double(VertConn > 0);
Aring = A; P = A;
for r = 2:opts.Nring
    P = double((P * A) > 0);
    Aring = double((Aring + P) > 0);
end

I = []; J = [];
for a = 1:numel(iSulci)
    gi = iSulci(a);
    for b = nb{a}(:)'
        gj = iSulci(b);
        if gj <= gi, continue; end                 % unordered, i<j, drop self
        if sum(Nv(gi,:) .* Nv(gj,:), 2) >= opts.NormalDot, continue; end   % (3)
        if Aring(gi, gj) ~= 0, continue; end        % (4) within N-ring
        I(end+1,1) = gi; J(end+1,1) = gj; %#ok<AGROW>
    end
end
pairs = sortrows([I J]);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_sulcal_walls.m`.
Expected: PASS `PASSED: test_sa_sulcal_walls (N pairs).` with N > 0.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_sulcal_walls.m dev/tests/test_sa_sulcal_walls.m
git commit -m "feat(sign-ambiguity): sa_sulcal_walls facing-wall pair detection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Component 1 metrics — `sa_smoothness`

Data-free: per-edge normal variation `δn` and Fiedler covariant variation `δf`, partitioned by SulciMap, plus singularity-energy concentration.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_smoothness.m`
- Test: `dev/tests/test_sa_smoothness.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_smoothness.m`:

```matlab
function test_sa_smoothness
% Component 1 decoupling claim: on real cortex, normal angular variation is
% ELEVATED on sulcal edges, while the Fiedler covariant variation is NOT
% (statistically flat across the sulcal/crown partition). SKIP if no cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

S = sa_smoothness(F, SulciMap);

% Normal variation must be higher on sulcal edges than crown edges.
assert(S.dn_sulci_median > S.dn_crown_median, ...
    'Normal variation should be elevated at sulci (sulci=%.3f, crown=%.3f).', ...
    S.dn_sulci_median, S.dn_crown_median);
% Fiedler covariant variation must NOT be elevated at sulci: the sulcal median
% should not exceed the crown median by more than the normal field does. Use a
% loose ratio test (decoupling): df elevation ratio << dn elevation ratio.
dnRatio = S.dn_sulci_median / max(S.dn_crown_median, eps);
dfRatio = S.df_sulci_median / max(S.df_crown_median, eps);
assert(dfRatio < dnRatio, ...
    'Fiedler variation should be less folding-coupled than the normal (dfRatio=%.2f, dnRatio=%.2f).', ...
    dfRatio, dnRatio);
% Singularity-energy fraction is a valid fraction.
assert(S.singEnergyFrac >= 0 && S.singEnergyFrac <= 1, 'singEnergyFrac out of range.');
fprintf('PASSED: test_sa_smoothness (dnRatio=%.2f, dfRatio=%.2f, singFrac=%.2f).\n', ...
    dnRatio, dfRatio, S.singEnergyFrac);
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sa_smoothness.m`.
Expected: FAIL — `Undefined function 'sa_smoothness'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_smoothness.m`:

```matlab
function S = sa_smoothness(F, SulciMap)
% SA_SMOOTHNESS: Component 1 (data-free) edge metrics for the sign-ambiguity
% continuity test, partitioned by SulciMap.
%
% USAGE:  S = sa_smoothness(F, SulciMap)
%   F        : struct from sa_frames (uses .VertConn, .Nv, .ConnEig, .R).
%   SulciMap : [nV x 1] binary (1 = sulcal vertex).
%
% Per mesh edge (i,j), upper triangle:
%   dn(i,j) = acos(clamp(n_i . n_j))                      -- normal angular variation
%   df(i,j) = |wrap( arg(zF_i) - arg(zF_j) - arg(R_ij) )| -- Fiedler covariant variation
%             with R_ij = -K(i,j)/|K(i,j)|, zF the per-component Fiedler eigenvector.
% An edge is "sulcal" if either endpoint is sulcal.
%
% OUTPUT S: dn/df medians per partition, the decoupling ratios, and the fraction
% of total df^2 energy on edges touching a singularity neighborhood.

VertConn = F.VertConn;
Nv = F.Nv;
K  = F.ConnEig.ConnLaplacian;
nV = size(Nv, 1);

% Build the combined Fiedler eigenvector zF (per-component Fiedler columns have
% disjoint support, so summing them yields one [nV x 1] complex field).
ConnEig = F.ConnEig;
comps = unique(ConnEig.Component(:))';
zF = zeros(nV, 1);
for c = comps
    col = find(ConnEig.Component(:) == c & ConnEig.CompRank(:) == 1, 1);
    if isempty(col), continue; end
    zc = ConnEig.Vectors(:, col);
    zF(zc ~= 0) = zc(zc ~= 0);
end

% Edge list (upper triangle of the mesh adjacency).
[ii, jj] = find(triu(VertConn > 0, 1));

% dn: normal angular variation.
nd = sum(Nv(ii,:) .* Nv(jj,:), 2);
nd = max(min(nd, 1), -1);
dn = acos(nd);

% df: covariant phase increment using the connection's own transport arg(R_ij).
lin = sub2ind(size(K), ii, jj);
Kij = full(K(lin));
Rij = -Kij ./ max(abs(Kij), eps);            % transport rotation (unit complex)
inc = angle(zF(ii)) - angle(zF(jj)) - angle(Rij);
df  = abs(atan2(sin(inc), cos(inc)));        % wrap to [0, pi], magnitude
% Undefined where either endpoint has no Fiedler support.
defined = (zF(ii) ~= 0) & (zF(jj) ~= 0);

% Partition edges by SulciMap (sulcal if either endpoint is sulcal).
edgeSulcal = (SulciMap(ii) > 0) | (SulciMap(jj) > 0);

S = struct();
S.dn_sulci_median = median(dn(edgeSulcal));
S.dn_crown_median = median(dn(~edgeSulcal));
S.df_sulci_median = median(df(defined & edgeSulcal));
S.df_crown_median = median(df(defined & ~edgeSulcal));
S.dn_dfRatio_sulci = S.dn_sulci_median / max(S.dn_crown_median, eps);

% Singularity-energy concentration: fraction of total df^2 on edges within the
% 5-ring of any singularity vertex.
sing = F.R.Singularities(:);
near = false(nV, 1);
if ~isempty(sing)
    A = double(VertConn > 0);
    reach = sparse(sing, 1, 1, nV, 1);
    for r = 1:5, reach = double((A * reach + reach) > 0); end
    near = reach > 0;
end
edgeNearSing = near(ii) | near(jj);
df2 = df(defined).^2;
sel = edgeNearSing(defined);
S.singEnergyFrac = sum(df2(sel)) / max(sum(df2), eps);

% Carry the raw edge data for figures.
S.edges = [ii jj];
S.dn = dn; S.df = df; S.defined = defined; S.edgeSulcal = edgeSulcal;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_smoothness.m`.
Expected: PASS, printing the ratios and singularity fraction (e.g. `dnRatio` clearly > `dfRatio`).

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_smoothness.m dev/tests/test_sa_smoothness.m
git commit -m "feat(sign-ambiguity): sa_smoothness Component 1 covariant metrics

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Component 2 metrics — `sa_continuity`

Applies the unconstrained kernel, picks the M100 peak, computes per-frame phase, conditions on physical continuity, and reports per-pair continuity metrics.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_continuity.m`
- Test: `dev/tests/test_sa_continuity.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_continuity.m`:

```matlab
function test_sa_continuity
% Component 2: end-to-end on the unconstrained kernel. Conditioned pair set is
% non-empty; constrained sign-flip rate exceeds the Fiedler phase-discontinuity
% rate; Fiedler phase discontinuity <= tess_tangents <= global-xyz (medians).
% SKIP if no cortex or no unconstrained kernel.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end
ResultsFile = find_unconstrained_kernel(SurfaceFile);
if isempty(ResultsFile)
    fprintf('SKIP: no unconstrained (nComponents=3) kernel on this surface.\n'); return;
end

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3, ...
              'Window',[0.05 0.15], 'MagRatio',0.5);
C = sa_continuity(F, SulciMap, ResultsFile, opts);

assert(C.nPairs > 0, 'Conditioned pair set must be non-empty.');
assert(C.signFlipRate > C.phaseDiscRate.fiedler, ...
    'Constrained sign-flip rate (%.2f) should exceed Fiedler phase-disc rate (%.2f).', ...
    C.signFlipRate, C.phaseDiscRate.fiedler);
assert(C.phaseDisc.fiedler_median <= C.phaseDisc.tang_median + 1e-6, ...
    'Fiedler phase discontinuity should be <= tess_tangents.');
assert(C.phaseDisc.tang_median <= C.phaseDisc.xyz_median + 1e-6, ...
    'tess_tangents phase discontinuity should be <= global-xyz.');
fprintf(['PASSED: test_sa_continuity (nPairs=%d, signFlip=%.2f, ', ...
    'phaseDisc fiedler/tang/xyz = %.2f/%.2f/%.2f).\n'], C.nPairs, C.signFlipRate, ...
    C.phaseDisc.fiedler_median, C.phaseDisc.tang_median, C.phaseDisc.xyz_median);
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end

function ResultsFile = find_unconstrained_kernel(SurfaceFile)
% Search all studies for a results file on this surface with nComponents==3.
ResultsFile = '';
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies), return; end
allStudy = [sStudies.Study];
for iS = 1:numel(allStudy)
    res = allStudy(iS).Result;
    for iR = 1:numel(res)
        if isempty(res(iR).FileName), continue; end
        try
            M = load(file_fullpath(res(iR).FileName), 'nComponents', 'SurfaceFile');
        catch, continue; end
        if isfield(M,'nComponents') && isequal(M.nComponents, 3)
            ResultsFile = res(iR).FileName; return;
        end
    end
end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sa_continuity.m`.
Expected: FAIL — `Undefined function 'sa_continuity'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_continuity.m`:

```matlab
function C = sa_continuity(F, SulciMap, ResultsFile, opts)
% SA_CONTINUITY: Component 2 of the sign-ambiguity continuity test. Applies an
% unconstrained MEG kernel, picks the peak sample, and compares phase
% continuity across facing sulcal-wall pairs under three frames.
%
% USAGE:  C = sa_continuity(F, SulciMap, ResultsFile, opts)
%   F           : struct from sa_frames.
%   SulciMap    : [nV x 1] binary.
%   ResultsFile : unconstrained (nComponents=3) results/kernel file.
%   opts        : MaxDist, NormalDot, Nring (pair detection); Window [t0 t1] s
%                 (peak search); MagRatio (continuity conditioning).
%
% OUTPUT C: J (peak current [nV x 3]), tPeak, pairs, conditioned mask, per-frame
% phase discontinuities (circular distance in [0,pi]) and the constrained
% sign-flip rate.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts,'MaxDist'),  opts.MaxDist  = 0.003; end
if ~isfield(opts,'NormalDot'),opts.NormalDot= -0.7;  end
if ~isfield(opts,'Nring'),    opts.Nring    = 3;     end
if ~isfield(opts,'Window'),   opts.Window   = [0.05 0.15]; end
if ~isfield(opts,'MagRatio'), opts.MagRatio = 0.5;   end

% ===== Apply kernel and pick the peak sample =====
Results = in_bst_results(ResultsFile, 1);   % LoadFull=1 -> Kernel*Data
assert(isequal(Results.nComponents, 3), 'sa_continuity requires nComponents=3.');
IGA  = Results.ImageGridAmp;                % [3nV x nTime]
Time = Results.Time(:)';
nV   = size(F.Nv, 1);
inWin = Time >= opts.Window(1) & Time <= opts.Window(2);
if ~any(inWin), inWin = true(size(Time)); end
gfp = sqrt(sum(IGA.^2, 1));                 % global field power over sources
gfp(~inWin) = -inf;
[~, ti] = max(gfp);
J = reshape(IGA(:, ti), 3, []).';           % [nV x 3] component-interleaved
C.tPeak = Time(ti);

% ===== Sulcal-wall pairs =====
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);
i = pairs(:,1); j = pairs(:,2);

% ===== Conditioning on physical continuity =====
Ji = J(i,:); Jj = J(j,:);
mi = sqrt(sum(Ji.^2,2)); mj = sqrt(sum(Jj.^2,2));
cosang = sum(Ji.*Jj,2) ./ max(mi.*mj, eps);
magOk  = min(mi,mj) ./ max(max(mi,mj), eps) > opts.MagRatio;
keep   = (cosang > 0) & magOk;
i = i(keep); j = j(keep);
C.pairs = [i j];
C.nPairs = numel(i);
if C.nPairs == 0
    C.signFlipRate = NaN;
    C.phaseDisc = struct('fiedler_median',NaN,'tang_median',NaN,'xyz_median',NaN);
    C.phaseDiscRate = struct('fiedler',NaN);
    C.J = J; return;
end

% ===== Constrained sign-flip rate =====
s = sum(J .* F.Nv, 2);                       % constrained scalar J . n
C.signFlipRate = mean(sign(s(i)) ~= sign(s(j)));

% ===== Per-frame phase discontinuity (circular distance) =====
phaseOf = @(fr) atan2(sum(J .* fr.e2, 2), sum(J .* fr.e1, 2));
circDist = @(a,b) abs(atan2(sin(a-b), cos(a-b)));
pf = phaseOf(F.fiedler); pt = phaseOf(F.tang); px = phaseOf(F.xyz);
df = circDist(pf(i), pf(j));
dt = circDist(pt(i), pt(j));
dx = circDist(px(i), px(j));
C.phaseDisc = struct('fiedler_median', median(df), ...
                     'tang_median',    median(dt), ...
                     'xyz_median',     median(dx));
% A pair is "phase-discontinuous" if the circular distance exceeds pi/2.
C.phaseDiscRate = struct('fiedler', mean(df > pi/2), ...
                         'tang',    mean(dt > pi/2), ...
                         'xyz',     mean(dx > pi/2));
% Carry per-pair vectors for figures.
C.df = df; C.dt = dt; C.dx = dx; C.J = J;
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_continuity.m`.
Expected: PASS (or `SKIP:` if no kernel). On the auditory data the constrained sign-flip rate should be substantially > 0 and the Fiedler median phase discontinuity the smallest of the three.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_continuity.m dev/tests/test_sa_continuity.m
git commit -m "feat(sign-ambiguity): sa_continuity Component 2 per-frame phase metrics

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: Sulcal-crossing profiles — `sa_crossing_profile`

Component 1 illustrative companion: surface shortest path between facing-wall pairs (descends through the fundus), with cumulative normal vs. Fiedler-frame turning along arc length.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_crossing_profile.m`
- Test: `dev/tests/test_sa_crossing_profile.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_crossing_profile.m`:

```matlab
function test_sa_crossing_profile
% sa_crossing_profile builds crossing paths between facing-wall pairs; along the
% path the normal turns substantially (~pi, descends fundus + back up) while the
% Fiedler frame turns less (smoother). SKIP if no cortex.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3);
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);

P = sa_crossing_profile(F, pairs, 3);
assert(~isempty(P), 'Expected at least one crossing profile.');
for k = 1:numel(P)
    assert(numel(P(k).path) >= 3, 'Path too short.');
    assert(isequal(size(P(k).arclen), size(P(k).nAngle)), 'arclen/nAngle size mismatch.');
    assert(P(k).nAngle(end) >= P(k).fAngle(end) - 1e-6, ...
        'Normal should turn at least as much as the Fiedler frame across a fold.');
end
fprintf('PASSED: test_sa_crossing_profile (%d profiles, max normal turn %.2f rad).\n', ...
    numel(P), max(arrayfun(@(x) x.nAngle(end), P)));
end

function SurfaceFile = find_cortex_20484V()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sa_crossing_profile.m`.
Expected: FAIL — `Undefined function 'sa_crossing_profile'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_crossing_profile.m`:

```matlab
function P = sa_crossing_profile(F, pairs, nProfiles, J)
% SA_CROSSING_PROFILE: Build sulcal-crossing profiles for Component 1's
% illustrative figure. For up to nProfiles facing-wall pairs, the surface mesh
% shortest path between the two walls descends through the fundus and back up;
% along it the normal turns ~pi while the Fiedler frame turns smoothly.
%
% USAGE:  P = sa_crossing_profile(F, pairs, nProfiles)
%         P = sa_crossing_profile(F, pairs, nProfiles, J)   % adds constrained scalar
%
%   F        : struct from sa_frames (uses Vtx, Nv, VertConn, fiedler).
%   pairs    : [nPairs x 2] facing-wall pairs (from sa_sulcal_walls).
%   nProfiles: max number of profiles to build.
%   J        : optional [nV x 3] current; if given, P(k).s is J.n along the path.
%
% OUTPUT P: struct array with fields path (vertex idx [L x 1]), arclen [L x 1],
% nAngle [L x 1] (cumulative normal turning, rad), fAngle [L x 1] (cumulative
% Fiedler-frame turning, rad), and s [L x 1] (constrained scalar, or []).

if nargin < 3 || isempty(nProfiles), nProfiles = 3; end
if nargin < 4, J = []; end
P = struct('path',{},'arclen',{},'nAngle',{},'fAngle',{},'s',{});
if isempty(pairs), return; end

Vtx = F.Vtx; Nv = F.Nv; fe1 = F.fiedler.e1;
% Weighted surface graph (Euclidean edge lengths).
[ei, ej] = find(triu(F.VertConn > 0, 1));
w = sqrt(sum((Vtx(ei,:) - Vtx(ej,:)).^2, 2));
G = graph(ei, ej, w, size(Vtx,1));

k = min(nProfiles, size(pairs,1));
for p = 1:k
    pth = shortestpath(G, pairs(p,1), pairs(p,2));
    if numel(pth) < 3, continue; end
    pth = pth(:);
    seg = sqrt(sum(diff(Vtx(pth,:),1,1).^2, 2));
    arclen = [0; cumsum(seg)];
    nAng = local_cum_turn(Nv(pth,:));
    fAng = local_cum_turn(fe1(pth,:));
    rec = struct('path', pth, 'arclen', arclen, 'nAngle', nAng, 'fAngle', fAng, 's', []);
    if ~isempty(J), rec.s = sum(J(pth,:) .* Nv(pth,:), 2); end
    P(end+1) = rec; %#ok<AGROW>
end
end

function ang = local_cum_turn(V)
% Cumulative absolute turning angle of a sequence of vectors (rad).
L = size(V,1);
ang = zeros(L,1);
for t = 2:L
    a = V(t-1,:); b = V(t,:);
    na = norm(a); nb = norm(b);
    if na < 1e-9 || nb < 1e-9, ang(t) = ang(t-1); continue; end
    c = max(min(dot(a,b)/(na*nb), 1), -1);
    ang(t) = ang(t-1) + acos(c);
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_crossing_profile.m`.
Expected: PASS, with the max normal turn near `pi` and `nAngle(end) >= fAngle(end)` for each profile.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_crossing_profile.m dev/tests/test_sa_crossing_profile.m
git commit -m "feat(sign-ambiguity): sa_crossing_profile Component 1 illustrative paths

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Figures — `sa_figures`

Renders the Component 1 partition figure, the Component 1 crossing profiles, and the Component 2 histogram + summary. No deep unit test (visual artifact); a minimal smoke test asserts files are written.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sa_figures.m`
- Test: `dev/tests/test_sa_figures.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sa_figures.m`:

```matlab
function test_sa_figures
% sa_figures writes the expected PNG files from precomputed S (Component 1),
% P (crossing profiles), and optional C (Component 2). Uses tiny synthetic
% structs so it runs headless and fast (no DB needed).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));

outDir = fullfile(tempdir, 'sa_fig_test');
if exist(outDir, 'dir'), rmdir(outDir, 's'); end

% Minimal synthetic Component-1 stats.
S = struct('dn_sulci_median',0.8,'dn_crown_median',0.3, ...
           'df_sulci_median',0.12,'df_crown_median',0.10,'singEnergyFrac',0.6, ...
           'dn',rand(100,1),'df',rand(100,1),'defined',true(100,1), ...
           'edgeSulcal',[true(50,1);false(50,1)]);
% Minimal synthetic crossing profiles.
al = linspace(0,0.01,12)';
P = struct('path',{(1:12)'}, 'arclen',{al}, ...
           'nAngle',{linspace(0,pi,12)'}, 'fAngle',{linspace(0,0.4,12)'}, ...
           's',{[ones(6,1); -ones(6,1)]});
% Minimal synthetic Component-2 stats.
C = struct('nPairs',40,'signFlipRate',0.7,'tPeak',0.1, ...
           'df',rand(40,1)*0.5,'dt',rand(40,1),'dx',rand(40,1)*2, ...
           'phaseDisc',struct('fiedler_median',0.2,'tang_median',0.6,'xyz_median',1.4));

files = sa_figures(S, P, C, outDir);
assert(~isempty(files), 'sa_figures returned no files.');
for k = 1:numel(files)
    assert(exist(files{k}, 'file') == 2, 'Expected figure file missing: %s', files{k});
end
rmdir(outDir, 's');
fprintf('PASSED: test_sa_figures (%d files).\n', numel(files));
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sa_figures.m`.
Expected: FAIL — `Undefined function 'sa_figures'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sa_figures.m`:

```matlab
function files = sa_figures(S, P, C, outDir)
% SA_FIGURES: Render the sign-ambiguity continuity figures.
%
% USAGE:  files = sa_figures(S, P, C, outDir)
%   S      : Component 1 stats (from sa_smoothness); [] to skip.
%   P      : crossing profiles (from sa_crossing_profile); [] to skip.
%   C      : Component 2 stats (from sa_continuity);  [] to skip.
%   outDir : output directory (created if absent).
% OUTPUT files: cell array of written PNG paths.

if ~exist(outDir, 'dir'), mkdir(outDir); end
files = {};

% ===== Component 1: decoupling bar chart + df/dn distributions =====
if ~isempty(S)
    h1 = figure('Visible','off','Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    bar([S.dn_sulci_median S.dn_crown_median; S.df_sulci_median S.df_crown_median]);
    set(gca,'XTickLabel',{'normal \delta n','Fiedler \delta f'});
    legend({'sulcal edges','crown edges'},'Location','northeast');
    ylabel('median angular variation (rad)');
    title(sprintf('Folding coupling (sing. energy frac = %.2f)', S.singEnergyFrac));
    subplot(1,2,2);
    edges = linspace(0, pi, 30);
    histogram(S.df(S.defined & S.edgeSulcal), edges, 'Normalization','probability'); hold on;
    histogram(S.df(S.defined & ~S.edgeSulcal), edges, 'Normalization','probability');
    legend({'sulcal','crown'}); xlabel('\delta f (rad)'); ylabel('prob');
    title('Fiedler covariant variation by partition');
    f1 = fullfile(outDir, 'c1_smoothness.png');
    saveas(h1, f1); close(h1); files{end+1} = f1;
end

% ===== Component 1: sulcal-crossing profiles =====
if ~isempty(P)
    nP = numel(P);
    hP = figure('Visible','off','Color','w','Position',[100 100 320*nP 320]);
    for k = 1:nP
        subplot(1, nP, k);
        x = 1000 * P(k).arclen;                 % mm
        plot(x, P(k).nAngle, '-o', 'DisplayName','normal turn'); hold on;
        plot(x, P(k).fAngle, '-s', 'DisplayName','Fiedler turn');
        if ~isempty(P(k).s)
            sn = sign(P(k).s);                       % constrained sign trace
            plot(x, (pi/2)*(1+sn), ':', 'DisplayName','constrained sign');
        end
        xlabel('arc length (mm)'); ylabel('cumulative turn (rad)');
        title(sprintf('crossing %d', k));
        if k == 1, legend('Location','northwest'); end
    end
    fP = fullfile(outDir, 'c1_profiles.png');
    saveas(hP, fP); close(hP); files{end+1} = fP;
end

% ===== Component 2: per-frame phase-discontinuity histogram + summary =====
if ~isempty(C)
    h2 = figure('Visible','off','Color','w','Position',[100 100 900 380]);
    subplot(1,2,1);
    edges = linspace(0, pi, 30);
    histogram(C.df, edges, 'Normalization','probability'); hold on;
    histogram(C.dt, edges, 'Normalization','probability');
    histogram(C.dx, edges, 'Normalization','probability');
    legend({'Fiedler','tess\_tangents','global-xyz'});
    xlabel('phase discontinuity across wall (rad)'); ylabel('prob');
    title(sprintf('Phase continuity, %d pairs, t=%.0f ms', C.nPairs, 1000*C.tPeak));
    subplot(1,2,2);
    bar([C.signFlipRate; C.phaseDisc.fiedler_median; C.phaseDisc.tang_median; C.phaseDisc.xyz_median]);
    set(gca,'XTickLabel',{'constrained sign-flip','Fiedler \theta','tang \theta','xyz \theta'});
    xtickangle(20); ylabel('rate / median rad');
    title('Constrained artifact vs. frame phase discontinuity');
    f2 = fullfile(outDir, 'c2_continuity.png');
    saveas(h2, f2); close(h2); files{end+1} = f2;
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sa_figures.m`.
Expected: PASS `PASSED: test_sa_figures (3 files).`

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sa_figures.m dev/tests/test_sa_figures.m
git commit -m "feat(sign-ambiguity): sa_figures rendering (partition, profiles, continuity)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 7: Driver — `sign_ambiguity_run`

Ties everything together: resolves surface + kernel, runs Component 1 always, Component 2 if a kernel is found, writes figures + `stats.md`/`stats.csv`, SKIPs cleanly.

**Files:**
- Create: `dev/benchmarks/sign_ambiguity/sign_ambiguity_run.m`
- Test: `dev/tests/test_sign_ambiguity_run.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_sign_ambiguity_run.m`:

```matlab
function test_sign_ambiguity_run
% Driver smoke test: runs end to end on the real DB, writes stats + figures.
% SKIP if no 20484V cortex (Component 1 needs at least the surface).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');

out = sign_ambiguity_run();
if isempty(out)
    fprintf('SKIP: driver reported no usable cortex.\n'); return;
end
assert(exist(out.statsFile, 'file') == 2, 'stats.md not written.');
assert(~isempty(out.figures), 'No figures written.');
assert(isfield(out, 'C1') && ~isempty(out.C1), 'Component 1 stats missing.');
fprintf('PASSED: test_sign_ambiguity_run (out=%s).\n', out.outDir);
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_sign_ambiguity_run.m`.
Expected: FAIL — `Undefined function 'sign_ambiguity_run'`.

- [ ] **Step 3: Write minimal implementation**

Create `dev/benchmarks/sign_ambiguity/sign_ambiguity_run.m`:

```matlab
function out = sign_ambiguity_run(SurfaceFile, ResultsFile)
% SIGN_AMBIGUITY_RUN: Driver for the connection-eigenmode sulcal-wall continuity
% test. Runs Component 1 (data-free smoothness) on a cortex surface, and
% Component 2 (real-data phase continuity) if an unconstrained kernel is found.
% Writes figures + stats under dev/benchmarks/sign_ambiguity/sign_ambiguity_run/.
%
% USAGE:  out = sign_ambiguity_run()                       % auto-resolve
%         out = sign_ambiguity_run(SurfaceFile)            % given surface
%         out = sign_ambiguity_run(SurfaceFile, ResultsFile)
%
% Returns [] if no usable cortex is found (SKIP). Otherwise a struct with
% outDir, statsFile, figures, C1 (Component 1 stats), C2 (Component 2 stats or []).

if nargin < 1 || isempty(SurfaceFile), SurfaceFile = local_find_cortex(); end
if isempty(SurfaceFile), out = []; return; end
if nargin < 2, ResultsFile = local_find_kernel(SurfaceFile); end

thisDir = fileparts(mfilename('fullpath'));
outDir  = fullfile(thisDir, 'sign_ambiguity_run');
if ~exist(outDir, 'dir'), mkdir(outDir); end

% Shared frames + SulciMap.
F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

% Component 1 (always).
C1 = sa_smoothness(F, SulciMap);

% Component 2 (if a kernel is available).
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3, ...
              'Window',[0.05 0.15], 'MagRatio',0.5);
C2 = [];
if ~isempty(ResultsFile)
    C2 = sa_continuity(F, SulciMap, ResultsFile, opts);
end

% Crossing profiles (Component 1 illustration; uses a few facing-wall pairs).
% Overlay the constrained scalar from the peak current when Component 2 ran.
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);
if ~isempty(C2)
    P = sa_crossing_profile(F, pairs, 3, C2.J);
else
    P = sa_crossing_profile(F, pairs, 3);
end

% Figures + stats.
figures = sa_figures(C1, P, C2, outDir);
statsFile = fullfile(outDir, 'stats.md');
local_write_stats(statsFile, SurfaceFile, ResultsFile, C1, C2);
local_write_csv(fullfile(outDir, 'stats.csv'), C1, C2);

out = struct('outDir', outDir, 'statsFile', statsFile, 'figures', {figures}, ...
             'C1', C1, 'C2', C2, 'P', {P}, 'SurfaceFile', SurfaceFile, 'ResultsFile', ResultsFile);
fprintf('sign_ambiguity_run: wrote %s\n', outDir);
end

% ---- local helpers ----
function SurfaceFile = local_find_cortex()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
end

function ResultsFile = local_find_kernel(SurfaceFile)
ResultsFile = '';
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies), return; end
allStudy = [sStudies.Study];
for iS = 1:numel(allStudy)
    res = allStudy(iS).Result;
    for iR = 1:numel(res)
        if isempty(res(iR).FileName), continue; end
        try, M = load(file_fullpath(res(iR).FileName), 'nComponents'); catch, continue; end
        if isfield(M,'nComponents') && isequal(M.nComponents, 3)
            ResultsFile = res(iR).FileName; return;
        end
    end
end
end

function local_write_stats(statsFile, SurfaceFile, ResultsFile, C1, C2)
fid = fopen(statsFile, 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, '# Sign-ambiguity continuity test\n\n');
fprintf(fid, '- Surface: `%s`\n', SurfaceFile);
fprintf(fid, '- Kernel: `%s`\n\n', iff(isempty(ResultsFile),'(none)',ResultsFile));
fprintf(fid, '## Component 1 (data-free smoothness)\n\n');
fprintf(fid, '| metric | sulcal | crown |\n|---|---|---|\n');
fprintf(fid, '| normal var dn (median rad) | %.4f | %.4f |\n', C1.dn_sulci_median, C1.dn_crown_median);
fprintf(fid, '| Fiedler var df (median rad) | %.4f | %.4f |\n', C1.df_sulci_median, C1.df_crown_median);
fprintf(fid, '\n- singularity energy fraction of df^2: %.3f\n', C1.singEnergyFrac);
fprintf(fid, '- decoupling: dn sulci/crown ratio = %.2f, df sulci/crown ratio = %.2f\n\n', ...
    C1.dn_sulci_median/max(C1.dn_crown_median,eps), C1.df_sulci_median/max(C1.df_crown_median,eps));
if ~isempty(C2)
    fprintf(fid, '## Component 2 (real-data phase continuity)\n\n');
    fprintf(fid, '- peak time: %.0f ms, conditioned pairs: %d\n', 1000*C2.tPeak, C2.nPairs);
    fprintf(fid, '- constrained sign-flip rate: %.3f\n\n', C2.signFlipRate);
    fprintf(fid, '| frame | median phase disc (rad) | disc rate (>pi/2) |\n|---|---|---|\n');
    fprintf(fid, '| Fiedler | %.3f | %.3f |\n', C2.phaseDisc.fiedler_median, C2.phaseDiscRate.fiedler);
    fprintf(fid, '| tess_tangents | %.3f | %.3f |\n', C2.phaseDisc.tang_median, C2.phaseDiscRate.tang);
    fprintf(fid, '| global-xyz | %.3f | %.3f |\n', C2.phaseDisc.xyz_median, C2.phaseDiscRate.xyz);
end
end

function local_write_csv(csvFile, C1, C2)
fid = fopen(csvFile, 'w');
c = onCleanup(@() fclose(fid));
fprintf(fid, 'component,metric,value\n');
fprintf(fid, 'C1,dn_sulci_median,%.6f\n', C1.dn_sulci_median);
fprintf(fid, 'C1,dn_crown_median,%.6f\n', C1.dn_crown_median);
fprintf(fid, 'C1,df_sulci_median,%.6f\n', C1.df_sulci_median);
fprintf(fid, 'C1,df_crown_median,%.6f\n', C1.df_crown_median);
fprintf(fid, 'C1,singEnergyFrac,%.6f\n', C1.singEnergyFrac);
if ~isempty(C2)
    fprintf(fid, 'C2,nPairs,%d\n', C2.nPairs);
    fprintf(fid, 'C2,signFlipRate,%.6f\n', C2.signFlipRate);
    fprintf(fid, 'C2,phaseDisc_fiedler,%.6f\n', C2.phaseDisc.fiedler_median);
    fprintf(fid, 'C2,phaseDisc_tang,%.6f\n', C2.phaseDisc.tang_median);
    fprintf(fid, 'C2,phaseDisc_xyz,%.6f\n', C2.phaseDisc.xyz_median);
end
end

function s = iff(cond, a, b)
if cond, s = a; else, s = b; end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_sign_ambiguity_run.m`.
Expected: PASS, writing `dev/benchmarks/sign_ambiguity/sign_ambiguity_run/{c1_smoothness.png, c2_continuity.png (if kernel), stats.md, stats.csv}`.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/sign_ambiguity/sign_ambiguity_run.m dev/tests/test_sign_ambiguity_run.m
git commit -m "feat(sign-ambiguity): sign_ambiguity_run driver + stats output

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Final verification

- [ ] Run all seven tests in sequence via the MATLAB MCP; confirm PASS or documented SKIP for each.
- [ ] Run `sign_ambiguity_run()` once interactively and visually inspect the two figures + `stats.md` to confirm the result statements in spec §5 hold (Fiedler covariant variation decoupled from SulciMap; constrained sign-flip rate high while Fiedler phase discontinuity lowest of the three frames).
- [ ] Dispatch a final code review over the whole `dev/benchmarks/sign_ambiguity/` directory.
- [ ] Use `superpowers:finishing-a-development-branch`.
```
