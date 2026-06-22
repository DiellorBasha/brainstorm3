# M3 Plan B — phase-readout math

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Turn a stored connection-Laplacian eigenmode into a usable geometric readout — a gauge-independent 3D tangent field, its magnitude, its singularities, and a between-subject (FreeSurfer-gauge) phase — plus the per-vertex FreeSurfer reference frame that the FS-gauge phase needs.

**Architecture:** Two pure-MATLAB functions on top of Plan A's nxr `vertexFrame` export and M2's `ConnEigenmodes`. `bst_tangent_face2vertex` transfers the per-face trivial-connection frame to vertices. `bst_conn_phase` decodes an eigenmode's complex coordinates into a 3D field via the nxr per-vertex frame (`w = Re(z)·e1 + Im(z)·e2`), then derives magnitude, singularities, and the FS-gauge phase (angle of the field in the per-vertex FS frame).

**Tech Stack:** MATLAB, Brainstorm (`in_tess_bst`, `tess_tangents`, `bst_conn_eigenmodes_ensure`, `in_tess_conn_eigenmodes`), the nxr `+nxr` wrappers (`measure.vertexFrame`, `interpolate.smoothVertex`).

**Reference spec:** `dev/connection_phase_readout_integration.md` (§3.2, §3.3).

**Depends on:** Plan A (the installed plugin now provides `nxr.manifold.measure.vertexFrame`); M2 (`ConnEigenmodes`).

**Scope note / refinement of the spec:** Plan B delivers the **gauge-independent Field** + the **FS-gauge (between-subject) Phase**. The *within-subject intrinsic* phase is, per spec §2/§9, best rendered **reference-free** in the viewer (Plan C) via stripes/iso-phase contours derived from `Field` — so it is NOT a committed scalar here. `bst_conn_phase` returns `Field`/`Magnitude`/`Singularities` always, and `Phase` (FS gauge) when a FS frame is supplied.

---

## File Structure

| Path | Responsibility |
|---|---|
| `toolbox/anatomy/bst_tangent_face2vertex.m` | Transfer a per-face tangent direction to a per-vertex orthonormal frame. |
| `toolbox/math/bst_conn_phase.m` | Decode a connection eigenmode → 3D field, magnitude, singularities, FS-gauge phase. |
| `dev/tests/test_tangent_face2vertex.m` | Test the face→vertex transfer (real cortex). |
| `dev/tests/test_conn_phase.m` | Gauge cross-check (decode vs `smoothVertex`) + `bst_conn_phase` outputs. |

**Test environment (both tests):** run via the MATLAB MCP against the live session (Brainstorm running, TutorialAuditory, real 20484-vertex cortex; nxr-compute loaded). Resolve the cortex by vertex count (same resolver idiom as prior tests); SKIP cleanly if absent. Never run bare `clear`.

---

## Task 1: `bst_tangent_face2vertex`

**Files:** Create `toolbox/anatomy/bst_tangent_face2vertex.m`; Test `dev/tests/test_tangent_face2vertex.m`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tangent_face2vertex.m`:

```matlab
function test_tangent_face2vertex
% Transfer the per-face trivial-connection frame (tess_tangents) to a per-vertex
% orthonormal frame and check it is unit, orthonormal, right-handed, and smooth.
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
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
TessMat = in_tess_bst(SurfaceFile);
F  = double(TessMat.Faces);
N  = TessMat.VertNormals;
nV = size(N, 1);

% Per-face trivial-connection frame (e1 = U).
[Uf, ~] = tess_tangents(SurfaceFile, 'NoSave', 1);

[Uv, Vv] = bst_tangent_face2vertex(F, Uf, N);

assert(isequal(size(Uv), [nV 3]) && isequal(size(Vv), [nV 3]), 'Uv/Vv must be nV x 3.');
assert(max(abs(sqrt(sum(Uv.^2,2)) - 1)) < 1e-6, 'Uv must be unit.');
assert(max(abs(sqrt(sum(Vv.^2,2)) - 1)) < 1e-6, 'Vv must be unit.');
assert(max(abs(sum(Uv .* Vv, 2))) < 1e-6, 'Uv . Vv must be ~0.');
Nu = N ./ max(sqrt(sum(N.^2,2)), eps);
assert(max(abs(sum(Uv .* Nu, 2))) < 1e-6, 'Uv must lie in the tangent plane (Uv . n ~0).');
cr = cross(Nu, Uv, 2);
assert(max(max(abs(cr - Vv))) < 1e-6, 'Vv must equal n x Uv (right-handed).');

% Smoothness: the per-vertex e1 should agree with most incident face e1's
% (the trivial-connection field is smooth away from the few singularities).
agree = abs(sum(Uv(F(:,1),:) .* Uf, 2));   % |cos angle| between vertex-1 frame and its face
frac = mean(agree > 0.9);
assert(frac > 0.9, 'Per-vertex frame should align with incident face frames for >90%% of faces (got %.2f).', frac);

fprintf('PASSED: per-vertex FS frame (nV=%d): orthonormal, right-handed, smooth (%.1f%% aligned).\n', nV, 100*frac);
fprintf('ALL TESTS PASSED: test_tangent_face2vertex\n');
end


function SurfaceFile = find_cortex_20484V()
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

Run `dev/tests/test_tangent_face2vertex.m` via the MATLAB MCP.
Expected: FAIL — `Unrecognized function or variable 'bst_tangent_face2vertex'`.

- [ ] **Step 3: Write the function**

Create `toolbox/anatomy/bst_tangent_face2vertex.m`:

```matlab
function [Uv, Vv] = bst_tangent_face2vertex(Faces, Uf, VertNormals)
% BST_TANGENT_FACE2VERTEX: Per-face tangent direction -> per-vertex orthonormal frame.
%
% USAGE:  [Uv, Vv] = bst_tangent_face2vertex(Faces, Uf, VertNormals)
%
% DESCRIPTION:
%     Transfers a per-FACE tangent direction Uf (e.g. the trivial-connection e1
%     from tess_tangents) to a per-VERTEX orthonormal tangent frame. Each vertex
%     averages the e1 vectors of its incident faces, projects the result into the
%     vertex tangent plane (removing the vertex-normal component) and renormalizes;
%     Vv = n x Uv completes a right-handed frame. The trivial-connection field is
%     smooth, so plain averaging is valid away from the (few) singularities, where
%     the incident directions cancel and Uv is left as an arbitrary in-plane unit
%     vector (those vertices are the field singularities).
%
% INPUTS:
%     Faces       : [nF x 3] 1-based triangle indices.
%     Uf          : [nF x 3] per-face tangent direction (unit).
%     VertNormals : [nV x 3] per-vertex normals.
%
% OUTPUTS:
%     Uv : [nV x 3] per-vertex first tangent (unit, in the tangent plane).
%     Vv : [nV x 3] per-vertex second tangent = n x Uv.
%
% SEE ALSO: tess_tangents, bst_conn_phase

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

Faces = double(Faces);
nV = size(VertNormals, 1);

% Accumulate incident-face e1 at each vertex (uniform average over incident faces).
ii   = [Faces(:,1); Faces(:,2); Faces(:,3)];   % [3nF x 1]
Urep = [Uf;          Uf;          Uf];          % [3nF x 3]
Usum = zeros(nV, 3);
for d = 1:3
    Usum(:, d) = accumarray(ii, Urep(:, d), [nV, 1]);
end

% Unit vertex normals.
N = VertNormals ./ max(sqrt(sum(VertNormals.^2, 2)), eps);

% Project the averaged direction into the tangent plane and renormalize.
Uv  = Usum - sum(Usum .* N, 2) .* N;
nrm = sqrt(sum(Uv.^2, 2));
% Guard singularities (cancellation): fall back to an arbitrary in-plane axis.
bad = nrm < 1e-9;
if any(bad)
    ref = repmat([1 0 0], sum(bad), 1);
    alt = abs(N(bad,1)) > 0.9;          % avoid degeneracy when n ~ x-axis
    ref(alt, :) = repmat([0 1 0], sum(alt), 1);
    proj = ref - sum(ref .* N(bad,:), 2) .* N(bad,:);
    Uv(bad, :) = proj;
    nrm(bad)   = sqrt(sum(proj.^2, 2));
end
Uv = Uv ./ max(nrm, eps);
Vv = cross(N, Uv, 2);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run `dev/tests/test_tangent_face2vertex.m` via the MATLAB MCP.
Expected: `ALL TESTS PASSED: test_tangent_face2vertex`. Paste the output. (A `SKIP:` is not acceptable.)

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/bst_tangent_face2vertex.m dev/tests/test_tangent_face2vertex.m
git commit -m "feat(conn-phase): bst_tangent_face2vertex (per-vertex FS reference frame)"
```

---

## Task 2: `bst_conn_phase`

**Files:** Create `toolbox/math/bst_conn_phase.m`; Test `dev/tests/test_conn_phase.m`.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_conn_phase.m`:

```matlab
function test_conn_phase
% (A) Gauge cross-check: decoding nxr's smoothVertex raw field via the exported
%     per-vertex frame must align with its world vectors -> the frame is the
%     correct gauge for vertex-domain complex coordinates (validates the decode
%     that bst_conn_phase relies on, incl. the connection eigenmodes which share
%     the same gauge).
% (B) bst_conn_phase outputs on the real cortex: gauge-independent 3D field
%     (tangent), magnitude = |z|, FS-gauge phase that winds, and singularities.
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
    fprintf('SKIP: no 20484-vertex cortex in the current protocol.\n');
    return;
end
% Work on a temp COPY so bst_conn_eigenmodes_ensure does not write ConnEigenmodes
% into the shared DB surface (named tess_cortex_*.mat so file_fullpath accepts it).
tmpFile = fullfile(tempdir, 'tess_cortex_connphase.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

TessMat = in_tess_bst(tmpFile);
V  = TessMat.Vertices;
F  = double(TessMat.Faces);
N  = TessMat.VertNormals;
nV = size(V, 1);

mctx   = nxr.manifold.context(V, F);
vFrame = nxr.manifold.measure.vertexFrame(mctx);

% ===== (A) Gauge cross-check via smoothVertex (nSym=1) =====
sm    = nxr.manifold.interpolate.smoothVertex(mctx, 1);   % nSym=1: true vector field
raw   = sm.vertexFieldRaw;     % [nV x 2] complex coords in the per-vertex frame
world = sm.vertexVectors;      % [nV x 3] world-space vectors
decoded = raw(:,1) .* vFrame.e1 + raw(:,2) .* vFrame.e2;   % [nV x 3]
dn = decoded ./ max(sqrt(sum(decoded.^2,2)), eps);
wn = world   ./ max(sqrt(sum(world.^2,2)),   eps);
dp = abs(sum(dn .* wn, 2));    % |cos angle|; allow global nRoSy sign
mask = (sqrt(sum(decoded.^2,2)) > 1e-6) & (sqrt(sum(world.^2,2)) > 1e-6);
assert(median(dp(mask)) > 0.99, ...
    'GAUGE: decoded smoothVertex raw field must align with its world vectors (median |cos|=%.4f).', median(dp(mask)));
fprintf('PASSED (A): vertex-frame gauge validated against smoothVertex (median |cos|=%.4f).\n', median(dp(mask)));

% ===== (B) bst_conn_phase =====
ConnEig = bst_conn_eigenmodes_ensure(tmpFile, 20);   % small, fast; explicit count skips scalar-axis ensure
[Uf, ~] = tess_tangents(tmpFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(F, Uf, N);
FsFrame = struct('e1', Uv, 'e2', Vv);

R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2);

assert(isequal(size(R.Field), [nV 3]), 'Field must be nV x 3.');
supp = find(any(R.Field ~= 0, 2));
assert(~isempty(supp), 'Field must be nonzero on the Fiedler support.');
% Field is tangent to the (nxr) vertex normal by construction.
tang = abs(sum(R.Field(supp,:) .* vFrame.normals(supp,:), 2));
assert(max(tang) < 1e-6, 'Field must be tangent (Field . n ~0).');
% Magnitude equals |z| (nonneg, positive on support).
assert(all(R.Magnitude(supp) > 0), 'Magnitude must be positive on support.');
% FS-gauge phase: finite on support, within [-pi, pi].
assert(all(isfinite(R.Phase(supp))), 'Phase must be finite on support.');
assert(max(abs(R.Phase(supp))) <= pi + 1e-9, 'Phase must lie in [-pi, pi].');
% Phase winds: spans more than half the circle on the support of one component.
comp1 = find(ConnEig.Component(:)==1 & ConnEig.CompRank(:)==1, 1);
idx1 = find(ConnEig.Vectors(:, comp1) ~= 0);
assert((max(R.Phase(idx1)) - min(R.Phase(idx1))) > pi, ...
    'FS-gauge phase should wind over a hemisphere (range > pi).');
% Singularities: nSing per component, with small magnitude (a field dip).
nComp = numel(unique(ConnEig.Component(:)));
assert(numel(R.Singularities) == 2 * nComp, 'Expected 2 singularities per component.');
assert(median(R.Magnitude(R.Singularities)) < median(R.Magnitude(supp)), ...
    'Singularity magnitude should dip below the component median.');

fprintf('PASSED (B): bst_conn_phase field/magnitude/phase/singularities on the cortex (nV=%d).\n', nV);
fprintf('ALL TESTS PASSED: test_conn_phase\n');
end


function SurfaceFile = find_cortex_20484V()
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

Run `dev/tests/test_conn_phase.m` via the MATLAB MCP.
Expected: FAIL — part (A) may pass (it uses only Plan A's `vertexFrame`), but the test fails at `bst_conn_phase` with `Unrecognized function or variable 'bst_conn_phase'`. (If part (A) itself fails, that is a real gauge problem — STOP and report it, do not proceed.)

- [ ] **Step 3: Write the function**

Create `toolbox/math/bst_conn_phase.m`:

```matlab
function R = bst_conn_phase(ConnEig, vFrame, varargin)
% BST_CONN_PHASE: Decode a connection-Laplacian eigenmode into a 3D field + phase.
%
% USAGE:  R = bst_conn_phase(ConnEig, vFrame)
%         R = bst_conn_phase(ConnEig, vFrame, 'Rank', 1, 'FsFrame', FsFrame, 'nSing', 2)
%
% DESCRIPTION:
%     Given a surface's connection eigenmodes (in_tess_conn_eigenmodes / M2) and
%     the nxr per-vertex tangent frame (nxr.manifold.measure.vertexFrame), decodes
%     the selected eigenmode's complex coordinates z into a gauge-independent 3D
%     tangent field  w = Re(z)*e1 + Im(z)*e2,  and derives its magnitude,
%     singularities, and the between-subject (FreeSurfer-gauge) phase. The eigen-
%     modes are block-structured per connected component (hemisphere); for each
%     component the column of within-component rank 'Rank' is used (Rank 1 = the
%     Fiedler / smoothest field). The intrinsic within-subject phase is NOT a
%     scalar here: it is rendered reference-free in the viewer (stripes / iso-phase
%     contours) from R.Field.
%
% INPUTS:
%     ConnEig : struct from in_tess_conn_eigenmodes (fields Vectors [nV x nModes]
%               complex, Component, CompRank).
%     vFrame  : struct with e1, e2, normals ([nV x 3]) from nxr vertexFrame.
%
% OPTIONS:
%     'Rank'    : within-component mode rank to use (default 1 = Fiedler).
%     'FsFrame' : struct with e1, e2 ([nV x 3]) per-vertex FreeSurfer frame; if
%                 given, R.Phase is the field angle in that frame. Default [] (NaN).
%     'nSing'   : singularities to report per component (default 2; Poincare-Hopf
%                 index sum on a hemisphere-sphere is 2).
%
% OUTPUTS:
%     R.Field        : [nV x 3] 3D tangent field (gauge-independent); 0 off-support.
%     R.Magnitude    : [nV x 1] |z|; 0 off-support.
%     R.Phase        : [nV x 1] FS-gauge phase in [-pi, pi]; NaN off-support / no FsFrame.
%     R.Singularities: [k x 1] vertex indices (the nSing smallest-|z| per component).
%     R.Rank         : the rank used.
%
% SEE ALSO: in_tess_conn_eigenmodes, bst_tangent_face2vertex, nxr.manifold.measure.vertexFrame

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
Rank    = 1;
FsFrame = [];
nSing   = 2;
for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'rank',    Rank = varargin{i+1};
        case 'fsframe', FsFrame = varargin{i+1};
        case 'nsing',   nSing = varargin{i+1};
    end
end

nV = size(vFrame.e1, 1);
if size(ConnEig.Vectors, 1) ~= nV
    error('bst_conn_phase:sizeMismatch', ...
        'ConnEig.Vectors (%d rows) must match vFrame (%d vertices).', size(ConnEig.Vectors,1), nV);
end

Field     = zeros(nV, 3);
Magnitude = zeros(nV, 1);
Phase     = nan(nV, 1);
Sing      = zeros(0, 1);

comps = unique(ConnEig.Component(:))';
for c = comps
    col = find(ConnEig.Component(:) == c & ConnEig.CompRank(:) == Rank, 1);
    if isempty(col)
        continue;
    end
    z   = ConnEig.Vectors(:, col);      % complex [nV x 1], nonzero only on component c
    idx = find(z ~= 0);                 % component support (block structure is exact-zero off it)
    zc  = z(idx);

    e1 = vFrame.e1(idx, :);
    e2 = vFrame.e2(idx, :);
    w  = real(zc) .* e1 + imag(zc) .* e2;   % 3D tangent field
    Field(idx, :)  = w;
    Magnitude(idx) = abs(zc);

    % Singularities: the nSing smallest-|z| vertices of this component.
    [~, ord] = sort(abs(zc), 'ascend');
    Sing = [Sing; idx(ord(1:min(nSing, numel(ord))))]; %#ok<AGROW>

    % Between-subject phase: the field angle in the per-vertex FreeSurfer frame.
    if ~isempty(FsFrame)
        U  = FsFrame.e1(idx, :);
        Vw = FsFrame.e2(idx, :);
        Phase(idx) = atan2(sum(w .* Vw, 2), sum(w .* U, 2));
    end
end

R = struct('Field', Field, 'Magnitude', Magnitude, 'Phase', Phase, ...
           'Singularities', Sing, 'Rank', Rank);
end
```

- [ ] **Step 4: Run the test to verify it PASSES**

Run `dev/tests/test_conn_phase.m` via the MATLAB MCP.
Expected: prints `PASSED (A)` (gauge), `PASSED (B)`, then `ALL TESTS PASSED: test_conn_phase`. Paste the output.

If part (B)'s phase-winding or singularity assertions fail on the real cortex (these probe genuine geometry, not just plumbing): do NOT weaken them blindly — STOP and report `DONE_WITH_CONCERNS` with the observed numbers (e.g. the actual phase range / singularity magnitudes), so the geometry can be reviewed. Do not change `bst_conn_phase`'s math to force a pass.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_conn_phase.m dev/tests/test_conn_phase.m
git commit -m "feat(conn-phase): bst_conn_phase (decode field, magnitude, singularities, FS-gauge phase)"
```

---

## Notes for the implementer

- **Gauge cross-check (Part A) is the key correctness gate** — it proves the nxr per-vertex frame correctly decodes vertex-domain complex coordinates, which is exactly what `bst_conn_phase` relies on (the connection eigenmodes share `smoothVertex`'s gauge). If it fails, the problem is upstream (Plan A), not in this code.
- **Intrinsic phase is intentionally not computed here** — it is a reference-free *rendering* (stripes / iso-phase contours from `R.Field`) done in Plan C. `bst_conn_phase` returns the gauge-independent `Field` that those renderings consume.
- **Singularity detection is a simple smallest-|z| heuristic** (sufficient for markers/validation); a robust discrete-index detector is a possible Plan C refinement.
- Never run bare `clear`; resolve the cortex by vertex count; SKIP cleanly if absent.
```
