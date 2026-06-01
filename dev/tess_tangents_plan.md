# tess_tangents Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `toolbox/anatomy/tess_tangents.m` — a globally consistent per-face tangent frame field computed per hemisphere via nxr's trivial connection, with singularities at the FreeSurfer registration-sphere poles, stored on the surface file.

**Architecture:** SurfaceFile-based function. Loads the cortex via `in_tess_bst`, requires `Reg.Sphere.Vertices`, guards on the nxr-compute plugin (SPM-style), splits into hemispheres with `tess_hemisplit`, and for each hemisphere extracts a local submesh, places two `+1` singularities at the max/min-z registration-sphere poles, and calls `nxr.manifold.interpolate.trivial`. The per-face direction/orthogonal vectors are scattered back to full-surface face indexing, oriented right-handed w.r.t. the face normal, and stored as `TessMat.TangentFrame` (struct, `Domain='face'`).

**Tech Stack:** MATLAB R2023b, Brainstorm (`in_tess_bst`, `tess_hemisplit`, `bst_plugin`, `bst_history`, `bst_save`, `file_fullpath`), nxr-compute MEX (`nxr.manifold.context`, `nxr.manifold.interpolate.trivial`). Tests run via the MATLAB MCP by direct function invocation (the repo idiom; `runtests` rejects function-style tests).

**Repo:** `/Users/diellorbasha/workspace/research/code/brainstorm3`, branch `feature/tess-tangents` (off `development`). Do NOT switch branches.

**Verified facts (do not re-derive):**
- `tess_hemisplit(sSurf)` → `[rH, lH, isConnected, iStruct, iRightScout, iLeftScout]` (right first, left second); uses the 'Structures' atlas, falls back to a y-coordinate split.
- `nxr.manifold.context(V, F)` takes 1-based faces (the mex subtracts 1 internally); returns `mctx`. `nxr.manifold.interpolate.trivial(mctx, singVerts, singValues)` takes 1-based `singVerts` and `singValues` that must sum to χ (=2 per hemisphere sphere); returns a struct with `directionVectors [nF×3]`, `orthogonalVectors [nF×3]`, `eulerCharacteristic`, `gaussBonnetSatisfied`.
- `TessMat.Reg.Sphere.Vertices [nV×3]` holds per-vertex registration-sphere coords; poles = max/min z. North/south per hemisphere are taken within that hemisphere's vertex set.
- `file_fullpath` returns an already-existing absolute path unchanged **iff** the filename matches a surface type — so temp test surfaces MUST be named `tess_cortex_*.mat` (same idiom as `dev/tests/test_io_eigenmodes_roundtrip.m`).
- `bst_history('add', FileMat, eventType, eventDesc)` returns the updated struct.
- The local DB has FreeSurfer-registered cortices (e.g. `sub-0002/tess_cortex_pial_low.mat`, 20484V/40960F, ico5 per hemisphere). The integration test discovers one at runtime and prefers a low-res one.
- nxr-compute is installable via `bst_plugin('Install','nxr-compute')` (release `plugin-v1.0.0` is live) and is already present in the user plugins dir from prior work.

---

## File Structure

- Create: `toolbox/anatomy/tess_tangents.m` — main function `tess_tangents` + local helper `solve_hemisphere`. One file, one responsibility (compute + store the tangent frame).
- Create: `dev/tests/test_tess_tangents_guard.m` — cheap guard test (errors on a surface without `Reg.Sphere`).
- Create: `dev/tests/test_tess_tangents.m` — integration test on a registered cortex (geometry, determinism, storage round-trip).

---

## Task 1: Implement `tess_tangents.m` + guard test

**Files:**
- Create: `toolbox/anatomy/tess_tangents.m`
- Create: `dev/tests/test_tess_tangents_guard.m`

- [ ] **Step 1: Write the guard test (failing)**

Create `dev/tests/test_tess_tangents_guard.m`:

```matlab
function test_tess_tangents_guard
% tess_tangents must error cleanly when the surface has no FreeSurfer
% registration sphere (Reg.Sphere). This path must NOT require nxr.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% Build a temp surface WITHOUT Reg.Sphere. Name it tess_cortex_*.mat so
% file_fullpath recognizes the type and returns the absolute path as-is.
[V, F] = tess_sphere(162);
tmpFile = fullfile(tempdir, 'tess_cortex_tangents_guard.mat');
bst_save(tmpFile, struct('Vertices', V, 'Faces', F, 'Comment', 'noreg'), 'v7');
cleanup = onCleanup(@() delete(tmpFile));

threw = false;
try
    tess_tangents(tmpFile, 'NoSave', 1);
catch ME
    threw = strcmp(ME.identifier, 'tess_tangents:noRegSphere');
    if ~threw
        rethrow(ME);
    end
end
assert(threw, 'tess_tangents did not raise tess_tangents:noRegSphere on an unregistered surface.');
fprintf('ALL TESTS PASSED: test_tess_tangents_guard\n');
end
```

- [ ] **Step 2: Run it to confirm it fails**

Load the MATLAB MCP eval tool via ToolSearch (`select:mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code`), then run:
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); addpath(genpath('dev/tests')); clear functions; test_tess_tangents_guard
```
Expected: FAIL — `tess_tangents` is undefined (`Undefined function or variable 'tess_tangents'`).

- [ ] **Step 3: Implement `tess_tangents.m`**

Create `toolbox/anatomy/tess_tangents.m`:

```matlab
function [U, V] = tess_tangents(SurfaceFile, varargin)
% TESS_TANGENTS: Globally consistent per-face tangent frame field.
%
% USAGE:  [U, V] = tess_tangents(SurfaceFile)              % compute, store, return
%         [U, V] = tess_tangents(SurfaceFile, 'NoSave', 1) % compute + return only
%
% DESCRIPTION:
%     Computes a globally consistent tangent frame field on a cortical surface
%     using nxr-compute's trivial connection. Two singularities are placed at the
%     north and south poles of each hemisphere's FreeSurfer registration sphere
%     (Reg.Sphere). Because the two hemispheres are disconnected genus-0 spheres,
%     the field is solved PER HEMISPHERE (each: Euler characteristic 2, two +1
%     singularities). The result is a per-FACE orthonormal tangent frame (U, V).
%
%     The frame is stored on the surface file as TessMat.TangentFrame (Domain
%     'face'). Transferring it to a per-vertex frame (via the Hodge star / DEC)
%     is future work; the storage format reserves a 'Domain' field for it.
%
% INPUT:
%     - SurfaceFile : Brainstorm cortex surface with a FreeSurfer registration
%                     sphere (TessMat.Reg.Sphere.Vertices).
% OPTIONS:
%     - NoSave : (logical) If true, do not write TangentFrame back to the file.
%                Default: false (store).
% OUTPUT:
%     - U : [nF x 3] per-face first tangent direction (e1).
%     - V : [nF x 3] per-face orthogonal tangent (e2); cross(U,V) is aligned to
%           the face normal (right-handed frame).
%
% Requires the nxr-compute plugin (no MATLAB fallback for the trivial connection).
%
% SEE ALSO: tess_normals, tess_hemisplit, tess_addsphere

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
NoSave = false;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'nosave', NoSave = varargin{i+1};
    end
end

%% ===== LOAD SURFACE =====
TessFile = file_fullpath(SurfaceFile);
TessMat  = in_tess_bst(SurfaceFile);
Vtx = TessMat.Vertices;
Fcs = double(TessMat.Faces);
nF  = size(Fcs, 1);

%% ===== REQUIRE FREESURFER REGISTRATION SPHERE =====
if ~isfield(TessMat, 'Reg') || ~isstruct(TessMat.Reg) || ...
   ~isfield(TessMat.Reg, 'Sphere') || ~isfield(TessMat.Reg.Sphere, 'Vertices') || ...
   isempty(TessMat.Reg.Sphere.Vertices)
    error('tess_tangents:noRegSphere', ...
        ['Surface has no FreeSurfer registration sphere (Reg.Sphere.Vertices). ' ...
         'Import with surface registration or run tess_addsphere first.']);
end
Sphere = TessMat.Reg.Sphere.Vertices;

%% ===== ENSURE NXR-COMPUTE PLUGIN =====
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
if ~isOk
    error('tess_tangents:nxrUnavailable', ...
        'tess_tangents requires the nxr-compute plugin: %s', errMsg);
end

%% ===== HEMISPHERE SPLIT =====
[rH, lH] = tess_hemisplit(TessMat);
hemis    = {lH(:)', rH(:)'};
hemiTags = {'L', 'R'};

%% ===== FACE NORMALS (for frame orientation) =====
fn = cross(Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:), Vtx(Fcs(:,3),:) - Vtx(Fcs(:,1),:));
fn = fn ./ max(sqrt(sum(fn.^2, 2)), eps);

%% ===== PER-HEMISPHERE TRIVIAL CONNECTION =====
U = zeros(nF, 3);
V = zeros(nF, 3);
assigned  = false(nF, 1);
singVerts = [];
singIdx   = [];
hemiOf    = {};
for h = 1:numel(hemis)
    vH = hemis{h};
    if isempty(vH)
        continue;
    end
    [Uh, Vh, fMask, polesGlobal] = solve_hemisphere(Vtx, Fcs, Sphere, vH);
    if any(assigned & fMask)
        error('tess_tangents:overlap', 'A face was assigned to both hemispheres.');
    end
    U(fMask, :)     = Uh;
    V(fMask, :)     = Vh;
    assigned(fMask) = true;
    singVerts = [singVerts; polesGlobal(:)];          %#ok<AGROW>
    singIdx   = [singIdx;   [1; 1]];                  %#ok<AGROW>
    hemiOf    = [hemiOf, hemiTags(h), hemiTags(h)];   %#ok<AGROW>
end

%% ===== COVERAGE CHECK =====
if ~all(assigned)
    error('tess_tangents:unassignedFaces', ...
        '%d of %d faces were not assigned to a hemisphere.', nnz(~assigned), nF);
end

%% ===== ORIENT FRAME (right-handed w.r.t. face normal) =====
flip = sum(cross(U, V, 2) .* fn, 2) < 0;
V(flip, :) = -V(flip, :);

%% ===== STORE =====
if ~NoSave
    Sing = struct();
    Sing.Vertices   = singVerts;
    Sing.Indices    = singIdx;
    Sing.Hemisphere = hemiOf;
    TF = struct();
    TF.Domain        = 'face';
    TF.U             = single(U);
    TF.V             = single(V);
    TF.Singularities = Sing;
    TF.Method        = 'nxr trivial-connection (FreeSurfer poles)';
    TessMat.TangentFrame = TF;
    TessMat = bst_history('add', TessMat, 'tangents', ...
        'Computed per-face tangent frame field (nxr trivial connection, FreeSurfer poles).');
    bst_save(TessFile, TessMat, 'v7');
end
end


%% ========================================================================
function [Uh, Vh, fMask, polesGlobal] = solve_hemisphere(Vtx, Fcs, Sphere, vH)
% Solve the trivial-connection frame on one hemisphere submesh.
%   Returns the per-face direction (Uh) and orthogonal (Vh) vectors for the
%   faces selected by fMask (into the full face array), plus the two pole
%   vertex indices (global) used as singularities.
nVtot = size(Vtx, 1);
isVH  = false(nVtot, 1);
isVH(vH) = true;
% Faces wholly inside this hemisphere (valid: hemispheres are disconnected).
fMask = all(isVH(Fcs), 2);
% Local re-indexing (global vertex -> 1-based local index).
map = zeros(nVtot, 1);
map(vH) = 1:numel(vH);
Floc = map(Fcs(fMask, :));
Vloc = Vtx(vH, :);
% Poles from the registration sphere restricted to this hemisphere.
sph = Sphere(vH, :);
[~, iN] = max(sph(:, 3));   % north pole (local index)
[~, iS] = min(sph(:, 3));   % south pole (local index)
% Trivial connection: two +1 singularities, sum = chi = 2.
mctx = nxr.manifold.context(Vloc, Floc);
r = nxr.manifold.interpolate.trivial(mctx, [iN; iS], [1; 1]);
if ~r.gaussBonnetSatisfied
    error('tess_tangents:gaussBonnet', ...
        'Gauss-Bonnet not satisfied for hemisphere (chi=%.3f, expected 2).', ...
        r.eulerCharacteristic);
end
Uh = r.directionVectors;
Vh = r.orthogonalVectors;
polesGlobal = vH([iN, iS]);
end
```

- [ ] **Step 4: Run the guard test to confirm it passes**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); addpath(genpath('dev/tests')); clear functions; test_tess_tangents_guard
```
Expected: PASS — `ALL TESTS PASSED: test_tess_tangents_guard`. (The Reg.Sphere check fires before the nxr install, so this is fast and needs no plugin.)

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_tangents.m dev/tests/test_tess_tangents_guard.m
git commit -m "tess_tangents: per-hemisphere trivial-connection tangent frame + guard test"
```

---

## Task 2: Integration test on a registered cortex

**Files:**
- Create: `dev/tests/test_tess_tangents.m`

- [ ] **Step 1: Write the integration test**

Create `dev/tests/test_tess_tangents.m`:

```matlab
function test_tess_tangents
% Integration test for tess_tangents on a real FreeSurfer-registered cortex:
% per-face orthonormal frame, right-handedness, determinism, and storage
% round-trip. Runs on a temp COPY so the DB surface is not mutated.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% --- Find a low-res FreeSurfer-registered cortex in the DB ---
srcFile = find_registered_cortex();
assert(~isempty(srcFile), 'No FreeSurfer-registered cortex found in the DB to test against.');
fprintf('Source cortex: %s\n', srcFile);

% --- Work on a temp copy (named tess_cortex_*.mat so file_fullpath accepts it) ---
tmpFile = fullfile(tempdir, 'tess_cortex_tangents_test.mat');
copyfile(file_fullpath(srcFile), tmpFile);
cleanup = onCleanup(@() delete(tmpFile));

% --- Compute + store ---
[U, V] = tess_tangents(tmpFile);
TessMat = in_tess_bst(tmpFile);
Vtx = TessMat.Vertices;  Fcs = double(TessMat.Faces);
nF  = size(Fcs, 1);
assert(isequal(size(U), [nF 3]), 'U must be nF x 3.');
assert(isequal(size(V), [nF 3]), 'V must be nF x 3.');

% --- Orthonormal, right-handed frame ---
assert(all(abs(sqrt(sum(U.^2,2)) - 1) < 1e-4), 'U not unit length.');
assert(all(abs(sqrt(sum(V.^2,2)) - 1) < 1e-4), 'V not unit length.');
assert(max(abs(sum(U.*V, 2))) < 1e-4, 'U and V not orthogonal.');
fn = cross(Vtx(Fcs(:,2),:) - Vtx(Fcs(:,1),:), Vtx(Fcs(:,3),:) - Vtx(Fcs(:,1),:));
fn = fn ./ max(sqrt(sum(fn.^2,2)), eps);
assert(all(sum(cross(U, V, 2) .* fn, 2) > -1e-6), 'Frame not right-handed w.r.t. face normal.');
fprintf('PASSED: per-face orthonormal, right-handed frame (%d faces).\n', nF);

% --- Storage round-trip ---
assert(isfield(TessMat, 'TangentFrame'), 'TangentFrame not stored.');
TF = TessMat.TangentFrame;
assert(strcmp(TF.Domain, 'face'), 'TangentFrame.Domain should be ''face''.');
assert(isequal(single(U), TF.U) && isequal(single(V), TF.V), 'Stored U/V mismatch.');
assert(numel(TF.Singularities.Vertices) == 4, 'Expected 4 singularities (2 per hemisphere).');
assert(all(TF.Singularities.Indices == 1), 'Singularity indices should all be +1.');
fprintf('PASSED: TangentFrame stored + reloaded (4 pole singularities).\n');

% --- Determinism ---
[U2, V2] = tess_tangents(tmpFile, 'NoSave', 1);
assert(isequal(U, U2) && isequal(V, V2), 'tess_tangents is not deterministic.');
fprintf('PASSED: deterministic across runs.\n');

fprintf('ALL TESTS PASSED: test_tess_tangents\n');
end


function SurfaceFile = find_registered_cortex()
% Return a low-res cortex FileName that has Reg.Sphere.Vertices, or '' if none.
SurfaceFile = '';
best = inf;
sSubjects = bst_get('ProtocolSubjects');
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try
            T = load(file_fullpath(surf(iF).FileName), 'Reg', 'Vertices');
        catch
            continue;
        end
        if isfield(T,'Reg') && isstruct(T.Reg) && isfield(T.Reg,'Sphere') ...
           && isfield(T.Reg.Sphere,'Vertices') && ~isempty(T.Reg.Sphere.Vertices)
            n = size(T.Vertices, 1);
            if n < best          % prefer the smallest (fastest to solve)
                best = n;
                SurfaceFile = surf(iF).FileName;
            end
        end
    end
end
end
```

- [ ] **Step 2: Run it**

Run (MATLAB MCP `evaluate_matlab_code`):
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); addpath(genpath('dev/tests')); clear functions; test_tess_tangents
```
Expected: PASS — three `PASSED` lines and `ALL TESTS PASSED: test_tess_tangents`. The run prints the chosen source cortex; the solve runs once per hemisphere.

- [ ] **Step 3: If a parity/Gauss-Bonnet assertion fails**

Do NOT loosen tolerances. Capture and report:
- If `tess_tangents:gaussBonnet` is raised: print the hemisphere's vertex/face counts and `r.eulerCharacteristic` — a χ ≠ 2 means the hemisphere submesh is not a clean genus-0 sphere (non-manifold or has boundary). Report for a planning decision.
- If orthogonality fails: print `max(abs(sum(U.*V,2)))` — nxr's `orthogonalVectors` should be exactly orthogonal per face; a large value indicates a marshalling/indexing problem.

- [ ] **Step 4: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_tess_tangents.m
git commit -m "tess_tangents: integration test (frame orthonormality, determinism, storage)"
```

---

## Task 3: Final integration — run both tests

- [ ] **Step 1: Run the full tess_tangents suite**

Run (MATLAB MCP `evaluate_matlab_code`), each expecting `ALL TESTS PASSED`:
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); addpath(genpath('dev/tests'));
for t = {'test_tess_tangents_guard','test_tess_tangents'}
    clear functions; fprintf('\n== %s ==\n', t{1}); feval(t{1});
end
```

- [ ] **Step 2: Confirm branch state**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git status -s
git log --oneline development..HEAD
```
Expected: clean working tree (no tess_tangents-related uncommitted changes); commits from Tasks 1–2 on `feature/tess-tangents`. Push only on explicit request.

---

## Self-Review

**Spec coverage:**
- §2 per-face frame, stored → Task 1 (`TangentFrame.Domain='face'`, `U/V`). ✓
- §3.3 per-hemisphere solve → Task 1 (`solve_hemisphere`, loop over `{lH, rH}`). ✓
- §4 SurfaceFile interface + `NoSave` → Task 1 (signature, arg parse). ✓
- §5 algorithm (load, Reg guard, nxr guard, hemisplit, submesh, poles, trivial, scatter, coverage, orient, store) → Task 1, all steps present. ✓
- §6 storage struct (`Domain/U/V/Singularities/Method` + History) → Task 1. ✓
- §7 nxr-only guarded + clean errors (`noRegSphere`, `nxrUnavailable`, `gaussBonnet`, `unassignedFaces`, `overlap`) → Task 1. ✓
- §8 tests: integration on registered cortex + missing-registration guard → Tasks 1 (guard) & 2 (integration); fixture resolved (DB cortex discovered at runtime, no synthetic fixture needed). ✓
- §9 out-of-scope (vertex frame, source-map analysis, GUI) → no tasks (correct). ✓

**Placeholder scan:** No TBD/TODO/"handle edge cases". All code blocks complete; the test discovers its fixture at runtime. ✓

**Type/name consistency:** `U`/`V` (`[nF×3]`), `solve_hemisphere` returns `[Uh, Vh, fMask, polesGlobal]`, `TangentFrame` fields `Domain/U/V/Singularities/Method`, error ids `tess_tangents:{noRegSphere,nxrUnavailable,gaussBonnet,unassignedFaces,overlap}` — used consistently across the function and both tests. `tess_hemisplit` return order `[rH, lH]` matches the verified signature (right, left), and the loop pairs `{lH, rH}` with tags `{'L','R'}` correctly. ✓
