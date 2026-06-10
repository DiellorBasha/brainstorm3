# `tess_frame` Facet Bundle (Phase 2) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Repurpose `tess_frame` into a self-contained function that computes the nxr facet bundle per hemisphere, stores it on the surface file as five top-level `1×2` fields (`Topology`/`Embedded`/`Intrinsic`/`Extrinsic`/`Gauge`), and returns the full-mesh 3D frame `{U,V,N}`.

**Architecture:** `tess_frame` loads the surface with Brainstorm I/O (`in_tess_bst`), splits hemispheres with `tess_hemisplit` (atlas L/R — never `conncomp`), runs `nxr_compute('facets', …)` on each hemisphere submesh with FreeSurfer-pole singularities for the trivial gauge, attaches scatter maps, stores the five groups via `bst_save('v7')`, and derives `{U,V,N}` from `Embedded.vertex.grid ⊙ Gauge.vertex.rotation`. No shared helper (`tess_store_perhemi` is intentionally NOT used). Phase 1 (the nxr `facets` command) is already merged and installed.

**Tech Stack:** MATLAB (Brainstorm toolbox), `nxr-compute` MEX (`facets` command), `matlab.unittest`, run via the MATLAB MCP.

**Prerequisite (already satisfied):** the installed `~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/nxr_compute.mexmaca64` exposes `nxr_compute('facets', h, gaugeType[, opts])` returning `{Topology, Embedded, Intrinsic, Extrinsic, Gauge}` (verified: smoke-tested on the installed binary).

**Reference facts (grounded against the code):**
- `tess_hemisplit` signature: `[rH, lH, isConnected, iStruct, iRightScout, iLeftScout] = tess_hemisplit(sSurf)` — `lH`=left, `rH`=right.
- Brainstorm I/O idiom (`tess_addsphere`): `TessMat = in_tess_bst(TessFile)` → modify → `bst_save(file_fullpath(TessFile), TessMat, 'v7')`, with `bst_history`.
- The `facets` Gauge sub-struct: `Gauge.type`, `Gauge.vertex.rotation` (`[nVh×1]` complex), `Gauge.face.rotation` (empty/deferred for trivial), `Gauge.singularity.{vertices,indices,source}`. `Embedded.vertex.grid` is the base (un-gauged) complex frame `e1+i·e2`; the resolved frame is `grid .* rotation`.
- No external code calls `tess_frame` (only its own help + its test) — safe to change its contract.

**MATLAB-session discipline:** before each test run, `rehash; clear test_tess_frame tess_frame;`. **Never** a bare `clear` (wipes Brainstorm `GlobalData`). Run via the MATLAB MCP (`mcp__plugin_brainstorm-dev_MATLAB__evaluate_matlab_code`). Tests that write to disk back up the `.mat` and restore via `onCleanup`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/anatomy/tess_frame.m` | **Replace.** Self-contained compute → store five facet groups → return `{U,V,N}`. |
| `dev/tests/test_tess_frame.m` | **Replace.** Property tests on the real 20484-vertex cortex (compute/store/scatter/derive/cache/NoSave/face-error). |

---

## Task 1: Repurpose `tess_frame` (compute + store + derive)

**Files:**
- Replace: `toolbox/anatomy/tess_frame.m`
- Replace: `dev/tests/test_tess_frame.m`

- [ ] **Step 1: Write the failing test (replace the whole file)**

Replace `dev/tests/test_tess_frame.m` with:

```matlab
function tests = test_tess_frame
% Tests for tess_frame: compute/store the nxr facet bundle + return {U,V,N}.
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

function test_store_five_groups_1x2(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for f = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'}
        verifyTrue(tc, isfield(T, f{1}) && isequal(size(T.(f{1})), [1 2]), ...
            sprintf('%s is 1x2', f{1}));
    end
    verifyEqual(tc, T.Embedded(1).Hemisphere, 'L');
    verifyEqual(tc, T.Embedded(2).Hemisphere, 'R');
    verifyEqual(tc, T.Gauge(1).Provenance.Backend, 'nxr');
    verifyEqual(tc, lower(T.Gauge(1).type), 'trivial');
end

function test_scatter_partition(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = sort([T.Embedded(1).GlobalVertices(:); T.Embedded(2).GlobalVertices(:)]);
    gf = sort([T.Embedded(1).GlobalFaces(:);    T.Embedded(2).GlobalFaces(:)]);
    verifyEqual(tc, gv, (1:size(T.Vertices,1))');
    verifyEqual(tc, gf, (1:size(T.Faces,1))');
end

function test_frame_orthonormal_fullmesh(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    [U,V,N] = tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    nV = size(T.Vertices,1);
    verifyEqual(tc, size(U), [nV 3]);
    verifyEqual(tc, size(V), [nV 3]);
    verifyEqual(tc, size(N), [nV 3]);
    verifyLessThan(tc, max(abs(sqrt(sum(U.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sqrt(sum(V.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sum(U.*V,2))), 1e-4);
    verifyLessThan(tc, max(abs(N - cross(U,V,2)), [], 'all'), 1e-4);
    verifyGreaterThan(tc, min(sqrt(sum(N.^2,2))), 0.9);   % every row filled
end

function test_frame_matches_grid_rotation(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    [U,V,~] = tess_frame(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        idx = T.Embedded(hh).GlobalVertices;
        c   = T.Embedded(hh).vertex.grid .* T.Gauge(hh).vertex.rotation;
        verifyLessThan(tc, max(abs(U(idx,:) - real(c)), [], 'all'), 1e-5);
        verifyLessThan(tc, max(abs(V(idx,:) - imag(c)), [], 'all'), 1e-5);
    end
end
```

- [ ] **Step 2: Run the test to verify it fails**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
runtests('dev/tests/test_tess_frame.m')
```
Expected: tests FAIL — the current `tess_frame` is a read-only accessor that errors `tess_frame:noBundle` (no stored bundle) and never writes `Embedded`/`Intrinsic`/`Extrinsic`. (`ForceRecompute` is not even an option yet.)

- [ ] **Step 3: Replace `toolbox/anatomy/tess_frame.m`**

Replace the ENTIRE file (keep a Brainstorm license header block like the current file's) with:

```matlab
function [U, V, N] = tess_frame(SurfaceFile, varargin)
% TESS_FRAME: Compute/store the nxr facet bundle per hemisphere; return the 3D frame.
%
% USAGE:  [U,V,N] = tess_frame(SurfaceFile)
%         [U,V,N] = tess_frame(SurfaceFile, 'Gauge','trivial', 'Domain','vertex', ...
%                              'ForceRecompute',1, 'NoSave',1)
%
% DESCRIPTION:
%     Loads the surface, splits hemispheres with tess_hemisplit (atlas L/R, never
%     conncomp), runs nxr_compute('facets', ...) on each hemisphere submesh (with
%     FreeSurfer-pole singularities for the trivial gauge), and stores the result
%     as five top-level 1x2 per-hemisphere struct arrays on the surface file:
%     TessMat.{Topology, Embedded, Intrinsic, Extrinsic, Gauge}. Each element is
%     the verbatim nxr facet struct (LOCAL indexing) plus GlobalVertices/
%     GlobalFaces/Hemisphere/Provenance scatter maps to the full-mesh order.
%
%     Returns the full-mesh intrinsic frame {U,V,N}: U=real(grid.*rot),
%     V=imag(grid.*rot), N=cross(U,V), where grid=Embedded.vertex.grid and
%     rot=Gauge.vertex.rotation (rot==1 for euclidean/levi-civita).
%
%     If the five fields are already present and ForceRecompute is false, the
%     frame is derived from the stored bundle without recomputing.
%
%     Vertex domain works for every gauge. Face domain works only for
%     euclidean/levi-civita: the trivial gauge's Gauge.face.rotation is deferred
%     in nxr, so face+trivial errors clearly.
%
% Requires the nxr-compute plugin; the trivial gauge needs a FreeSurfer reg sphere.
%
% SEE ALSO: tess_hemisplit, tess_tangents, tess_normals

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
    Gauge='trivial'; Domain='vertex'; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge=lower(varargin{i+1});
            case 'domain',         Domain=lower(varargin{i+1});
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end
    if ~ismember(Domain, {'vertex','face'})
        error('tess_frame:badDomain', 'Domain must be ''vertex'' or ''face''.');
    end

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    Groups = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'};
    haveAll = all(cellfun(@(f) isfield(TessMat,f) && ~isempty(TessMat.(f)), Groups));

    if ForceRecompute || ~haveAll
        TessMat = local_compute_store(TessFile, TessMat, Gauge, NoSave);
    end

    [U, V, N] = local_derive_frame(TessMat, Domain);
end

% ----------------------------------------------------------------------------
function TessMat = local_compute_store(TessFile, TessMat, Gauge, NoSave)
    Groups = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'};

    % require nxr-compute (no MATLAB fallback)
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_frame:nxrUnavailable', 'tess_frame requires nxr-compute: %s', errMsg);
    end

    % trivial gauge needs a FreeSurfer registration sphere
    if strcmpi(Gauge,'trivial')
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_frame:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % hemisphere split from atlas labels (never conncomp)
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_frame:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','facets', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    Arr = struct();
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_frame:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        opts = struct();
        if strcmpi(Gauge,'trivial')
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));   % north pole (local index)
            [~, iS] = min(sph(:,3));   % south pole (local index)
            opts.singVerts  = [iN; iS];
            opts.singValues = [1; 1];
        end

        h = nxr_compute('create', Vloc, Floc);
        S = nxr_compute('facets', h, Gauge, opts);
        nxr_compute('destroy', h);

        for f = 1:numel(Groups)
            s = S.(Groups{f});
            s.GlobalVertices = vH;
            s.GlobalFaces    = find(fMask);
            s.Hemisphere     = tags{hh};
            s.Provenance     = prov;
            Arr.(Groups{f})(hh) = s;   %#ok<AGROW>
        end
    end

    for f = 1:numel(Groups), TessMat.(Groups{f}) = Arr.(Groups{f}); end

    if ~NoSave
        TessMat_full = load(TessFile);
        for f = 1:numel(Groups), TessMat_full.(Groups{f}) = TessMat.(Groups{f}); end
        TessMat_full = bst_history('add', TessMat_full, 'facets', ...
            sprintf(['Stored nxr facet bundle (gauge=%s) as per-hemisphere ' ...
                     'Topology/Embedded/Intrinsic/Extrinsic/Gauge.'], Gauge));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end

% ----------------------------------------------------------------------------
function [U, V, N] = local_derive_frame(TessMat, Domain)
    Emb = TessMat.Embedded; Ga = TessMat.Gauge;
    if strcmp(Domain,'vertex')
        nElem = size(TessMat.Vertices,1);
    else
        nElem = size(TessMat.Faces,1);
    end
    U = zeros(nElem,3); V = zeros(nElem,3);

    for hh = 1:numel(Emb)
        gridH = Emb(hh).(Domain).grid;            % nElemH x 3 complex
        if strcmpi(Ga(hh).type, 'trivial')
            if strcmp(Domain,'face')
                if ~isfield(Ga(hh).face,'rotation') || isempty(Ga(hh).face.rotation)
                    error('tess_frame:faceTrivialDeferred', ...
                        'Face-domain trivial frame needs Gauge.face.rotation (empty/deferred in nxr).');
                end
                rot = Ga(hh).face.rotation;
            else
                rot = Ga(hh).vertex.rotation;
            end
        else
            rot = ones(size(gridH,1),1);
        end
        cRot = gridH .* rot;                       % broadcast over the 3 columns
        if strcmp(Domain,'vertex')
            idx = Emb(hh).GlobalVertices;
        else
            idx = Emb(hh).GlobalFaces;
        end
        U(idx,:) = real(cRot);
        V(idx,:) = imag(cRot);
    end

    N = cross(U, V, 2);
end
```

- [ ] **Step 4: Run the test to verify it passes**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
runtests('dev/tests/test_tess_frame.m')
```
Expected: 4 tests PASS.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_frame.m dev/tests/test_tess_frame.m
git commit -m "feat(tess-frame): self-contained facet-bundle compute/store + {U,V,N}"
```

---

## Task 2: Cache-return, NoSave, and face-domain trivial error

**Files:**
- Modify: `dev/tests/test_tess_frame.m` (add tests; implementation already supports these)

- [ ] **Step 1: Add the tests**

Append these three test functions to `dev/tests/test_tess_frame.m` (before the final `end` of the file is not needed — these are top-level `local` test functions; place them after `test_frame_matches_grid_rotation`):

```matlab
function test_cache_return_no_recompute(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);
    % stamp a sentinel into the stored provenance, then call without ForceRecompute
    TF = load(TessFile);
    TF.Embedded(1).Provenance.ComputeDate = 'SENTINEL';
    bst_save(TessFile, TF, 'v7');

    [U,~,~] = tess_frame(SurfaceFile);          % must cache-return (no recompute)
    T = in_tess_bst(SurfaceFile, 0);
    verifyEqual(tc, T.Embedded(1).Provenance.ComputeDate, 'SENTINEL');   % unchanged => not recomputed
    verifyEqual(tc, size(U), [size(T.Vertices,1) 3]);
end

function test_nosave_returns_without_writing(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    % clean slate: strip the five fields if present
    TF = load(TessFile);
    for f = {'Topology','Embedded','Intrinsic','Extrinsic','Gauge'}
        if isfield(TF, f{1}), TF = rmfield(TF, f{1}); end
    end
    bst_save(TessFile, TF, 'v7');

    [U,~,~] = tess_frame(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyEqual(tc, size(U,2), 3);
    verifyFalse(tc, isfield(T,'Embedded') && ~isempty(T.Embedded));   % not written
end

function test_face_trivial_errors(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>

    tess_frame(SurfaceFile, 'ForceRecompute', 1);   % ensure bundle present (vertex)
    verifyError(tc, @() tess_frame(SurfaceFile, 'Domain','face'), 'tess_frame:faceTrivialDeferred');
end
```

- [ ] **Step 2: Run the full test file**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
runtests('dev/tests/test_tess_frame.m')
```
Expected: 7 tests PASS (the implementation from Task 1 already supports cache-return, NoSave, and the face-trivial error — these tests lock that behavior in). If `test_cache_return_no_recompute` fails, confirm the cache guard checks all five fields and returns before `local_compute_store`.

- [ ] **Step 3: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add dev/tests/test_tess_frame.m
git commit -m "test(tess-frame): cache-return, NoSave, face-trivial error"
```

---

## Final verification

- [ ] **Run the test file; expect 7 passing.**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
r = runtests('dev/tests/test_tess_frame.m');
fprintf('\nTOTAL: %d passed, %d failed\n', sum([r.Passed]), sum([r.Failed]));
assert(all([r.Passed]) && ~any([r.Failed]), 'Some tests failed.');
```

- [ ] **Then complete via superpowers:finishing-a-development-branch.**

---

## Notes for the implementer

- **Single source of truth:** the five stored groups are the verbatim nxr `facets` output per hemisphere plus the four scatter/provenance fields. Do not post-process the nxr structs.
- **Gauge default is `trivial`** with FreeSurfer-pole singularities (`singValues=[1;1]`, Gauss-Bonnet χ=2 per genus-0 hemisphere). For `euclidean`/`levi-civita`, `opts` stays empty and `Gauge.vertex.rotation` is identity, so the frame is just `grid`.
- **Face domain** intentionally errors under the trivial gauge (`Gauge.face.rotation` is deferred/empty in nxr) — this is expected behavior, not a bug to work around.
- **Isolation:** every test backs up the full `.mat` and restores via `onCleanup`. Do not weaken the storage tests to `NoSave`.
- **Out of scope:** deleting `tess_store_perhemi` and the four `tess_topology`/`tess_geometry`/`tess_gauge`/`tess_bundle` writers (separate cleanup). They still produce the older `Topology`/`Geometry`/`Gauge` fields; `tess_frame` now additionally produces `Embedded`/`Intrinsic`/`Extrinsic` and overwrites `Topology`/`Gauge` with the (compatible) `facets` versions on recompute.
```
