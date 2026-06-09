# Surface-File Bundle Storage Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Store the nxr-compute coordinate bundle (`Topology`/`Geometry`/`Gauge`) on the Brainstorm surface file as per-hemisphere `1×2` struct arrays, with a derived `{U,V,N}` frame accessor.

**Architecture:** One writer `tess_bundle.m` splits the cortex by hemisphere, calls `nxr_compute('bundle', …)` on each submesh (nxr never integrates hemispheres), and stores each result verbatim (local indexing) plus `GlobalVertices`/`GlobalFaces` scatter maps into `TessMat.{Topology,Geometry,Gauge}(h)`. Light by default; heavy `.operators` opt-in. A pure accessor `tess_frame.m` rebuilds the full-mesh `{U,V,N}` frame from `grid ⊙ gauge-rotation`.

**Tech Stack:** MATLAB (Brainstorm toolbox), `nxr-compute` MEX (geometry-central), `matlab.unittest`, run via the MATLAB MCP.

**Prerequisite (must be true before executing):** the **new** `nxr-compute` MEX (with the `bundle` command) is built and installed into `~/.brainstorm/plugins/nxr-compute/…/nxr_compute.mexmaca64`, and `nxr_compute('version')` / a `bundle` smoke call succeeds in the running session. The stale MEX rejects `bundle` ("Unknown command"), which would make every test below fail at dispatch, not at our logic.

**Reference patterns (read before starting):**
- `toolbox/anatomy/tess_coordinates.m` — option parsing, `Reg.Sphere` guard, hemisphere split, FreeSurfer-pole singularities, nxr handle lifecycle, load-full→set→`bst_history`→`bst_save 'v7'`.
- `toolbox/anatomy/tess_operators.m` — `NoSave`/`ForceRecompute`/cache-return shape, nxr-required guard.
- `dev/tests/test_tess_coordinates.m` — `local_find_cortex(20484)`, `brainstorm nogui`, backup/restore file isolation via `onCleanup`.
- nxr `CLAUDE.md` "MATLAB coordinate-system bundle" — the `bundle` command contract.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/anatomy/tess_bundle.m` | **Create.** Per-hemisphere nxr `bundle` → `TessMat.{Topology,Geometry,Gauge}` (1×2). Light default, heavy opt-in. |
| `toolbox/anatomy/tess_frame.m` | **Create.** Pure accessor: full-mesh `{U,V,N}` from stored `Geometry.grid ⊙ Gauge.rotation`. |
| `dev/tests/test_tess_bundle.m` | **Create.** Property tests (light + heavy) on the real cortex. |
| `dev/tests/test_tess_frame.m` | **Create.** Frame-derivation tests (vertex domain; face-trivial error). |

Each test file carries its own `local_find_cortex(20484)` helper (copied verbatim from `test_tess_coordinates.m`), so the files are independent.

---

## Task 1: `tess_bundle` — light, per-hemisphere store

**Files:**
- Create: `toolbox/anatomy/tess_bundle.m`
- Test: `dev/tests/test_tess_bundle.m`

- [ ] **Step 1: Write the failing test file (light-path tests)**

Create `dev/tests/test_tess_bundle.m`:

```matlab
function tests = test_tess_bundle
% Property tests for tess_bundle (per-hemisphere nxr bundle store) on real cortex.
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

function test_three_struct_arrays_1x2(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(B.Topology), [1 2]);
    verifyEqual(tc, size(B.Geometry), [1 2]);
    verifyEqual(tc, size(B.Gauge),    [1 2]);
    verifyEqual(tc, B.Gauge(1).Provenance.Backend, 'nxr');   % writer-added provenance
end

function test_scatter_maps_partition_mesh(tc)
    SurfaceFile = local_cortex();
    B = tess_bundle(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = [B.Geometry(1).GlobalVertices(:); B.Geometry(2).GlobalVertices(:)];
    gf = [B.Geometry(1).GlobalFaces(:);    B.Geometry(2).GlobalFaces(:)];
    % disjoint + complete partition of vertices and faces
    verifyEqual(tc, sort(gv), (1:size(T.Vertices,1))');
    verifyEqual(tc, sort(gf), (1:size(T.Faces,1))');
    verifyEqual(tc, B.Geometry(1).Hemisphere, 'L');
    verifyEqual(tc, B.Geometry(2).Hemisphere, 'R');
end

function test_no_operators_in_light(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyFalse(tc, isfield(B.Topology(1), 'operators') && ~isempty(B.Topology(1).operators));
    verifyFalse(tc, isfield(B.Geometry(1), 'operators') && ~isempty(B.Geometry(1).operators));
    verifyFalse(tc, isfield(B.Gauge(1),    'operators') && ~isempty(B.Gauge(1).operators));
end

function test_vertex_grid_orthonormal(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        c  = B.Geometry(hh).vertex.grid;        % nVh x 3 complex
        e1 = real(c); e2 = imag(c);
        n1 = sqrt(sum(e1.^2,2)); n2 = sqrt(sum(e2.^2,2));
        verifyLessThan(tc, max(abs(n1-1)), 1e-4);
        verifyLessThan(tc, max(abs(n2-1)), 1e-4);
        verifyLessThan(tc, max(abs(sum(e1.*e2,2))), 1e-4);
    end
end

function test_trivial_gauge_gauss_bonnet(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);   % default gauge='trivial'
    for hh = 1:2
        verifyEqual(tc, lower(B.Gauge(hh).type), 'trivial');
        verifyEqual(tc, sum(B.Gauge(hh).singularity.indices), 2, 'AbsTol', 1e-6);
    end
end

function test_field_is_stored(tc)
    % Real save path, then restore the file (isolation).
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_bundle(SurfaceFile, 'ForceRecompute', 1);   % default: save
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Topology') && isequal(size(T.Topology),[1 2]));
    verifyTrue(tc, isfield(T,'Geometry') && isequal(size(T.Geometry),[1 2]));
    verifyTrue(tc, isfield(T,'Gauge')    && isequal(size(T.Gauge),[1 2]));
end
```

- [ ] **Step 2: Run the tests to verify they fail**

In the MATLAB MCP session:
```matlab
rehash; clear test_tess_bundle tess_bundle; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: all tests ERROR with `Undefined function 'tess_bundle'` (not yet created).

- [ ] **Step 3: Write the minimal `tess_bundle.m` (light path)**

Create `toolbox/anatomy/tess_bundle.m`:

```matlab
function B = tess_bundle(SurfaceFile, varargin)
% TESS_BUNDLE: Store the nxr-compute coordinate bundle as per-hemisphere fields.
%
% USAGE:  B = tess_bundle(SurfaceFile)
%         B = tess_bundle(SurfaceFile, 'Gauge','trivial', 'NoSave',1, 'ForceRecompute',1)
%
% DESCRIPTION:
%     nxr-compute never integrates the two hemispheres, so each cortex is bundled
%     one connected genus-0 hemisphere at a time. The three results are stored as
%     1x2 struct arrays ((1)=left, (2)=right). Each element is the verbatim nxr
%     submesh bundle (LOCAL indexing) plus GlobalVertices/GlobalFaces/Hemisphere
%     scatter maps to TessMat.Vertices/.Faces global order.
%
% OPTIONS:
%     'Gauge'           'trivial' (default) | 'levi-civita' | 'euclidean'
%     'Operators'       false (default) | true   - heavy .operators (Task 2)
%     'Coupling'        'ambient' (default) | 'product'   - covariantLaplacian (heavy)
%     'Mass'            'lumped' (default) | 'galerkin'    - mass variant (heavy)
%     'NoSave'          false (default) | true
%     'ForceRecompute'  false (default) | true
%
% OUTPUT (also stored as TessMat.Topology/.Geometry/.Gauge, each 1x2):
%     B.Topology, B.Geometry, B.Gauge   - 1x2 struct arrays
%
% Requires the nxr-compute plugin; trivial gauge needs a FreeSurfer reg sphere.
%
% Authors: Diellor Basha, 2026

    % --- options ---
    Gauge='trivial'; Operators=false; Coupling='ambient'; Mass='lumped';
    NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          Gauge=varargin{i+1};
            case 'operators',      Operators=varargin{i+1};
            case 'coupling',       Coupling=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return (light only; heavy is opt-in and always recomputes) ---
    if ~ForceRecompute && ~Operators && isfield(TessMat,'Topology') && ~isempty(TessMat.Topology)
        B = struct();                          % explicit assignment: struct() ctor
        B.Topology = TessMat.Topology;         % mishandles struct-array field values
        B.Geometry = TessMat.Geometry;
        B.Gauge    = TessMat.Gauge;
        return;
    end

    % --- trivial gauge needs a FreeSurfer registration sphere ---
    if strcmpi(Gauge,'trivial')
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_bundle:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % --- require nxr-compute (no MATLAB fallback) ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_bundle:nxrUnavailable', 'tess_bundle requires nxr-compute: %s', errMsg);
    end
    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    % --- hemisphere split from import labels (not geometry) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_bundle:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};

    Vtx = double(TessMat.Vertices);
    Fcs = double(TessMat.Faces);
    nVtot = size(Vtx,1);

    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_bundle:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        % nxr submesh bundle
        h = nxr_compute('create', Vloc, Floc);
        opts = struct();
        if strcmpi(Gauge,'trivial')
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));   % north pole (local index)
            [~, iS] = min(sph(:,3));   % south pole (local index)
            opts.singVerts  = [iN; iS];
            opts.singValues = [1; 1];
        end
        Bh = nxr_compute('bundle', h, Gauge, opts);
        nxr_compute('destroy', h);

        % attach scatter maps; build the 1x2 struct arrays
        sT = Bh.Topology; sT.GlobalVertices = vH; sT.GlobalFaces = find(fMask); sT.Hemisphere = tags{hh}; sT.Provenance = prov;
        sG = Bh.Geometry; sG.GlobalVertices = vH; sG.GlobalFaces = find(fMask); sG.Hemisphere = tags{hh}; sG.Provenance = prov;
        sA = Bh.Gauge;    sA.GlobalVertices = vH; sA.GlobalFaces = find(fMask); sA.Hemisphere = tags{hh}; sA.Provenance = prov;
        TopoArr(hh) = sT;  GeoArr(hh) = sG;  GaArr(hh) = sA;  %#ok<AGROW>
    end

    B = struct();                  % explicit assignment (struct() ctor mishandles
    B.Topology = TopoArr;          % struct-array field values)
    B.Geometry = GeoArr;
    B.Gauge    = GaArr;

    % --- save ---
    if ~NoSave
        TessMat_full = load(TessFile);
        TessMat_full.Topology = B.Topology;
        TessMat_full.Geometry = B.Geometry;
        TessMat_full.Gauge    = B.Gauge;
        TessMat_full = bst_history('add', TessMat_full, 'bundle', ...
            sprintf('Stored nxr bundle (gauge=%s) as per-hemisphere Topology/Geometry/Gauge.', Gauge));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end
```

- [ ] **Step 4: Run the tests to verify they pass**

```matlab
rehash; clear test_tess_bundle tess_bundle; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: all 6 tests PASS. If `test_scatter_maps_partition_mesh` fails on the union check, confirm `tess_hemisplit` returns disjoint vertex sets covering all vertices for this surface (it does for FreeSurfer-imported cortex).

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_bundle.m dev/tests/test_tess_bundle.m
git commit -m "feat(tess-bundle): per-hemisphere nxr bundle store (light)"
```

---

## Task 2: `tess_bundle` — heavy `.operators` opt-in

**Files:**
- Modify: `toolbox/anatomy/tess_bundle.m` (per-hemisphere loop: pass `operators`/`coupling`/`mass` to nxr)
- Test: `dev/tests/test_tess_bundle.m` (add heavy-path tests)

- [ ] **Step 1: Add the failing heavy-path tests**

Append to `dev/tests/test_tess_bundle.m`:

```matlab
function test_heavy_dec_identity(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1);
    for hh = 1:2
        op = B.Topology(hh).operators;
        Z = op.dec.d1 * op.dec.d0;            % d o d = 0
        verifyLessThan(tc, full(max(abs(Z(:)))), 1e-9);
    end
end

function test_heavy_connection_laplacian_hermitian(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1);
    for hh = 1:2
        K = B.Gauge(hh).operators.laplacian;   % complex nVh x nVh
        verifyTrue(tc, issparse(K) && ~isreal(K));            % complex sparse
        verifyLessThan(tc, full(max(abs(K - K'), [], 'all')), 1e-6);   % Hermitian
    end
end

function test_heavy_covariant_laplacian_shape_sym(tc)
    B = tess_bundle(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1);
    for hh = 1:2
        nVh = numel(B.Gauge(hh).GlobalVertices);
        L3  = B.Gauge(hh).operators.covariantLaplacian;       % 3N x 3N real
        verifyEqual(tc, size(L3), [3*nVh, 3*nVh]);
        verifyTrue(tc, isreal(L3));
        verifyLessThan(tc, full(max(abs(L3 - L3.'), [], 'all')), 1e-9);   % symmetric
    end
end
```

- [ ] **Step 2: Run to verify the new tests fail**

```matlab
rehash; clear test_tess_bundle tess_bundle; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: the 3 new tests FAIL — the light writer does not request `.operators`, so `B.Topology(hh).operators` is absent (`Reference to non-existent field 'operators'`). The 6 light tests still PASS.

- [ ] **Step 3: Pass operators options through to nxr**

In `toolbox/anatomy/tess_bundle.m`, inside the per-hemisphere loop, extend the `opts` build (immediately after the trivial-gauge `if` block, before the `nxr_compute('bundle', …)` call):

```matlab
        if Operators
            opts.operators = true;
            opts.coupling  = Coupling;
            opts.mass      = Mass;
        end
```

- [ ] **Step 4: Run the full file to verify all pass**

```matlab
rehash; clear test_tess_bundle tess_bundle; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: all 9 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_bundle.m dev/tests/test_tess_bundle.m
git commit -m "feat(tess-bundle): heavy .operators opt-in (DEC/K/covariantLaplacian)"
```

---

## Task 3: `tess_frame` — derived `{U,V,N}` accessor

**Files:**
- Create: `toolbox/anatomy/tess_frame.m`
- Test: `dev/tests/test_tess_frame.m`

- [ ] **Step 1: Write the failing test file**

Create `dev/tests/test_tess_frame.m`:

```matlab
function tests = test_tess_frame
% Tests for tess_frame (derived {U,V,N} from the stored bundle) on real cortex.
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

function local_ensure_bundle(SurfaceFile)
    % Make sure a trivial-gauge bundle is on the file; restore is the caller's job.
    T = in_tess_bst(SurfaceFile, 0);
    if ~isfield(T,'Geometry') || isempty(T.Geometry)
        tess_bundle(SurfaceFile, 'ForceRecompute', 1);
    end
end

function test_frame_fullmesh_orthonormal(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);

    [U,V,N] = tess_frame(SurfaceFile);          % default vertex domain
    T = in_tess_bst(SurfaceFile, 0);
    nV = size(T.Vertices,1);
    verifyEqual(tc, size(U), [nV 3]);
    verifyEqual(tc, size(V), [nV 3]);
    verifyEqual(tc, size(N), [nV 3]);
    % orthonormal, right-handed, no unfilled rows
    verifyLessThan(tc, max(abs(sqrt(sum(U.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sqrt(sum(V.^2,2))-1)), 1e-4);
    verifyLessThan(tc, max(abs(sum(U.*V,2))), 1e-4);
    verifyLessThan(tc, max(abs(N - cross(U,V,2)), [], 'all'), 1e-4);
    verifyGreaterThan(tc, min(sqrt(sum(N.^2,2))), 0.9);    % every row filled
end

function test_frame_matches_grid_rotation(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);

    [U,V,~] = tess_frame(SurfaceFile);
    T = in_tess_bst(SurfaceFile, 0);
    for hh = 1:2
        idx = T.Geometry(hh).GlobalVertices;
        c   = T.Geometry(hh).vertex.grid .* T.Gauge(hh).vertex.rotation;   % grid (x) rotation
        verifyLessThan(tc, max(abs(U(idx,:) - real(c)), [], 'all'), 1e-5);
        verifyLessThan(tc, max(abs(V(idx,:) - imag(c)), [], 'all'), 1e-5);
    end
end

function test_face_trivial_errors(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    local_ensure_bundle(SurfaceFile);
    verifyError(tc, @() tess_frame(SurfaceFile, 'Domain','face'), 'tess_frame:faceTrivialDeferred');
end
```

- [ ] **Step 2: Run to verify failure**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
runtests('dev/tests/test_tess_frame.m')
```
Expected: tests ERROR with `Undefined function 'tess_frame'`.

- [ ] **Step 3: Write `tess_frame.m`**

Create `toolbox/anatomy/tess_frame.m`:

```matlab
function [U, V, N] = tess_frame(SurfaceFile, varargin)
% TESS_FRAME: Derived per-element intrinsic frame {U,V,N} from the stored bundle.
%
% USAGE:  [U,V,N] = tess_frame(SurfaceFile)                    % 'vertex' domain
%         [U,V,N] = tess_frame(SurfaceFile, 'Domain','face')
%
% DESCRIPTION:
%     Reads the per-hemisphere TessMat.Geometry/.Gauge (1x2), applies the gauge
%     rotation to the complex grid (c = e1 + i*e2), and scatters to full-mesh
%     order via GlobalVertices/GlobalFaces. Returns U=real, V=imag, N=cross(U,V).
%     Computes nothing new; persists nothing.
%
%     Vertex domain works for every gauge. Face domain works only for
%     euclidean/levi-civita: the trivial gauge's Gauge.face.rotation is deferred
%     in nxr, so face+trivial errors clearly.
%
% Authors: Diellor Basha, 2026

    Domain = 'vertex';
    for i = 1:2:numel(varargin)
        if strcmpi(varargin{i}, 'domain'); Domain = lower(varargin{i+1}); end
    end
    if ~ismember(Domain, {'vertex','face'})
        error('tess_frame:badDomain', 'Domain must be ''vertex'' or ''face''.');
    end

    TessMat = in_tess_bst(SurfaceFile, 0);
    if ~isfield(TessMat,'Geometry') || isempty(TessMat.Geometry) ...
       || ~isfield(TessMat,'Gauge') || isempty(TessMat.Gauge)
        error('tess_frame:noBundle', 'Surface has no stored bundle; run tess_bundle first.');
    end
    Geo = TessMat.Geometry; Ga = TessMat.Gauge;

    if strcmp(Domain,'vertex')
        nElem = size(TessMat.Vertices,1);
    else
        nElem = size(TessMat.Faces,1);
    end
    U = zeros(nElem,3); V = zeros(nElem,3);

    for hh = 1:numel(Geo)
        gridH = Geo(hh).(Domain).grid;          % nElemH x 3 complex
        if strcmpi(Ga(hh).type, 'trivial')
            if strcmp(Domain,'face')
                % nxr ships Gauge.face.rotation as a deferred (empty) placeholder.
                % Future-proof: support it the moment the build populates it.
                if ~isfield(Ga(hh).face,'rotation') || isempty(Ga(hh).face.rotation)
                    error('tess_frame:faceTrivialDeferred', ...
                        'Face-domain trivial frame needs Gauge.face.rotation (empty/deferred in nxr).');
                end
                rot = Ga(hh).face.rotation;      % nElemH x 1 complex (when populated)
            else
                rot = Ga(hh).vertex.rotation;    % nElemH x 1 complex
            end
        else
            rot = ones(size(gridH,1),1);
        end
        cRot = gridH .* rot;                      % broadcast over the 3 columns
        if strcmp(Domain,'vertex')
            idx = Geo(hh).GlobalVertices;
        else
            idx = Geo(hh).GlobalFaces;
        end
        U(idx,:) = real(cRot);
        V(idx,:) = imag(cRot);
    end

    N = cross(U, V, 2);
end
```

- [ ] **Step 4: Run to verify pass**

```matlab
rehash; clear test_tess_frame tess_frame; disp('rehashed');
runtests('dev/tests/test_tess_frame.m')
```
Expected: all 3 tests PASS.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_frame.m dev/tests/test_tess_frame.m
git commit -m "feat(tess-frame): derived {U,V,N} from stored bundle (vertex domain)"
```

---

## Final verification

- [ ] **Run both test files together; expect 12 passing.**

```matlab
rehash; clear test_tess_bundle test_tess_frame tess_bundle tess_frame; disp('rehashed');
results = [runtests('dev/tests/test_tess_bundle.m'), runtests('dev/tests/test_tess_frame.m')];
disp(table([results.Passed]', [results.Failed]', 'VariableNames', {'Passed','Failed'}));
assert(all([results.Passed]) && ~any([results.Failed]), 'Some tests failed.');
```

- [ ] **Then complete via superpowers:finishing-a-development-branch.**

---

## Notes for the implementer

- **MATLAB function caching:** in the persistent MCP session, new/edited functions and newly-added local test functions are not reliably auto-picked-up. Before each `runtests`, run `rehash; clear <specific function names>;`. **Never** a bare `clear` (it wipes Brainstorm's `GlobalData`); only `clear <name1> <name2> …`.
- **nxr field-name contract:** the tests read `B.Geometry(h).vertex.grid`, `B.Gauge(h).type`, `B.Gauge(h).vertex.rotation`, `B.Gauge(h).singularity.indices`, `B.Topology(h).operators.dec.{d0,d1}`, `B.Gauge(h).operators.{laplacian,covariantLaplacian}`. These match nxr's documented `bundle` schema. If a sub-field name differs in the freshly built MEX, fix the **test's** field path to the real name (the writer stores the struct verbatim, so it does not need to change) and note the discrepancy.
- **Struct-array assignment** (`TopoArr(hh) = sT`) requires both hemispheres' elements to carry identical field names — they do, because both flow through the same nxr call plus the same three appended maps.
- **Isolation:** every test that writes to disk backs up the full `.mat` and restores it via `onCleanup`. Do not weaken this to `NoSave` on the storage tests — they must exercise the real save path.
