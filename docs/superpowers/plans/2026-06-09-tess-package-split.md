# Topology / Geometry / Gauge Package Split Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Split the single `tess_bundle` writer into three per-package writers `tess_topology` / `tess_geometry` / `tess_gauge` (each mapping to nxr's standalone command), keeping `tess_bundle` as a light wrapper that produces all three in one nxr `bundle` pass.

**Architecture:** A shared private helper `tess_store_perhemi.m` owns the per-hemisphere orchestration common to all four writers — guards (nxr / atlas labels / connected-hemisphere / reg-sphere), hemisphere split, per-hemisphere submesh build + nxr handle lifecycle + FreeSurfer-pole computation, scatter-map attachment (`GlobalVertices`/`GlobalFaces`/`Hemisphere`/`Provenance`), `1×2` assembly, and save. Each writer is a thin shell: it parses its own options, builds a `Provenance` struct and a `ComputeFn` callback that calls the matching nxr command, and delegates to the helper. `tess_bundle` uses nxr's `bundle` command (one pass, all three fields); the three standalone writers use `topology`/`geometry`/`gauge`. `tess_frame.m` is unchanged.

**Tech Stack:** MATLAB (Brainstorm toolbox), `nxr-compute` MEX, `matlab.unittest`, run via the MATLAB MCP.

**Prerequisite:** the installed nxr-compute MEX supports the standalone commands `nxr_compute('topology', h[, opts])`, `('geometry', h[, opts])`, `('gauge', h, type[, opts])`, and `('bundle', h, type[, opts])` (all verified live to return `schemaVersion`+per-element groups, and `.operators` when `opts.operators=true`).

**Reference patterns (read first):**
- `toolbox/anatomy/tess_bundle.m` (current, merged) — the per-hemisphere loop, guards, scatter maps, provenance, save pattern being factored out.
- `dev/tests/test_tess_bundle.m` — the test conventions (`local_find_cortex(20484)`, `brainstorm nogui`, backup/restore via `onCleanup`).

**MATLAB-session discipline (all tasks):** before each test run, `rehash; clear <names>;` (e.g. `rehash; clear test_tess_topology tess_topology tess_store_perhemi;`). **Never** a bare `clear`/`clear all` (wipes Brainstorm `GlobalData`). Tests that exercise the save path must back up the `.mat` and restore via `onCleanup`.

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/anatomy/tess_store_perhemi.m` | **Create.** Shared private orchestration: guards, split, per-hemisphere nxr call via a `ComputeFn`, scatter maps, `1×2` assembly, save. |
| `toolbox/anatomy/tess_bundle.m` | **Modify.** Refactor to a thin wrapper over the helper using nxr `bundle` (behavior unchanged). |
| `toolbox/anatomy/tess_topology.m` | **Create.** Thin writer → `TessMat.Topology` via nxr `topology`. |
| `toolbox/anatomy/tess_geometry.m` | **Create.** Thin writer → `TessMat.Geometry` via nxr `geometry` (+`Mass`). |
| `toolbox/anatomy/tess_gauge.m` | **Create.** Thin writer → `TessMat.Gauge` via nxr `gauge` (+`Gauge` type, `Coupling`, `Mass`). |
| `dev/tests/test_tess_topology.m` | **Create.** |
| `dev/tests/test_tess_geometry.m` | **Create.** |
| `dev/tests/test_tess_gauge.m` | **Create.** |
| `dev/tests/test_tess_bundle.m` | **Unchanged** — runs as a regression check that the refactor preserved behavior. |

---

## Task 1: Shared helper + refactor `tess_bundle`

**Files:**
- Create: `toolbox/anatomy/tess_store_perhemi.m`
- Modify: `toolbox/anatomy/tess_bundle.m`
- Test: `dev/tests/test_tess_bundle.m` (unchanged; used as regression)

- [ ] **Step 1: Run the existing bundle tests to capture the green baseline**

```matlab
rehash; clear test_tess_bundle tess_bundle; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: 10 Passed (baseline to preserve through the refactor).

- [ ] **Step 2: Write the shared helper**

Create `toolbox/anatomy/tess_store_perhemi.m`:

```matlab
function B = tess_store_perhemi(SurfaceFile, FieldNames, NeedSphere, Provenance, UseCache, NoSave, ComputeFn)
% TESS_STORE_PERHEMI: Shared per-hemisphere nxr-bundle orchestration.
%
% Splits the cortex by hemisphere, runs ComputeFn on each hemisphere submesh,
% attaches scatter maps, assembles 1x2 struct arrays, and stores the requested
% fields on the surface file. Backs the tess_topology/geometry/gauge/bundle writers.
%
% INPUTS:
%   SurfaceFile  Brainstorm surface file name.
%   FieldNames   cellstr subset of {'Topology','Geometry','Gauge'} to store.
%   NeedSphere   logical; require a FreeSurfer reg sphere + compute poles (gauge).
%   Provenance   struct attached as .Provenance on every stored element.
%   UseCache     logical; if true and all FieldNames already present, return cached.
%   NoSave       logical; skip writing to disk.
%   ComputeFn    @(h, polesLocal) -> struct whose fields are (a superset of)
%                FieldNames, each the nxr struct for that hemisphere submesh.
%                polesLocal is [iN; iS] local indices (empty when ~NeedSphere).
%
% OUTPUT: B, struct with the FieldNames as 1x2 per-hemisphere struct arrays.
%
% Authors: Diellor Basha, 2026

    TessFile = file_fullpath(SurfaceFile);
    TessMat  = in_tess_bst(SurfaceFile, 0);

    % --- cache return: all requested fields already present ---
    if UseCache
        haveAll = true;
        for f = 1:numel(FieldNames)
            if ~isfield(TessMat, FieldNames{f}) || isempty(TessMat.(FieldNames{f}))
                haveAll = false; break;
            end
        end
        if haveAll
            B = struct();
            for f = 1:numel(FieldNames), B.(FieldNames{f}) = TessMat.(FieldNames{f}); end
            return;
        end
    end

    % --- require nxr-compute ---
    [isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
    if ~isOk
        error('tess_store_perhemi:nxrUnavailable', 'requires nxr-compute: %s', errMsg);
    end

    % --- require FreeSurfer registration sphere (gauge only) ---
    if NeedSphere
        if ~isfield(TessMat,'Reg') || ~isstruct(TessMat.Reg) || ~isfield(TessMat.Reg,'Sphere') ...
           || ~isfield(TessMat.Reg.Sphere,'Vertices') || isempty(TessMat.Reg.Sphere.Vertices)
            error('tess_store_perhemi:noRegSphere', ...
                'Trivial gauge needs a FreeSurfer registration sphere (Reg.Sphere.Vertices).');
        end
    end

    % --- require Structures atlas with L/R labels ---
    hasLabels = false;
    if isfield(TessMat,'Atlas') && ~isempty(TessMat.Atlas)
        iStruct = find(strcmpi({TessMat.Atlas.Name}, 'Structures'), 1);
        if ~isempty(iStruct) && ~isempty(TessMat.Atlas(iStruct).Scouts)
            scouts = TessMat.Atlas(iStruct).Scouts;
            labels = {scouts.Label}; regions = {scouts.Region};
            reg1 = cellfun(@(c) c(1), regions(~cellfun(@isempty, regions)), 'UniformOutput', false);
            hasL = any(strcmpi(labels,'lh')) || any(strcmpi(reg1,'L'));
            hasR = any(strcmpi(labels,'rh')) || any(strcmpi(reg1,'R'));
            hasLabels = hasL && hasR;
        end
    end
    if ~hasLabels
        error('tess_store_perhemi:noHemisphereLabels', ...
            'Surface has no Structures atlas with left/right hemisphere labels.');
    end

    % --- hemisphere split (import labels) ---
    [rH, lH, isConn] = tess_hemisplit(TessMat);
    if isConn
        error('tess_store_perhemi:connectedHemispheres', ...
            'Hemispheres are connected; nxr bundles each as an independent component.');
    end
    hemis = {lH(:), rH(:)}; tags = {'L','R'};
    Vtx = double(TessMat.Vertices); Fcs = double(TessMat.Faces); nVtot = size(Vtx,1);

    Arr = struct();
    for hh = 1:2
        vH = hemis{hh};
        if isempty(vH)
            error('tess_store_perhemi:emptyHemisphere', 'Hemisphere %s has no vertices.', tags{hh});
        end
        isV = false(nVtot,1); isV(vH) = true;
        fMask = all(isV(Fcs), 2);
        map = zeros(nVtot,1); map(vH) = 1:numel(vH);
        Vloc = Vtx(vH,:);
        Floc = map(Fcs(fMask,:));

        poles = [];
        if NeedSphere
            sph = TessMat.Reg.Sphere.Vertices(vH,:);
            [~, iN] = max(sph(:,3));
            [~, iSp] = min(sph(:,3));
            poles = [iN; iSp];
        end

        h = nxr_compute('create', Vloc, Floc);
        S = ComputeFn(h, poles);
        nxr_compute('destroy', h);

        for f = 1:numel(FieldNames)
            s = S.(FieldNames{f});
            s.GlobalVertices = vH;
            s.GlobalFaces    = find(fMask);
            s.Hemisphere     = tags{hh};
            s.Provenance     = Provenance;
            Arr.(FieldNames{f})(hh) = s;   %#ok<AGROW>
        end
    end

    B = struct();
    for f = 1:numel(FieldNames), B.(FieldNames{f}) = Arr.(FieldNames{f}); end

    % --- save ---
    if ~NoSave
        TessMat_full = load(TessFile);
        for f = 1:numel(FieldNames), TessMat_full.(FieldNames{f}) = B.(FieldNames{f}); end
        TessMat_full = bst_history('add', TessMat_full, 'bundle', ...
            sprintf('Stored %s (per-hemisphere nxr).', strjoin(FieldNames, '/')));
        bst_save(TessFile, TessMat_full, 'v7');
    end
end
```

- [ ] **Step 3: Refactor `tess_bundle.m` onto the helper**

Replace the entire body of `toolbox/anatomy/tess_bundle.m` (keep the leading help comment + license block; replace the option-parse-through-save implementation) with:

```matlab
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

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','bundle', 'NxrVersion',nxrVer, 'Gauge',Gauge, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_bundle(h, poles, Gauge, Operators, Coupling, Mass);
    B = tess_store_perhemi(SurfaceFile, {'Topology','Geometry','Gauge'}, ...
            strcmpi(Gauge,'trivial'), prov, ~ForceRecompute && ~Operators, NoSave, computeFn);
end

function S = local_compute_bundle(h, poles, Gauge, Operators, Coupling, Mass)
    opts = struct();
    if strcmpi(Gauge,'trivial')
        opts.singVerts = poles; opts.singValues = [1; 1];
    end
    if Operators
        opts.operators = true; opts.coupling = Coupling; opts.mass = Mass;
    end
    S = nxr_compute('bundle', h, Gauge, opts);   % returns {Topology,Geometry,Gauge}
```

(The function signature line `function B = tess_bundle(SurfaceFile, varargin)` and the help/license header stay; everything from the old `% --- options ---` onward is replaced by the above, and the new `local_compute_bundle` subfunction is appended.)

- [ ] **Step 4: Run the bundle regression tests**

```matlab
rehash; clear test_tess_bundle tess_bundle tess_store_perhemi; disp('rehashed');
runtests('dev/tests/test_tess_bundle.m')
```
Expected: 10 Passed (identical to the Step-1 baseline). If `test_three_struct_arrays_1x2` fails on `Provenance.Backend`, confirm the helper attaches `Provenance` (it does) — the value is still `'nxr'`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_store_perhemi.m toolbox/anatomy/tess_bundle.m
git commit -m "refactor(tess-bundle): extract tess_store_perhemi shared helper"
```

---

## Task 2: `tess_topology`

**Files:**
- Create: `toolbox/anatomy/tess_topology.m`
- Test: `dev/tests/test_tess_topology.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_topology.m`:

```matlab
function tests = test_tess_topology
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

function test_topology_struct_array(tc)
    Topo = tess_topology(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(Topo), [1 2]);
    verifyTrue(tc, isfield(Topo(1), 'schemaVersion'));
    verifyEqual(tc, Topo(1).Hemisphere, 'L');
    verifyEqual(tc, Topo(2).Hemisphere, 'R');
    verifyEqual(tc, Topo(1).Provenance.Package, 'topology');
end

function test_topology_partition(tc)
    SurfaceFile = local_cortex();
    Topo = tess_topology(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    gv = sort([Topo(1).GlobalVertices(:); Topo(2).GlobalVertices(:)]);
    verifyEqual(tc, gv, (1:size(T.Vertices,1))');
end

function test_topology_operators(tc)
    Topo = tess_topology(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1);
    for hh = 1:2
        op = Topo(hh).operators;
        verifyEqual(tc, nnz(op.dec.d1 * op.dec.d0), 0);   % d o d = 0 (exact)
    end
end

function test_topology_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_topology(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Topology') && isequal(size(T.Topology),[1 2]));
end
```

- [ ] **Step 2: Run to verify failure**

```matlab
rehash; clear test_tess_topology tess_topology; disp('rehashed');
runtests('dev/tests/test_tess_topology.m')
```
Expected: ERROR `Undefined function 'tess_topology'`.

- [ ] **Step 3: Write `tess_topology.m`**

Create `toolbox/anatomy/tess_topology.m`:

```matlab
function Topology = tess_topology(SurfaceFile, varargin)
% TESS_TOPOLOGY: Store the nxr-compute halfedge Topology package per hemisphere.
%
% USAGE:  Topology = tess_topology(SurfaceFile)
%         Topology = tess_topology(SurfaceFile, 'Operators',1, 'NoSave',1, 'ForceRecompute',1)
%
% Stores TessMat.Topology as a 1x2 per-hemisphere struct array ((1)=left,(2)=right).
% With 'Operators',1 the nxr graph laplacian + DEC operators are attached.
%
% Authors: Diellor Basha, 2026

    Operators=false; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'operators',      Operators=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','topology', 'NxrVersion',nxrVer, ...
                  'Operators',logical(Operators), 'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_topology(h, Operators);
    B = tess_store_perhemi(SurfaceFile, {'Topology'}, false, prov, ...
            ~ForceRecompute && ~Operators, NoSave, computeFn);
    Topology = B.Topology;
end

function S = local_compute_topology(h, Operators)
    opts = struct();
    if Operators, opts.operators = true; end
    S = struct('Topology', nxr_compute('topology', h, opts));
end
```

- [ ] **Step 4: Run to verify pass**

```matlab
rehash; clear test_tess_topology tess_topology tess_store_perhemi; disp('rehashed');
runtests('dev/tests/test_tess_topology.m')
```
Expected: 4 Passed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_topology.m dev/tests/test_tess_topology.m
git commit -m "feat(tess-topology): per-hemisphere nxr topology writer"
```

---

## Task 3: `tess_geometry`

**Files:**
- Create: `toolbox/anatomy/tess_geometry.m`
- Test: `dev/tests/test_tess_geometry.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_geometry.m` (copy the `local_cortex`/`local_find_cortex` helpers verbatim from `test_tess_topology.m`, then):

```matlab
function test_geometry_struct_array(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    verifyEqual(tc, size(Geo), [1 2]);
    verifyEqual(tc, Geo(1).Provenance.Package, 'geometry');
end

function test_geometry_grid_orthonormal(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        c = Geo(hh).vertex.grid; e1 = real(c); e2 = imag(c);
        verifyLessThan(tc, max(abs(sqrt(sum(e1.^2,2))-1)), 1e-4);
        verifyLessThan(tc, max(abs(sqrt(sum(e2.^2,2))-1)), 1e-4);
        verifyLessThan(tc, max(abs(sum(e1.*e2,2))), 1e-4);
    end
end

function test_geometry_operators(tc)
    Geo = tess_geometry(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1, 'Mass', 'lumped');
    for hh = 1:2
        op = Geo(hh).operators;
        verifyTrue(tc, isfield(op,'laplacian') && isfield(op,'mass') && isfield(op,'hodge'));
        verifyTrue(tc, issparse(op.mass.lumped));
    end
end

function test_geometry_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_geometry(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Geometry') && isequal(size(T.Geometry),[1 2]));
end
```

(Header line: `function tests = test_tess_geometry` + `tests = functiontests(localfunctions); end`, same shape as `test_tess_topology.m`.)

- [ ] **Step 2: Run to verify failure**

```matlab
rehash; clear test_tess_geometry tess_geometry; disp('rehashed');
runtests('dev/tests/test_tess_geometry.m')
```
Expected: ERROR `Undefined function 'tess_geometry'`.

- [ ] **Step 3: Write `tess_geometry.m`**

Create `toolbox/anatomy/tess_geometry.m`:

```matlab
function Geometry = tess_geometry(SurfaceFile, varargin)
% TESS_GEOMETRY: Store the nxr-compute Geometry package per hemisphere.
%
% USAGE:  Geometry = tess_geometry(SurfaceFile)
%         Geometry = tess_geometry(SurfaceFile, 'Operators',1, 'Mass','lumped', ...)
%
% Stores TessMat.Geometry as a 1x2 per-hemisphere struct array. With 'Operators',1
% the nxr cotan laplacian + mass + hodge operators are attached.
%
% Authors: Diellor Basha, 2026

    Operators=false; Mass='lumped'; NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'operators',      Operators=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','geometry', 'NxrVersion',nxrVer, ...
                  'Operators',logical(Operators), 'Mass',Mass, 'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_geometry(h, Operators, Mass);
    B = tess_store_perhemi(SurfaceFile, {'Geometry'}, false, prov, ...
            ~ForceRecompute && ~Operators, NoSave, computeFn);
    Geometry = B.Geometry;
end

function S = local_compute_geometry(h, Operators, Mass)
    opts = struct();
    if Operators, opts.operators = true; opts.mass = Mass; end
    S = struct('Geometry', nxr_compute('geometry', h, opts));
end
```

- [ ] **Step 4: Run to verify pass**

```matlab
rehash; clear test_tess_geometry tess_geometry tess_store_perhemi; disp('rehashed');
runtests('dev/tests/test_tess_geometry.m')
```
Expected: 4 Passed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_geometry.m dev/tests/test_tess_geometry.m
git commit -m "feat(tess-geometry): per-hemisphere nxr geometry writer"
```

---

## Task 4: `tess_gauge`

**Files:**
- Create: `toolbox/anatomy/tess_gauge.m`
- Test: `dev/tests/test_tess_gauge.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_tess_gauge.m` (copy the `local_cortex`/`local_find_cortex` helpers verbatim from `test_tess_topology.m`, then):

```matlab
function test_gauge_struct_array(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);   % default trivial
    verifyEqual(tc, size(Ga), [1 2]);
    verifyEqual(tc, lower(Ga(1).type), 'trivial');
    verifyEqual(tc, Ga(1).Provenance.Package, 'gauge');
end

function test_gauge_gauss_bonnet(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1);
    for hh = 1:2
        verifyEqual(tc, sum(Ga(hh).singularity.indices), 2, 'AbsTol', 1e-6);
    end
end

function test_gauge_operators(tc)
    Ga = tess_gauge(local_cortex(), 'NoSave', 1, 'ForceRecompute', 1, 'Operators', 1, 'Coupling', 'ambient');
    for hh = 1:2
        nVh = numel(Ga(hh).GlobalVertices);
        K = Ga(hh).operators.laplacian;
        verifyTrue(tc, issparse(K));
        verifyGreaterThan(tc, full(max(abs(imag(K)), [], 'all')), 0);     % genuinely complex
        verifyLessThan(tc, full(max(abs(K - K'), [], 'all')), 1e-6);      % Hermitian
        L3 = Ga(hh).operators.covariantLaplacian;
        verifyEqual(tc, size(L3), [3*nVh, 3*nVh]);
    end
end

function test_gauge_is_stored(tc)
    SurfaceFile = local_cortex();
    TessFile = file_fullpath(SurfaceFile);
    backup = load(TessFile);
    restorer = onCleanup(@() bst_save(TessFile, backup, 'v7'));  %#ok<NASGU>
    tess_gauge(SurfaceFile, 'ForceRecompute', 1);
    T = in_tess_bst(SurfaceFile, 0);
    verifyTrue(tc, isfield(T,'Gauge') && isequal(size(T.Gauge),[1 2]));
end
```

(Header line: `function tests = test_tess_gauge` + `tests = functiontests(localfunctions); end`.)

- [ ] **Step 2: Run to verify failure**

```matlab
rehash; clear test_tess_gauge tess_gauge; disp('rehashed');
runtests('dev/tests/test_tess_gauge.m')
```
Expected: ERROR `Undefined function 'tess_gauge'`.

- [ ] **Step 3: Write `tess_gauge.m`**

Create `toolbox/anatomy/tess_gauge.m`:

```matlab
function Gauge = tess_gauge(SurfaceFile, varargin)
% TESS_GAUGE: Store the nxr-compute Gauge package per hemisphere.
%
% USAGE:  Gauge = tess_gauge(SurfaceFile)
%         Gauge = tess_gauge(SurfaceFile, 'Gauge','trivial', 'Operators',1, 'Coupling','ambient', 'Mass','lumped')
%
% Stores TessMat.Gauge as a 1x2 per-hemisphere struct array. The trivial gauge
% places singularities at the FreeSurfer sphere poles. With 'Operators',1 the
% nxr connection Laplacian + covariant Laplacian are attached.
%
% Authors: Diellor Basha, 2026

    GaugeType='trivial'; Operators=false; Coupling='ambient'; Mass='lumped';
    NoSave=false; ForceRecompute=false;
    for i = 1:2:numel(varargin)
        switch lower(varargin{i})
            case 'gauge',          GaugeType=varargin{i+1};
            case 'operators',      Operators=varargin{i+1};
            case 'coupling',       Coupling=varargin{i+1};
            case 'mass',           Mass=varargin{i+1};
            case 'nosave',         NoSave=varargin{i+1};
            case 'forcerecompute', ForceRecompute=varargin{i+1};
        end
    end

    nxrVer = '';
    try, nxrVer = nxr_compute('version'); catch, end   %#ok<CTCH>
    prov = struct('Backend','nxr', 'Package','gauge', 'NxrVersion',nxrVer, 'Gauge',GaugeType, ...
                  'Operators',logical(Operators), 'Coupling',Coupling, 'Mass',Mass, ...
                  'ComputeDate',datestr(now,'yyyy-mm-dd HH:MM:SS'));

    computeFn = @(h, poles) local_compute_gauge(h, poles, GaugeType, Operators, Coupling, Mass);
    B = tess_store_perhemi(SurfaceFile, {'Gauge'}, strcmpi(GaugeType,'trivial'), prov, ...
            ~ForceRecompute && ~Operators, NoSave, computeFn);
    Gauge = B.Gauge;
end

function S = local_compute_gauge(h, poles, GaugeType, Operators, Coupling, Mass)
    opts = struct();
    if strcmpi(GaugeType,'trivial')
        opts.singVerts = poles; opts.singValues = [1; 1];
    end
    if Operators
        opts.operators = true; opts.coupling = Coupling; opts.mass = Mass;
    end
    S = struct('Gauge', nxr_compute('gauge', h, GaugeType, opts));
end
```

- [ ] **Step 4: Run to verify pass**

```matlab
rehash; clear test_tess_gauge tess_gauge tess_store_perhemi; disp('rehashed');
runtests('dev/tests/test_tess_gauge.m')
```
Expected: 4 Passed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/tess_gauge.m dev/tests/test_tess_gauge.m
git commit -m "feat(tess-gauge): per-hemisphere nxr gauge writer"
```

---

## Final verification

- [ ] **Run all four test files; expect 22 passing (10 bundle + 4 + 4 + 4).**

```matlab
rehash; clear test_tess_bundle test_tess_topology test_tess_geometry test_tess_gauge ...
      tess_bundle tess_topology tess_geometry tess_gauge tess_store_perhemi tess_frame; disp('rehashed');
r = [ runtests('dev/tests/test_tess_bundle.m'), ...
      runtests('dev/tests/test_tess_topology.m'), ...
      runtests('dev/tests/test_tess_geometry.m'), ...
      runtests('dev/tests/test_tess_gauge.m') ];
fprintf('\nTOTAL: %d passed, %d failed (of %d)\n', sum([r.Passed]), sum([r.Failed]), numel(r));
assert(all([r.Passed]) && ~any([r.Failed]), 'Some tests failed.');
```

- [ ] **Then complete via superpowers:finishing-a-development-branch.**

---

## Notes for the implementer

- The four writers are intentionally near-identical thin shells over `tess_store_perhemi`; the only differences are (a) which nxr command the `ComputeFn` calls, (b) which options each accepts, and (c) the `FieldNames`/`NeedSphere` passed to the helper. Do not add logic to the writers beyond option parsing + provenance + the `ComputeFn`.
- `tess_bundle` must remain behaviorally identical — Task 1's success criterion is the **unchanged** `test_tess_bundle.m` passing 10/10 after the refactor.
- nxr command opts: `nxr_compute('topology'/'geometry', h, opts)` and `nxr_compute('gauge', h, type, opts)` accept an empty `struct()` for the light path (verified). If a light call rejects an empty struct, omit the third/fourth arg when no options are set and note it.
- Provenance now carries a `Package` field (`'topology'/'geometry'/'gauge'/'bundle'`) so a stored field records which writer produced it.
- This refactor is independent of the in-flight nxr connection-Laplacian PSD fix; it does not touch operator math.
