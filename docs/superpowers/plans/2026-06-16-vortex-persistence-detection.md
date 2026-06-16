# Vortex-core Persistence Detection — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the "all 1-ring extrema + amplitude gate" vortex-core detector with topological-persistence-ranked detection, and rewire the Helmholtz GUI gate to persistence.

**Architecture:** A new geometry-free function `bst_vortex_persistence` computes 0-D persistence of the superlevel-set filtration of a scalar field (union-find / elder rule), separately for `+field` and `−field`. `bst_dirac_helmholtz` calls it **per hemisphere** (using cached per-hemisphere 1-ring adjacency + vertex coords), localizes each core to sub-vertex precision, and assembles a back-compatible struct array carrying `persistence`. `view_helmholtz`/`panel_helmholtz` gate markers by persistence instead of `|ω|`.

**Tech Stack:** MATLAB (Brainstorm). Tests are plain functions with a `chk(label,cond)` failure counter, run via the MATLAB MCP (`run_matlab_file`) or `matlab -batch`.

**Commits:** This repo's `development` branch is user-managed — commit steps below are OPTIONAL and gated by the user. Stage logically, but do not push/commit unless the user asks.

**Preconditions for running tests:** MATLAB up, Brainstorm started, `TutorialAuditory` protocol loaded (Subject01 cortex = `Surface(5)`, 20484V, with Dirac + Laplace-Beltrami operators).

---

### Task 1: `bst_vortex_persistence` — topological persistence detector

**Files:**
- Create: `toolbox/math/bst_vortex_persistence.m`
- Test: `dev/tests/test_vortex_persistence.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_vortex_persistence.m`:

```matlab
function test_vortex_persistence()
% Pure unit tests for bst_vortex_persistence on hand-computed tiny graphs.
% Author: Diellor Basha, 2026
    nFail = 0;

    % Chain 1-2-3-4-5, field = [5 2 4 1 3].
    % +field: global max v1 (Inf); v3 dies at saddle v2 -> pers 4-2=2;
    %         v5 dies at saddle v4 -> pers 3-1=2.
    % -field: global max v4 (Inf, min of field); v2 dies at saddle v3 -> pers 2.
    nb = {2; [1;3]; [2;4]; [3;5]; 4};
    f  = [5;2;4;1;3];
    C  = bst_vortex_persistence(f, [], 'Neighbors', nb);

    nFail = nFail + chk('5 features total', numel(C.vertex) == 5);
    nFail = nFail + chk('v1 is +global',  any(C.vertex==1 & C.chirality==1 & isinf(C.persistence)));
    nFail = nFail + chk('v4 is -global',  any(C.vertex==4 & C.chirality==-1 & isinf(C.persistence)));
    nFail = nFail + chk('v3 +core pers=2', any(C.vertex==3 & C.chirality==1  & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('v5 +core pers=2', any(C.vertex==5 & C.chirality==1  & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('v2 -core pers=2', any(C.vertex==2 & C.chirality==-1 & abs(C.persistence-2)<1e-12));
    nFail = nFail + chk('sorted by persistence desc', issorted(C.persistence,'descend'));

    % MinPersistence=3 keeps only the two Inf globals.
    C2 = bst_vortex_persistence(f, [], 'Neighbors', nb, 'MinPersistence', 3);
    nFail = nFail + chk('MinPersistence keeps only globals', numel(C2.vertex)==2 && all(isinf(C2.persistence)));

    % Flat field -> one global per sign, no spurious finite-persistence cores.
    Cf = bst_vortex_persistence(zeros(5,1), [], 'Neighbors', nb);
    nFail = nFail + chk('flat field -> 2 globals only', numel(Cf.vertex)==2 && all(isinf(Cf.persistence)));

    % Faces path builds 1-ring: single triangle, peak at the highest vertex.
    Ct = bst_vortex_persistence([3;1;2], [1 2 3]);
    nFail = nFail + chk('triangle global is vertex 1', any(Ct.vertex==1 & isinf(Ct.persistence) & Ct.chirality==1));

    fprintf('\n==== test_vortex_persistence: %d failed ====\n', nFail);
    if nFail > 0, error('test_vortex_persistence FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run the test, verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_vortex_persistence.m`
Expected: error `Undefined function 'bst_vortex_persistence'` (or similar).

- [ ] **Step 3: Implement `bst_vortex_persistence.m`**

Create `toolbox/math/bst_vortex_persistence.m`:

```matlab
function Cores = bst_vortex_persistence(field, faces, varargin)
% BST_VORTEX_PERSISTENCE  Vortex cores as persistent extrema of a scalar field.
%
% Cores are local extrema of FIELD on a triangular mesh, ranked by topological
% persistence (peak-minus-saddle contrast) from one union-find pass over the
% superlevel-set filtration (the elder rule). Maxima of +FIELD are positive-
% chirality cores; maxima of -FIELD (minima of FIELD) are negative-chirality.
% Every connected component contributes one never-merging (Inf-persistence)
% global core, so the result is correct on disconnected meshes (e.g. 2 hemispheres).
%
% USAGE:
%   Cores = bst_vortex_persistence(field, faces)
%   Cores = bst_vortex_persistence(field, faces, 'MinPersistence', p, 'Neighbors', nb)
%
% INPUT:
%   field  : [nV x 1] scalar field on mesh vertices (e.g. the stream function psi)
%   faces  : [nF x 3] 1-based triangles (used only to build the 1-ring if Neighbors
%            is not supplied; pass [] when supplying Neighbors)
%   'MinPersistence' : keep cores with persistence >= this (default 0; Inf always kept)
%   'Neighbors'      : {nV x 1} precomputed 1-ring adjacency (reuse across frames)
%
% OUTPUT: Cores (columnar struct, sorted by persistence descending):
%   .vertex .value .persistence .chirality(+1/-1) .isGlobal .birth .death
%
% Author: Diellor Basha, 2026

    field = field(:);
    nV = numel(field);

    MinPersistence = 0;  neighbors = {};
    for k = 1:2:numel(varargin)
        switch lower(varargin{k})
            case 'minpersistence', MinPersistence = varargin{k+1};
            case 'neighbors',      neighbors      = varargin{k+1};
            otherwise, error('bst_vortex_persistence: unknown option %s', varargin{k});
        end
    end
    if isempty(neighbors)
        neighbors = i_one_ring(faces, nV);
    elseif numel(neighbors) ~= nV
        error('bst_vortex_persistence: Neighbors must have one cell per vertex (%d vs %d).', ...
              numel(neighbors), nV);
    end

    pos = i_superlevel(field,  neighbors, nV);     % +field
    neg = i_superlevel(-field, neighbors, nV);     % -field

    vertex      = [pos.peak;        neg.peak];
    persistence = [pos.persistence; neg.persistence];
    birth       = [pos.birth;       neg.birth];
    death       = [pos.death;       neg.death];
    chirality   = [ ones(numel(pos.peak),1); -ones(numel(neg.peak),1)];
    isGlobal    = isinf(persistence);
    value       = field(vertex);

    keep = persistence >= MinPersistence;          % Inf always passes
    idx  = find(keep);
    [persistence, ord] = sort(persistence(idx), 'descend');
    idx  = idx(ord);

    Cores = struct('vertex',vertex(idx), 'value',value(idx), 'persistence',persistence, ...
                   'chirality',chirality(idx), 'isGlobal',isGlobal(idx), ...
                   'birth',birth(idx), 'death',death(idx));
end

% ---- 0-D persistence of the superlevel-set filtration of FIELD ----
function F = i_superlevel(field, neighbors, nV)
    [~, order] = sort(field, 'descend');
    parent     = zeros(nV,1);     % union-find parent; 0 = not yet added
    peakVertex = zeros(nV,1);     % peak vertex carried by each component root
    isAdded    = false(nV,1);
    featPeak   = zeros(nV,1);
    featDeath  = zeros(nV,1);
    nFeat      = 0;

    for k = 1:nV
        v = order(k);  isAdded(v) = true;
        nb = neighbors{v};  nb = nb(isAdded(nb));
        if isempty(nb)
            parent(v) = v;  peakVertex(v) = v;            % new component born (a maximum)
        else
            roots = unique(arrayfun(@findRoot, nb(:)));
            [~, hi]  = max(field(peakVertex(roots)));      % elder = highest peak
            survivor = roots(hi);
            survivorPeak = peakVertex(survivor);
            parent(v) = survivor;
            for r = roots(:)'
                if r ~= survivor
                    nFeat = nFeat + 1;
                    featPeak(nFeat)  = peakVertex(r);       % younger dies here
                    featDeath(nFeat) = field(v);            % at saddle value
                    parent(r) = survivor;
                end
            end
            peakVertex(survivor) = survivorPeak;
        end
    end
    % every surviving root never merges -> a global (Inf-persistence) feature
    for v = 1:nV
        if parent(v) == v
            nFeat = nFeat + 1;
            featPeak(nFeat)  = peakVertex(v);
            featDeath(nFeat) = -Inf;
        end
    end

    featPeak  = featPeak(1:nFeat);
    featDeath = featDeath(1:nFeat);
    F.peak        = featPeak;
    F.birth       = field(featPeak);
    F.death       = featDeath;
    F.persistence = F.birth - featDeath;                    % Inf for globals

    function r = findRoot(x)
        r = x;
        while parent(r) ~= r, r = parent(r); end
        while parent(x) ~= r, nx = parent(x); parent(x) = r; x = nx; end
    end
end

% ---- 1-ring vertex adjacency from triangle faces ----
function neighbors = i_one_ring(faces, nV)
    if isempty(faces)
        error('bst_vortex_persistence: need faces or Neighbors.');
    end
    e = [faces(:,[1 2]); faces(:,[2 3]); faces(:,[3 1])];
    A = sparse([e(:,1);e(:,2)], [e(:,2);e(:,1)], true, nV, nV);
    neighbors = cell(nV,1);
    for v = 1:nV, neighbors{v} = find(A(:,v)); end
end
```

- [ ] **Step 4: Run the test, verify it passes**

Run: `dev/tests/test_vortex_persistence.m`
Expected: `==== test_vortex_persistence: 0 failed ====`

- [ ] **Step 5: Lint**

Run MCP `check_matlab_code` on `toolbox/math/bst_vortex_persistence.m`; resolve anything beyond benign Brainstorm idioms.

- [ ] **Step 6: Commit (optional, user-gated)**

```bash
git add toolbox/math/bst_vortex_persistence.m dev/tests/test_vortex_persistence.m
git commit -m "feat(vortex): topological-persistence core detector"
```

---

### Task 2: Cache per-hemisphere adjacency + geometry in `Prepare`

**Files:**
- Modify: `toolbox/math/bst_dirac_helmholtz.m` (`Prepare`, ~lines 55–96)

- [ ] **Step 1: Add cache fields to the `Op` allocation**

In `Prepare`, change the `deal` line (currently allocating `Op.D … Op.Gz`) to also allocate the new per-hemisphere cells. Replace:

```matlab
    [Op.D, Op.vH, Op.Nf, Op.Wfv, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz] = deal(cell(1,nH));
```
with:
```matlab
    [Op.D, Op.vH, Op.Nf, Op.Wfv, Op.M, Op.cholK, Op.free, Op.totMass, ...
     Op.Gx, Op.Gy, Op.Gz, Op.NbH, Op.VtxH, Op.VnH] = deal(cell(1,nH));
    Op.Vtx = Vtx;
```

- [ ] **Step 2: Populate the cache inside the hemisphere loop**

Inside the `for hh = 1:nH` loop, after `Floc`/`Vloc` are computed (just after line 65 `Floc = mapV(...); Vloc = Vtx(vH, :);`), add:

```matlab
        % per-hemisphere 1-ring (LOCAL indexing) + geometry for core detection
        eLoc = [Floc(:,[1 2]); Floc(:,[2 3]); Floc(:,[3 1])];
        Aloc = sparse([eLoc(:,1);eLoc(:,2)], [eLoc(:,2);eLoc(:,1)], true, nVh, nVh);
        nbLoc = cell(nVh,1);
        for vv = 1:nVh, nbLoc{vv} = find(Aloc(:,vv)); end
        Op.NbH{hh}  = nbLoc;
        Op.VtxH{hh} = Vloc;
        Op.VnH{hh}  = Surf.VertNormals(vH, :);
```

- [ ] **Step 3: Smoke-check Prepare still runs**

Run (MCP `evaluate_matlab_code`):
```matlab
SurfaceFile = bst_get('Subject',1).Surface(5).FileName; Surf = in_tess_bst(SurfaceFile,0);
D = tess_operators(SurfaceFile,'Dirac'); L = tess_operators(SurfaceFile,'Laplace-Beltrami');
Op = bst_dirac_helmholtz('Prepare', D, L, Surf);
fprintf('NbH hemis=%d, VtxH{1}=%dx%d, has Vtx=%d\n', numel(Op.NbH), size(Op.VtxH{1},1), size(Op.VtxH{1},2), isfield(Op,'Vtx'));
```
Expected: `NbH hemis=2, VtxH{1}=Nx3, has Vtx=1` (N = hemisphere vertex count).

*Note:* `tess_operators(SurfaceFile, variant)` returns the operator struct; if it returns empty, load via the `i_load_op` helper already in `test_dirac_helmholtz.m`.

- [ ] **Step 4: Commit (optional, user-gated)**

```bash
git add toolbox/math/bst_dirac_helmholtz.m
git commit -m "feat(helmholtz): cache per-hemisphere 1-ring + geometry in Prepare"
```

---

### Task 3: Rewrite `FindCores` (per-hemisphere persistence + sub-vertex), preserve charge convention

**Files:**
- Modify: `toolbox/math/bst_dirac_helmholtz.m` (`Frame` lines 136–137; `FindCores` lines 171–192; add helpers)
- Test: `dev/tests/test_dirac_helmholtz.m` (extend)

- [ ] **Step 1: Write the failing test additions**

Append these checks to `test_dirac_helmholtz.m` inside the function, just before the final `fprintf('\n==== ...')` line:

```matlab
    % --- (6) persistence schema + per-hemisphere globals ---
    Op6 = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    % plant one positive psi bump per hemisphere; check each hemisphere yields a global core
    vL = Op6.vH{1}(round(numel(Op6.vH{1})/2));
    vR = Op6.vH{2}(round(numel(Op6.vH{2})/2));
    psi6 = zeros(size(V,1),1);
    psi6 = psi6 + exp(-sum((V-V(vL,:)).^2,2)/(2*0.01^2));
    psi6 = psi6 + exp(-sum((V-V(vR,:)).^2,2)/(2*0.01^2));
    cr = bst_dirac_helmholtz('FindCoresOp', psi6, Op6, zeros(size(V,1),1));
    nFail = nFail + chk('cores carry persistence field', isfield(cr,'persistence') && isfield(cr,'pos') && isfield(cr,'isGlobal'));
    glob = cr(logical([cr.isGlobal]));
    nFail = nFail + chk('one +global per hemisphere (>=2 globals)', sum([glob.chirality]==1) >= 2);
    nFail = nFail + chk('both planted bumps are cores', ismember(vL,[cr.iVertex]) && ismember(vR,[cr.iVertex]));
    nFail = nFail + chk('cores sorted by persistence desc', issorted([cr.persistence],'descend'));
    cL = cr([cr.iVertex]==vL);
    nFail = nFail + chk('sub-vertex pos near planted peak', ~isempty(cL) && norm(cL(1).pos - V(vL,:)) < 0.004);
```

Note this introduces a new public entry `'FindCoresOp'` (field + Op + omega → struct array) used by Frame and tests; the legacy `'FindCores'` (field, VertConn, omega) stays for back-compat.

- [ ] **Step 2: Run the test, verify it fails**

Run: `dev/tests/test_dirac_helmholtz.m`
Expected: FAIL — `Undefined ... 'FindCoresOp'` or missing fields.

- [ ] **Step 3: Rewrite `Frame`'s core calls and `FindCores`, add helpers**

In `Frame`, replace lines 136–137:
```matlab
    Ht.Cores    = FindCores(Ht.Psi, Op.VertConn, Ht.Curl);    % vortex cores (sign = vorticity)
    Ht.Sources  = FindCores(Ht.Phi, Op.VertConn, Ht.Div);     % sources/sinks (sign = divergence)
```
with:
```matlab
    Ht.Cores    = FindCoresOp(Ht.Psi, Op, Ht.Curl);           % vortex cores (sign = vorticity)
    Ht.Sources  = FindCoresOp(Ht.Phi, Op, Ht.Div);            % sources/sinks (sign = divergence)
```

Replace the whole legacy `FindCores` function (lines 171–192) with the following three functions:

```matlab
%% ===== core detection: persistence-ranked extrema, per hemisphere =====
function cores = FindCoresOp(field, Op, omega) %#ok<DEFNU>
% Per-hemisphere persistence detection + sub-vertex localization. Returns a struct
% array (sorted by persistence desc) carrying both the legacy fields (iVertex,
% charge, omega) and the new ones (persistence, isGlobal, birth, death, pos).
    cores = i_empty_cores();
    for hh = 1:numel(Op.vH)
        vH = Op.vH{hh};  nb = Op.NbH{hh};  Vloc = Op.VtxH{hh};  Vn = Op.VnH{hh};
        fl = field(vH);
        C  = bst_vortex_persistence(fl, [], 'Neighbors', nb);
        for k = 1:numel(C.vertex)
            vloc = C.vertex(k);  vg = vH(vloc);
            cores(end+1) = i_make_core(vg, omega(vg), C.chirality(k), C.persistence(k), ...
                C.isGlobal(k), C.birth(k), C.death(k), ...
                i_subvertex(vloc, fl, nb, Vloc, Vn)); %#ok<AGROW>
        end
    end
    if ~isempty(cores)
        [~, ord] = sort([cores.persistence], 'descend');
        cores = cores(ord);
    end
end

% Legacy entry (field + VertConn): runs on the full mesh; pos = NaN (no geometry).
function cores = FindCores(field, VertConn, omega) %#ok<DEFNU>
    nV = numel(field);
    [ii, jj] = find(VertConn);
    nb = accumarray(ii, jj, [nV 1], @(x){x}, {zeros(0,1)});
    C  = bst_vortex_persistence(field, [], 'Neighbors', nb);
    cores = i_empty_cores();
    for k = 1:numel(C.vertex)
        v = C.vertex(k);
        cores(end+1) = i_make_core(v, omega(v), C.chirality(k), C.persistence(k), ...
            C.isGlobal(k), C.birth(k), C.death(k), nan(1,3)); %#ok<AGROW>
    end
end

function s = i_empty_cores()
    s = struct('iVertex',{},'charge',{},'omega',{},'persistence',{}, ...
               'isGlobal',{},'birth',{},'death',{},'pos',{});
end

function s = i_make_core(vg, om, chirality, persistence, isGlobal, birth, death, pos)
    if om ~= 0, ch = sign(om); else, ch = -chirality; end   % preserve legacy fallback
    s = struct('iVertex',vg, 'charge',ch, 'omega',om, 'persistence',persistence, ...
               'isGlobal',logical(isGlobal), 'birth',birth, 'death',death, 'pos',pos);
end

% Sub-vertex localization: quadratic fit of FIELD over the 1-ring in a tangent chart.
function p = i_subvertex(vloc, field, nb, Vloc, Vn)
    v0 = Vloc(vloc,:);  p = v0;
    ns = nb{vloc};
    if numel(ns) < 5, return; end
    n = Vn(vloc,:);  n = n / max(norm(n), eps);
    e0 = [1 0 0]; if abs(e0*n') > 0.9, e0 = [0 1 0]; end
    t1 = e0 - (e0*n')*n; t1 = t1/norm(t1);  t2 = cross(n, t1);
    off = Vloc(ns,:) - v0;
    u = off*t1';  w = off*t2';  g = field(ns) - field(vloc);
    A = [u, w, 0.5*u.^2, u.*w, 0.5*w.^2];
    if rcond(A'*A) < 1e-10, return; end
    c = A \ g;                       % [b; c; d; e; f]
    H = [c(3) c(4); c(4) c(5)];
    if rcond(H) < 1e-8, return; end
    uw = -H \ [c(1); c(2)];
    if norm(uw) > max(sqrt(u.^2 + w.^2)), return; end   % reject runaway -> keep vertex
    p = v0 + uw(1)*t1 + uw(2)*t2;
end
```

Also update the USAGE comment block near the top (line ~28-29) to mention `FindCoresOp` alongside the legacy `FindCores`.

- [ ] **Step 4: Run the test, verify it passes**

Run: `dev/tests/test_dirac_helmholtz.m`
Expected: `==== test_dirac_helmholtz: 0 failed ====` (all original checks + the 5 new ones).

- [ ] **Step 5: Lint** `toolbox/math/bst_dirac_helmholtz.m` (MCP `check_matlab_code`).

- [ ] **Step 6: Commit (optional, user-gated)**

```bash
git add toolbox/math/bst_dirac_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "feat(helmholtz): persistence-ranked cores + sub-vertex localization"
```

---

### Task 4: Gate the GUI by persistence (`view_helmholtz`)

**Files:**
- Modify: `toolbox/gui/view_helmholtz.m` (`UpdateFrame`, lines 126–131; `i_readout`, lines 212–226)

- [ ] **Step 1: Replace the amplitude gate with a persistence gate**

Replace lines 126–131:
```matlab
    % --- component markers, pruned by the magnitude gate (fraction of frame max |omega|) ---
    mk = comp.Markers;
    if ~isempty(mk) && St.GateFrac > 0
        om = abs([mk.omega]);  mx = max(om);
        if mx > 0; mk = mk(om >= St.GateFrac * mx); end
    end
```
with:
```matlab
    % --- component markers, pruned by the persistence gate (fraction of frame max
    %     persistence); the never-merging global core of each hemisphere is always kept ---
    mk = comp.Markers;
    if ~isempty(mk) && St.GateFrac > 0
        pr  = [mk.persistence];
        mxf = max([pr(isfinite(pr)), 0]);
        if mxf > 0
            mk = mk(isinf(pr) | (pr >= St.GateFrac * mxf));
        end
    end
```

- [ ] **Step 2: Add top-persistence to the readout**

In `i_readout`, for the `'vortex'` and `'source'` cases, append the top persistence. Replace the `'vortex'` case body:
```matlab
        case 'vortex'
            if isempty(mk); txt = '0 vortices, 0 antivortices';
            else; np=sum([mk.charge]>0); nn=sum([mk.charge]<0); txt=sprintf('%d vortices (+), %d antivortices (-), net %+d', np, nn, np-nn); end
```
with:
```matlab
        case 'vortex'
            if isempty(mk); txt = '0 vortices, 0 antivortices';
            else
                np=sum([mk.charge]>0); nn=sum([mk.charge]<0);
                pr=[mk.persistence]; tp=max(pr(isfinite(pr)));
                if isempty(tp); tp=0; end
                txt=sprintf('%d vortices (+), %d antivortices (-), net %+d; top persistence %.2g', np, nn, np-nn, tp);
            end
```
(Leave the `'source'` and other cases unchanged.)

- [ ] **Step 3: Verify the GUI path constructs without error on a real frame**

Run (MCP `evaluate_matlab_code`) — exercises the exact `Frame` + gate logic the GUI uses, without opening a figure:
```matlab
HMos = in_bst_headmodel('Subject01/S01_AEF_20131218_01_notch/headmodel_surf_os_meg.mat',0);
Surf = in_tess_bst(HMos.SurfaceFile,0);
D = tess_operators(HMos.SurfaceFile,'Dirac'); L = tess_operators(HMos.SurfaceFile,'Laplace-Beltrami');
Op = bst_dirac_helmholtz('Prepare', D, L, Surf);
G = double(HMos.Gain); nV=size(Surf.Vertices,1);
% any nonzero unconstrained frame is fine for a structural check:
rng(0); Jt = randn(3*nV,1)*1e-9;
Ht = bst_dirac_helmholtz('Frame', Op, Jt);
mk = Ht.Cores; pr=[mk.persistence]; mxf=max(pr(isfinite(pr)));
kept50 = mk(isinf(pr) | (pr >= 0.5*mxf));
fprintf('cores=%d  globals=%d  kept@0.5=%d  hasPos=%d\n', numel(mk), sum([mk.isGlobal]), numel(kept50), all(arrayfun(@(c)numel(c.pos)==3,mk)));
```
Expected: `cores=… globals>=2 kept@0.5<cores hasPos=1` (gate reduces the count; ≥2 hemispheric globals; every core has a 1×3 pos).

- [ ] **Step 4: Commit (optional, user-gated)**

```bash
git add toolbox/gui/view_helmholtz.m
git commit -m "feat(helmholtz): gate vortex markers by persistence, not amplitude"
```

---

### Task 5: Relabel the gate slider (`panel_helmholtz`)

**Files:**
- Modify: `toolbox/gui/panel_helmholtz.m` (~line 46, the "Marker threshold (magnitude gate)" label; line 6 header comment)

- [ ] **Step 1: Update the label text**

At the gate-slider construction (~line 46), change the displayed label from the magnitude-gate wording to persistence. Find the `gui_component('label', ...)` (or equivalent) that reads "Marker threshold (magnitude gate)" / "magnitude" and change its string to `'Marker threshold (persistence)'`. Update the panel header comment on line 6 ("slider prunes weak singular points") to "slider prunes low-persistence singular points".

- [ ] **Step 2: Verify the panel builds**

Run (MCP `evaluate_matlab_code`):
```matlab
[bstPanel,~] = panel_helmholtz('CreatePanel', [], (1:50)');
disp(class(bstPanel));   % expect a BstPanel object, no error
```
Expected: prints a panel class name, no error.

- [ ] **Step 3: Commit (optional, user-gated)**

```bash
git add toolbox/gui/panel_helmholtz.m
git commit -m "ui(helmholtz): relabel gate slider as persistence"
```

---

### Task 6: Full regression + integration verification

**Files:** none (verification only)

- [ ] **Step 1: Run both test suites**

Run (MCP `run_matlab_file`): `dev/tests/test_vortex_persistence.m` then `dev/tests/test_dirac_helmholtz.m`.
Expected: both print `0 failed`.

- [ ] **Step 2: Integration on the real S01 alpha vortex**

Run (MCP `evaluate_matlab_code`):
```matlab
gui_brainstorm('SetCurrentProtocol', 2);   % TutorialAuditory (if not already)
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
[sS,~]=bst_get('DataFile',df); cm=in_bst_channel(sS.Channel(1).FileName); ty={cm.Channel.Type};
HMos=in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'],0);
G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(ty(:),'MEG');
NC=load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
HMf=HMos; HMf.Gain=G(iMEG,:);
OPT=struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','dspm2018');
OPT.NoiseCovMat.NoiseCov=Cn; OPT.ChannelTypes=ty(iMEG);
Rd=bst_inverse_dirac(HMf,OPT);
DM=in_bst_data(df); [~,iT]=min(abs(DM.Time-22.6)); Jt=Rd.ImagingKernel*double(DM.F(iMEG,iT));
Surf=in_tess_bst(HMos.SurfaceFile,0);
D=tess_operators(HMos.SurfaceFile,'Dirac'); L=tess_operators(HMos.SurfaceFile,'Laplace-Beltrami');
Op=bst_dirac_helmholtz('Prepare',D,L,Surf); Ht=bst_dirac_helmholtz('Frame',Op,Jt);
mk=Ht.Cores; pr=[mk.persistence]; mxf=max(pr(isfinite(pr)));
for gf=[0 .25 .5 .75]
    n=sum(isinf(pr)|(pr>=gf*mxf));
    fprintf('gate %.2f -> %d cores\n', gf, n);
end
top=mk(1); fprintf('top core: vertex %d @ [%.0f %.0f %.0f]mm pers=%.2g chir=%+d\n', top.iVertex, top.pos*1e3, top.persistence, top.charge);
```
Expected: core count **decreases monotonically** as the gate rises (e.g. ~34 → a handful); the top-persistence core sits in superior parietal (≈[-45 -12 74] mm region) — the known alpha vortex.

- [ ] **Step 3: Final lint sweep** of all three modified/created source files; confirm no new non-idiom warnings.

- [ ] **Step 4: Commit (optional, user-gated)**

```bash
git add -A
git commit -m "test(helmholtz): persistence detector regression + integration checks"
```

---

## Self-review notes

- **Spec coverage:** detector (Task 1) ✓; per-hemisphere + cache (Tasks 2–3) ✓; sub-vertex localization (Task 3, `i_subvertex`) ✓; output schema with `persistence/isGlobal/birth/death/pos` + preserved `iVertex/charge/omega` (Task 3, `i_make_core`) ✓; persistence gate (Task 4) ✓; panel relabel (Task 5) ✓; tests at unit/per-hemisphere/integration levels (Tasks 1,3,6) ✓.
- **Charge convention preserved:** `i_make_core` keeps `charge = sign(omega)`, falling back to `-chirality` when `omega==0` — so the existing test (2) (`charge<0` for a planted +ψ max with ω=0) stays green, and marker colors keep their meaning.
- **Naming consistency:** new public entry `FindCoresOp`; `i_make_core` / `i_empty_cores` / `i_subvertex` used consistently in Task 3; `Op.NbH/VtxH/VnH/Vtx` defined in Task 2 and consumed in Task 3.
- **Equivalence note:** because hemispheres are disconnected mesh components, "record every surviving root" inside `bst_vortex_persistence` yields one global per hemisphere; calling it per-hemisphere (Task 3) additionally reuses cached adjacency and is robust to any within-hemisphere disconnection. Both honor the per-hemisphere decision.
