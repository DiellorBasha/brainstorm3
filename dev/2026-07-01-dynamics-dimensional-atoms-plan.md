# Dimensionality-aware atom construction — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make realised atoms match the operator's fiber dimension (scalar/tangent/quaternion) and carry the direction their fiber requires, fixing the `bst_eigenfilter('Atom')` mis-seed and adding a default direction, `SeedDir` persistence, a GUI direction picker, and a quiver+norm display.

**Architecture:** The correct dimensionality-aware, direction-carrying seeding **already exists** in `bst_eigenwavelet('Atom')` (builds the delta in the source layout per variant, maps to the eigenbasis via `bst_eigenfilter`'s `RowMap`). This plan **consolidates onto it**: fix the parallel, buggy `bst_eigenfilter('Atom')` (the one the panel calls) to seed the same way with a `seedDir` argument, then reuse `ToVec`/`Op.Frame` for display decode. Refines spec §3.1 (no new `i_operator_fiber` resolver — `RowMap` is the embed, `ToVec` is the decode; only a tiny `field_type`→C helper is added).

**Tech Stack:** MATLAB, Brainstorm; `bst_eigenfilter`/`bst_eigenwavelet`/`manifold_ft`/`RowMap`, the nxr operator registry (`Op.Registry.Primary.field_type`), `figure_3d` `QuiverVectorOverride`.

## Global Constraints

- **Dimensionality source of truth:** `ax.Operator.Registry.Primary.field_type` ∈ {`real`,`complex`,`quaternion`} → C ∈ {1,2,4}; fallback `C = size(ax.Phi{1},1)/numel(ax.GlobalVertices{1})`.
- **Interleaved quaternion layout** `[w,x,y,z]` per vertex; physical current = imaginary slots `(v-1)*4+[2,3,4]`, `w=0`. Carry the FULL quaternion through project→filter→reconstruct; extract the imaginary 3-vector only at the end.
- **Default direction:** Dirac → seed vertex surface normal (`VertNormals`); Connection Laplacian → complex `1` (frame e1, angle 0); scalar → `1`.
- **Backward compatibility:** `bst_eigenfilter('Atom')`'s new `seedDir` is trailing/optional; the scalar path and the B/C cached-projection Apply path must be unaffected.
- **Reuse, don't duplicate:** seeding mirrors `bst_eigenwavelet('Atom')` (lines 314-333) + `RowMap`; Dirac decode uses `bst_eigenwavelet('ToVec')`; never re-derive the layout.
- **No `clear` while Brainstorm runs** (use `rehash`); prefer `matlab -batch` for headless tests (Apple-silicon GUI session drops `GlobalData`).
- **Commit trailers** (development branch): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>` + `Claude-Session: https://claude.ai/code/session_01EXBb3MmD4g8QcFi8NCCByL`.

---

### Task 1: `field_type`→dimensionality helper in `bst_eigenfilter`

**Files:**
- Modify: `toolbox/eigen/bst_eigenfilter.m` (add local `Fiber` verb + `i_fiber` local)
- Test: `dev/test_atom_fiber.m` (new)

**Interfaces:**
- Produces: `[C, kind] = bst_eigenfilter('Fiber', ax)` — `C` ∈ {1,2,4}, `kind` ∈ {`scalar`,`tangent`,`quaternion`}. Reads `ax.Operator.Registry.Primary.field_type` (`real`/`complex`/`quaternion`); falls back to `size(ax.Phi{1},1)/numel(ax.GlobalVertices{1})`.

- [ ] **Step 1: Write the failing test** (`dev/test_atom_fiber.m`)

```matlab
function tests = test_atom_fiber
tests = functiontests(localfunctions);
end
function test_fieldtype_map(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    map = {'Laplace-Beltrami',1,'scalar'; 'Connection Laplacian',2,'tangent'; 'Dirac',4,'quaternion'};
    for i=1:size(map,1)
        ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant',map{i,1},'nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
        [C,kind] = bst_eigenfilter('Fiber', ax);
        verifyEqual(tc, C, map{i,2}, sprintf('%s C', map{i,1}));
        verifyEqual(tc, kind, map{i,3}, sprintf('%s kind', map{i,1}));
    end
end
function test_layout_fallback(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    ax.Operator = [];                          % strip Registry -> must fall back to Phi layout
    [C,kind] = bst_eigenfilter('Fiber', ax);
    verifyEqual(tc, C, 4); verifyEqual(tc, kind, 'quaternion');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (headless): `matlab -batch "addpath(genpath('/Users/diellorbasha/workspace/research/code/brainstorm3/toolbox')); brainstorm server; results = runtests('dev/test_atom_fiber.m'); disp(results); exit"`
Expected: FAIL — `bst_eigenfilter('Fiber', …)` unknown verb.

- [ ] **Step 3: Implement `Fiber` verb + `i_fiber` local** in `bst_eigenfilter.m`

Add the verb to the `eval(macro_method)` set (it already dispatches locals), and add:

```matlab
% Fiber dimensionality of the operator: C components/vertex + kind tag. Source of truth is the nxr
% registry field_type (real/complex/quaternion); falls back to the Phi row/vertex ratio.
function [C, kind] = Fiber(ax) %#ok<DEFNU>
    [C, kind] = i_fiber(ax);
end
function [C, kind] = i_fiber(ax)
    C = [];  ft = '';
    if isfield(ax,'Operator') && isstruct(ax.Operator) && isfield(ax.Operator,'Registry') ...
            && ~isempty(ax.Operator.Registry) && isfield(ax.Operator.Registry,'Primary') ...
            && ~isempty(ax.Operator.Registry.Primary) && isfield(ax.Operator.Registry.Primary,'field_type')
        ft = ax.Operator.Registry.Primary.field_type;
    end
    switch lower(ft)
        case 'real',       C = 1;
        case 'complex',    C = 2;
        case 'quaternion', C = 4;
    end
    if isempty(C)                                   % pre-registry binary -> derive from the Phi layout
        nV = numel(ax.GlobalVertices{1});
        C  = round(size(ax.Phi{1},1) / max(nV,1));
    end
    switch C
        case 1, kind = 'scalar';
        case 2, kind = 'tangent';
        case 4, kind = 'quaternion';
        otherwise, kind = 'scalar';
    end
end
```

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS.
- [ ] **Step 5: Commit** — `git add toolbox/eigen/bst_eigenfilter.m dev/test_atom_fiber.m && git commit` (message: `feat(eigen): bst_eigenfilter Fiber verb — operator dimensionality from field_type`).

---

### Task 2: Fix `bst_eigenfilter('Atom')` — source-layout seed + `seedDir` + `RowMap`

**Files:**
- Modify: `toolbox/eigen/bst_eigenfilter.m` (the `Atom` function, ~line 66)
- Test: `dev/test_atom_seed.m` (new)

**Interfaces:**
- Consumes: `bst_eigenfilter('Fiber', ax)` (Task 1); `RowMap(F, ax, h)` (existing local); `manifold_ft` (existing).
- Produces: `[W, gv] = bst_eigenfilter('Atom', ax, kernelName, kernelParams, seedVert, seedDir)` — `seedDir` optional. Scalar: amplitude (default 1). Tangent: complex (default 1). Dirac: 3-vector (default `[1;0;0]`). `W` is the realised field in the eigenbasis reconstruction (`[C·n_block × nT]`).

- [ ] **Step 1: Write the failing test** (`dev/test_atom_seed.m`)

```matlab
function tests = test_atom_seed
tests = functiontests(localfunctions);
end
function test_dirac_seed_locates_at_vertex(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    gv = ax.GlobalVertices{1};  seed = gv(round(numel(gv)/2));       % a vertex that is NOT local index 1
    kp = struct('lmax', max(ax.Lambda{1}(:)));
    [W, gvb] = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [0;0;1]);  % +Z dipole
    n = numel(gvb);
    imag3 = [W(2:4:end,1) W(3:4:end,1) W(4:4:end,1)];  nrm = sqrt(sum(imag3.^2,2));
    [~,ipk] = max(nrm);
    Vxyz = getfield(in_tess_bst(surf,0),'Vertices'); %#ok
    dmm = norm(Vxyz(gvb(ipk),:) - Vxyz(seed,:))*1000;
    verifyLessThan(tc, dmm, 15, 'Dirac impulse must peak within 15mm of the seed (mis-seed regression)');
end
function test_dirac_seed_direction(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    gv = ax.GlobalVertices{1};  seed = gv(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    Wx = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [1;0;0]);
    Wz = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed, [0;0;1]);
    verifyGreaterThan(tc, norm(Wx - Wz), 1e-9, 'different seed directions must give different fields');
end
function test_scalar_unchanged(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    [W, gvb] = bst_eigenfilter('Atom', ax, 'diffusion', kp, seed);      % no seedDir -> default amplitude 1
    verifyEqual(tc, size(W,1), numel(gvb));                             % scalar: one row per vertex
    verifyGreaterThan(tc, W(1), max(W)*0.5);                            % peak at the seed (local idx 1)
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (current `Atom` has no `seedDir` arg; `test_dirac_seed_locates_at_vertex` fails because `sparse(loc,…)` mis-seeds for `seed ≠ vertex 1`).

- [ ] **Step 3: Replace the seed construction in `Atom`.** Current buggy lines:

```matlab
    loc = find(gv == seedVert, 1);
    c0  = manifold_ft(Phi, M, full(sparse(loc,1,1,size(Phi,1),1)));    % seed in the eigenbasis [K x 1]
```

Change the signature to `function [W, gv] = Atom(ax, KernelName, KernelParams, seedVert, seedDir)` and replace those two lines with a source-layout seed embedded via `RowMap` (mirrors `bst_eigenwavelet('Atom')`):

```matlab
    [C, kind] = i_fiber(ax);
    nSrc = 0; for hh=1:numel(ax.GlobalVertices), nSrc = max(nSrc, max(ax.GlobalVertices{hh}(:))); end
    switch kind
        case 'scalar'
            if nargin < 5 || isempty(seedDir), seedDir = 1; end
            F = zeros(nSrc, 1);       F(seedVert) = seedDir(1);
        case 'tangent'
            if nargin < 5 || isempty(seedDir), seedDir = 1; end
            F = complex(zeros(nSrc,1)); F(seedVert) = seedDir(1);           % complex tangent in the frame
        case 'quaternion'
            if nargin < 5 || isempty(seedDir), seedDir = [1;0;0]; end
            seedDir = seedDir(:);
            F = zeros(3*nSrc, 1);     F((seedVert-1)*3 + (1:3)) = seedDir(1:3);   % 3-vector dipole
    end
    [srcRows, dstRows, nrows, mapMsg] = RowMap(F, ax, blk);
    if ~isempty(mapMsg), error('bst_eigenfilter:Atom', '%s', mapMsg); end
    U = zeros(nrows, 1);  if ~isreal(F), U = complex(U); end
    U(dstRows) = F(srcRows);
    c0 = manifold_ft(Phi, M, U);                                          % seed in the eigenbasis [K x 1]
```

(`blk`, `Phi`, `Lam`, `M`, `gv` are already resolved above in `Atom`; the existing `g(λ)` domain switch that follows is unchanged.)

- [ ] **Step 4: Run test to verify it passes** — Expected: 3/3 PASS (`Fiber`/`i_fiber` from Task 1 is in the same file).
- [ ] **Step 5: Commit** — `feat(eigen): dimensionality-aware Atom seed (source-layout + seedDir via RowMap); fixes multi-component mis-seed`.

---

### Task 3: Default-direction helper in the panel

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (add local `i_atom_default_dir`)
- Test: `dev/test_atom_default_dir.m` (new)

**Interfaces:**
- Consumes: `bst_eigenfilter('Fiber', ax)` (Task 1); `in_tess_bst` (existing).
- Produces: `dir = panel_bst_dynamics('i_atom_default_dir', ax, seedVert)` — Dirac: unit `[nx ny nz]` surface normal at `seedVert`; tangent/scalar: `1`.

- [ ] **Step 1: Write the failing test** (`dev/test_atom_default_dir.m`)

```matlab
function tests = test_atom_default_dir
tests = functiontests(localfunctions);
end
function test_dirac_normal(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(10);
    d = panel_bst_dynamics('i_atom_default_dir', ax, seed);
    verifyEqual(tc, numel(d), 3);  verifyEqual(tc, norm(d), 1, 'AbsTol', 1e-6);
    S = in_tess_bst(surf, 0);
    verifyGreaterThan(tc, abs(dot(d(:), S.VertNormals(seed,:)')), 0.99, 'must be the seed surface normal');
end
function test_scalar_amplitude(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',6,'TimeWindow',[0 .04],'SampleRate',100));
    d = panel_bst_dynamics('i_atom_default_dir', ax, ax.GlobalVertices{1}(1));
    verifyEqual(tc, d, 1);
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (`i_atom_default_dir` undefined).
- [ ] **Step 3: Implement `i_atom_default_dir`** in `panel_bst_dynamics.m`

```matlab
% Default impulse direction per operator dimensionality: Dirac -> seed surface normal (unit 3-vector);
% tangent/scalar -> 1 (frame e1 / amplitude). See atom-operator-applicability.
function dir = i_atom_default_dir(ax, seedVert) %#ok<DEFNU>
    [~, kind] = bst_eigenfilter('Fiber', ax);
    if ~strcmp(kind, 'quaternion'), dir = 1; return; end
    dir = [0 0 1];                                       % fallback if normals are missing
    try
        S = in_tess_bst(ax.SurfaceFile, 0);
        if isfield(S,'VertNormals') && ~isempty(S.VertNormals) && seedVert <= size(S.VertNormals,1)
            n = S.VertNormals(seedVert, :);  if norm(n) > 0, dir = n / norm(n); end
        end
    catch %#ok<CTCH>
    end
    dir = dir(:)';
end
```

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS.
- [ ] **Step 5: Commit** — `feat(dynamics): default impulse direction (Dirac normal / scalar amplitude)`.

---

### Task 4: `G.SeedDir` in the atom data model

**Files:**
- Modify: `toolbox/db/db_template.m:360` (add `SeedDir` to the `atomgroup` template)
- Test: `dev/test_atom_seeddir_model.m` (new)

**Interfaces:**
- Produces: `db_template('atomgroup')` has field `SeedDir` (default `[]`); persists through `bst_dynamics('Save'/'Load')` and `NewGroup`/`AddGroup`.

- [ ] **Step 1: Write the failing test** (`dev/test_atom_seeddir_model.m`)

```matlab
function tests = test_atom_seeddir_model
tests = functiontests(localfunctions);
end
function test_template_has_seeddir(tc)
    G = db_template('atomgroup');
    verifyTrue(tc, isfield(G, 'SeedDir'));
    verifyEmpty(tc, G.SeedDir);
end
function test_roundtrip(tc)
    T = bst_dynamics('New', 'test');
    G = bst_dynamics('NewGroup', 'a'); G.vertices = 5; G.Operator = 'Dirac'; G.SeedDir = [0 0 1];
    T = bst_dynamics('AddGroup', T, G);
    f = [tempname '.mat'];  bst_dynamics('Save', f, T);  T2 = bst_dynamics('Load', f);  delete(f);
    verifyEqual(tc, T2.Groups(1).SeedDir, [0 0 1]);
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (`SeedDir` not in template).
- [ ] **Step 3: Add the field** at `db_template.m` (after `KernelParams`, line 360):

```matlab
            'KernelParams', [], ...     % GENERATOR: kernel param struct (bst_eigfilter_controls('ToKernel'))
            'SeedDir',      [], ...     % GENERATOR: impulse direction -- [] scalar | complex tangent | [nx ny nz] Dirac
```

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS (`bst_dynamics('Load')` template-fills via `struct_copy_fields`, so old files get `SeedDir=[]`).
- [ ] **Step 5: Commit** — `feat(dynamics): persist atom SeedDir (impulse direction) in the group model`.

---

### Task 5: Panel realiser — pass `SeedDir`, return decoded `V3`

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` — `i_atom_realise` (line ~595), `i_atom_preview_impulse` (~635), `OnCreateAtom` (~1035)
- Test: `dev/test_atom_realise_v3.m` (new)

**Interfaces:**
- Consumes: `bst_eigenfilter('Atom', ax, k, kp, seed, dir)` (Task 2); `bst_eigenwavelet('ToVec', W, EigenMat)` (existing, Dirac decode); `i_atom_default_dir` (Task 3); `Op.Frame` for tangent decode.
- Produces: `[W, gv, isSigned, V3] = i_atom_realise(st, kernel, kp, seed, variant, seedDir)` — `V3` is `[nV × 3]` full-surface ambient vectors (Dirac: imag slots; tangent: `a·e1+b·e2`; scalar: `[]`). `OnCreateAtom` stores `G.SeedDir = i_atom_default_dir(ax, seed)`.

- [ ] **Step 1: Write the failing test** (`dev/test_atom_realise_v3.m`)

```matlab
function tests = test_atom_realise_v3
tests = functiontests(localfunctions);
end
function test_dirac_v3_shape_and_dir(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    st.SurfaceFileOverride = surf;                                  % i_atom_axes reads ax.SurfaceFile from st
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Dirac','nModes',80,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    % call the pure realise-core directly (exposed for tests): [W,gv,V3] = i_atom_realise_core(ax,k,kp,seed,dir)
    [~, gv, V3] = panel_bst_dynamics('i_atom_realise_core', ax, 'diffusion', kp, seed, [0 0 1]); %#ok
    nV = 0; for h=1:numel(ax.GlobalVertices), nV=max(nV,max(ax.GlobalVertices{h}(:))); end
    verifyEqual(tc, size(V3), [nV 3]);
    [~,ipk] = max(sqrt(sum(V3.^2,2)));
    verifyGreaterThan(tc, abs(V3(ipk,3)), max(abs(V3(ipk,1)), abs(V3(ipk,2))), 'peak dipole ~ +Z');
end
function test_scalar_v3_empty(tc)
    surf = getenv('BST_TEST_SURF'); if isempty(surf), surf='sub-MTL0002/tess_cortex_pial_low.mat'; end
    ax = bst_eigen('Axes', struct('SurfaceFile',surf,'Variant','Laplace-Beltrami','nModes',40,'TimeWindow',[0 .04],'SampleRate',100));
    seed = ax.GlobalVertices{1}(1);  kp = struct('lmax', max(ax.Lambda{1}(:)));
    [~,~,V3] = panel_bst_dynamics('i_atom_realise_core', ax, 'diffusion', kp, seed, 1);
    verifyEmpty(tc, V3);
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (`i_atom_realise_core` undefined).

- [ ] **Step 3: Add `i_atom_realise_core` + rewire `i_atom_realise`.** Add a pure core that both the test and `i_atom_realise` call:

```matlab
% Pure realise-core: run the atom on ax at (seed,dir); return the raw field W [C*n x nT], its global
% vertices gv, and the decoded full-surface ambient vectors V3 [nV x 3] ([] for scalar).
function [W, gv, V3] = i_atom_realise_core(ax, kernel, kp, seed, seedDir) %#ok<DEFNU>
    [W, gv] = bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir);
    [C, kind] = bst_eigenfilter('Fiber', ax);
    nV = 0; for h=1:numel(ax.GlobalVertices), nV = max(nV, max(ax.GlobalVertices{h}(:))); end
    V3 = [];
    switch kind
        case 'quaternion'
            n = numel(gv);  im = [W(2:4:end,1) W(3:4:end,1) W(4:4:end,1)];   % imag 3-vector at frame 1
            V3 = zeros(nV,3);  V3(gv,:) = im;
        case 'tangent'
            % complex (a+bi) in the operator frame -> ambient 3-vector a*e1 + b*e2 (Op.Frame per hemi)
            blk = 1; for h=1:numel(ax.GlobalVertices), if any(ax.GlobalVertices{h}==seed), blk=h; break; end, end
            Fr = ax.Operator.Frame{blk};  a = real(W(:,1));  b = imag(W(:,1));
            V3 = zeros(nV,3);  V3(gv,:) = a.*Fr.e1 + b.*Fr.e2;
    end
    if C == 1, V3 = []; end
end
```

Then in `i_atom_realise` (which currently returns `[W, gv, isSigned]`), add `seedDir` input + `V3` output and delegate the field build to the core (keep the existing `i_paintable_scalar` + normalisation for `W`):

```matlab
function [W, gv, isSigned, V3] = i_atom_realise(st, kernel, kp, seed, variant, seedDir)
    if (nargin < 5) || isempty(variant), variant = 'Laplace-Beltrami'; end
    W = []; gv = []; isSigned = false; V3 = [];
    ax = i_atom_axes(st, variant);  if isempty(ax), return; end
    if ~isstruct(kp), kp = struct(); end
    if ~isfield(kp,'lmax') || isempty(kp.lmax), kp.lmax = max(ax.Lambda{1}(:)); end
    if (nargin < 6) || isempty(seedDir), seedDir = i_atom_default_dir(ax, seed); end
    try
        [Wraw, gv, V3] = i_atom_realise_core(ax, kernel, kp, seed, seedDir);
    catch %#ok<CTCH>
        W = []; gv = []; V3 = []; return;
    end
    nGv = numel(gv);
    W = i_paintable_scalar(Wraw, nGv);
    if size(W,1) ~= nGv, W = []; V3 = []; return; end
    if any(strcmp(variant, {'Laplace-Beltrami','LB-Connectome'}))
        [W, isSigned] = i_atom_normalize(W, ax.Mass{i_seed_block(ax, seed)});
    else
        pk = max(abs(W(:)));  if pk > 0, W = W / pk; end;  isSigned = false;
    end
end
```

Update `i_atom_preview_impulse` to capture `V3` and pass it on (Task 6 consumes it):

```matlab
    [W, gv, isSigned, V3] = i_atom_realise(st, k, kp, seed, variant, i_field(st,'atomSeedDir',[]));
    ...
        view_dynamics('SetAtomField', st.hFig, W, gv, isSigned, V3);
```

In `OnCreateAtom`, store the default direction on the new atom (after `ax`/`seed` are known, before building `G`):

```matlab
    sdir = i_atom_default_dir(ax, seed);
    ... G = i_default_atom('diffusion', kp, seed, ax.SurfaceFile, sprintf('atom%d', numel(st.T.Groups)+1), op);
    G.SeedDir = sdir;
```

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS.
- [ ] **Step 5: Commit** — `feat(dynamics): realise atoms with direction + decode ambient V3 for display`.

---

### Task 6: Quiver+norm display (`view_dynamics` / `figure_3d`)

**Files:**
- Modify: `toolbox/gui/view_dynamics.m` — `SetAtomField` (~278), `ClearAtomField` (~302)
- Test: `dev/test_atom_quiver_overlay.m` (new; headless appdata check)

**Interfaces:**
- Consumes: `V3` from `i_atom_realise` (Task 5); `figure_3d('SetShowSourceVectors', hFig, iTess, onoff)` + `QuiverVectorOverride` appdata (existing).
- Produces: `view_dynamics('SetAtomField', hFig, W, gv, isSigned, V3)` — when `V3` non-empty sets `QuiverVectorOverride` + source vectors on; empty/absent clears them. `ClearAtomField` also clears the override + disables source vectors.

- [ ] **Step 1: Write the failing test** (`dev/test_atom_quiver_overlay.m`) — pure appdata contract on a bare figure:

```matlab
function tests = test_atom_quiver_overlay
tests = functiontests(localfunctions);
end
function test_v3_sets_override(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig));
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',false));
    V3 = [0 0 1; 0 0 1; 0 0 0; 0 0 0];
    view_dynamics('SetAtomField', hFig, [1;1;0;0], (1:4)', false, V3);
    verifyEqual(tc, getappdata(hFig,'QuiverVectorOverride'), V3);
end
function test_empty_v3_clears(tc)
    hFig = figure('Visible','off');  c = onCleanup(@() close(hFig));
    setappdata(hFig, 'DynamicsOverlay', struct('Op','none','iTess',1,'nV',4,'srcDS',[],'srcResult',[], ...
        'AtomField',[],'AtomGV',[],'AtomSigned',false,'AtomWin',[]));
    setappdata(hFig, 'Surface', struct('hPatch',[],'DataSource',struct('Type','Source','FileName','x'), ...
        'ShowSourceVectors',true));
    setappdata(hFig, 'QuiverVectorOverride', [0 0 1]);
    view_dynamics('SetAtomField', hFig, [1;1;0;0], (1:4)', false, []);
    verifyEmpty(tc, getappdata(hFig,'QuiverVectorOverride'));
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (`SetAtomField` has no `V3` arg / ignores it).

- [ ] **Step 3: Extend `SetAtomField` + `ClearAtomField`.** Change the signature to `SetAtomField(hFig, W, gv, isSigned, V3)` (`if nargin<5, V3=[]; end`). After the existing scalar-paint body (which sets `D.AtomField`, colormap, `i_dynamics_overlay(hFig)`), append:

```matlab
    if (nargin >= 5) && ~isempty(V3) && (size(V3,2) == 3)
        setappdata(hFig, 'QuiverVectorOverride', V3);
        try, figure_3d('SetShowSourceVectors', hFig, D.iTess, 1); catch, end %#ok<CTCH>
    else
        setappdata(hFig, 'QuiverVectorOverride', []);
        try, figure_3d('SetShowSourceVectors', hFig, D.iTess, 0); catch, end %#ok<CTCH>
    end
```

In `ClearAtomField`, after the existing clear, add:

```matlab
    setappdata(hFig, 'QuiverVectorOverride', []);
    try, D = getappdata(hFig,'DynamicsOverlay'); if ~isempty(D) && isfield(D,'iTess'), figure_3d('SetShowSourceVectors', hFig, D.iTess, 0); end, catch, end %#ok<CTCH>
```

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS. (Note: `SetShowSourceVectors` on the bare test figure is guarded by `try`; `PlotSourceVectors` no-ops without a real patch — the test only asserts the appdata contract.)
- [ ] **Step 5: Commit** — `feat(dynamics): draw Dirac/tangent atoms as quivers over the norm colormap`.

---

### Task 7: GUI direction picker (Atom section)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` — Atom-section build (near the operator radios / atom params), `i_select_atom_load` (~1087), a new `OnPickSeedDir` + `i_build_dir_control`.
- Test: `dev/test_dir_control.m` (new; the control-model builder is pure) + **live** verification (MCP).

**Interfaces:**
- Consumes: `bst_eigenfilter('Fiber', ax)` (Task 1); `i_atom_default_dir` (Task 3); the existing cortex vertex pick (`figure_3d`).
- Produces: `spec = panel_bst_dynamics('i_dir_control_spec', kind)` → struct `{show, type}` (`scalar`→`show=false`; `tangent`→`type='angle'`; `quaternion`→`type='preset'` with presets `{Normal,+X,+Y,+Z,Pick-on-surface}`). Setting a direction updates `st.atomSeedDir` + the selected atom's `G.SeedDir` and re-previews via `i_atom_preview`.

- [ ] **Step 1: Write the failing test** (`dev/test_dir_control.m`)

```matlab
function tests = test_dir_control
tests = functiontests(localfunctions);
end
function test_control_spec(tc)
    s = panel_bst_dynamics('i_dir_control_spec', 'scalar');     verifyFalse(tc, s.show);
    a = panel_bst_dynamics('i_dir_control_spec', 'tangent');    verifyTrue(tc, a.show); verifyEqual(tc, a.type, 'angle');
    q = panel_bst_dynamics('i_dir_control_spec', 'quaternion'); verifyTrue(tc, q.show); verifyEqual(tc, q.type, 'preset');
    verifyTrue(tc, ismember('Normal', q.presets));
    verifyTrue(tc, ismember('Pick-on-surface', q.presets));
end
function test_preset_to_dir(tc)
    verifyEqual(tc, panel_bst_dynamics('i_preset_dir', '+Z'), [0 0 1]);
    verifyEqual(tc, panel_bst_dynamics('i_preset_dir', '+X'), [1 0 0]);
end
```

- [ ] **Step 2: Run test to verify it fails** — Expected: FAIL (`i_dir_control_spec`/`i_preset_dir` undefined).

- [ ] **Step 3: Implement the control model + wiring.** Add the pure helpers:

```matlab
function s = i_dir_control_spec(kind) %#ok<DEFNU>
    switch kind
        case 'tangent',    s = struct('show',true,'type','angle','presets',{{}});
        case 'quaternion', s = struct('show',true,'type','preset','presets',{{'Normal','+X','+Y','+Z','Pick-on-surface'}});
        otherwise,         s = struct('show',false,'type','none','presets',{{}});
    end
end
function d = i_preset_dir(name) %#ok<DEFNU>
    switch name
        case 'Normal', d = [];                 % resolved to the seed normal by the caller
        case '+X', d = [1 0 0];  case '+Y', d = [0 1 0];  case '+Z', d = [0 0 1];
        otherwise, d = [];
    end
end
```

Then, in the Atom-section GUI build, add a direction control (a `JComboBox` for `quaternion`, a numeric angle field for `tangent`, hidden for `scalar`) built from `i_dir_control_spec(kind)` where `kind = second output of bst_eigenfilter('Fiber', i_atom_axes(st, i_atom_op(st)))`. On change:
- preset `Normal` → `i_atom_default_dir(ax, seed)`; `+X/+Y/+Z` → `i_preset_dir(name)`; `Pick-on-surface` → arm a one-shot pick (`figure_3d` returns the clicked vertex; `dir = unit(Vxyz(pick,:) - Vxyz(seed,:))`); `angle θ` → `complex(cos θ, sin θ)`.
- store on `st.atomSeedDir` and the selected atom `st.T.Groups(curAtom).SeedDir`, then `i_atom_preview()`.
- rebuild/refresh the control in `i_select_atom_load` and `OnSetOperator` so it matches the current operator's `kind` (hidden for scalar).

- [ ] **Step 4: Run test to verify it passes** — Expected: 2/2 PASS.
- [ ] **Step 5: Live verify (MCP/headless GUI):** launch a Dirac dSPM Dynamics session → Create atom → quivers along the seed normal over the norm colormap; set preset `+Z` → quivers reorient to +Z; `Pick-on-surface` → click a vertex → quivers point seed→target; a Laplace-Beltrami atom shows no direction control and no quivers. Screenshot each.
- [ ] **Step 6: Commit** — `feat(dynamics): GUI impulse-direction picker (preset/angle/pick-on-surface)`.

---

## Self-Review

**Spec coverage:** resolver/dimensionality → Task 1 (`Fiber`, refines spec §3.1 by reusing `RowMap`); seeding fix + direction → Task 2; default direction → Task 3; data model `SeedDir` → Task 4; realise + decode `V3` → Task 5; quiver+norm display → Task 6; GUI picker → Task 7. Out-of-scope items (face domains, real/conformal companion display, Apply-path direction) correctly absent.

**Placeholder scan:** none — every code step has concrete code; the headless run command is spelled out; the tangent `field_type` string (`complex`) is confirmed, not deferred.

**Type consistency:** `bst_eigenfilter('Fiber', ax)→[C,kind]` used identically in Tasks 2/3/5/7; `bst_eigenfilter('Atom', ax, kernel, kp, seed, seedDir)` signature consistent Tasks 2/5; `i_atom_realise_core(ax,kernel,kp,seed,seedDir)→[W,gv,V3]` consistent Tasks 5/(test); `SetAtomField(hFig,W,gv,isSigned,V3)` consistent Tasks 5/6; `G.SeedDir` consistent Tasks 4/5/7; decode uses interleaved imag slots `(2:4:end,3:4:end,4:4:end)` matching `RowMap`/`ToVec`.
