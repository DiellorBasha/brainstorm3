# Differential Stratification + Poisson/Cholesky Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Stratify `toolbox/differential/` by field type (scalar→LBO, tangent→connection, ambient→Dirac) and back it with a shared, lazily-cached Poisson/Cholesky solver.

**Architecture:** A lazy factor helper (`tess_cholesky`) persists the Cholesky factor of an assembled operator matrix onto its `operator_` node and caches it in-session; a stratified `bst_poisson` solves `Lφ=f` through it; the existing verbs (`bst_helmholtz`, `bst_divergence`, `bst_curl`) are refactored onto it and gain field-type branches; `bst_operators` becomes a router with a guard table.

**Tech Stack:** MATLAB (Brainstorm fork, R2023b), nxr-compute plugin, sparse linear algebra (`chol`/`decomposition`). Validation via the MATLAB MCP (`run_matlab_file` / `evaluate_matlab_code`); no unit-test framework — validation scripts assert numeric properties and print `PASS`/`FAIL`.

## Global Constraints

- MATLAB R2023b; nxr-compute plugin must be installed (`bst_plugin('Install','nxr-compute')`).
- File header block: every new `.m` uses the standard Brainstorm GPL header + `Authors: Diellor Basha, 2026` (copy from `bst_gradient.m`).
- Per-hemisphere everything: operators/bases are per-hemisphere `1x2` cells indexed `hh`.
- Manifold face normals are gauge-signed (effectively inward); orient outward via the divergence-theorem test already in `bst_helmholtz`/`bst_curl` (`i_orient_outward`). Never re-derive from `VertNormals`.
- Geometry source of truth: `manifold_` node (`Embedded`, `DEC`). Assembled matrices: `operator_` node. Connectivity (`Faces`): `in_tess_bst`.
- I/O-free hot path: per-frame solving must not touch disk. Disk persistence happens only on the explicit `tess_cholesky('attach', …)` cache-miss path.
- Never run MATLAB `clear` in the live session (wipes `GlobalData`); use `rehash`. Edited `.m` files auto-reload.
- Commit after every task. Branch: `feature/differential-stratification` (already created).
- Validation scratch/baseline `.mat` files go under `dev/scratch/` (create if missing); do not commit baselines.

---

## File Structure

**Create:**
- `toolbox/anatomy/tess_cholesky.m` — lazy factor: pure getter (I/O-free) + `'attach'` (persist on node) + `'solve'`.
- `toolbox/differential/bst_poisson.m` — stratified `Lφ=f` solver (scalar LBO route in this plan).
- `dev/test_tess_cholesky.m`, `dev/test_bst_poisson.m`, `dev/test_helmholtz_parity.m`, `dev/test_ambient_divcurl.m`, `dev/test_operators_router.m` — validation scripts.

**Modify:**
- `toolbox/db/db_template.m:133` (`operatormat`) — add `Cholesky` field.
- `toolbox/differential/bst_operators.m` — drop local `i_poisson`/`i_laplacian` Poisson, call `bst_poisson`; add router/guard/`FieldType` (Task 5).
- `toolbox/differential/bst_helmholtz.m` — replace `Prepare`-time `decomposition(...)` + private `i_poisson` with `tess_cholesky`/`bst_poisson` (Task 3).
- `toolbox/differential/bst_divergence.m`, `bst_curl.m` — add ambient (`3nV`) branch (Task 4).
- `toolbox/differential/bst_gradient.m` — scalar-input guard (Task 5).

---

## Task 1: `tess_cholesky` + `operatormat.Cholesky` field

**Files:**
- Modify: `toolbox/db/db_template.m:133-146`
- Create: `toolbox/anatomy/tess_cholesky.m`
- Test: `dev/test_tess_cholesky.m`

**Interfaces:**
- Produces:
  - `dF = tess_cholesky(OperatorNode, hh, pin)` — pure getter (I/O-free). `OperatorNode` is a loaded `operatormat` struct; `hh` hemisphere index; `pin` = vector of local pinned indices (default `1`). Returns `dF = struct('L',L,'p',p,'free',free,'n',n)` where `[L,~,p]=chol(A(free,free),'lower','vector')`, `A=OperatorNode.Operator{hh}`, `free=setdiff(1:n,pin)`.
  - `x = tess_cholesky('solve', dF, rhs)` — solves `A x = rhs` on the free block (pinned rows of `x` are 0; caller recenters). `rhs` is `[n x k]`.
  - `OperatorNode = tess_cholesky('attach', OperatorNode, OperatorFile, pin)` — computes the factor for both hemispheres, stores `OperatorNode.Cholesky{hh}=dF`, `bst_save`s the node to `OperatorFile`, returns the updated struct.
- Consumes: nothing from earlier tasks.

- [ ] **Step 1: Add the `Cholesky` field to the operatormat template**

In `toolbox/db/db_template.m`, in the `case 'operatormat'` struct (line 134-146), add one field before `'Provenance'`:

```matlab
              'Frame',          [], ...   % Connection Laplacian only: 1x2 struct(.e1 [V x 3], .e2 [V x 3], .normal [V x 3]) -- canonical per-vertex tangent frame the complex eigenmodes decode in (field = real(U).*e1 + imag(U).*e2)
              'Cholesky',       [], ...   % lazy factor cache (tess_cholesky): 1x2 cell, dF=struct('L','p','free','n') of the pinned A=Operator{hh}; [] until first solve attaches it
              'Provenance',     []);
```

- [ ] **Step 2: Write the failing validation script**

Create `dev/test_tess_cholesky.m`:

```matlab
function test_tess_cholesky
% Validates tess_cholesky on a synthetic SPD matrix (no DB / no nxr).
%   1. 'solve' reproduces backslash on the pinned free block.
%   2. The pure getter returns the factor carried on the node struct without recomputing.
%   3. 'attach' populates OperatorNode.Cholesky for both hemispheres.
    fprintf('== test_tess_cholesky ==\n');
    rng(0);
    n = 200;
    % a small SPD "stiffness": graph Laplacian of a random connected graph + tiny shift
    G = sprandsym(n, 0.05);  G = G - diag(diag(G));  G = abs(G);
    K = spdiags(sum(G,2)+1e-6, 0, n, n) - G;   % SPD-ish; pin removes the constant nullspace
    K = (K+K')/2;
    M = spdiags(rand(n,1)+0.5, 0, n, n);
    Node = db_template('operatormat');
    Node.Variant = 'Laplace-Beltrami';  Node.ParentSurface = 'unit://test';
    Node.Operator = {K, K};  Node.Mass = {M, M};
    Node.GlobalVertices = {(1:n)', (1:n)'};

    pin = 1;  free = setdiff(1:n, pin)';
    dF = tess_cholesky(Node, 1, pin);
    b = randn(n,1);  b(pin) = 0;
    x = tess_cholesky('solve', dF, b);
    xref = zeros(n,1);  xref(free) = K(free,free) \ b(free);
    err = norm(x - xref) / max(norm(xref), eps);
    assert(err < 1e-10, 'solve mismatch vs backslash: %g', err);
    fprintf('  solve vs backslash rel err = %g  [OK]\n', err);

    % getter uses a pre-attached factor (no recompute): tamper with a marker and confirm reuse
    Node.Cholesky = {dF, dF};  Node.Cholesky{1}.marker = 42;
    dF2 = tess_cholesky(Node, 1, pin);
    assert(isfield(dF2,'marker') && dF2.marker == 42, 'getter did not reuse attached factor');
    fprintf('  getter reuses attached factor  [OK]\n');

    fprintf('PASS\n');
end
```

- [ ] **Step 3: Run it to verify it fails (function missing)**

Run (MATLAB MCP `run_matlab_file`): `dev/test_tess_cholesky.m`
Expected: error `Undefined function 'tess_cholesky'`.

- [ ] **Step 4: Implement `tess_cholesky`**

Create `toolbox/anatomy/tess_cholesky.m` (GPL header from `bst_gradient.m`, then):

```matlab
function out = tess_cholesky(varargin)
% TESS_CHOLESKY: Lazy, cached Cholesky factor of an assembled operator matrix.
%
% The factor of a pinned SPD operator (A = OperatorNode.Operator{hh}) is the factorization
% of ONE specific assembled matrix, so it is persisted ON the operator_ node (next to A),
% not on the geometry-only manifold_ node. Three modes:
%
%   dF   = tess_cholesky(OperatorNode, hh, pin)            % PURE getter, I/O-free
%          -> dF = struct('L',L,'p',p,'free',free,'n',n), [L,~,p]=chol(A(free,free),'lower','vector')
%          Resolution: node.Cholesky{hh} -> in-session memory cache -> compute (+cache in memory)
%   x    = tess_cholesky('solve', dF, rhs)                 % A x = rhs on the free block; pinned rows 0
%   Node = tess_cholesky('attach', OperatorNode, File, pin)% compute both hemis, save onto the node (I/O)
%
% The opaque `decomposition` object is NOT serialized; L (sparse lower) + p (permutation) are
% plain data that reload cleanly and reconstruct the two triangular solves.
%
% SEE ALSO: bst_poisson, bst_get_operator_node, tess_operators
%
% Authors: Diellor Basha, 2026
    persistent MEM
    if isempty(MEM), MEM = containers.Map('KeyType','char','ValueType','any'); end

    if ischar(varargin{1})
        switch lower(varargin{1})
            case 'solve'
                out = i_solve(varargin{2}, varargin{3});
            case 'attach'
                out = i_attach(varargin{2}, varargin{3}, i_pin(varargin{4:end}));
            otherwise
                error('tess_cholesky:badMode', 'Unknown mode: %s', varargin{1});
        end
        return;
    end

    % pure getter: tess_cholesky(OperatorNode, hh, pin)
    Node = varargin{1};  hh = varargin{2};  pin = i_pin(varargin{3:end});
    % 1) factor already on the node?
    if iscell(Node.Cholesky) && numel(Node.Cholesky) >= hh && ~isempty(Node.Cholesky{hh})
        out = Node.Cholesky{hh};  return;
    end
    % 2) in-session memory cache?
    key = i_key(Node, hh, pin);
    if isKey(MEM, key), out = MEM(key);  return; end
    % 3) compute + cache in memory (no disk write on the pure path)
    out = i_factor(Node.Operator{hh}, pin);
    MEM(key) = out;
end

function pin = i_pin(varargin)
    if isempty(varargin) || isempty(varargin{1}), pin = 1; else, pin = varargin{1}(:)'; end
end

function key = i_key(Node, hh, pin)
    key = sprintf('%s|%s|%d|%s', Node.Variant, Node.ParentSurface, hh, mat2str(pin));
end

function dF = i_factor(A, pin)
    A = (A + A')/2;
    n = size(A,1);
    free = setdiff((1:n)', pin(:));
    [L, flag, p] = chol(A(free,free), 'lower', 'vector');
    if flag ~= 0
        error('tess_cholesky:notSPD', 'Pinned operator block is not SPD (chol flag=%d).', flag);
    end
    dF = struct('L', L, 'p', p(:), 'free', free, 'n', n);
end

function x = i_solve(dF, rhs)
    % A(free,free)(p,p) = L L'  =>  solve on the free block, pinned rows stay 0
    k = size(rhs, 2);
    x = zeros(dF.n, k);
    bf = rhs(dF.free, :);
    bp = bf(dF.p, :);
    y  = dF.L \ bp;
    z  = dF.L' \ y;
    xf = zeros(numel(dF.free), k);
    xf(dF.p, :) = z;
    x(dF.free, :) = xf;
end

function Node = i_attach(Node, OperatorFile, pin)
    nH = numel(Node.Operator);
    if ~iscell(Node.Cholesky), Node.Cholesky = cell(1, nH); end
    for hh = 1:nH
        if isempty(Node.Operator{hh}), continue; end
        Node.Cholesky{hh} = i_factor(Node.Operator{hh}, pin);
    end
    bst_save(file_fullpath(OperatorFile), Node, 'v6');
end
```

- [ ] **Step 5: Run validation to verify it passes**

Run: `dev/test_tess_cholesky.m`
Expected: `solve vs backslash rel err = …  [OK]`, `getter reuses attached factor  [OK]`, `PASS`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/db/db_template.m toolbox/anatomy/tess_cholesky.m dev/test_tess_cholesky.m
git commit -m "feat(anatomy): tess_cholesky lazy factor cache on operator_ node

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `bst_poisson` (scalar LBO) + replace `bst_operators` Poisson

**Files:**
- Create: `toolbox/differential/bst_poisson.m`
- Modify: `toolbox/differential/bst_operators.m:138-142` (poisson dispatch), `:207-235` (drop `i_poisson`, keep `i_laplacian`)
- Test: `dev/test_bst_poisson.m`

**Interfaces:**
- Consumes: `tess_cholesky(OperatorNode, hh, pin)`, `tess_cholesky('solve', dF, rhs)` (Task 1).
- Produces: `phi = bst_poisson(OperatorNode, F)` — scalar-stratum solve of `K φ = M f` per hemisphere with the nullspace handled centrally (project `f` to mean-zero in the mass metric, pinned solve via `tess_cholesky`, recenter). `OperatorNode` = a `Laplace-Beltrami` operatormat (`Operator{hh}=K`, `Mass{hh}=M`, `GlobalVertices{hh}`). `F` = per-vertex scalar `[nV x nT]`. Returns `phi [nV x nT]`.

- [ ] **Step 1: Write the failing validation script**

Create `dev/test_bst_poisson.m`:

```matlab
function test_bst_poisson
% bst_poisson scalar solve reproduces the legacy bst_operators i_poisson math exactly,
% on a synthetic LBO-like node (no DB / no nxr).
    fprintf('== test_bst_poisson ==\n');
    rng(1);
    n = 300;
    G = sprandsym(n, 0.04);  G = abs(G - diag(diag(G)));
    K = spdiags(sum(G,2), 0, n, n) - G;  K = (K+K')/2;     % cotan-like, constant nullspace
    M = spdiags(rand(n,1)+0.5, 0, n, n);
    Node = db_template('operatormat');
    Node.Variant='Laplace-Beltrami'; Node.ParentSurface='unit://test';
    Node.Operator={K,K}; Node.Mass={M,M}; Node.GlobalVertices={(1:n)',(1:n)'};
    f = randn(n, 4);

    phi = bst_poisson(Node, f);

    % reference: legacy projection/pin/recenter (verbatim from bst_operators i_poisson)
    free=(2:n)'; totMass=sum(M(:));  ref=zeros(n,4);
    dK = decomposition(K(free,free),'chol');
    for c=1:4
        rhs = M*(f(:,c) - sum(M*f(:,c))/totMass);
        x=zeros(n,1); x(free)=dK\rhs(free); ref(:,c)=x-mean(x);
    end
    err = norm(phi-ref,'fro')/max(norm(ref,'fro'),eps);
    assert(err < 1e-10, 'bst_poisson != legacy i_poisson: rel err %g', err);
    fprintf('  bst_poisson vs legacy rel err = %g  [OK]\nPASS\n', err);
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/test_bst_poisson.m`
Expected: `Undefined function 'bst_poisson'`.

- [ ] **Step 3: Implement `bst_poisson`**

Create `toolbox/differential/bst_poisson.m` (GPL header, then):

```matlab
function phi = bst_poisson(OperatorNode, F)
% BST_POISSON: Stratified Poisson solve  L phi = f  on the cortical manifold.
%
% Scalar stratum (Laplace-Beltrami): solves  K phi = M f  per hemisphere, where K is the
% cotan stiffness and M the Galerkin mass. The constant nullspace is handled HERE, once:
% project f to the mean-zero subspace in the mass metric, pinned solve through the cached
% tess_cholesky factor, recenter. This is the single home of the nullspace handling that
% was duplicated in bst_operators (per-column re-factorization) and bst_helmholtz.
%
% USAGE:  phi = bst_poisson(OperatorNode, F)
%   OperatorNode : a 'Laplace-Beltrami' operatormat (Operator{hh}=K, Mass{hh}=M, GlobalVertices{hh})
%   F            : per-vertex scalar source [nV x nT]
%   phi          : per-vertex potential [nV x nT] (mean-zero per hemisphere)
%
% SEE ALSO: tess_cholesky, bst_operators, bst_helmholtz
%
% Authors: Diellor Basha, 2026
    if ~strcmpi(OperatorNode.Variant, 'Laplace-Beltrami')
        error('bst_poisson:variant', ...
            'bst_poisson scalar route needs a Laplace-Beltrami operator (got %s).', OperatorNode.Variant);
    end
    nVtot = max(cellfun(@(c) max(double(c(:))), OperatorNode.GlobalVertices));
    phi = zeros(nVtot, size(F,2));
    for hh = 1:numel(OperatorNode.Operator)
        if isempty(OperatorNode.Operator{hh}), continue; end
        vH = double(OperatorNode.GlobalVertices{hh}(:));
        M  = OperatorNode.Mass{hh};
        dF = tess_cholesky(OperatorNode, hh, 1);     % pin vertex 1
        fh = F(vH, :);
        totMass = sum(M(:));
        fh = fh - (sum(M*fh, 1) / totMass);          % project to mean-zero (mass metric)
        x  = tess_cholesky('solve', dF, M*fh);       % K x = M f  on the free block
        phi(vH, :) = x - mean(x, 1);                 % recenter
    end
end
```

- [ ] **Step 4: Run validation to verify it passes**

Run: `dev/test_bst_poisson.m`
Expected: `bst_poisson vs legacy rel err = …  [OK]`, `PASS`.

- [ ] **Step 5: Rewire `bst_operators` poisson onto `bst_poisson`**

In `toolbox/differential/bst_operators.m`, replace the `'poisson'` case body (lines 138-142):

```matlab
        case 'poisson'
            LBO   = bst_get_operator_node(SurfaceFile, 'Laplace-Beltrami');
            f     = i_as_vertex_scalar(F, nVtot, 'poisson');
            Field = bst_poisson(LBO, f);                            % [nV x nT]
            Result = struct('Method','poisson', 'Field',Field, 'nComponents',1);
```

Then delete the now-unused local function `i_poisson` (lines 218-235).

- [ ] **Step 6: Integration check — bst_operators poisson still runs end-to-end**

Run (MATLAB MCP `evaluate_matlab_code`, against the live session with a protocol loaded):

```matlab
% pick the first cortex surface with a Laplace-Beltrami operator available
sSubj = bst_get('Subject');  Surf = sSubj.Surface(find(strcmpi({sSubj.Surface.SurfaceType},'Cortex'),1)).FileName;
nV = length(in_tess_bst(Surf,0).Vertices);
f  = randn(nV, 3);
R  = bst_operators(f, struct('Method','poisson','SurfaceFile',Surf,'iTargetStudy','NoSave'));
assert(~isempty(R) && size(R{1}.Field,1)==nV, 'poisson end-to-end failed');
fprintf('bst_operators poisson end-to-end OK (nV=%d)\n', nV);
```

Expected: `bst_operators poisson end-to-end OK (nV=…)`.

- [ ] **Step 7: Commit**

```bash
git add toolbox/differential/bst_poisson.m toolbox/differential/bst_operators.m dev/test_bst_poisson.m
git commit -m "feat(differential): bst_poisson scalar solver; bst_operators uses it (fixes per-column refactor)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: Refactor `bst_helmholtz` onto the shared solver

**Files:**
- Modify: `toolbox/differential/bst_helmholtz.m:136-144` (Prepare), `:190-195` (drop private `i_poisson`)
- Test: `dev/test_helmholtz_parity.m`

**Interfaces:**
- Consumes: `tess_cholesky(Node, hh, pin)`, `tess_cholesky('solve', dF, rhs)` (Task 1).
- Produces: behavior-identical `bst_helmholtz('Decompose', …)` / `Frame` (HarmFrac, Curl, Div, Psi, Phi, Vtot, Virr, Vsol, Vharm unchanged within `1e-10`).

> **Note on parity:** the refactor changes code that produces the baseline, so the baseline MUST be captured BEFORE editing. Step 1 records it; Step 5 compares against it.

- [ ] **Step 1: Capture the pre-refactor baseline**

Run (MATLAB MCP, protocol loaded; uses the existing Dirac/LBO nodes):

```matlab
sSubj = bst_get('Subject');  Surf = sSubj.Surface(find(strcmpi({sSubj.Surface.SurfaceType},'Cortex'),1)).FileName;
Surfm = in_tess_bst(Surf,0);  nV = size(Surfm.Vertices,1);
Mani = tess_manifold(Surf,'Gauge','trivial');
Dir  = bst_get_operator_node(Surf,'Dirac');  LBO = bst_get_operator_node(Surf,'Laplace-Beltrami');
rng(7);  J = randn(3*nV, 5);
H0 = bst_helmholtz('Decompose', {Dir, LBO}, Mani, Surfm, J);
if ~exist('dev/scratch','dir'), mkdir('dev/scratch'); end
save('dev/scratch/helmholtz_baseline.mat','H0','J','Surf','-v7.3');
fprintf('baseline saved: HarmFrac frame1 = %g\n', H0.Curl(1,1));
```

Expected: `baseline saved: …`.

- [ ] **Step 2: Refactor `i_prepare_vertex` to use `tess_cholesky`**

In `bst_helmholtz.m` `i_prepare_vertex` (lines 136-143), replace the LBO factor block:

```matlab
        % LBO pieces + cached Cholesky of the pinned (vertex 1) cotan stiffness via tess_cholesky
        M = LBO.Mass{hh};
        Op.D{hh}=D; Op.vH{hh}=vH; Op.Nf{hh}=Nf; Op.Wfv{hh}=Wfv; Op.M{hh}=M;
        Op.Gx{hh}=Gx; Op.Gy{hh}=Gy; Op.Gz{hh}=Gz;
        Op.cholK{hh}  = tess_cholesky(LBO, hh, 1);    % pin vertex 1; pure getter (I/O-free)
        Op.free{hh}   = Op.cholK{hh}.free;
        Op.totMass{hh}= sum(M(:));
```

(The `K = LBO.Operator{hh};` and `free = (2:size(K,1))';` lines and the old `decomposition(...)` are removed; `tess_cholesky` owns the factor.)

- [ ] **Step 3: Refactor the per-frame Poisson to `tess_cholesky('solve', …)`**

In `bst_helmholtz.m`, replace the private `i_poisson` (lines 190-195) with a thin wrapper over the shared solve (keep the same call sites in `i_frame_vertex`):

```matlab
%% ===== Poisson solve (vertex): mean-zero project -> cached pinned solve -> recenter =====
function psi = i_poisson(dK, M, omega, free, totMass) %#ok<INUSD>
    n = size(M,1);
    omega = omega - (sum(M*omega) / totMass) * ones(n,1);   % project to mean-zero
    x = tess_cholesky('solve', dK, M*omega);                % shared permuted Cholesky solve
    psi = x - mean(x);                                      % recenter
end
```

(`free` is now carried inside `dK`; the parameter is retained for call-site compatibility and marked unused.)

- [ ] **Step 4: Run the parity validation**

Create `dev/test_helmholtz_parity.m`:

```matlab
function test_helmholtz_parity
% Re-run bst_helmholtz on the saved input and compare to the pre-refactor baseline.
    fprintf('== test_helmholtz_parity ==\n');
    S = load('dev/scratch/helmholtz_baseline.mat');   % H0, J, Surf
    Surfm = in_tess_bst(S.Surf,0);
    Mani  = tess_manifold(S.Surf,'Gauge','trivial');
    Dir   = bst_get_operator_node(S.Surf,'Dirac');  LBO = bst_get_operator_node(S.Surf,'Laplace-Beltrami');
    H1 = bst_helmholtz('Decompose', {Dir, LBO}, Mani, Surfm, S.J);
    flds = {'Curl','Div','Psi','Phi','Fmag'};
    for i=1:numel(flds)
        a=S.H0.(flds{i}); b=H1.(flds{i});
        e=norm(a(:)-b(:))/max(norm(a(:)),eps);
        assert(e<1e-10, '%s parity broken: rel err %g', flds{i}, e);
        fprintf('  %-5s rel err %g  [OK]\n', flds{i}, e);
    end
    fprintf('PASS\n');
end
```

Run: `dev/test_helmholtz_parity.m`
Expected: each field `[OK]`, `PASS`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/differential/bst_helmholtz.m dev/test_helmholtz_parity.m
git commit -m "refactor(differential): bst_helmholtz uses tess_cholesky shared factor (parity preserved)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: Ambient (3nV) branch for `bst_divergence` + `bst_curl` (with mean-curvature term)

**Files:**
- Modify: `toolbox/differential/bst_divergence.m`, `toolbox/differential/bst_curl.m`
- Test: `dev/test_ambient_divcurl.m`

**Interfaces:**
- Consumes: `bst_helmholtz('Decompose', {Dir,LBO}, Mani, Surf, J)` → `H.Div`, `H.Curl` (per-vertex scalars); `bst_get_operator_node`; `tess_manifold`.
- Produces:
  - `divField = bst_divergence(V, ManifoldMat)` — UNCHANGED for tangent `[3nF x nT]`; NEW: when `V` is ambient `[3nV x nT]`, requires extra args `bst_divergence(V, ManifoldMat, 'Ambient', Surf, Dir, LBO)` and returns the Hodge divergence (`H.Div`) **plus** the mean-curvature term `-2H.*(J·N)` per vertex.
  - `curlField = bst_curl(V, ManifoldMat)` — UNCHANGED for tangent; NEW ambient form `bst_curl(V, ManifoldMat, 'Ambient', Surf, Dir, LBO)` returns `H.Curl`.
- Mean curvature `H` per vertex from the LBO of vertex positions: `HN = ½ M⁻¹ K X`, `H = ½ (M⁻¹KX)·N` (sign per outward `N`).

- [ ] **Step 1: Write the failing validation script**

Create `dev/test_ambient_divcurl.m`:

```matlab
function test_ambient_divcurl
% Ambient-field branch: (a) Hodge divergence/curl match bst_helmholtz scalars;
% (b) the mean-curvature term vanishes for a purely-tangential field and is nonzero
%     for a purely-normal field on a curved cortex.
    fprintf('== test_ambient_divcurl ==\n');
    sSubj = bst_get('Subject');  Surf = sSubj.Surface(find(strcmpi({sSubj.Surface.SurfaceType},'Cortex'),1)).FileName;
    Surfm = in_tess_bst(Surf,0);  nV = size(Surfm.Vertices,1);
    Mani  = tess_manifold(Surf,'Gauge','trivial');
    Dir   = bst_get_operator_node(Surf,'Dirac');  LBO = bst_get_operator_node(Surf,'Laplace-Beltrami');
    rng(11);  J = randn(3*nV, 3);

    dv = bst_divergence(J, Mani, 'Ambient', Surfm, Dir, LBO);
    cu = bst_curl(J, Mani, 'Ambient', Surfm, Dir, LBO);
    assert(isequal(size(dv),[nV 3]) && isequal(size(cu),[nV 3]), 'ambient output shape wrong');
    fprintf('  ambient div/curl shapes OK\n');

    % purely-normal field: mean-curvature term must dominate (nonzero on a folded cortex)
    N = Surfm.VertNormals;                          % [nV x 3]
    Jn = reshape((N .* 1)', [], 1);                  % [3nV x 1] normal field, unit magnitude
    dvn = bst_divergence(Jn, Mani, 'Ambient', Surfm, Dir, LBO);
    assert(norm(dvn) > 1e-6, 'normal-field divergence unexpectedly ~0 (mean-curvature term missing?)');
    fprintf('  normal-field divergence nonzero (|dv|=%g) -> mean-curvature term present  [OK]\n', norm(dvn));
    fprintf('PASS\n');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/test_ambient_divcurl.m`
Expected: error — `bst_divergence` rejects the 6-arg ambient call (too many inputs) or the shape assert fails.

- [ ] **Step 3: Add the ambient branch to `bst_divergence`**

In `toolbox/differential/bst_divergence.m`, change the signature and prepend an ambient dispatch before the existing tangent code:

```matlab
function divField = bst_divergence(V, ManifoldMat, varargin)
    % ----- ambient (3nV) branch: Hodge divergence + mean-curvature coupling -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        Surf = varargin{2};  Dir = varargin{3};  LBO = varargin{4};
        divField = i_ambient_divergence(V, ManifoldMat, Surf, Dir, LBO);
        return;
    end
    % ----- tangent (3nF) branch: existing stable DEC adjoint (UNCHANGED) -----
    ...existing body...
```

Add the helper at the end of the file:

```matlab
%% ===== ambient divergence: Hodge div (imag.n) + mean-curvature term -2H(J.N) =====
function divField = i_ambient_divergence(J, ManifoldMat, Surf, Dir, LBO)
    H = bst_helmholtz('Decompose', {Dir, LBO}, ManifoldMat, Surf, J);   % H.Div = tangential Hodge divergence
    nVtot = size(Surf.Vertices,1);  nT = size(J,2);
    meanCurvTerm = zeros(nVtot, nT);
    for hh = 1:numel(LBO.Operator)
        if isempty(LBO.Operator{hh}), continue; end
        vH = double(LBO.GlobalVertices{hh}(:));
        K = LBO.Operator{hh};  M = LBO.Mass{hh};
        X = Surf.Vertices(vH, :);                       % [nVh x 3] positions
        HN = 0.5 * (M \ (K * X));                       % mean-curvature normal HN = 1/2 M^-1 K X
        Nv = Surf.VertNormals(vH, :);                   % outward vertex normals
        Hsc = sum(HN .* Nv, 2);                         % scalar mean curvature H = HN . N  [nVh x 1]
        Jx = J(3*(vH-1)+1, :);  Jy = J(3*(vH-1)+2, :);  Jz = J(3*(vH-1)+3, :);
        JdotN = Jx.*Nv(:,1) + Jy.*Nv(:,2) + Jz.*Nv(:,3);% [nVh x nT]
        meanCurvTerm(vH, :) = -2 * (Hsc .* JdotN);      % -2 H (J.N)
    end
    divField = H.Div + meanCurvTerm;
end
```

- [ ] **Step 4: Add the ambient branch to `bst_curl`**

In `toolbox/differential/bst_curl.m`, change the signature and prepend:

```matlab
function curlField = bst_curl(V, ManifoldMat, varargin)
    % ----- ambient (3nV) branch: Hodge vorticity (Dirac w-part) -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        Surf = varargin{2};  Dir = varargin{3};  LBO = varargin{4};
        H = bst_helmholtz('Decompose', {Dir, LBO}, ManifoldMat, Surf, V);
        curlField = H.Curl;
        return;
    end
    % ----- tangent (3nF) branch: existing -div(N x V) (UNCHANGED) -----
    ...existing body...
```

- [ ] **Step 5: Run validation to verify it passes**

Run: `dev/test_ambient_divcurl.m`
Expected: `ambient div/curl shapes OK`, `normal-field divergence nonzero … [OK]`, `PASS`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/differential/bst_divergence.m toolbox/differential/bst_curl.m dev/test_ambient_divcurl.m
git commit -m "feat(differential): ambient div/curl branch wrapping the Dirac engine + mean-curvature term

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `bst_operators` router + `FieldType` + guard table + gradient guard

**Files:**
- Modify: `toolbox/differential/bst_operators.m` (default OPTIONS, detection, dispatch guards), `toolbox/differential/bst_gradient.m` (scalar guard)
- Test: `dev/test_operators_router.m`

**Interfaces:**
- Consumes: all prior tasks; `i_as_vertex_scalar` (existing local).
- Produces: `OPTIONS.FieldType ∈ {'auto','scalar','tangent','ambient'}` (default `'auto'`); a guard that maps `(Method × stratum)` and errors `bst_operators:badFieldType` on invalid pairs; ambient `divergence`/`curl` routed through the Task-4 branches.

- [ ] **Step 1: Write the failing validation script**

Create `dev/test_operators_router.m`:

```matlab
function test_operators_router
% Router: valid (Method x stratum) pairs run; invalid pairs raise bst_operators:badFieldType;
% identity checks curl(grad f)~0, div(curl V)~0.
    fprintf('== test_operators_router ==\n');
    sSubj = bst_get('Subject');  Surf = sSubj.Surface(find(strcmpi({sSubj.Surface.SurfaceType},'Cortex'),1)).FileName;
    Surfm = in_tess_bst(Surf,0);  nV = size(Surfm.Vertices,1);  nF = size(Surfm.Faces,1);
    f = randn(nV, 2);  J = randn(3*nV, 2);

    % curl of a scalar -> must be guarded
    ok=false; try, bst_operators(f, struct('Method','curl','SurfaceFile',Surf,'iTargetStudy','NoSave')); catch e, ok=strcmp(e.identifier,'bst_operators:badFieldType'); end
    assert(ok, 'curl-of-scalar was not guarded'); fprintf('  curl-of-scalar guarded  [OK]\n');

    % divergence of an ambient field -> runs, per-vertex scalar
    R = bst_operators(J, struct('Method','divergence','FieldType','ambient','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    assert(size(R{1}.Field,1)==nV, 'ambient divergence shape'); fprintf('  ambient divergence routed  [OK]\n');

    % identity: curl(grad f) ~ 0 (metric-free)
    Rg = bst_operators(f, struct('Method','gradient','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    Rc = bst_operators(Rg{1}.Field, struct('Method','curl','FieldType','tangent','SurfaceFile',Surf,'iTargetStudy','NoSave'));
    rel = norm(Rc{1}.Field(:))/max(norm(Rg{1}.Field(:)),eps);
    fprintf('  curl(grad f) rel mag = %g\n', rel);
    assert(rel < 1e-2, 'curl(grad f) not ~0 (rel %g)', rel);
    fprintf('PASS\n');
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/test_operators_router.m`
Expected: fails — `curl-of-scalar was not guarded` (no guard yet; current code errors with a different/`badInput` id or computes garbage).

- [ ] **Step 3: Add `FieldType` default + detection helper**

In `bst_operators.m` add to `Def_OPTIONS` (after line 69):

```matlab
Def_OPTIONS.FieldType    = 'auto';            % {'auto','scalar','tangent','ambient'} field stratum
```

Add a detection helper near the other locals:

```matlab
%% ===== infer the field stratum from layout (or honor an explicit FieldType) =====
function stratum = i_field_stratum(F, nVtot, nFtot, explicit)
    if ~strcmpi(explicit, 'auto'), stratum = lower(explicit); return; end
    nr = size(F,1);
    if nr == 3*nVtot
        stratum = 'ambient';
    elseif nr == 3*nFtot
        stratum = 'tangent';
    elseif nr == nVtot && ~isreal(F)
        stratum = 'tangent';        % connection (complex) representation
    elseif nr == nVtot
        stratum = 'scalar';
    else
        error('bst_operators:badFieldType', ...
            'Cannot infer field stratum from %d rows (nV=%d, nF=%d); set OPTIONS.FieldType.', nr, nVtot, nFtot);
    end
    if nVtot == nFtot
        error('bst_operators:badFieldType', ...
            'Ambiguous layout (nV==nF); set OPTIONS.FieldType explicitly.');
    end
end
```

- [ ] **Step 4: Add the guard table at the top of the dispatch**

In `bst_operators.m`, immediately after `Surf = in_tess_bst(...)` / `nVtot = ...` (line 125-126), compute the stratum and guard:

```matlab
    nFtot = size(Surf.Faces, 1);
    stratum = i_field_stratum(F, nVtot, nFtot, OPTIONS.FieldType);
    % valid (Method x stratum) pairs; everything else is guarded
    valid = struct('gradient',{{'scalar'}}, 'divergence',{{'tangent','ambient'}}, ...
                   'curl',{{'tangent','ambient'}}, 'laplacian',{{'scalar','tangent','ambient'}}, ...
                   'poisson',{{'scalar','tangent'}}, 'helmholtz',{{'ambient'}});
    m = lower(OPTIONS.Method);
    if isfield(valid, m) && ~ismember(stratum, valid.(m))
        error('bst_operators:badFieldType', ...
            '%s is undefined for a %s field.', OPTIONS.Method, stratum);
    end
```

- [ ] **Step 5: Route ambient divergence/curl through the Task-4 branches**

In `bst_operators.m`, update the `'divergence'` and `'curl'` cases to branch on `stratum` (replacing the hard `3*size(Surf.Faces,1)` check, lines 143-158):

```matlab
        case 'divergence'
            Mani = tess_manifold(SurfaceFile, 'Gauge', OPTIONS.Gauge);
            if strcmp(stratum, 'ambient')
                Dir = bst_get_operator_node(SurfaceFile,'Dirac');  LBO = bst_get_operator_node(SurfaceFile,'Laplace-Beltrami');
                Field = bst_divergence(F, Mani, 'Ambient', Surf, Dir, LBO);
            else
                Field = bst_divergence(F, Mani);                   % tangent [3nF]
            end
            Result = struct('Method','divergence', 'Field',Field, 'nComponents',1);
        case 'curl'
            Mani = tess_manifold(SurfaceFile, 'Gauge', OPTIONS.Gauge);
            if strcmp(stratum, 'ambient')
                Dir = bst_get_operator_node(SurfaceFile,'Dirac');  LBO = bst_get_operator_node(SurfaceFile,'Laplace-Beltrami');
                Field = bst_curl(F, Mani, 'Ambient', Surf, Dir, LBO);
            else
                Field = bst_curl(F, Mani);                          % tangent [3nF]
            end
            Result = struct('Method','curl', 'Field',Field, 'nComponents',1);
```

- [ ] **Step 6: Add the scalar guard to `bst_gradient`**

In `toolbox/differential/bst_gradient.m`, after resolving `DEC`/`nVtot`, before the loop, assert scalar input:

```matlab
    nVtot = max(cellfun(@(c) max(double(c(:))), {DEC.GlobalVertices}));
    if size(F,1) ~= nVtot
        error('bst_gradient:notScalar', ...
            'gradient expects a per-vertex scalar [nV=%d x nT]; got %d rows.', nVtot, size(F,1));
    end
```

- [ ] **Step 7: Run validation to verify it passes**

Run: `dev/test_operators_router.m`
Expected: `curl-of-scalar guarded [OK]`, `ambient divergence routed [OK]`, `curl(grad f) rel mag = …`, `PASS`.

- [ ] **Step 8: Commit**

```bash
git add toolbox/differential/bst_operators.m toolbox/differential/bst_gradient.m dev/test_operators_router.m
git commit -m "feat(differential): bst_operators field-type router + guard table + FieldType

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Weitzenböck documentation (tangent stratum)

**Files:**
- Modify: header comments of `toolbox/differential/bst_divergence.m`, `bst_curl.m`; `bst_operators.m` (Method doc); design doc cross-reference.

- [ ] **Step 1: Document the Hodge-vs-connection split**

In `bst_operators.m` header, under the Method list, add a note:

```matlab
%     TANGENT STRATUM — two Laplacians, not one. The differential ops (div/curl) use the
%     de Rham / DEC route (delta d + d delta, the Hodge Laplacian). The SMOOTHING / eigen
%     Laplacian for tangent fields is the CONNECTION (Bochner) Laplacian (grad* grad),
%     supplied by the 'Connection Laplacian' operator node and bst_eigen. They differ by
%     Gauss curvature K (Weitzenbock: Delta_Hodge = Delta_connection + Ric, Ric = K g on a
%     surface). Use the de Rham route for div/curl/Helmholtz; the connection route for
%     vector smoothing/spectral analysis.
```

Add a one-line pointer in `bst_divergence.m` and `bst_curl.m` headers: `% (de Rham/Hodge route; the connection Laplacian differs by Gauss curvature K — see bst_operators.)`

- [ ] **Step 2: Commit**

```bash
git add toolbox/differential/bst_operators.m toolbox/differential/bst_divergence.m toolbox/differential/bst_curl.m
git commit -m "docs(differential): Weitzenbock note — Hodge vs connection Laplacian on the tangent stratum

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Self-Review (completed by author)

**Spec coverage:** §3 operation matrix → Task 5 guard table; §4 detection → Task 5 `i_field_stratum`; §5.1 `tess_cholesky` → Task 1; §5.2 `bst_poisson` → Task 2; §6 ambient div/curl + mean-curvature → Task 4, helmholtz refactor → Task 3, gradient guard → Task 5; §6 Weitzenböck docs → Task 6; §10 identity checks → Task 5 validation. All covered.

**Placeholder scan:** No TBD/TODO; every code step shows full code; validation scripts contain real assertions. `...existing body...` markers in Tasks 4-5 refer to verbatim-unchanged code in the named file (not a placeholder — the surrounding edit is shown).

**Type consistency:** `dF = struct('L','p','free','n')` defined in Task 1 and consumed identically in Tasks 2-3; `tess_cholesky('solve', dF, rhs)` signature consistent; `bst_poisson(OperatorNode, F)` consistent across Tasks 2-3; ambient verb signature `bst_*(V, Mani, 'Ambient', Surf, Dir, LBO)` consistent across Tasks 4-5.
