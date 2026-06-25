# Differential Helmholtz `Compute` Refactor — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `process_helmholtz('Compute', J, Cov)` the single stateless flat-covariant Helmholtz algorithm (div/curl + potentials), with `bst_operators` resolving operators and `bst_divergence`/`bst_curl` owning ambient div/curl, then delete `bst_helmholtz`.

**Architecture:** `bst_operators` (resolver/router) → `process_helmholtz('Compute')` (pure per-frame algorithm) sequences `bst_divergence` + `bst_curl` (ambient, from the Covariant operator) + `bst_poisson` (potentials, factor cached by `tess_cholesky`). All three workflow states (navigate / detect / save) call the same stateless `Compute`, mirroring `process_bandpass('Compute')` + `VisualizationFilters`. Cores and the vortex-detection stack are retired to `dev/experimental`. `view_helmholtz`'s compute is repointed (its deletion is Spec 2).

**Tech Stack:** MATLAB, Brainstorm process framework (`macro_method` dispatch), nxr-compute operators (`Covariant` node), `tess_cholesky` cached factorization.

**Spec:** `docs/superpowers/specs/2026-06-25-differential-helmholtz-compute-refactor-design.md`

## Global Constraints

- **Algorithm content is ported VERBATIM** from the current `bst_helmholtz` flat-covariant
  *vertex* path (`i_prepare_vertex` + `i_frame_vertex`). Do not re-derive the formulas; copy
  them. The layering changes (div→`bst_divergence`, curl→`bst_curl`, recipe→`process_helmholtz`),
  the math does not.
- **No invented process hooks.** Only standard ones: `GetDescription`, `Run`, `Compute`,
  `FormatComment`, `GetOptions`, `GetFileTag`. The expensive Cholesky is cached by
  `tess_cholesky` (node → `persistent MEM` → compute), never a `'Prepare'` verb.
- **Sign convention:** `s = +1` (the flat-covariant divergence is the true divergence; no flip).
- **`HarmFrac` guard:** keep `if harmDen > 0 ... else 0`. Never `max(harmDen, eps)` — source
  amplitudes ~1e-11 make `harmDen ~1e-22 < eps`; the clamp corrupts the ratio.
- **Cores are removed, not relocated.** Do not call `bst_vortex_persistence`.
- **MATLAB env bootstrap** (per session, headless): `addpath('<repo root>'); brainstorm nogui;`
  then select the protocol holding the prepared cortex (`bst_set('iProtocol', N)`).
- **Test cortex:** a cortex with `manifold_` + `Covariant` operator nodes prepared. The dev
  harness uses `Subject01/tess_cortex_pial_low.mat` (20484 V). Each test names it in one place.
- **License header:** every new `.m` file carries the standard Brainstorm GPL header block
  (copy from `process_bandpass.m`), authored "Diellor Basha, 2026".

---

## File structure

| File | Action | Responsibility |
|---|---|---|
| `dev/tests/baselines/helmholtz_baseline.m` | Create | Capture old `bst_helmholtz` outputs as the parity oracle (run ONCE before edits) |
| `toolbox/math/bst_face2vertex.m` | Create | Shared area-weighted face→vertex averaging map `Wfv` (3 cross-module callers → promoted per the no-scatter rule) |
| `toolbox/differential/bst_poisson.m` | Modify | Accept the `'Covariant'` node (its `Operator` is the cotan stiffness) |
| `toolbox/differential/bst_divergence.m` | Modify | Ambient branch computes divergence directly from the Covariant node; drop `LBO` |
| `toolbox/differential/bst_curl.m` | Modify | Ambient branch computes vorticity directly from the Covariant node; drop `LBO` |
| `toolbox/process/functions/process_helmholtz.m` | Create | `Compute` = flat-covariant recipe; `Run` = save; standard skeleton |
| `toolbox/differential/bst_operators.m` | Modify | `'helmholtz'`/`'divergence'`/`'curl'` Methods route to the new engines; drop `LBO` |
| `toolbox/process/functions/process_helmholtz_events.m` | Modify | Use `process_helmholtz('Compute')` instead of `bst_helmholtz` |
| `toolbox/gui/view_helmholtz.m` | Modify | Per-frame compute → `process_helmholtz('Compute')`; resolve `Cov` once. NOT deleted. |
| `dev/experimental/` | Create + move | Receives the retired vortex-detection stack |
| `toolbox/differential/bst_helmholtz.m` | Delete | Fully superseded |
| `dev/test_ambient_divcurl.m`, `dev/test_helmholtz_covariant.m`, `dev/test_helmholtz_parity.m` | Modify | Update to the new API |

---

### Task 1: Capture the parity baseline (oracle)

The old `bst_helmholtz` is the source of truth. Capture its outputs on a fixed synthetic
field **before any edit**, so every later task asserts equality against it.

**Files:**
- Create: `dev/tests/baselines/helmholtz_baseline.m`

**Interfaces:**
- Produces: `dev/tests/baselines/helmholtz_baseline.mat` with struct `B` holding fields
  `Surf` (char), `J` `[3nV x 1]`, and per-vertex `Div, Curl, Phi, Psi, Fmag, Hmag`
  `[nV x 1]`, `Virr, Vsol, Vtot, Hresid` `[nV x 3]`, scalar `HarmFrac`.

- [ ] **Step 1: Write the capture script**

```matlab
function helmholtz_baseline
% Capture old bst_helmholtz('Frame') outputs as the refactor parity oracle.
% Run ONCE on the unmodified tree, BEFORE editing any differential/ file.
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';   % cortex w/ manifold_ + Covariant ready
    Surf = in_tess_bst(SurfaceFile, 0);
    nV   = size(Surf.Vertices, 1);
    Cov  = tess_operators(SurfaceFile, 'Covariant');
    LBO  = tess_operators(SurfaceFile, 'Laplace-Beltrami');
    Mani = tess_manifold(SurfaceFile);
    rng(7);  J = randn(3*nV, 1);                          % fixed-seed synthetic source frame

    Op = bst_helmholtz('Prepare', {Cov, LBO}, Mani, Surf, 'Domain', 'vertex');
    Ht = bst_helmholtz('Frame', Op, J, false);            % cores off

    B = struct('Surf', SurfaceFile, 'J', J, ...
        'Div', Ht.Div, 'Curl', Ht.Curl, 'Phi', Ht.Phi, 'Psi', Ht.Psi, ...
        'Fmag', Ht.Fmag, 'Hmag', Ht.Hmag, 'Virr', Ht.Virr, 'Vsol', Ht.Vsol, ...
        'Vtot', Ht.Vtot, 'Hresid', Ht.Hresid, 'HarmFrac', Ht.HarmFrac); %#ok<NASGU>
    outDir = fileparts(mfilename('fullpath'));
    save(fullfile(outDir, 'helmholtz_baseline.mat'), 'B');
    fprintf('helmholtz_baseline: saved %d-vertex baseline (HarmFrac=%.3e)\n', nV, Ht.HarmFrac);
end
```

- [ ] **Step 2: Run it on the unmodified tree**

Run (MATLAB, with Brainstorm started + protocol selected):
`run_matlab_file dev/tests/baselines/helmholtz_baseline.m`
Expected: prints `helmholtz_baseline: saved 20484-vertex baseline (HarmFrac=…)`, creates
`dev/tests/baselines/helmholtz_baseline.mat`.

- [ ] **Step 3: Commit the baseline**

```bash
git add dev/tests/baselines/helmholtz_baseline.m dev/tests/baselines/helmholtz_baseline.mat
git commit -m "test(differential): capture bst_helmholtz parity baseline before refactor"
```

---

### Task 2: `bst_poisson` accepts the `'Covariant'` node

**Files:**
- Modify: `toolbox/differential/bst_poisson.m:39-42` (the variant guard) + the header (`:11`).
- Test: `dev/tests/test_bst_poisson_covariant.m` (create)

**Interfaces:**
- Consumes: a `Covariant` operatormat (`Operator{hh}` = cotan stiffness, `Mass{hh}` = M,
  `GlobalVertices{hh}`).
- Produces: `phi = bst_poisson(OperatorNode, F)` works for `Variant ∈ {Laplace-Beltrami, Covariant}`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_bst_poisson_covariant
    SurfaceFile = 'Subject01/tess_cortex_pial_low.mat';
    Cov = tess_operators(SurfaceFile, 'Covariant');
    nV  = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    rng(3);  f = randn(nV, 1);
    phi = bst_poisson(Cov, f);                       % must not error on 'Covariant'
    assert(isequal(size(phi), [nV 1]), 'phi shape wrong');
    % K phi = M f residual on the free block, per hemisphere
    for hh = 1:numel(Cov.Operator)
        vH = double(Cov.GlobalVertices{hh}(:));  K = (Cov.Operator{hh}+Cov.Operator{hh}')/2;
        M = Cov.Mass{hh};  free = 2:numel(vH);  ph = phi(vH);  fh = f(vH);
        fh = fh - sum(M*fh)/sum(M(:));  r = K(free,free)*ph(free) - M(free,:)*fh;
        assert(norm(r)/max(norm(M(free,:)*fh),eps) < 1e-8, 'Poisson residual too large (hemi %d)', hh);
    end
    fprintf('PASS test_bst_poisson_covariant\n');
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `run_matlab_file dev/tests/test_bst_poisson_covariant.m`
Expected: FAIL with `bst_poisson:variant` ("needs a Laplace-Beltrami operator (got Covariant)").

- [ ] **Step 3: Widen the variant guard**

In `bst_poisson.m` replace the guard:

```matlab
    if ~any(strcmpi(OperatorNode.Variant, {'Laplace-Beltrami', 'Covariant'}))
        error('bst_poisson:variant', ...
            'bst_poisson scalar route needs a Laplace-Beltrami or Covariant operator (got %s).', OperatorNode.Variant);
    end
```

And update the header line documenting the operand:
```matlab
%   OperatorNode : a 'Laplace-Beltrami' or 'Covariant' operatormat (Operator{hh}=K cotan stiffness, Mass{hh}=M, GlobalVertices{hh})
```

- [ ] **Step 4: Run to verify it passes**

Run: `run_matlab_file dev/tests/test_bst_poisson_covariant.m`
Expected: `PASS test_bst_poisson_covariant`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/differential/bst_poisson.m dev/tests/test_bst_poisson_covariant.m
git commit -m "feat(differential): bst_poisson accepts the Covariant operator node"
```

---

### Task 3: `bst_divergence` ambient branch computes divergence directly

Replace the delegation to `bst_helmholtz('Decompose')` with the flat-covariant strong
divergence ported from `bst_helmholtz` `i_prepare_vertex`/`i_frame_vertex`. Drop the `LBO` operand.

**Files:**
- Modify: `toolbox/differential/bst_divergence.m` (ambient branch + `i_ambient_divergence` +
  header signature lines).
- Test: `dev/tests/test_bst_divergence_ambient.m` (create)

**Interfaces:**
- Consumes: `Cov` (Covariant node), `J` `[3nV x nT]`, `helmholtz_baseline.mat`.
- Produces: `divField = bst_divergence(J, ManifoldMat, 'Ambient', Surf, Cov)` `[nV x nT]`
  (`ManifoldMat`/`Surf` accepted for signature symmetry with the tangent branch but unused;
  `nVtot` derived from `Cov.GlobalVertices`).

- [ ] **Step 1: Write the failing test**

```matlab
function test_bst_divergence_ambient
    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    div = bst_divergence(B.J, [], 'Ambient', [], Cov);
    assert(isequal(size(div), size(B.Div)), 'div shape mismatch');
    rel = norm(div - B.Div) / max(norm(B.Div), eps);
    assert(rel < 1e-10, 'ambient divergence differs from baseline (rel=%.2e)', rel);
    fprintf('PASS test_bst_divergence_ambient (rel=%.2e)\n', rel);
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `run_matlab_file dev/tests/test_bst_divergence_ambient.m`
Expected: FAIL (old ambient branch needs `Surf`/`Dir`/`LBO`; called with `[]` it errors or
returns wrong values).

- [ ] **Step 3: Implement the ambient branch**

In `bst_divergence.m`, replace the ambient dispatch block:

```matlab
    % ----- ambient (3nV) branch: flat-covariant surface divergence (Covariant node) -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        % USAGE: bst_divergence(J, ManifoldMat, 'Ambient', Surf, Cov)
        % ManifoldMat/Surf are accepted for signature symmetry with the tangent branch
        % but are unused here: the flat-covariant divergence needs only the Covariant node.
        Cov = varargin{end};
        divField = i_ambient_divergence(V, Cov);
        return;
    end
```

Replace the old `i_ambient_divergence` with the ported strong divergence:

```matlab
%% ===== ambient divergence: flat-covariant surface divergence (incl. -2H(J.N) coupling) =====
% Ported from bst_helmholtz i_prepare_vertex/i_frame_vertex. Per hemisphere, from the
% Covariant node: strong per-face divergence Gx*Jx+Gy*Jy+Gz*Jz, area-weighted to vertices by
% Wfv. s=+1 (calibrated: this IS the true surface divergence, already includes the
% mean-curvature coupling -2H(J.N); a constant ambient field gives 0 even on folds).
function divField = i_ambient_divergence(J, Cov)
    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    nT = size(J, 2);
    divField = zeros(nVtot, nT);
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Wfv = bst_face2vertex(C.Faces, C.FaceArea);   % shared math helper
        Jx = J(3*(vH-1)+1, :);  Jy = J(3*(vH-1)+2, :);  Jz = J(3*(vH-1)+3, :);
        divF = Gx*Jx + Gy*Jy + Gz*Jz;                 % [nFh x nT] per-face surface divergence
        divField(vH, :) = s * (Wfv * divF);
    end
end
```

First create the shared helper (one home, three callers — `bst_divergence`, `bst_curl`,
`process_helmholtz`). Create `toolbox/math/bst_face2vertex.m` (with the GPL header):

```matlab
function Wfv = bst_face2vertex(Faces, FaceArea)
% BST_FACE2VERTEX: Area-weighted face->vertex averaging map Wfv [nVh x nFh].
% Wfv * x_face gives, per vertex, the area-weighted mean of its incident faces' values.
% Ported from the bst_helmholtz flat-covariant decomposition; shared by bst_divergence,
% bst_curl, and process_helmholtz.
%
% USAGE:  Wfv = bst_face2vertex(Faces, FaceArea)
%   Faces    : [nFh x 3] LOCAL vertex indices (1..nVh) per face
%   FaceArea : [nFh x 1] face areas
%
% Authors: Diellor Basha, 2026
% <GPL header block here>
    nFh = size(Faces, 1);  nVh = max(Faces(:));
    I3 = [Faces(:,1); Faces(:,2); Faces(:,3)];  J3 = [(1:nFh)'; (1:nFh)'; (1:nFh)'];
    Wfv = sparse(I3, J3, repmat(FaceArea, 3, 1), nVh, nFh);
    Wfv = spdiags(1 ./ max(sum(Wfv, 2), eps), 0, nVh, nVh) * Wfv;
end
```

Update the `bst_divergence` header USAGE/INPUTS lines: drop `LBO`; `'Ambient'` now takes
`…, Surf, Cov`.

- [ ] **Step 4: Run to verify it passes**

Run: `run_matlab_file dev/tests/test_bst_divergence_ambient.m`
Expected: `PASS test_bst_divergence_ambient (rel=…)` with rel < 1e-10.

- [ ] **Step 5: Commit**

```bash
git add toolbox/differential/bst_divergence.m dev/tests/test_bst_divergence_ambient.m
git commit -m "feat(differential): bst_divergence ambient branch computes divergence from the Covariant node"
```

---

### Task 4: `bst_curl` ambient branch computes vorticity directly

**Files:**
- Modify: `toolbox/differential/bst_curl.m` (ambient branch + header).
- Test: `dev/tests/test_bst_curl_ambient.m` (create)

**Interfaces:**
- Consumes: `Cov`, `J`, `helmholtz_baseline.mat`, and `i_face_to_vertex` (local copy — `bst_curl`
  is a separate file; duplicate the small helper as a local function).
- Produces: `curlField = bst_curl(J, ManifoldMat, 'Ambient', Surf, Cov)` `[nV x nT]`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_bst_curl_ambient
    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    cu = bst_curl(B.J, [], 'Ambient', [], Cov);
    assert(isequal(size(cu), size(B.Curl)), 'curl shape mismatch');
    rel = norm(cu - B.Curl) / max(norm(B.Curl), eps);
    assert(rel < 1e-10, 'ambient vorticity differs from baseline (rel=%.2e)', rel);
    fprintf('PASS test_bst_curl_ambient (rel=%.2e)\n', rel);
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `run_matlab_file dev/tests/test_bst_curl_ambient.m`
Expected: FAIL (old ambient branch delegates to `bst_helmholtz` with `Dir`/`LBO`).

- [ ] **Step 3: Implement the ambient branch**

In `bst_curl.m`, replace the ambient dispatch block:

```matlab
    % ----- ambient (3nV) branch: flat-covariant vorticity (Covariant node) -----
    if ~isempty(varargin) && strcmpi(varargin{1}, 'Ambient')
        % USAGE: bst_curl(J, ManifoldMat, 'Ambient', Surf, Cov)  (ManifoldMat/Surf unused here)
        Cov = varargin{end};
        curlField = i_ambient_curl(V, Cov);
        return;
    end
```

Add the ported vorticity + the local `i_face_to_vertex` helper:

```matlab
%% ===== ambient vorticity: flat-covariant curl.n (Covariant node) =====
% Ported from bst_helmholtz i_frame_vertex: per-face curl vector, projected on the face normal,
% area-weighted to vertices. s=+1. Vorticity is intrinsic (no mean-curvature term).
function curlField = i_ambient_curl(J, Cov)
    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    nT = size(J, 2);
    curlField = zeros(nVtot, nT);
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);  Nf = C.FaceNormal;
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Wfv = bst_face2vertex(C.Faces, C.FaceArea);   % shared math helper (Task 3)
        Jx = J(3*(vH-1)+1, :);  Jy = J(3*(vH-1)+2, :);  Jz = J(3*(vH-1)+3, :);
        omF = zeros(nFh, nT);
        for t = 1:nT
            cv = [Gy*Jz(:,t) - Gz*Jy(:,t), Gz*Jx(:,t) - Gx*Jz(:,t), Gx*Jy(:,t) - Gy*Jx(:,t)];
            omF(:, t) = sum(cv .* Nf, 2);             % vorticity = curl . n
        end
        curlField(vH, :) = s * (Wfv * omF);
    end
end
```

Update the `bst_curl` header USAGE/INPUTS lines: drop `LBO`; `'Ambient'` takes `…, Surf, Cov`.

- [ ] **Step 4: Run to verify it passes**

Run: `run_matlab_file dev/tests/test_bst_curl_ambient.m`
Expected: `PASS test_bst_curl_ambient (rel=…)` with rel < 1e-10.

- [ ] **Step 5: Commit**

```bash
git add toolbox/differential/bst_curl.m dev/tests/test_bst_curl_ambient.m
git commit -m "feat(differential): bst_curl ambient branch computes vorticity from the Covariant node"
```

---

### Task 5: `process_helmholtz('Compute')` — the flat-covariant recipe

The pure stateless algorithm. Sequences `bst_divergence` + `bst_curl` for the strong Div/Curl,
then ports the weak-Hodge potential solve + reconstruction from `bst_helmholtz` `i_frame_vertex`,
using `bst_poisson`/`tess_cholesky` for the cached factor.

**Files:**
- Create: `toolbox/process/functions/process_helmholtz.m` (skeleton + `Compute`; `Run`/options
  in Task 6).
- Test: `dev/tests/test_process_helmholtz_compute.m` (create)

**Interfaces:**
- Consumes: `Cov`, `J` `[3nV x nT]`, the Task 3/4 engines, `bst_poisson`, `tess_cholesky`.
- Produces: `Ht = process_helmholtz('Compute', J, Cov)` — struct with `Div, Curl, Phi, Psi,
  Fmag, Hmag` `[nV x nT]`, `Virr, Vsol, Vtot, Hresid` `[nV x 3]` (single-frame), `HarmFrac` scalar.
  Multi-column `J`: scalar/vector fields are `[nV x nT]`; `Virr/Vsol/Vtot/Hresid` and `HarmFrac`
  are returned for the LAST column (matches the single-frame `Frame` contract used by callers).

- [ ] **Step 1: Write the failing test**

```matlab
function test_process_helmholtz_compute
    d = load(fullfile(fileparts(mfilename('fullpath')), 'baselines', 'helmholtz_baseline.mat'));
    B = d.B;  Cov = tess_operators(B.Surf, 'Covariant');
    Ht = process_helmholtz('Compute', B.J, Cov);
    chk = @(f) norm(Ht.(f)(:) - B.(f)(:)) / max(norm(B.(f)(:)), eps);
    for f = {'Div','Curl','Phi','Psi','Fmag','Hmag','Virr','Vsol','Vtot','Hresid'}
        rel = chk(f{1});
        assert(rel < 1e-9, '%s differs from baseline (rel=%.2e)', f{1}, rel);
    end
    assert(abs(Ht.HarmFrac - B.HarmFrac) < 1e-9*max(abs(B.HarmFrac),eps), 'HarmFrac differs');
    fprintf('PASS test_process_helmholtz_compute\n');
end
```

- [ ] **Step 2: Run to verify it fails**

Run: `run_matlab_file dev/tests/test_process_helmholtz_compute.m`
Expected: FAIL ("Undefined function 'process_helmholtz'").

- [ ] **Step 3: Create the process skeleton + `Compute`**

Create `toolbox/process/functions/process_helmholtz.m` (standard `process_bandpass` skeleton;
`Run`/options are stubs filled in Task 6). Include the GPL header block. Body:

```matlab
function varargout = process_helmholtz( varargin )
% PROCESS_HELMHOLTZ: Helmholtz-Hodge decomposition of a 3-D cortical source vector field.
%
% Compute (pure, stateless, I/O-free) is the flat-covariant algorithm: it sequences
% bst_divergence + bst_curl (ambient div/curl from the Covariant operator) and recovers the
% scalar potentials phi/psi via the weak Hodge solve (bst_poisson, factor cached by
% tess_cholesky), plus the irrotational/solenoidal/harmonic reconstruction. Run loops Compute
% over a source series and saves results maps. Same algorithm for the GUI ephemeral feedback
% (view_helmholtz / panel_bst_dynamics) and the on-file save.
%
% USAGE:  Ht = process_helmholtz('Compute', J, Cov)   % J [3nV x nT], Cov = 'Covariant' node
%         OutputFiles = process_helmholtz('Run', sProcess, sInputs)
%
% Authors: Diellor Basha, 2026
% <GPL header block here>
eval(macro_method);
end

%% ===== EXTERNAL CALL: pure flat-covariant Helmholtz of a 3-D source frame =====
% Ported verbatim from the old bst_helmholtz i_prepare_vertex/i_frame_vertex (vertex domain).
function Ht = Compute(J, Cov)
    s = +1;
    nVtot = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
    nT = size(J, 2);
    z1 = zeros(nVtot, nT);  z3 = zeros(nVtot, 3);
    Ht = struct('Div',z1, 'Curl',z1, 'Phi',z1, 'Psi',z1, 'Fmag',z1, 'Hmag',z1, ...
                'Vtot',z3, 'Virr',z3, 'Vsol',z3, 'Hresid',z3, 'HarmFrac',0);
    % Strong fields from the dedicated engines (single home of div/curl).
    Ht.Div  = bst_divergence(J, [], 'Ambient', [], Cov);
    Ht.Curl = bst_curl(J, [], 'Ambient', [], Cov);
    harmNum = 0;  harmDen = 0;
    for hh = 1:numel(Cov.Covariant)
        C = Cov.Covariant{hh};  vH = double(Cov.GlobalVertices{hh}(:));
        nFh = size(C.Faces, 1);  nVh = numel(vH);
        Gx = C.ScalarGrad(1:nFh,:);  Gy = C.ScalarGrad(nFh+1:2*nFh,:);  Gz = C.ScalarGrad(2*nFh+1:3*nFh,:);
        Nf = C.FaceNormal;  Af = C.FaceArea;  nv = C.Frame.normal;
        W   = spdiags(Af, 0, nFh, nFh);
        Wfv = bst_face2vertex(C.Faces, Af);           % shared math helper (Task 3)
        Fvf = sparse([(1:nFh)';(1:nFh)';(1:nFh)'], [C.Faces(:,1);C.Faces(:,2);C.Faces(:,3)], 1/3, nFh, nVh);
        nx=Nf(:,1); ny=Nf(:,2); nz=Nf(:,3);
        Sx = spdiags(ny,0,nFh,nFh)*Gz - spdiags(nz,0,nFh,nFh)*Gy;
        Sy = spdiags(nz,0,nFh,nFh)*Gx - spdiags(nx,0,nFh,nFh)*Gz;
        Sz = spdiags(nx,0,nFh,nFh)*Gy - spdiags(ny,0,nFh,nFh)*Gx;
        dF = tess_cholesky(Cov, hh, 2);                 % cached pinned factor (pin vertex 1, free=2:nVh)
        for t = 1:nT
            Jx = J(3*(vH-1)+1, t);  Jy = J(3*(vH-1)+2, t);  Jz = J(3*(vH-1)+3, t);  Jv = [Jx Jy Jz];
            Jf = [Fvf*Jx, Fvf*Jy, Fvf*Jz];
            divw  = s * (Gx'*W*Jf(:,1) + Gy'*W*Jf(:,2) + Gz'*W*Jf(:,3));   % weak divergence source
            vortw = s * (Sx'*W*Jf(:,1) + Sy'*W*Jf(:,2) + Sz'*W*Jf(:,3));   % weak vorticity source
            phi = tess_cholesky('solve', dF, divw);   phi = phi - mean(phi);
            psi = tess_cholesky('solve', dF, vortw);  psi = psi - mean(psi);
            Virr = Wfv * (s * [Gx*phi, Gy*phi, Gz*phi]);
            Vsol = Wfv * (s * cross(Nf, [Gx*psi, Gy*psi, Gz*psi], 2));
            Jn   = Jx.*nv(:,1) + Jy.*nv(:,2) + Jz.*nv(:,3);
            Hres = Jv - Virr - Vsol - Jn.*nv;
            Ht.Phi(vH,t)=phi;  Ht.Psi(vH,t)=psi;  Ht.Fmag(vH,t)=sqrt(Jx.^2+Jy.^2+Jz.^2);
            Ht.Hmag(vH,t)=sqrt(sum(Hres.^2,2));
            if t == nT     % last-frame vector fields (single-frame contract)
                Ht.Vtot(vH,:)=Jv;  Ht.Virr(vH,:)=Virr;  Ht.Vsol(vH,:)=Vsol;  Ht.Hresid(vH,:)=Hres;
                av = full(sum(Cov.Mass{hh}, 2));
                harmNum = harmNum + sum(av .* sum(Hres.^2,2));
                harmDen = harmDen + sum(av .* sum(Jv.^2,2));
            end
        end
    end
    if harmDen > 0, Ht.HarmFrac = harmNum / harmDen; else, Ht.HarmFrac = 0; end
end
```

> NB: `tess_cholesky(Cov, hh, 2)` pins vertex 1 (free block `2:nVh`), matching the old
> `i_prepare_vertex` `free = 2:nVh`. The cached factor makes per-frame `Compute` cheap.

- [ ] **Step 4: Run to verify it passes**

Run: `run_matlab_file dev/tests/test_process_helmholtz_compute.m`
Expected: `PASS test_process_helmholtz_compute`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_helmholtz.m dev/tests/test_process_helmholtz_compute.m
git commit -m "feat(differential): process_helmholtz('Compute') flat-covariant Helmholtz algorithm"
```

---

### Task 6: `process_helmholtz` Run / GetDescription / FormatComment (save state)

Add the on-file "Save" path. Reuse the existing event-maps behavior (kernel-link →
band-passed sensors → source frames → maps), now driven by `Compute`.

**Files:**
- Modify: `toolbox/process/functions/process_helmholtz.m` (add `GetDescription`, `Run`,
  `FormatComment`, `GetOptions`, local save/load helpers — port from the current
  `process_helmholtz_events.m`, replacing the `bst_helmholtz` Prepare/Frame calls with
  `process_helmholtz('Compute')`).
- Test: `dev/tests/test_process_helmholtz_run.m` (create — smoke test on the alpha block)

**Interfaces:**
- Consumes: `Compute` (Task 5); `tess_operators(Surf,'Covariant')`.
- Produces: `OutputFiles = process_helmholtz('Run', sProcess, sInputs)` writing the
  J / |J| / Phi / Psi maps (same four as `process_helmholtz_events` today).

- [ ] **Step 1: Write `GetDescription` / `FormatComment` / `GetOptions`**

Model on `process_bandpass` conventions; `Category='Custom'` (3-vector source → scalar maps):

```matlab
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Helmholtz-Hodge decomposition';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 338;
    sProcess.Description = 'https://neuroimage.usc.edu/brainstorm';
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'results'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    % Event group (timepoints) + sensor band (mirror the maps process)
    sProcess.options.eventname.Comment = 'Phase-marker event (timepoints): ';
    sProcess.options.eventname.Type    = 'text';
    sProcess.options.eventname.Value   = 'alpha_peak';
    sProcess.options.freqband.Comment = {'Delta (2-4 Hz)','Theta (4-8 Hz)','Alpha (8-13 Hz)', ...
        'Beta (13-30 Hz)','Gamma (30-60 Hz)','Custom (below)','<B>Bandpass applied to sensors:</B>'; ...
        'delta','theta','alpha','beta','gamma','custom',''};
    sProcess.options.freqband.Type    = 'radio_linelabel';
    sProcess.options.freqband.Value   = 'alpha';
    sProcess.options.freqrange.Comment = 'Custom frequency range:';
    sProcess.options.freqrange.Type    = 'freqrange';
    sProcess.options.freqrange.Value   = {[8, 13], 'Hz', []};
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = ['Helmholtz-Hodge: ', sProcess.options.eventname.Value];
end
```

- [ ] **Step 2: Write `Run` (port from `process_helmholtz_events`, swap engine)**

Copy `process_helmholtz_events`'s `Run` + `i_load_recording` + `i_save_results` into
`process_helmholtz.m`, changing ONLY the operator/engine block:

```matlab
        % ===== RESOLVE OPERATOR (once) + DECOMPOSE AT EACH EVENT =====
        Cov = bst_operators_resolve(R.SurfaceFile);     % see helper below (find-or-create Covariant)
        nV  = max(cellfun(@(c) max(double(c(:))), Cov.GlobalVertices));
        Jvec = zeros(3*nV, nE);  Fmag = zeros(nV, nE);  Phi = zeros(nV, nE);  Psi = zeros(nV, nE);
        for k = 1:nE
            bst_progress('text', sprintf('File %d/%d: Helmholtz frame %d/%d...', iFile, length(sInputs), k, nE));
            J  = R.ImagingKernel * Fbp(:, iT(k));        % [3nV x 1] unconstrained source vector
            Ht = process_helmholtz('Compute', J, Cov);   % the pure algorithm
            Fmag(:,k)=Ht.Fmag;  Phi(:,k)=Ht.Phi;  Psi(:,k)=Ht.Psi;
            Jvec(:,k)=reshape(Ht.Vtot', [], 1);
        end
```

Add the resolver helper (kept local so the process owns one resolution call):

```matlab
%% ===== resolve the Covariant operator for a surface (find-or-create) =====
function Cov = bst_operators_resolve(SurfaceFile)
    Cov = tess_operators(SurfaceFile, 'Covariant');   % find-or-create (cached node)
end
```

> `i_save_results`/`i_load_recording`: copy verbatim from `process_helmholtz_events.m`, change
> the `Function` string and history text to `process_helmholtz`.

- [ ] **Step 3: Write the Run smoke test**

```matlab
function test_process_helmholtz_run
    % Smoke: Run on an alpha kernel-link produces 4 results maps with correct shapes.
    KernelLink = 'Subject01/.../results_...alpha_link.mat';   % a kernel link w/ alpha_peak events
    sProcess = process_helmholtz('GetDescription');
    sProcess.options.eventname.Value = 'alpha_peak';
    sProcess.options.freqband.Value  = 'alpha';
    sInputs = bst_process('GetInputStruct', KernelLink);
    out = process_helmholtz('Run', sProcess, sInputs);
    assert(numel(out) == 4, 'expected 4 maps, got %d', numel(out));
    fprintf('PASS test_process_helmholtz_run (%d maps)\n', numel(out));
end
```

- [ ] **Step 4: Run the smoke test**

Run: `run_matlab_file dev/tests/test_process_helmholtz_run.m`
Expected: `PASS test_process_helmholtz_run (4 maps)`. (If no alpha link exists in the test
protocol, mark this step skipped and rely on Task 7's `process_helmholtz_events` parity instead.)

- [ ] **Step 5: Commit**

```bash
git add toolbox/process/functions/process_helmholtz.m dev/tests/test_process_helmholtz_run.m
git commit -m "feat(differential): process_helmholtz Run/GetDescription save path on Compute"
```

---

### Task 7: Repoint consumers to `process_helmholtz('Compute')`; drop dead `LBO`

**Files:**
- Modify: `toolbox/differential/bst_operators.m:166,175,183-189` (divergence/curl/helmholtz Methods).
- Modify: `toolbox/process/functions/process_helmholtz_events.m:172-193` (engine swap).
- Modify: `toolbox/gui/view_helmholtz.m:51-57,116` (compute repoint; do NOT delete).
- Test: rerun `dev/tests/test_bst_operators.m` (update its helmholtz/div/curl calls).

**Interfaces:**
- Consumes: `process_helmholtz('Compute')`, `bst_divergence`/`bst_curl` ambient (Cov-only).
- Produces: no toolbox/ reference to `bst_helmholtz`.

- [ ] **Step 1: `bst_operators` — drop LBO, route helmholtz to Compute**

`divergence`/`curl` ambient Methods:
```matlab
                Dir = tess_operators(SurfaceFile,'Covariant');
                Field = bst_divergence(F, Mani, 'Ambient', Surf, Dir);     % curl: bst_curl(...)
```
`helmholtz` Method:
```matlab
        case 'helmholtz'
            Cov = tess_operators(SurfaceFile, 'Covariant');
            if size(F,1) ~= 3*nVtot
                Messages = sprintf('bst_operators: helmholtz needs a [3nV x nT] vector field (got %d rows, 3nV=%d).', size(F,1), 3*nVtot);
                isError = 1; break;
            end
            H = process_helmholtz('Compute', F, Cov);
            Result = struct('Method','helmholtz', 'Field',H.Curl, 'nComponents',1, 'Helmholtz',H);
```

- [ ] **Step 2: `process_helmholtz_events` — swap engine**

Replace lines 170-193 (the Prepare/Frame block) with the Task 6 `Cov`-resolve + per-event
`process_helmholtz('Compute', J, Cov)` block (same as Task 6 Step 2). Remove the `LBO`/`Mani`/
`Surf`/`bst_helmholtz` lines.

- [ ] **Step 3: `view_helmholtz` — repoint compute (keep the viewer)**

At resolve time (≈line 51-57) replace the `bst_helmholtz('Prepare', …)` with resolving and
stashing the Covariant node:
```matlab
    Cov = tess_operators(SurfaceFile, 'Covariant');
    setappdata(hFig, 'HelmholtzCov', Cov);    % held for per-frame Compute (no Prepare verb)
```
At the per-frame call (≈line 116) replace `bst_helmholtz('Frame', St.Op, Jt, false)` with:
```matlab
    Cov = getappdata(St.hFig, 'HelmholtzCov');
    Ht = process_helmholtz('Compute', Jt, Cov);  St.Cache(iT) = Ht;
```
Remove the now-unused `LBO`/`Op`/`Mani`/`Surf` plumbing that only fed `bst_helmholtz`.
(Display verbs `SetComponent`/`SetSmoothing`/`UpdateFrame` are untouched. Cores: `view_helmholtz`
no longer gets `Ht.Cores`; gate any core readout behind `isfield(Ht,'Cores')` — Spec 2 restores
detection. If a core-readout call would now error, stub it to empty.)

- [ ] **Step 4: Verify no `bst_helmholtz` references remain in toolbox/**

Run: `grep -rn "bst_helmholtz" toolbox/`
Expected: only doc-comment mentions (SEE ALSO), no live calls. Fix any live call found.

- [ ] **Step 5: Run the operator + events regressions**

Run: `run_matlab_file dev/tests/test_bst_operators.m` (update its ambient div/curl/helmholtz
calls to drop `LBO` / use Covariant first).
Expected: PASS. Also re-run Tasks 3–5 tests (all still PASS).

- [ ] **Step 6: Commit**

```bash
git add toolbox/differential/bst_operators.m toolbox/process/functions/process_helmholtz_events.m toolbox/gui/view_helmholtz.m dev/tests/test_bst_operators.m
git commit -m "refactor(differential): route helmholtz consumers through process_helmholtz('Compute'); drop dead LBO"
```

---

### Task 8: Retire the vortex-detection stack to `dev/experimental/`

**Files:**
- Create dir: `dev/experimental/`
- Move: `toolbox/process/functions/process_vortex_track.m`, `toolbox/math/bst_vortex_track.m`,
  `toolbox/math/bst_vortex_link_step.m`, `toolbox/math/bst_vortex_persistence.m` → `dev/experimental/`.

**Interfaces:**
- Produces: no `toolbox/` reference to `bst_vortex_*` or `process_vortex_track`.

- [ ] **Step 1: Confirm the file locations**

Run: `find toolbox -name "bst_vortex_*.m" -o -name "process_vortex_track.m"`
Expected: the four files above (adjust paths in Step 2 to match).

- [ ] **Step 2: Move them out of the toolbox path**

```bash
mkdir -p dev/experimental
git mv toolbox/process/functions/process_vortex_track.m dev/experimental/
git mv toolbox/math/bst_vortex_track.m               dev/experimental/
git mv toolbox/math/bst_vortex_link_step.m           dev/experimental/
git mv toolbox/math/bst_vortex_persistence.m         dev/experimental/
```

- [ ] **Step 3: Verify nothing in toolbox/ references them**

Run: `grep -rn "bst_vortex_\|process_vortex_track" toolbox/`
Expected: no matches. If any remain (e.g. a tree menu entry for the vortex-track process),
remove that entry.

- [ ] **Step 4: Commit**

```bash
git add -A dev/experimental toolbox/
git commit -m "chore(differential): retire vortex-detection stack to dev/experimental (rebuild later)"
```

---

### Task 9: Delete `bst_helmholtz` and update remaining dev tests

**Files:**
- Delete: `toolbox/differential/bst_helmholtz.m`
- Modify: `dev/test_ambient_divcurl.m`, `dev/test_helmholtz_covariant.m`, `dev/test_helmholtz_parity.m`

**Interfaces:**
- Consumes: nothing references `bst_helmholtz` after Task 7.

- [ ] **Step 1: Delete the engine**

```bash
git rm toolbox/differential/bst_helmholtz.m
```

- [ ] **Step 2: Update the dev tests to the new API**

In each test, replace `bst_helmholtz('Decompose'/'Prepare'/'Frame', …)` with the new calls:
- div/curl: `bst_divergence(J, [], 'Ambient', [], Cov)` / `bst_curl(...)` (Cov = `tess_operators(Surf,'Covariant')`).
- full decomposition: `process_helmholtz('Compute', J, Cov)`.
`dev/test_ambient_divcurl.m`: change `Dir = tess_operators(Surf,'Dirac')` → `Cov =
tess_operators(Surf,'Covariant')` and the two ambient calls accordingly (drop `LBO`).

- [ ] **Step 3: Verify the tree has no dangling opener and the codebase is clean**

Run: `grep -rn "bst_helmholtz" toolbox/ dev/`
Expected: no live calls anywhere (doc mentions only). Fix any straggler.

- [ ] **Step 4: Run the full differential test sweep**

Run each: `test_bst_poisson_covariant`, `test_bst_divergence_ambient`, `test_bst_curl_ambient`,
`test_process_helmholtz_compute`, `test_bst_operators`, `test_ambient_divcurl`,
`test_helmholtz_covariant`.
Expected: all PASS.

- [ ] **Step 5: Commit**

```bash
git add -A
git commit -m "refactor(differential): delete bst_helmholtz; algorithm now lives in process_helmholtz/bst_divergence/bst_curl"
```

---

## Self-review

**Spec coverage:**
- §6.1 ambient divergence → Task 3 ✅ | §6.2 ambient vorticity → Task 4 ✅ | §6.3 drop LBO →
  Tasks 3,4,7 ✅ | §6.4 bst_poisson Covariant → Task 2 ✅ | §6.5 process_helmholtz Compute →
  Task 5 ✅ | §6.5 Run/save → Task 6 ✅ | §6.6 bst_operators route → Task 7 ✅ | §6.7 cores
  removed (not relocated) → no task creates a detector; Task 7 Step 3 stubs the core readout;
  Task 9 deletes the cores code with bst_helmholtz ✅ | §6.8 retirements: bst_helmholtz delete
  → Task 9; view_helmholtz compute-repoint (not deleted) → Task 7; events on Compute → Task 7;
  vortex stack out → Task 8 ✅.
- Parity discipline (capture-before-edit) → Task 1 ✅.

**Placeholder scan:** the only deliberately non-literal item is the `KernelLink` path in Task 6
Step 3 (depends on the live protocol) — flagged with a skip fallback. No `TODO`/"handle edge
cases"/uncoded steps elsewhere.

**Type consistency:** `Compute(J, Cov)` returns the field names used by every consumer
(`Div/Curl/Phi/Psi/Fmag/Hmag/Virr/Vsol/Vtot/Hresid/HarmFrac`), matching the baseline struct
`B` (Task 1) and the old `Ht`. The area-weighted face→vertex map is the single shared
`bst_face2vertex(Faces, FaceArea)` (Task 3), called by `bst_divergence`, `bst_curl`, and
`process_helmholtz` — promoted to `math/` per the no-scatter rule (3 cross-module callers).
`tess_cholesky` pin argument is `2` everywhere (free block `2:nVh`), matching old `free = 2:nVh`.
```
