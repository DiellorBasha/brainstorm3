# bst_helmholtz manifold-geometry consolidation — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make the `manifold_` node the single source of truth for per-face geometry (add `Embedded.face.area` via nxr-compute), then rename `bst_dirac_helmholtz` → `bst_helmholtz` with vertex/face `switch` dispatch reading geometry from the manifold instead of recomputing it.

**Architecture:** Two phases. Phase A adds a face-area accessor to nxr's `facets` Embedded struct (C++ + MEX rebuild + manifold schemaVersion gate). Phase B renames and merges the two Helmholtz files into one domain-dispatching, I/O-free orchestrator that consumes manifold geometry + operator-node operators. Phase A gates Phase B.

**Tech Stack:** C++ (nxr-compute, geometry-central, CMake/ninja), MATLAB MEX, Brainstorm MATLAB toolbox.

## Global Constraints

- nxr-compute repo: `/Users/diellorbasha/workspace/research/code/nxr-compute`. Build: `bash scripts/build.sh Release`. MEX artifact: `build/Release/nxr_compute.mexmaca64`.
- Managed-plugin install dir (BST loads from here): `~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/nxr_compute.mexmaca64`. Back up the existing MEX (`.bak-YYYYMMDD-HHMM`) before overwriting — stale-binary trap.
- Brainstorm dev repo: `/Users/diellorbasha/workspace/research/code/brainstorm3`. Work on branch `feature/bst-helmholtz-manifold-consolidation`.
- `bst_helmholtz` MUST stay I/O-free: no `load`/`bst_get`/`in_bst_*`/`tess_*`/GUI calls. Callers resolve and pass `ManifoldMat` + operator node + `'Domain'`.
- MATLAB execution is via the live MATLAB MCP session; tests that need a protocol use Subject 1 of the loaded protocol (TutorialAuditory).
- nxr rule (from nxr CLAUDE.md): use geometry-central for geometry; do NOT reimplement areas.
- Geometry (normals/areas/centroids/positions) ← manifold node; operators (Dirac `D`, LBO `K`/`M`, `gradFace`) ← operator node.

---

## Phase A — Add `Embedded.face.area` to the manifold node

### Task A1: Expose `face.area` in nxr `facets` (C++ source + MEX)

**Files:**
- Modify: `nxr-compute/include/nxr/facets.h` (FaceView struct, ~line 25-29)
- Modify: `nxr-compute/src/facets.cpp` (~line 409, next to `centroid()`)
- Modify: `nxr-compute/bindings/mex/src/nxr_compute_mex.cpp` (`buildEmbeddedStruct`, ~line 1545-1570)

**Interfaces:**
- Produces: `EmbeddedFacet::FaceView::area() -> Eigen::VectorXd [nF]`; MEX `Embedded.face.area [nF×1] double`; `Embedded.schemaVersion == 2`.

- [ ] **Step 1: Declare the accessor in `facets.h`.** In `class EmbeddedFacet`, `struct FaceView`, after the `centroid()` line add:

```cpp
        Eigen::VectorXd  area()     const;   // [nF]
```

- [ ] **Step 2: Implement it in `facets.cpp`** next to the other `EmbeddedFacet::FaceView` accessors (after the `centroid()` definition, ~line 409):

```cpp
Eigen::VectorXd  EmbeddedFacet::FaceView::area()       const { return m.lightGeometry().faceAreas; }
```

(`MeshGeometry::faceAreas` is a `VectorXd [nF]` already filled by `lightGeometry()`; confirmed in `include/nxr/compute.h`.)

- [ ] **Step 3: Marshal it + bump schema in the MEX `buildEmbeddedStruct`.** Replace the face-struct block:

```cpp
    { const char* f[] = {"normal","grid","centroid"};
      mxArray* g = mxCreateStructMatrix(1,1,3,f);
      mxSetField(g,0,"normal",   eigenMatrixToMx(fv.normal()));
      mxSetField(g,0,"grid",     eigenComplexMatrixToMx(fv.grid()));
      mxSetField(g,0,"centroid", eigenMatrixToMx(fv.centroid()));
      mxSetField(s,0,"face",g); }
```

with (add `"area"`, bump count 3→4, add the `mxSetField`):

```cpp
    { const char* f[] = {"normal","grid","centroid","area"};
      mxArray* g = mxCreateStructMatrix(1,1,4,f);
      mxSetField(g,0,"normal",   eigenMatrixToMx(fv.normal()));
      mxSetField(g,0,"grid",     eigenComplexMatrixToMx(fv.grid()));
      mxSetField(g,0,"centroid", eigenMatrixToMx(fv.centroid()));
      mxSetField(g,0,"area",     eigenVectorToMx(fv.area()));
      mxSetField(s,0,"face",g); }
```

And bump the Embedded schema line `mxSetField(s,0,"schemaVersion",scalarToMx(1));` → `scalarToMx(2)`. Update the header comment block to list `E.face.area [nF x 1] double` and `schemaVersion == 2`.

- [ ] **Step 4: Commit the source change.**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
git add include/nxr/facets.h src/facets.cpp bindings/mex/src/nxr_compute_mex.cpp
git commit -m "feat(facets): expose Embedded.face.area [nF] (schemaVersion 2)"
```

### Task A2: Rebuild the MEX, install it, verify `face.area` from MATLAB

**Files:**
- Build artifact: `nxr-compute/build/Release/nxr_compute.mexmaca64`
- Install target: `~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/nxr_compute.mexmaca64`

- [ ] **Step 1: Build.**

```bash
cd /Users/diellorbasha/workspace/research/code/nxr-compute
bash scripts/build.sh Release 2>&1 | tail -20
```

Expected: build completes; `build/Release/nxr_compute.mexmaca64` is updated. If the toolchain fails, STOP and report (do not fall back silently).

- [ ] **Step 2: Back up the installed MEX and install the new one.**

```bash
DST=~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b
cp "$DST/nxr_compute.mexmaca64" "$DST/nxr_compute.mexmaca64.bak-$(date +%Y%m%d-%H%M)"
cp /Users/diellorbasha/workspace/research/code/nxr-compute/build/Release/nxr_compute.mexmaca64 "$DST/nxr_compute.mexmaca64"
```

- [ ] **Step 3: Verify from MATLAB (MCP).** Restart any stale MEX by clearing the function, then create a small mesh context and check the field. Run via MATLAB MCP `evaluate_matlab_code`:

```matlab
clear nxr_compute            % drop any loaded stale MEX (function only, NOT workspace)
rehash path
[V,F] = bst_canonical_cortex(2562);   % small canonical mesh (never hand-build)
h = nxr_compute('create', V, double(F));
S = nxr_compute('facets', h, 'trivial', struct());
nxr_compute('destroy', h);
assert(isfield(S.Embedded.face,'area'), 'no area field');
a = S.Embedded.face.area;
fprintf('schemaVersion=%d  nF=%d  min=%.3e sum=%.4f\n', S.Embedded.schemaVersion, numel(a), min(a), sum(a));
assert(S.Embedded.schemaVersion==2 && numel(a)==size(F,1) && all(a>0), 'area/schema check failed');
disp('PASS: Embedded.face.area present, positive, schemaVersion 2');
```

Expected: `PASS` line; `min(a) > 0`; `numel(a) == nFaces`.

- [ ] **Step 4: Commit nothing (build artifacts are not tracked in brainstorm3).** Record the install in the task log only.

### Task A3: `tess_manifold` recompute-on-stale-schema gate

**Files:**
- Modify: `brainstorm3/toolbox/anatomy/tess_manifold.m` (reuse path, around the `bst_get('ManifoldFileForSurface', …)` block ~line 86-100)
- Test: `brainstorm3/dev/tests/test_manifold_face_area.m` (create)

**Interfaces:**
- Consumes: `Embedded(1).schemaVersion` from the loaded manifold node.
- Produces: a constant `REQUIRED_EMBEDDED_SCHEMA = 2` gate; stale nodes recompute.

- [ ] **Step 1: Write the failing test** `dev/tests/test_manifold_face_area.m`:

```matlab
function test_manifold_face_area()
% tess_manifold yields Embedded.face.area (post nxr schemaVersion 2), positive, summing to area.
    SurfaceFile = bst_get('Subject',1).Surface(bst_get('Subject',1).iCortex).FileName;
    M = tess_manifold(SurfaceFile, 'ForceRecompute', true);   % force a fresh nxr facets
    for hh = 1:numel(M.Embedded)
        a = M.Embedded(hh).face.area;
        nF = numel(M.Embedded(hh).GlobalFaces);
        assert(M.Embedded(hh).schemaVersion>=2, 'h%d schemaVersion<2', hh);
        assert(numel(a)==nF && all(a>0), 'h%d area size/positivity', hh);
        fprintf('  h%d nF=%d sum(area)=%.4f m^2\n', hh, nF, sum(a));
    end
    disp('PASS: test_manifold_face_area');
end
```

- [ ] **Step 2: Run it — expect FAIL** if `tess_manifold` reuses a stale (schema-1) node or lacks `ForceRecompute`. Via MCP `run_matlab_file`. Expected: error on missing `area` or schema<2 (proves the gate is needed).

- [ ] **Step 3: Add the schema gate to `tess_manifold.m`.** Near the option parsing add the constant; in the reuse branch (where an existing manifold of the requested Gauge is loaded and returned), add a schema check that forces recompute. Concretely, after loading `existEntry`/the cached `ManifoldMat` in the non-ForceRecompute reuse path, before returning it:

```matlab
    REQUIRED_EMBEDDED_SCHEMA = 2;   % face.area added in nxr facets schemaVersion 2
    ...
    % inside the reuse branch, after loading the candidate node `Mcand`:
    emb = Mcand.Embedded;
    isStale = isempty(emb) || ~isfield(emb,'schemaVersion') || isempty(emb(1).schemaVersion) ...
              || (emb(1).schemaVersion < REQUIRED_EMBEDDED_SCHEMA) || ~isfield(emb(1).face,'area');
    if ~isStale
        ManifoldMat = Mcand;   % reuse
        return;
    end
    % else: fall through to recompute (optionally db_delete the stale node first, as the overwrite path does)
```

(Match the file's actual reuse-branch variable names; the logic is: stale schema OR missing `face.area` ⇒ do not reuse ⇒ recompute.)

- [ ] **Step 4: Run the test — expect PASS.** Via MCP `run_matlab_file dev/tests/test_manifold_face_area.m`. Expected: two `h# nF=… sum(area)=…` lines + `PASS`.

- [ ] **Step 5: Commit.**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/anatomy/tess_manifold.m dev/tests/test_manifold_face_area.m
git commit -m "feat(manifold): face.area + recompute-on-stale-schema gate (REQUIRED_EMBEDDED_SCHEMA=2)"
```

### Task A4 (scoped cleanup): `bst_face_leadfield` reads `Embedded.face.area`

**Files:** Modify `brainstorm3/toolbox/forward/bst_face_leadfield.m` (~line 193-201)

- [ ] **Step 1: Replace local area recompute with a manifold read.** Where `A_f` is computed (the comment "tess_manifold carries no face-area field"), instead read it from the manifold loop already present:

```matlab
        A_f(gf,:)   = M.Embedded(hh).face.area;        % canonical face area (schemaVersion>=2)
```

Remove the local `face_areas`/cross-product area computation and update the stale comment.

- [ ] **Step 2: Run the face-leadfield test** (`dev/tests/` whichever exercises `bst_face_leadfield`, e.g. `test_face_leadfield.m`) via MCP. Expected: unchanged pass (areas identical to within fp tolerance).

- [ ] **Step 3: Commit.**

```bash
git add toolbox/forward/bst_face_leadfield.m
git commit -m "refactor(leadfield): read face area from manifold node (single source of truth)"
```

---

## Phase B — `bst_dirac_helmholtz` → `bst_helmholtz`

### Task B1: Rename + Domain dispatch + vertex geometry from manifold

**Files:**
- Rename: `toolbox/math/bst_dirac_helmholtz.m` → `toolbox/math/bst_helmholtz.m`
- Modify: the renamed file (dispatch + vertex `Prepare`/`Frame` geometry source)
- Modify test: `dev/tests/test_dirac_helmholtz.m` (re-point + new signature)

**Interfaces:**
- Produces:
  - `bst_helmholtz('Prepare', OperatorNode, ManifoldMat, 'Domain','vertex') -> Op` (sets `Op.Domain='vertex'`). `OperatorNode` carries `FirstOrder.Intrinsic{hh}` (Dirac D), `Operator{hh}`/`Mass{hh}` (LBO K/M), `GlobalVertices{hh}`. `ManifoldMat` is `db_template('manifoldmat')` with `Embedded(hh).{face.normal,face.area,face.centroid,vertex.position,vertex.normal,GlobalVertices,GlobalFaces}`.
  - `bst_helmholtz('Frame', Op, Jt [,withCores]) -> Ht` (routes on `Op.Domain`).
  - Shared helper `i_pinned_solve(chol, rhs, free)` (defined in B3).

- [ ] **Step 1: Rename the file (preserve history).**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git mv toolbox/math/bst_dirac_helmholtz.m toolbox/math/bst_helmholtz.m
```

- [ ] **Step 2: Rename the function + add Domain dispatch.** In `bst_helmholtz.m`, rename `function varargout = bst_dirac_helmholtz(...)` → `bst_helmholtz`. Change the dispatcher so `Prepare` records the domain and `Frame` routes:

```matlab
function varargout = bst_helmholtz(varargin)
    [varargout{1:nargout}] = feval(varargin{:});
end

function Op = Prepare(OperatorNode, ManifoldMat, varargin) %#ok<DEFNU>
    Domain = 'vertex';
    for i = 1:2:numel(varargin)
        if strcmpi(varargin{i},'Domain'), Domain = lower(varargin{i+1}); end
    end
    switch Domain
        case 'vertex', Op = i_prepare_vertex(OperatorNode, ManifoldMat);
        case 'face',   Op = i_prepare_face(OperatorNode, ManifoldMat);
        otherwise, error('bst_helmholtz:badDomain','Domain must be vertex|face');
    end
    Op.Domain = Domain;
end

function Ht = Frame(Op, Jt, withCores) %#ok<DEFNU>
    if nargin<3 || isempty(withCores), withCores = true; end
    switch Op.Domain
        case 'vertex', Ht = i_frame_vertex(Op, Jt, withCores);
        case 'face',   Ht = i_frame_face(Op, Jt, withCores);
    end
end
```

- [ ] **Step 3: Move the current vertex `Prepare` body into `i_prepare_vertex`, sourcing geometry from `ManifoldMat`.** Replace the per-hemisphere geometry recompute. Old (delete):

```matlab
        % face normals (oriented outward via the surface vertex normals)
        e1 = Vloc(Floc(:,2),:) - Vloc(Floc(:,1),:);
        e2 = Vloc(Floc(:,3),:) - Vloc(Floc(:,1),:);
        Nf = cross(e1, e2, 2);  Af = sqrt(sum(Nf.^2,2));
        Nf = Nf ./ max(Af, eps);
        vn = Surf.VertNormals(vH(Floc(:,1)), :);
        flip = sum(Nf .* vn, 2) < 0;  Nf(flip,:) = -Nf(flip,:);
```

New (read canonical geometry; keep `Floc`/`Vloc` only for the FEM gradient + incidence assembly, sourced from manifold vertex positions):

```matlab
        E   = ManifoldMat.Embedded(hh);
        vH  = double(E.GlobalVertices(:));
        Vloc = E.vertex.position;                 % [nVh x 3] canonical positions
        Nf   = E.face.normal;                     % [nFh x 3] gauge-consistent outward
        Aface = E.face.area;                      % [nFh x 1] true area (no twoA/Af hazard)
        % local face connectivity from the operator/topology scatter (GlobalFaces -> local)
        Floc = i_local_faces(ManifoldMat, hh);    % [nFh x 3] local vertex indices (helper below)
```

Use `Aface` (true area) for the face→vertex incidence `Wfv` weights (replaces the old `Af`=2·area entries; harmless either way but now unambiguous). Use `Vloc`, `Nf` for the FEM gradient blocks `Gx/Gy/Gz` exactly as before (they consume geometry, now manifold-sourced). Set `Op.VtxH{hh}=Vloc`, `Op.VnH{hh}=E.vertex.normal`, `Op.NbH{hh}` from `Floc` 1-ring.

Add the `i_local_faces` helper (reconstruct local face→vertex indices from the manifold scatter maps — the manifold stores `GlobalFaces` and `GlobalVertices`; map global faces' vertex triples into local indices):

```matlab
function Floc = i_local_faces(ManifoldMat, hh)
    E  = ManifoldMat.Embedded(hh);
    Fg = double(ManifoldMat.Topology(hh).face.vertexIndices); % [nFh x 3] GLOBAL vertex ids (see note)
    vH = double(E.GlobalVertices(:));
    mapV = zeros(max(vH),1); mapV(vH) = 1:numel(vH);
    Floc = mapV(Fg);
end
```

NOTE: confirm the manifold exposes per-face vertex indices (Topology group). If not present, the caller already has `Surf.Faces`; in that case pass `Surf` through `Prepare` ONLY for connectivity (topology, not geometry) and derive `Floc` from `Surf.Faces(GlobalFaces,:)`. Decide during execution based on the actual Topology fields; prefer manifold, fall back to `Surf.Faces` connectivity (still no geometry recompute).

- [ ] **Step 4: Move the vertex `Frame` body into `i_frame_vertex`** unchanged except it consumes `Op` fields now populated from the manifold (no code change beyond the rename — the math is identical). Keep `Decompose` (whole-series) routed to vertex.

- [ ] **Step 5: Re-point + update `dev/tests/test_dirac_helmholtz.m`.** Change calls from `bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf)` to resolve a manifold and use the new signature:

```matlab
    Op = bst_helmholtz('Prepare', Dirac, ManifoldMat, 'Domain','vertex');
    Ht = bst_helmholtz('Frame', Op, J(:,t));
```

where `Dirac` is the operator node (carrying LBO too, or pass a combined operator struct — see B-note) and `ManifoldMat = tess_manifold(SurfaceFile)`. The `PoissonSolve`-verb assertion (line 15) is handled in B3.

- [ ] **Step 6: Run the vertex Helmholtz test** via MCP. Expected: HarmFrac→~0 round-trip and orthonormality assertions PASS, matching pre-refactor values.

- [ ] **Step 7: Commit.**

```bash
git add toolbox/math/bst_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "refactor(helmholtz): rename to bst_helmholtz, Domain dispatch, vertex geometry from manifold"
```

**B-note (operator inputs):** the vertex path needs both Dirac `D` and LBO `K`/`M`. Either pass a single operator node that carries both (preferred — extend the resolver), or pass `{Dirac, LBO}` as a cell in the `OperatorNode` slot. Choose the cell form `OperatorNode = {DiracNode, LBONode}` to avoid changing the operator-node schema; `i_prepare_vertex` unpacks `D=OperatorNode{1}.FirstOrder.Intrinsic{hh}`, `K=OperatorNode{2}.Operator{hh}`, `M=OperatorNode{2}.Mass{hh}`. The face path takes a single face operator node.

### Task B2: Fold in the face branch; gradFace from operator node; delete `_face`

**Files:**
- Modify: `toolbox/math/bst_helmholtz.m` (add `i_prepare_face`, `i_frame_face`, face helpers)
- Delete: `toolbox/math/bst_dirac_helmholtz_face.m`
- Modify test: `dev/tests/test_dirac_helmholtz_face.m`

**Interfaces:**
- Consumes: face operator node with `FaceAux.GradFace [3F×F]`, `GlobalVertices`/`GlobalFaces`; `ManifoldMat.Embedded(hh).face.{normal,area,centroid}`.
- Produces: `i_prepare_face`, `i_frame_face` (same output struct as vertex `Frame`).

- [ ] **Step 1: Port `_face.m`'s `Prepare` into `i_prepare_face`, replacing the fresh-nxr build with operator-node + manifold reads.** Delete the in-Prepare nxr build:

```matlab
        h = nxr_compute('create', Vloc, Floc);
        G = nxr_compute('operators', h, 'gradFace');
        nxr_compute('destroy', h);
        e1 = Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:);  e2 = Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:);
        Nf = cross(e1,e2,2);  twoA = sqrt(sum(Nf.^2,2));  Nf = Nf./max(twoA,eps);  Af = twoA/2;
        vn = Surf.VertNormals(vH(Floc(:,1)),:);  flip = sum(Nf.*vn,2)<0;  Nf(flip,:)=-Nf(flip,:);
        Cf = (Vloc(Floc(:,1),:)+Vloc(Floc(:,2),:)+Vloc(Floc(:,3),:))/3;
```

Replace with:

```matlab
        E   = ManifoldMat.Embedded(hh);
        G   = OperatorNode.FaceAux{hh}.GradFace;   % [3F x F] from operator node
        Nf  = E.face.normal;                       % gauge-consistent outward
        Af  = E.face.area;                         % true area
        Cf  = E.face.centroid;                     % barycentric centroid
        nFh = numel(E.GlobalFaces);
```

The remaining face assembly (`SkewG = i_block_rotation(Nf,nFh)*G`, coupled `A`, `cholA`, `freeA`, `NbF` dual adjacency) is unchanged. Keep `i_block_rotation`, `i_dual_adjacency` as shared/face helpers.

- [ ] **Step 2: Port `_face.m`'s `Frame` into `i_frame_face`** unchanged (math identical; consumes the now-manifold-sourced `Op` fields). Port `i_find_cores` (face) as the face core helper.

- [ ] **Step 3: Delete the old face file.**

```bash
git rm toolbox/math/bst_dirac_helmholtz_face.m
```

- [ ] **Step 4: Re-point + update the face test.** In `dev/tests/test_dirac_helmholtz_face.m`, build/resolve a face operator node (Hodge-Face) + manifold, and call:

```matlab
    FaceOp = tess_operators(SurfaceFile, 'Hodge-Face', 'NoSave', true);   % carries FaceAux.GradFace
    Mani   = tess_manifold(SurfaceFile);
    Op = bst_helmholtz('Prepare', FaceOp, Mani, 'Domain','face');
    Ht = bst_helmholtz('Frame', Op, Jf);
```

(Replace the prior `bst_dirac_helmholtz_face('Prepare', DiracOp, LBO, Surf)` call.) Keep the existing B-orthonormality / HarmFrac assertions.

- [ ] **Step 5: Run the face test** via MCP. Expected: HarmFrac→~0, W_F-orthonormality PASS.

- [ ] **Step 6: Commit.**

```bash
git add toolbox/math/bst_helmholtz.m dev/tests/test_dirac_helmholtz_face.m
git commit -m "refactor(helmholtz): fold face branch into bst_helmholtz (gradFace from operator node, geometry from manifold); delete _face"
```

### Task B3: Shared `i_pinned_solve`; delete dead `PoissonSolve` verb

**Files:** Modify `toolbox/math/bst_helmholtz.m`; modify `dev/tests/test_dirac_helmholtz.m`

- [ ] **Step 1: Add the shared solver helper.**

```matlab
function x = i_pinned_solve(dChol, rhs, free)
% Solve a pinned SPD system using a cached Cholesky factor: x(free) = dChol\rhs(free),
% remaining (pinned) entries 0, then recenter to mean-zero. Used by both domains.
    x = zeros(size(rhs));
    x(free) = dChol \ rhs(free);
    x = x - mean(x);
end
```

Refactor `i_poisson` (vertex) to call `i_pinned_solve` after the mean-zero RHS projection; refactor the face coupled solve in `i_frame_face` to call `i_pinned_solve(Op.cholA{hh}, [divS;curlS], Op.freeA{hh})` then split `phi`/`psi`.

- [ ] **Step 2: Delete the `PoissonSolve` verb** (the `function psi = PoissonSolve(...)` block) and its USAGE-comment line.

- [ ] **Step 3: Update the test that used `PoissonSolve`** (`test_dirac_helmholtz.m:15`). Replace the `bst_dirac_helmholtz('PoissonSolve', K, M, ...)` round-trip assertion with an equivalent check driven through `bst_helmholtz('Frame', …)` (verify `Curl`/`Psi` reconstruct a known vorticity field), or assert via a small direct `i_pinned_solve`-equivalent computed inline in the test (since the helper is now internal). Concrete: synthesize `J = n × grad(psi0)` for a known `psi0`, decompose, assert `corr(Ht.Psi, psi0) > 0.99` on the dominant hemisphere.

- [ ] **Step 4: Run the test** via MCP. Expected: PASS.

- [ ] **Step 5: Commit.**

```bash
git add toolbox/math/bst_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "refactor(helmholtz): unify cached-factor solve (i_pinned_solve); drop dead PoissonSolve verb"
```

### Task B4: Update callers + extract shared operator loader

**Files:**
- Modify: `toolbox/gui/view_helmholtz.m`, `toolbox/process/functions/process_vortex_track.m`
- Create: `toolbox/anatomy/bst_get_operator_node.m` (shared resolver) — OR a local in a shared util; prefer a small dedicated function next to `tess_operators`.

- [ ] **Step 1: Extract the duplicated `i_op` into a shared function** `toolbox/anatomy/bst_get_operator_node.m`:

```matlab
function Op = bst_get_operator_node(SurfaceFile, variant)
% Find-or-create an operator node of `variant` under `SurfaceFile`, return the loaded node.
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = in_bst_operator(sSubject.Surface(iSurf).Operator(k).FileName);
            if strcmpi(S.Variant, variant), Op = S; break; end
        end
    end
    if isempty(Op), tess_operators(SurfaceFile, variant); Op = bst_get_operator_node(SurfaceFile, variant); end
end
```

- [ ] **Step 2: Update `view_helmholtz.m`.** Replace the local `i_op` calls + the `Surf` arg with the new resolver + manifold:

```matlab
    Dirac = bst_get_operator_node(SurfaceFile, 'Dirac');
    LBO   = bst_get_operator_node(SurfaceFile, 'Laplace-Beltrami');
    Mani  = tess_manifold(SurfaceFile);
    Op = bst_helmholtz('Prepare', {Dirac, LBO}, Mani, 'Domain','vertex');
    ...
    Ht = bst_helmholtz('Frame', St.Op, Jt, needCores);
```

Delete the local `i_op` function from `view_helmholtz.m`.

- [ ] **Step 3: Update `process_vortex_track.m`** identically (replace `i_op` + `bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf)` with the `{Dirac,LBO}` + `Mani` form; delete its local `i_op`).

- [ ] **Step 4: Run the view + track tests/benchmarks** via MCP: `dev/tests/test_helmholtz_view.m`. Expected: PASS (decomposition + markers unchanged).

- [ ] **Step 5: Commit.**

```bash
git add toolbox/anatomy/bst_get_operator_node.m toolbox/gui/view_helmholtz.m toolbox/process/functions/process_vortex_track.m
git commit -m "refactor(helmholtz): callers pass manifold + shared bst_get_operator_node loader"
```

### Task B5: Parity test + re-point remaining tests/benchmarks

**Files:**
- Create: `dev/tests/test_helmholtz_manifold_parity.m`
- Modify: `dev/benchmarks/bench_dirac_face_helmholtz.m`, `dev/benchmarks/bench_dirac_streamribbon_real.m` (re-point names)

- [ ] **Step 1: Write the parity test** comparing manifold-sourced geometry against the legacy cross/flip recompute, asserting the decomposition is unchanged:

```matlab
function test_helmholtz_manifold_parity()
% Manifold-sourced geometry must reproduce the legacy cross/flip decomposition.
    SurfaceFile = bst_get('Subject',1).Surface(bst_get('Subject',1).iCortex).FileName;
    Surf = in_tess_bst(SurfaceFile, 0);
    Mani = tess_manifold(SurfaceFile);
    % legacy normals (cross + VertNormals flip) per hemisphere vs manifold normals
    for hh = 1:numel(Mani.Embedded)
        E = Mani.Embedded(hh); vH = double(E.GlobalVertices(:));
        Floc = i_legacy_local_faces(Surf, E);     % helper: Surf.Faces(GlobalFaces) -> local
        Vloc = Surf.Vertices(vH,:);
        e1=Vloc(Floc(:,2),:)-Vloc(Floc(:,1),:); e2=Vloc(Floc(:,3),:)-Vloc(Floc(:,1),:);
        Nfl = cross(e1,e2,2); Nfl = Nfl./max(sqrt(sum(Nfl.^2,2)),eps);
        vn = Surf.VertNormals(vH(Floc(:,1)),:); fl = sum(Nfl.*vn,2)<0; Nfl(fl,:)=-Nfl(fl,:);
        dotp = sum(Nfl .* E.face.normal, 2);      % should be ~+1 if signs agree
        fprintf('  h%d normal agreement: min dot = %.4f  (frac<0: %.3f)\n', hh, min(dotp), mean(dotp<0));
        assert(min(dotp) > 0.99, 'h%d face-normal sign/dir mismatch vs manifold', hh);
        a_legacy = sqrt(sum(cross(e1,e2,2).^2,2))/2;
        assert(max(abs(a_legacy - E.face.area)) < 1e-9*max(E.face.area), 'h%d area mismatch', hh);
    end
    disp('PASS: test_helmholtz_manifold_parity');
end
```

(If `min dot` is negative for some faces, that is the sign-convention risk flagged in the spec — STOP and reconcile before relying on the decomposition.)

- [ ] **Step 2: Run the parity test** via MCP. Expected: per-hemisphere agreement `min dot > 0.99`, area match, PASS. If it fails on sign, report and reconcile (do not proceed).

- [ ] **Step 3: Re-point benchmark filenames** (`bst_dirac_helmholtz` → `bst_helmholtz`, signatures updated to `{Dirac,LBO}`/`FaceOp` + `Mani` + `'Domain'`). These are dev scripts; update call sites only.

- [ ] **Step 4: Grep to confirm no stale references remain.**

```bash
grep -rn "bst_dirac_helmholtz" toolbox dev | grep -v "\.bak"
```

Expected: no hits in `toolbox/`; only historical doc/spec mentions remain.

- [ ] **Step 5: Commit.**

```bash
git add dev/tests/test_helmholtz_manifold_parity.m dev/benchmarks/bench_dirac_face_helmholtz.m dev/benchmarks/bench_dirac_streamribbon_real.m
git commit -m "test(helmholtz): manifold-geometry parity guard; re-point benchmarks to bst_helmholtz"
```

---

## Self-review notes

- **Spec coverage:** Phase A (face.area nxr + schema gate + leadfield cleanup) = A1–A4. Phase B (rename, dispatch, vertex+face geometry rewire, i_pinned_solve, drop PoissonSolve, caller+i_op, parity test, re-point) = B1–B5. All spec sections mapped.
- **Sign-convention risk** (spec "Top risk") is covered by B5 parity test as a hard gate.
- **I/O-free constraint:** `bst_helmholtz` takes `OperatorNode` + `ManifoldMat` structs; all `tess_*`/`bst_get`/`in_*` calls live in callers/tests (B1 B-note, B4).
- **Open verification during execution:** whether the manifold `Topology` group exposes per-face vertex indices (B1 Step 3 NOTE) — fall back to `Surf.Faces` connectivity (not geometry) if absent.
