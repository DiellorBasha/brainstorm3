# Derived-Anatomy Nodes — SP2 (`operator_` node) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add the `operator_` derived-anatomy node (surface child holding `{Operator A, Mass B}` per hemisphere) + the `tess_operators` assembler for the three operators (Laplace-Beltrami, Connection Laplacian, Dirac), closing the loop for SP3's `tess_eigen`.

**Architecture:** Replicate SP1's merged `manifold_` node pattern for `operator`. SP1 already provides the surface `.Operator` child list (backfilled by the `db_update` 5.04 migration), `node_create_subject` nesting, `tree_callbacks` menu plumbing, the `'manifold'`/`'operator'` node types render via the unmodified Java `BstNode`, and `db_add_manifold`/`bst_get('ManifoldFile')` as direct templates.

**Tech Stack:** MATLAB (Brainstorm DB/tree), nxr-compute MEX (`operators` command), MATLAB MCP.

**Templates to mirror (all merged on `development`):**
- `toolbox/db/db_add_manifold.m` → `db_add_operator.m`
- `toolbox/core/bst_get.m` `case 'ManifoldFile'` → `case 'OperatorFile'`
- `toolbox/anatomy/tess_manifold.m` (assembler shape: guards, `tess_hemisplit`, per-hemi `nxr_compute('create'/…/'destroy')`, scatter maps, `db_add_*`) → `tess_operators.m`
- `toolbox/db/db_template.m` `manifoldmat` → `operatormat`; `surface` already has `.Operator`.
- `toolbox/io/file_gettype.m`/`file_short.m`/`file_fullpath.m` `manifold` branch → `operator`.
- `toolbox/tree/node_create_subject.m` `.Manifold` nesting loop → `.Operator`.
- `toolbox/tree/tree_callbacks.m` `case 'manifold'` + "Compute manifold" + `ManifoldView_/ManifoldDelete_Callback` → `operator`.

**Session discipline:** run via MATLAB MCP; only ever `clear <named>`, never bare `clear`; DB-mutating verifications restore the protocol (`file_delete` + `db_reload_database('current')`).

**⚑ CANONICAL MESH RULE (non-negotiable):** Never hand-build a mesh (octahedron, sphere, tetrahedron, …) as input to `nxr_compute('create', …)`. geometry-central's `ManifoldSurfaceMesh` **segfaults and crashes MATLAB** on malformed input. Every test/probe that reaches `create` MUST get its mesh from `bst_canonical_cortex(20484)` (the real cortex), split per hemisphere with `tess_hemisplit`, and reindex locally — exactly as `tess_manifold` does. Production/test code must call `nxr_safe_create(V, F)` (validates → clean error) instead of `nxr_compute('create', V, F)`. Synthetic data is allowed ONLY for pure-algorithm tests that never call `create`.

---

## Task 1: DB data types + file-type wiring + `AnyFile` fix

**Files:** `db_template.m`, `file_gettype.m`, `file_short.m`, `file_fullpath.m`, `bst_get.m`.

- [ ] **Step 1: `db_template('operatormat')`** — add (mirror `manifoldmat`):
```matlab
    case 'operatormat'
        template = struct(...
              'Comment',        '', ...
              'ParentSurface',  '', ...
              'Variant',        '', ...   % 'Laplace-Beltrami' | 'Connection Laplacian' | 'Dirac'
              'Operator',       [], ...   % 1x2 per-hemisphere sparse A
              'Mass',           [], ...   % 1x2 per-hemisphere sparse B
              'GlobalVertices', [], ...   % 1x2 per-hemisphere scatter maps
              'Provenance',     []);
```

- [ ] **Step 2: `operator` file type** — in `file_gettype.m` add a `_operator_` branch (mirror `_manifold_`, before `_tess`) returning `'operator'`; in `file_short.m` and `file_fullpath.m` add `'operator'` to the anatomy-type case lists (mirror the `'manifold'` additions).

- [ ] **Step 3: `bst_get('AnyFile')` — recognize manifold + operator** (fixes SP1's "View file contents" error). Read the `case 'AnyFile'` block in `bst_get.m`; add routing so a `manifold`/`operator` file type resolves via the corresponding `ManifoldFile`/`OperatorFile` logic (or returns the subject/surface context `view_struct` needs). Mirror how `AnyFile` dispatches other anatomy types.

- [ ] **Step 4: Verify (MATLAB MCP, DB-free where possible)**
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); rehash;
o = db_template('operatormat'); assert(all(isfield(o,{'Variant','Operator','Mass','GlobalVertices','ParentSurface'})),'operatormat');
assert(strcmpi(file_gettype('operator_x_250101.mat'),'operator'),'file_gettype operator');
assert(strcmpi(file_gettype('manifold_x.mat'),'manifold') && strcmpi(file_gettype('tess_cortex_low.mat'),'cortex'),'others intact');
disp('T1 OK');
```

- [ ] **Step 5: Commit** — `feat(db): operatormat schema + operator file type + bst_get AnyFile for manifold/operator`.

---

## Task 2: `db_add_operator` + `bst_get('OperatorFile')`

**Files:** create `toolbox/db/db_add_operator.m`; modify `bst_get.m`.

- [ ] **Step 1:** Create `db_add_operator(iSubject, ParentSurfaceFile, OperatorMat, Comment)` by copying `db_add_manifold.m` verbatim and substituting `Manifold`→`Operator`, `manifold_`→`operator_`. Keep the surface-array normalization loop (the fresh-rebuild that avoids the dissimilar-structures error) and the `iSubject==0` round-trip. License header included.

- [ ] **Step 2:** Add `case 'OperatorFile'` to `bst_get.m` by mirroring `case 'ManifoldFile'` (scan `Surface(*).Operator`, return `[sSubject,iSubject,iSurface,iOperator]`).

- [ ] **Step 3: Live verify + cleanup** (mirror SP1 Task 3's verify): on a subject with a cortex, `db_add_operator(iSubj, surfFile, struct-with-stub-Operator, 'TEST op')`, confirm `.Operator` list gains the entry, `bst_get('OperatorFile', newFile)` resolves, then delete file + remove entry + `db_reload_database('current')`. Expect `T2 OK`.

- [ ] **Step 4: Commit** — `feat(db): db_add_operator + bst_get OperatorFile`.

---

## Task 3: `tess_operators` assembler (+ test)

**Files:** create `toolbox/anatomy/tess_operators.m`; test `dev/tests/test_tess_operators.m`.

- [ ] **Step 0 (verify nxr calls):** Before writing, confirm the exact nxr `operators` calls **on the canonical cortex** (MCP) — NEVER a hand-built mesh (see the Canonical Mesh Rule above). Get the mesh with `SurfaceFile = bst_canonical_cortex(20484)`, `tess_hemisplit` it, build a local per-hemisphere submesh (mask + reindex, as `tess_manifold` does), then `h = nxr_safe_create(Vloc, Floc)`. On that handle confirm `nxr_compute('operators',h,'laplacian','cotan')`, `('mass','galerkin')`, `('dirac',0.5)` all return as expected; and determine the **Connection Laplacian** call — `nxr_compute('operators',h,'laplacian','connection')` builds in the ACTIVE gauge; if the trivial gauge with FS-pole singularities is required, find how to set it (e.g. an opts struct with `singVerts/singValues`, or a `setGauge`/`gauge` call first). Document the exact connection call. If the connection operator can't be fetched in the trivial gauge via `operators`, report and we adjust (LBO + Dirac can land first).

- [ ] **Step 1: Write `tess_operators(SurfaceFile, OperatorName, varargin)`** — mirror `tess_manifold` structure:
  - Options: `Tau` (0.5, for Dirac/connection trivial gauge), `NoSave`, `ForceRecompute`.
  - Map `OperatorName` → `Variant`: `'Laplace-Beltrami'`|`'Connection Laplacian'`|`'Dirac'` (accept case-insensitive; error on unknown with a clear list).
  - Guards (nxr install; Structures L/R atlas). **NOTE (Step-0 finding):** unlike `tess_manifold`'s `facets` call, the three *operators* do NOT need the trivial gauge or FS-pole singularities — `laplacian/connection` builds the intrinsic Levi-Civita connection internally (verified bit-identical with/without a `facets` call) and `dirac`/`cotan`/`mass` act directly on the discrete submesh. So **no `Reg.Sphere` guard and no `singVerts/singValues` opts** here. `tess_hemisplit`; per hemisphere build local submesh, `h=nxr_safe_create(Vloc,Floc)` (validating wrapper — NOT raw `nxr_compute('create')`):
    - **Laplace-Beltrami:** `A = nxr_compute('operators',h,'laplacian','cotan')`; `B = nxr_compute('operators',h,'mass','galerkin')`.
    - **Connection Laplacian:** `A = nxr_compute('operators',h,'laplacian','connection')` (complex Hermitian; intrinsic Levi-Civita connection — no gauge/FS-pole config needed, per Step-0 finding); `B = nxr_compute('operators',h,'mass','galerkin')`.
    - **Dirac:** `A = nxr_compute('operators',h,'dirac',Tau)`; `B = kron(nxr_compute('operators',h,'mass','galerkin'), speye(4))`.
    - `nxr_compute('destroy',h)`.
  - Assemble `OperatorMat = db_template('operatormat')`: 1×2 `Operator`/`Mass`/`GlobalVertices`, `Variant`, `Provenance` (Backend='nxr', NxrVersion, Variant, Tau if applicable, ComputeDate).
  - `iSubject` via `bst_get('SurfaceFile', SurfaceFile)`; `db_add_operator(...)` unless `NoSave`. Comment e.g. `sprintf('%s operator', Variant)`. Return `OperatorMat`. License header.

- [ ] **Step 2: Test `dev/tests/test_tess_operators.m`** (real 20484 cortex via `bst_canonical_cortex(20484)` — NEVER a hand-built mesh; backup/restore via `onCleanup`):
  - For each of the 3 variants: `tess_operators(cortexFile, variant, 'NoSave',1, 'ForceRecompute',1)` returns 1×2 `Operator/Mass`.
  - Shapes/symmetry: LBO `A` real `[nVₕ×nVₕ]` symmetric, `B` `[nVₕ×nVₕ]`; Dirac `A` `[4nVₕ×4nVₕ]` symmetric, `B == kron(Mg,I4)` `[4nVₕ×4nVₕ]`; Connection `A` `[nVₕ×nVₕ]` (complex), Hermitian.
  - Smoke eigensolve: `d = eigs(A, B, 6, 'smallestabs')` returns finite real(ish) eigenvalues per variant (Dirac → 4-fold clusters).
  - Save path for one variant: `tess_operators(cortexFile,'Dirac','ForceRecompute',1)` creates an `operator_` node (`.Operator` non-empty; `bst_get('OperatorFile',…)` resolves); cleanup.

- [ ] **Step 3: Run / commit** — `feat(tess-operators): per-hemisphere {operator,mass} nodes for LBO/connection/Dirac`.

---

## Task 4: Tree nesting + context menus

**Files:** `node_create_subject.m`, `tree_callbacks.m`.

- [ ] **Step 1:** In `node_create_subject.m`, add an `.Operator` child-nesting loop beside the `.Manifold` loop (same `CreateNode('operator', …)` + `isfield` guard, `nodeSurface.add`).

- [ ] **Step 2:** In `tree_callbacks.m`:
  - Surface popup: a **"Compute operator"** submenu with three items → `bst_call(@tess_operators, filenameRelative, 'Laplace-Beltrami'|'Connection Laplacian'|'Dirac')`. (Use a `gui_component('Menu', …)` parent + three MenuItems; place near "Compute manifold", read-only guarded.)
  - `case 'operator'` node popup: **View** (`OperatorView_Callback`, field-name inspector mirroring `ManifoldView_Callback`) + **Delete** (`OperatorDelete_Callback`, mirror `ManifoldDelete_Callback`, single-selection guarded).
  - Add `'operator'` to the `filenameFull` SUBJECTS switch.

- [ ] **Step 3: Verify** (MCP): create an operator node via `tess_operators`, confirm `.Operator` non-empty + `CreateNode('operator',…)` returns created; exercise the Delete callback path; cleanup. GUI tree rendering = user's visual check.

- [ ] **Step 4: Commit** — `feat(tree): nest operator nodes under surface + Compute operator / View / Delete menus`.

---

## Final verification
- [ ] All three operator variants compute into nodes on the real cortex; nodes nest under the surface; View/Delete work; `bst_get('AnyFile')` resolves manifold + operator ("View file contents" works). Restore protocol.
- [ ] Then **superpowers:finishing-a-development-branch**; proceed to SP3 (`tess_eigen` reading operator nodes → eigen nodes).

## Notes for the implementer
- This is a faithful replication of merged SP1 — read each SP1 template and match its idioms; the only genuinely new logic is `operatormat`, the per-variant nxr fetch in `tess_operators`, and the `AnyFile` fix.
- **Connection-Laplacian nxr call is the one unknown** (Step 0) — verify it before implementing that branch; LBO + Dirac can land even if connection needs follow-up.
- No Java build; `'operator'` renders with the default icon (custom icon = optional later [[bst-java-fork]] work).
- DB-mutating tests must restore the protocol (`file_delete` + reload).
