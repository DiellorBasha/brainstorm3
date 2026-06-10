# Derived-Anatomy Nodes — SP3 (`eigen_` node) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Add the `eigen_` derived-anatomy node (one file/node per eigensolve, holding an operator's eigenbasis) + the `tess_eigen` assembler (find-or-create the operator, `eigs(A,B)`, Rayleigh–Ritz for Dirac), closing the derived-anatomy trio (manifold/operator/eigen).

**Architecture:** Replicate SP2's merged `operator_` node pattern for `eigen`. SP1+SP2 already provide the surface `.Eigen` child list (backfilled by the `db_update` 5.04 migration), `node_create_subject` nesting, `tree_callbacks` menu plumbing, the `'eigen'` node type rendering via the unmodified Java `BstNode`, and `db_add_operator`/`bst_get('OperatorFile')`/`tess_operators` as direct templates. The only genuinely new logic is the `eigenmat` schema, the find-or-create-operator + eigensolve in `tess_eigen`, and the `local_ritz_basis` nested function.

**Tech Stack:** MATLAB (Brainstorm DB/tree), nxr-compute MEX (via `tess_operators`), MATLAB MCP.

**Templates to mirror (all merged on `development`):**
- `toolbox/db/db_add_operator.m` → `db_add_eigen.m`
- `toolbox/core/bst_get.m` `case 'OperatorFile'` → `case 'EigenFile'`; `case 'operator'` in `AnyFile` → `case 'eigen'`
- `toolbox/anatomy/tess_operators.m` (assembler shape: option parsing, guards, return, license) → `tess_eigen.m`
- `toolbox/anatomy/tess_dirac_eigenmodes.m` (the `eigs(A,B,...)` + `local_ritz_basis` Rayleigh–Ritz logic, and the Dirac multiplet over-fetch cap) → copy `local_ritz_basis` as a nested function in `tess_eigen.m`
- `toolbox/db/db_template.m` `operatormat` → `eigenmat`; `surface` already has `.Eigen`.
- `toolbox/io/file_gettype.m`/`file_short.m`/`file_fullpath.m` `operator` branch → `eigen`.
- `toolbox/tree/node_create_subject.m` `.Operator` nesting loop → `.Eigen`.
- `toolbox/tree/tree_callbacks.m` `case 'operator'` + "Compute operator" + `OperatorView_/OperatorDelete_Callback` → `eigen`.

**⚑ CANONICAL MESH RULE (non-negotiable):** `tess_eigen` reaches `nxr_compute('create')` transitively via `tess_operators` (find-or-create). All tests/probes get their mesh ONLY from `bst_canonical_cortex(20484)`; never hand-build a mesh. `tess_operators` already routes through `nxr_safe_create`.

**Session discipline:** run via MATLAB MCP; start `if ~brainstorm('status'); brainstorm nogui; end` + `[isOk,~]=bst_plugin('Install','nxr-compute')`; only ever `clear <named>`, never bare `clear`; `rehash` after editing .m files; DB-mutating verifications restore the protocol (`file_delete` + `db_reload_database('current')`).

---

## Task 1: DB data type + file-type wiring + `AnyFile`

**Files:** `db_template.m`, `file_gettype.m`, `file_short.m`, `file_fullpath.m`, `bst_get.m`.

- [ ] **Step 1: `db_template('eigenmat')`** — add (mirror `operatormat`):
```matlab
    case 'eigenmat'
        template = struct(...
              'Comment',        '', ...
              'ParentSurface',  '', ...
              'OperatorFile',   '', ...   % path ref to the operator node solved
              'Variant',        '', ...   % 'Laplace-Beltrami' | 'Connection Laplacian' | 'Dirac'
              'Phi',            [], ...   % 1x2 per-hemi cell: eigenvectors
              'Lambda',         [], ...   % 1x2 per-hemi cell: eigenvalues [K x 1]
              'K',              [], ...   % modes requested
              'GlobalVertices', [], ...   % 1x2 per-hemi cell: scatter maps
              'Provenance',     []);
```

- [ ] **Step 2: `eigen` file type** — in `file_gettype.m` add an `_eigen_` branch (mirror `_operator_`, before `_tess`) returning `'eigen'`; in `file_short.m` and `file_fullpath.m` add `'eigen'` to the anatomy-type case lists (mirror the `'operator'` additions).

- [ ] **Step 3: `bst_get('AnyFile')` — recognize eigen** — add a `case 'eigen'` to the `AnyFile` switch in `bst_get.m`, mirroring `case 'operator'`: resolve via `EigenFile`, and when `nargout>=5` set `sItem = sStudy.Surface(iSurf).Eigen(iItem)`.

- [ ] **Step 4: Verify (MATLAB MCP)**
```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3'); rehash;
e = db_template('eigenmat'); assert(all(isfield(e,{'OperatorFile','Variant','Phi','Lambda','K','GlobalVertices','ParentSurface'})),'eigenmat');
assert(strcmpi(file_gettype('eigen_x_250101.mat'),'eigen'),'file_gettype eigen');
assert(strcmpi(file_gettype('operator_x.mat'),'operator') && strcmpi(file_gettype('tess_cortex_low.mat'),'cortex'),'others intact');
[a,b,c,d,ee] = bst_get('AnyFile','eigen_x_250101.mat'); %#ok<ASGLU>  % must NOT error (empty match OK)
disp('T1 OK');
```

- [ ] **Step 5: Commit** — `feat(db): eigenmat schema + eigen file type + bst_get AnyFile for eigen`.

---

## Task 2: `db_add_eigen` + `bst_get('EigenFile')`

**Files:** create `toolbox/db/db_add_eigen.m`; modify `bst_get.m`.

- [ ] **Step 1:** Create `db_add_eigen(iSubject, ParentSurfaceFile, EigenMat, Comment)` by copying `db_add_operator.m` verbatim and substituting `Operator`→`Eigen`, `operator_`→`eigen_`. Keep the surface-array normalization loop (fresh-rebuild avoiding the dissimilar-structures error) and the `iSubject==0` round-trip. License header included.

- [ ] **Step 2:** Add `case 'EigenFile'` to `bst_get.m` by mirroring `case 'OperatorFile'` (scan `Surface(*).Eigen`, return `[sSubject,iSubject,iSurface,iEigen]`). Add the usage line to the header comment block (mirror the `OperatorFile` line).

- [ ] **Step 3: Live verify + cleanup** (mirror SP2 Task 2): on a subject with a cortex, `db_add_eigen(iSubj, surfFile, struct-with-stub-Phi, 'TEST eig')`, confirm `.Eigen` list gains the entry, `bst_get('EigenFile', newFile)` resolves, then delete file + remove entry + `db_reload_database('current')`. Expect `T2 OK`.
```matlab
if ~brainstorm('status'); brainstorm nogui; end
P = bst_get('ProtocolSubjects'); iSubj = 1;
surfFile = ''; for i=1:numel(P.Subject(iSubj).Surface); if strcmpi(P.Subject(iSubj).Surface(i).SurfaceType,'Cortex'); surfFile = P.Subject(iSubj).Surface(i).FileName; break; end; end
em = db_template('eigenmat'); em.Variant='Dirac'; em.Phi={sparse(4,2),sparse(4,2)}; em.ParentSurface=surfFile;
newFile = db_add_eigen(iSubj, surfFile, em, 'TEST eig');
[s,iS,iSurf,iE] = bst_get('EigenFile', newFile); assert(~isempty(s) && ~isempty(iE), 'EigenFile resolves');
file_delete(file_fullpath(newFile), 1);
PS = bst_get('ProtocolSubjects'); PS.Subject(iS).Surface(iSurf).Eigen(iE) = []; bst_set('ProtocolSubjects', PS); db_save();
db_reload_database('current'); disp('T2 OK');
```

- [ ] **Step 4: Commit** — `feat(db): db_add_eigen + bst_get EigenFile`.

---

## Task 3: `tess_eigen` assembler (+ test)

**Files:** create `toolbox/anatomy/tess_eigen.m`; test `dev/tests/test_tess_eigen.m`.

- [ ] **Step 1: Write `tess_eigen(SurfaceFile, OperatorName, varargin)`** — mirror `tess_operators` structure:
  - Options: `K` (default 400), `Tau` (0.5, forwarded to Dirac operator creation), `NoSave`, `ForceRecompute`.
  - Map `OperatorName` → `Variant` (case-insensitive): `'Laplace-Beltrami'`|`'Connection Laplacian'`|`'Dirac'` (error on unknown with the list).
  - Guards: nxr install; resolve `[~,iSubject,iSurface]` via `bst_get('SurfaceFile', SurfaceFile)`.
  - **Find-or-create operator:** read `sSubject.Surface(iSurface).Operator`; find the first entry with `strcmpi(.Variant, Variant)`. If found, `OperatorFile = that.FileName`. Else `OperatorMat = tess_operators(SurfaceFile, OperatorName, 'Tau', Tau)` and re-resolve its `OperatorFile` via `bst_get('OperatorFile', ...)` (or have `tess_operators` return the saved path; resolve from the refreshed `.Operator` list by Variant). Load the operator node with `Op = load(file_fullpath(OperatorFile))` (it is a plain `.mat`, not a standard bst file type) → per-hemi cells `A = Op.Operator{hh}`, `B = Op.Mass{hh}`, `gv = Op.GlobalVertices{hh}`.
  - **Per hemisphere `hh = 1:2`:**
    - `nVh` from `size(A,1)` (LBO/Connection `nVh`; Dirac `4*nVh`).
    - Dirac multiplet handling: round `K` to a multiple of 4 and over-fetch (reuse the `nRequest = min(K + max(8, ceil(0.3*K)), <cap>)` logic from `tess_dirac_eigenmodes`); LBO/Connection: `nRequest = min(K+8, size(A,1)-2)`.
    - `[V, D] = eigs(A, B, nRequest, 'smallestabs');` `lam = real(diag(D));` sort ascending.
    - **B-orthonormalize.** Dirac AND Connection → `V = local_ritz_basis(A, B, V);` (nested fn, see Step 1b) — BOTH have degenerate multiplets, so plain normalization leaves the set B-non-orthonormal (Connection residual ~0.16 observed). LBO → standard B-orthonormalization (`V = V ./ sqrt(real(diag(V' * B * V))).'` with a small rank guard). Keep the first `K` columns; `Phi{hh}=V(:,1:K); Lambda{hh}=lam(1:K);`
    - `GlobalVertices{hh} = gv;`
  - Assemble `EigenMat = db_template('eigenmat')` with the 1×2 cells, `OperatorFile`, `Variant`, `ParentSurface=SurfaceFile`, `K`, and `Provenance` (Backend='nxr', NxrVersion from `nxr_compute('version')`, Variant, K, Tau when Dirac, Ortho='Rayleigh-Ritz' for Dirac else 'B-orthonormal', ComputeDate=`datestr(now,'yyyy-mm-dd HH:MM:SS')`).
  - `db_add_eigen(iSubject, SurfaceFile, EigenMat, sprintf('%s eigenmodes (K=%d)', Variant, K))` unless `NoSave`. Return `EigenMat`. License header.

- [ ] **Step 1b: `local_ritz_basis` nested function** — adapt the Rayleigh–Ritz B-orthonormalization from `tess_dirac_eigenmodes.m` into `tess_eigen.m` as a `local_`-prefixed nested/subfunction (rank-revealing B-orthonormalization + Rayleigh–Ritz; throws `:rankDeficient` if the fetched set is too small). Must be **complex-safe** (use `'` Hermitian transpose, not `.'`) since the Connection variant is complex Hermitian. Do NOT create a standalone shared file; do NOT modify `tess_dirac_eigenmodes.m`.

- [ ] **Step 2: Test `dev/tests/test_tess_eigen.m`** (real 20484 cortex via `bst_canonical_cortex(20484)` — NEVER a hand-built mesh; backup/restore via `onCleanup`):
  - Use a small `K` (e.g. `K=12`, Dirac rounds to 12) to keep eigs fast.
  - For each of the 3 variants: `tess_eigen(cortexFile, variant, 'K',12, 'NoSave',1, 'ForceRecompute',1)` returns 1×2 `Phi`/`Lambda`.
  - Shapes: LBO/Connection `Phi{hh}` `[nVh×K]`, Dirac `Phi{hh}` `[4nVh×K]`; `Lambda{hh}` `[K×1]`.
  - **B-orthonormality:** load the matching operator (`bst_get('OperatorFile')` on `EigenMat.OperatorFile`), per hemi assert `norm(Phi' * B * Phi - I) < 1e-6` (Dirac uses the Rayleigh–Ritz path).
  - `Lambda` ascending and real-ish; Dirac → near-4-fold clusters.
  - **Find-or-create:** on a surface with no operator node of that Variant, `tess_eigen(..., 'NoSave',1)` creates the operator node (assert `.Operator` gains a matching-Variant entry); a second call reuses it (no new operator entry). Clean up created operator/eigen nodes.
  - Save path: `tess_eigen(cortexFile,'Dirac','K',12)` creates an `eigen_` node (`.Eigen` non-empty; `bst_get('EigenFile',…)` resolves; stored `OperatorFile` resolves via `bst_get('OperatorFile')`); cleanup (delete eigen + operator files, remove entries, reload).

- [ ] **Step 3: Run / commit** — static-check; run the test live (iterate to green); restore protocol. Commit `feat(tess-eigen): per-hemisphere eigenbasis nodes (find-or-create operator; Rayleigh-Ritz for Dirac)`.

---

## Task 4: Tree nesting + context menus

**Files:** `node_create_subject.m`, `tree_callbacks.m`.

- [ ] **Step 1:** In `node_create_subject.m`, add an `.Eigen` child-nesting loop beside the `.Operator` loop (same `CreateNode('eigen', …)` + `isfield` guard, `nodeSurface.add`).

- [ ] **Step 2:** In `tree_callbacks.m`:
  - Surface popup: a **"Compute eigenmodes"** submenu with three items → `bst_call(@tess_eigen, filenameRelative, 'Laplace-Beltrami'|'Connection Laplacian'|'Dirac')`, placed near "Compute operator", under the same `cortex && ~ReadOnly` guard.
  - `case 'eigen'` node popup: **View** (`EigenView_Callback`, field-name inspector mirroring `OperatorView_Callback`) + **Delete** (`EigenDelete_Callback`, mirror `OperatorDelete_Callback`, single-selection guarded, file_delete + remove `.Eigen` entry + `bst_set` + `db_save` + `UpdateNode`).
  - Add `'eigen'` to the `filenameFull` SUBJECTS switch.

- [ ] **Step 3: Verify** (MCP): create an eigen node via `tess_eigen` (small `K`), confirm `.Eigen` non-empty + a programmatic `node_create_subject` rebuild adds an `'eigen'` node nested under the cortex; exercise the Delete callback path (file gone, entry removed, `bst_get('EigenFile')` no longer resolves); cleanup (also remove the operator node created by find-or-create). GUI tree rendering = user's visual check.

- [ ] **Step 4: Commit** — `feat(tree): nest eigen nodes under surface + Compute eigenmodes / View / Delete menus`.

---

## Final verification
- [ ] All three eigen variants compute into nodes on the canonical cortex; nodes nest under the surface; View/Delete work; `bst_get('AnyFile')` resolves manifold + operator + eigen; find-or-create builds/reuses the operator correctly; B-orthonormality holds per variant. Restore protocol.
- [ ] Run `test_tess_eigen`, `test_tess_operators`, `test_tess_manifold` — all green.
- [ ] Then **superpowers:finishing-a-development-branch**; SP4 (retire old TessMat eigenmode fields/IO + migrate leadfield/inverse consumers) is the next sub-project.

## Notes for the implementer
- Faithful replication of merged SP2 — read each SP2 template and match its idioms; the only genuinely new logic is `eigenmat`, the find-or-create + eigensolve in `tess_eigen`, the `local_ritz_basis` nested fn, and the `AnyFile`/`EigenFile` wiring.
- Per-hemi `Phi`/`Lambda`/`GlobalVertices` are 1×2 **cell** arrays (like `operatormat`'s `Operator`/`Mass`).
- No Java build; `'eigen'` renders with the default icon.
- DB-mutating tests must restore the protocol (`file_delete` + reload), and must also clean up any operator node created by find-or-create.
