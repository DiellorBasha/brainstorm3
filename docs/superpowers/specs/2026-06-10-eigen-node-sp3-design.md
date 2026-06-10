# Derived-Anatomy Nodes — SP3 (`eigen_` node) — Design

**Date:** 2026-06-10
**Status:** Design (approved) — pending spec review → plan
**Repo:** `brainstorm3` (MATLAB). Builds directly on SP1 (`manifold_`, merged) and SP2 (`operator_`, merged on `development`).

## Goal

Add the third derived-anatomy node type — **`eigen_`**, a surface child holding the eigenbasis of an operator — and the assembler `tess_eigen`. Each eigensolve is its **own file/node** (one `eigen_*.mat` per solve), with the variant carried as a label + comment + filename, and a **path reference to the operator node** it was computed from. This replaces the originally-sketched single nested `Eigen.{Scalar,Vector2,Vector3}` struct: the file-based DB makes one node per eigensolve the natural unit (parallel to how `operator_` nodes carry their `Variant`).

## Decisions (from brainstorming)

- **One eigen file/node per eigensolve** (not a single node with nested Scalar/Vector2/Vector3 fields). Variant by label + comment + file ending + `OperatorFile` path reference.
- **Scope = eigen node only.** Migrating the eigenmode leadfield/inverse consumers and retiring the old TessMat eigenmode fields is deferred to SP4.
- **Find-or-create operator.** `tess_eigen(SurfaceFile, OperatorName, ...)` reuses an existing operator node of the matching `Variant` on the surface, or builds one via `tess_operators` if absent, then eigensolves it and stores the resulting `OperatorFile` reference.
- **Variant labels mirror the operator names** — `'Laplace-Beltrami'` | `'Connection Laplacian'` | `'Dirac'` (the eigenvector field dimensionality follows from the variant).
- **Rayleigh–Ritz for Dirac AND Connection** degenerate multiplets, included as a `local_`-prefixed **nested function inside `tess_eigen.m`** (Brainstorm convention — no standalone shared helper; the existing copy in `tess_dirac_eigenmodes.m` stays as-is). The nested copy is **complex-safe** (Hermitian transpose) for the complex Connection operator. Only LBO uses plain B-orthonormalization.
- **Per-hemisphere**, atlas split via `tess_hemisplit` (never `conncomp`) — consistent with SP2.
- **Canonical Mesh Rule** applies: any `nxr_compute('create')` reached transitively (via `tess_operators` find-or-create) goes through `nxr_safe_create` on a `bst_canonical_cortex` mesh; tests never hand-build meshes.

## Architecture — mirror SP2's `operator_` node

SP1/SP2 established the derived-node infrastructure (surface child lists incl. `.Eigen`, backfilled by the `db_update` 5.04 migration; `node_create_subject` nesting; `tree_callbacks` menus; `db_add_*`/`bst_get(*File)` resolution; the Java node type rendering with a default icon; file-type wiring; `bst_get('AnyFile')` support). SP3 replicates that pattern for `eigen`.

### New / changed pieces

**1. `db_template('eigenmat')`** (file content):
```
Comment, ParentSurface,            % back-reference (path)
OperatorFile,                      % path ref to the operator node it solved   ← KEY
Variant,                           % 'Laplace-Beltrami' | 'Connection Laplacian' | 'Dirac'
Phi,            % 1x2 per-hemi cell: eigenvectors
                %   LBO        [nVh x K]  real
                %   Connection [nVh x K]  complex
                %   Dirac      [4nVh x K]
Lambda,         % 1x2 per-hemi cell: eigenvalues [K x 1]
K,                                 % modes requested
GlobalVertices, % 1x2 per-hemi cell: scatter maps (copied from the operator node)
Provenance      % Backend='nxr', NxrVersion, Variant, K, Tau (Dirac), Ortho='Rayleigh-Ritz'|'B-orthonormal', ComputeDate
```
`Phi`/`Lambda`/`GlobalVertices` use 1×2 **cell** arrays (per-hemi matrices of differing dimensions can't share a struct-array field — same choice as `operatormat`). `db_template('surface').Eigen` already exists (wired in SP1) and is populated now.

**2. File type `eigen`** — add to `file_gettype` (`eigen_` prefix → `'eigen'`, before `_tess`), `file_short`, `file_fullpath` (anat-relative), and `bst_get('AnyFile')` (route via `EigenFile`), mirroring `operator`.

**3. `db_add_eigen(iSubject, ParentSurfaceFile, EigenMat, Comment)`** + **`bst_get('EigenFile', FileName)`** → `[sSubject, iSubject, iSurface, iEigen]` — exact mirror of `db_add_operator` / `bst_get('OperatorFile')` (ProtocolSubjects round-trip incl. `iSubject==0`, append to the parent surface's `.Eigen` list, `file_unique` name, `db_save`, `UpdateNode`).

**4. `node_create_subject`** — add an `.Eigen` child-nesting loop next to the `.Operator` loop (`CreateNode('eigen', …)` + `isfield` guard).

**5. `tree_callbacks`** — surface popup: **"Compute eigenmodes ▸ {Laplace-Beltrami | Connection Laplacian | Dirac}"** submenu (each → `tess_eigen(filenameRelative, <variant>)`), placed near "Compute operator", cortex + write guarded. `case 'eigen'` node popup: **View** (`EigenView_Callback`, field inspector) + **Delete** (`EigenDelete_Callback`, mirror `OperatorDelete_Callback`). Add `'eigen'` to the `filenameFull` SUBJECTS switch.

**6. `tess_eigen(SurfaceFile, OperatorName, ...)`** (new):
- Map `OperatorName` → `Variant` (case-insensitive; clear error on unknown).
- Guards: nxr-compute installed; Structures L/R atlas (needed by the find-or-create path → `tess_operators` → `tess_hemisplit`).
- **Find-or-create operator:** scan `sSubject.Surface(iSurface).Operator` for an entry with matching `Variant`; if found, resolve its file; else `tess_operators(SurfaceFile, OperatorName [, 'Tau',Tau])`. Load the operator node's per-hemi `(Operator A, Mass B)` cells and its `GlobalVertices`.
- **Per hemisphere:** `nModes = K` (Dirac: round to a multiple of 4 and over-fetch for the 4-fold multiplets, reusing the cap logic from `tess_dirac_eigenmodes`); `[V, D] = eigs(A, B, nRequest, 'smallestabs')`; B-orthonormalize. **Dirac and Connection → `local_ritz_basis`** (nested function) for their degenerate multiplets; LBO → plain B-orthonormalization. Keep the first `K` modes; `Phi{hh}=V`, `Lambda{hh}=diag(D)`.
- Assemble `EigenMat = db_template('eigenmat')` (1×2 `Phi`/`Lambda`/`GlobalVertices`, `OperatorFile`, `Variant`, `ParentSurface`, `Provenance`). `db_add_eigen` unless `NoSave`. Comment e.g. `sprintf('%s eigenmodes (K=%d)', Variant, K)`. Return `EigenMat`.
- Options: `K` (default **400**), `Tau` (0.5, forwarded to operator creation for Dirac), `NoSave`, `ForceRecompute`. License header.

## Tests / validation
- `db_template('eigenmat')` schema; `file_gettype('eigen_*')='eigen'`; existing types intact; `bst_get('AnyFile')` resolves manifold/operator/eigen.
- `db_add_eigen` + `bst_get('EigenFile')` live round-trip (create/resolve/cleanup).
- `tess_eigen` per variant on the canonical cortex (`bst_canonical_cortex(20484)` — never a hand-built mesh): produces 1×2 `Phi`/`Lambda` of the right shapes; **B-orthonormality** `Phiᵀ B Phi ≈ I` per hemi; `Lambda` ascending & real-ish (Dirac → 4-fold clusters); the stored `OperatorFile` resolves via `bst_get('OperatorFile')`; **find-or-create** builds the operator node when absent and reuses it when present; save-path creates an `eigen_` node nested under the surface. Restore the protocol (`file_delete` + `db_reload_database('current')`).
- Tree rendering / menus: user visual check.

## Out of scope (SP3)
- Migrating leadfield/inverse consumers (`bst_dirac_eigenmode_leadfield`, `bst_inverse_eigenmodes`, …) to read eigen nodes — SP4.
- Retiring old TessMat `Eigenmodes`/`ConnEigenmodes`/`DiracEigen` fields + `in_/out_tess_*eigenmodes` IO — SP4.
- Rich eigenbasis browser/visualization — later.
