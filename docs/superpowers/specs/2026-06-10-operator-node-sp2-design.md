# Derived-Anatomy Nodes — SP2 (`operator_` node) — Design

**Date:** 2026-06-10
**Status:** Design (approved) — pending spec review → plan
**Repo:** `brainstorm3` (MATLAB). Builds directly on SP1 (`2026-06-10-derived-anatomy-nodes-design.md`, merged).

## Goal

Add the second derived-anatomy node type — **`operator_`**, a surface child holding an operator + its mass matrix — and the assembler `tess_operators`. This closes the dependency loop for SP3: `tess_eigen` will read an operator node and eigensolve it. Scope is exactly the three operators the eigensolves need: **Laplace-Beltrami, Connection Laplacian, Dirac** (Mass/DEC/other variants deferred).

## Decisions (from brainstorming)

- **Order:** operator node (SP2) before eigen node (SP3) — eigen consumes the operator.
- **Bundled mass:** each `operator_` node stores the `{Operator A, Mass B}` pair for its generalized eigenproblem, so `tess_eigen` runs `eigs(A,B)` directly. A standalone Mass node is deferred.
- **Relationship:** `operator_` and `eigen_` are both **children of the surface**; the eigen node (SP3) stores an `OperatorFile` **path reference** to the operator it was computed from (Brainstorm's reference-by-path style, like results→head model).
- **Per-hemisphere:** operators (and thus LBO too) are assembled per hemisphere via `tess_hemisplit` (atlas L/R, never `conncomp`) — consistent with connection/Dirac and the project rule. This is a deliberate change from the old whole-mesh/`conncomp` `tess_eigenmodes`.

## Architecture — mirror SP1's `manifold_` node

SP1 already established and merged the derived-node infrastructure: surface child lists (`.Operator` exists and is backfilled by the `db_update` 5.04 migration), `node_create_subject` nesting, `tree_callbacks` menus, `db_add_*`/`bst_get(*File)` resolution, the `'manifold'` Java node type rendering with a default icon, and the `manifold` file-type wiring. SP2 replicates that pattern for `operator`.

### New / changed pieces

**1. `db_template('operatormat')`** (file content):
```
Comment, ParentSurface,
Variant,        % 'Laplace-Beltrami' | 'Connection Laplacian' | 'Dirac'
Operator,       % 1x2 per-hemisphere sparse A
Mass,           % 1x2 per-hemisphere sparse B (matches A's eigenproblem)
GlobalVertices, % 1x2 per-hemisphere scatter maps
Provenance      % Backend='nxr', NxrVersion, Variant, params (e.g. Tau), ComputeDate
```
The bundled `B` per variant: LBO → scalar vertex mass `M`; Connection → the connection-Laplacian's mass `M`; Dirac → `kron(Mg, I4)`. So a consumer reads `(Operator(hh), Mass(hh))` and solves the generalized problem with no extra assembly.

**2. File type `operator`** — add to `file_gettype` (`operator_` prefix → `'operator'`), `file_short`, `file_fullpath` (anat-relative), mirroring `manifold`. **Also fold in the deferred `bst_get('AnyFile')` support for BOTH `manifold` and `operator`** so the generic "View file contents" menu works (SP1 left this erroring "File type not recognized").

**3. `db_add_operator(iSubject, ParentSurfaceFile, OperatorMat, Comment)`** + **`bst_get('OperatorFile', FileName)`** — mirror `db_add_manifold` / `bst_get('ManifoldFile')` exactly (ProtocolSubjects round-trip incl. `iSubject==0`, append to the parent surface's `.Operator` list, `file_unique` name in the anat dir, `db_save`, `UpdateNode`).

**4. `node_create_subject`** — add an `.Operator` child-nesting loop next to the `.Manifold` loop (same `CreateNode('operator', …)` + `isfield` guard pattern).

**5. `tree_callbacks`** — surface popup: **"Compute operator ▸ {Laplace-Beltrami | Connection Laplacian | Dirac}"** submenu (each calls `tess_operators(filenameRelative, <variant>)`). `case 'operator'` node popup: **View** (field-name inspector stub, like the manifold View) + **Delete** (mirror `ManifoldDelete_Callback`). Add `'operator'` to the `filenameFull` SUBJECTS switch.

**6. `tess_operators(SurfaceFile, OperatorName, ...)`** (new):
- Switch `OperatorName` → variant. Guards (nxr-compute installed; Structures L/R atlas; `Reg.Sphere` for the trivial gauge — needed by the Dirac/connection trivial gauge).
- Per hemisphere (`tess_hemisplit`; local submesh; `nxr_compute('create', Vloc, Floc)`):
  - **Laplace-Beltrami:** `A = operators(h,'laplacian','cotan')`, `B = operators(h,'mass','galerkin')`.
  - **Connection Laplacian:** `A = operators(h,'laplacian','connection')` (complex K, in the trivial gauge with FS-pole singularities), `B = operators(h,'mass','galerkin')`.
  - **Dirac:** `A = operators(h,'dirac',Tau)`, `B = kron(operators(h,'mass','galerkin'), speye(4))`.
  - `destroy`; attach `GlobalVertices`.
- Assemble `OperatorMat = db_template('operatormat')` (1×2 `{Operator, Mass}`, `Variant`, `Provenance`); `db_add_operator` unless `NoSave`. Options: `Tau` (0.5, Dirac/connection trivial gauge), `NoSave`, `ForceRecompute`.
- Return `OperatorMat`.

## Sets up SP3
`tess_eigen(SurfaceFile, OperatorName, ...)` will: find-or-create the operator node (`tess_operators`), read its `(Operator, Mass)` per hemi, `eigs(A,B)` + B-orthonormalize (Rayleigh-Ritz for Dirac's degenerate multiplets — see [[dirac-eigenmode-leadfield]]), and write an `eigen_` node (`Variant` Scalar/Vector2/Vector3) carrying `OperatorFile` = the operator node's path.

## Tests / validation
- `db_template('operatormat')` schema; `file_gettype('operator_*')='operator'`; existing types intact.
- `db_add_operator` + `bst_get('OperatorFile')` live round-trip (create/resolve/cleanup), incl. `bst_get('AnyFile')` now resolving manifold+operator.
- `tess_operators` per variant on the real cortex: produces a 1×2 `{Operator,Mass}` node with correct shapes/symmetry per variant (LBO `A=Aᵀ` real `[nVₕ²]`; Dirac `A` `[4nVₕ²]` symmetric, `B=kron(Mg,I4)`); node nests under the surface; `eigs(A,B,6,'smallestabs')` returns real eigenvalues (smoke). Cleanup after.
- Tree rendering / menus: user visual check in the GUI.

## Out of scope (SP2)
- `tess_eigen` / `eigen_` node (SP3).
- Mass/DEC/other operator variants (later).
- Retiring old fields/functions/IO and frame/tangent consumers (SP4).
