# Cortical-Flow Distillation — Milestone 1 Design

**Date:** 2026-08-07
**Status:** Approved design, pending implementation plan
**Scope:** Milestone 1 of distilling the proof-of-concept cortical-flow methods into
clean, self-contained Brainstorm commits: icosphere prerequisite, operator
primitives, `tess_eigen`. Later phases (filtering, visualization/detection, Dirac
kernel projection) get their own specs after the Milestone 1 gate.

## 1. Goal

Replicate the validated PoC methods (currently built on the nxr-compute MEX
plugin, custom derived-anatomy DB nodes, and the bst-java fork) as pure-MATLAB,
plugin-free, self-contained commits on a clean feature branch, verified in
stages against the PoC implementations.

## 2. Decisions (settled during brainstorming)

| Question | Decision |
|---|---|
| Destination | New clean branch (e.g. `feature/cortical-flow-core`) off upstream `master`, PET-deliverable style. Commits read as Diellor's own: no Claude trailers, no plan/design docs on the branch. Specs + verification harness stay on `development`. |
| nxr-compute | Retired. Pure-MATLAB reimplementation; nxr used only as an off-branch validation oracle during development. |
| Milestone 1 operator scope | LBO + Dirac. Connectome variants and Connection Laplacian deferred. |
| Storage | Embedded in the surface `.mat` file (accepting larger files / slower first load). No new node types, no bst-java changes, no `db_update` migration. |
| Compute vs. store | Operators are cheap → computed on demand, never stored. Eigenmodes are expensive → the only stored artifact (`SurfaceMat.Eigen`). |
| Verification gates | Per phase: (1) numerical parity vs nxr-built results on a real cortex, (2) analytic checks on a clean icosphere, (3) live Brainstorm demo script the user runs. |
| Structure | Approach A: layered phases with a user-reviewed gate after each (Phase 0 icosphere → Phase 1 operators → Phase 2 eigen = Milestone 1). |
| Dirac inverse-kernel track | Separate later phase on the same branch. It is a *projection* of the existing unconstrained MNE ImagingKernel (vector form in Dirac eigenmodes, norm form in LBO eigenmodes), NOT a port of `bst_inverse_dirac`. |

## 3. Architecture

### 3.1 Primitives — pure functions, no I/O

- `tess_massmatrix(Vertices, Faces)` → sparse Galerkin (consistent FEM) mass
  matrix B. Parity target: nxr `mass/galerkin`, near machine precision.
- `tess_laplacian(Vertices, Faces)` → sparse cotangent stiffness A.
  Parity target: nxr `laplacian/cotan`.
- `tess_dirac(...)` → Dirac operator. **Gated on Open Question 2** (exact
  variant); not implemented until that flag is resolved.

Vectorized sparse assembly (3-corner-angle accumulation). Negative cotan
weights from obtuse triangles are kept (matches nxr). Inputs validated
(NaN/Inf, degenerate faces).

### 3.2 Pencil layer — `tess_operators(SurfaceFile, OperatorName)`, compute-only

- Splits hemispheres from the **Structures atlas** via `tess_hemisplit`
  (never `conncomp`); errors clearly if the atlas is missing or a hemisphere
  is not a closed 2-manifold (pointing at the icosphere import path).
- Calls the primitives per hemisphere; returns the pencil `(A, B){1x2}`.
- Also accepts the `Operator` recipe struct stored in an `Eigen` entry (see
  3.4) so a stored basis's pencil can be reassembled exactly.
- No saving, no caching, no DB interaction, no `NoSave`/`ForceRecompute`/
  `Interactive` options. Assembly is sub-second.

### 3.3 `tess_eigen(SurfaceFile, OperatorName, ...)` — the stored artifact

- Find-or-create against `SurfaceMat.Eigen.<Variant>`: reuse when the stored
  `Operator` recipe matches the request and stored `nModes` >= requested
  (truncate); otherwise recompute and **replace** that variant's slot.
- Solves `A*phi = lambda*B*phi` per hemisphere; B-orthonormalizes;
  Rayleigh-Ritz re-orthogonalization within degenerate Dirac multiplets.
  Solver strategy is **provisional pending Open Question 1**.
- Errors on non-convergence; never silently returns fewer modes.

### 3.4 Data model — embedded `Eigen` field, named variants

One new optional field in the surface file. Variant names map to valid MATLAB
identifiers through a single canonical helper
(`'Laplace-Beltrami'` → `LaplaceBeltrami`, `'Dirac'` → `Dirac`; later
`LBConnectome`, `DiracConnectome`, `ConnectionLaplacian`).

```matlab
SurfaceMat.Eigen
  .LaplaceBeltrami
      .Phi      {1x2}   % {lh, rh} per-hemisphere B-orthonormal eigenvectors
      .Lambda   {1x2}   % ascending; lambda_1 ~ 0
      .nModes           % modes per hemisphere (default 400)
      .Operator         % machine-readable recipe: everything tess_operators
                        %   needs to reassemble the SAME pencil:
                        %   .Name ('Laplace-Beltrami'), .Tau ([]),
                        %   .Assembly (function versions, mass type)
      .Solver           % how the eigensolve ran: .Method, shift/tolerances
                        %   (per Open Question 1), .Date
  .Dirac
      .Phi      {1x2}   % 4*nVh rows (quaternion components, vertex-major);
                        %   numeric class TBD per Open Question 2
      .Operator         % .Name='Dirac', .Tau, .Assembly (per Open Question 2)
      ...
```

Invariants and conventions:

- **Row convention (never stored):** hemi-local row k = k-th *sorted* global
  vertex index of that hemisphere per the Structures atlas — exactly
  `tess_hemisplit` output order. Documented in both function headers. No
  `GlobalVertices` field: hemisphere membership is derived at run time.
- **One slot per variant, by construction.** Recomputing with different
  parameters (e.g. a different Dirac Tau) replaces the slot; parameter-distinct
  bases of the same variant cannot coexist. Deliberate cache semantics —
  prevents silent accumulation of stale bases.
- **Reproducibility invariant:**
  `tess_operators(SurfaceFile, S.Eigen.<v>.Operator)` reproduces the pencil;
  `S.Eigen.<v>.Solver` pins the eigensolve. Traceability never depends on
  parsing History prose.
- **Consistency guard at reuse:** per-hemisphere `Phi` row counts must sum to
  `length(Vertices)` (times 4 for Dirac) and the Structures atlas must be
  present; mismatch forces recompute. (Geometry edits create new surface files per Brainstorm
  convention, so in-place staleness is an edge case, not the norm.)

### 3.5 File writes and DB — the mrimask pattern

Follows the established in-place field-addition precedent
(`bst_memory` mrimask cache, `tess_addsphere`):

1. `s.Eigen = ...; bst_save(SurfaceFile, s, 'v7', 1)` — append/update mode,
   only the `Eigen` variable is rewritten.
2. `bst_history('add', SurfaceFile, 'eigen', desc)` — narrative row, e.g.
   `'Laplace-Beltrami: 400 modes/hemisphere, solver=..., tess_eigen v1.0'`.
3. If the surface is loaded, patch `GlobalData.Surface(iSurf)` in memory
   directly — no unload, no view disruption.
4. **No DB calls.** FileName/Comment/SurfaceType are unchanged; the database
   and tree see the same surface. Embedded data survives `db_reload_subjects`
   by construction.

`in_tess_bst` uses plain `load()` (no field whitelist), so the new field flows
through every existing consumer untouched.

Size budget: ico5 LBO K=400/hemi ≈ 66 MB; Dirac ≈ 262 MB if real double
(class TBD per Open Question 2). `-v7` has a 2 GB per-variable ceiling and
`Eigen` is one variable — comfortable margin at Milestone 1 scope; fallback if
ever crossed is splitting per-variant top-level fields. Cost accepted: 3D
views pay the extra read on first load of the surface file.

## 4. Phases and gates

### Phase 0 — Icosphere prerequisite (extraction, not rewrite)

Re-derive as clean commits from the `development` branch (already pure MATLAB):

1. `tess_downsize` `'icosphere'` method — FreeSurfer/MNE-style per-hemisphere
   uniform resampling via the registration sphere.
2. `tess_manifold` consolidated manifold check/repair — guarantees each
   hemisphere is a closed 2-manifold (the property the operators depend on).
3. Import wiring — FreeSurfer import method option + BIDS import defaults.

**Gate 0:** import an OMEGA/tutorial subject; confirm ico5 cortex and
per-hemisphere manifold pass.

### Phase 1 — Operator primitives + pencil layer

`tess_massmatrix`, `tess_laplacian`, compute-only `tess_operators` (LBO path).
The Dirac commit within this phase waits on Open Question 2.

**Gate 1:** parity vs nxr-built A, B on a real cortex (near machine
precision); analytic icosphere checks (mass row sums = vertex areas, total =
sphere area; `L*ones = 0`; symmetry).

### Phase 2 — `tess_eigen` (= Milestone 1 gate)

Embedded-storage eigensolve per 3.3/3.4/3.5. Solver settings provisional
pending Open Question 1.

**Gate 2:** analytic sphere spectrum `lambda = l(l+1)/R^2` with
multiplicities 2l+1 on a clean icosphere; subspace correlation vs the nxr-era
eigenbasis on a tutorial cortex (compared as subspaces, not columns — signs
and multiplet mixing are gauge); live GUI script displaying eigenmodes on the
cortex.

## 5. Verification harness

Lives on `development` (never on the clean branch): one driver script per
gate producing a saved report the user reviews. nxr-compute remains installed
on the development side purely as the comparison oracle.

## 6. Open questions (pinned; must be resolved at the stated boundary)

### Open Question 1 — Eigensolver strategy (blocks Phase 2 implementation)

The known conditioning bug (`eigs 'smallestabs'` with sigma=0 returns wrong
eigenbases on ill-conditioned pencils) means the clean `tess_eigen` must not
blindly copy the old code *or* an assumed shift-invert recipe. Resolve by
numerical experiment before Phase 2: exact shift choice; whether the LBO
pencil (Galerkin B) needs it at all or only Dirac/ill-conditioned pencils;
run-to-run reproducibility. Parity vs the nxr-era basis is done as subspace
correlation, since pre-fix results may themselves need re-validation.

### Open Question 2 — Exact Dirac variant (blocks the Phase 1 Dirac commit)

Two different operators are both called "Dirac" in the PoC lineage:

- the **intrinsic quaternionic relative-Dirac operator** `dirac(tau)`
  (Crane-style, assembled from face geometry) — what nxr builds and what the
  PoC `tess_operators` stores;
- the **ambient-flat lift** `cotanL (x) I4` — the flow work's validated
  anchor, found equal-or-better for MEG source vectors, and nearly free once
  the LBO pencil exists.

Different spectra, different assembly cost (quaternionic block assembly vs a
`kron`), different parity targets, different `Phi` numeric class/size.
Decide before the Phase 1 Dirac commit which is *the* Dirac operator on the
clean branch (or whether tau interpolates with ambient-flat as a special
case). The Dirac kernel-projection phase depends on this choice.

## 7. Later phases (outline only — each gets its own brainstorm + spec)

1. **Filtering framework** — eigenfilter/wavelet machinery operating on the
   embedded basis.
2. **Visualization / detection** — view/panel layer using stock node types
   only.
3. **Dirac kernel projection** — existing unconstrained MNE ImagingKernel
   expressed in eigenmode coordinates: vector kernel in Dirac eigenmodes,
   norm (scalar) form in LBO eigenmodes. Depends on Open Questions 1 and 2.

## 8. Non-goals (Milestone 1)

- No connectome operator variants, no Connection Laplacian.
- No new tree node types, java changes, or DB schema migrations.
- No GUI/panel work beyond the Gate 2 demo script.
- No port of `bst_inverse_dirac` (the Dirac-basis inverse solver).
- No changes to `in_tess_bst` / `bst_memory` load paths (no selective loading).
