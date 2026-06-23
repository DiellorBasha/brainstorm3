# bst_helmholtz — manifold geometry consolidation (design)

**Date:** 2026-06-23
**Status:** Approved (brainstorming) → ready for implementation plan
**Author:** Diellor Basha (with Claude)

## Goal

Consolidate the cortical Helmholtz–Hodge decomposition around two canonical
derived-anatomy nodes:

- the **`manifold_` node** as the single source of truth for *geometry* (face
  normals, areas, centroids, vertex positions/normals, frames), and
- the **`operator_` node** as the source of *differential operators* (Dirac `D`,
  LBO `K`/`M`, `gradFace`).

Concretely:

1. **Phase A** — add a per-face **area** field to the manifold node, computed by
   nxr-compute (so `tess_manifold` fetches *all* geometry, leaving no geometric
   primitive recomputed by consumers).
2. **Phase B** — rename `bst_dirac_helmholtz` → **`bst_helmholtz`** with a
   `switch`-based vertex/face dispatch, and rewire it to read geometry from the
   manifold node instead of recomputing it from raw vertices. Fold
   `bst_dirac_helmholtz_face` in as the `face` branch and delete it.

Phase A gates Phase B (the face area must exist in the manifold before the
Helmholtz rewire can consume it).

## Background / current state

`toolbox/math/bst_dirac_helmholtz.m` (vertex, 259 lines) and
`bst_dirac_helmholtz_face.m` (face, 150 lines) are a matched pair implementing the
same verb contract (`Prepare`/`Frame`, shared output struct) over two domains.
Both **recompute per-face geometry from raw vertices** in `Prepare`:

- face normals via `cross(e1,e2)` then a **flip heuristic** (`Nf·VertNormals < 0`)
  to disambiguate winding — a workaround for not having canonical normals;
- face areas via `‖cross‖`, stored under a variable named `Af` that means **2×area**
  in the vertex file but **true area** in the face file (a latent naming hazard,
  currently harmless only because every consumer is scale-invariant);
- centroids (face file) via the mean of the three vertices.

Meanwhile the `manifold_` node (`db_template('manifoldmat')`, built by
`tess_manifold` from `nxr_compute('facets', …)`) already stores canonical,
gauge-consistent per-face geometry, and `bst_face_leadfield` already consumes it
(`M.Embedded(hh).face.normal/.centroid`, "not tess_tangents"). The Helmholtz pair
predates / bypasses that backbone.

**Confirmed gap:** the manifold node has **no face-area field**. A live-node dump
shows the only area in the entire tree is `Intrinsic.vertex.dualArea` (per-vertex
lumped mass, not per-face). `Embedded.face` = `{normal, grid, centroid}`. nxr
*internally* computes `faceAreas` (geometry-central `requireFaceAreas`;
`geometry_bundle.cpp` fills `g.faceAreas`) — it is simply not exposed on the
`Embedded` struct.

## Design decisions (with rationale)

| Decision | Choice | Rationale |
|---|---|---|
| Where face area lives | `Embedded.face.area` in the manifold node, emitted by nxr | Area is pure geometry → belongs in the manifold, not borrowed from operator-node masses. Matches the "create the missing nxr operator, no workarounds" rigor directive. |
| Stale-schema manifold nodes | **Recompute on access** (schemaVersion gate in `tess_manifold`) | Keeps the on-disk node genuinely canonical; cost is one nxr recompute per surface on first access. |
| Merge shape | Domain-dispatching orchestrator (one `switch`, separate vertex/face bodies) | Keeps the two different algorithms readable; avoids `if isFace` spaghetti and avoids two files. |
| I/O boundary | `bst_helmholtz` stays **I/O-free**; caller resolves + passes `ManifoldMat` + operator node | Consistent with the math/ purity established for this module (callers own `in_*`/`bst_get`/`tess_manifold`). |
| Geometry vs operators | Geometry (normals/areas/centroids/positions) ← manifold; operators (`D`, `K`/`M`, `gradFace`) ← operator node | The operator/geometry split the rest of the pipeline already uses. |
| Rename strategy | Hard rename `bst_dirac_helmholtz` → `bst_helmholtz`; delete `_face`; update all callers | Matches prior clean renames in this codebase (no deprecated shims). |

## Phase A — `Embedded.face.area` via nxr-compute

### nxr-compute (C++), repo `~/workspace/research/code/nxr-compute`

1. `include/nxr/facets.h` — add to `EmbeddedFacet::FaceView`:
   `Eigen::VectorXd area() const;  // [nF]`
2. `src/facets.cpp` — implement alongside the existing accessors
   (`FaceView::normal` → `geometry::frames(m).normals`, `FaceView::centroid` →
   `lightGeometry().faceCentroids`):
   `Eigen::VectorXd EmbeddedFacet::FaceView::area() const { return m.lightGeometry().faceAreas; }`
   (Confirm `MeshGeometry`/`lightGeometry()` exposes `faceAreas`; `geometry_bundle.cpp`
   already resizes/fills `g.faceAreas`. If not surfaced on `lightGeometry`, add it
   there or use the `requireFaceAreas()` pattern used by the operators.)
3. `bindings/mex/src/nxr_compute_mex.cpp` — in `buildEmbeddedStruct` (~line 1561),
   add `"area"` to the face field list `{"normal","grid","centroid"}` → `{…,"area"}`
   and `mxSetField(g, 0, "area", eigenVectorToMx(fv.area()));` (use/observe the
   existing Eigen→mx helper for a `VectorXd`). **Bump the `Embedded` `schemaVersion`.**
4. Rebuild the MEX and **copy the rebuilt binary into**
   `~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/` (managed-plugin
   stale-binary trap — BST otherwise loads the old MEX).
5. nxr C++ test (`test/test_facets.cpp`): assert `area()` has `nF` entries, all
   positive, and `sum(area)` equals the mesh total area within tolerance.

### Brainstorm (MATLAB)

6. `tess_manifold.m` — define the required `Embedded.schemaVersion` constant; in the
   find-or-load reuse path (around the `bst_get('ManifoldFileForSurface', …)` block),
   treat a cached/loaded node whose `Embedded.schemaVersion` is below the required
   version as a **cache miss → recompute** (do not reuse). No change to the wholesale
   `S.Embedded` store — the new `area` field flows through automatically.
7. MATLAB test: a recomputed manifold node has `Embedded.face.area` `[nF×1]`, all
   positive, summing to the hemisphere area; a stale-schema node triggers recompute.

### Scoped cleanup (Phase A)

8. `bst_face_leadfield.m:196` — stop recomputing `A_f` locally; read
   `M.Embedded(hh).face.area`. (Reinforces single-source-of-truth; small, optional —
   drop if scope tightens.)

## Phase B — `bst_dirac_helmholtz` → `bst_helmholtz`

### File-level

- Rename `toolbox/math/bst_dirac_helmholtz.m` → `toolbox/math/bst_helmholtz.m`.
- Fold `bst_dirac_helmholtz_face.m` in as the `face` branch; **delete** `_face.m`.

### Public API

```
Op = bst_helmholtz('Prepare', OperatorNode, ManifoldMat, 'Domain', 'vertex'|'face')
Ht = bst_helmholtz('Frame',   Op, J [, withCores])      % routes on Op.Domain
H  = bst_helmholtz('Decompose', OperatorNode, ManifoldMat, J, 'Domain', …)  % vertex whole-series convenience
```

- `Prepare` tags `Op.Domain`; `Frame` dispatches on it (no shape-sniffing of `J`).
- The vertex-only whole-series `Decompose` convenience is retained under the vertex
  branch. The dead `PoissonSolve` verb is **removed** (only `test_dirac_helmholtz.m`
  used it; that assertion moves onto the `Prepare`/`Frame` path or the shared solver
  helper).

### `Op`-field mapping — vertex branch

| `Op` field | New source | Was |
|---|---|---|
| `D` | operator node `FirstOrder.Intrinsic{hh}` | (same) |
| `Nf` (face normal) | `Manifold.Embedded(hh).face.normal` | `cross`+`VertNormals` flip — **deleted** |
| face area (for `Wfv`) | `Manifold.Embedded(hh).face.area` | `‖cross‖` (`Af`=2A) — **deleted** |
| `M`, `cholK`, `free`, `totMass` | operator-node LBO `K`/`M` + local factorization | (same; factor stays local) |
| `Gx`/`Gy`/`Gz` (FEM P1 gradient) | assembled locally **from manifold geometry** (normals, areas, vertex positions) | assembled from recomputed geometry |
| `Wfv` (face→vertex incidence) | assembled locally from topology + manifold area | local (recomputed area) |
| `VtxH` (vertex positions) | `Manifold.Embedded(hh).vertex.position` | `Vtx(vH,:)` |
| `VnH` (vertex normals, sub-vertex chart) | `Manifold.Embedded(hh).vertex.normal` | `Surf.VertNormals(vH,:)` |
| `NbH` (1-ring, core detection) | local from topology | (same) |

Net: the submesh-reconstruction prologue, the cross-product/flip, and all geometry
recomputation are removed; `Surf.VertNormals` dependency is eliminated. The FEM
gradient operator is still *assembled* in `Prepare` (it is an operator, not raw
geometry) but consumes only manifold geometry.

### `Op`-field mapping — face branch

| `Op` field | New source | Was |
|---|---|---|
| `G` (`gradFace`) | operator node `FaceAux.GradFace` (Hodge-Face/Dirac-Face) | fresh `nxr_compute('create'/'operators')` — **deleted** |
| `Nf` | `Manifold.Embedded(hh).face.normal` | `cross`+flip |
| `Af`/`WF` (true area) | `Manifold.Embedded(hh).face.area` | `‖cross‖/2` |
| `Cf` (centroid) | `Manifold.Embedded(hh).face.centroid` | mean of 3 verts |
| `SkewG` | local: `i_block_rotation(Nf)·G` (manifold `Nf`) | (same, manifold normal) |
| `cholA`, `freeA` | local factorization of coupled `A` | (same; factor stays local) |
| `NbF` (dual adjacency, core detection) | local from topology | (same) |

The face branch now requires a face operator node providing `GradFace`; the caller
resolves it. **Implication:** face tests/benchmarks must create-or-resolve a face
operator node instead of relying on the in-`Prepare` nxr build.

### Internal cleanup

- One shared `i_pinned_solve(chol, rhs, free)` helper unifies the duplicated
  cached-factor reuse (vertex dual-Poisson + face coupled solve).
- Normalize `i_empty_cores` and the `bst_vortex_persistence` core-detection loop into
  shared helpers where the two domains genuinely coincide; keep domain-specific
  geometry (sub-vertex quadratic fit vs centroid positions) separate.

### Callers

- `toolbox/gui/view_helmholtz.m` and `toolbox/process/functions/process_vortex_track.m`:
  resolve the manifold via `tess_manifold(SurfaceFile)`, pass `ManifoldMat` + the
  operator node + `'Domain'` into `bst_helmholtz`.
- Extract the duplicated local `i_op(SurfaceFile, variant)` operator-loader (currently
  byte-identical in both callers) into one shared helper. (Scoped cleanup — drop if
  scope tightens.)
- Production callers use the **vertex** domain; the **face** domain is exercised by
  tests/benchmarks only (face core detection is still deferred from the GUI).

### Tests/benchmarks to re-point

`dev/tests/test_dirac_helmholtz.m`, `dev/tests/test_dirac_helmholtz_face.m`,
`dev/tests/test_helmholtz_view.m`, `dev/benchmarks/bench_dirac_face_helmholtz.m`,
`dev/benchmarks/bench_dirac_streamribbon_real.m` → `bst_helmholtz`.

## Testing strategy

- **Phase A:** nxr C++ area test; MATLAB manifold-area test; stale-schema-recompute
  test.
- **Phase B:** both existing Helmholtz tests (re-pointed) stay green — the key
  invariants are the **HarmFrac→0 round-trip** and **B-/W_F-orthonormality**, plus
  core/source detection. Add a **parity test**: with geometry sourced from the manifold
  vs the old cross/flip recompute, assert per-face normals agree (up to the documented
  sign convention) and `HarmFrac`, `Curl`, `Div` are unchanged within tolerance on a
  real surface.

## Top risk

The manifold `face.normal` is "outward, trivial-connection sign"; the old code derived
the sign via the `VertNormals` flip. Divergence (`imag·n`) and the solenoidal field
(`n×grad`) depend on that sign, so a per-face sign mismatch would silently corrupt the
decomposition. **The parity test is the correctness guard** for the geometry swap; if
signs disagree, reconcile the convention (or apply a one-time sign alignment) before
proceeding.

## Sequencing

1. Phase A complete and verified (nxr rebuilt + installed, manifold recompute yields
   `area`, tests green).
2. Then Phase B (rename + dispatch + geometry rewire + caller/test updates), gated on
   the Phase B parity test.

## Out of scope / follow-ups

- Moving the vertex FEM P1 gradient into the operator node (it stays assembled in
  `Prepare`, sourced from manifold geometry).
- Bringing the face domain into the GUI (face core detection remains deferred).
- Any change to the inverse/eigen pipelines.
