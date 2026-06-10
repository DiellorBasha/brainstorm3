# nxr Facet Bundle → TessMat via `tess_frame` — Design

**Date:** 2026-06-09
**Status:** Design (pending review)
**Repos:** `nxr-compute` (C++/MEX) + `brainstorm3` (MATLAB)

## Goal

Bring nxr-compute's per-element manifold quantities into the Brainstorm surface
file, organized to mirror nxr's **representation-facet** taxonomy
(topology / embedded / intrinsic / extrinsic) plus the gauge, so Brainstorm can
use the precomputed frames and geometry/topology quantities — in particular the
ambient 3D vertex frame that lives in the same world coordinates as the leadfield.

A single self-contained `tess_frame` is the Brainstorm entry point: it loads the
surface, splits hemispheres, runs nxr per hemisphere, stores the facet groups on
the TessMat file, and returns the full-mesh `{U,V,N}` frame.

## Background & key findings (verified live, 2026-06-09)

- **There is no `embedded` MEX command** — not in the installed binary, the
  freshly built `main` binary, or any branch/commit. Verified by probing both
  binaries and `git grep` across all revisions. The dispatched data commands are
  `topology`, `geometry`, `gauge`, `bundle` (+ `operators` on `main`).
- The C++ **facet layer already exists** (`include/nxr/facets.h`,
  `src/facets.cpp`): typed views `topology()`, `embedded()`, `intrinsic()`,
  `extrinsic()`, `gauge()`, `operators()`. The work is marshalling these views to
  MATLAB structs — the data is already computed and cached.
- The flat `geometry` command **mixes** embedded/intrinsic/extrinsic per element
  and **omits** `vertex.position` (= input `V`), `vertex/face.normal`, and
  `principalDir`. Grouping by representation in C++ (not a fragile MATLAB
  field-name regroup) is therefore the better design and exposes the dropped
  members.
- The ambient frame is `embedded().vertex().grid()` = `e1 + i·e2` in world
  coordinates; `N = real(grid) × imag(grid)`. The trivial gauge realizes it as
  `exp(i·φ_v) · grid` via `gauge().vertex().rotation`.
- **MEX marshalling pattern** (`bindings/mex/src/nxr_compute_mex.cpp`):
  `buildTopologyStruct` / `buildGeometryStruct` / `buildGaugeStruct` use
  `mxCreateStructMatrix` + typed converters (`scalarToMx`,
  `indexVectorToMx1Based`, real/complex matrix converters); `cmdBundle` composes
  them. New facet builders follow the same shape.
- **Brainstorm native I/O** (e.g. `tess_addsphere`): `in_tess_bst(TessFile)` to
  load → modify → `bst_save(file_fullpath(TessFile), TessMat, 'v7')`, with
  `bst_history` for provenance. Hemisphere split via `tess_hemisplit` (atlas
  L/R) — **never `conncomp`**.

## Decomposition

Two phases; Phase 1 is a hard dependency of Phase 2 and ships first. Each phase
becomes its own implementation plan.

---

## Phase 1 — nxr-compute: facet data commands

**Repo:** `~/workspace/research/code/nxr-compute` (branch `main`).

### New struct builders

Source the typed facet views; marshal with existing converters. Each carries a
`schemaVersion` (== 1), matching `topology`/`geometry`/`gauge`.

**`buildEmbeddedStruct(h)`**

| Field | Type / size |
|---|---|
| `vertex.position` | `[nV×3]` double (raw input `V`) |
| `vertex.normal` | `[nV×3]` double |
| `vertex.grid` | `[nV×3]` complex (`e1 + i·e2`) |
| `face.normal` | `[nF×3]` double |
| `face.grid` | `[nF×3]` complex |
| `face.centroid` | `[nF×3]` double |

**`buildIntrinsicStruct(h)`**

| Field | Type / size |
|---|---|
| `vertex.dualArea` | `[nV×1]` double |
| `vertex.angleSum` | `[nV×1]` double |
| `edge.length` | `[nE×1]` double |
| `edge.cotanWeight` | `[nE×1]` double |
| `halfedge.cotanWeight` | `[nH×1]` double |
| `halfedge.transportAlong` | `[nH×1]` complex |
| `halfedge.transportAcross` | `[nH×1]` complex |

**`buildExtrinsicStruct(h)`**

| Field | Type / size |
|---|---|
| `vertex.curvature2RoSy` | `[nV×1]` complex (deviatoric `q`) |
| `vertex.meanCurvature` | `[nV×1]` double |
| `vertex.principalDir` | `[nV×3]` double (max principal direction) |
| `edge.dihedralAngle` | `[nE×1]` double |

### New dispatch commands (additive — `geometry`/`bundle` untouched)

- Standalone: `nxr_compute('embedded', h)`, `('intrinsic', h)`, `('extrinsic', h)`.
- Grouped: **`nxr_compute('facets', h, gaugeType[, opts])`** →
  `{Topology, Embedded, Intrinsic, Extrinsic, Gauge}`.
  - `Topology` / `Gauge` reuse `buildTopologyStruct` / `buildGaugeStruct`.
  - `opts` carries `singVerts` / `singValues` for the trivial gauge (same
    contract as `bundle`).
  - **Data-only** — no `.operators` (those remain on the `operators` command;
    out of scope here).

Register all four in the dispatcher and the "Available:" help string.

### Build / install

1. TDD in `bindings/mex/test/` (see below), build with `scripts/build.sh Release`.
2. Back up the installed binary
   (`~/.brainstorm/plugins/nxr-compute/nxr-compute-mex-r2023b/nxr_compute.mexmaca64`
   → `.bak.<date>`), copy the new `build/Release/nxr_compute.mexmaca64` in.
3. Smoke-check `nxr_compute('facets', h, 'trivial', opts)` in the live session.
4. (Optional) bump `cmdVersion` to `0.2.0` so Phase 2 can assert a minimum.

### Phase 1 tests (`bindings/mex/test/test_facets.m`)

- `facets` returns exactly `{Topology, Embedded, Intrinsic, Extrinsic, Gauge}`.
- `embedded.vertex.position == V` (input), element counts consistent
  (`nV/nE/nF/nH`).
- `embedded.vertex.grid` orthonormal: `|e1|=|e2|=1`, `e1·e2≈0`; and
  `embedded.vertex.normal ≈ e1×e2`.
- `extrinsic.vertex.principalDir` present and unit-norm (proves the facet path
  exposes what flat `geometry` drops).
- Trivial gauge on a genus-0 mesh: `sum(Gauge.singularity.indices) == 2`
  (Gauss-Bonnet).
- Standalone `embedded`/`intrinsic`/`extrinsic` are byte-identical to the
  corresponding sub-structs of `facets`.

**Acceptance:** new tests pass; the existing `topology`/`geometry`/`gauge`/
`bundle` tests still pass (back-compat); MEX builds and installs cleanly.

---

## Phase 2 — brainstorm3: self-contained `tess_frame`

**Repo:** `~/workspace/research/code/brainstorm3`. No shared helper
(`tess_store_perhemi` is **not** used and **not** deleted here).

### Contract

```
[U, V, N] = tess_frame(SurfaceFile, ...)
```

**Options:** `'Gauge'` (default `'trivial'`), `'Domain'` (`'vertex'` default;
`'face'` errors under trivial — `Gauge.face.rotation` is deferred/empty),
`'NoSave'`, `'ForceRecompute'`.

### Algorithm

1. `TessMat = in_tess_bst(SurfaceFile)`; `TessFile = file_fullpath(SurfaceFile)`.
2. **Cache return:** if not `ForceRecompute` and all five facet fields present,
   derive and return `{U,V,N}` without recompute.
3. Guards: require nxr-compute (`bst_plugin('Install','nxr-compute')`); require
   `Reg.Sphere.Vertices` for the trivial gauge; require a `Structures` atlas with
   L/R labels.
4. `[rH, lH, isConn] = tess_hemisplit(TessMat)`; error if `isConn`.
5. **Per hemisphere** (L then R):
   - Build local submesh `Vloc/Floc` (remap face indices to local).
   - FS-pole singularities: `sph = Reg.Sphere.Vertices(vH,:)`; north = `argmax z`,
     south = `argmin z` (local indices).
   - `h = nxr_compute('create', Vloc, Floc)`;
     `S = nxr_compute('facets', h, 'trivial', struct('singVerts',[iN;iS],'singValues',[1;1]))`;
     `nxr_compute('destroy', h)`.
   - Attach `GlobalVertices = vH`, `GlobalFaces = find(fMask)`, `Hemisphere`,
     `Provenance` to each of the five groups.
6. **Store** five top-level `1×2` per-hemisphere struct arrays
   `TessMat.{Topology, Embedded, Intrinsic, Extrinsic, Gauge}` ((1)=L, (2)=R);
   `bst_history` + `bst_save(TessFile, TessMat, 'v7')` (unless `NoSave`).
7. **Return** full-mesh `{U,V,N}`: per hemisphere
   `c = Embedded.vertex.grid .* Gauge.vertex.rotation`, scatter
   `U(GlobalVertices,:) = real(c)`, `V = imag(c)`, then `N = cross(U,V,2)`.

### Storage layout

Five top-level fields, each `1×2`, each element = the verbatim nxr facet struct
for that hemisphere + `GlobalVertices` / `GlobalFaces` / `Hemisphere` /
`Provenance`. `Provenance`: `Backend='nxr'`, `NxrVersion`, `Gauge`, `ComputeDate`.

### Phase 2 tests (`dev/tests/test_tess_frame.m`, real 20484-vertex cortex)

- `{U,V,N}` are `[nV×3]`; orthonormal; `N == U×V`; every row filled (no gaps).
- The five fields are stored, each `1×2`; `GlobalVertices`/`GlobalFaces`
  disjoint-and-complete partitions of the mesh; `Hemisphere` is `L`/`R`.
- `{U,V,N}` matches `Embedded.vertex.grid ⊙ Gauge.vertex.rotation` per hemisphere.
- Trivial face domain errors (`tess_frame:faceTrivialDeferred`).
- Cache-return path returns without recompute; `ForceRecompute` recomputes.
- Save isolation: back up the `.mat`, restore via `onCleanup`.

**Acceptance:** all Phase 2 tests pass against the installed (Phase-1-rebuilt)
MEX.

---

## Edge cases / risks

- **Face-domain trivial gauge:** `Gauge.face.rotation` is empty/deferred in nxr →
  `tess_frame(...,'Domain','face')` errors clearly (unchanged behavior).
- **Intrinsic-Delaunay caveat:** facet *data* slices the embedded-sourced light
  geometry (per nxr CLAUDE.md), identical to the flat command — no new behavior.
- **Connected hemispheres:** hard error (nxr requires a single genus-0
  component per call); covered by guard.
- **MEX/version drift:** Phase 2 captures `nxr_compute('version')` in
  `Provenance`; optional version bump in Phase 1 lets Phase 2 assert a minimum.

## Out of scope (explicit)

- Deleting `tess_store_perhemi` and fixing/removing the four
  `tess_topology`/`tess_geometry`/`tess_gauge`/`tess_bundle` writers — separate
  cleanup.
- A leadfield→cortical transform accessor (uses the stored frame; later work).
- `.operators` in the `facets` command (available via `operators`).
- A synthesized resolved-frame `Embedded` field (the `Embedded` group is the nxr
  facet: `position`/`normal`/`grid`).
