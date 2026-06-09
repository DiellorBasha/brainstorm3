# Surface-File Bundle Storage — Design

**Date:** 2026-06-09
**Status:** Draft for review
**Branch:** `geom/surface-bundle`
**Supersedes:** `2026-06-08-surface-operators-coordinates-design.md` (the `Operators`/`Coordinates` field approach; that branch is unmerged and is replaced by this one).

## Purpose

Make the Brainstorm surface file a **direct projection of the nxr-compute
coordinate-system bundle**. `nxr_compute('bundle', h, gaugeType, opts)` returns
three element-grouped structures — `Topology`, `Geometry`, `Gauge` — assembled
from geometry-central. This design stores them verbatim as three new canonical
fields:

- `TessMat.Topology`
- `TessMat.Geometry`
- `TessMat.Gauge`

These replace the monolithic `TessMat.nxr` bundle (written by
`tess_nxr_populate.m`) and the just-superseded `Operators`/`Coordinates` fields
with a single, schema-versioned, nxr-native store. The surface file becomes the
sole geometric authority; the head model keeps its Cartesian `Gain` + a
`SurfaceFile` reference, and any consumer (frame re-coordinatization, eigenmode
leadfield) reads these fields.

nxr is canonical. Brainstorm's own `VertNormals`/`VertConn`/face normals are
**not** reconciled into these fields — geometry-central's intrinsic frames (the
flattened one-ring vertex tangent space, the per-triangle face tangent space)
are the truth. The only shared quantity with Brainstorm is the **vertex
positions** (Brainstorm SCS), passed into `nxr_compute('create', V, F)`, which
is what embeds the nxr frames into the same Cartesian space as the leadfield.

## Decisions (settled)

- **Three fields mirror the bundle.** `TessMat.{Topology,Geometry,Gauge}` store
  the nxr `bundle` structs as-is (1-based indices, complex `grid`, etc.).
- **Per-hemisphere `1×2` struct arrays.** nxr never integrates the two
  hemispheres, so each field is a `1×2` struct array (`(1)`=left, `(2)`=right),
  each element the verbatim nxr submesh bundle in local indexing plus
  `GlobalVertices`/`GlobalFaces`/`Hemisphere` scatter maps. Internal
  connectivity stays local; no halfedge re-indexing across hemispheres.
- **Light by default; heavy `.operators` opt-in.** The default write stores the
  *light* bundle (frames, areas, curvature, rotations, singularities, halfedge
  connectivity). The heavy `.operators` sub-structs (DEC `d0/d1/hodge*`, cotan
  `L`, mass, complex connection `K`, `covariantLaplacian`) are written **only**
  when explicitly requested, because they dominate file size and are
  deterministic from `Topology`+`Geometry` (rehydratable on demand).
- **One writer.** `tess_bundle.m` replaces `tess_operators.m` and
  `tess_coordinates.m`. It calls `nxr_compute('bundle', ...)` and stores the
  result.
- **The `{U,V,N}` frame is a derived view, not a stored field.** It is
  `Geometry.<domain>.grid ⊙ Gauge.<domain>.rotation`, embedded. A small helper
  (`tess_frame.m`) returns it on demand; nothing persists it.
- **nxr is canonical; Brainstorm geometry is ignored** (see Purpose). No
  alignment-to-Brainstorm-normal assertion.
- **Additive / non-breaking.** Legacy `TessMat.nxr`, `TessMat.TangentFrame`,
  `tess_nxr_populate.m`, `tess_tangents.m`, and `bst_face_leadfield.m` are
  untouched. Brief on-disk duplication with the legacy `nxr` bundle is accepted.
- **`tess_tangents.m` stays deprecated** (header + soft warning already added)
  and functional this step; migrating its callers is a future step.

## Scope

**In scope**
- `tess_bundle.m` → `TessMat.{Topology,Geometry,Gauge}` (light default, heavy opt-in).
- `tess_frame.m` → derived `{U,V,N}` view from `Geometry.grid` ⊙ `Gauge.rotation`.
- Property-based validation of both, run on a real cortical surface.

**Out of scope (future steps / separate specs)**
- The vertex connection-eigenmode leadfield (`Gauge.Eigen`, `L̃`) — separate
  spec, consumes this one.
- Face-native source space (the DEC dual / connection-Laplacian face domain).
- Migrating consumers off `TessMat.nxr` / `TessMat.TangentFrame`; retiring
  `tess_nxr_populate.m` / `tess_tangents.m`.
- `MeshData` and `Gauge.face.rotation` (deferred on the nxr side).

## `tess_bundle.m` → `TessMat.{Topology,Geometry,Gauge}`

Self-contained compute+store+return, mirroring `tess_tangents.m`'s shape.

```
B = tess_bundle(SurfaceFile)                                   % light, trivial gauge, save
B = tess_bundle(SurfaceFile, 'Gauge', 'trivial')              % gauge type
B = tess_bundle(SurfaceFile, 'Operators', 1)                  % include heavy .operators
B = tess_bundle(SurfaceFile, 'NoSave', 1, 'ForceRecompute', 1)
```

**Options**
- `Gauge` — `'euclidean' | 'levi-civita' | 'trivial'` (default `'trivial'`,
  FreeSurfer-pole singularities, the gauge the eigenmode leadfield needs).
- `Operators` — `false` (default, light) | `true` (attach `.operators` to each
  of the three structs via `opts.operators=true`).
- `Coupling` — `'ambient' | 'product'` (default `'ambient'`), only consulted
  when `Operators` is true (controls `Gauge.operators.covariantLaplacian`).
- `Mass` — `'lumped' | 'galerkin'` (default `'lumped'`), only when `Operators`
  is true.
- `NoSave`, `ForceRecompute` — as in `tess_operators.m`/`tess_coordinates.m`.

**Returned / stored value — per-hemisphere struct arrays.** Because each
hemisphere is bundled independently, each of the three fields is a **`1×2`
struct array** indexed by hemisphere (`(1)`=left, `(2)`=right):

- `TessMat.Topology`  `1×2` — `TessMat.Topology(h)` is hemisphere `h`'s nxr Topology.
- `TessMat.Geometry`  `1×2` — `TessMat.Geometry(h).vertex.grid`, etc.
- `TessMat.Gauge`     `1×2` — `TessMat.Gauge(h).vertex.rotation`, etc.

Each element is the **verbatim nxr submesh bundle in local indexing**, plus two
maps the writer adds so consumers can scatter to global `TessMat.Vertices`/
`.Faces` order:
- `GlobalVertices` `[nVh×1]` — local→global vertex indices (`= lH` / `rH`).
- `GlobalFaces`    `[nFh×1]` — local→global face indices (the hemisphere face mask).
- `Hemisphere` — `'L'` / `'R'`.

These maps live on each of the three structs (small, cheap). All connectivity
cross-references inside a hemisphere's `Topology`/`.operators` stay **local and
self-consistent** — never remapped — so the fragile halfedge re-indexing of a
stitched mesh is avoided entirely. Global assembly happens only where it is
trivial (scatter of per-element vertex/face values) or where the math is
naturally block-diagonal (operators, eigenmodes).

**Field schema (mirrors the nxr bundle; representative fields):**
- `Topology` — 1-based halfedge struct-of-arrays (`vertex/edge/face/corner/
  halfedge` groups; `0` = none/boundary). `+.operators` (heavy): `laplacian`
  (graph `d0ᵀd0`), `dec.{d0,d1}`.
- `Geometry` — `totalArea`; `vertex.{dualAreas, angleSums, curvature(2-RoSy q +
  meanCurvature), grid(V×3 complex, c=e1+i·e2)}`; `edge.{lengths, cotanWeights,
  dihedralAngles}`; `face.{areas, centroids, grid(F×3 complex)}`;
  `halfedge.{cotanWeights, vectorsInVertex, vectorsInFace, transportAlong,
  transportAcross}`; `corner.{angles, scaledAngles}`. `+.operators` (heavy):
  `laplacian` (cotan), `mass.{lumped,galerkin}`, `hodge.{h0,h1,h2,h1inv}`.
- `Gauge` — `type`; `vertex.rotation` (V×1 complex `exp(iθ)` vs Levi-Civita
  grid; trivial only); `singularity.{vertices(1-based), indices(Σ=χ), source}`.
  `+.operators` (heavy): `laplacian` (complex `V×V` connection `K` in the
  current gauge), `covariantLaplacian` (`3N×3N` real, `coupling` ambient/product).
- Each struct carries its nxr `schemaVersion`; the writer adds provenance
  (`Backend='nxr'`, `nxr_version`, `ComputeDate`, gauge/coupling/mass options used).

**Computation (per hemisphere)**
```matlab
bst_plugin('Install','nxr-compute');                 % nxr-required, no fallback
[rH, lH, isConn] = tess_hemisplit(TessMat);          % import-label split
assert(~isConn, 'tess_bundle:connectedHemispheres'); % each must be its own component
hemis = {lH(:), rH(:)}; tags = {'L','R'};
for hh = 1:2
    vH   = hemis{hh};                                % global vertex indices
    fMask = all(ismember(TessMat.Faces, vH), 2);     % faces fully inside this hemi
    map  = zeros(size(TessMat.Vertices,1),1); map(vH) = 1:numel(vH);
    Vloc = double(TessMat.Vertices(vH,:));
    Floc = map(TessMat.Faces(fMask,:));              % local 1-based faces
    [iN,iS] = local_poles(TessMat.Reg.Sphere.Vertices(vH,:));  % N/S sphere poles (local)
    h = nxr_compute('create', Vloc, Floc);
    opts = struct('singVerts',[iN;iS], 'singValues',[1;1]);    % trivial-gauge poles
    if Operators, opts.operators = true; opts.coupling = Coupling; opts.mass = Mass; end
    Bh = nxr_compute('bundle', h, Gauge, opts);
    nxr_compute('destroy', h);
    % attach scatter maps, accumulate into the 1x2 struct arrays
    for f = {'Topology','Geometry','Gauge'}
        s = Bh.(f{1});
        s.GlobalVertices = vH;  s.GlobalFaces = find(fMask);  s.Hemisphere = tags{hh};
        B.(f{1})(hh) = s;
    end
end
```
(`local_poles` = max/min sphere-`z`; the singularity/Gauss–Bonnet handling
mirrors `tess_coordinates.m`'s `solve_hemisphere`.) For `euclidean`/`levi-civita`
gauges, `singVerts`/`singValues` are omitted.
Singularity placement for the `trivial` gauge reuses the hemisphere-split +
FreeSurfer-pole logic currently in `tess_coordinates.m`/`tess_tangents.m`
(import-label hemisphere split, north/south sphere poles, Gauss–Bonnet check).

**Per-hemisphere by construction (settled).** nxr does **not** integrate the two
hemispheres — every `bundle` call operates on a single connected genus-0
submesh. So `tess_bundle` splits the cortex (`tess_hemisplit` → `lH`, `rH`,
asserting the two are not connected), and for each hemisphere: builds the
submesh in **local** indexing, creates an nxr handle, calls `bundle` with that
hemisphere's two FreeSurfer poles as the trivial-gauge singularities, and
destroys the handle. This is also the *correct* representation downstream: the
connection-Laplacian eigenbasis (Spec 2) is per connected component — there are
no cross-hemisphere modes — so the natural store is per hemisphere, not a
stitched single mesh.

**Save** load full file → set the three fields → `bst_history('add', …)` →
`bst_save(..., 'v7')` (the `out_tess_eigenmodes.m` / `tess_operators.m` pattern).

**Caveat (record in code):** `Geometry.operators.mass.lumped` is the lumped
**barycentric** mass (`vertexLumpedMassMatrix`), not Voronoi/circumcentric —
per nxr's own `CLAUDE.md`. Document on the field; don't trust a "Voronoi" label.

## `tess_frame.m` → derived `{U,V,N}`

```
[U,V,N] = tess_frame(SurfaceFile)                 % default 'vertex' domain
[U,V,N] = tess_frame(SurfaceFile, 'Domain','face')
```
Iterates the two hemispheres, and for each reads
`TessMat.Geometry(h).<domain>.grid` and `TessMat.Gauge(h).<domain>.rotation`,
applies `c' = c .* rotation` (gauge rotation of the grid), then **scatters** into
full-mesh rows via `GlobalVertices` (vertex domain) or `GlobalFaces` (face
domain): `U(global,:) = real(c')`, `V(global,:) = imag(c')`, `N = cross(U,V)`
(`nVert×3` / `nFace×3` full-mesh). Pure function over the stored bundle —
computes nothing new, persists nothing. If `Gauge(h).<domain>.rotation` is absent
(euclidean/levi-civita: grid is already the gauge), the grid is used directly
(`rotation ≡ 1`). The full-mesh `{U,V,N}` is the per-source `SO(3)` transform a
consumer applies to re-express a Cartesian `Gain` block in intrinsic
coordinates.

**Domain default is `vertex`** because that is the source space we are
targeting and the only domain for which the **trivial** gauge currently ships a
rotation: nxr exposes `Gauge.vertex.rotation` but `Gauge.face.rotation` is
**deferred** (see nxr `CLAUDE.md`). So `tess_frame` supports: vertex domain for
all gauges; face domain only for `euclidean`/`levi-civita` (rotation ≡ 1). A
face-domain **trivial** frame must wait for `Gauge.face.rotation` and is part of
the future face-native work — `tess_frame` errors clearly if asked for it.

## Backward compatibility

- All three fields are **new and additive**; no `db_template` change needed for
  load. Legacy readers (`tess_nxr_populate`, `tess_tangents`,
  `bst_face_leadfield`, scalar `Eigenmodes` consumers) are untouched.
- The morning's `Operators`/`Coordinates` fields are **not** introduced (that
  branch is superseded, unmerged). Any surface already carrying them is
  unaffected — `tess_bundle` writes different field names.

## Validation (property tests, real surface, MATLAB MCP)

Surface: the loaded protocol's `tess_cortex_*_low` (≈20484 verts; needs
nxr-compute + FreeSurfer registration sphere). Follow the demo/`tess_operators`
TDD pattern (`matlab.unittest`, `local_find_cortex(20484)`, backup/restore file
isolation via `onCleanup`).

**`tess_bundle` (light):**
- Each of `Topology`, `Geometry`, `Gauge` is a `1×2` struct array; each element
  carries `schemaVersion` and `GlobalVertices`/`GlobalFaces`/`Hemisphere`.
- The two hemispheres' `GlobalVertices` are **disjoint** and their **union is
  `1:nVert`** (every vertex assigned once); likewise `GlobalFaces` over `1:nFace`.
- Light write has **no** `.operators` on any struct.
- `Geometry(h).vertex.grid` is `nVh×3` complex; `.face.grid` is `nFh×3` complex.
  Per element the embedded frame is orthonormal: with `e1=real(c)`,
  `e2=imag(c)`, `‖e1‖=‖e2‖=1`, `e1·e2≈0` (the `SO(3)` precondition the ambient
  operator relies on).
- `Gauge(h).type` matches the requested gauge; trivial carries `singularity`
  with `Σ indices == χ == 2` for the hemisphere (Gauss–Bonnet), `source` recorded.
- Stored on disk (real save path) and reloads via `in_tess_bst` as a `1×2` array.

**`tess_bundle` (heavy, `Operators=1`):** (per hemisphere `h`)
- `Topology(h).operators.dec.{d0,d1}` present; `d1*d0 == 0` (`d∘d=0`).
- `Geometry(h).operators.{laplacian,mass,hodge}` present; cotan
  `K=d0'*hodge.h1*d0` symmetric, near-zero row sums.
- `Gauge(h).operators.laplacian` is complex `nVh×nVh`, Hermitian (`‖K−K'‖≈0`).
- `Gauge(h).operators.covariantLaplacian` is `3nVh×3nVh` real, symmetric;
  default coupling `ambient`; world-form sanity
  `blockdiag(F)·L3·blockdiag(F)' ≈ kron(I3, cotanL)` (the defining ambient identity).

**`tess_frame`:** (vertex domain, full-mesh output)
- `[U,V,N]` are `nVert×3`, orthonormal per row, `N==cross(U,V)`, no zero rows
  (every vertex filled by exactly one hemisphere's scatter).
- For each hemisphere, the scattered rows equal `Geometry(h).vertex.grid`
  decomposition rotated by `Gauge(h).vertex.rotation` (consistency with the
  stored bundle).
- Face domain under `trivial` gauge errors clearly (deferred
  `Gauge.face.rotation`); face domain under `levi-civita` returns the grid.

## Resolved decisions

1. **Field names:** `TessMat.Topology`, `TessMat.Geometry`, `TessMat.Gauge`
   (top-level, mirror the nxr `bundle`).
2. **Persistence:** light by default; heavy `.operators` opt-in via
   `'Operators',1`. Frame `{U,V,N}` is derived (`tess_frame`), never stored.
3. **Writer:** single `tess_bundle.m`; supersedes `tess_operators.m` /
   `tess_coordinates.m`.
4. **Canonical authority:** nxr frames; Brainstorm normals ignored; vertex
   positions (SCS) are the only shared quantity.
5. **Per-hemisphere shape:** each field is a `1×2` struct array (nxr never
   integrates hemispheres; the connection-Laplacian basis is per component).
   Local connectivity is never remapped; `GlobalVertices`/`GlobalFaces` carry
   the scatter to global order.

## Future steps (built on this)

- **Vertex connection-eigenmode leadfield** (next spec): `eigs(K,M)` over
  `Gauge.operators.laplacian` → `TessMat.Gauge.Eigen.{Vectors,Values}`;
  `L̃ = (G·cᵀ)·Vectors`; `ambient covariantLaplacian` as inverse regularizer.
- Face-native source space via the DEC dual (connection-Laplacian face domain).
- Migrate `bst_face_leadfield.m` and other `tess_tangents`/`nxr` consumers onto
  the bundle; retire `tess_tangents.m` and `tess_nxr_populate.m`.
