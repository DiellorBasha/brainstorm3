# view_manifold — Consolidated Manifold Viewer — Design

- **Date:** 2026-06-12
- **Status:** Draft for review
- **Author:** Diellor Basha (with Claude)
- **Related:** `view_tangents.m` (replaced), `tess_frame.m`, `tess_manifold.m`, `tess_tangents.m`
  (deprecated), `view_connection_phase.m`, `tree_callbacks.m`, `db_template('manifoldmat')`

## 1. Goal

Provide a single consolidated viewer, `view_manifold`, for a manifold DB node (`manifold_*.mat`).
It reads the canonical frame data stored in the node and renders it on the parent cortex. The first
view is the **per-vertex tangent-basis + normal frame**; the function is structured so later, more
complex views (curvature RoSy, parallel transport, gauge) are added as new view modes.

This replaces the deprecated `view_tangents` (which recomputed a per-face frame via `tess_tangents`)
and makes the manifold node the source of frame information for the GUI.

## 2. Background and data

`tess_manifold` stores a `manifold_*.mat` node (`db_template('manifoldmat')`) with 1×2
per-hemisphere structs: `Topology, Embedded, Intrinsic, Extrinsic, Gauge`, plus `ParentSurface`.

The canonical per-vertex frame is derived exactly as in `tess_frame.m`'s `local_derive_frame`:

```
grid = Embedded(hh).vertex.grid     % [nVh x 3] COMPLEX  (encodes U + iV)
rot  = Gauge(hh).vertex.rotation    % [nVh x 1] complex unit (gauge rotation)
cRot = grid .* rot
U = real(cRot);  V = imag(cRot);  N = cross(U, V)
```

Verified live on the canonical cortex: `U,V` are unit, mutually orthogonal, both ⟂ `N`, and
`cross(U,V)` reproduces `Embedded.vertex.normal` to ~1e-16.

Other node data available for **future** views (out of scope here): `Extrinsic.vertex.principalDir`
(unit curvature direction), `meanCurvature`, `Intrinsic.halfedge.transport*`, `Gauge.singularity`.

**Per-face limitation:** the per-face frame needs `Gauge.face.rotation`, which nxr leaves
**empty/deferred for the trivial gauge** (documented in `tess_frame`). So `view_manifold` displays
the **per-vertex** frame only. Per-face is deferred until nxr populates `Gauge.face.rotation`.

## 3. Scope

**In scope (all unblocked):**
1. New `toolbox/gui/view_manifold.m` — per-vertex frame view from the manifold node.
2. Delete `toolbox/gui/view_tangents.m` and `dev/tests/test_view_tangents.m`.
3. Manifold node menu → single **"View manifold"** item → `view_manifold`; remove "Display tangent
   basis", the interim `ManifoldViewTangents_Callback`, and the field-dump `ManifoldView_Callback`.
4. Repoint `view_connection_phase` off `tess_tangents` → the manifold/`tess_frame` per-vertex frame.
5. Mark `tess_tangents.m` `@deprecated` in its header (keep the function).

**Out of scope (deferred / blocked):**
- Physically removing `tess_tangents.m`: blocked — `bst_wavefront_track` and `tess_nxr_populate`
  consume **per-face** frames the manifold cannot yet supply (nxr `Gauge.face.rotation` deferral).
  A later task removes it once per-face frames exist (or those callers are reworked).
- Per-face frame view, curvature/transport/gauge views — future `view_manifold` view modes.
- The deeper `view_connection_phase` behavior rework (a separate later step).

## 4. Components

### 4.1 `toolbox/gui/view_manifold.m` (new)

Public: `hFig = view_manifold(ManifoldFile, 'View','frame', 'MaxArrows', N)`.

- **Macro dispatch:** `view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)` routes to the pure
  subfunction for headless tests; otherwise routes to the GUI entry.
- **Load + validate:** `file_exist` guard → `bst_error`; `load` the node; require `ParentSurface`,
  `Embedded` (1×2), `Gauge` (1×2).
- **View dispatch:** `switch lower(View)` — `'frame'` implemented; unknown → `bst_error`.
- **`DeriveVertexFrame(Embedded, Gauge, nVert)`** (pure) → struct with `P` (anchors, `[nVert x 3]`),
  `U, V, N` (`[nVert x 3]`, zeros off-support), and `Sing` (global singularity vertex ids). Per
  hemisphere: `grid=Embedded(hh).vertex.grid`, `rot=Gauge(hh).vertex.rotation`,
  `cRot=grid.*rot`, `U=real(cRot)`, `V=imag(cRot)`, `N=cross(U,V)`, scattered to global indices via
  `Embedded(hh).GlobalVertices`; `P=Embedded(hh).vertex.position`;
  `Sing=GlobalVertices(Gauge(hh).singularity.vertices)`. Errors on shape mismatch / missing fields.
- **GUI (`ViewFigure`)**: `view_surface(Surface, 0, [.5 .5 .5], 'NewFigure', 0)`; opaque, no smooth,
  wireframe edges on (`panel_surface('SetSurfaceSmooth'/'SetSurfaceEdges')`); resolve `Axes3D` and
  **`hold(hAxes,'on')` before any `quiver3`** (axes-reset trap); hijack `KeyPressFcn`; status label.
- **`DrawArrows` (nested):** delete prior `Tag^='tangent'`; subsample faces→here vertices via
  `ArrowSubsample(nVert, nArrows)`; build glyphs via `ArrowField(P, N, U, V, idx, len, offset)`
  (global glyph length `quiverSize * 0.5 * meanEdge`, base offset along `N`); headless `quiver3`
  for `U` and `V` (shared color, `ShowArrowHead off`, `Tag 'tangentU'/'tangentV'`), `N` (toggle,
  `Tag 'tangentN'`), singularity lollipops (toggle, radial lift, `Tag 'tangentSing*'`); legend +
  status label.
- **Keyboard** (same idiom as `view_tangents`): Left/Right density; Shift+Up/Down glyph length;
  Ctrl+Up/Down width; `N` toggle normals (off by default); `P` toggle singularities; `H` help;
  otherwise → stashed callback.
- **Pure helpers** `ArrowSubsample`, `ArrowField` copied from `view_tangents` (the established local
  idiom), anchored at vertices.

### 4.2 `toolbox/tree/tree_callbacks.m`

- Manifold popup (`case 'manifold'`): one item `'View manifold'` → `view_manifold(filenameFull)`.
- Remove the `'Display tangent basis'` item and the `ManifoldViewTangents_Callback` function
  (added in the previous relocation commit) and the `ManifoldView_Callback` field-dump stub.
- `Delete` item unchanged.

### 4.3 `toolbox/gui/view_connection_phase.m`

Replace (≈ lines 77–79):
```matlab
[Uf, ~]  = tess_tangents(SurfaceFile, 'NoSave', 1);
[Uv, Vv] = bst_tangent_face2vertex(Fcs, Uf, Nv);
FsFrame  = struct('e1', Uv, 'e2', Vv);
```
with:
```matlab
[Uv, Vv] = tess_frame(SurfaceFile);          % manifold gauge frame (vertex domain)
FsFrame  = struct('e1', Uv, 'e2', Vv);
```
`tess_frame` returns the full-mesh per-vertex `[U,V,N]` from the manifold/facets bundle (computing
and storing on the surface `TessMat` if absent). This drops both `tess_tangents` and
`bst_tangent_face2vertex` from this path. Frame semantics shift to the manifold gauge frame
(intended); the deeper connection-phase rework remains a later step.

### 4.4 `toolbox/anatomy/tess_tangents.m`

Add an `@deprecated` note to the header: superseded by the manifold frame (`view_manifold` /
`tess_frame`); retained only for the per-face callers (`bst_wavefront_track`, `tess_nxr_populate`)
until per-face manifold frames are available. No behavior change.

### 4.5 Deletions

- `toolbox/gui/view_tangents.m`
- `dev/tests/test_view_tangents.m`

## 5. Data flow

```
view_manifold(ManifoldFile)
  → load node; Surface = ManifoldMat.ParentSurface; TessMat = in_tess_bst(Surface); nVert
  → DeriveVertexFrame(Embedded, Gauge, nVert) → {P, U, V, N, Sing}
  → view_surface(Surface, opaque) + wireframe + hold(Axes3D)
  → DrawArrows(): quiver3 U,V (+N toggle) at P; singularity lollipops
  → KeyPressFcn: density / length / width / N / P / H
```

## 6. Edge cases & error handling

- Missing file / no `ParentSurface` / `Embedded` not 1×2 / `Gauge` not 1×2 → `bst_error` + return.
- `Embedded(hh).vertex.grid` row count ≠ `numel(GlobalVertices(hh))` → `bst_error`.
- Empty `Gauge(hh).singularity` → no lollipops (non-fatal).
- Near-zero frame vectors (degenerate) → ε-guarded unit normalization in `ArrowField` (no NaN).
- Unknown `'View'` value → `bst_error`.

## 7. Testing

### 7.1 Headless — `dev/tests/test_manifold_frame.m`

Via `view_manifold('DeriveVertexFrame', Embedded, Gauge, nVert)` (no figure):
- Synthetic 2-hemi `Embedded`/`Gauge`: `grid` = `(e1 + i e2)` for chosen orthonormal `e1,e2`, `rot=1`.
  Assert `U=e1`, `V=e2`, `N=cross(e1,e2)`, all unit, `U⟂V`, scattered to the right global vertices,
  off-support rows zero.
- A non-trivial `rot=exp(iθ)` rotates `U,V` in-plane by `θ` (frame stays orthonormal).
- `Gauge.singularity.vertices` map to the correct global ids in `Sing`.
- Shape-mismatch / missing-field inputs error.

### 7.2 Live — `dev/tests/test_view_manifold.m`

Requires Brainstorm + the registered manifold node:
- `view_manifold(ManifoldFile)` returns a valid figure; `Axes3D` and the cortex patch survive the
  glyph draw; `tangentU`/`tangentV` quivers exist with non-empty `UData`.
- `N` key adds a `tangentN` quiver; `P` toggles `tangentSing`; density key changes glyph count.
- Close cleanly.

### 7.3 Repoint guard

- `grep -n "tess_tangents" toolbox/gui/view_connection_phase.m` → empty.
- `tess_frame(SurfaceFile)` returns a `[nVert x 3]` `U` on the canonical cortex (sanity).

## 8. Future work (out of scope)

- Per-face frame view (needs nxr `Gauge.face.rotation`).
- Curvature-RoSy, parallel-transport, gauge-rotation view modes in `view_manifold`.
- Physical removal of `tess_tangents` after the per-face callers are repointed.
- The deeper `view_connection_phase` behavior rework.
