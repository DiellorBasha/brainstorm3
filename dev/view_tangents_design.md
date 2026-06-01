# view_tangents — Design Spec

**Date:** 2026-06-01
**Branch:** `feature/tess-tangents`
**Status:** Approved (brainstorming) → ready for implementation plan

## Goal

Provide `view_tangents.m`, an interactive 3D viewer that displays the per-face
tangent frame field computed and stored by `tess_tangents` (`TessMat.TangentFrame`),
so the frame can be visually inspected and verified — in particular the global
consistency of the field and the four singularities at the FreeSurfer registration
poles.

## Approach

Standalone GUI function `toolbox/gui/view_tangents.m`, modeled on the existing
`toolbox/gui/view_leadfield_vectors.m`. It opens a normal Brainstorm surface
figure (`view_surface`), `hold on`, and overlays the frame with MATLAB's built-in
`quiver3` — exactly the established Brainstorm vector-field idiom (also used by
`panel_opticalflow.m`). Compute stays in `tess_tangents`; display stays here
(clean separation). No new rendering machinery.

Rejected alternatives:
- A `'Display'` mode inside `tess_tangents.m` — welds GUI onto the compute
  function, breaks the `view_*` naming convention.
- A generic `view_surface_vectors` shared utility — a refactor of working code
  (`view_leadfield_vectors`) for no immediate gain. Possible future consolidation,
  not now.

## Interface

```matlab
hFig = view_tangents(SurfaceFile)                     % default density
hFig = view_tangents(SurfaceFile, 'MaxArrows', 3000)  % override initial arrow count
```

- `SurfaceFile` — Brainstorm cortex surface (relative or full path).
- Option `'MaxArrows'` (default `2000`) — initial number of face frames drawn.
- Returns the figure handle `hFig`.

## Data flow

1. `TessFile = file_fullpath(SurfaceFile)`; `TessMat = in_tess_bst(SurfaceFile)`.
2. **Get the frame.** If `TessMat.TangentFrame` is missing/empty, call
   `tess_tangents(SurfaceFile)` (computes + stores; this is where the
   nxr / Reg.Sphere / hemisphere-label requirements live), then re-load the
   surface so `TangentFrame.Singularities` is available.
3. Validate `TangentFrame.Domain == 'face'`; otherwise
   `error('view_tangents:unsupportedDomain', ...)` (vertex frames are future work).
4. Pull `U`, `V` (`Nf×3`, cast to `double`) and the `Singularities` struct.
5. Display geometry (plain MATLAB — viewer stays plugin-free for stored frames):
   - **Face normals** via `tess_normals(Vertices, Faces)` (reuse; returns
     `[VertNormals, FaceNormals]`).
   - **Face centroids** = `(V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:)) / 3`
     (one-liner; no dedicated helper exists, matching `tess_curvature` / `fem_resect`).
   - **meanEdge** = mean triangle edge length, used as the base arrow size and the
     surface offset distance.

## Geometry & rendering

- Base figure: `view_surface(SurfaceFile, 0.5, [.5 .5 .5], 'NewFigure')`
  (semi-transparent grey cortex), grab the `'Axes3D'` handle, `hold on`.
- **Offset arrow bases off the surface** along the face normal:
  `base = centroid + 0.25*meanEdge * faceNormal` (the `panel_opticalflow` trick),
  so arrows sit on the cortex instead of sinking in.
- **Equal-length arrows:** normalize U, V, N to unit and scale by
  `quiverSize * meanEdge`; call `quiver3` with autoscaling **off** (scale arg `0`)
  so every frame arrow has the same length (a frame field has no magnitude to encode).
- Three `quiver3` sets over the current face subset:
  | Arrow | Color | Tag |
  |-------|-------|-----|
  | U | yellow `[1 1 0]` | `tangentU` |
  | V | yellow `[1 1 0]` (same as U) | `tangentV` |
  | normal | magenta `[1 0 1]` | `tangentN` |

  U and V share a color (the in-plane tangent cross); the normal is distinct.

## Singularity markers

`TangentFrame.Singularities` holds the global pole vertex indices and hemisphere
tags (4 poles — north/south per hemisphere). Draw with `plot3` as large filled
red markers (`MarkerSize` ~10, Tag `tangentSing`) at `Vertices(Singularities.Vertices,:)`.
These are the points the field rotates around — the primary visual check.

## Interaction

Hijack the figure `KeyPressFcn` (same idiom as `view_leadfield_vectors`): unhandled
keys fall through to the saved callback so normal rotation/standard views still work.
Every handled key triggers a single `DrawArrows()` redraw (which first deletes the
`tangent*`-tagged objects).

| Key | Action |
|-----|--------|
| Left / Right | fewer / more arrows (density ÷ / × 1.5, clamped `[~200, nFaces]`) |
| Shift + Up/Down | arrow length × / ÷ 1.2 (`quiverSize`) |
| Ctrl + Up/Down | line width × / ÷ 1.2 (`quiverWidth`) |
| N | toggle normal arrows |
| P | toggle singularity (pole) markers |
| H | help dialog listing the shortcuts |

A bottom `uicontrol` text label shows the legend + live state:
*"U,V (yellow) · normal (magenta) · `<n>`/`<total>` faces shown · H for help"*.

## Decomposition

- `local_arrow_field(Vtx, Fcs, U, V, N, idx, scale, offset)` — **pure** local helper.
  Given the subsample index list and scaling, returns the per-arrow base XYZ and
  arrow UVW for U, V, N. No figure handles inside, so it is unit-testable on its
  own. `DrawArrows()` just feeds its output to `quiver3`.
- Subsample indices are **deterministic**: `idx = unique(round(linspace(1, nF, nArrows)))`
  (uniform stride over face index; reproducible for testing).

## Tree integration

In `toolbox/tree/tree_callbacks.m`, the `'cortex'` popup block (alongside the other
cortex-only items, ~line 1263), single-node only:

```matlab
gui_component('MenuItem', jPopup, [], 'Display tangent basis', IconLoader.ICON_DISPLAY, [], ...
    @(h,ev)bst_call(@view_tangents, filenameRelative));
```

`bst_call` ensures any error surfaces as a clean Brainstorm error dialog rather than
a stack trace.

## Error handling

- Auto-compute errors propagate **with their identifiers**
  (`tess_tangents:noRegSphere`, `:noHemisphereLabels`, `:nxrUnavailable`,
  `:connectedHemispheres`, `:gaussBonnet`) — the viewer never swallows them.
- `TangentFrame.Domain ~= 'face'` → `error('view_tangents:unsupportedDomain', ...)`.
- **Known limitation:** a read-only protocol with no stored frame will fail at the
  `tess_tangents` save step (error propagates). No NoSave display fallback for now
  (YAGNI); revisit if it becomes a real workflow need.

## Testing

`dev/tests/test_view_tangents.m`, run via the MATLAB MCP `evaluate_matlab_code`
(project pattern; function-style tests are invoked directly, not via `runtests`).

1. **Integration / auto-compute + render.** Discover a registered, labeled cortex
   (reuse the discovery approach from `test_tess_tangents.m`), copy it to a temp file
   with **no** stored `TangentFrame`. Then:
   - `hFig = view_tangents(tmpFile, 'MaxArrows', 500)`.
   - Assert auto-compute fired: the temp file now has `TangentFrame` stored.
   - Assert `ishandle(hFig)`.
   - Assert `tangentU` / `tangentV` / `tangentN` quiver objects exist with **equal**
     arrow counts (≤ 500).
   - Assert U and V share a color that **differs** from N.
   - Assert `tangentSing` marker count == `numel(Singularities.Vertices)` (4 for two
     hemispheres).
   - `close(hFig)`.
2. **Unit test of `local_arrow_field`.** Deterministic counts, equal arrow lengths,
   and that bases are offset off the surface along the normal — independent of any
   figure.

## Files

- Create: `toolbox/gui/view_tangents.m`
- Modify: `toolbox/tree/tree_callbacks.m` (cortex popup — one menu item)
- Create: `dev/tests/test_view_tangents.m`
- Reuse (no change): `toolbox/anatomy/tess_tangents.m`, `toolbox/anatomy/tess_normals.m`,
  `toolbox/gui/view_surface.m`
