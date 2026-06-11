# Source-Vector Quiver Overlay for the Cortical Activations View — Design

- **Date:** 2026-06-11
- **Status:** Draft for review
- **Author:** Diellor Basha (with Claude)
- **Related:** `bst_dirac.m`, `bst_inverse_dirac.m`, `figure_3d.m`, `panel_surface.m`, `view_leadfield_vectors.m`, `view_connection_phase.m`

## 1. Goal

Add a **quiver (arrow) overlay** to Brainstorm's existing "Cortical Activations → Display on cortex" 3D view so the user can see the **direction** of an unconstrained (3-component) source field — and **time-step through it** (e.g. watch the alpha-cycle field rotate/evolve) — exactly as the colormap display already does, with the arrows updating on every time-cursor change.

The motivating use case is the Dirac-eigenmode vector inverse field (the posterior alpha burst in the `data_block001_band` test segment), but the overlay reads the generic unconstrained 3-vector of whatever source is displayed.

## 2. Background and conceptual framing

The recovered source is an **unconstrained vector field**: three Cartesian components per cortical vertex, in the SCS/head frame shared by the displayed surface. A 3D vector anchored at a surface vertex is a **section of the ambient bundle restricted to the surface — the pullback bundle `f*Tℝ³`** — which over each point is a full copy of ℝ³, *not* the 2D tangent plane. Display is therefore purely **ambient and Cartesian**: the arrow components are the source-vector components, drawn directly.

The **Dirac operator and its eigenmodes are out of scope for the renderer.** They enter only when we later do differential analysis / filtering of the field, which must respect the cortical surface geometry. This spec is the display layer that makes the field visible; the analysis layer is a separate, future spec.

## 3. Non-goals

- **No changes to the inverse.** `bst_inverse_dirac` already returns the full vector kernel (`ImagingKernel` is `[3·nVert × nCh]`) and the mode-coefficient kernel (`ImagingKernelMode`). Nothing about the solver changes.
- **No Dirac/mode reconstruction at display time.** Arrows are the ambient 3-vector, not a band-limited reconstruction from selected modes (a possible future extension).
- **No differential analysis** (divergence, curl, flow, dispersion) in this spec.
- **No new top-level viewer window.** The overlay lives inside the existing 3D figure.

## 4. Existing code touchpoints

From a survey of the rendering pipeline:

- **Time-cursor → color redraw:** `figure_3d.m:ColormapChangedCallback` → `panel_surface.m:UpdateSurfaceColormap` → `figure_3d.m:UpdateSurfaceColor` (sets `FaceVertexCData` on the cortex patch). The overlay update hooks into this same path.
- **Un-collapsed 3-vector availability:** `bst_memory.m:GetResultsValues` loads the full 3-component values (`ImageGridAmp(:,iTime)` or kernel×data) *before* `bst_source_orient` collapses them to an RMS scalar. The overlay needs the values at the pre-collapse point.
- **Quiver precedent:** `view_leadfield_vectors.m:DrawArrows` uses `quiver3` anchored at grid locations; `view_connection_phase.m` decimates a cortex field for arrow display (`step = ceil(nSupra / MaxArrows)`, arrow-key density scaling). These are the models for `PlotSourceVectors`.
- **Sparse-overlay precedent:** `panel_surface.m:PlotGrid` is the structural analog for a new `PlotSourceVectors(hFig, iTess)` helper.
- **Toggle location:** `figure_3d.m:DisplayFigurePopup` (right-click menu) — add a `CheckBoxMenuItem`.

## 5. Design

### 5.1 Components and interfaces

- **`PlotSourceVectors(hFig, iTess)`** *(new; in `panel_surface.m` or `figure_3d.m`, analogous to `PlotGrid`)* — creates the `quiver3` handle for surface `iTess`, anchored at the (optionally decimated) vertex set, offset slightly along vertex normals. Stores the handle and the chosen vertex index set on `TessInfo(iTess)`. Idempotent: removes any prior quiver before recreating.
- **`UpdateSourceVectors(hFig, iTess)`** *(new)* — fetches the current-time 3-vector for the stored vertex set, normalizes to unit length (with ε-guard), and updates the existing quiver handle's `UData/VData/WData` (anchors fixed → no recreation). No-op if the overlay is off or the source is not unconstrained.
- **Un-oriented fetch** — obtain the per-vertex 3-vector at the current time without the RMS collapse. Preferred: a small option/path through `GetResultsValues` (e.g. an `ApplyOrient=0` fetch) returning `[3·nVert × 1]` for the current time; fall back to multiplying the stored `ImagingKernel` rows by the current data column if a clean fetch is not exposed.
- **Toggle** — `CheckBoxMenuItem` "Show source vectors (quiver)" in `DisplayFigurePopup`; calls a setter that flips `TessInfo(iTess).ShowSourceVectors`, then `PlotSourceVectors` (on) / deletes the handle (off), and a redraw.
- **State on `TessInfo(iTess)`:** `ShowSourceVectors` (bool), `SourceVectorHandle` (graphics handle), `SourceVectorIdx` (vertex indices drawn), `SourceVectorScale` (length in mm), `SourceVectorMaxArrows` / threshold (the "threshold the quiver number" control).

### 5.2 Glyph semantics (defaults)

| Property | Default | Notes |
|---|---|---|
| Encodes | **Direction only** | Amplitude is read from the existing colormap on the surface. |
| Length | **Unit-normalized**, scaled by a global `SourceVectorScale` (mm) | ε-guard: vertices with magnitude < ε are not normalized (skip / zero-length) to avoid drawing random directions from numerical noise. |
| Color | Single neutral fixed color (configurable) | Not magnitude — that would duplicate the colormap. |
| Anchor | Vertex position + small offset along vertex normal | Avoids depth-buffer occlusion on the opaque cortex. |
| Surface | **Opaque by default** | Transparency left to the existing Data/Surface sliders. |
| Vertex set | **Entire field by default** | Control to threshold/decimate the arrow count (by amplitude threshold and/or `MaxArrows` decimation). |
| Frame | Ambient SCS Cartesian | Components → `U,V,W` directly; no projection, no transform. |

### 5.3 Data flow (per time step)

```
global time cursor changes
  → figure_3d:ColormapChangedCallback
    → panel_surface:UpdateSurfaceColormap        (existing: colormap = RMS norm)
      → figure_3d:UpdateSurfaceColor             (existing: FaceVertexCData)
      → UpdateSourceVectors(hFig, iTess)         (NEW: ~1 added call)
          • fetch 3-vector at iTime for SourceVectorIdx (un-oriented)
          • normalize (ε-guard) → unit directions
          • set quiver UData/VData/WData (× SourceVectorScale)
```

### 5.4 Default arrow count / performance

"Entire field" means an arrow per vertex by default (~20k on the standard cortex). Because `UpdateSourceVectors` only rewrites the quiver's `U/V/W` data each step (anchors fixed), per-step cost is low. If interactive time-stepping is sluggish at full density, the `MaxArrows`/threshold control decimates; we will measure on the alpha block and set a sensible default cap only if needed (documented, not silent).

## 6. Edge cases and error handling

- **Constrained / scalar source (`nComponents == 1`):** no direction to show → the menu item is disabled (or the overlay silently no-ops) for that result.
- **Kernel vs full-matrix results:** support both (kernel×data on the fly, or `ImageGridAmp` slice).
- **Near-zero magnitude:** ε-guard prevents NaN/normalization noise; such vertices draw no arrow.
- **Decimation reproducibility:** the drawn vertex index set is fixed when the overlay is created (stored), so arrows don't "swim" between frames.
- **Multiple surfaces / overlays:** state is per-`TessInfo`, so a figure with several layers behaves predictably.

## 7. Validation (manual, GUI)

1. Load the Dirac unconstrained source for `data_block001_band` (posterior alpha), open Display on cortex.
2. Toggle "Show source vectors (quiver)" → arrows appear over the active posterior region; amplitude still legible via the colormap.
3. Time-step across the burst (≈21.5–24 s): arrows update smoothly and the field direction evolves cycle-by-cycle over the right parieto-occipital generator.
4. Confirm arrows are visible against the opaque cortex (offset works) and that toggling/transparency sliders behave.
5. Threshold/decimation control reduces arrow count as expected; full-density time-stepping remains interactive (or a documented cap is applied).

## 8. Future work (out of scope here)

- Dirac-operator differential analysis/filtering of the field (divergence, curl, Helmholtz, flow, dispersion) in the eigen-coefficient domain.
- Optional band-limited display: reconstruct the arrow field from a chosen Dirac mode range via `ImagingKernelMode` (Transform/Reconstruct).
- Direction-encoded arrow coloring as an alternative to neutral color.
