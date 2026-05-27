# Spec: Icosphere (FreeSurfer/MNE-style) cortex downsampling

## Goal

Add FreeSurfer/MNE-style **icosahedral** downsampling to Brainstorm so a cortex is resampled
**per hemisphere** to a standard ico level — ico3/4/5/6 = **642 / 2562 / 10242 / 40962 vertices per
hemisphere** — producing clean, closed, **2-manifold** meshes with uniform vertex spacing. This
replaces `reducepatch`'s inhomogeneous, non-manifold decimation (the defects `tess_manifold`
exists to repair) for users who want fsaverage-style cortices.

| FreeSurfer/MNE | per hemisphere | total cortex |
|----------------|----------------|--------------|
| ico3 | 642   | 1284   |
| ico4 | 2562  | 5124   |
| ico5 | 10242 | **20484** |
| ico6 | 40962 | 81924  |

## Key architectural finding (why per-hemisphere, and why at import)

A Brainstorm cortex is **two hemispheres concatenated** (`tess_concatenate({lh,rh})`), and that
function **stacks the registered spheres**: `Reg.Sphere.Vertices = [lh_sphere ; rh_sphere]`. Each
hemisphere's `?h.sphere.reg` is its **own full unit sphere centered at the origin**, so a combined
cortex's `Reg.Sphere` is **two coincident/overlapping spheres** (verified: both centroids ≈ origin,
both radius 0.1, both x ∈ [-0.10, 0.10]).

Consequence: running one ico grid against a *combined* cortex's sphere pulls vertices from
**whichever hemisphere is nearest**, yielding ~`nIco` *total* (not per-hemisphere) and ico faces
that **fuse the hemispheres** (verified: result had 1 connected component vs the original's 2).

Therefore icosphere must operate on **one hemisphere (one sphere)** at a time. FreeSurfer import
**already** downsamples each hemisphere separately (`tess_downsize(?hFile, nVertHemi, ...)`) *before*
`tess_concatenate`, so the per-hemisphere downsize step is the correct integration point. No
`tess_hemisplit` needed.

## Component 1 — `tess_downsize.m`: the `'icosphere'` method (single surface)

Operates on a **single** surface that has a FreeSurfer `Reg.Sphere` (i.e. one hemisphere).

1. Require `Reg.Sphere.Vertices` (matching vertex count); error clearly if absent.
2. **Guard**: compute connected components; if **> 1** (e.g. a merged L+R cortex), `bst_error`
   ("use icosphere via FreeSurfer import, per hemisphere") and abort — prevents the silent fusion.
3. Normalize the registered sphere to the unit sphere.
4. `[IcoVert, IcoFaces] = tess_sphere(newNbVertices)` — snaps to the nearest ico count.
5. Map each ico node to the nearest subject vertex (max dot product, chunked) → `sel`.
6. `resolve_ico_collisions` keeps the mapping injective (reassign duplicates to nearest unused).
7. Vertices = subject coords at `sel`; **faces = `IcoFaces`** (the ico grid's connectivity).
8. **Match the source surface's winding sign** (Brainstorm convention — required for forward
   modeling, normals, and correct rendering). Compute `sign(Σ v1·(v2×v3))` for both the source
   hemisphere and the ico mesh (both closed → valid winding indicator) and flip the ico faces if
   the signs differ. Do **not** force an absolute "outward" rule: Brainstorm's convention is the
   *opposite* sign, so forcing outward inverts the normals (renders dark / bumpy).
9. `I`/`J`: sort `sel` ascending, remap faces — exact reducepatch convention, so atlases / scouts /
   `Reg.Sphere` remap downstream untouched.
10. Override the comment/filename count to the actual `nIco`.

## Component 2 — `import_anatomy_fs.m`: per-hemisphere integration (primary path)

- New optional arg `Method` (default `'reducepatch'`), backward-compatible (arg #9).
- Interactive ASK block: a **method** radio (`reducepatch` / `icosphere`); if icosphere, an
  **ico-level** radio (ico3/4/5/6) → `nVertHemi = [642 2562 10242 40962]`. Non-interactive: derive
  `nVertHemi = round(nVertices/2)` (tess_downsize snaps).
- `Method` is threaded into the **6** per-hemisphere downsize calls (pial / white / mid × L/R). The
  **cerebellum** call stays `reducepatch` (no sphere registration).
- The mid surface (`tess_average({pial,white})`) inherits the pial sphere, so icosphere works on it.
- `tess_concatenate` offsets faces, so the combined `cortex_NV` has the correct `2 × nIco` vertices
  and **two separate components**. Naming (`cortex_20484V`) derives automatically from actual counts.

## Verification (done, live)

- `checkcode` on both files: no errors (only pre-existing style lint).
- Guard: combined cortex → **2 components (refused)**; single hemisphere → **1 component (allowed)**.
- Real FreeSurfer hemisphere (`fsaverage` lh.white) → `tess_downsize(...,'icosphere',ico5)` →
  **10242 vertices, 20480 faces, 1 component, `tess_check_manifold` = 1, Reg.Sphere carried, I/J = 10242**.

## Non-goals / notes

- **Select** nearest existing vertices (MNE-style), not barycentric interpolation — preserves the
  `I` index map for atlas/scout/registration transfer.
- Prefer **pure-ico** counts (the dialog enforces this); arbitrary counts snap to mixed grids.
- The same per-hemisphere `reducepatch` pattern exists in the CAT12/SimNIBS/BrainVISA/BrainSuite
  importers — generalizing icosphere to them is deferred (out of scope; FreeSurfer only for now).
- Existing buggy `cortex_ico_*` surfaces (produced by running icosphere on a combined cortex) are
  invalid (fused) and should be deleted/regenerated.
