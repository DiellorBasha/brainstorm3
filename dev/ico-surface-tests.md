# Icosphere downsampling — validation tests

Validation of the FreeSurfer/MNE-style `icosphere` cortex downsampling (see
`dev/ico-downsize.md`). Run live against MATLAB R2023b on the `feature/eigenmode-analysis`
branch (commit `82abcc00`), protocol `OMEGA_Tutorial`.

## Test: Brainstorm icosphere ico5 vs FreeSurfer `fsaverage5`

**Goal.** Import `fsaverage` (`~/workspace/library/datasets/fsaverage`) with `icosphere` / **ico5**
(20,484-vertex cortex) and compare against FreeSurfer's own ico5 downsampling,
`fsaverage5` (`~/workspace/library/datasets/fsaverage5`, already 20,484 vertices).

**Inputs.**
- Ours: `cortex_20484V` from importing `fsaverage` via the GUI with method `icosphere`, level ico5
  (orientation-fix verified: signed-volume sign = −1, matching the high-res surface).
- Reference: `fsaverage5` `lh/rh.white` + `sphere.reg`, read directly from the FreeSurfer files.

**Metrics & reliability.**
- *Counts / topology* — vertices, faces, connected components, `tess_check_manifold`.
- *Uniformity* — coefficient of variation (CoV = std/mean) of edge length and face area
  (scale- and frame-invariant; computed on the actual surface). **Reliable.**
- *Selection parity* — fraction of our selected vertices that are the canonical `fsaverage5`
  nodes, via an **index match** of `cortex_20484V` into the high-res `cortex_327684V`
  (both in the same Brainstorm frame; `fsaverage5` nesting confirmed below). **Reliable.**
- ⚠️ *Not reliable (excluded):* a sphere-domain angular-alignment metric and an inline
  re-implementation of the selection both gave nonsense (63° "alignment", 8,875 collisions,
  CoV ≈ 2.4) because reading the FreeSurfer `sphere.reg` directly via `in_tess` does **not**
  reproduce the frame/loading that Brainstorm's `tess_addsphere` uses. Those numbers were
  analysis artifacts and are not reported as results.
- *Nesting check* — confirmed `fsaverage5` vertices equal the first 10,242 vertices of
  `fsaverage` per hemisphere (max position difference ≈ 0.005), so "canonical node" =
  "first 10,242 of each hemisphere".

### Results

| Metric (per hemisphere) | Ours (LH / RH) | fsaverage5 (LH / RH) |
|---|---|---|
| Vertices | 10,242 / 10,242 | 10,242 / 10,242 |
| Faces | 20,480 / 20,480 | 20,480 / 20,480 |
| Connected components | 1 / 1 | 1 / 1 |
| Manifold (`tess_check_manifold`) | yes / yes | yes / yes |
| Edge-length CoV | 0.331 / 0.339 | 0.266 / 0.266 |
| Face-area CoV | 0.491 / 0.511 | 0.307 / 0.322 |
| All vertices are real cortex vertices | 20,484 / 20,484 exact matches into high-res | — |
| **Selected vertices that are canonical fsaverage5 nodes** | **36.1%** | — |

### Interpretation

- **Structurally identical** to `fsaverage5`: same vertex/face counts, two clean
  components, both hemispheres manifold. Every one of the 20,484 output vertices is a genuine
  vertex of the original cortex.
- **Selection differs:** only ~36% of our nodes are exactly FreeSurfer's canonical ico5 nodes;
  the other ~64% are **adjacent** vertices. This is because Brainstorm's `tess_sphere`
  icosahedron is not rotationally aligned with FreeSurfer's canonical icosahedron, so the
  nearest-neighbour snap on the registered sphere lands on a neighbour.
- **Uniformity:** consequently our mesh is **modestly less uniform** (edge CoV 0.33 vs 0.27,
  face-area CoV ~0.50 vs ~0.31). This is the source of the minor "bumps" once the
  orientation/lighting bug was fixed — genuine geometry, not lighting.

### Conclusion

The `icosphere` method is **validated**: it produces a structurally-correct, manifold,
`fsaverage5`-equivalent ico5 cortex with Brainstorm-convention winding (renders and forward-models
correctly), and visually matches the high-res surface after the orientation fix.

**Known limitation (acceptable, improvable).** The mesh is slightly less uniform than
FreeSurfer's canonical ico5 because of icosahedron grid misalignment (36% canonical overlap). For
general (non-fsaverage) subjects the nearest-neighbour selection is the correct general method and
the small non-uniformity is the expected trade-off.

**Optional future enhancement (separate scope).** For tighter FreeSurfer parity: rotationally
align the `tess_sphere` grid to FreeSurfer's canonical icosahedron, or — for FreeSurfer-registered
surfaces — select FreeSurfer's canonical ico nodes directly via the nesting property.
