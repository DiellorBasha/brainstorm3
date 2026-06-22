# The Brainstorm Scouts System — reference & expansion notes

Reference map of Brainstorm's cortical region-of-interest ("Scouts") subsystem, written so we
can come back to it when expanding scouts with the file-based `bst_eigen` differential-geometry
system. Line numbers are indicative (large files drift) — treat them as starting points.

Status: 2026-06 reference. The "Expansion" section (§7) is the forward-looking design, not yet built.

---

## 1. What a Scout *is* (one paragraph)

A **Scout** is a flat, unordered **set of vertex indices** into a cortical surface mesh, plus a
**Seed** vertex and a scalar **reduction Function** (mean / max / PCA / power / std). It is a
boolean mask whose only geometric structure is mesh adjacency (the surface's sparse `VertConn`).
There is no basis, no spectrum, no orientation/flow, no metric beyond 0/1 connectivity, and no
multi-scale structure. Everything the system does — grow, edit, display, reduce — is built on that
mask. The single place it brushes against spectral geometry is `Function = 'pca'`, which projects
the patch onto the dominant eigenvector of the **data** covariance (a degenerate one-mode version of
what an operator eigenbasis projection does).

---

## 2. File ownership map

| Concern | File |
|---|---|
| Data model | `toolbox/db/db_template.m` — `db_template('Scout')`, `db_template('Atlas')`, `surfacemat` |
| GUI / lifecycle | `toolbox/gui/panel_scout.m` (the Scout tab; create/grow/edit/save/display, ~5700 lines) |
| Grow (region swell) | `toolbox/anatomy/tess_scout_swell.m` (one-hop neighbor grow over `VertConn`) |
| Reduction (core op) | `toolbox/math/bst_scout_value.m` |
| Apply scouts to data | `toolbox/process/functions/process_extract_scout.m` |
| Import atlases | `toolbox/io/import_label.m` (FreeSurfer `.annot`/`.label`, `.gii`, `_scout.mat`) |
| Display | `panel_scout.m` `PlotScouts` (colored `patch` + boundary contour + seed marker + label) |

---

## 3. Data model — scouts live *inside* the surface file

A `tess_*.mat` surface carries an **array of atlases**, each holding a `.Scouts` array; `iAtlas`
selects the active one. So one surface can carry many parcellations simultaneously (e.g.
'User scouts', 'Structures', Desikan-Killiany).

```matlab
% db_template('surfacemat')  (toolbox/db/db_template.m ~105-119)
surfacemat.Atlas   = [Atlas ...]   % 1 x nAtlas
surfacemat.iAtlas  = 1             % 1-based active atlas
surfacemat.VertConn= []            % sparse [nV x nV] logical adjacency (shared-face) -> grow/connectivity
surfacemat.Vertices= []            % [nV x 3] mm

% db_template('atlas')  (~902-905)
Atlas = struct('Name','User scouts', 'Scouts', [Scout ...])

% db_template('scout')  (~912-926)
Scout = struct( ...
  'Vertices', [], ...   % indices into the surface mesh (the mask)
  'Seed',     [], ...   % initial/centroid vertex (region grow origin, display center)
  'Color',    [], ...   % [R G B] in [0,1]
  'Label',    '', ...   % name, e.g. 'M1 R'
  'Function', 'Mean', ...% reduction: Mean|Max|Power|Std|PCA|FastPCA|All|None
  'Region',   'UU', ... % 2-char: {L,R,U} x {F,P,T,O,C,U} (hemisphere x lobe) for auto-color/metadata
  'Handles',  [] );     % cached graphics handles (hPatch,hContour,hLabel,hScout) -- stripped before save
```

Key takeaway: **a scout = `Vertices` (mask) + `Seed` + `Function`**. The rest is cosmetic/metadata.

---

## 4. Lifecycle (create / grow / edit / save)

- **Create** (`panel_scout.m` CreatePanel + add-scout callback): "Create scout" toggle; click a
  vertex -> `Seed`, `Vertices=[Seed]`; shift-click to add/remove vertices.
- **Grow / Shrink**: `tess_scout_swell(verts, VertConn)` returns the **one ring of neighbors**:
  ```matlab
  newverts = find(max(vconn(iverts,:), [], 1));   % all vertices one hop from the scout
  newverts = setdiff(newverts, iverts);
  ```
  Grow/Shrink iterate this. Pure mesh adjacency, **unweighted** (every edge counts as distance 1).
  -> This is exactly the limitation §7 fixes: swell by *geodesic diameter*, not hop-count.
- **Edit**: rename (`Label`), recolor (`Color`), set reduction (`Function`), merge/split.
- **Persist** (`SetScouts` -> `GlobalData.Surface(iSurf).Atlas(iAtlas).Scouts`, flag
  `isAtlasModified`; `SaveScouts` strips `Handles` and `bst_save`s into the surface `.mat`, or
  exports to FreeSurfer `.annot`/`.label`).

---

## 5. The reduction — the heart of the system (`bst_scout_value.m`)

Given source data `F` (`[nVert x nTime]` constrained, or `[3*nVert x nTime]` unconstrained) for a
scout's rows, collapse to `[1 x nTime]` (or `[k x nTime]`) by `Function`:

| Function | Operation |
|---|---|
| `mean` | `mean(F,1)` |
| `power` | `mean(F.^2,1)` |
| `std` | `std(F,[],1)` |
| `max` | signed max-abs vertex |
| `pca` / `pca2023` | project onto top eigenvector of the patch **data covariance** (`eig(Covar)`), sign-handled |
| `all` / `none` | no reduction (return all vertex series) |

Unconstrained (3-vector) sources are flattened first via `XyzFunction` (`norm` = `mean(‖·‖)`, or
`pca`), with optional normal-based **sign-flip** for orientation consistency (`isSignFlip`,
constrained case). `process_extract_scout.m` is the wrapper: load results -> select the scout's
source rows -> `bst_scout_value` per scout -> one matrix row per scout (a `matrixmat`).

---

## 6. Display (brief)

`PlotScouts` (`panel_scout.m`): build a MATLAB `patch` from the faces whose 3 vertices are all in
the scout (renumbered to local indices), `FaceColor=Color`, `FaceAlpha` user-set, `Tag='ScoutPatch'`;
plus a boundary contour, a seed-marker sphere, and a text label. Handles cached in `Scout.Handles(iFig)`.

---

## 7. Expansion via `bst_eigen` (the forward-looking design)

A scout is `{vertex mask} + {scalar reduction}`. The `bst_eigen` system replaces the **data**
covariance with the **operator** eigenbasis and adds differential-geometry structure scouts can't
represent. The mapping:

| Scout concept | Eigen / differential-geometry analogue |
|---|---|
| Region = vertex mask | Region = **spectral subspace** (a band of operator modes), or a mask lifted to its restricted eigenbasis |
| Reduce = mean/max/pca | Reduce = **basis projection** `C = Phi' * B * F` (`manifold_ft`) — coefficients, not an average |
| `pca` = top covariance mode | top **operator** mode (LBO / Connection / Dirac) — geometry-aware, data-independent |
| Grow via `VertConn` (0/1 hops) | **metric-aware grow**: heat-kernel / geodesic distance, isolines (see below) |
| Scalar field, no orientation | **vector/flow aware** (Dirac, Hodge, connection Laplacian) |
| Single scale | **multi-scale** via eigenvalue bands / spectral wavelets (the `eigfilter` tiles) |

Question a scout answers: *"what is the average activity in this patch?"*
Question the eigen version answers: *"what is this region's spectral signature / flow / multi-scale
content?"* — and the projection/reduction machinery already exists (`manifold_ft`/`ift`,
`bst_eigenspectrum`, `bst_eigenfilter`, the `eigfilter` bank).

### 7a. Principled cortical distance (heat method) — the key idea

Replace hop-count growth with a **real geodesic distance field**, via the geometry-central / nxr
heat method (Crane et al.):

1. **Heat filter localized to a vertex**: a delta at the seed, low-passed by a heat kernel
   `exp(-t*lambda)` in the operator eigenbasis (the `eigfilter` 'heat' kernel) — i.e. short-time
   heat diffusion from the seed. (`bst_eigenfilter` applied to a seed delta = a localized heat bump.)
2. **Gradient** of that diffused field on the surface (`nxr` face gradient / `gradFace`).
3. **Normalize** the gradient vector field to unit magnitude (the heat-method trick: the normalized
   gradient approximates the *direction* of the geodesic distance gradient).
4. **Poisson solve** (scalar LBO) for the scalar field whose gradient best matches that unit field:
   `lap phi = div(normalized grad)` -> `phi` = the **geodesic distance** from the seed.

This gives a principled distance-on-the-cortex function (curvature/metric aware), unlike `VertConn`
hops. nxr-compute already exposes the pieces (heat, gradFace, lapFace/Poisson, divergence) used by
geometry-central's `HeatMethodDistanceSolver`.

### 7b. Isolines for diameter-based swelling

With a true distance field `phi(seed)`, **swell a scout by geodesic diameter** rather than
connectivity: the scout of radius `r` = `{ phi <= r }`, and its boundary = the **isoline**
`phi = r`. This replaces `tess_scout_swell`'s ring-by-ring hop growth with a smooth, isotropic,
metric-correct dilation — "grow by 5 mm" instead of "grow by 1 ring". `nxr` isoline extraction
(`isolines`) gives the boundary contour directly.

### 7c. An "eigen-scout" (first concrete step)

A region defined by `(eigen_ node + mode band)` (or a vertex mask lifted to its restricted
eigenbasis), whose reduction is `bst_eigenspectrum` / `manifold_ft` over the masked vertices — i.e.
`process_extract_scout`'s analogue driven by `bst_eigen` instead of `bst_scout_value`. Output: the
region's eigen-coefficients / spectrum, not a single mean time series.

---

## 8. Pointers

- Operators/eigenbasis: `tess_operators`, `tess_eigen` (LBO / Connection Laplacian / Dirac),
  `in_bst_operator` / `in_bst_eigen` (the file-based `operator_`/`eigen_` nodes).
- Projection primitives: `toolbox/math/manifold_ft.m`, `manifold_ift.m`.
- Eigen analysis: `toolbox/eigen/bst_eigen.m` (orchestrator), `bst_eigenspectrum.m`,
  `bst_eigenfilter.m`, `eigfilter/` (kernel bank incl. `heat`).
- Differential geometry backend: the nxr-compute plugin (`nxr_compute('operators', ...,
  'gradFace' | 'lapFace' | ...)`), geometry-central heat-method distance.
- Helmholtz/Hodge (vector-field decomposition, already file-based): `bst_dirac_helmholtz`,
  `bst_dirac_helmholtz_face`.
