# tess_connectome + whole-brain connectome operators — design

**Goal.** Promote the prototype connectome/combined Laplacians into proper Brainstorm operators:
`tess_connectome` (find-or-create connectome from anatomy fibers) and two named operators in
`tess_operators` (`Connectome Laplacian`, `Combined Laplacian`), with a registration-based fallback that
lets any subject "borrow" the default HCP-1065 connectome when it has no DWI of its own.

## Reuse (do NOT duplicate Brainstorm)
- **Fibers are anatomy already** — `fibers_*.mat`, the subject `Fibers` field, `db_set_template`. The
  HCP-1065 fibers are bundled on `@default_subject` (Surface "HCP-1065"; see
  `build_hcp1065_default_fibers.m`).
- **`fibers_helper('AssignToScouts', FibMat, ConnectFile, ScoutCentroids)`** — assigns each fiber's two
  endpoints to the nearest centroid; returns `Scouts.Assignment` `[nFib x 2]`. Pass **scout centroids**
  for an ROI connectome or **all surface vertices** for a vertex connectome. REUSED as-is.
- **`tess_interp_tess2tess`** — surface→surface projection via FreeSurfer registration spheres (the
  group-analysis mechanism). REUSED for the template→subject mapping.
- **Operator find-or-create pattern** — `tess_operators` + `db_add_operator` + `in_bst_operator` +
  `bst_get('OperatorFileForSurface')` (mirrors `tess_eigen`). REUSED for storage/caching.

## New code (small)
1. Assemble the connectome matrix from `AssignToScouts` output (`C[a,b] += 1` per fiber).
2. Vertex-mode sulcal smoothing `W = S^p Wraw S^p'` (from `proto_fiber_connectome`: raw endpoints leave
   ~46% of vertices isolated; smoothing connects 99%).
3. The two Laplacian operators.

## tess_connectome(SurfaceFile, Opts)
- `Opts.Resolution` = `'roi'` (atlas scouts; needs `Opts.Atlas`, default Desikan) | `'vertex'`.
- `Opts.SmoothHops` (vertex mode, default 3), `Opts.ForceRecompute` (default false → reuse cached).
- **Fiber resolution:** subject's own `Fibers` if present; else the **default HCP-1065** fibers,
  endpoints projected onto this surface via `tess_interp_tess2tess(templateCortex, SurfaceFile)`
  (add registration spheres with `tess_addsphere` if missing). This is the "assume the HCP connectome"
  path — spherical surface co-registration, NOT streamline warping.
- Output: a cached connectome node (sparse `W`, the assignment, provenance: which fibers, which atlas,
  resolution, registration source).

## tess_operators: two named operators (stored as operator_ files)
- **`Connectome Laplacian`** — graph Laplacian of the `tess_connectome` `W` (symmetric normalized
  `I - D^-1/2 W D^-1/2`, whole-brain). Stored like `Laplace-Beltrami`.
- **`Combined Laplacian`** — `L = L_LBO (block-diag per-hemi) + gamma * L_connectome`. `gamma` balances
  geodesic vs network (default from the fiber-length / spectrum scale; see the eigenwavelet mm
  calibration). FIRST operator that returns a SINGLE whole-brain basis (not `{LH, RH}`): its eigenmodes
  (via `tess_eigen`) span both hemispheres, and the DEC (`bst_gradient/divergence/curl`) flows across
  the callosum. Keep per-hemisphere operators unchanged (geometry-pure) so existing callers don't break.

## Registration (the template→subject question) — RESOLVED
Possible via the registration sphere. `tess_interp_tess2tess` gives template-vertex ↔ subject-vertex
correspondence; map HCP fiber endpoints onto the subject cortex, then build the connectome on the
subject's OWN vertices. Caveats recorded: it's a group-average substrate (correct for default anatomy /
subjects lacking DWI; subject DWI supersedes when available), and the projection inherits the spherical
registration's accuracy.

## Build order
1. `tess_connectome` (ROI + vertex, own-fibers path) + validate vs `proto_fiber_connectome`.
2. The template→subject registration fallback (`tess_interp_tess2tess`).
3. `Connectome Laplacian` + `Combined Laplacian` in `tess_operators` (+ `tess_eigen` whole-brain path).
4. Operator-node storage + `db_update` migration (reuse the generic operator node).

## Open questions
- Operator-node tree type: reuse the generic `operator_` node (bst-java fork) vs a connectome-specific
  node — start generic.
- `gamma` default + whether to expose fiber-length-band stratification now or later (later).
