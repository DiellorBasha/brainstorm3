# Helmholtz view: per-hemisphere bookkeeping + live trajectory tracking

**Date:** 2026-06-16
**Status:** approved (design)
**Builds on:** persistence detection (731d5a3e), tracker (dca6509e), marker sizing (e532b304).

## Motivation

Each FreeSurfer cortical hemisphere is a closed genus-0 surface (chi=2) and the two are
disconnected components, so singular-point bookkeeping (Morse/Poincare-Hopf counts of
psi/phi extrema) is a per-hemisphere invariant. Today `view_helmholtz` detects cores
per hemisphere but **gates and counts them whole-brain**, so a strong vortex on one side
can suppress real cores on the other and the readout lumps both. Separately, the viewer
shows only per-frame cores with no sense of motion. This change makes bookkeeping
per-hemisphere and adds a live, accumulate-on-play trajectory overlay with intrinsic-mean
(Karcher) centroids.

## Decisions (settled with user)

- Trajectory = **temporal, accumulated over contiguous forward play/scrub** (honors the
  active-frame rule; nothing precomputed over the series).
- Geodesics drawn by **mesh-edge Dijkstra** (nxr's edgeFlip/`query.line` is a stub; only
  heat-method distance is wrapped, and we don't need it here).
- **Karcher mean centroid** per trajectory via `nxr.manifold.query.center` (per hemisphere;
  capped to top-N tracks for responsiveness).
- Per-hemisphere **counts AND gate normalization**.

## A. Per-hemisphere bookkeeping

- `bst_dirac_helmholtz` `FindCoresOp`: add `.hemi` (= hemisphere loop index `hh`) to each
  core via `i_make_core`. Legacy full-mesh `FindCores` sets `hemi = 1`.
- `view_helmholtz` gate (UpdateFrame): keep core if
  `isGlobal OR persistence >= GateFrac * max(finite persistence within the SAME hemisphere)`.
- `i_readout` (vortex/source kinds): per-hemisphere line, e.g.
  `LH: 2 vortices(+), 1 antivortex(-) | RH: 1(+), 2(-)`. Falls back gracefully if a
  hemisphere has none.

## B. "Track trajectory" toggle (temporal, accumulated)

- `panel_helmholtz`: add a "Track trajectory" checkbox -> `view_helmholtz('SetTrack', hFig, isOn)`.
- `HelmholtzState` gains: `Track` (bool, default false), `LastIT`, `Tracks` (struct array:
  `.coreVerts` [1xL] ordered core vertices, `.path` [1xP] concatenated geodesic vertex
  indices, `.chirality`, `.hemi`, `.centroid` (1x3 or [])), `Graph` (cached mesh edge graph),
  `Ctx` (1xnH cached per-hemisphere nxr contexts), `V2H` (vertex->hemi map).
- On `UpdateFrame` with `Track` on:
  1. Reset rule: if toggled on, `|iT-LastIT|~=1`, iT decreased, or Component/smoothing
     changed since last -> clear `Tracks`, reseed each current gated core as a new track
     (path = its vertex, no segment). Set `LastIT=iT`; return after drawing.
  2. Contiguous step (`iT==LastIT+1`): take current gated cores (the displayed markers `mk`).
     Greedy-match each previous track head to a current core with **same chirality AND same
     hemi**, nearest Euclidean within `MaxJump` (default 0.012 m), one-to-one (sorted by
     distance). For a match: `path = [path, i_geodesic(Graph, prevVert, curVert)]`,
     append `curVert` to `coreVerts`. Unmatched current cores -> new tracks (birth).
     Unmatched heads -> stop extending (track stays drawn, dropped from heads).
  3. Draw: clear `Tag 'HelmholtzTrack'`; draw each track's polyline through `V(path,:)`
     colored by chirality (red `+`, blue `-`), `LineWidth` ~2.
- `SetTrack(false)`: clear `Tracks`, delete `'HelmholtzTrack'` + `'HelmholtzCentroid'` objects.

## C. Karcher-mean centroid

- Build `Ctx{hh} = nxr.manifold.context(Vloc, Floc)` once per hemisphere (lazy, cached in
  state); `Vloc/Floc` are the hemisphere submesh (reindexed local), same construction as
  `bst_dirac_helmholtz` Prepare / `tess_tangents`.
- After updating tracks, for the **top-N tracks** (default N=5) by current length (>=3 frames)
  by persistence: map `coreVerts` (global) -> local indices for that track's hemi, call
  `c = nxr.manifold.query.center(Ctx{hemi}, localVerts)`. If `c` is a vertex index ->
  `centroid = Vloc(c,:)`; if a 3x1 coord -> use directly (map local->ambient). Draw centroids
  as `Tag 'HelmholtzCentroid'` markers (e.g. filled diamond, chirality color, black edge).
- Cost guard: vector-heat solve per track per frame; capping to top-N (len>=3) bounds it.
  Context factorizations built once. Note in code: centroids may lag during rapid play and
  settle when the cursor stops.

## D. Components / files

- `toolbox/math/bst_dirac_helmholtz.m` — `.hemi` field (A).
- `toolbox/gui/view_helmholtz.m` — per-hemi gate + readout (A); Track state, UpdateFrame
  accumulation, geodesic draw, centroids, SetTrack (B,C); helper subfunctions
  `i_gate_per_hemi`, `i_link_step`, `i_geodesic`, `i_build_graph`, `i_hemi_context`.
- `toolbox/gui/panel_helmholtz.m` — "Track trajectory" checkbox + `OnTrack` callback (B).

## E. Testing (TDD)

Pure helpers unit-tested in a new `dev/tests/test_helmholtz_track.m`:
- `i_gate_per_hemi(mk, frac)`: two hemispheres with different max persistence; a core that
  would be cut by a global threshold but survives its own hemisphere's threshold is kept;
  globals always kept.
- `i_link_step(prevHeads, curCores, MaxJump)`: same-chirality+hemi match within radius;
  chirality mismatch and hemi mismatch block links; jump>radius -> birth+death.
- `i_geodesic(G, a, b)`: shortest path on a hand-built weighted line/grid graph equals the
  known path; same-vertex returns that vertex; disconnected returns empty.
Plus `bst_dirac_helmholtz` test: `Ht.Cores` carry `.hemi` consistent with `Op.vH`.
Integration: open `view_helmholtz` on the S01 Dirac source, enable Track, step frames ->
per-hemisphere polylines + centroids appear; readout splits LH/RH.

## Non-goals

Cross-frame identity persistence beyond contiguous play; trajectory export from the GUI
(the batch `process_vortex_track` already produces saved dipoles); smoothing the geodesic
beyond mesh edges; a MaxJump slider (fixed default for now).
