# Design: Connection Laplacian onto the canonical operator/eigen path

**Date:** 2026-06-22
**Status:** Approved (design + "go ahead"); executing.
**Branch:** `refactor/file-based-dirac-consolidation` (Increment 5).
**Relates to:** `dev/2026-06-22-eigen-spectral-consolidation.md` (Increment 5).

## Goal

The connection Laplacian is one operator *variant* among Dirac / Laplace–Beltrami. Its
operator and its (complex, per-hemisphere) eigenmodes belong in the canonical `operator_`
and `eigen_` nodes via `tess_operators`/`tess_eigen` — which already support the
`'Connection Laplacian'` variant. Retire the parallel standalone implementation (a third
storage mechanism: surface-embedded `ConnEigenmodes`), and **stash** the phase /
sign-correction / wavefront research cluster — it has no canonical equivalent yet and will
be rewired onto the `eigen_` node in a later increment.

## Background (verified)

- **Canonical already supports it.** `tess_operators('Connection Laplacian')` →
  `A = nxr 'laplacian','connection'` (complex Hermitian, per-hemisphere) + `B = mass/galerkin`;
  `tess_eigen('Connection Laplacian')` → complex-safe Rayleigh-Ritz into per-hemisphere
  `Phi{1,2}` (complex, B-orthonormal) with real-positive `Lambda{1,2}`. GUI menus
  "Compute operator/eigenmodes → Connection Laplacian" already create these nodes.
- **Legacy standalone** computes a WHOLE-MESH `ConnEig` (`.Vectors [nV×nModes]` complex +
  `Component/CompRank/Order` block metadata), split by connected component, stored
  surface-embedded as `TessMat.ConnEigenmodes`.
- **No live consumer.** Every consumer of the legacy axis is the phase viewer/decoder, the
  experimental sign-correction/wavefront/cwt-fiedler cluster, or tests/demos. No `process_*`
  / `panel_*` / live `view_*` references it except the one tree menu item (removed here).

## Decisions (agreed)

- Stash the sign-correction/wavefront/cwt-fiedler experiments **with** the phase work (they
  share the connection-Fiedler dependency and are experimental/test-only).
- Stash mechanism: **move to `dev/stash/connection-phase/`** (off the Brainstorm path,
  inactive, git-preserved).

## ① Validate the canonical path (gate, no code change)

Compute on the cortex and assert:
- `tess_operators('Connection Laplacian')`: `Operator{hh}` complex Hermitian, `Mass{hh}` real SPD.
- `tess_eigen('Connection Laplacian','nModes',K)`: `Phi{hh}` complex `[nVh×K]`,
  `Phi{hh}'·B{hh}·Phi{hh} ≈ I` (B-orthonormal), `Lambda{hh}` real and ≥ 0.
- Optional pre-deletion snapshot: compare the per-hemisphere Fiedler (`Phi{hh}(:,1)`) magnitude
  to the legacy whole-mesh Fiedler for sanity (confidence only; legacy is being deleted).

## ② Delete (canonical replaces)

- `toolbox/anatomy/tess_conn_eigenmodes.m`
- `toolbox/anatomy/bst_conn_eigenmodes_ensure.m`
- `toolbox/anatomy/tess_connection_laplacian.m`  (legacy whole-mesh builder; callers = the above + its tests only)
- `toolbox/io/in_tess_conn_eigenmodes.m`, `toolbox/io/out_tess_conn_eigenmodes.m`
- Legacy tests: `dev/tests/test_conn_eigenmodes_{compute,ensure,roundtrip,manifold_guard}.m`,
  `dev/tests/test_connection_laplacian_{guard,nsym,smoke,spectrum}.m`

The surface-embedded `ConnEigenmodes` field is orphaned (nothing reads it after deletion) —
no DB migration needed; it is an unused extra field on any surface that has it.

## ③ Stash → `dev/stash/connection-phase/`

Move (git mv) out of `toolbox/`, off the path, with a `README.md` documenting the rewire path:

- `toolbox/math/bst_conn_phase.m`
- `toolbox/gui/view_connection_phase.m`
- `toolbox/inverse/bst_source_sign_correct.m`, `bst_face_sign_correct.m`
- `toolbox/inverse/bst_wavefront_track.m`, `bst_face_wavefront_track.m`
- `toolbox/inverse/bst_cwt_fiedler_pipeline.m`
- dev plans: `dev/connection_phase_readout_plan_{A,B,C}.md`, `dev/connection_phase_readout_integration.md`
- cluster tests/demos/benchmarks: `dev/tests/test_conn_phase.m`, `test_view_connection_phase.m`,
  `test_wavefront_pipeline.m`, `dev/tests/run_cwt_fiedler.m`,
  `dev/demo/connection_laplacian_demo.m`, `dev/demo/alpha_wave_phase_demo.m`,
  `dev/benchmarks/sign_ambiguity/` (whole dir)

These keep their now-dangling legacy calls (e.g. `bst_conn_eigenmodes_ensure`,
`in_tess_conn_eigenmodes`). The README states the rewire: consume the canonical
`eigen_` node via `tess_eigen('Connection Laplacian')`, reassemble per-hemisphere
`Phi{1,2}` to whole-mesh via `GlobalVertices`, and take the **Fiedler** as `Phi{hh}(:,1)`
(smallest-eigenvalue mode per hemisphere) — replacing the legacy `Component==1 & CompRank==1`
selection.

## ④ GUI disconnect

Remove the "View connection phase" tree menu item (`tree_callbacks.m:2563`) and its
`OperatorViewConnPhase_Callback` (`~4054–4072`), since `view_connection_phase` moves to the
stash. (The menu hung off the `operator_` node context menu.)

## Keep in place (not stashed, not deleted)

- `tess_tangents` — live (the constrained/loose face-leadfield frame in `bst_face_leadfield`).
- `bst_tangent_face2vertex` — generic per-face→per-vertex frame helper; orphaned after the
  stash but harmless and general; leave with its unit test.
- `bst_face_eigenmode_leadfield` — kept in Increment 4 (separate face-leadfield experiment).

## Validation

- 5a: live assertions above (canonical connection operator/eigen correctness).
- After 5b/5c: `grep` shows no `toolbox/` reference to the deleted legacy symbols or the
  stashed files; `checkcode` clean on the edited `tree_callbacks.m`; the Dirac/LBO paths
  untouched (regression: `tess_eigen('Dirac')` still reuses).
- Compute "Connection Laplacian" operator + eigen via the functions end-to-end and confirm
  the nodes register and load via the canonical `in_bst_operator`/`in_bst_eigen`.

## Out of scope

- Rewiring the phase/sign-correction work onto the `eigen_` node (a later increment).
- The `bst_eigen_spectral_*` umbrella (Increment 6).
- The legacy scalar-LBO Tier-1/2 retirements.

## Sub-commits

1. **5a** — validate the canonical connection path (live; no code).
2. **5b** — stash move (git mv) + `dev/stash/connection-phase/README.md` + GUI disconnect.
3. **5c** — delete the legacy standalone compute/IO + its tests.
