# Stashed: connection-phase / sign-correction / wavefront research

**Stashed:** 2026-06-22 (Increment 5 of the file-based Dirac consolidation).
**Why:** these consumed the *legacy* standalone connection-Laplacian implementation
(`bst_conn_eigenmodes_ensure`, `in/out_tess_conn_eigenmodes`, `tess_conn_eigenmodes`,
`tess_connection_laplacian`, surface-embedded `ConnEigenmodes`), which has been **retired**.
The connection Laplacian is now a first-class operator variant, computed and stored via the
canonical `tess_operators('Connection Laplacian')` + `tess_eigen('Connection Laplacian')`
(`operator_` / `eigen_` DB nodes). This research cluster has no canonical equivalent yet, so
it is set aside here (off the Brainstorm path = inactive) to be **rewired** onto the `eigen_`
node in a later increment.

These files keep their now-dangling calls to the deleted legacy functions. They will NOT run
as-is; that is expected.

## Contents

- `bst_conn_phase.m` — decode a connection eigenmode into 3D tangent field + phase.
- `view_connection_phase.m` — the phase viewer (was launched from the operator_ node menu).
- `bst_source_sign_correct.m`, `bst_face_sign_correct.m` — sulcal sign-flip correction using
  the connection **Fiedler** magnitude as voting weights.
- `bst_wavefront_track.m`, `bst_face_wavefront_track.m` — wavefront tracking.
- `bst_cwt_fiedler_pipeline.m` — experimental CWT + Fiedler inverse pipeline.
- `connection_phase_readout_*.md` — the design/plan docs.
- tests/demos/benchmarks for the above.

## How to rewire onto the canonical eigen_ node

The legacy stored a **whole-mesh** `ConnEig` with `.Vectors [nV×nModes]` (complex) and
`Component`/`CompRank`/`Order` metadata. The canonical `eigen_` node (Variant
`'Connection Laplacian'`) stores **per-hemisphere** `Phi{1,2}` (complex `[nVh×K]`,
M-orthonormal) + `Lambda{1,2}` + `GlobalVertices{1,2}`.

Mapping:

```matlab
EM = tess_eigen(SurfaceFile, 'Connection Laplacian', 'nModes', K);   % or load the eigen_ node
% reassemble whole-mesh:
nV = ...; Phi_whole = complex(zeros(nV, K));
for hh = 1:2
    Phi_whole(EM.GlobalVertices{hh}, :) = EM.Phi{hh};
end
% the per-hemisphere Fiedler (smallest-eigenvalue mode) replaces the legacy
% Component==1 & CompRank==1 selection:
fiedler_left  = EM.Phi{1}(:, 1);
fiedler_right = EM.Phi{2}(:, 1);
```

Notes:
- The legacy split by connected *component*; the canonical splits by *hemisphere*
  (atlas L/R, never conncomp). For a standard cortex these coincide (2 components = 2 hemis).
- The new connection eigensolve (`tess_eigen` → `local_solve_connection`) drops the
  spurious non-positive discretization modes and returns only the genuine positive
  low-frequency modes (the legacy kept the spurious ones). When rewiring, expect a slightly
  cleaner mode set.
- The connection Laplacian operator is non-PSD by discretization; if you need the legacy's
  whole-mesh regularized assembly for reference, see git history for
  `toolbox/anatomy/tess_connection_laplacian.m` (deleted in Increment 5).
