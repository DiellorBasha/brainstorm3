# Gate 2 = Milestone 1 report — Laplace-Beltrami eigenmodes (tess_eigen)

Clean-branch checkout: `~/workspace/research/code/brainstorm3-clean`
HEAD: `5840d76c` — "Anatomy: Add tess_eigen, embedded per-hemisphere Laplace-Beltrami
eigenbases" (8 commits off `master`; full list at the bottom).

## 1. Pinned solver parameters (Open Question 1 — RESOLVED)

From `experiment_shift_invert.m` / `shift_experiment_results.md` (K=400, real
cortex pencil, both hemispheres of sub-0002):

**PINNED: `TauRel=1e-4`, `StartVector=deterministic-ones`**

Decision-rule application:

- **Shift-invert arms** (TauRel in {1e-6, 1e-4, 1e-2}, both hemispheres):
  `maxres` in **1.14e-10 – 3.69e-10**, tau-invariant (same order of magnitude
  regardless of tau); orthonormality `orth` in 5.19e-14 – 6.09e-14.
  Deterministic-start reproducibility (`det-rep`) and cross-arm subspace
  disagreement are both at **machine epsilon** (0, ±2.22e-16) across all three
  tau values on both hemispheres — the choice of tau is immaterial to the
  result, so per the decision rule's tie-break, the middle value 1e-4 is
  pinned as default.
- **sigma=0 (`smallestabs`) arm — doubly disqualified:**
  - *Combined run (sigma=0 attempted first, single background invocation):*
    process-fatal. MATLAB emitted `RCOND = 1.650491e-16` ("First input matrix
    is close to singular or badly scaled") on Hemisphere 1's factorization and
    the batch process died mid-run — no second arm was ever reached, no
    results file survived beyond the Hemisphere-1 header line.
  - *Isolated rerun (sigma=0 alone, after shift-invert results already safely
    on disk):* same RCOND=1.65e-16 warning fired but the process survived
    this time — and returned numerically catastrophic results: **maxres=0.992
    (Hemi 1) / 0.991 (Hemi 2)**, ~10 orders of magnitude worse than the
    shift-invert arms' ~1e-10 floor. `rand-rep` ~0.07 also shows the returned
    subspace is sensitive to the random start vector, unlike shift-invert's
    rand-rep ~0.
  - Conclusion: sigma=0 is either **fatal** (kills the process,
    non-deterministically) or **silently wrong** (fails the residual test by
    ~10 orders of magnitude). Both outcomes are disqualifying independent of
    the exact residual threshold — this is the resolution of Open Question 1.
    `tess_eigen` uses shift-invert exclusively.

## 2. Sphere spectrum verification (`test_tess_eigen_sphere.m`)

Analytic check against l(l+1) spherical-harmonic eigenvalues (K=50, unit
two-sphere mesh, 2562 verts/hemi) plus storage round-trip (embed / reuse+
truncate / ForceRecompute replace / History row) and a residual+orthonormality
check on the recomputed basis.

- **Max relative error vs analytic spectrum: 0.021** (assert threshold 0.05).
- Pencil residual < 1e-10, B-orthonormality (Frobenius) < 1e-10 — both met.
- Rerun result: **PASSED** (2026-08-08 rerun, ~10s wall time via
  `dev/verify/phase0/run_matlab.sh`). Log: `test_tess_eigen_sphere.log`.
- Non-fatal MATLAB warnings noted in the log ("Ignoring issym field...") are
  benign — `eigs` receiving a matrix-free/struct options combination it
  chooses to ignore a redundant symmetry hint on; does not affect results.

## 3. Cortex verification vs the independent nxr oracle pencil (`test_tess_eigen_cortex.m`)

K=400/hemisphere on the real ico5 cortex (sub-0002, `tess_cortex_pial_low.mat`,
10242 verts/hemi), checked against the **independent** nxr-built oracle pencil
(`dev/verify/phase1/oracle_lbo_sub0002.mat`) — a cross-implementation check,
not a self-consistency check.

- Residual vs independent nxr pencil: **7.7e-10 (hemi 1) / 1.5e-10 (hemi 2)**
  (assert threshold 1e-8).
- Orthonormality: **6.7e-14** (assert threshold 1e-10).
- Zero mode recovered (`lambda1` ≈ 0 relative to `lambda400`) on both
  hemispheres.
- Solve vs reuse timing (Task 3 reference run): **solve 5.7s, reuse 0.24s**
  (reuse >10x faster, as required).
- File growth on embedding K=400 basis for both hemispheres: **+63.6 MB**
  (2×10242×400 doubles ≈ 65.5 MB raw; `-v7` gzip keeps most of it since
  eigenvectors are smooth).
- `in_tess_bst` passes the embedded `Eigen` field through untouched.

**Rerun deviation (2026-08-08):** re-running `test_tess_eigen_cortex.m` this
time reproduced the correctness numbers essentially exactly (residual
7.69e-10 / 1.48e-10, orth 6.73e-14 / 6.81e-14, matching the reference run
above to 3 significant figures) but the storage-growth assertion **failed**:
`solve time: 0.2 s`, file size unchanged (65.6 MB → 65.6 MB), tripping
`assert(info1.bytes > info0.bytes + 40e6, ...)`. Root cause: Task 4's
`make_eigenmode_movie.m` runs against the isolated protocol's live
`sub-0002/tess_cortex_pial_low.mat` and *intentionally* embeds the K=400
`Eigen` field permanently into that shared source file (by design — "this
becomes the standing Gate-2 demo protocol", per that script's own header
comment). `test_tess_eigen_cortex.m` copies from that same source file
(`Orig`) into a scratch working copy, so once Task 4 has run, the scratch
copy already carries the cached `Eigen` field before `tess_eigen` is ever
called — the "solve" call in the rerun is actually a cache **reuse** (hence
0.2s, matching the reference run's 0.24s reuse timing almost exactly) rather
than a fresh solve, and there is nothing left to grow. This is a test-harness
sequencing artifact (Task 4 mutating shared state that Task 3's test assumes
is pristine), not a regression in `tess_eigen` correctness — the numeric
checks that matter (residual vs. independent oracle, orthonormality, zero
mode) all reproduced. Full log: `test_tess_eigen_cortex.log`.

## 4. Eigenmode demo artifacts (Task 4)

- `make_eigenmode_movie.m` embeds the first 100 (of the K=400 computed) LBO
  eigenmodes as a Brainstorm results file (mode index = "time" axis) into the
  isolated protocol under `sub-0002 > @intra`, comment "LBO eigenmodes 1-100
  (mode=frame)".
- Results file path (isolated protocol, protocol-relative under sub-0002's
  `@intra` analysis study): `results_lbo_eigenmodes*.mat` — created via
  `db_add_data` / `bst_process('GetNewFilename', ...)`; open via
  `dev/verify/phase0/launch_gui_clean.sh`, sub-0002 > @intra > "LBO
  eigenmodes 1-100 (mode=frame)".
- Snapshot PNGs captured for modes 2, 10, 50 (headless `out_figure_image`
  succeeded — no skip needed):
  - `dev/verify/phase2/eigenmode_002.png`
  - `dev/verify/phase2/eigenmode_010.png`
  - `dev/verify/phase2/eigenmode_050.png`
- **Viewing note:** Brainstorm's default amplitude-threshold slider hides
  sub-peak structure on these maps — lower the threshold slider in the GUI
  to see the full harmonic pattern (higher-order modes have many small
  lobes that sit under the default threshold).

## 5. Milestone 1 status: COMPLETE

Gate 2 = Milestone 1 (Laplace-Beltrami eigenmodes: massmatrix → laplacian →
operators → eigen, cotangent stiffness through to a validated, cached,
GUI-visible per-hemisphere eigenbasis) is complete. 8 clean commits on the
lab clean-branch checkout (`~/workspace/research/code/brainstorm3-clean`,
`master..HEAD`, oldest first):

1. `523cc323` — Anatomy: Add icosphere (FreeSurfer/MNE-style) per-hemisphere cortex downsampling
2. `acd1f5cd` — IO: FreeSurfer import cortex downsampling method option (icosphere, default for auto-import)
3. `64f6d56d` — Bugfix: Comment out unused sun.misc.BASE64Decoder import (unresolvable on Java 11+ / MATLAB R2023b+)
4. `bc480ca6` — IO: BIDS import cortex downsampling method and resolution options (default icosphere ico5)
5. `d190fbfa` — Anatomy: Add tess_massmatrix, Galerkin mass matrix for triangular surfaces
6. `6fed101c` — Anatomy: Add tess_laplacian, cotangent stiffness matrix for triangular surfaces
7. `b08b2a42` — Anatomy: Add tess_operators, per-hemisphere Laplace-Beltrami operator pencil
8. `5840d76c` — Anatomy: Add tess_eigen, embedded per-hemisphere Laplace-Beltrami eigenbases

## 6. Open follow-on decisions (not resolved here)

- **Open Question 2:** Dirac variant of the eigen path (deferred to Phase 3+).
- **Phase 3 scoping:** filtering framework vs. Dirac kernel projection first —
  still to be decided by the user.

## 7. Deviations from the task brief

- **Step 1 rerun of `test_tess_eigen_cortex.m` failed its storage-growth
  assertion** (see §3 above) due to Task 4 having permanently embedded the
  `Eigen` field into the shared isolated-protocol source file that this test
  copies from. The correctness assertions (residual vs. independent oracle,
  orthonormality, zero-mode recovery) all reproduced successfully; only the
  "growth implies a fresh solve happened" assertion fired, because in this
  post-Task-4 environment the copy source is no longer pristine. Reported
  here as a deviation per the brief's request rather than silently patched;
  the authoritative solve/reuse timing and file-growth numbers reported in
  §3 are from the original Task 3 reference run, captured before Task 4
  mutated the shared source file. No code or test changes were made as a
  result — this is a harness/state-sequencing note, not a `tess_eigen` defect.
- `test_tess_eigen_sphere.m` reran cleanly with no deviation.
- No other deviations from the brief. Step 3 commit stages only
  `dev/verify/phase2/` (force-added, git-ignored dir) and
  `dev/verify/phase0/run_matlab.sh`; pre-existing dirty files in the
  working tree (`dev/demo/_figures/*.png` deletions, `toolbox/dynamics/
  bst_dynamics.m`, `source_script_run.png`) are untouched, per instruction.

## Fix wave (post final review)

The Milestone-1 final review requested one hardening fix, amended in place
(tess_eigen commit c2ef53f5 -> 5840d76c; prior 7 SHAs unchanged): the
cache-reuse guard now verifies Structures-atlas presence (L/R labels) and
per-hemisphere row counts via tess_hemisplit before serving a cached basis
(protects against the known atlas-wipe gotcha), plus a forward-comment on the
stored Tau recipe field for the future Dirac variant. Sphere test extended
with an atlas-strip guard case (asserts tess_operators:noHemisphereLabels);
full test re-passed. Deferred to Phase 3+: VariantToField collision edge,
issym/datestr cosmetics, per-surface write serialization if batch concurrency
arrives.
