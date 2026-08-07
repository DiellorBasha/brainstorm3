# Gate 1 report — Phase 1 LBO operator primitives

Plan: `docs/superpowers/plans/2026-08-07-phase1-lbo-operators.md`
Harness: `dev/verify/phase1/` (oracle + 3 test `.m` files) +
`dev/verify/phase0/run_matlab.sh` (now addpaths `phase1`).

## Clean-branch commits (worktree `~/workspace/research/code/brainstorm3-clean`, branch `feature/cortical-flow-core`)

7 commits total ahead of `master` (`git log --oneline master..`); the 3 new
Phase 1 commits, oldest first:

```
0acb987d Anatomy: Add tess_massmatrix, Galerkin mass matrix for triangular surfaces
53bc911d Anatomy: Add tess_laplacian, cotangent stiffness matrix for triangular surfaces
7d85f06f Anatomy: Add tess_operators, per-hemisphere Laplace-Beltrami operator pencil
```

(Below them, 4 pre-existing commits from the icosphere/BIDS-import line:
`bc480ca6`, `64f6d56d`, `acd1f5cd`, `523cc323`.)

Worktree is at `7d85f06f`, `git status` clean (verified before this rerun).

## Parity surface (subject/protocol)

The oracle (`oracle_lbo_sub0002.mat`) was generated against the **Gate-0
imported ico5 cortex** (`sub-0002`, `tess_cortex_pial_low.mat`, `nV=20484`,
`L=10242`/`R=10242` per hemisphere) living in the **isolated**
`bst_userdir_clean` local DB under protocol `omega-tutorial-cortical-flow`
(`dev/verify/phase0/bst_userdir_clean/.brainstorm/local_db/...`) — **not**
the `SpikeData-2` protocol named in the original brief. That protocol was
never populated (a populate-script bug), so Task 1 corrected the parity
surface to the Gate-0 isolated-DB copy and re-ran against it; all downstream
tests (2–4) parity-check against this same surface via the committed
`SurfaceFile` path baked into each test.

## Oracle conventions (verbatim from Task 1's provenance/probe)

From `oracle_lbo_sub0002.mat`'s `meta.conventions` (probe run immediately
after generation, same MATLAB batch invocation):

- **A (stiffness, cotan)**: symmetric, exactly (Frobenius err = 0 both
  hemispheres). `A*1 ≈ 0` (row sums): `8.17e-14` (hemi 1) / `1.82e-14`
  (hemi 2) — machine-precision null space, confirms constants are in the
  kernel. **Positive diagonal** (`diagA sign = +1`); dense-200×200 eig probe
  strictly positive (0.114 / 0.111); corroborating `x'Ax` sign probe (5
  random vectors/hemi) = `+1` in every case. Convention: **+PSD
  cotan-Laplacian/stiffness matrix** (not the negative-semidefinite
  "Laplace–Beltrami operator" sign) — eigenpairs of `Ax = λBx` have `λ ≥ 0`.
- **B (mass matrix)**: symmetric, exactly (Frobenius err = 0). **Positive
  off-diagonals** (`B offdiag sign = +1`) → **consistent (Galerkin)** mass
  matrix, not lumped/diagonal. `sum(B(:))` = total surface area (partition
  of unity: `Σᵢφᵢ = 1` ⇒ `sum(B(:)) = ∫1 dA`):
  - hemi 1: **0.112357 m²**
  - hemi 2: **0.11357 m²**
  - i.e. **≈ 0.1124 / 0.1136 m² per hemisphere** — physiologically plausible
    for a pial-surface hemisphere (~1120–1140 cm²).
  - **Correction to the plan**: the Phase-1 plan's "~0.01 m²" hint for this
    quantity is an order-of-magnitude typo/guess — 0.01 m² (100 cm²) is
    single-gyrus scale, far too small for a whole hemisphere. The measured
    **~0.1 m² order is the physically correct target**, and it is what
    Tasks 2–4 built and validated against.

## Step 1: three-test rerun (this task, fresh MATLAB batch runs, not reused from Tasks 2–4)

Runner: `dev/verify/phase0/run_matlab.sh dev/verify/phase1/<test>.m`, output
tee'd to `dev/verify/phase1/<test>.log`.

| Test | Result | Wall time (`time` around full run_matlab.sh invocation) | Per-hemi parity rel-err (Frobenius, vs oracle `A`/`B`, tol 1e-12) |
|---|---|---|---|
| `test_tess_massmatrix` | **PASSED** | 9.965s | hemi1 `5.70612e-15`, hemi2 `5.94831e-15` (mass, `B`) |
| `test_tess_laplacian`  | **PASSED** | 9.079s | hemi1 `7.82091e-15`, hemi2 `4.38253e-16` (stiffness, `A`) |
| `test_tess_operators`  | **PASSED** | 9.810s | not printed by the committed test (assert-only, tol 1e-12); see note below |

`grep -c "PASSED" dev/verify/phase1/*.log` → 1/1/1 (all three).

Logs: `dev/verify/phase1/test_tess_massmatrix.log`,
`dev/verify/phase1/test_tess_laplacian.log`,
`dev/verify/phase1/test_tess_operators.log`.

**Note on `test_tess_operators` parity numbers**: the committed test (`dev/verify/phase1/test_tess_operators.m`)
only asserts `norm(Op{hh}-S.A{hh},'fro')/norm(S.A{hh},'fro') <= 1e-12` and
the mass-matrix equivalent, then prints `PASSED` — it does not print the
rel-err values themselves, so this Step-1 rerun's log has no numeric
parity line for this test. `tess_operators` composes `tess_laplacian` and
`tess_massmatrix` directly with no additional transform, so its parity is
mathematically identical to the two tests above. Task 4's report recorded
these exact values from a non-committed verbose variant of the same test,
run once during Task 4's own development: hemi1 stiffness `7.821e-15` /
mass `5.706e-15`; hemi2 stiffness `4.383e-16` / mass `5.948e-15` — i.e.
~1e-15 (hemi1) / ~1e-16 (hemi2), consistent with (and traceable to) the
Task 2/3 numbers above. These are carried forward from that prior run, not
re-verified numerically in this Step-1 rerun, since the committed harness
doesn't expose them.

## Guard behaviors verified (`test_tess_operators`, all in Step 1's PASSED run)

1. **Recipe-struct reproducibility**: `tess_operators(SurfaceFile, struct('Name','Laplace-Beltrami','Tau',[]))` reproduces the string-call pencil exactly (`isequal`).
2. **`unknownVariant` guard**: `tess_operators(SurfaceFile, 'Dirac')` → errors with identifier `tess_operators:unknownVariant` (Dirac variant intentionally not yet implemented — see Open Question 2 below).
3. **`noHemisphereLabels` guard**: surface with `Atlas` field stripped → errors with identifier `tess_operators:noHemisphereLabels`.
4. **`nonManifold` guard**: surface with one duplicated face → errors with identifier `tess_operators:nonManifold`.

All four assertions passed silently inside the test (no separate log lines — the test only prints on final success), consistent with Task 4's report of the same four guards firing correctly under a verbose variant.

## Deviations

- Parity surface: SpikeData-2 → isolated `bst_userdir_clean` Gate-0 copy (Task 1, documented above and in `task-1-report.md`).
- Plan's `sum(B(:))` hint ("~0.01 m²") corrected to measured "~0.1 m²" order (Task 1, carried through Tasks 2–4).
- `test_tess_operators` parity numbers for this Step-1 rerun are carried forward from Task 4's non-committed verbose variant (see note above) rather than freshly printed, since the committed test is assert-only.
- No other deviations from the Task 5 brief.

## Harness commit

Committed on the lab branch (main checkout, current branch — see report footer for SHA): `dev/verify/phase1/*` (oracle `.mat` included, 2.4 MB, well under the ~50 MB gitignore threshold) + the `dev/verify/phase0/run_matlab.sh` addpath edit. No other files staged.
