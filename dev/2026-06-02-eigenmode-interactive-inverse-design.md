# Eigenmode source mapping in the interactive "Compute sources" dialog

**Date:** 2026-06-02
**Author:** Diellor Basha
**Status:** Design approved, pending implementation plan

## Motivation

Eigenmode source mapping currently runs only through the Process box
(`process_eigenmodes_inverse`). Users expect to reach it the same way as every
other inverse method: right-click a recording → **Compute sources [2018]** →
pick a method alongside Minimum norm / dSPM / sLORETA / LCMV / Dipoles. This
work wires the eigenmode solver into that interactive dialog and refactors the
Process-box entry into a thin wrapper, following Brainstorm's standard
`Compute()`-core pattern.

## Background: how the inverse pipeline is wired today

The interactive and batch paths share one core, `process_inverse_2018('Compute', …)`:

```
right-click → "Compute sources [2018]"        tree_callbacks.m:1391
   → selectHeadmodelAndComputeSources()        tree_callbacks.m:3546
   → panel_protocols('TreeInverse', …)          panel_protocols.m:1355
   → process_inverse_2018('Compute', …)          process_inverse_2018.m:130
        ├─ Run(sProcess, sInputs)  → Compute(…)   (batch / process box)
        ├─ gui_show_dialog(@panel_inverse_2018)   process_inverse_2018.m:245
        └─ switch(OPTIONS.InverseMethod)          process_inverse_2018.m:665
              case 'minnorm','gls','lcmv' → bst_inverse_linear_2018
              case 'mem'                  → be_main
```

- `panel_inverse_2018.m:137-170` defines the Method radio buttons.
- `panel_inverse_2018.m:155-167` already gates methods by `HeadModelType`
  (`'mixed'` disables beamformer/dipoles, `'volume'` constrains orientation).
- `panel_inverse_2018.m:626-747` (`GetPanelContents` / `GetSelectedMethod`)
  returns the chosen `InverseMethod` / `InverseMeasure`.
- Generic results-node creation lives at `process_inverse_2018.m:764-810` and is
  already `HeadModelType`-aware (`'surface'`/`'volume'`/`'mixed'`).

The eigenmode head model differs structurally: `Gain` is `[nCh × K]` (composed
`L·Φ`) with an `isEigenmode = 1` flag and an `Eigenvalues` field. The standard
linear/beamformer/dipole/MEM solvers cannot run on it, and conversely
`bst_inverse_eigenmodes` (`toolbox/inverse/bst_inverse_eigenmodes.m`) rejects a
normal head model (`:60-64`).

## Goals

- Add an **"Eigenmode source mapping"** method to the interactive dialog,
  available only when the active head model is an eigenmode leadfield.
- Reuse the shared `Compute()` core — no parallel inverse stack.
- Keep `process_eigenmodes_inverse` working for batch/scripted pipelines, as a
  thin wrapper over the shared core.
- Preserve the existing outputs: cortex sources node (`Φ·M̃` kernel) and an
  optional eigenmode-coefficients matrix node.

## Non-goals

- No change to the eigenmode forward model (`bst_eigenmode_leadfield`) or the
  "Cortex surface harmonics" head-model option.
- No change to the mode-space solver math in `bst_inverse_eigenmodes`.
- No new spectral-prior formulas (reuse `bst_eigenmode_prior`).

## Design decisions (resolved in brainstorming)

1. **Dialog UI** — a dedicated "Eigenmode source mapping" radio plus its own
   Spectral-prior sub-panel; standard methods are disabled when an eigenmode
   head model is active (mutually exclusive gating).
2. **Outputs** — cortex sources by default, with a checkbox to also save the
   eigenmode-coefficients matrix node.
3. **Code structure** — solve + reconstruct live in shared helpers called from
   the `Compute()` core; no duplicated math.
4. **Process box** — `process_eigenmodes_inverse` becomes a thin wrapper over
   `process_inverse_2018('Compute', …)`, matching how Brainstorm processes wrap
   their interactive counterparts.
5. **Whitener** — reuse Brainstorm's per-modality whitener by copying it
   verbatim into a shared helper; do not modify the proven
   `bst_inverse_linear_2018`. **Whitening is ON by default** so the eigenmode
   inverse replicates the standard sensor-data procedure (channels, SSP, bad
   channels, average reference, whitener) used by every other method. A toggle
   lets the user drop the whitener (identity `iW`, "pure" solve) when explicitly
   chosen; SSP/bad-channels/avg-ref still apply. See "Whitening".

## Architecture & data flow

```
right-click → "Compute sources [2018]"  ──┐
                                          ├─→ process_inverse_2018('Compute', …)
process_eigenmodes_inverse (batch) ───────┘        │  (shared core)
   (thin wrapper: packs OPTIONS,                     │
    calls Compute)                                   │  switch InverseMethod
                                                     ├ minnorm/gls/lcmv → bst_inverse_linear_2018
                                                     ├ mem              → be_main
                                                     └ eigenmode        → bst_inverse_eigenmodes
                                                                            → reconstruct Φ·M̃
                                                                            → (opt) coefficients M̃·d
```

Three helpers carry the eigenmode-specific work:

- **`bst_inverse_eigenmodes('SolvePure', …)`** (exists) — mode-space MAP solve;
  given the compressed leadfield `L_tilde`, eigenvalues, an externally-supplied
  whitener `iW`, and SSP projector `Proj`, returns the mode-space kernel `M̃`
  `[K × nGoodCh]`. The `'SolvePure'` entry already accepts injected `iW`/`Proj`
  (`bst_inverse_eigenmodes.m:6,33`), so the interactive path reuses Brainstorm's
  whitener instead of the solver's internal one.
- **`bst_noise_whitener`** (new) — a **verbatim** extraction of the per-modality
  whitener from `bst_inverse_linear_2018.m` (the block at `:347-449` plus the
  local subfunctions `truncate_and_regularize_covariance` `:1098` and
  `cov1para_local` `:1236`). See "Whitening" below.
- **Reconstruction helper** `bst_eigenmode_reconstruct` (new, standalone file) —
  loads Φ from the surface (`in_tess_eigenmodes`), returns the cortex kernel
  `Φ·M̃` `[nVert × nGoodCh]` and, when requested, the mode-space kernel `M̃` for
  the coefficients node. Kept as its own file so the interactive case and the
  batch wrapper share one code path.

## Component changes

### `panel_inverse_2018.m` (dialog)

- **Method radio** — add "Eigenmode source mapping" to the Method ButtonGroup
  (`:137-170`).
- **Spectral-prior sub-panel** — Log (2026) / Flat (none) / Power (1/f),
  default Log. The existing Measure panel (Current density / dSPM / sLORETA) is
  reused for the eigenmode measure (`mne`/`dspm`/`sloreta`).
- **Coefficients checkbox** — "Also save eigenmode coefficients", default off.
- **Whitening checkbox** — "Apply noise whitening (recommended)", **default on**.
  When unchecked, the solver receives identity `iW` (pure sensor-space solve);
  SSP, bad channels, and average reference still apply.
- **Mutually exclusive gating** — `CreatePanel` gains an `isEigenmode`
  argument. When true: select the eigenmode radio, disable the standard
  methods, and hide irrelevant source-orientation options. When false: hide the
  eigenmode radio. Mirrors the existing `HeadModelType` gating at `:155-167`.
- **`GetPanelContents` / `GetSelectedMethod`** (`:626-747`) — return
  `InverseMethod='eigenmode'`, `InverseMeasure∈{mne,dspm,sloreta}`, plus new
  `EigenmodePrior∈{log,flat,power}` and `SaveCoefficients` (bool) fields.

### `process_inverse_2018.m` (shared core)

- **Detect `isEigenmode`** in the study loop (`:184-227`) from
  `HeadModel(iHeadModel)` and pass it to `gui_show_dialog(@panel_inverse_2018, …)`
  (`:245`).
- **New `case 'eigenmode'`** in the solver switch (`:665`): build the whitener
  with `bst_noise_whitener(OPTIONS.NoiseCovMat, OPTIONS.ChannelTypes,
  OPTIONS.NoiseMethod, OPTIONS.NoiseReg)` and reuse the SSP `Proj` that
  `Compute()` already constructs (`:555`); load `L_tilde = HM.Gain(GoodChannel,:)`
  and `lambdas = HM.Eigenvalues`; call
  `bst_inverse_eigenmodes('SolvePure', L_tilde, lambdas, iW, Proj, Measure,
  EigenmodePrior, Alpha, SNR, Unreg)`; then `bst_eigenmode_reconstruct` to set
  `Results.ImagingKernel = Φ·M̃`. Set `OPTIONS.FunctionName = ['eigenmode_' Measure]`
  so the saved node label matches the existing convention.
- **Coefficients node** — after the generic results node is saved (`:764-810`),
  if `SaveCoefficients` and the input is non-raw, build the `matrix_eigencoeffs`
  node (per-mode time series `θ = M̃·d`, rows labeled `Mode k (lam=…)`), exactly
  as `process_eigenmodes_inverse` does today.

### `process_eigenmodes_inverse.m` (batch wrapper)

- `GetDescription` keeps Index 339 / SubGroup "Sources" and exposes the
  eigenmode options (method / prior / SNR / nModes / save-coefficients).
- `Run()` packs those into `OPTIONS` (with `InverseMethod='eigenmode'`,
  `DisplayMessages=0`) and calls `process_inverse_2018('Compute', iStudies,
  iDatas, OPTIONS)`. The ~120 lines of solve/reconstruct/node-building are
  removed in favor of the shared core.

## Whitening

The eigenmode solver currently ships a simple single-global-block whitener
(`bst_inverse_eigenmodes.m:85-92`, equivalent to the legacy `bst_whitener.m`).
That is correct for single-modality data but biased for combined modalities
(MEG mag+grad, MEG+EEG), where the regularization is dominated by the
highest-variance modality. The correct per-modality whitener exists **only**
inside `bst_inverse_linear_2018.m` — there is no standalone function for it.

**Decision:** copy the per-modality whitener **verbatim** into a new shared
helper `bst_noise_whitener.m`:

- Source: the block at `bst_inverse_linear_2018.m:347-449` (per-modality
  truncate/regularize loop, zeroing cross-modality covariance) plus the local
  subfunctions `truncate_and_regularize_covariance` (`:1098`) and
  `cov1para_local` (`:1236`, Ledoit-Wolf shrinkage). Copied literally — no
  reinterpretation of the math.
- **`bst_inverse_linear_2018.m` is NOT modified** in this change. The standard
  MNE/dSPM/sLORETA/LCMV/dipole solves stay byte-for-byte identical, so the
  proven path cannot regress.
- The `'eigenmode'` case calls `bst_noise_whitener` to produce `iW`, reuses
  `Compute()`'s existing SSP `Proj`, and feeds both into
  `bst_inverse_eigenmodes('SolvePure', …)`. No double-whitening: the solver's
  internal whitener/projector are bypassed on this path.

Because the dialog reuses `panel_inverse_2018`, the existing **noise-covariance
regularization** controls (`NoiseMethod` = reg/shrink/diag, `NoiseReg`) apply to
the eigenmode method unchanged — same units balancing as every other method.

**Default vs. option:** whitening is **on by default** — the eigenmode case
builds `iW` with `bst_noise_whitener` so the full sensor-data procedure (channel
selection, SSP, bad channels, average reference, per-modality whitener) matches
every other Brainstorm inverse. The dialog's "Apply noise whitening" checkbox
(default checked) lets the user opt out, in which case `iW = eye(nGoodCh)` is
passed to `'SolvePure'` — a pure sensor-space solve. SSP/bad-channels/avg-ref are
folded into `L_tilde` and the data path regardless, so only the whitening step is
removed. The batch wrapper exposes the same toggle (default on).

**Tech debt (flagged, not done now):** this leaves two copies of the whitener
logic. A later change can refactor `bst_inverse_linear_2018` to call
`bst_noise_whitener`, but only with full regression testing of the standard
methods — out of scope here.

Source-orientation options (`fixed`/`loose`/`free`) are irrelevant to eigenmode
and are hidden/ignored in the dialog and the solver branch.

## Edge cases

- **Raw input** — cortex sources kernel is produced (Time left empty, viewer
  fetches from the data file); the coefficients node is skipped with a warning
  (coefficients require imported data), matching current behavior.
- **No noise covariance** — fall back to identity whitening with a warning
  (existing `bst_inverse_eigenmodes` behavior).
- **Non-eigenmode head model selected with eigenmode forced** — the solver
  already errors clearly (`:60-64`); the gating prevents this in the GUI.
- **`nModes` cap** — 0 means all modes in the head model (unchanged semantics).

## Testing

- **Pure** — extend `dev/tests/test_inverse_eigenmodes_pure.m` to cover the
  reconstruction helper (`Φ·M̃` shape and values) independent of the DB.
- **Whitener parity** — assert `bst_noise_whitener` returns the same `iW` as the
  inline whitener in `bst_inverse_linear_2018` for representative single- and
  multi-modality covariances (guards the verbatim copy against drift).
- **e2e** — extend `dev/tests/test_eigenmodes_inverse_e2e.m` to drive the new
  path via `process_inverse_2018('Compute', iStudies, iDatas, OPTIONS)` with
  `InverseMethod='eigenmode'`; assert the cortex results node is created and, with
  the checkbox on, the `matrix_eigencoeffs` node appears.
- **Manual** — right-click a recording whose active head model is an eigenmode
  leadfield → confirm only the eigenmode method is enabled, run it, and verify
  the cortex source map and scout time series display normally.

## Files touched

| File | Change |
|------|--------|
| `toolbox/inverse/panel_inverse_2018.m` | Eigenmode radio, prior sub-panel, coefficients checkbox, `isEigenmode` gating, `GetPanelContents`/`GetSelectedMethod` fields |
| `toolbox/process/functions/process_inverse_2018.m` | `isEigenmode` detection, dialog arg, `case 'eigenmode'` (build `iW` + reuse `Proj`, call `SolvePure`, reconstruct), coefficients node |
| `toolbox/inverse/bst_noise_whitener.m` (new) | Verbatim copy of the per-modality whitener from `bst_inverse_linear_2018.m:347-449` + subfunctions `:1098`,`:1236` |
| `toolbox/inverse/bst_inverse_eigenmodes.m` | Unchanged math; driven via `'SolvePure'` with injected `iW`/`Proj` |
| `toolbox/inverse/bst_inverse_linear_2018.m` | **Not modified** (proven solve kept intact) |
| `toolbox/inverse/bst_eigenmode_reconstruct.m` (new, standalone) | Load Φ, return `Φ·M̃` and optional `M̃` |
| `toolbox/process/functions/process_eigenmodes_inverse.m` | Reduce to thin wrapper over `Compute()` |
| `dev/tests/test_inverse_eigenmodes_pure.m`, `dev/tests/test_eigenmodes_inverse_e2e.m` | Extend coverage (incl. `bst_noise_whitener` parity vs. the 2018 inline whitener) |

## Fix: across-hemisphere mode selection (one-hemisphere bug)

`tess_eigenmodes` solves each hemisphere as a separate connected component and
stores the eigenvectors **grouped by component** (`[LH modes…, RH modes…]`,
not globally sorted by eigenvalue). The forward composer
(`bst_eigenmode_leadfield`) and the inverse reconstruction
(`bst_eigenmode_reconstruct`) originally took the **first K columns**
(`Phi(:,1:K)`), so when "Number of modes" `K` was ≤ the per-hemisphere count
(e.g. `K=300` with 300 modes per hemisphere), the leadfield and reconstruction
spanned **only the left hemisphere** — the right hemisphere received exactly zero
sources.

**Fix:** select the **K globally-lowest-eigenvalue** modes across all components
(both hemispheres), recorded in `headmodel.ModeIndices` so the leadfield and the
reconstruction use the exact same columns in the same order. "Number of modes"
now means a **total** budget (≈ K/2 per hemisphere for near-symmetric
hemispheres). Files: `bst_eigenmode_leadfield.m` (selection + `ModeIndices`),
`bst_eigenmode_reconstruct.m` (optional `ModeIndices` arg; falls back to
`Phi(:,1:K)` for old head models), `process_inverse_2018.m` (loads/passes
`ModeIndices`). Regression: `dev/tests/test_eigenmode_hemisphere_pure.m`.

**ACTION REQUIRED:** existing eigenmode head models have the one-hemisphere
selection baked into their `Gain` and lack `ModeIndices`; they must be
**recomputed** ("Compute head model" → Cortex surface harmonics) to get
both-hemisphere coverage.
