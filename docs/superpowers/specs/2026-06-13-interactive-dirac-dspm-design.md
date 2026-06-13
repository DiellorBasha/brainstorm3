# Interactive Dirac-dSPM — replace Eigen-dSPM in the GUI

**Date:** 2026-06-13
**Author:** Diellor Basha (with Claude)
**Status:** Approved (design forks settled interactively)

## Problem

The Dirac eigenmode inverse (`bst_inverse_dirac`, wrapped by the pipeline process
`process_inverse_dirac`) is fully implemented and validated, but it is reachable
**only** from the Process tab (*Sources → "Compute sources: Dirac eigenmodes"*).
The interactive right-click **"Compute sources"** flow instead exposes the older
scalar **Eigen-dSPM** method (the "Eigenmode source mapping" radio inside
`panel_inverse_2018`). We are deprecating Eigen-dSPM and want the Dirac method to
be the interactive source-mapping option.

## Decisions (settled)

1. **Wiring:** Dirac gets its **own interactive flow** — a new
   `process_inverse_dirac('ComputeInteractive', ...)` entry plus a small dedicated
   options dialog, reached from a new tree menu item. We do **not** fold Dirac into
   the shared `panel_inverse_2018` (its tau/nModes options don't fit there, and
   reusing the proven `process_inverse_dirac.Run` is lower-risk).
2. **Deprecation depth:** **Remove Eigen-dSPM from the GUI, keep the code dormant.**
   The `bst_inverse_eigenmodes` solver and the `'eigenmode'` branch in
   `process_inverse_2018` stay on disk (still used by the accuracy-benchmark
   harness/tests); only the GUI entry point is removed.

## Architecture

Three units, each with a single clear purpose:

### 1. `toolbox/inverse/panel_inverse_dirac.m` (new)

A compact options panel mirroring the option set already defined in
`process_inverse_dirac.GetDescription`:

- **Measure** — radio: *Current density map* / *dSPM* / *sLORETA*
  (values `amplitude` / `dspm2018` / `sloreta`; default `dspm2018`).
- **Signal-to-noise ratio (SNR)** — value (default 3).
- **Dirac modes per hemisphere** — value (default 400).
- **Dirac τ (intrinsic↔extrinsic mix)** — value (default 0.5).
- **Noise covariance regularization (fraction)** — value (default 0.1).
- **Sensor types or names** — text (default `MEG`).

Standard Brainstorm panel contract: `CreatePanel()` builds the Swing form;
`GetPanelContents()` returns a struct `s` with fields
`Measure, SNR, nModes, Tau, NoiseReg, SensorTypes`. Modeled on the structure of
`panel_inverse_2018` but minimal (no source-orientation/output sections — the
Dirac source space is always unconstrained 3-vector cortex).

### 2. `process_inverse_dirac('ComputeInteractive', bstNodes)` (new subfunction)

- Resolve `bstNodes` → channel studies, mirroring `panel_protocols('TreeInverse')`
  node resolution: `data`/`rawdata` → `tree_dependencies`; otherwise
  `tree_channel_studies(bstNodes, 'NoIntra')`. Bail out cleanly on empty / `-10`.
- Pop the dialog via `gui_show_dialog('Compute sources', @panel_inverse_dirac, ...)`.
  Cancel → return `{}`.
- Map the returned struct into a synthetic `sProcess` (clone of
  `GetDescription()` with `options.*.Value` overwritten) and a minimal `sInputs`
  array (`sInputs(k).iStudy = iStudies(k)`; `Run` builds shared kernels per channel
  study and never reads the data samples, so no file payload is needed).
- Call the existing `Run(sProcess, sInputs)` unchanged. Surface errors via
  `bst_error` if no outputs and an error message was produced.

`Run` itself is **not modified** — the interactive path feeds it the same
`sProcess`/`sInputs` shape it already consumes from the pipeline.

### 3. `toolbox/tree/tree_callbacks.m` (edits)

- **Add** a menu item *"Compute sources: Dirac eigenmodes"* in
  `fcnPopupComputeSources()` (line ~2825), right under the existing
  *"Compute sources [2018]"* item. This single edit covers all six call sites
  (data, condition, subject, etc.). Callback:
  `@(h,ev)bst_call(@process_inverse_dirac, 'ComputeInteractive', bstNodes)`.
- **Add** the same item at the head-model-node menu (line ~1399), via a thin
  wrapper that first double-clicks the node to set it active (mirroring
  `selectHeadmodelAndComputeSources`) then calls `ComputeInteractive`.

### 4. `toolbox/inverse/panel_inverse_2018.m` (edit — remove Eigen-dSPM from GUI)

Make the **"Eigenmode source mapping"** radio permanently hidden and disabled,
regardless of the `isEigenmode` argument. Concretely, replace the conditional
visibility block (lines ~173–188) with an unconditional
`jRadioMethodEig.setVisible(0); jRadioMethodEig.setEnabled(0);`, and drop the
`isEigenmode`-forcing branch that disabled the standard methods. All other
eigenmode plumbing in the panel (`jPanelEig`, getters, `ctrl` fields,
`GetSelectedMethod`/`GetMethodComment` `'eigenmode'` cases) stays in place but is
unreachable — the sub-panel is only shown when the now-never-selected radio is
selected. This keeps the code dormant per the deprecation decision with a minimal
diff.

## Data flow

```
right-click recording/condition/headmodel
  └─ "Compute sources: Dirac eigenmodes"
       └─ process_inverse_dirac('ComputeInteractive', bstNodes)
            ├─ resolve nodes → iStudies (channel studies)
            ├─ gui_show_dialog(@panel_inverse_dirac) → options struct
            ├─ build sProcess (options) + sInputs (.iStudy)
            └─ process_inverse_dirac('Run', sProcess, sInputs)
                 └─ bst_inverse_dirac(HeadModel, OPT)  [unchanged]
                      → shared results_DiracEig_KERNEL_*.mat + tree refresh
```

## Error handling

- No channel / no head model / no noise covariance in a study → per-study
  `bst_report('Error', ...)` and skip (already handled in `Run`).
- Head model not a surface / not unconstrained / scalar eigenmode model →
  existing `Run` guards produce a clear error.
- Dialog cancelled → `ComputeInteractive` returns `{}` silently.
- Empty/invalid node selection → early return / `bst_error`.

## Testing / validation

- **Static:** `checkcode` (M-Lint) on the three touched/new files.
- **Interactive smoke test** (MATLAB MCP, headless Brainstorm server):
  on a tutorial protocol with an unconstrained surface head model + noise cov,
  call `process_inverse_dirac('ComputeInteractive', <data nodes>)` and confirm a
  `results_DiracEig_KERNEL_*.mat` shared kernel appears and loads.
- **Regression:** confirm the standard `panel_inverse_2018` panel still opens and
  offers MN/dSPM/sLORETA/beamformer/dipole/MEM with the Eigenmode radio gone, and
  that the dormant `process_inverse_2018` `'eigenmode'` branch + benchmark harness
  are untouched.

## Out of scope (YAGNI)

- No new head-model type or storage format — uses the existing unconstrained
  surface head model; the Dirac basis is found-or-created internally.
- No removal/retirement of `bst_inverse_eigenmodes` or the process-level eigenmode
  branch (kept dormant for benchmarks).
- No changes to `bst_inverse_dirac` or `process_inverse_dirac.Run`.
```
