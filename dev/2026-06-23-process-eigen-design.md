# process_eigen orchestrator + process_eigenspectrum — Design

**Date:** 2026-06-23
**Author:** Diellor Basha
**Status:** Approved design → implementation plan next

## Goal

Expose the eigen-domain analysis engine (`bst_eigen`) in the Brainstorm Process
panel so eigenspectra can be batch-computed across many recordings/subjects and
stored as `timefreq_` nodes in the tree — the spatial-spectral twin of
`process_psd` (Welch PSD).

## Background / why this is small

The compute and database layers already exist and mirror the timefreq stack:

| Temporal (time-frequency) | Spatial-spectral (eigen) | Status |
|---|---|---|
| `bst_psd` (Welch windowed FFT → power) | `bst_eigenspectrum` (windowed manifold FT `C = Φᵀ(B·U)` → mode power) | exists |
| `bst_timefreq` (read → dispatch → save timefreq node) | `bst_eigen` (read → dispatch `spectrum`/`filter`/`wavelet` → save node) | exists |
| `process_timefreq` (shared Run + dispatch by caller name) | **`process_eigen`** | **build** |
| `process_psd` / `process_hilbert` / `process_fft` … (thin wrappers) | **`process_eigenspectrum`** (+ future filter/wavelet) | **build** |

`bst_eigen` already owns the full read → dispatch → save pipeline for the
`spectrum` method (`toolbox/eigen/bst_eigen.m`):

- reads `results` files via `in_bst_results(InitFile, 1, 'ImageGridAmp','Time','SurfaceFile','nComponents')`
  — the `1` (full load) means a kernel-only result is expanded to a source map
  automatically, so kernel-applied-on-the-fly is already handled;
- restricts to `OPTIONS.TimeWindow`;
- computes the windowed eigenspectrum **per hemisphere** via `ComputeEigenspectrum`
  → `bst_eigenspectrum` (Welch-style time windowing: `WinLength`/`WinOverlap`/`WinFunc`,
  identical scheme to `bst_psd`);
- writes a `timefreq_eigenspectrum` node (`Freqs = sqrt(Lambda)`, `Method='spectrum'`)
  via `bst_save` + `db_add_data` when `OPTIONS.iTargetStudy` is set, or returns the
  struct when it is `'NoSave'`.

There is exactly **one functional gap** inside `bst_eigen` for the panel use case:
`GetEigenBasis` (bst_eigen.m ~lines 305–309) still **errors** on the implicit
resolution path instead of calling the resolver that now exists in `bst_get`:

```matlab
% bst_get('EigenFileForSurface', SurfaceFile, Variant, nModes, Tau)  (bst_get.m:1381–1411)
```

`bst_get('EigenFileForSurface', ...)` already loops a surface's `Eigen` cache,
filters by `Variant` (and `Tau` for Dirac-type), and returns the best-fit node
(smallest `nModes ≥` requested). The SurfaceFile is implicit from the results
file; the **Variant** is the one thing the caller must supply, because a single
surface can host several `eigen_` nodes at once (Laplace-Beltrami, Connection
Laplacian, Dirac, Dirac-Face, Hodge-Face, and/or different `nModes`).

## Decisions (locked)

1. **`process_eigen` is a pure orchestrator** — no `GetDescription`, never shown
   in the panel. It is only the shared `Run()` that siblings delegate to, mirroring
   how `process_psd`/`process_hilbert` delegate to `process_timefreq('Run',…)`.
   (Unlike `process_timefreq`, which also doubles as the Morlet member;
   `process_eigen` deliberately does not double as a member.)
2. **Eigen file resolved implicitly from `SurfaceFile`**, via
   `bst_get('EigenFileForSurface', SurfaceFile, Variant, nModes, Tau)`.
3. **Variant chosen by a dropdown**, default **Laplace-Beltrami**. `nModes` and
   `Tau` are NOT exposed; resolution uses defaults (largest available basis for
   the variant; `Tau` only consulted for Dirac-type).
4. **v1 scope = orchestrator + spectrum only.** `process_eigenfilter` /
   `process_eigenwavelet` are *named* in the dispatch map but not created here.
5. **Welch-style time windowing is a first-class option** (`win_length`,
   `win_overlap`, `win_std`), exactly like `process_psd`. Leaving `win_length`
   empty computes a single full-window spectrum.
6. **Spectrum node shape is unchanged** — reuse `bst_eigen`'s existing
   `BuildSpectrumTimefreq` output (`Freqs = sqrt(Lambda)`).

## Architecture

```
process_eigenspectrum (registered)  ──┐
   [future: process_eigenfilter,      ├─→ process_eigen('Run',…) ──→ bst_eigen(files, OPTIONS)
            process_eigenwavelet]   ──┘   (pure: dispatch by         (existing compute + save)
                                           func2str(sProcess.Function))
```

## Components

### 1. `toolbox/eigen/bst_eigen.m` — close the resolver gap

- Add `Def_OPTIONS.Tau = []`.
- In `GetEigenBasis`, replace the `error('bst_eigen:AutoResolveTODO', …)` block
  with:
  - call `bst_get('EigenFileForSurface', SurfaceFile, OPTIONS.Variant, OPTIONS.nModes, OPTIONS.Tau)`;
  - if it returns empty, raise a clear error:
    `"No '<Variant>' eigenbasis found for this surface — compute it first."`;
  - otherwise read the relative path from the returned cache entry
    (`sSubject.Surface(iSurface).Eigen(iEigen).FileName`) and continue to
    `in_bst_eigen` / `in_bst_operator` as today.
- No change to the dispatch, compute, or save paths.

### 2. `toolbox/process/functions/process_eigen.m` (NEW) — shared Run + dispatch

- **No `GetDescription`** (pure orchestrator).
- `Run(sProcess, sInputs)`:
  1. Map `func2str(sProcess.Function)` → Method:
     `process_eigenspectrum`→`'spectrum'`, `process_eigenfilter`→`'filter'`,
     `process_eigenwavelet`→`'wavelet'`. Unknown caller → error.
  2. Build `OPTIONS = bst_eigen()` defaults, then fill from `sProcess.options`
     (Method, Variant, Measure, WinLength, WinOverlap, WinFunc, TimeWindow).
     Leave `iTargetStudy = []` so each output lands in its own input's study.
  3. `OutputFiles = bst_eigen({sInputs.FileName}, OPTIONS)` — `bst_eigen` already
     loops inputs and resolves each input's study via `bst_get('AnyFile', …)`.

### 3. `toolbox/process/functions/process_eigenspectrum.m` (NEW) — registered sibling

- `GetDescription`:
  - Comment: `'Eigenspectrum (spatial PSD)'`
  - Category `'Custom'`, SubGroup `'Frequency'`, an unused Index in that band.
  - **InputTypes `{'results'}`**, OutputTypes `{'timefreq'}`.
  - Options:
    - `variant` — dropdown {Laplace-Beltrami (default), Connection Laplacian,
      Dirac, Dirac-Face, Hodge-Face}
    - `measure` — radio {power (default), magnitude}
    - `timewindow` — time window
    - `win_length` — value, seconds (default e.g. 1 s; empty → single window)
    - `win_overlap` — value, percent (default 50)
    - `win_std` — checkbox "Save standard deviation across windows"
      (maps `WinFunc = 'mean'` vs `'mean+std'`), mirroring `process_psd`.
- `Run`: one line — `OutputFiles = process_eigen('Run', sProcess, sInputs)`.

## Data flow (one input)

`results` file → `bst_eigen` reads `ImageGridAmp` (kernel auto-applied) +
`SurfaceFile` → `GetEigenBasis` resolves the `eigen_` node from
`(SurfaceFile, Variant)` → `ComputeEigenspectrum` per hemisphere (Welch windowing)
→ `timefreq_eigenspectrum` node (`Freqs = sqrt(λ)`, `Method='spectrum'`), parented
to the input. Batch = the Process panel looping inputs/subjects — no extra code.

## Error handling

- No matching `eigen_` node for `(surface, variant)` → clear "compute the basis
  first" error; reported per input, does not abort the rest of the batch.
- Volume head model / empty `SurfaceFile` → skipped with a message (eigenspectrum
  needs a surface basis).
- Vertex-count / variant-layout mismatch → surfaced from `ComputeEigenspectrum`.

## Testing (MATLAB MCP, live S01 protocol)

1. **Resolver** — `bst_get('EigenFileForSurface', surf, 'Laplace-Beltrami')`
   returns the LBO node for the input's surface.
2. **Dispatch** — caller-name → method map (`process_eigenspectrum`→`'spectrum'`),
   unknown caller errors.
3. **End-to-end implicit** — run `process_eigenspectrum` on a real `results` file
   with **no** EigenFile passed; assert a `timefreq_eigenspectrum` node is created,
   `Freqs == sqrt(λ)`, `Method == 'spectrum'`, TF shape matches the mode count.
4. **Windowed** — with `win_length` set, assert multiple Welch windows are used
   (`Nwin > 1`) and the result differs from the single-window spectrum; with
   `win_std`, assert `Std` is populated.
5. **Negative** — a variant with no basis on that surface → expected clear error,
   batch continues with the remaining inputs.

## Out of scope (fast-follows)

- `process_eigenfilter`, `process_eigenwavelet` (reuse this orchestrator).
- Exposing `nModes` / `Tau` as user options.
- Group/across-subject averaging of eigenspectra.
- Any change to `BuildSpectrumTimefreq` node shape or the spectrum viewer.
