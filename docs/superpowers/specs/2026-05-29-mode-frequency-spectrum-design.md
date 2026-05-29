# Mode-Frequency Spectrum — Design Spec

**Date:** 2026-05-29
**Author:** Diellor Basha (with Claude)
**Branch:** `feature/eigenmode-spectrum`
**Status:** Approved design → ready for implementation plan

---

## Goal

Compute the **mode-frequency spectrum**: for each Laplace–Beltrami eigenmode, the
temporal-frequency content of its **mode coefficient** over time. We do this by
running Brainstorm's existing spectral engine (`bst_timefreq`, methods PSD and FFT)
on **mode coefficient time series** that are built in memory on the fly and never
saved to disk.

This is the temporal-frequency counterpart to the (already shipped) spatial
eigenspectrum: the spatial spectrum asks *which spatial patterns are present*; the
mode-frequency spectrum asks *how each spatial pattern oscillates in time*.

## Architecture

Two new batch processes under **Sources → Frequency**, mirroring the option panels
and UI of the existing `process_psd` / `process_fft`, both delegating to one shared
engine. The engine builds the mode coefficient time series in memory (via a **mode
kernel** that folds the eigenmode projection into the imaging kernel), then hands
the matrix to `bst_timefreq` exactly the way the existing scout/cluster path does
(`DataToProcess = {matrix}` + `RowNames`/`TimeVector`/`nComponents`/`SurfaceFile`).
The output is a standard timefreq/spectrum file that opens in Brainstorm's existing
`figure_spectrum` line viewer with no custom display code this round.

**Self-contained — no edits to shared/core files** (`bst_timefreq`, `bst_psd`,
`process_psd`, `process_fft` are reused unchanged). This keeps the upstream-merge
surface minimal, consistent with the codebase's "add a new `process_eigenmodes_*`"
pattern.

### Why this reuse is mathematically exact (verified in code)

For a *source* input, `bst_timefreq` already does read-window → (kernel) → Welch/FFT
(`bst_timefreq.m:346‑380`, calling `bst_psd` at 532/543). Crucially, `bst_psd`
applies the imaging kernel **after** the FFT, in the frequency domain
(`bst_psd.m:157‑158`: `TFwin = ImagingKernel * TFwin`, *then* power at 163). This is
valid because the FFT is linear: `FFT(K·sensors) = K·FFT(sensors)`.

The mode coefficient is the same kind of linear map of the sensors,
`coeff(t) = (Φ'·M·K)·F(t) = ModeKernel·F(t)`, so `FFT(coeff) = ModeKernel·FFT(sensors)`.
Therefore **building the mode coefficient time series first and running PSD/FFT on it
yields the identical result** to substituting `ModeKernel` for `ImagingKernel` inside
`bst_psd`. Our self-contained path (build `coeff`, feed the in-memory matrix to
`bst_timefreq`) is exact, not an approximation. We chose it over a guarded
`bst_timefreq` edit specifically to avoid touching core files.

## Tech Stack

MATLAB; Brainstorm process framework (`macro_method` dispatch); `bst_timefreq`
(timefreq engine); `bst_eigenmodes_project` / `tess_laplacian` / `in_tess_eigenmodes`
(eigenmode math); MATLAB MCP for testing.

---

## Terminology (used everywhere — code, GUI labels, file comments)

| Math | Plain name | Meaning |
|------|------------|---------|
| `u(t)` | **cortical map** | source value at every vertex at one instant |
| `φ_k` | **eigenmode** / **mode** | one fixed spatial pattern on the cortex |
| `λ_k` | **spatial frequency** | the eigenvalue (smooth → low, wiggly → high) |
| `c_k(t)` | **mode coefficient** | how strongly mode `k` is present in `u(t)` |
| `c_k(t)` over all `t` | **mode coefficient time series** | `[nModes × nTime]`, one line per mode |
| PSD/FFT of that | **mode frequency spectrum** (per mode) | power/amplitude vs temporal frequency |
| stacked over modes | **mode-frequency spectrum** | the joint object (λ axis × ω axis) |
| `Φ'·M·Kernel` | **mode kernel** | `[nModes × nChannels]`, sensors → mode coefficients |

---

## Scope

**In scope (v1):**
- Input: surface source **results** files (kernel-linked, including raw-linked; or
  full `ImageGridAmp`).
- **Constrained orientation only** (`nComponents == 1`).
- Default **300 modes** (user-settable "Number of modes" option).
- Methods **PSD (Welch)** and **FFT**.
- Process a **bounded time window** (default: full time for imported/epoched data;
  a user-selected window for raw — never the whole recording).
- Visualization: the **existing `figure_spectrum` lines view** (one spectrum line
  per mode), available for free on the output file.

**Out of scope (deferred):**
- The joint 2-D **heatmap** viewer (spatial freq × temporal freq, color = power) —
  next round.
- Morlet / Hilbert methods (3-D mode×time×freq output; need phase handling).
- Unconstrained sources (needs the planned tangent-frame / vector-heat work).
- Full-raw streaming Welch (v1 processes a bounded window in one read).

---

## File Structure

**New files:**

- `toolbox/math/bst_eigenmodes_modekernel.m` — *pure helper.*
  `ModeKernel = bst_eigenmodes_modekernel(Eigenmodes, MassMatrix, ImagingKernel, nModes)`
  returns `Φ(:,1:nModes)' · M · ImagingKernel` → `[nModes × nChannels]`. With an
  empty `ImagingKernel`, returns the projector `Φ(:,1:nModes)' · M`
  (`[nModes × nVertices]`) to apply directly to a stored `ImageGridAmp`.
  No DB access, no I/O — independently unit-testable.

- `toolbox/process/functions/process_eigenmodes_freq.m` — *shared engine, not listed
  in the menu (no `GetDescription`).* `Run(sProcess, sInputs, Method)`:
  validates inputs, builds the mode coefficient time series in memory, assembles
  `tfOPTIONS`, calls `bst_timefreq`, tags the output. `Method ∈ {'psd','fft'}`.

- `toolbox/process/functions/process_eigenmodes_psd.m` — *listed process.*
  `GetDescription` mirrors `process_psd` options (time window, window length,
  overlap, units, win_std) plus the shared "Number of modes" option.
  `Run` → `process_eigenmodes_freq('Run', sProcess, sInputs, 'psd')`.

- `toolbox/process/functions/process_eigenmodes_fft.m` — *listed process.*
  `GetDescription` mirrors `process_fft` options (time window) plus "Number of
  modes". `Run` → `process_eigenmodes_freq('Run', sProcess, sInputs, 'fft')`.

**Reused unchanged:** `bst_eigenmodes_project`, `in_tess_eigenmodes`,
`tess_laplacian`, `bst_timefreq`, and the saved `Eigenmodes.MassMatrix` / `MassType`
(skip recomputing the mass matrix when present).

---

## Data Flow (per input result file, inside `process_eigenmodes_freq`)

1. **Load metadata only** —
   `in_bst_results(File, 0, 'ImagingKernel','ImageGridAmp','GoodChannel',
   'nComponents','DataFile','Time','SurfaceFile','HeadModelType','Atlas')`.
   `isLoadFull = 0` so the kernel is never expanded here.

2. **Validate** (`bst_report('Error', ...)` + `continue` on failure):
   surface head model; no atlas; `nComponents == 1`; eigenmodes present on
   `SurfaceFile` (`in_tess_eigenmodes`); eigenmode vertex count == source vertex
   count.

3. **Mass matrix** — reuse `Eigenmodes.MassMatrix` if present; else
   `[~, M] = tess_laplacian(V, F, 'MassType', Eigenmodes.MassType)`.

4. **Time window** — from the `timewindow` option; default whole time for
   imported/epoched, user-selected for raw.

5. **Build mode coefficient time series (in memory, capped at nModes):**
   - *Kernel result* (`ImagingKernel` non-empty): bounded **windowed** read of the
     sensor data `F` for the selected window (raw-safe — never the whole recording);
     `ModeKernel = bst_eigenmodes_modekernel(Eig, M, ImagingKernel, nModes)`;
     `coeff = ModeKernel · F(GoodChannel, window)`.
   - *Full result* (`ImageGridAmp` present, no kernel): load `ImageGridAmp` for the
     window; `P = bst_eigenmodes_modekernel(Eig, M, [], nModes)` (= `Φ'·M`);
     `coeff = P · ImageGridAmp(:, window)`.
   - `coeff` is `[nModes × nTimeWindow]`.

6. **Run the engine:**
   - `tfOPTIONS = bst_timefreq();` (defaults), then set `Method`, and map options
     (PSD: `WinLength`, `WinOverlap`, `PowerUnits`, `WinFunc`, `TimeWindow`;
     FFT: `TimeWindow`).
   - `DataToProcess = {coeff}`; `tfOPTIONS.TimeVector = windowTimeVector`;
     `tfOPTIONS.ListFiles = {File}`; `tfOPTIONS.nComponents = 1`;
     `tfOPTIONS.SurfaceFile = {SurfaceFile}`;
     `tfOPTIONS.RowNames = { {'Mode 1 (λ=…)', …, 'Mode K (λ=…)'} }`;
     output as a single file (`Output = 'all'`).
   - Mark the result as a generic named-signal spectrum (`DataType = 'matrix'`) so
     rows are treated as named modes (not sources) in the viewers.
   - `[OutputFiles, Messages, isError] = bst_timefreq(DataToProcess, tfOPTIONS);`
     surface messages via `bst_report`.

7. **History / comment** — record method, #modes, mass type, eigenvalue range,
   source file, and time window on the output file.

---

## Options / UI

- **PSD** (`process_eigenmodes_psd`): Time window; Window length (s); Window overlap
  (%); Units (physical/normalized); Save std across windows — identical wording to
  `process_psd`. Plus **Number of modes** (default 300).
- **FFT** (`process_eigenmodes_fft`): Time window. Plus **Number of modes** (300).
- The "Number of modes" option caps the projection to the first K modes; if fewer
  modes exist than requested, use all available and note it in the report.

## Output & Visualization

- A standard timefreq/spectrum file, rows = modes (labeled `Mode k (λ=…)`),
  x = temporal frequency. PSD → power vs frequency; FFT → amplitude vs frequency.
- Opens via Brainstorm's existing **Power spectrum** display (`view_spectrum`,
  `figure_spectrum`): one line per mode. No custom viewer this round.
- The joint 2-D heatmap is explicitly deferred to a follow-up round.

---

## Testing

- **Pure unit test** — `dev/tests/test_bst_eigenmodes_modekernel.m`: on a tiny
  synthetic mesh + random kernel, assert `ModeKernel · F == Φ'·M·(Kernel·F)` to
  numerical tolerance (the folded route equals project-after-source); assert the
  no-kernel branch returns `Φ'·M` and that mode-capping selects the first K rows.
  Prints `ALL TESTS PASSED`.
- **Integration test** — `dev/tests/test_process_eigenmodes_freq.m`: build a small
  synthetic source result + eigenmodes in a temp protocol, run PSD and FFT, assert
  the output timefreq file has shape `[nModes × 1 × nFreq]`, correct RowNames
  (mode labels), a sane ascending positive frequency vector, and finite power.
  Prints `ALL TESTS PASSED`.
- **Manual GUI validation** by the user afterward (the usual loop): compute
  eigenmodes → run "Mode-frequency spectrum (PSD/FFT)" on a source result →
  open the power-spectrum lines view.

---

## Risks / Edge Cases

- **Raw windowed read:** the engine reads sensor data for the selected window the
  same way `bst_timefreq` does — load the raw `sFile` descriptor, compute
  `SamplesBounds` from the time window, single `in_fread` call (replicating the
  ~5-line `ReadRawRecordings` pattern in `bst_timefreq.m:960`). Memory/cost are
  identical to a standard PSD/FFT run.
- **Raw default time window:** mirror standard PSD/FFT — the `timewindow` option is
  empty by default and an empty window means the whole file is read in one bounded
  `in_fread`, exactly like `process_psd`/`process_fft`. Imported/epoched data uses
  the full epoch. (No special guard; behavior matches the existing frequency tools.)
- **Bad-segment exclusion (raw):** because we feed the mode coefficient time series
  to `bst_timefreq` via the in-memory matrix path (rather than letting it read the
  raw file itself), automatic bad-segment exclusion is **not** inherited in v1. This
  matches the "light" scope; a later round can pass `BadSegments` through if needed.
- **Mode capping vs available modes:** request 300 but a surface may have fewer;
  clamp and report.
- **Mass-matrix provenance:** prefer the saved `Eigenmodes.MassMatrix`; only
  recompute when absent, and use the recorded `MassType` so the projection matches
  how the eigenmodes were built.
- **DataType for in-memory matrices:** ensure `bst_timefreq` tags the output so it
  opens in the spectrum viewer with per-mode rows (mirror the scout/cluster path).
- **Frequency resolution / window length** come from the standard PSD panel; no new
  validation beyond what `process_psd` already does.
