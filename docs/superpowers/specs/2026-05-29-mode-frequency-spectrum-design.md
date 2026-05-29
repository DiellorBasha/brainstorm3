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
engine. The engine reproduces `bst_timefreq`'s source-input branch (read the
windowed sensor data + bad segments exactly as it does), then makes the **identical
`bst_psd` call** standard source PSD makes — substituting only the **mode kernel**
(`Φ'·M·ImagingKernel`) for the imaging kernel and mode labels for the row names — and
packages the result into a standard timefreq/spectrum file. That file opens in
Brainstorm's existing `figure_spectrum` line viewer with no custom display code this
round.

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
Therefore **substituting `ModeKernel` for `ImagingKernel` in the `bst_psd` call gives
the mode-coefficient spectrum directly** — exact, not an approximation.

We realize this by calling `bst_psd` ourselves with the same arguments
`bst_timefreq` uses (same windowed `F`, same `BadSegments`, same Welch options),
only swapping the kernel. This is **self-contained** (no edits to `bst_timefreq` /
`bst_psd` / `process_psd` / `process_fft`) **and** preserves Brainstorm's exact
behavior — crucially the bad-segment exclusion — so the mode-frequency spectrum is
guaranteed consistent with the source time-frequency spectrum on the same data.
(We deliberately do **not** use the simpler in-memory-matrix path through
`bst_timefreq`, because that path drops raw bad-segment handling and would diverge.)

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
  validates inputs, reads the windowed sensor data + bad segments the same way
  `bst_timefreq`'s source branch does, builds the mode kernel, calls `bst_psd`
  (substituting the mode kernel for the imaging kernel), and saves the result as a
  timefreq file. `Method ∈ {'psd','fft'}`.

- `toolbox/process/functions/process_eigenmodes_psd.m` — *listed process.*
  `GetDescription` mirrors `process_psd` options (time window, window length,
  overlap, units, win_std) plus the shared "Number of modes" option.
  `Run` → `process_eigenmodes_freq('Run', sProcess, sInputs, 'psd')`.

- `toolbox/process/functions/process_eigenmodes_fft.m` — *listed process.*
  `GetDescription` mirrors `process_fft` options (time window) plus "Number of
  modes". `Run` → `process_eigenmodes_freq('Run', sProcess, sInputs, 'fft')`.

**Reused unchanged (called, not edited):** `bst_psd` (the actual spectral compute +
bad-segment skipping + kernel application), `in_tess_eigenmodes`, `tess_laplacian`,
`in_fread` / `panel_record('GetBadSegments')` (the windowed-read + bad-segment
pattern mirrored from `bst_timefreq`'s `ReadRawRecordings`), and the saved
`Eigenmodes.MassMatrix` / `MassType` (skip recomputing the mass matrix when present).

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

5. **Read sensor data + bad segments for the window** (mirroring
   `bst_timefreq.m:346‑380` + `ReadRawRecordings`):
   - *Kernel result* (`ImagingKernel` non-empty): from the associated `DataFile`,
     load `sFile`/`ChannelMat`; compute `SamplesBounds` from the time window;
     `[F, TimeVector] = in_fread(...)` (raw — bounded read, never the whole file) or
     `F = sMat.F` (imported); restrict `F = F(GoodChannel, :)`; obtain
     `BadSegments` via `panel_record('GetBadSegments', ...)` exactly as the standard
     path does. Build `ModeKernel = bst_eigenmodes_modekernel(Eig, M, ImagingKernel, nModes)`.
   - *Full result* (`ImageGridAmp` present, no kernel): `F = ImageGridAmp(:, window)`
     (the sources themselves); `ModeKernel = bst_eigenmodes_modekernel(Eig, M, [], nModes)`
     (= `Φ'·M`); `BadSegments` per the standard imported-data path.

6. **Compute the spectrum via the identical `bst_psd` call** standard source PSD/FFT
   makes, swapping only the kernel:
   - PSD: `[TF, Freqs, Nwin] = bst_psd(F, sfreq, WinLength, WinOverlap, BadSegments, ModeKernel, WinFunc, PowerUnits, IsRelative);`
   - FFT: `[TF, Freqs] = bst_psd(F, sfreq, [], 0, BadSegments, ModeKernel, [], PowerUnits);`
   - Because `bst_psd` skips bad-segment windows *before* applying the kernel
     (`bst_psd.m:157`), bad-segment treatment is bit-for-bit identical to a source
     PSD. `TF` is `[nModes × 1 × nFreq]`.

7. **Package + save** as a standard timefreq file (model the save on
   `process_eigenmodes_spectrum`'s output + Brainstorm's `timefreq` template):
   `RowNames = {'Mode 1 (λ=…)', …}`, `Freqs`, `Method`, `nComponents = 1`,
   `SurfaceFile`, `DataType = 'matrix'` (rows treated as named modes, not sources),
   `DataFile = ''`; `db_add_data` into the input's study. History records method,
   #modes, mass type, eigenvalue range, source file, and time window.

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
- **Parity test (the key correctness check)** — the linear relation
  `mode = Φ'·M · source` holds on the **complex** spectrum (power is nonlinear), so
  run both with **FFT, measure = 'none'** (complex): standard source FFT
  (`process_fft`) → `X_src`, mode FFT → `X_mode`; assert
  `X_mode == (Φ'·M) · X_src` to numerical tolerance. Repeat on a raw file **with a
  marked bad segment** — equality there proves the bad-segment windows excluded are
  identical. Prints `ALL TESTS PASSED`.
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
- **Bad-segment exclusion (raw):** must match standard Brainstorm exactly —
  otherwise the mode-frequency spectrum would diverge from the source
  time-frequency spectrum on the same data. Guaranteed here by passing the same
  `BadSegments` into the same `bst_psd` call the source path uses (only the kernel
  differs). The plan must verify parity against a source PSD on a raw file with a
  marked bad segment.
- **Mode capping vs available modes:** request 300 but a surface may have fewer;
  clamp and report.
- **Mass-matrix provenance:** prefer the saved `Eigenmodes.MassMatrix`; only
  recompute when absent, and use the recorded `MassType` so the projection matches
  how the eigenmodes were built.
- **Output tagging:** since we build/save the timefreq file ourselves, set the
  fields (`DataType = 'matrix'`, `RowNames` = mode labels, `Method`, `Freqs`) so the
  file opens in the spectrum viewer with per-mode rows (model on a standard PSD
  timefreq struct + `process_eigenmodes_spectrum`'s save).
- **Frequency resolution / window length** come from the standard PSD panel; no new
  validation beyond what `process_psd` already does.
