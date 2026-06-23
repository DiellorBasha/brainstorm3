# Band-Power Period Detector — Design

**Date:** 2026-06-23
**Author:** Diellor Basha
**File:** `toolbox/process/functions/process_evt_refphase.m`

## Goal

Detect time windows where band-limited oscillatory **power** (amplitude) is high
on a continuous recording, and store them as an extended-event group
(onset/offset). First temporal marker for the joint spatiotemporal marker table.

## Motivation — why a new detector

The existing burst detectors (`process_evt_detect_burst` CFAR,
`process_evt_detect_burst_wavelet`) are **spectral-contrast** methods: they ask
"is there a narrowband peak sticking out of the broadband spectrum?" and are
deliberately **amplitude-invariant**. On an alpha-dominated segment alpha is a
sustained spectral peak, so they fire continuously — on the 38 s test block
`Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat` the CFAR default
produced **1 window covering 100%** of the file. That is the wrong tool for
"when is the band *strong*".

This detector measures **amplitude** with a **relative threshold on the signal's
own distribution**, ported from the artifact-detector architecture
(`process_evt_detect`: bandpass → rectify → `k·std` → crossings). It is
self-calibrating, so by construction it never marks everything. On the same
block, percentile thresholding gives selective windows (p85 → 8 windows, ~13%).

## Method — `Compute(F, TimeVector, OPTIONS, validMask) -> [evt(2×nP), markers, stats]`

**Amplitude / period detection:**
1. (optional) **Normalize** each channel by its robust scale over valid samples
   (mixed sensor types — equalize units before the GFP).
2. **Bandpass** to target band — `[Fbp, FiltSpec] = process_bandpass(... 'bst-hfilter-2019', 0)`;
   exclude `FiltSpec.transient` edge samples from the valid mask.
3. **Global Field Power** across sensors — `gfp = std(Fbp, 1, 1)` (the exact
   measure drawn by *Extra → Show GFP*, `figure_timeseries.m:3698`; Lehmann &
   Skrandies 1980). Single-channel fallback: `gfp = abs(Fbp)`.
4. **Smooth/lowpass** the GFP — `movmean` over the smoothing window. For a
   narrowband signal ~82% of GFP energy IS the slow envelope, so smoothing
   yields a clean amplitude trace directly (no Hilbert needed; verified on the
   test block: GFP peaks at 20.8 Hz, 0% energy at 8–13 Hz, 82% < 4 Hz).
5. **Relative threshold** by mode → enter `Thi`, exit `Tlo`, **estimated over
   valid samples only** (bad segments excluded so artifacts don't skew it):
   - `percentile`: `prctile(sig, enter/exit)` (0–100; defaults 85/75)
   - `mad`: `median + k·(1.4826·MAD)` (k; e.g. 2/1.5)
   - `std`: `median + k·std`; `Tlo` clamped ≤ `Thi`.
6. **Hysteresis scan** — enter when `sig ≥ Thi`, stay until `sig < Tlo`; an
   invalid (bad/edge) sample closes any open window. **Merge** gaps < `minGap`.

**Per-cycle phase markers (from GFP at 2f):** a band-limited GFP is sign-blind,
so it oscillates at **2f** with no fundamental. Bandpass the GFP at
`[2·fLow, 2·fHigh]` (clamped < Nyquist) for a near-sinusoid; its **local maxima**
(φ=0) are field-magnitude maxima = alpha extrema (consecutive peaks alternate
alpha+/−; the source field's own sign recovers polarity downstream), and its
**local minima** (φ=π) are alpha zero-crossings.
7. Within each candidate period, require ≥ 2 GFP peaks (genuine oscillation),
   **snap onset/offset to the first/last peak** (sign-consistent, integer cycles),
   then enforce `minDuration` on the snapped bounds.
8. Emit: extended `evt` (periods) + `markers.peak` / `markers.trough` (point
   events, only those inside kept periods). `markerMode ∈ {both,peak,trough,none}`;
   `none` ⇒ amplitude bounds, no snapping/markers.

`Run` stores three groups: `<name>` (extended), `<name>_peak`, `<name>_trough`
(simple). RMS-of-per-channel-envelope is near-equivalent to smoothed GFP
(r = 0.97) — GFP chosen for display fidelity and the free 2f phase reference.

## Process options (GUI) — deliberately minimal

| Option | Type | Default |
|---|---|---|
| `eventname` | text | `power_alpha` |
| `sensortypes` | text | `MEG` |
| `normalize` | checkbox (mixed sensor types) | off |
| `timewindow` | timewindow | [] (whole file) |
| `freqband` | radio_linelabel {delta,theta,alpha,beta,gamma,custom} | `alpha` |
| `freqrange` | freqrange | `[8 13] Hz` |
| `thresholdmode` | radio_linelabel {percentile, mad, std} | `percentile` |
| `enterthresh` | value | `85` |
| `exitthresh` | value | `75` |

**Auto-derived from the band center `fc = mean(freqRange)` (no GUI knob):**
- `smoothing = 0.5/fc` — a boxcar of ½ a period nulls the GFP's 2f ripple
  (empirically −23 dB suppression, 99% envelope kept at r=0.5; the non-monotonic
  dip at r=0.75 confirms the boxcar-null mechanism).
- `minDuration = 3/fc` (≥ 3 cycles), `minGap = 2/fc` (merge if closer).
- Phase markers are **always** emitted (peak + trough); the per-cycle marker is
  unambiguous, so no mode selector. (`Compute` still accepts explicit
  smoothing/minDuration/minGap overrides — `[]` means auto — for testing.)
- Menu label: **"Detect bursts (phase polarity)"** (canonical; CFAR/Wavelet
  retired to `dev/experimental/`).

`Run` reads raw broadband sensor data (continuous-only; empty-sensor and
empty-detection guards; mixed-sensor-type warning), masks bad segments, and
stores three event groups: `<name>` (extended), `<name>_peak`, `<name>_trough`
(simple). `bst_report` summary line counts periods + markers.

## Testing — `dev/test_evt_refphase.m` (7/7)

- **Synthetic topographic burst (MAD):** one localized period + peak/trough
  markers, all markers inside the period.
- **Real (`data_block001_02`):** 3–15 periods, markers within periods, monotonic.
- **2f marker frequency:** peak ISI ≈ 1/(2·f) (≈0.048 s for alpha) — confirms the
  GFP frequency-doubling.
- **Flat signal:** no events, no markers.
- **Bad-segment masking:** injected artifact excluded, real alpha still found.
- **markerMode=none:** amplitude bounds, no marker groups.
- **Normalization:** mixed-scale channels handled.

## Scope / notes

- Bad-segment masking, filter-transient trim, and optional per-channel
  normalization are implemented (ported from `process_evt_detect`).
- Input is **raw broadband** sensor data — the process bandpasses internally
  (target band for GFP/amplitude; `[2·fLow,2·fHigh]` on the GFP for phase markers).
  Do not pre-filter.
- Auto-discovered by Brainstorm (no registration needed).
- Not committed yet (untracked on `development`).
- Future: split `_peak` into polarity-consistent odd/even sub-trains (the source
  field sign recovers polarity downstream, so deferred).
