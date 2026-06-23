# Helmholtz Maps at Events — Design

**Date:** 2026-06-23
**Author:** Diellor Basha
**File:** `toolbox/process/functions/process_helmholtz_events.m`

## Goal

Batch (process_) version of the `view_helmholtz` GUI: at each timepoint of a
phase-marker event group (e.g. `alpha_peak` from `process_evt_refphase`),
reconstruct the unconstrained source vector field and run the Helmholtz-Hodge
decomposition, storing the static source maps as results ("sources") files. First
step of the static source-map workflow (extrema detection / thresholding come
later, with the differential + eigenwavelet tools).

## Input

A **kernel-link results file** (`InputTypes = {'results'}`), UNCONSTRAINED
(`nComponents == 3`) — the (sensor data + inverse kernel) pairing the user drags
from the tree, exactly as `view_helmholtz` consumes. From it: `ImagingKernel`
`[3nV×nGood]`, `GoodChannel`, `SurfaceFile`, `DataFile`.

## Options

| Option | Type | Default |
|---|---|---|
| `eventname` | text (a simple phase-marker group) | `alpha_peak` |
| `freqband` | radio {delta,theta,alpha,beta,gamma,custom} | `alpha` |
| `freqrange` | freqrange | `[8 13] Hz` |

The band must match the band the events were detected in (the bandpass applied to
the sensors before reconstruction).

## Pipeline (mirrors view_helmholtz)

1. Load the kernel link → `ImagingKernel`, `GoodChannel`, `SurfaceFile`, `DataFile`.
   Guard: `nComponents==3` and `ImagingKernel` non-empty.
2. Load the linked recording (`data` or `raw`) → sensors `F`, `TimeVector`, events.
   Find the `eventname` group; `iT = bst_closest(eventTimes, TimeVector)`.
3. Bandpass the good channels: `Fbp = process_bandpass(F(GoodChannel,:), band)`.
4. Resolve operators once and prepare:
   `Dirac/LBO = bst_get_operator_node(SurfaceFile,…)`, `Mani = tess_manifold`,
   `Surf = in_tess_bst`, `Op = bst_helmholtz('Prepare', {Dirac,LBO}, Mani, Surf, 'Domain','vertex')`.
5. Per event: `J = ImagingKernel * Fbp(:,iT(k))` `[3nV×1]`;
   `Ht = bst_helmholtz('Frame', Op, J, false)` (cores off — detection is a later step).
   Collect `Ht.Fmag`, `Ht.Phi`, `Ht.Psi`, and `reshape(Ht.Vtot',[],1)`.
6. Save 4 results files, `Time` = the event times (one frame per event):

| File | `ImageGridAmp` | nComp | Colormap |
|---|---|---|---|
| `<event> \| Source vector J` | `Vtot` (interleaved) | 3 | — |
| `<event> \| Total field \|J\|` | `Fmag` | 1 | source |
| `<event> \| Potential Φ` | `Phi` | 1 | stat2 (signed) |
| `<event> \| Stream Ψ` | `Psi` | 1 | stat2 (signed) |

## Gotchas handled

- **Output folder**: derived from `R.DataFile`, NOT `bst_fileparts(sInput.FileName)`
  — the input is a `link|kernel|data` filename whose `bst_fileparts` is a bogus
  `link/` folder (caused a spurious "disk full" write error).
- **Single event**: a one-frame map is duplicated so the viewer has a time axis.
- **Vector layout**: `Ht.Vtot` is `[nV×3]`; stored as `[3nV×1]` interleaved
  (`reshape(Vtot',[],1)`), matching Brainstorm's unconstrained convention.
- **Raw vs data**: `i_load_recording` handles both (`in_fread` for raw links).

## Verification (MCP, TutorialAuditory)

Input: `MN: MEG(Unconstr) 2018` kernel link `[61452×272]` on
`data_block001_02` (alpha_peak = 64 events). Output: 4 maps, 64 frames, **1.6 s**
(cached Cholesky in `Prepare`). Sanity: `|J| ≥ 0`; Φ, Ψ signed; vector field
`3×20484` rows. ✓

## Scope / next

- No eigen, no frequency stacking, no extrema detection yet (deliberately simple).
- Next: isolate source/sink (Φ extrema) and vortex/antivortex (Ψ extrema) with
  source-domain amplitude thresholding (a spatial analogue of the GFP gate), then
  fold time × source × frequency × eigenmode into one descriptor table.
- Not committed yet (untracked on `development`).
