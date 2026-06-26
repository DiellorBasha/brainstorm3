# PREVENT-AD MEG import + preprocessing + source mapping

**Date:** 2026-06-25
**Author:** Diellor Basha (with Claude)
**Status:** Approved — implementation

## Goal

Adapt `toolbox/script/tutorial_omega.m` into a custom script that imports,
preprocesses, and source-maps the PREVENT-AD BIDS MEG dataset into the existing
`preventad` Brainstorm protocol. Two non-negotiable customizations versus the
tutorial:

1. **ico5 anatomy downsampling** (20484 vertices total) is a hard requirement —
   the differential/spectral analysis engine depends on the FreeSurfer
   sphere-registered icosphere correspondence; reducepatch breaks it.
2. **Dirac dSPM inverse** in addition to the standard dSPM inverse, both built
   on a single unconstrained overlapping-spheres head model.

Test on a single subject (`sub-MTL0002`) before any batch run.

## Inputs / signature

```
preventad_import(BidsDir, SubjectNames)
```

- `BidsDir` (default `/Volumes/SpikeData-2/workspace/library/datasets/preventad/meg`)
- `SubjectNames` — comma-string or cellarray; empty = all subjects. The
  single-subject test passes `'sub-MTL0002'`; the batch run passes `''` (all) or
  an explicit list.

Location: `dev/preventad_import.m` (dev scratch script, not a shipped tutorial).

## Key facts established from the codebase

- `process_import_bids` exposes `downsamplemethod` (`'icosphere'`/`'reducepatch'`)
  and `icolevel` (`'ico3'..'ico6'`). `icolevel='ico5'` → 20484 verts; icosphere is
  FreeSurfer-only and auto-detects `derivatives/freesurfer/sub-*`. On the icosphere
  path `nvertices` is ignored.
- `process_import_bids` exposes `selectsubj` (comma-separated names, empty=all) →
  single-subject test and batch use the same script.
- `process_inverse_dirac` accepts a **standard unconstrained surface head model**
  (`Gain = [nChan x 3*nVert]`, i.e. overlapping spheres) and performs the
  Dirac-basis transform + reconstruction internally. No separate Dirac head-model
  node is required. Defaults: 400 modes/hemi, tau=0.5, SNR=3, noisereg=0.1,
  measure `dspm2018`.
- A single overlapping-spheres head model feeds both inverses (orientation
  constraint is applied at inverse time, not at head-model time).
- `sub-MTL0002` has MEG under `ses-02/meg` (task-rest run-01/02 + task-noise);
  FreeSurfer anatomy under `derivatives/freesurfer/sub-MTL0002`.

## Pipeline (deltas from tutorial_omega in **bold**)

1. **Import BIDS** — `process_import_bids` with **`downsamplemethod='icosphere'`,
   `icolevel='ico5'`**, `channelalign=0`, **`selectsubj=SubjectNames`**.
2. Remove head points → refine registration → convert to continuous CTF.
3. **Resample 2400 → 1200 Hz** (`process_resample`, freq=1200) done FIRST so all
   filtering runs on half the samples and the surviving continuous file is
   half-size. New Nyquist (600 Hz) keeps every notched harmonic (60-300 Hz) and
   the full analysis band (<=90 Hz) → lossless for this analysis.
4. Preprocessing — notch 60/120/180/240/300 Hz, high-pass 0.3 Hz, PSD-after +
   spectrum snapshot, delete intermediate (raw, resample, notch) folders.
5. Artifact cleaning — select `task-rest`; detect cardiac (ECG), blinks (VEOG),
   saccades (HEOG); remove cardiac events coinciding with blinks (±0.25s); SSP
   ECG (cardiac) + SSP EOG (blink) + SSP EOG (saccade), all MEG, usessp/select=1.
   Then AUTOMATIC bad-segment detection (`process_evt_detect_badsegment`,
   sensitivity 3, 1-7Hz movement + 40-240Hz muscle/SQUID) renamed to `bad_*` so
   the segments are excluded downstream. Registration + SSP snapshots.
   (PREVENT-AD rest runs carry ECG, VEOG, HEOG; the noise scan has none.)
   NOTE on bad CHANNELS: Brainstorm has no recommended fully-automatic MEG
   bad-channel detector (`process_detectbad` carries an explicit "not
   recommended / manual inspection only" warning; `process_detectbad_mad` is
   segment-based). Not included by default to avoid killing good channels.
6. Source estimation:
   - Noise covariance from `task-noise`.
   - ONE overlapping-spheres head model (unconstrained).
   - **NEW: Manifold** via `process_tess_manifold` (gauge=trivial) on each
     subject's cortex, BEFORE the inverses (find-or-create manifold_ node for
     later differential/spectral analysis).
   - Standard **dSPM** via `process_inverse_2018`, **UNCONSTRAINED** (`SourceOrient
     'free'`, nComponents=3) — matches the Dirac dSPM source space.
   - **NEW: Dirac dSPM** via `process_inverse_dirac` on the same head model
     (process defaults; also unconstrained 3-vector).
7. Power maps — **standard dSPM AND Dirac dSPM**: PSD (6 bands) → relative norm →
   project to default anatomy → spatial smooth 3 → average → contact-sheet snapshot.
   Standard runs on the per-file dSPM results directly; Dirac (a SHARED kernel,
   DataFile='') requires resolving its per-recording LINK nodes first (process_psd
   can't run on the bare shared kernel). Dirac PSD comment tagged `,Dirac` to keep
   the two distinguishable (esp. in Group_analysis). [2026-06-26: Dirac added; the
   first 66 imported subjects predate this and lack Dirac power maps -> backfill
   them via a GUI process job.]
8. Save report (display wrapped in try/catch: `bst_report('Open')` hits a removed
   JDK class `sun.misc.BASE64Encoder` on Apple silicon; the report .mat is still
   saved).

## New process built for this pipeline

`toolbox/process/functions/process_tess_manifold.m` — thin pipeline wrapper
around `tess_manifold`. Resolves the unique subjects in the input files, computes
(find-or-create) a `manifold_` node on each subject's cortex (gauge option,
ForceRecompute option), passes inputs through unchanged. Anatomy-level, idempotent.

## Storage (measured)

- Per subject @ 2400 Hz: ~3.05 GB (anat 0.34 + data 2.7).
- Per subject @ 1200 Hz + manifold: ~2.02 GB (continuous halves; anat/headmodel/
  kernels are sampling-independent so total drops ~34%, not 50%; manifold +19 MB).
- Full ico5-importable set (111 subjects): ~237 GB. SpikeData-2 free: 920 GB.

## Protocol handling

The script imports into the **existing** `preventad` protocol — it does NOT
delete/recreate the protocol (unlike `tutorial_omega`, which wipes `TutorialOmega`).
This makes per-subject re-runs and incremental batch additions safe.

## Validation (single-subject test)

Run on `sub-MTL0002`, then confirm:
1. Imported cortex surface has 20484 vertices (ico5).
2. Both inverse kernels are present in the DB (standard dSPM + Dirac dSPM).
3. Report contains no errors; spectrum/registration/SSP/power snapshots present.

Once green, the same script runs the batch by passing the full subject list (or
empty for all).

## Decisions (from brainstorming)

- Dirac basis: one head model, internal projection (no separate node).
- Power maps: standard dSPM only.
- Dirac params: process defaults.
- Script location: `dev/`.
- Power-map tail (project-to-default-anatomy + average across the subject's two
  rest runs): kept as-is for parity with the tutorial.
