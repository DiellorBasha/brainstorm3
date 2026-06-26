# PET metadata in the Brainstorm data model — design

**Date:** 2026-06-26
**Branch:** `experimental/pet`
**Status:** Design approved — pending implementation plan

## Goal

Make tracer / timing / scanner metadata a first-class part of the Brainstorm PET
data model. Unlike structural MRI, a PET volume is uninterpretable without its
tracer type, frame timing, injection timing, and scanner/recon context. That
metadata must be captured at import (from the BIDS JSON sidecar when available,
the NIfTI header as fallback), carried through processing, and inspectable from
the database tree. Every downstream PET step (temporal windowing, SUVR, PVC,
surface projection) will read from this single authoritative place.

This is the **first** piece of a larger effort to produce a robust static SUVR
map on the cortical surface. It deliberately does **not** implement SUVR,
temporal windowing, surface projection, or partial-volume correction.

## Context — what already exists

- `import_mri.m` already detects PET (sets `volType='PET'`, tags filename
  `_volpet`) purely from the `'PET'` substring in the import comment, and
  `tree_callbacks.m` types a node as `volpet` from that `_volpet` filename tag.
- `sMri.Header` already stores the **full decoded NIfTI header** (`hdr.dim`,
  `hdr.nifti`, ...). So header data is present but raw/unstructured.
- The BIDS JSON sidecar (tracer, injection, frame timing, scanner) is **not**
  captured anywhere — it is discarded at import.
- `preventad_pet_import.m` imports the raw 4D PET, motion-corrects (`mri_realign`),
  coregisters + reslices to T1 (`mri_coregister`), deletes intermediates, and
  saves a single `PET <tracer>` base node. The original BIDS path (and its JSON
  sidecar) is available at import time but not propagated.
- The PET tree context menu is built in `tree_callbacks.m` `fcnPetProcessing`
  (currently: Realign frames / Compute SUVR / Project volume to surface).

## Grounding facts (PREVENT-AD, verified from the data)

Two tracers per subject, HRRT scanner, `Units = nCi/ml`, `ImageDecayCorrected =
true`, `ReconFilterType = none` (near-intrinsic resolution, ~2.5–3 mm — raises
the importance of PVC in a later piece). Frame timing decoded as
`post-injection = FrameTimesStart − InjectionStart`:

| Tracer (`trc-`)      | Role    | InjectionStart | Frames   | Coverage post-injection | Standard window |
|----------------------|---------|----------------|----------|-------------------------|-----------------|
| `18FNAV4694`         | amyloid | −2400 s        | 6 × 5 min| **40–70 min**           | ~50–70 min      |
| `18Fflortaucipir`    | tau     | −4800 s        | 4 × 5 min| **80–100 min**          | 80–100 min      |

Both acquisitions are already late-window scans; the short frames exist mainly to
enable inter-frame motion correction. (This is why the later windowing step is
about validation + provenance, not window *selection* — out of scope here.)

## Data model — `sMri.PET`

A curated/derived summary plus the verbatim decoded JSON, attached to the PET
volume's MRI struct:

```
sMri.PET
  .Source        'bids-json' | 'nifti-header' | 'none'   provenance of the metadata
  .Tracer        .Name .Radionuclide .Units
  .Injection     .InjectedRadioactivity .InjectedRadioactivityUnits
                 .Mode .TimeZero .InjectionStart .ScanStart
  .Frames        .TimesStart [1xN] .Duration [1xN] .N
                 .MidTimes [1xN]            (derived: TimesStart + Duration/2)
                 .CoverageMinPI [t0 t1]     (derived, minutes post-injection:
                                             ([min TimesStart, max TimesStart+Duration]
                                              − InjectionStart)/60)
  .Decay         .ImageDecayCorrected .ImageDecayCorrectionTime
  .Scanner       .Manufacturer .Model .ReconMethod
                 .ReconFilterType .ReconFilterSize .AttenuationCorrection
  .Json          <raw decoded JSON struct, verbatim; [] if header-only>
```

- Fields not available from the source are set to `[]` / `'n/a'` and `Source`
  reflects what was found. Downstream code detects missing data via `Source` and
  empty fields.
- Derived fields are computed once at import.

## Components

### 1. `pet_read_metadata.m` (new — shared helper)

`PET = pet_read_metadata(PetFile, sMri)`

- **JSON (primary):** locate the BIDS sidecar next to `PetFile`
  (`*_pet.nii`/`*_pet.nii.gz` → `*_pet.json`, skipping `._*` AppleDouble files),
  decode it (`bst_jsondecode`), map fields into the curated schema, keep the raw
  decode in `PET.Json`, set `Source='bids-json'`.
- **NIfTI-header fallback:** if no JSON, fill what `sMri.Header` provides (frame
  count from `dim`, durations if encoded, voxel size, units), mark
  tracer/injection/scanner unavailable, set `Source='nifti-header'`.
- **None:** neither available → `Source='none'`, empty schema.
- Pure read/parse — no DB side effects. Directly unit-testable against the BIDS
  JSON files.

### 2. Core wiring — `import_mri.m`

Inside the existing PET branch (sets `volType='PET'` / `_volpet`), call
`pet_read_metadata(MriFile, sMri)` and assign `sMri.PET`. Confined to the
already-PET-specific path; every PET import inherits metadata capture (GUI BIDS
import, tutorial, and `preventad_pet_import`).

### 3. Carry-through (correctness requirement)

Realign / coregister / reslice preserve frame count, so `.PET` timing stays valid
through them. `preventad_pet_import.m` re-attaches `sBase.PET = sImp.PET` before
the final `out_mri_bst` save, so the surviving `PET <tracer>` base node carries
the metadata (the raw import it derives from is deleted). No change to
`mri_realign` / `mri_coregister`.

### 4. GUI inspectability — `tree_callbacks.m`

Add a **"PET information"** item to the `fcnPetProcessing` submenu. On click, load
the volume and show a **formatted read-only dialog** (`java_dialog('msgbox', …)`)
containing: tracer + units; injection summary; a **frames table** (start /
duration / mid-time / minutes-post-injection); the coverage window; decay-
correction status; scanner / recon. If `.PET` is absent (volumes imported before
this change) the dialog says so plainly. A simple text dialog (not a custom
panel) for v1 — lowest risk, immediately auditable.

## Validation (first 2 subjects: `sub-MTL0002`, `sub-MTL0005`)

1. **Unit test** `pet_read_metadata` directly against the 2 subjects' JSON
   sidecars (both tracers): assert tracer name, units (`nCi/ml`), frame counts
   (6 amyloid / 4 tau), and derived coverage (40–70 / 80–100 min).
2. **NIfTI fallback test:** call `pet_read_metadata` with the sidecar hidden/absent
   and assert `Source='nifti-header'` with frame count still populated and
   tracer/injection marked unavailable.
3. **End-to-end:** on the `experimental/pet` protocol, delete the existing PET
   bases for the 2 test subjects (imported before this change) and re-run
   `preventad_pet_import` for them; confirm each surviving base node carries
   `.PET` and that "PET information" renders correctly for both tracers.

## Scope guard (explicitly out of scope)

No SUVR, no temporal windowing/aggregation, no surface projection, no PVC. This
piece only **captures, carries, and displays** the metadata those later pieces
will depend on.

## Files touched

- `toolbox/anatomy/pet_read_metadata.m` — new helper.
- `toolbox/io/import_mri.m` — call helper in PET branch.
- `dev/preventad_pet_import.m` — re-attach `.PET` before save.
- `toolbox/tree/tree_callbacks.m` — "PET information" menu + display callback.
- `dev/test_pet_metadata.m` — new 2-subject validation script (unit + fallback +
  end-to-end checks).
