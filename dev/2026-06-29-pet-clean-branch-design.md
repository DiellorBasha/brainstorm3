# Clean upstream-compatible PET pipeline branch — Design

**Date:** 2026-06-29
**Author:** Diellor Basha (with Claude Code)
**Status:** Approved design — proceeding to implementation plan

## Goal

Extract the validated, standard PET processing pipeline from the heavily
interleaved `development` branch into a clean branch that depends only on
official upstream Brainstorm — sequestering the experimental
differential-geometry work that depends on the development infrastructure.

Two deliverables, both off a freshened `master`:

- **Branch A — `feature/pet-pvc`**: toolbox-only standard PET pipeline
  (metadata capture, PVC, SUVR refinement, surface projection). Self-contained
  against official upstream. Pushed to origin.
- **Branch B — `benchmark/pet-pvc`**: branches off A, adds the standard-pipeline
  validation benchmarks, with a **tag** (`pet-pvc-v1`) to reference in the
  eventual pull request to `brainstorm-tools/brainstorm3`. Pushed to origin
  (branch + tag).

## Context / current state

- `master` is `0` ahead / `53` behind `upstream/master` → bringing it current is
  a clean **fast-forward**.
- `development` is `1057` ahead / `197` behind `master`. The ~40 standard-PET
  commits are interleaved among ~1000 unrelated commits (eigen, dynamics, dirac,
  connectome, atom designer), which is why cherry-picking history was rejected in
  favor of a clean reconstruct.
- Merge-base for the 3-way reconstruction:
  `4983e2f0` (`merge-base development upstream/master`).
- **Upstream already contains the base PET scaffold**: `pet_process.m`,
  `panel_process_pet.m`, `panel_import_pet.m`, `mri_interp_vol2tess.m`,
  `tutorial_pet_processing.m`. The work to extract is the **enhancement layer**
  (PVC / SUVR / GTM / metadata) on top.
- The standard pipeline is **self-contained**: a grep of the standard files for
  experimental/dev-only symbols (`pet_epicenter`, `pet_spread*`, `connectome`,
  `bst_eigen`, `bst_dirac`, `nxr_*`, `manifold`, `tess_operators`, `bst_dynamics`,
  `hodge`, `helmholtz`, `pet_pvc_surface`) returned nothing.

## Commit hygiene & attribution (hard requirement)

The deliverable branches (`feature/pet-pvc`, `benchmark/pet-pvc`) must read as
the developer's own work:

- **No AI traces.** Commits carry **no** `Co-Authored-By: Claude` and **no**
  `Claude-Session` trailers. Author is the developer (`Diellor Basha`).
- **No plan/design `.md` files** committed to the clean branches. This design
  doc and any implementation plan live **only on `development`** (the dev
  branch), never on `feature/pet-pvc` or `benchmark/pet-pvc`.
- Commit messages are written in the project's normal style, as if authored by
  hand.

## Strategy: reconstruct clean (not cherry-pick)

Take the final working state of the standard PET files from `development` and
re-apply them onto the freshened `master`, grouped into logical commits.
Robust against the interleaving; guarantees upstream-compatibility.

## Step 0 — Freshen master

1. `git fetch upstream` (done).
2. Fast-forward local `master` to `upstream/master`.
3. Push `origin master`.

## Step 1 — Build Branch A (`feature/pet-pvc`, off updated master)

### 1a. Handle existing branch name collision
`feature/pet-pvc` already exists (local + origin, the old PVC work,
tip `69fc3b87`). **Archive it**: rename local `feature/pet-pvc` →
`archive/pet-pvc-old` (kept locally for safety). The new clean branch reuses
`feature/pet-pvc`; the eventual `git push --force-with-lease origin feature/pet-pvc`
replaces `origin/feature/pet-pvc`.

### 1b. File inventory

**New files** (copy verbatim from `development`; absent upstream):
- `toolbox/anatomy/pet_pvc.m` — PETPVE12 Müller-Gärtner volume PVC
- `toolbox/anatomy/pet_gtm.m` — Rousset GTM volume PVC
- `toolbox/anatomy/pet_suvr.m` — robust SUVR reference (erosion + trimmed mean)
- `toolbox/anatomy/pet_scanner_fwhm.m` — PSF FWHM from scanner metadata
  (helper called by both `pet_pvc` and `pet_gtm`)
- `toolbox/anatomy/pet_read_metadata.m` — BIDS-JSON / NIfTI metadata capture

**Modified PET files** (re-apply development changes onto upstream's current
version via 3-way against merge-base `4983e2f0`; all changes in these files are
PET, resolve toward upstream structure + the PET additions):
- `toolbox/anatomy/pet_process.m`
- `toolbox/gui/panel_process_pet.m`
- `toolbox/gui/panel_import_pet.m`
- `toolbox/anatomy/mri_interp_vol2tess.m`
- `toolbox/script/tutorial_pet_processing.m`

**Shared files** (port **only the PET hunks** onto upstream's version — these
carry ~1000 commits of unrelated dev changes, so no full-file take):
- `toolbox/io/import_mri.m` — PET import branch + `pet_read_metadata` call
- `toolbox/tree/tree_callbacks.m` — PET-processing menu + "PET information"

**Explicitly excluded / sequestered** (remain on `development`):
- `toolbox/anatomy/pet_pvc_surface.m` (surface-native PVC — experimental)
- `toolbox/anatomy/pet_epicenter.m` (Morse-Smale / persistence)
- `toolbox/anatomy/pet_spread_simulate.m`, `toolbox/anatomy/pet_spread_invert.m`
  (reaction-diffusion coupling)
- all `dev/` scripts and benchmarks

### 1c. Commit grouping (logical, each builds on the last)
1. **PET metadata capture** — `pet_read_metadata` + `import_mri` PET branch +
   `panel_import_pet` + tree "PET information"
2. **Robust SUVR reference** — `pet_suvr` + `pet_process` SUVR wiring
3. **Volume PVC (Müller-Gärtner)** — `pet_pvc` + `pet_scanner_fwhm` +
   pipeline/panel wiring
4. **Volume PVC (Rousset GTM)** — `pet_gtm` + pipeline/panel wiring
5. **Surface projection refinement** — `mri_interp_vol2tess` + `pet_process`
   mid-centered profile
6. **Tutorial up to date** — `tutorial_pet_processing` + tree PET-processing menu

`pet_process.m` is touched across several commits — acceptable. If hunks are not
cleanly separable, adjacent commits are merged rather than forcing artificial
splits.

### 1d. Validation
Launch Brainstorm from Branch A in MATLAB against the clean tree; run
`tutorial_pet_processing.m` (import → PVC → SUVR → surface projection). Confirm
no missing-function / dev-dependency errors and that none of the sequestered
functions are referenced.

## Step 2 — Build Branch B (`benchmark/pet-pvc`, off A)

1. Add standard-pipeline validation benchmarks:
   - `dev/benchmarks/compare_pvc_petsurfer.m` (PVC vs PETsurfer parity)
   - `dev/benchmarks/bench_pet_surface_recovery.m` (synthetic vol→surface recovery)
   - Exclude `validate_surface_pvc_vs_petsurfer.m` (surface PVC sequestered) and
     the spread/epicenter benchmarks.
2. Verify these benchmarks reference no sequestered/dev-only infra. If any do,
   note it (benchmarks may reference cohort data paths and may not run headless;
   included as reference artifacts, validated only for dependency-cleanliness).
3. Commit, **tag `pet-pvc-v1`**, push origin (branch + tag).

## Risks

- Hunk-extraction from `import_mri.m` / `tree_callbacks.m` is the delicate part
  (large files, mixed PET + non-PET changes) — port carefully, diff against
  upstream.
- Upstream's 53-commit evolution may have touched `pet_process.m` / panels →
  3-way conflicts resolved preserving upstream structure + adding features.
- Benchmarks may not run headless (cohort data paths / plugins); treated as
  reference artifacts.
