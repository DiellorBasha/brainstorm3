# Clean upstream-compatible PET pipeline branch — Implementation Plan

> **For agentic workers:** This is git + MATLAB reconstruction work, not codebase
> TDD. Steps use checkbox (`- [ ]`) syntax for tracking. Each task ends with a
> verification gate (git/MATLAB command + expected output) before committing.
> Recommended execution: **inline in the current session** — several steps need
> human/Claude judgment (isolating PET hunks, 3-way conflict resolution) and the
> MATLAB MCP for validation.

**Goal:** Reconstruct the validated standard PET pipeline (metadata, PVC, SUVR,
surface projection) from `development` onto a freshened `master` as a clean
`feature/pet-pvc` branch that depends only on official upstream Brainstorm, plus
a `benchmark/pet-pvc` branch (tagged) carrying the validation benchmarks.

**Architecture:** Upstream already has the PET scaffold (`pet_process`,
`panel_process_pet`, `panel_import_pet`, `mri_interp_vol2tess`, tutorial). We add
the enhancement layer (new files) and port the PET-only deltas onto upstream's
current versions of the existing files, then carve into logical, working commits.

**Tech Stack:** git (origin = DiellorBasha/brainstorm3, upstream =
brainstorm-tools/brainstorm3), MATLAB R2023b + Brainstorm, MATLAB MCP for
validation.

## Global Constraints

- **No AI traces on deliverable branches.** Commits on `feature/pet-pvc` and
  `benchmark/pet-pvc` carry NO `Co-Authored-By: Claude` and NO `Claude-Session`
  trailers. Author = `Diellor Basha <diellorbasha@gmail.com>` (current git
  config — verify, don't override).
- **No plan/design `.md` on deliverable branches.** Never `git add` any
  `dev/*.md` or design/plan doc onto `feature/pet-pvc` / `benchmark/pet-pvc`.
  These docs live on `development` only.
- **Upstream-only dependencies.** No standard PET file may reference
  `pet_epicenter`, `pet_spread_*`, `pet_pvc_surface`, `connectome`, `bst_eigen`,
  `bst_dirac`, `nxr_*`, `manifold`, `tess_operators`, `bst_dynamics`, `hodge`,
  `helmholtz`.
- **Sequestered (never copy to clean branches):** `pet_pvc_surface.m`,
  `pet_epicenter.m`, `pet_spread_simulate.m`, `pet_spread_invert.m`, all `dev/`
  scripts (benchmarks handled separately on `benchmark/pet-pvc`).
- **Key refs:** merge-base = `4983e2f0`; old PVC branch tip = `69fc3b87`;
  `master` is 0 ahead / 53 behind `upstream/master` (clean fast-forward).

---

## Phase 0 — Pre-flight

### Task 0.1: Stash development working-tree changes

The working tree has uncommitted changes (`dev/preventad_import.m`, benchmark
PNGs) that belong on `development`. Stash so branch operations start clean.

- [ ] **Step 1: Confirm current branch and dirty state**

Run:
```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git rev-parse --abbrev-ref HEAD   # expect: development
git status --short | head
```
Expected: branch `development`; a list of `M`/`??` entries.

- [ ] **Step 2: Stash including untracked, with a clear label**

Run:
```bash
git stash push -u -m "WIP dev working tree before pet-pvc clean branch (2026-06-29)"
git stash list | head -1
```
Expected: one stash `stash@{0}: ... pet-pvc clean branch (2026-06-29)`.
Record this stash ref — it is restored in Task 7.1.

- [ ] **Step 3: Verify clean tree**

Run: `git status --short`
Expected: empty output.

### Task 0.2: Freshen `master` to `upstream/master` and push

- [ ] **Step 1: Confirm upstream is fetched and the delta is a fast-forward**

Run:
```bash
git fetch upstream
git rev-list --left-right --count master...upstream/master
```
Expected: `0	53` (master 0 ahead, 53 behind → fast-forward safe).

- [ ] **Step 2: Fast-forward local master**

Run:
```bash
git checkout master
git merge --ff-only upstream/master
git rev-list --left-right --count master...upstream/master
```
Expected: `0	0` after merge.

- [ ] **Step 3: Push freshened master to origin**

Run:
```bash
git push origin master
```
Expected: origin master updated (fast-forward, no force).

---

## Phase 1 — Branch setup

### Task 1.1: Archive the existing `feature/pet-pvc`

- [ ] **Step 1: Confirm the old branch tip before renaming**

Run: `git log --oneline -1 feature/pet-pvc`
Expected: `69fc3b87 Fix pet_pvc.m: use bst_plugin Install+Load ...`

- [ ] **Step 2: Rename local branch to archive**

Run:
```bash
git branch -m feature/pet-pvc archive/pet-pvc-old
git rev-parse --abbrev-ref archive/pet-pvc-old >/dev/null && echo "archived"
```
Expected: `archived`. (origin/feature/pet-pvc still points at `69fc3b87`; it is
replaced in Task 5.2 via `--force-with-lease`.)

### Task 1.2: Create the clean `feature/pet-pvc` off updated master

- [ ] **Step 1: Branch from master**

Run:
```bash
git checkout master
git checkout -b feature/pet-pvc
git log --oneline -1
```
Expected: HEAD == current `master` tip (the freshened upstream tip).

- [ ] **Step 2: Confirm author identity (no override)**

Run: `git config user.name && git config user.email`
Expected: `Diellor Basha` / `diellorbasha@gmail.com`.

---

## Phase 2 — Diff audit (drives all porting)

### Task 2.1: Audit the per-file PET delta vs upstream

This determines exactly what to port. The "modified" PET files may already be
near-identical to upstream (upstream advanced 53 commits); the real deltas may be
small.

- [ ] **Step 1: Diffstat of every standard PET file, upstream → development**

Run:
```bash
for f in toolbox/anatomy/pet_process.m toolbox/gui/panel_process_pet.m \
         toolbox/gui/panel_import_pet.m toolbox/anatomy/mri_interp_vol2tess.m \
         toolbox/script/tutorial_pet_processing.m \
         toolbox/io/import_mri.m toolbox/tree/tree_callbacks.m; do
  echo "=== $f ==="
  git diff --stat upstream/master..development -- "$f"
done
```
Expected: a diffstat per file. Record which files have large vs small deltas.

- [ ] **Step 2: Full PET-delta for the two SHARED files (manual hunk selection)**

Run:
```bash
git diff upstream/master..development -- toolbox/io/import_mri.m   > /tmp/pet_import_mri.diff
git diff upstream/master..development -- toolbox/tree/tree_callbacks.m > /tmp/pet_tree_callbacks.diff
wc -l /tmp/pet_import_mri.diff /tmp/pet_tree_callbacks.diff
```
Then read both diffs and classify each hunk as **PET** (keep) or **non-PET dev**
(drop). PET hunks to keep:
- `import_mri.m`: the `pet_read_metadata` call + storing PET metadata on `sMri`
  (and any PET import branch additions not already upstream).
- `tree_callbacks.m`: the "PET information" context-menu entry and any
  PVC/SUVR-related PET-processing menu additions not already upstream.
Expected output of this step: a written list of the exact PET hunks (line
ranges) to port. This list is the input to Task 3.3.

- [ ] **Step 3: Confirm the NEW files are absent upstream**

Run:
```bash
for f in pet_pvc pet_gtm pet_suvr pet_scanner_fwhm pet_read_metadata; do
  git cat-file -e upstream/master:toolbox/anatomy/$f.m 2>/dev/null \
    && echo "PRESENT $f" || echo "absent  $f"
done
```
Expected: all five `absent` (they are new files to add verbatim).

---

## Phase 3 — Reconstruct the working tree (materialize + validate before commit)

Build the full final state first, validate the whole pipeline, then carve commits
in Phase 4.

### Task 3.1: Add the new standard PET files verbatim

**Files:** Create on branch (from `development`):
`toolbox/anatomy/pet_read_metadata.m`, `pet_suvr.m`, `pet_pvc.m`,
`pet_scanner_fwhm.m`, `pet_gtm.m`.

- [ ] **Step 1: Check the files out of development into the working tree**

Run:
```bash
git checkout development -- \
  toolbox/anatomy/pet_read_metadata.m \
  toolbox/anatomy/pet_suvr.m \
  toolbox/anatomy/pet_pvc.m \
  toolbox/anatomy/pet_scanner_fwhm.m \
  toolbox/anatomy/pet_gtm.m
git status --short
```
Expected: five `A`/`M` entries staged (checkout stages them).

- [ ] **Step 2: Verify they carry no sequestered/dev deps**

Run:
```bash
grep -nE 'pet_epicenter|pet_spread|pet_pvc_surface|connectome|bst_eigen|bst_dirac|nxr_|manifold|tess_operators|bst_dynamics|hodge|helmholtz' \
  toolbox/anatomy/pet_read_metadata.m toolbox/anatomy/pet_suvr.m \
  toolbox/anatomy/pet_pvc.m toolbox/anatomy/pet_scanner_fwhm.m \
  toolbox/anatomy/pet_gtm.m || echo "CLEAN"
```
Expected: `CLEAN`.

### Task 3.2: Materialize the modified PET-only files (3-way onto upstream)

**Files:** Modify `toolbox/anatomy/pet_process.m`,
`toolbox/gui/panel_process_pet.m`, `toolbox/anatomy/mri_interp_vol2tess.m`,
`toolbox/gui/panel_import_pet.m`, `toolbox/script/tutorial_pet_processing.m`.

All changes in these files are PET, so apply the development delta over the
merge-base via 3-way; conflicts mean upstream also touched the region — resolve
keeping upstream structure + the PET additions.

- [ ] **Step 1: For each file, apply the dev delta with 3-way**

Run (per file; example shown for one, repeat for all five):
```bash
for f in toolbox/anatomy/pet_process.m toolbox/gui/panel_process_pet.m \
         toolbox/anatomy/mri_interp_vol2tess.m toolbox/gui/panel_import_pet.m \
         toolbox/script/tutorial_pet_processing.m; do
  echo "=== $f ==="
  git diff 4983e2f0..development -- "$f" | git apply --3way --whitespace=nowarn - \
    && echo "OK $f" || echo "CONFLICT $f"
done
git status --short
```
Expected: `OK` for each, or `CONFLICT` flagged.

- [ ] **Step 2: Resolve any conflicts**

For any file marked `CONFLICT`, open it, find `<<<<<<<`/`=======`/`>>>>>>>`
markers, and resolve by keeping upstream's structural changes while preserving
the PET feature additions. Then:
```bash
grep -rl '^<<<<<<<' toolbox/ || echo "NO MARKERS LEFT"
```
Expected: `NO MARKERS LEFT`.

- [ ] **Step 3: Confirm these files reference only upstream + new PET fns**

Run:
```bash
grep -nE 'pet_epicenter|pet_spread|pet_pvc_surface|connectome|bst_eigen|bst_dirac|nxr_|manifold|tess_operators|bst_dynamics|hodge|helmholtz' \
  toolbox/anatomy/pet_process.m toolbox/gui/panel_process_pet.m \
  toolbox/anatomy/mri_interp_vol2tess.m toolbox/gui/panel_import_pet.m \
  toolbox/script/tutorial_pet_processing.m || echo "CLEAN"
```
Expected: `CLEAN`.

### Task 3.3: Port PET hunks into the shared files

**Files:** Modify `toolbox/io/import_mri.m`, `toolbox/tree/tree_callbacks.m`.
Start from upstream's version (already on the branch); add ONLY the PET hunks
identified in Task 2.1 Step 2.

- [ ] **Step 1: Apply the selected PET hunks**

Using the hunk list from Task 2.1, apply only those hunks. Preferred mechanism —
build a reduced patch containing only the PET hunks and apply with 3-way:
```bash
# After hand-editing the two .diff files down to PET-only hunks:
git apply --3way --whitespace=nowarn /tmp/pet_import_mri.diff
git apply --3way --whitespace=nowarn /tmp/pet_tree_callbacks.diff
grep -rl '^<<<<<<<' toolbox/io/import_mri.m toolbox/tree/tree_callbacks.m || echo "NO MARKERS"
```
Expected: `NO MARKERS`. (If `git apply` rejects a reduced patch due to context,
fall back to hand-editing the two files to insert the PET blocks at the matching
locations.)

- [ ] **Step 2: Verify no non-PET dev changes leaked in**

Run:
```bash
git diff --stat upstream/master -- toolbox/io/import_mri.m toolbox/tree/tree_callbacks.m
git diff upstream/master -- toolbox/io/import_mri.m toolbox/tree/tree_callbacks.m \
  | grep -E '^\+' | grep -iE 'pet|metadata' | head
```
Expected: diffstat shows only modest additions; the added `+` lines are
PET/metadata-related only. Eyeball the full diff to confirm no unrelated dev code.

### Task 3.4: Static validation of the whole working tree

- [ ] **Step 1: Whole-tree grep for sequestered symbols (all changed files)**

Run:
```bash
git diff --name-only upstream/master | while read f; do
  grep -nE 'pet_epicenter|pet_spread|pet_pvc_surface|connectome|bst_eigen|bst_dirac|nxr_|manifold|tess_operators|bst_dynamics|hodge|helmholtz' "$f" && echo "  ^ in $f"
done || true
echo "scan done"
```
Expected: no matches before `scan done`.

- [ ] **Step 2: MATLAB syntax-check every changed .m file**

Via MATLAB MCP (or `matlab -batch`), run `checkcode` on each changed file and
confirm no parse errors (Brainstorm-idiom warnings are fine):
```matlab
files = strsplit(strtrim(evalc("!git diff --name-only upstream/master -- '*.m'")));
for i=1:numel(files); if isempty(files{i}); continue; end; checkcode(files{i}); end
```
Expected: no `Parse error` / `Syntax error` lines.

### Task 3.5: End-to-end MATLAB validation of the pipeline

- [ ] **Step 1: Launch Brainstorm from this branch**

Use the `brainstorm-dev:brainstorm-start` skill (adds this checkout to path and
starts Brainstorm). Confirm it boots without errors about missing PET functions.

- [ ] **Step 2: Run the PET tutorial / pipeline**

Run `tutorial_pet_processing.m` (or drive `pet_process` on a known PET subject):
import → static mean → PVC (MG and GTM) → robust SUVR → surface projection.
Expected: completes with no `Undefined function` errors and produces a
surface SUVR map. Record the run output.

- [ ] **Step 3: Gate**

If anything fails, fix in the working tree and re-run 3.4–3.5 before proceeding
to commits. Do not commit a non-working tree.

---

## Phase 4 — Carve feature-grouped commits

Working tree is fully materialized and validated. Now commit subsets in logical
order. New function files are standalone (safe at every intermediate commit);
`pet_process.m` wiring lands only after all referenced functions are committed.
**Every commit: plain message, no AI trailers.** First unstage everything:

- [ ] **Step 0: Reset the index (keep working tree)**

Run: `git reset -q` then `git status --short`
Expected: all reconstructed changes shown as unstaged/untracked, working tree
intact.

### Task 4.1: Commit 1 — PET metadata capture

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/anatomy/pet_read_metadata.m toolbox/io/import_mri.m
git commit -m "Capture BIDS/NIfTI PET metadata into the data model

Add pet_read_metadata (BIDS-JSON primary, NIfTI fallback) and wire it into
the MRI/PET import path so scanner and tracer metadata are stored on sMri.PET."
git log --oneline -1
```
Expected: commit created; `git show --stat HEAD` lists exactly those two files.

- [ ] **Step 2: Verify no trailers**

Run: `git log -1 --format='%B' | grep -iE 'claude|co-authored' || echo "CLEAN"`
Expected: `CLEAN`.

### Task 4.2: Commit 2 — Robust SUVR reference

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/anatomy/pet_suvr.m
git commit -m "Add robust SUVR reference (pet_suvr)

Reference-region SUVR with erosion + trimmed-mean for a stable reference
uptake estimate."
git show --stat HEAD | tail -3
```
Expected: commit lists only `pet_suvr.m`. Verify no trailers (as 4.1 Step 2).

### Task 4.3: Commit 3 — Volume PVC (Müller-Gärtner) + scanner FWHM

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/anatomy/pet_pvc.m toolbox/anatomy/pet_scanner_fwhm.m
git commit -m "Add Muller-Gartner volume PVC (pet_pvc) + scanner FWHM

PETPVE12-based Muller-Gartner partial-volume correction; pet_scanner_fwhm
derives the PSF FWHM from scanner metadata when not supplied."
git show --stat HEAD | tail -3
```
Expected: lists only those two files. Verify no trailers.

### Task 4.4: Commit 4 — Volume PVC (Rousset GTM)

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/anatomy/pet_gtm.m
git commit -m "Add Rousset GTM regional volume PVC (pet_gtm)

Geometric transfer matrix regional PVC; extracerebral signal modelled as a
nuisance region."
git show --stat HEAD | tail -3
```
Expected: lists only `pet_gtm.m`. Verify no trailers.

### Task 4.5: Commit 5 — Wire PVC/SUVR/GTM + surface projection into the pipeline

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/anatomy/pet_process.m toolbox/gui/panel_process_pet.m \
        toolbox/gui/panel_import_pet.m toolbox/anatomy/mri_interp_vol2tess.m
git commit -m "Wire PVC/SUVR into pet_process pipeline + GUI

pet_process performs average -> PVC (MG or GTM, in volume) -> robust SUVR
-> surface projection (mid-centered profile via mri_interp_vol2tess);
panel_process_pet exposes the PVC/SUVR options."
git show --stat HEAD | tail -6
```
Expected: lists those four files. Verify no trailers.

### Task 4.6: Commit 6 — Tutorial + tree menus

- [ ] **Step 1: Stage + commit**

Run:
```bash
git add toolbox/script/tutorial_pet_processing.m toolbox/tree/tree_callbacks.m
git commit -m "Update PET tutorial + tree menus for PVC/SUVR pipeline

Bring tutorial_pet_processing up to date with the PVC/SUVR steps; add the
PET information / PET processing context-menu entries."
git show --stat HEAD | tail -3
```
Expected: lists those two files. Verify no trailers.

### Task 4.7: Branch-level verification

- [ ] **Step 1: Working tree must be clean (everything committed)**

Run: `git status --short`
Expected: empty (all reconstructed content is now committed; nothing stray).

- [ ] **Step 2: Branch contains only PET toolbox files, no docs, no trailers**

Run:
```bash
echo "--- files changed vs master ---"
git diff --name-only master..feature/pet-pvc
echo "--- any dev/docs/plan files? (should be none) ---"
git diff --name-only master..feature/pet-pvc | grep -E '^(dev/|docs/)' || echo "NONE"
echo "--- any AI trailers across the branch? (should be none) ---"
git log master..feature/pet-pvc --format='%B' | grep -iE 'claude|co-authored' || echo "NONE"
echo "--- author check ---"
git log master..feature/pet-pvc --format='%an <%ae>' | sort -u
```
Expected: file list is the standard PET toolbox set only; `NONE` for docs;
`NONE` for trailers; author is `Diellor Basha <diellorbasha@gmail.com>` only.

- [ ] **Step 3: Re-run MATLAB validation at branch HEAD**

Repeat Task 3.5 (Brainstorm boot + tutorial run) at the final commit to confirm
the carved history produces the same working pipeline.

---

## Phase 5 — Push `feature/pet-pvc`

### Task 5.1: Push, replacing the old origin branch

- [ ] **Step 1: Push with lease (origin still holds the old PVC tip 69fc3b87)**

Run:
```bash
git push --force-with-lease=feature/pet-pvc:69fc3b87 origin feature/pet-pvc
git rev-parse --abbrev-ref feature/pet-pvc@{upstream}
```
Expected: push succeeds; upstream tracking set to `origin/feature/pet-pvc`.
(If lease fails, the remote moved since fetch — stop and re-inspect, do not blind
force.)

---

## Phase 6 — `benchmark/pet-pvc` (validation + tag)

### Task 6.1: Branch off the clean pipeline

- [ ] **Step 1: Create the branch**

Run:
```bash
git checkout -b benchmark/pet-pvc
git log --oneline -1
```
Expected: HEAD == `feature/pet-pvc` tip.

### Task 6.2: Add the standard-pipeline validation benchmarks

**Files:** Create `dev/benchmarks/compare_pvc_petsurfer.m`,
`dev/benchmarks/bench_pet_surface_recovery.m` (volume-pipeline validation only;
exclude `validate_surface_pvc_vs_petsurfer.m`, spread/epicenter benchmarks).

- [ ] **Step 1: Check the benchmarks out of development**

Run:
```bash
git checkout development -- \
  dev/benchmarks/compare_pvc_petsurfer.m \
  dev/benchmarks/bench_pet_surface_recovery.m
```

- [ ] **Step 2: Verify they reference no sequestered/dev-only infra**

Run:
```bash
grep -nE 'pet_epicenter|pet_spread|pet_pvc_surface|connectome|bst_eigen|bst_dirac|nxr_|manifold|tess_operators|bst_dynamics|hodge|helmholtz' \
  dev/benchmarks/compare_pvc_petsurfer.m dev/benchmarks/bench_pet_surface_recovery.m \
  || echo "CLEAN"
```
Expected: `CLEAN`. If a benchmark references sequestered infra, drop it and note
the exclusion in the commit message.

- [ ] **Step 3: Commit (no trailers)**

Run:
```bash
git add dev/benchmarks/compare_pvc_petsurfer.m dev/benchmarks/bench_pet_surface_recovery.m
git commit -m "Add PET PVC/SUVR validation benchmarks

PVC-vs-PETsurfer parity comparison and synthetic volume-to-surface recovery
benchmark for the standard PET pipeline."
git log -1 --format='%B' | grep -iE 'claude|co-authored' || echo "CLEAN"
```
Expected: `CLEAN`.

### Task 6.3: Tag and push

- [ ] **Step 1: Create the tag**

Run:
```bash
git tag -a pet-pvc-v1 -m "Standard PET pipeline + validation benchmarks (PR reference)"
git tag -n1 pet-pvc-v1
```
Expected: tag `pet-pvc-v1` listed.

- [ ] **Step 2: Push branch + tag**

Run:
```bash
git push origin benchmark/pet-pvc
git push origin pet-pvc-v1
```
Expected: both pushed (new branch + new tag, no force).

---

## Phase 7 — Cleanup

### Task 7.1: Return to development and restore stashed WIP

- [ ] **Step 1: Switch back and restore the stash**

Run:
```bash
git checkout development
git stash list | head -1
git stash pop
git status --short
```
Expected: back on `development`; the Task 0.1 stash restored (preventad_import.m,
PNGs) with no conflicts.

- [ ] **Step 2: Final sanity of the deliverables**

Run:
```bash
git branch --list 'feature/pet-pvc' 'benchmark/pet-pvc' 'archive/pet-pvc-old'
git ls-remote --heads origin feature/pet-pvc benchmark/pet-pvc
git ls-remote --tags origin pet-pvc-v1
```
Expected: all three local branches present; both remote branches and the tag
present on origin.

---

## Self-Review notes

- **Spec coverage:** Step 0 freshen master ✓; archive old branch ✓ (1.1);
  new files ✓ (3.1); modified PET files 3-way ✓ (3.2); shared-file PET hunks ✓
  (3.3); sequester enforced ✓ (constraints + 3.2/3.4 greps); 6 logical commits ✓
  (4.1–4.6); validation ✓ (3.5/4.7); benchmarks + tag branch ✓ (Phase 6);
  no-AI-trace + no-docs constraints ✓ (4.1–4.7 checks).
- **Judgment steps (need inline execution):** 2.1 Step 2 hunk classification,
  3.2 Step 2 conflict resolution, 3.3 hunk porting, all MATLAB validation.
- **Risk:** if upstream's 53 commits substantially rewrote `pet_process.m`, the
  3-way in 3.2 may need careful manual merge; the validation gate (3.5) catches
  regressions before any commit.
