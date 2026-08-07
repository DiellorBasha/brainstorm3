# Gate 0 Report — Icosphere Clean-Branch Extraction

Date: 2026-08-07
Clean worktree: `~/workspace/research/code/brainstorm3-clean` (branch `feature/cortical-flow-core`)
Harness (this checkout): `~/workspace/research/code/brainstorm3/dev/verify/phase0/` (branch `feature/import-fibers-pipeline`)
Dataset: `/Volumes/SpikeData-2/workspace/library/datasets/omega-tutorial` (OMEGA tutorial BIDS dataset, sub-0002)

## 1. Clean branch commits

Base (merge-base with `master`): `172a5ac3` — "Bugfix: Process export 1 file, do not add Condition|Subject to filename"

Four commits on `feature/cortical-flow-core` ahead of `master` (oldest first):

| SHA | Subject |
|---|---|
| `523cc323` | Anatomy: Add icosphere (FreeSurfer/MNE-style) per-hemisphere cortex downsampling |
| `acd1f5cd` | IO: FreeSurfer import cortex downsampling method option (icosphere, default for auto-import) |
| `64f6d56d` | Bugfix: Comment out unused sun.misc.BASE64Decoder import (unresolvable on Java 11+ / MATLAB R2023b+) |
| `bc480ca6` | IO: BIDS import cortex downsampling method and resolution options (default icosphere ico5) |

**Why the 4th commit (`64f6d56d`) was needed:** `toolbox/io/import_label.m` had a top-level `import sun.misc.BASE64Decoder;` statement left over from older Java runtimes. `sun.misc.BASE64Decoder` is an internal Sun class removed from the public API in Java 9+; on MATLAB R2023b (bundled JRE is Java 11+) the `import` statement itself throws, making `import_label.m` fail to parse at all. FreeSurfer anatomy import calls into the atlas/label-import path, so without this fix `test_import_fs_ico` (Task 4) and the FS-anatomy stage of `test_import_bids_ico` (Task 5) are fatal on R2023b, independent of the icosphere work. This mirrors the same fix already applied on `dev` at `fa2d8b2e`; commenting out the unused import (the decoder itself is not called anywhere in the file) is the minimal, no-behavior-change fix. Not an icosphere-feature commit, but required for the harness (and any FS-atlas import) to run at all on this MATLAB version.

## 2. Test re-run results (this session, clean re-run of the full harness)

Runner: `dev/verify/phase0/run_matlab.sh` — MATLAB R2023b `-batch`, isolated `user.home` override (`bst_userdir_clean/`) so `brainstorm('server','local')` creates/uses its own `.brainstorm/local_db`, fully isolated from the developer's real `~/.brainstorm`. Each test executed against the clean worktree (`brainstorm3-clean`) via `cd($WT)` inside the MATLAB session, with the harness (`tess_repair.m` oracle) added to path.

| # | Test | Result | Wall time | Log |
|---|---|---|---|---|
| 1 | `test_tess_repair_unit` | PASSED | ~12s (16:10:12–16:10:24) | `test_tess_repair_unit.log` |
| 2 | `test_tess_downsize_ico` | PASSED | ~13s (16:10:29–16:10:42) | `test_tess_downsize_ico.log` |
| 3 | `test_import_fs_ico` | PASSED | ~70s (16:10:46–16:11:56) | `test_import_fs_ico.log` |
| 4 | `test_import_bids_ico` | PASSED | ~67s (16:12:01–16:13:08) | `test_import_bids_ico.log` |

`grep -l "PASSED" dev/verify/phase0/*.log` matches all four logs.

### Test 1 — `test_tess_repair_unit` (oracle self-check)
Synthetic closed sphere (642 verts): clean mesh validates manifold; a duplicated flipped-winding face is correctly detected as non-manifold; `tess_repair(...,'Repair',1)` removes exactly the spurious face and restores manifoldness; repair on an already-clean mesh is a no-op (`isequal` on vertices/faces). All four assertions passed.

### Test 2 — `test_tess_downsize_ico` (icosphere downsampler, synthetic)
Synthetic registered sphere (40962 verts, self-registration) downsampled via `tess_downsize(..., 10242, 'icosphere')`. Log line: `Icosphere resampling: requested 10242, snapped to ico grid of 10242 vertices.` Result: 10242 vertices (exact ico grid snap), injective + sorted vertex mapping `I`, `tess_repair` confirms closed 2-manifold, winding (signed volume sign) matches source.

### Test 3 — `test_import_fs_ico` (FreeSurfer import, live dataset)
Non-interactive `import_anatomy_fs` on `sub-0002` FreeSurfer output (`derivatives/freesurfer/sub-0002/ses-mri/anat`, lh/rh pial+white+sphere.reg, ~153k verts/hemi native). Default method resolves to `'icosphere'`. Resulting cortex: **20484 vertices** (`ico5`, asserted exactly — `assert(size(Cortex.Vertices,1) == 20484, ...)` did not fire). Per-hemisphere 2-manifold check via `tess_hemisplit` + `tess_repair` oracle passed for both hemispheres (no assertion failure). Protocol `FsIcoUnit` deleted at end of test (scratch protocol, not the kept one).

### Test 4 — `test_import_bids_ico` (BIDS import, live dataset, DEFAULT options)
`process_import_bids` run on the OMEGA tutorial BIDS root with **no explicit `downsamplemethod`/`icolevel`** — exercising the `GetDescription` default (icosphere / ico5). Protocol `omega-tutorial-cortical-flow`, subject `sub-0002`. Result (from log):

```
Cortex: sub-0002
  vertices=20484 (expect 20484), faces=40960, components=2 (expect 2), manifold=1 (expect 1)
ALL TESTS PASSED: test_import_bids_ico
```

**ico5 vertex count: 20484** (10242/hemisphere), 40960 faces, 2 connected components, both hemispheres pass the `tess_hemisplit` + `tess_repair` 2-manifold oracle check. Protocol intentionally **kept** (not deleted) for Gate 0 GUI inspection, per the test's own design.

## 3. Known pre-existing limitation (not introduced by this branch)

During the BIDS import run, the log shows:

```
Couldn't open /Volumes/SpikeData-2/workspace/library/datasets/omega-tutorial/sub-0002/ses-01/anat/._sub-0002_ses-01_T1w.nii.gz.
java.util.zip.ZipException: Not in GZIP format
	at java.base/java.util.zip.GZIPInputStream.readHeader(GZIPInputStream.java:166)
	...
	at org.brainstorm.file.Unpack.gunzip(Unpack.java:19)
```

The dataset directory contains macOS AppleDouble sidecar files (e.g. `._sub-0002_ses-01_T1w.nii.gz`, visible via `ls -la` on the dataset root — `._sub-0002`, `._participants.tsv`, etc.). `process_import_bids`'s upstream directory scan picks up the AppleDouble file alongside the real `.nii.gz`, tries to gunzip it, and the scan aborts with this error **after anatomy import completes but before MEG recording linking**. Consequence: the kept `omega-tutorial-cortical-flow` protocol contains **anatomy only** for `sub-0002` — confirmed on disk: `anat/sub-0002` (surfaces present) but `data/sub-0002` contains only `@intra`/`@default_study` (no raw/recording session folders). The importer nonetheless reports `nVertices=20484`, 2 components, both hemispheres manifold — i.e. the icosphere cortex import itself completed correctly; the abort happens in a later, unrelated stage (MEG linking) of the same BIDS scan.

This is **pre-existing upstream `process_import_bids` behavior** (AppleDouble files are a macOS filesystem artifact of the dataset copy, not something this branch's icosphere changes touch) and is **not introduced by the icosphere clean-branch work**. It does not block Gate 0 — the vertex/manifold checks that Gate 0 cares about (ico5, per-hemisphere 2-manifold cortex from both the FreeSurfer path and the BIDS-default path) both passed. Flagging it here so Phase 1 planning is not surprised that the kept protocol has no MEG recordings.

## 4. Deviations from the brief

- Test filenames: the brief referred to a "live" BIDS test; the actual script is `test_import_bids_ico.m` (not `_live`), per the task instructions' correction.
- No rejected hunks or manual conflict resolution was needed this run — reran the harness clean end-to-end, no code changes to the clean worktree were required (worktree stayed at `bc480ca6`, confirmed clean before and after).
- No standalone-FS-folder skip needed — `test_import_fs_ico` ran directly against the live dataset's FreeSurfer derivatives folder and passed.
- Report additionally documents the AppleDouble/BIDS-scan-abort limitation (see §3) and the 4th commit's rationale (see §1), per corrected task instructions.
- A `launch_gui_clean.sh` script was added (not in the original brief, which had two inaccuracies in its GUI launch instructions — see the accompanying task report) to give a working, isolation-safe way to inspect the kept protocol in the desktop GUI. **Untested, mirrors tested runner isolation**: this session did not launch it live. `run_matlab.sh`'s isolation mechanism (java `user.home` override -> isolated `.brainstorm/local_db`) is confirmed working by all four headless test runs above; `launch_gui_clean.sh` applies the identical override with `brainstorm('start','local')` (the GUI-mode analog of the `brainstorm('server','local')` call `run_matlab.sh` uses) and `matlab -desktop -r` instead of `-batch`. Spawning a real interactive MATLAB desktop + Brainstorm GUI window on the active console session from an unattended agent run was judged out of scope / disruptive for this task, so only the shell script's syntax was checked (`bash -n`), not an actual GUI launch. Recommend a human runs it once before relying on it.

## 5. Files in this harness commit

- `run_matlab.sh` — isolated headless MATLAB runner (`-batch`) against the clean worktree
- `launch_gui_clean.sh` — isolated desktop MATLAB launcher (`-desktop`) for GUI inspection of the kept protocol
- `tess_repair.m` — manifold-validation/repair oracle (harness-local, added to path, not part of the clean branch)
- `test_tess_repair_unit.m`, `test_tess_downsize_ico.m`, `test_import_fs_ico.m`, `test_import_bids_ico.m` — the four Gate 0 tests
- `*.log` — this run's captured output for all four tests
- `.gitignore` — excludes `bst_userdir_clean/` (isolated Brainstorm user dir / DB, not checked in)
