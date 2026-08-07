# Phase 0 — Icosphere Clean Branch Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create clean branch `feature/cortical-flow-core` off upstream master and re-derive the icosphere discretization prerequisite (tess_repair, tess_downsize 'icosphere', FreeSurfer/BIDS import wiring) as self-contained commits, verified by Gate 0.

**Architecture:** Extraction, not new development — apply the hunks of specific
`development`-branch commits onto fresh upstream files in a dedicated git
worktree. The math/code already exists and is validated; the work is clean
re-derivation plus per-task verification.

**Tech Stack:** MATLAB R2023b (`matlab -batch`), Brainstorm (no plugins), git worktree.

**Spec:** `docs/superpowers/specs/2026-08-07-cortical-flow-distillation-design.md`

## Global Constraints

- Clean-branch commits carry NO Claude co-author trailers and NO `dev/*.md`
  docs/plans. Commit messages follow upstream Brainstorm style
  (`Anatomy: ...`, `IO: ...`).
- NEVER switch branches in the main checkout
  (`~/workspace/research/code/brainstorm3` — the user's active work lives
  there). All clean-branch work happens in the worktree
  `~/workspace/research/code/brainstorm3-clean` (`$WT` below).
- Pure MATLAB only on the clean branch: no nxr-compute, no bst-java, no DB
  schema changes.
- Verification harness lives in the MAIN checkout under `dev/verify/phase0/`
  (committed to the lab branch, never to the clean branch).
- MATLAB runs headless: `/Applications/MATLAB_R2023b.app/bin/matlab -batch`
  (adjust path if `which matlab` resolves elsewhere). Never use `clear` in a
  live shared MATLAB session; batch runs are isolated so it doesn't apply.
- All Brainstorm tests use an ISOLATED database dir
  (`bst_set('BrainstormDbDir', <scratch>/bst_db_clean)`) and disposable
  protocols deleted on success — never the user's real brainstorm_db.
- Source-of-truth commits on `development` (main checkout):
  `82abcc00` (tess_downsize icosphere + import_anatomy_fs Method arg),
  `894f8c74` (import_anatomy auto-import wiring),
  `tess_repair.m` current file on `development` (consolidation lineage
  `6c0f5631` → rename),
  BIDS wiring: `99f54332`, `5462916e`, `0524f046`, `4c3f9b1e`, `ec26e190`,
  `00f2ec90` (file `process_import_bids.m` only).
- Gate 0 dataset: `/Users/diellorbasha/workspace/library/datasets/omega-tutorial`
  (BIDS, subject `sub-0002`).

---

### Task 1: Branch + worktree + harness scaffold

**Files:**
- Create: worktree `~/workspace/research/code/brainstorm3-clean` on new branch `feature/cortical-flow-core`
- Create (main checkout): `dev/verify/phase0/run_matlab.sh`

**Interfaces:**
- Produces: `$WT` = worktree path used by every later task; `run_matlab.sh <script.m>` runs a harness script against the worktree's Brainstorm in an isolated DB.

- [ ] **Step 1: Update local master to upstream (no checkout switch)**

```bash
cd ~/workspace/research/code/brainstorm3
git fetch upstream
git fetch upstream master:master     # ff-only update of the ref; fails loudly if not ff
git log --oneline -1 master          # note the new base SHA
```

- [ ] **Step 2: Create the worktree + branch**

```bash
git worktree add ~/workspace/research/code/brainstorm3-clean -b feature/cortical-flow-core master
ls ~/workspace/research/code/brainstorm3-clean/brainstorm.m   # sanity: worktree populated
```

- [ ] **Step 3: Create the MATLAB runner for harness scripts**

Create `dev/verify/phase0/run_matlab.sh` (main checkout):

```bash
#!/bin/bash
# Run a harness .m script against the CLEAN worktree's Brainstorm, isolated DB.
# Usage: ./run_matlab.sh /abs/path/to/script.m
set -euo pipefail
WT="$HOME/workspace/research/code/brainstorm3-clean"
SCRIPT="$1"
DBDIR="${BST_DB_CLEAN:-$HOME/workspace/research/code/brainstorm3/dev/verify/phase0/bst_db_clean}"
mkdir -p "$DBDIR"
MATLAB="/Applications/MATLAB_R2023b.app/bin/matlab"
[ -x "$MATLAB" ] || MATLAB="$(command -v matlab)"
"$MATLAB" -batch "cd('$WT'); brainstorm server; bst_set('BrainstormDbDir','$DBDIR'); db_import('$DBDIR'); run('$SCRIPT'); brainstorm stop; exit(0);"
```

`chmod +x dev/verify/phase0/run_matlab.sh`. Note: if `db_import` errors on an
empty dir on first run, drop that call — `bst_set('BrainstormDbDir',...)`
followed by protocol creation is sufficient; verify on first use in Task 2.

- [ ] **Step 4: Verify MATLAB launches and Brainstorm server starts in the worktree**

```bash
echo "disp(['BST OK, ver: ' bst_get('Version').Version])" > /tmp/bst_smoke.m  # any scratch path
dev/verify/phase0/run_matlab.sh /tmp/bst_smoke.m
```

Expected: prints `BST OK, ver: ...` with no dialog/hang. (No commit this task; the harness is committed in Task 6.)

---

### Task 2: `tess_repair` on the clean branch

**Files:**
- Create (worktree): `toolbox/anatomy/tess_repair.m`
- Create (main checkout): `dev/verify/phase0/test_tess_repair_unit.m`

**Interfaces:**
- Produces: `[Vertices, Faces, isManifold, report] = tess_repair(Vertices, Faces, 'Repair',0/1, 'RequireClosed',0/1)` — used by Tasks 3–6 to assert manifoldness.

- [ ] **Step 1: Write the unit test (main checkout)**

Create `dev/verify/phase0/test_tess_repair_unit.m`:

```matlab
% TEST_TESS_REPAIR_UNIT: validation + repair on a synthetic closed mesh.
[V, F] = tess_sphere(642);
% 1) clean closed sphere validates
[~, ~, isM, rep] = tess_repair(V, F);
assert(isM, 'clean sphere must validate as manifold');
% 2) duplicated face (flipped winding) -> non-manifold edges detected
F2 = [F; F(1, [2 1 3])];
[~, ~, isM2, rep2] = tess_repair(V, F2);
assert(~isM2, 'corrupted mesh must fail validation');
% 3) repair removes the spurious face and restores manifoldness
[Vr, Fr, isM3, rep3] = tess_repair(V, F2, 'Repair', 1);
assert(isM3, 'repair must restore manifoldness');
assert(size(Fr,1) == size(F,1), 'repair should remove exactly the spurious face');
% 4) repairing an already-clean mesh is a no-op
[Vn, Fn, isM4] = tess_repair(V, F, 'Repair', 1);
assert(isM4 && isequal(Fn, F) && isequal(Vn, V), 'repair of clean mesh must be a no-op');
disp('test_tess_repair_unit PASSED');
```

- [ ] **Step 2: Run it — must FAIL (function absent in worktree)**

```bash
dev/verify/phase0/run_matlab.sh ~/workspace/research/code/brainstorm3/dev/verify/phase0/test_tess_repair_unit.m
```

Expected: error `Unrecognized function ... 'tess_repair'`.

- [ ] **Step 3: Extract the file from development**

```bash
cd ~/workspace/research/code/brainstorm3
git show development:toolbox/anatomy/tess_repair.m > ~/workspace/research/code/brainstorm3-clean/toolbox/anatomy/tess_repair.m
```

Then open the worktree copy and confirm: header comment mentions no nxr, no
DB, no tess_manifold; authorship line reads Diellor Basha. It is a pure
function — no edits expected.

- [ ] **Step 4: Run the test — must PASS**

Same command as Step 2. Expected: `test_tess_repair_unit PASSED`. If
assertion 3 fails on face count, inspect `rep3` — the repair criterion is
"remove the face whose normal deviates most from neighbors"; the duplicated
face qualifies. Debug before proceeding (systematic-debugging), do not relax
the assertion.

- [ ] **Step 5: Commit (worktree, clean style, no trailer)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_repair.m
git commit -m "Anatomy: Add tess_repair, 2-manifold surface validation and repair"
```

---

### Task 3: `tess_downsize` 'icosphere' method

**Files:**
- Modify (worktree): `toolbox/anatomy/tess_downsize.m`
- Create (main checkout): `dev/verify/phase0/test_tess_downsize_ico.m`

**Interfaces:**
- Consumes: `tess_repair` (Task 2), upstream `tess_sphere` (unchanged — the dev-branch tess_sphere diff is upstream drift, NOT ours; do not touch it).
- Produces: `tess_downsize(TessFile, N, 'icosphere')` accepting the new method string; ico grid snap (642/2562/10242/40962), injective vertex mapping, local helper `resolve_ico_collisions`.

- [ ] **Step 1: Write the unit test (main checkout)**

Create `dev/verify/phase0/test_tess_downsize_ico.m`:

```matlab
% TEST_TESS_DOWNSIZE_ICO: icosphere downsampling on a synthetic registered sphere.
Protocol = 'IcoDownsizeUnit';
gui_brainstorm('DeleteProtocol', Protocol);
gui_brainstorm('CreateProtocol', Protocol, 0, 0);
[~, iSubject] = db_add_subject('TestSubj', [], 0, 0);
% Synthetic closed surface WITH registration sphere (the sphere is its own reg)
[V, F] = tess_sphere(40962);
TessMat = db_template('surfacemat');
TessMat.Comment  = 'sphere_40962';
TessMat.Vertices = V;
TessMat.Faces    = F;
TessMat.Reg.Sphere.Vertices = V;
ProtocolInfo = bst_get('ProtocolInfo');
TessFile = bst_fullfile(ProtocolInfo.SUBJECTS, 'TestSubj', 'tess_sphere_40962.mat');
bst_save(TessFile, TessMat, 'v7');
db_add_surface(iSubject, TessFile, TessMat.Comment);
% Downsample to the 10242 ico grid
[NewFile, iSurf, I, J] = tess_downsize(TessFile, 10242, 'icosphere');
NewMat = in_tess_bst(NewFile, 0);
assert(size(NewMat.Vertices,1) == 10242, 'expected snap to ico grid 10242');
assert(numel(unique(I)) == numel(I), 'vertex mapping must be injective');
assert(issorted(I), 'kept-vertex indices must be sorted (reducepatch convention)');
[~, ~, isM] = tess_repair(NewMat.Vertices, NewMat.Faces);
assert(isM, 'icosphere output must be a closed 2-manifold');
% Winding must match the source (signed volume same sign)
sv = @(Vv,Ff) sum(sum(Vv(Ff(:,1),:) .* cross(Vv(Ff(:,2),:), Vv(Ff(:,3),:), 2)));
assert(sign(sv(NewMat.Vertices,NewMat.Faces)) == sign(sv(V,F)), 'winding must match source');
gui_brainstorm('DeleteProtocol', Protocol);
disp('test_tess_downsize_ico PASSED');
```

- [ ] **Step 2: Run it — must FAIL**

```bash
dev/verify/phase0/run_matlab.sh ~/workspace/research/code/brainstorm3/dev/verify/phase0/test_tess_downsize_ico.m
```

Expected: failure at the `tess_downsize(..., 'icosphere')` call (unknown
method / case not handled). If it fails EARLIER (protocol/API calls), fix the
test's Brainstorm calls first — the API idioms above follow
`dev/tests/test_import_bids_ico_live.m` on development.

- [ ] **Step 3: Apply the icosphere hunks from commit 82abcc00**

```bash
cd ~/workspace/research/code/brainstorm3
git show 82abcc00 -- toolbox/anatomy/tess_downsize.m > /tmp/ico_downsize.patch   # scratch path
cd ~/workspace/research/code/brainstorm3-clean
git apply --3way /tmp/ico_downsize.patch || git apply --reject /tmp/ico_downsize.patch
```

If hunks are rejected (upstream drift), apply them manually: the commit adds
(1) `'icosphere'` to the Method doc-list and interactive method dialog,
(2) the `case 'icosphere'` block, (3) the `resolve_ico_collisions` local
function. Compare against `git show development:toolbox/anatomy/tess_downsize.m`
for the authoritative current form of these three pieces. Do NOT copy the
whole dev file (it may contain unrelated drift); add only these pieces.

- [ ] **Step 4: Run the test — must PASS**

Same command as Step 2. Expected: `test_tess_downsize_ico PASSED`.

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/anatomy/tess_downsize.m
git commit -m "Anatomy: Add icosphere (FreeSurfer/MNE-style) per-hemisphere cortex downsampling"
```

---

### Task 4: FreeSurfer import wiring (`import_anatomy_fs` + `import_anatomy`)

**Files:**
- Modify (worktree): `toolbox/io/import_anatomy_fs.m`, `toolbox/io/import_anatomy.m`
- Create (main checkout): `dev/verify/phase0/test_import_fs_ico.m`

**Interfaces:**
- Consumes: `tess_downsize(..., 'icosphere')` (Task 3).
- Produces: `import_anatomy_fs(iSubject, FsDir, nVertices, isInteractive, sFid, isExtraMaps, isVolumeAtlas, isKeepMri, Method)` — 9th arg `Method` (`'icosphere'`|`'reducepatch'`|`[]`=ask); auto-import defaults `Method='icosphere'`, `nVertices=[]` → ico5 (10242/hemi, 20484 total). Task 5's BIDS path calls this signature.

- [ ] **Step 1: Locate a FreeSurfer folder for the live check**

```bash
find /Users/diellorbasha/workspace/library/datasets/omega-tutorial -maxdepth 4 -name "surf" -type d | head -3
```

Expected: a `derivatives/freesurfer*/sub-0002*/surf` path; note it as `$FSDIR`
(the folder CONTAINING surf/, i.e. the subject's FS directory). If none
exists, skip this task's live check (Steps 2/5 assertions) and rely on Task
5's BIDS live check, which exercises the same code path — note the skip in
the Gate 0 report.

- [ ] **Step 2: Write the live check (main checkout)**

Create `dev/verify/phase0/test_import_fs_ico.m` (replace `<FSDIR>` with the
path found in Step 1):

```matlab
% TEST_IMPORT_FS_ICO: non-interactive FS import must yield an ico5 manifold cortex.
FsDir = '<FSDIR>';
Protocol = 'FsIcoUnit';
gui_brainstorm('DeleteProtocol', Protocol);
gui_brainstorm('CreateProtocol', Protocol, 0, 0);
[~, iSubject] = db_add_subject('FsSubj', [], 0, 0);
errorMsg = import_anatomy_fs(iSubject, FsDir, [], 0, [], 0, 1, 0, 'icosphere');
assert(isempty(errorMsg), 'import_anatomy_fs error: %s', errorMsg);
sSubject = bst_get('Subject', iSubject);
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
Cortex = in_tess_bst(CortexFile, 0);
assert(size(Cortex.Vertices,1) == 20484, 'expected ico5 cortex (20484), got %d', size(Cortex.Vertices,1));
% Per-hemisphere manifold check: split on the Structures atlas, never conncomp
[rH, lH] = tess_hemisplit(in_tess_bst(CortexFile));
for hemi = {rH, lH}
    iV = hemi{1};
    iF = all(ismember(Cortex.Faces, iV), 2);
    lut = zeros(size(Cortex.Vertices,1),1); lut(iV) = 1:numel(iV);
    [~, ~, isM] = tess_repair(Cortex.Vertices(iV,:), lut(Cortex.Faces(iF,:)));
    assert(isM, 'hemisphere fails 2-manifold check');
end
gui_brainstorm('DeleteProtocol', Protocol);
disp('test_import_fs_ico PASSED');
```

Note: verify `tess_hemisplit`'s exact return signature on the worktree
(`help tess_hemisplit`) before running — it is an upstream function; if it
returns `[iRH, iLH, isConnected]` adjust the destructuring accordingly.

- [ ] **Step 3: Apply the wiring hunks**

```bash
cd ~/workspace/research/code/brainstorm3
git show 82abcc00 -- toolbox/io/import_anatomy_fs.m > /tmp/fs_method.patch
git show 894f8c74 -- toolbox/io/import_anatomy.m   > /tmp/auto_import.patch
cd ~/workspace/research/code/brainstorm3-clean
git apply --3way /tmp/fs_method.patch  || git apply --reject /tmp/fs_method.patch
git apply --3way /tmp/auto_import.patch || git apply --reject /tmp/auto_import.patch
```

On rejects, apply manually with `git show development:<file>` as the
authoritative form. The pieces: `Method` 9th argument + interactive
method/resolution dialog + icosphere branch in `import_anatomy_fs`;
`FsMethod='icosphere'` auto-import defaults + updated `import_anatomy_fs`
call sites in `import_anatomy`. Do NOT bring dev-only drift (e.g. the removed
overwrite-confirmation dialog in dev's import_anatomy is unrelated upstream
divergence — keep upstream's version of everything outside these pieces).

- [ ] **Step 4: Run the live check — must PASS**

```bash
dev/verify/phase0/run_matlab.sh ~/workspace/research/code/brainstorm3/dev/verify/phase0/test_import_fs_ico.m
```

Expected: `test_import_fs_ico PASSED` (runtime: several minutes — full FS import).

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/io/import_anatomy_fs.m toolbox/io/import_anatomy.m
git commit -m "IO: FreeSurfer import cortex downsampling method option (icosphere, default for auto-import)"
```

---

### Task 5: BIDS import wiring (`process_import_bids`)

**Files:**
- Modify (worktree): `toolbox/process/functions/process_import_bids.m`
- Create (main checkout): `dev/verify/phase0/test_import_bids_ico.m`

**Interfaces:**
- Consumes: `import_anatomy_fs` 9-arg signature (Task 4).
- Produces: `process_import_bids` options `downsamplemethod` (default `'icosphere'`) + ico-level selector (default ico5) + `ResolveAnatDownsample` helper threading method/level to the anatomy import.

- [ ] **Step 1: Copy the proven live test (main checkout)**

```bash
cd ~/workspace/research/code/brainstorm3
git show development:dev/tests/test_import_bids_ico_live.m > dev/verify/phase0/test_import_bids_ico.m
```

Open it and adjust ONLY: (a) function name to `test_import_bids_ico`; (b) it
already targets `/Users/diellorbasha/workspace/library/datasets/omega-tutorial`
sub-0002 with DEFAULT options and asserts an ico5 (20484) manifold cortex —
keep those assertions verbatim; (c) if it calls dev-only helpers, inline the
`bst_get('Subject')`/`iCortex` idiom from Task 4's test instead.

- [ ] **Step 2: Run it — must FAIL**

```bash
dev/verify/phase0/run_matlab.sh ~/workspace/research/code/brainstorm3/dev/verify/phase0/test_import_bids_ico.m
```

Expected: failure because `process_import_bids` has no `downsamplemethod`
option yet (or default reducepatch cortex ≠ 20484).

- [ ] **Step 3: Apply the combined BIDS wiring diff**

```bash
cd ~/workspace/research/code/brainstorm3
git diff 894f8c74..00f2ec90 -- toolbox/process/functions/process_import_bids.m > /tmp/bids_wiring.patch
cd ~/workspace/research/code/brainstorm3-clean
git apply --3way /tmp/bids_wiring.patch || git apply --reject /tmp/bids_wiring.patch
```

Upstream has drifted on this file — expect partial rejects. Resolve manually;
the semantic pieces (authoritative form:
`git show development:toolbox/process/functions/process_import_bids.m`):
(1) `GetDescription` options `downsamplemethod` (combobox, default icosphere)
and ico-level (default ico5, self-documenting labels); (2) local function
`ResolveAnatDownsample` (FreeSurfer-only icosphere, reducepatch fallback with
warning, method-aware vertex validation); (3) `Run`/call-site threading of
method+level into `import_anatomy_fs(..., Method)`. Bring nothing else from
dev's version of this file.

- [ ] **Step 4: Run the live test — must PASS (this is Gate 0's core check)**

Same command as Step 2. Expected: final `PASSED` line with ico5 20484 +
per-hemisphere manifold assertions green. Runtime: several minutes.

- [ ] **Step 5: Commit (worktree)**

```bash
cd ~/workspace/research/code/brainstorm3-clean
git add toolbox/process/functions/process_import_bids.m
git commit -m "IO: BIDS import cortex downsampling method and resolution options (default icosphere ico5)"
```

---

### Task 6: Gate 0 report + harness commit + user review

**Files:**
- Create (main checkout): `dev/verify/phase0/gate0_report.md`
- Commit (main checkout, lab branch): `dev/verify/phase0/*`

**Interfaces:**
- Consumes: all Task 2–5 tests.
- Produces: the Gate 0 evidence document the user reviews before Phase 1 begins.

- [ ] **Step 1: Re-run the full harness clean, capturing output**

```bash
cd ~/workspace/research/code/brainstorm3
for t in test_tess_repair_unit test_tess_downsize_ico test_import_fs_ico test_import_bids_ico; do
  dev/verify/phase0/run_matlab.sh "$PWD/dev/verify/phase0/$t.m" 2>&1 | tee "dev/verify/phase0/$t.log"
done
grep -l "PASSED" dev/verify/phase0/*.log
```

Expected: all four logs contain their `PASSED` line (or three, with the Task 4
skip documented if no standalone FS folder existed).

- [ ] **Step 2: Write `dev/verify/phase0/gate0_report.md`**

Contents (concrete, from the logs — no prose padding): clean-branch base SHA
and the 4 commit SHAs with subjects; per-test one-line result + runtime;
the ico5 vertex count and per-hemisphere manifold verdicts from the live
imports; any deviations (rejected hunks resolved manually, skipped FS check).

- [ ] **Step 3: Commit the harness on the lab branch (main checkout)**

```bash
cd ~/workspace/research/code/brainstorm3
git add dev/verify/phase0/
git commit -m "verify(phase0): Gate 0 harness + report for icosphere clean-branch extraction

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

(Trailer is fine here — this is the lab branch, not the clean branch.)

- [ ] **Step 4: Present Gate 0 to the user**

Show: the report, `git -C ~/workspace/research/code/brainstorm3-clean log --oneline master..` (expect 4 clean commits), and ask for the Gate 0 verdict. Phase 1 planning starts only after approval. Do NOT push the clean branch anywhere until the user says so.
