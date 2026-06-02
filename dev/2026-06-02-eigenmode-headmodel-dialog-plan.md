# Eigenmode Head-Model Dialog Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add the eigenmode leadfield to the "Compute head model" dialog as a Source-space option ("Cortex surface harmonics"), composing `L̃ = L·Φ` on an in-memory base leadfield and saving only the harmonic head-model node.

**Architecture:** All changes live in `toolbox/forward/panel_headmodel.m`. The `ComputeHeadModel` branch is the testable logic (computes the base via `bst_headmodeler` in-memory, composes via the merged `bst_eigenmode_leadfield`, saves only the harmonic node) — covered by an e2e that calls `ComputeHeadModel` with a hand-built `sMethod`, bypassing the GUI. The GUI controls (a 4th Source-space radio + a modes field) are wired separately and validated by `checkcode` + manual launch.

**Tech Stack:** MATLAB / Brainstorm. Tests are `dev/tests/*.m` functions printing `ALL TESTS PASSED`, run via the MATLAB MCP (`run_matlab_file`). Spec: `dev/2026-06-02-eigenmode-headmodel-dialog-design.md`.

---

## File structure

| File | Action | Responsibility |
|---|---|---|
| `toolbox/forward/panel_headmodel.m` | Modify | Source-space radio + modes field (GUI); `GetPanelContents` mapping; `ComputeHeadModel` eigenspace branch (compose + save-only-harmonic). |
| `toolbox/forward/bst_eigenmode_leadfield.m` | Reuse unchanged | Composition engine. |
| `toolbox/io/in_tess_eigenmodes.m` | Reuse unchanged | Eigenmode load + availability check. |
| `dev/tests/test_headmodel_dialog_eigenmode_e2e.m` | Create | e2e for the `ComputeHeadModel` eigenspace branch. |

---

## Task 1: `ComputeHeadModel` eigenspace branch (logic + e2e)

This is the real logic and is testable without the GUI by passing `sMethod` directly.

**Files:**
- Modify: `toolbox/forward/panel_headmodel.m` (`ComputeHeadModel`)
- Test: `dev/tests/test_headmodel_dialog_eigenmode_e2e.m`

- [ ] **Step 1: Write the failing e2e test** — create `dev/tests/test_headmodel_dialog_eigenmode_e2e.m`:

```matlab
function test_headmodel_dialog_eigenmode_e2e
% Smoke: ComputeHeadModel with the "Cortex surface harmonics" source space
% (passed as sMethod, bypassing the GUI) produces exactly ONE harmonic head
% model node (Gain=[nCh x K], isEigenmode, Comment '... | harmonic') and no
% extra base node. Skips cleanly without a suitable study.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study: MEG channel + subject cortex with precomputed eigenmodes
iStudyTarget = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.Channel) || isempty(s.Channel.FileName); continue; end
    if ~any(strcmpi(s.Channel.DisplayableSensorTypes, 'MEG')); continue; end
    sSubj = bst_get('Subject', s.BrainStormSubject);
    if isempty(sSubj) || isempty(sSubj.iCortex); continue; end
    [~, isEig] = in_tess_eigenmodes(sSubj.Surface(sSubj.iCortex).FileName);
    if isEig; iStudyTarget = iS; break; end
end
if isempty(iStudyTarget)
    disp('SKIP: no study with MEG channel + cortex eigenmodes.'); return;
end

sBefore   = bst_get('Study', iStudyTarget);
nHmBefore = numel(sBefore.HeadModel);

% Build a harmonic sMethod (bypass the GUI)
sMethod = struct('Comment', 'Overlapping spheres | harmonic', ...
    'HeadModelType', 'surface', 'SourceCompression', 'eigenmode', 'nModes', 0, ...
    'MEGMethod', 'os_meg', 'EEGMethod', '', 'ECOGMethod', '', 'SEEGMethod', '', ...
    'NIRSMethod', '', 'SaveFile', 1, 'Interactive', 0);
OutFiles = panel_headmodel('ComputeHeadModel', iStudyTarget, sMethod);
assert(~isempty(OutFiles), 'Harmonic head model must be produced.');

% Exactly one new node, and it is the harmonic one
sAfter = bst_get('Study', iStudyTarget);
assert(numel(sAfter.HeadModel) == nHmBefore + 1, 'Exactly one new head model node expected (no base node).');
HM = in_bst_headmodel(sAfter.HeadModel(sAfter.iHeadModel).FileName, 0);
assert(isfield(HM,'isEigenmode') && HM.isEigenmode==1, 'Active node must be the harmonic head model.');
assert(size(HM.Gain,2) == HM.nModes, 'Gain must be [nCh x nModes].');
assert(~isempty(HM.Eigenvalues), 'Eigenvalues must be stored.');
assert(~isempty(regexp(HM.Comment, '\| harmonic$', 'once')), 'Comment must end with "| harmonic".');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `dev/tests/test_headmodel_dialog_eigenmode_e2e.m`. Expected: an error/assert failure because `ComputeHeadModel` does not yet honor `SourceCompression='eigenmode'` (it will try to save a normal base node and won't set `isEigenmode`). If it SKIPs, develop against OMEGA sub-0002 (study with cortex eigenmodes + MEG).

- [ ] **Step 3a: Add the eigenspace flag.** In `toolbox/forward/panel_headmodel.m`, find:

```matlab
    isOpenMEEG = any(strcmpi(allMethods, 'openmeeg'));
    isDuneuro = any(strcmpi(allMethods, 'duneuro'));
```
and insert immediately after:
```matlab
    % "Cortex surface harmonics" source space: compose L*Phi on an in-memory base
    % leadfield and save ONLY the harmonic node (not a separate base node).
    isEigenSpace = isfield(sMethod, 'SourceCompression') && strcmpi(sMethod.SourceCompression, 'eigenmode');
```

- [ ] **Step 3b: Validate eigenmodes early + force in-memory base.** Find (inside the study loop):

```matlab
        % Load channel description
        ChannelFile = file_fullpath(sStudy.Channel.FileName);
        ChannelMat = in_bst_channel(ChannelFile);
```
and insert immediately BEFORE it:
```matlab
        % Eigenmode source space requires precomputed eigenmodes on the cortex
        if isEigenSpace
            [~, isEig] = in_tess_eigenmodes(sSubject.Surface(sSubject.iCortex).FileName);
            if ~isEig
                errMessage = 'No eigenmodes on the cortex surface. Compute eigenmodes first.';
                continue;
            end
        end
```
Then find:
```matlab
        if sMethod.SaveFile
            OPTIONS.HeadModelFile = bst_fileparts(ChannelFile);
        else
            OPTIONS.HeadModelFile = '';
        end
```
and replace with:
```matlab
        if sMethod.SaveFile && ~isEigenSpace
            OPTIONS.HeadModelFile = bst_fileparts(ChannelFile);
        else
            OPTIONS.HeadModelFile = '';   % eigenspace: base computed in-memory; only the harmonic node is saved
        end
```

- [ ] **Step 3c: Compose + save only the harmonic node.** Find:

```matlab
        % ===== Add new HeadModel in Brainstorm Database =====
        % If a file was saved
        if ~isempty(OPTIONS.HeadModelFile)
```
and insert immediately BEFORE that block:
```matlab
        % ===== EIGENMODE COMPOSITION ("Cortex surface harmonics") =====
        if isEigenSpace
            baseHM = OPTIONS.HeadModelMat;            % in-memory base (no file written)
            [Eig, isEig] = in_tess_eigenmodes(baseHM.SurfaceFile);
            if ~isEig
                errMessage = 'No eigenmodes on the cortex surface. Compute eigenmodes first.';
                continue;
            end
            nModesReq = 0;
            if isfield(sMethod, 'nModes') && ~isempty(sMethod.nModes)
                nModesReq = sMethod.nModes;
            end
            CompHM = bst_eigenmode_leadfield(baseHM, Eig, 'nModes', nModesReq);
            CompHM.Comment = OPTIONS.Comment;         % already "<base> | harmonic" (set by UpdateComment)
            CompHM = bst_history('add', CompHM, 'eigenmode_leadfield', ...
                sprintf('Eigenmode leadfield (%d modes) composed in Compute head model', CompHM.nModes));
            % Save only the harmonic node
            StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
            OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_eigenmode');
            bst_save(OutputFile, CompHM, 'v7');
            newHeadModel = db_template('HeadModel');
            newHeadModel.FileName      = file_short(OutputFile);
            newHeadModel.Comment       = CompHM.Comment;
            newHeadModel.HeadModelType = 'surface';
            newHeadModel.MEGMethod     = OPTIONS.MEGMethod;
            newHeadModel.EEGMethod     = OPTIONS.EEGMethod;
            newHeadModel.ECOGMethod    = OPTIONS.ECOGMethod;
            newHeadModel.SEEGMethod    = OPTIONS.SEEGMethod;
            newHeadModel.NIRSMethod    = OPTIONS.NIRSMethod;
            iHeadModel = length(sStudy.HeadModel) + 1;
            sStudy.HeadModel(iHeadModel) = newHeadModel;
            sStudy.iHeadModel = iHeadModel;
            bst_set('Study', iStudy, sStudy);
            panel_protocols('UpdateNode', 'Study', iStudy);
            OutputFiles{end+1} = file_short(OutputFile);
            continue;     % skip the normal base-node save block
        end

```
(Note: the harmonic comment relies on `OPTIONS.Comment` already ending in `| harmonic`. Until Task 2 wires `UpdateComment`, the e2e passes `sMethod.Comment = 'Overlapping spheres | harmonic'` directly, which flows into `OPTIONS.Comment` via `struct_copy_fields` — so this task is testable on its own.)

- [ ] **Step 4: Run the test; verify it PASSES** — `dev/tests/test_headmodel_dialog_eigenmode_e2e.m` → `ALL TESTS PASSED`. If `bst_eigenmode_leadfield` rejects the in-memory base, inspect `OPTIONS.HeadModelMat` field names (`.Gain` unconstrained `[nCh×3·nVert]`, `.GridOrient`, `.SurfaceFile`, `.HeadModelType`) — the engine needs those. Report any field-name fix.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/forward/panel_headmodel.m dev/tests/test_headmodel_dialog_eigenmode_e2e.m
git commit -m "$(printf 'Head model dialog: eigenmode (Cortex surface harmonics) compose+save branch\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: GUI controls (radio + modes field + mapping)

Wires the dialog so a user can pick the source space. GUI Swing code isn't unit-testable headlessly; verify with `checkcode` + a manual launch.

**Files:**
- Modify: `toolbox/forward/panel_headmodel.m` (panel creation, `ctrl` struct, `UpdateComment`, `GetPanelContents`)

- [ ] **Step 1: Add the radio + modes field.** Find:

```matlab
    jRadioGridMixed   = gui_component('Radio', jPanelSourceSpace, 'br hfill', 'Custom source model', jButtonGroupGridType,  'Estimates sources based on the atlas "Source model"', @UpdateComment);
    if ~isMixed
        jRadioGridMixed.setEnabled(0);
    end
    % Default: surface
    jRadioGridSurface.setSelected(1);
```
and replace with:
```matlab
    jRadioGridMixed   = gui_component('Radio', jPanelSourceSpace, 'br hfill', 'Custom source model', jButtonGroupGridType,  'Estimates sources based on the atlas "Source model"', @UpdateComment);
    if ~isMixed
        jRadioGridMixed.setEnabled(0);
    end
    % Cortex surface harmonics (LBO eigenmode basis of the cortical source space)
    jRadioGridHarmonics = gui_component('Radio', jPanelSourceSpace, 'br hfill', 'Cortex surface harmonics', jButtonGroupGridType, 'Re-expresses the cortical source space in the LBO eigenmode basis (L*Phi)', @UpdateComment);
    gui_component('label', jPanelSourceSpace, 'br', 'Number of modes (0=all): ');
    jTextNModes = gui_component('text', jPanelSourceSpace, 'tab hfill', '0');
    jTextNModes.setEnabled(0);
    % Default: surface
    jRadioGridSurface.setSelected(1);
```

- [ ] **Step 2: Add controls to the `ctrl` struct.** Find:

```matlab
                  'jRadioGridMixed',     jRadioGridMixed, ...
```
and replace with:
```matlab
                  'jRadioGridMixed',     jRadioGridMixed, ...
                  'jRadioGridHarmonics', jRadioGridHarmonics, ...
                  'jTextNModes',         jTextNModes, ...
```

- [ ] **Step 3: Update `UpdateComment`** for the modes-field enable, NIRS gating, and the label. Find:

```matlab
        % Disable NIRS for other than surface head model
        if isNirs
            if jRadioGridVolume.isSelected() || jRadioGridMixed.isSelected()
                jCheckMethodNIRS.setSelected(0);
                jCheckMethodNIRS.setEnabled(0);
            else
                jCheckMethodNIRS.setEnabled(1);
            end
        end
```
and replace with:
```matlab
        % Enable the "Number of modes" field only for the harmonics source space
        jTextNModes.setEnabled(jRadioGridHarmonics.isSelected());
        % Disable NIRS for anything other than a plain cortex-surface head model
        if isNirs
            if jRadioGridVolume.isSelected() || jRadioGridMixed.isSelected() || jRadioGridHarmonics.isSelected()
                jCheckMethodNIRS.setSelected(0);
                jCheckMethodNIRS.setEnabled(0);
            else
                jCheckMethodNIRS.setEnabled(1);
            end
        end
```
Then find:
```matlab
        % Grid type
        if jRadioGridVolume.isSelected()
            Comment = [Comment ' (volume)'];
        elseif jRadioGridMixed.isSelected()
            Comment = [Comment ' (mixed)'];
        else
            %Comment = [Comment ' (cortex)'];
        end
```
and replace with:
```matlab
        % Grid type
        if jRadioGridVolume.isSelected()
            Comment = [Comment ' (volume)'];
        elseif jRadioGridMixed.isSelected()
            Comment = [Comment ' (mixed)'];
        elseif jRadioGridHarmonics.isSelected()
            Comment = [Comment ' | harmonic'];
        else
            %Comment = [Comment ' (cortex)'];
        end
```

- [ ] **Step 4: Map the radio in `GetPanelContents`.** Find:

```matlab
    % Get source space
    if ctrl.jRadioGridSurface.isSelected()
        s.HeadModelType = 'surface';
    elseif ctrl.jRadioGridVolume.isSelected()
        s.HeadModelType = 'volume';
    elseif ctrl.jRadioGridMixed.isSelected()
        s.HeadModelType = 'mixed';
    end
```
and replace with:
```matlab
    % Get source space
    s.SourceCompression = 'none';
    if ctrl.jRadioGridSurface.isSelected()
        s.HeadModelType = 'surface';
    elseif ctrl.jRadioGridVolume.isSelected()
        s.HeadModelType = 'volume';
    elseif ctrl.jRadioGridMixed.isSelected()
        s.HeadModelType = 'mixed';
    elseif ctrl.jRadioGridHarmonics.isSelected()
        s.HeadModelType = 'surface';            % base physics is a surface run
        s.SourceCompression = 'eigenmode';
        nModes = str2double(char(ctrl.jTextNModes.getText()));
        if isnan(nModes) || nModes < 0; nModes = 0; end
        s.nModes = round(nModes);
    end
```

- [ ] **Step 5: Static check.** Load `mcp__plugin_brainstorm-dev_MATLAB__check_matlab_code` via ToolSearch and run `checkcode` on `toolbox/forward/panel_headmodel.m`. Expected: no ERROR-level diagnostics introduced (ignore pre-existing style/info and the standard `varargout`/`eval` idiom warnings; the new `jRadioGridHarmonics`/`jTextNModes` are used in `ctrl`, `UpdateComment`, and `GetPanelContents`).

- [ ] **Step 6: Re-run the Task-1 e2e** — `dev/tests/test_headmodel_dialog_eigenmode_e2e.m` must still print `ALL TESTS PASSED` (the GUI changes must not break the direct-`sMethod` path).

- [ ] **Step 7: Manual GUI verification (report, do not automate).** In a Brainstorm GUI session: right-click a study → *Compute head model* → confirm the **"Cortex surface harmonics"** radio appears in the Source space group, the **Number of modes** field enables only when it's selected, the Comment shows `… | harmonic`, and running it produces one `… | harmonic` head model node. Note this in the report as manually verified (or flag if a GUI session isn't available).

- [ ] **Step 8: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/forward/panel_headmodel.m
git commit -m "$(printf 'Head model dialog: add Cortex surface harmonics source-space control\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Self-review (done)

- **Spec coverage:** radio in Source space (T2 Step 1); modes field (T2 Step 1/3); `GetPanelContents` mapping (T2 Step 4); in-memory base + compose + save-only-harmonic (T1 Step 3a–3c); `<base> | harmonic` label (T2 Step 3 + carried via `OPTIONS.Comment`); eigenmode prerequisite validation (T1 Step 3b/3c); reuse engine (T1 Step 3c). Testing: e2e (T1) + checkcode/manual (T2). All covered.
- **Placeholders:** none — every edit shows exact old→new code.
- **Consistency:** `SourceCompression`/`nModes` fields set in `GetPanelContents` (T2) match what `ComputeHeadModel` reads (T1: `isEigenSpace`, `sMethod.nModes`). The harmonic label is produced once (`UpdateComment`) and consumed via `OPTIONS.Comment`; the e2e supplies it directly so T1 is independently testable.
