# Dirac Eigenmode Leadfield Process (Phase 1 — save the artifact) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A process `process_dirac_eigenmode_leadfield` that composes the base unconstrained head model into the Dirac eigenbasis (via the existing `bst_dirac_eigenmode_leadfield`) and saves it as a `headmodel_dirac_eigenmode_*.mat` DB node — the first-class "Fourier-transformed head model" the viewer (Phase 2) will load.

**Architecture:** Mirror `process_eigenmode_leadfield` exactly: pick a base non-eigenmode **unconstrained surface** head model, get the cortex `DiracEigen` (compute via `tess_dirac_eigenmodes` if absent / different Tau-K), compose, `bst_save('headmodel_dirac_eigenmode')`, register the DB node. The math (`bst_dirac_eigenmode_leadfield`) is already built + tested; this is thin DB glue.

**Tech Stack:** MATLAB (Brainstorm process), `matlab.unittest`-free script test for the static dispatch path, live validation via the MATLAB MCP.

**Reference (read first):** `toolbox/process/functions/process_eigenmode_leadfield.m` — the exact base-model selection, save (`bst_process('GetNewFilename', StudyDir, 'headmodel_eigenmode')` + `bst_save(...,'v7')`), and DB registration (`db_template('headmodel')`, `bst_set('Study',...)`, `panel_protocols('UpdateNode'/'SelectNode')`, `db_save`) pattern this mirrors.

**Session discipline:** run via the MATLAB MCP; before each test, `rehash; clear <names>;` — **never** a bare `clear` (wipes Brainstorm `GlobalData`).

---

## File Structure

| File | Responsibility |
|---|---|
| `toolbox/process/functions/process_dirac_eigenmode_leadfield.m` | **Create.** Compose + save the Dirac eigenmode head model. |
| `dev/tests/test_process_dirac_eigenmode_leadfield.m` | **Create.** Static dispatch sanity (GetDescription/FormatComment). |

---

## Task 1: `process_dirac_eigenmode_leadfield` + sanity test

**Files:**
- Create: `toolbox/process/functions/process_dirac_eigenmode_leadfield.m`
- Test: `dev/tests/test_process_dirac_eigenmode_leadfield.m`

- [ ] **Step 1: Write the failing sanity test**

Create `dev/tests/test_process_dirac_eigenmode_leadfield.m`:

```matlab
function test_process_dirac_eigenmode_leadfield
% Static dispatch sanity (no DB): GetDescription + FormatComment.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);

sProcess = process_dirac_eigenmode_leadfield('GetDescription');
assert(strcmp(sProcess.Comment, 'Compute Dirac eigenmode leadfield'), 'Comment');
assert(isfield(sProcess.options,'nmodes') && isfield(sProcess.options,'tau'), 'options nmodes+tau');
assert(strcmpi(sProcess.SubGroup, 'Sources'), 'SubGroup');

sProcess.options.nmodes.Value = {0, '', 0};
c0 = process_dirac_eigenmode_leadfield('FormatComment', sProcess);
assert(contains(lower(c0), 'default'), 'FormatComment default');

sProcess.options.nmodes.Value = {50, '', 0};
c1 = process_dirac_eigenmode_leadfield('FormatComment', sProcess);
assert(contains(c1, '50'), 'FormatComment count');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run to verify it fails**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_process_dirac_eigenmode_leadfield process_dirac_eigenmode_leadfield; disp('rehashed');
test_process_dirac_eigenmode_leadfield
```
Expected: `Undefined function 'process_dirac_eigenmode_leadfield'`.

- [ ] **Step 3: Write `toolbox/process/functions/process_dirac_eigenmode_leadfield.m`**

```matlab
function varargout = process_dirac_eigenmode_leadfield( varargin )
% PROCESS_DIRAC_EIGENMODE_LEADFIELD: Compose the unconstrained leadfield into the Dirac eigenbasis.
%
% Produces a headmodel_dirac_eigenmode_*.mat node whose Gain [nCh x 2K] expresses
% the UNCONSTRAINED leadfield in the Dirac (curvature-aware vector) eigenbasis
% (Gain = Psi' * B * Phi_D). Requires a surface unconstrained head model; computes
% the Dirac eigenbasis on the cortex (tess_dirac_eigenmodes) if absent.
%
% Authors: Diellor Basha, 2026
eval(macro_method);
end

function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Compute Dirac eigenmode leadfield';
    sProcess.Category    = 'Custom';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 339;
    sProcess.Description = '';
    sProcess.InputTypes  = {'data', 'raw'};
    sProcess.OutputTypes = {'data', 'raw'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.isSeparator = 0;
    sProcess.options.nmodes.Comment = 'Number of Dirac eigenmodes per hemisphere (0 = default 400): ';
    sProcess.options.nmodes.Type    = 'value';
    sProcess.options.nmodes.Value   = {0, '', 0};
    sProcess.options.tau.Comment    = 'Curvature weighting tau [0-1]: ';
    sProcess.options.tau.Type       = 'value';
    sProcess.options.tau.Value      = {0.5, '', 3};
    sProcess.options.label_info.Comment = ['<FONT color="#777777">Expresses the UNCONSTRAINED leadfield in the ' ...
        'Dirac (curvature-aware vector) eigenbasis.<BR>Requires a surface head model; computes the Dirac ' ...
        'eigenbasis on the cortex if absent.</FONT>'];
    sProcess.options.label_info.Type = 'label';
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    n = sProcess.options.nmodes.Value{1};
    if n > 0; Comment = sprintf('Dirac eigenmode leadfield (%d modes)', n);
    else;     Comment = 'Dirac eigenmode leadfield (default modes)'; end
end

function OutputFiles = Run(sProcess, sInputs) %#ok<DEFNU>
    OutputFiles = {};
    nModes = sProcess.options.nmodes.Value{1};
    Tau    = sProcess.options.tau.Value{1};
    Kbasis = nModes; if Kbasis <= 0, Kbasis = 400; end

    [sStudy, iStudy] = bst_get('Study', sInputs(1).iStudy);
    if isempty(sStudy.HeadModel)
        bst_report('Error', sProcess, sInputs, 'No head model available for this study.'); return;
    end

    % Choose a BASE (non-eigenmode) UNCONSTRAINED surface head model.
    iBase = []; HeadModel = [];
    candidates = [sStudy.iHeadModel, setdiff(1:numel(sStudy.HeadModel), sStudy.iHeadModel)];
    for ic = candidates
        if ic < 1 || ic > numel(sStudy.HeadModel); continue; end
        try, hmC = in_bst_headmodel(sStudy.HeadModel(ic).FileName, 0); catch, continue; end
        isEig  = isfield(hmC,'isEigenmode')      && ~isempty(hmC.isEigenmode)      && hmC.isEigenmode;
        isDEig = isfield(hmC,'isDiracEigenmode') && ~isempty(hmC.isDiracEigenmode) && hmC.isDiracEigenmode;
        if ~isEig && ~isDEig && strcmpi(hmC.HeadModelType, 'surface')
            iBase = ic; HeadModel = hmC; break;
        end
    end
    if isempty(iBase)
        bst_report('Error', sProcess, sInputs, 'No base surface head model found.'); return;
    end
    HeadModelFile = sStudy.HeadModel(iBase).FileName;

    if mod(size(HeadModel.Gain,2), 3) ~= 0
        bst_report('Error', sProcess, sInputs, ...
            'Base head model must be unconstrained (Gain [nCh x 3*nVert]).'); return;
    end

    % Dirac eigenbasis on the surface (compute if absent / different Tau or K)
    T = in_tess_bst(HeadModel.SurfaceFile, 0);
    needCompute = ~isfield(T,'DiracEigen') || isempty(T.DiracEigen) ...
        || ~isequal([T.DiracEigen.Tau],[Tau Tau]) || any([T.DiracEigen.nModes] ~= Kbasis);
    if needCompute
        bst_progress('text', 'Computing Dirac eigenbasis...');
        DiracEigen = tess_dirac_eigenmodes(HeadModel.SurfaceFile, 'Tau', Tau, 'K', Kbasis);
    else
        DiracEigen = T.DiracEigen;
    end

    nVertHM = size(HeadModel.Gain, 2) / 3;
    nVertDE = sum(arrayfun(@(d) numel(d.GlobalVertices), DiracEigen));
    if nVertHM ~= nVertDE
        bst_report('Error', sProcess, sInputs, sprintf( ...
            'Head model has %d vertices but Dirac eigenbasis covers %d.', nVertHM, nVertDE)); return;
    end

    % Compose (all Kbasis modes) + history
    CompHM = bst_dirac_eigenmode_leadfield(HeadModel, DiracEigen);
    CompHM = bst_history('add', CompHM, 'dirac_eigenmode_leadfield', ...
        sprintf('Composed Dirac eigenmode leadfield: %d modes (tau=%.3g) from %s', ...
        CompHM.nModes, Tau, HeadModelFile));

    % Save + register (mirrors process_eigenmode_leadfield)
    StudyDir = bst_fileparts(file_fullpath(sStudy.FileName));
    OutputFile = bst_process('GetNewFilename', StudyDir, 'headmodel_dirac_eigenmode');
    bst_save(OutputFile, CompHM, 'v7');

    sHeadModel = db_template('headmodel');
    sHeadModel.FileName      = file_short(OutputFile);
    sHeadModel.Comment       = CompHM.Comment;
    sHeadModel.HeadModelType = CompHM.HeadModelType;
    iHM = length(sStudy.HeadModel) + 1;
    sStudy.HeadModel(iHM) = sHeadModel;
    sStudy.iHeadModel     = iHM;
    bst_set('Study', iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', iStudy);
    panel_protocols('SelectNode', [], file_short(OutputFile));
    db_save();

    % DB side-effect; pass inputs through unchanged.
    OutputFiles = {sInputs.FileName};
end
```

- [ ] **Step 4: Run the sanity test to verify it passes**

```matlab
cd('/Users/diellorbasha/workspace/research/code/brainstorm3');
rehash; clear test_process_dirac_eigenmode_leadfield process_dirac_eigenmode_leadfield; disp('rehashed');
test_process_dirac_eigenmode_leadfield
```
Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/process/functions/process_dirac_eigenmode_leadfield.m dev/tests/test_process_dirac_eigenmode_leadfield.m
git commit -m "feat(process-dirac-eigenmode-leadfield): save unconstrained leadfield in Dirac eigenbasis as headmodel node"
```

---

## Task 2: Live validation on a real study

**Files:** none (validation only; the process's real exercise is producing a loadable DB node).

- [ ] **Step 1: Find a study with a base unconstrained surface head model**

In the MATLAB MCP session (Brainstorm running):
```matlab
if ~brainstorm('status'); brainstorm nogui; end
ProtocolStudies = bst_get('ProtocolStudies');
% Scan studies for a non-eigenmode surface head model with an unconstrained Gain.
found = '';
for is = 1:numel(ProtocolStudies.Study)
    st = ProtocolStudies.Study(is);
    for ih = 1:numel(st.HeadModel)
        try, hm = in_bst_headmodel(st.HeadModel(ih).FileName, 0, 'HeadModelType','Gain','isEigenmode'); catch, continue; end
        isEig = isfield(hm,'isEigenmode') && ~isempty(hm.isEigenmode) && hm.isEigenmode;
        if ~isEig && strcmpi(hm.HeadModelType,'surface') && mod(size(hm.Gain,2),3)==0
            found = sprintf('study %d (%s), headmodel %d', is, st.Name, ih); fprintf('FOUND: %s\n', found); break;
        end
    end
    if ~isempty(found); break; end
end
if isempty(found); fprintf('NO suitable base head model found in this protocol.\n'); end
```
If none is found, document that live validation requires a protocol with a surface head model and STOP here (report the process is implemented + sanity-tested; full DB validation deferred to a protocol that has one — Phase 2's viewer will also exercise it).

- [ ] **Step 2: Run the process on a data file in that study and confirm the node**

(Only if Step 1 found a study; use a `data`/`raw` input file from it.)
```matlab
% sFiles = a data file in the found study (bst_process('CallProcess', ...) or build sInputs).
sProcess = process_dirac_eigenmode_leadfield('GetDescription');
sProcess.options.nmodes.Value = {40, '', 0};   % small K for a fast validation run
sProcess.options.tau.Value    = {0.5, '', 3};
% Build sInputs for one data file in the found study, then:
process_dirac_eigenmode_leadfield('Run', sProcess, sInputs);

% Verify a headmodel_dirac_eigenmode node was created and is well-formed:
st = bst_get('Study', sInputs(1).iStudy);
hmFile = st.HeadModel(st.iHeadModel).FileName;
HM = in_bst_headmodel(hmFile, 0);
assert(isfield(HM,'isDiracEigenmode') && HM.isDiracEigenmode==1, 'isDiracEigenmode set');
assert(size(HM.Gain,2) == HM.nModes, 'Gain [nCh x 2K] matches nModes');
assert(~isempty(strfind(hmFile,'headmodel_dirac_eigenmode')), 'node filename');
fprintf('LIVE VALIDATION OK: %s, Gain %s\n', hmFile, mat2str(size(HM.Gain)));
```
Expected: `LIVE VALIDATION OK` with `Gain [nCh x 2K]`. Restore the protocol afterward if the validation should not persist the node (delete the test node + `db_save`), or keep it to drive Phase 2.

- [ ] **Step 3: Commit any doc notes (no source change expected)**

If Step 1 found nothing, add a one-line note to the plan's status; otherwise no commit needed in this task.

---

## Final verification

- [ ] **Sanity test green** (`test_process_dirac_eigenmode_leadfield` → `ALL TESTS PASSED`).
- [ ] **If a suitable study exists,** live validation produced a `headmodel_dirac_eigenmode` node with `isDiracEigenmode` and `Gain [nCh×2K]`.
- [ ] **Then complete via superpowers:finishing-a-development-branch**, and proceed to Plan 2 (Phase 2 — viewer extension).

---

## Notes for the implementer

- This is **thin DB glue** mirroring `process_eigenmode_leadfield`; the math is the already-tested `bst_dirac_eigenmode_leadfield`. Do not re-implement composition.
- The base head model must be **unconstrained** (`Gain [nCh×3·nVert]`); `in_bst_headmodel(file, 0)` returns it unconstrained. Guarded.
- `nmodes` is the per-hemisphere Dirac basis size `K` (0 → default 400); `bst_dirac_eigenmode_leadfield` then composes all `K` (so the saved `nModes = 2K`).
- Computing the Dirac eigenbasis fresh triggers the `eigs` solve (minutes at K=400) — that's expected for "Compute Dirac eigenmode leadfield"; `tess_dirac_eigenmodes` caches it on the surface.
- Process index 339 is adjacent to the scalar `process_eigenmode_leadfield` (338); confirm 339 is free, bump if taken.
