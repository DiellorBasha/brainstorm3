# Upstream Database Schema Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** The fork stops changing Brainstorm's database format: it runs on an upstream database (DbVersion 5.03), and a script (written here, run only with Diellor's OK) takes the real database back from the fork's 5.07 to 5.03.

**Architecture:** The fork keeps the in-memory shape it uses today: each surface entry carries `Manifold`, `Operator`, `Eigen` and `Wavelet` lists, read by about 100 call sites (bst_get, tree, tess_*, tests). Those lists become memory-only. A new module, `db_derived_nodes`, rebuilds them from the files in the subject folder whenever subjects are loaded or parsed (`Attach`). Each derived file already names its parent (`ParentSurface`, or `ParentEigen` for wavelets). The module also removes them from every copy `db_save` writes (`Strip`). The fork's version bumps 5.04–5.07 are deleted, so `CurrentDbVersion` returns to upstream's 5.03.

**Tech Stack:** MATLAB R2023b, Brainstorm fork (`development`), tests are `functiontests` in `dev/tests/`, run with `matlab -batch "results = runtests('dev/tests/<file>.m'); disp(table(results))"`.

**Spec:** The decision and its reasons are in the conversation of 2026-10-26 and in this plan's header (no separate spec). Diellor: "We need it to be only compatible with the upstream and the most recent up to date Brainstorm"; "yes, do 1 before 2". Field map of the fork's schema: section *Background* below.

## Global Constraints

- Work on branch `refactor/upstream-db-schema` (from `development` 8ed040b0) in the worktree `~/workspace/research/code/brainstorm3@refactor-upstream-db-schema`. Never switch branches in the main checkout `~/workspace/research/code/brainstorm3`.
- Upstream format means: `CurrentDbVersion = 5.03`, and no `Manifold`/`Operator`/`Eigen`/`Wavelet`/`Filterbank` field in any `Surface` struct written to `protocol.mat` or `brainstorm.mat`.
- Do not change the in-memory shape (`db_template('Surface')` keeps its four lists with their current entry fields). Readers stay untouched.
- Do not change derived files on disk, their names, or their content.
- Never delete data. Tests create and delete only their own protocol (`TestUpstreamSchema`) and temp folders (`tempname`).
- **The downgrade script (Task 6) is written and dry-run tested on a copy only. Running it on the real database (`/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db`) needs Diellor's explicit OK.**
- MATLAB: one agent MATLAB batch run at a time on the hub (`operations/scripts/mem-guard.sh`). Never `clear` in a live Brainstorm session.
- Self-describing variable names in new code (`derivedNodes`, not `D`). Stage files by explicit path, never `git add -A`.
- Commits on this branch end with the session's attribution lines (it is not a clean deliverable branch).

## Background (field map, from `development`)

- `toolbox/core/bst_startup.m:208` `CurrentDbVersion = 5.07`. `toolbox/db/db_update.m` blocks `< 5.04` (378–403), `< 5.05` (405–430), `< 5.06` (432–457) widen `Surface` with `db_template('Surface')` through the local `NormalizeSurfaceArray` (613). Block `< 5.07` (468–499) fills spec fields through `BackfillDerivedSpecs` (650) and the helpers `i_backfill`, `i_fill_spec`, `i_loadmeta`, `i_provtau` (676–755).
- Entry fields (`db_template.m:36-39`):
  - `Manifold{FileName,Comment,Gauge}`
  - `Eigen{FileName,Comment,Variant,nModes,Tau,OperatorFile}`
  - `Operator{FileName,Comment,Variant,Tau}`
  - `Wavelet{FileName,Comment,ParentEigen}`
- Derived files sit in the subject folder as `manifold_*.mat`, `operator_*.mat`, `eigen_*.mat`, `wavelet_*.mat`. `file_gettype` returns `manifold`/`operator`/`eigen`/`wavelet`. The files carry these variables:
  - Manifold, operator and eigen files: `Comment` and `ParentSurface` (`file_short` form).
  - Eigen files: `Variant`, `nModes` (or `K` in old files), `Provenance.Tau`, `OperatorFile`.
  - Operator files: `Variant`, `Provenance.Tau`.
  - Manifold files: `Provenance.Gauge`.
  - Wavelet files: `Comment`, `ParentEigen` (no `ParentSurface`).
- Load paths: `bst_startup.m:242` (from brainstorm.mat), `db_load_protocol.m:100` (from protocol.mat), `toolbox/db/private/db_parse_subject.m` (reload, fix, import). Upstream's parser ignores unknown file types, so it skips derived files. Today a reload therefore drops every derived node; this plan fixes that as a side effect.
- Save path: `db_save.m:50` (`BstMat.ProtocolsListSubjects`) and `db_save.m:87` (`ProtocolMat.ProtocolSubjects`).
- Every delete path (tree callbacks 4221–4360, `db_delete_surface_node`, `db_add_*` supersede) deletes the file as well as the entry, so a disk rebuild cannot resurrect a deleted node.
- Local branches still at 5.07: `development`, `fix/preventad-import-inputstruct`, `feat/preventad-multimodal`, `fix/fibers-acpc-registration`, `feature/import-fibers-pipeline`, `feature/source-flow-processes`, `refactor/pet-metadata`. The production copy `library/software/brainstorm3` is at 5.03.

## Review Focus

1. **macOS AppleDouble files** (`._eigen_*.mat`, present on the external drive SpikeData-2) must be ignored, not loaded or listed. Test: Task 1 `testIgnoresAppleDouble`.
2. **Orphan derived files**: the subject or surface was renamed, or the parent eigen was deleted, so `ParentSurface`/`ParentEigen` matches nothing. Expected: skipped with one `BST> Derived file not attached` line naming the file, never an error, never silently lost from the dry-run report. Test: Task 1 `testOrphanIsReportedNotAttached`, Task 6 dry-run lists them.
3. **Huge derived files** (operator and eigen matrices, hundreds of MB). `Attach` must read only the small variables. Test: Task 1 `testScanReadsOnlySmallFields`.
4. **Repeated attach** (startup, then protocol load, then reload) must not duplicate entries. Test: Task 1 `testAttachIsIdempotent`.
5. **Protocols not loaded, or read-only**: `db_save` skips them, and `brainstorm.mat` holds empty templates for them. `Strip` must accept empty and missing `DefaultSubject` and `Subject`. Test: Task 1 `testStripHandlesEmptyProtocol`.

## File Structure

- Create `toolbox/db/db_derived_nodes.m`: owns the memory-only derived lists (`Scan`, `AttachSubject`, `Attach`, `Strip`, `NormalizeSurfaceArray`).
- Modify `toolbox/db/db_save.m`: strip before writing.
- Modify `toolbox/core/bst_startup.m`, `toolbox/db/db_load_protocol.m`, `toolbox/db/private/db_parse_subject.m`: attach after loading.
- Modify `toolbox/db/db_update.m`: delete the 5.04–5.07 blocks and their helpers (`NormalizeSurfaceArray` moves to the new module).
- Modify `toolbox/core/bst_startup.m:208`: `CurrentDbVersion = 5.03`.
- Create `dev/tests/test_db_derived_nodes.m` (unit, no database) and `dev/tests/test_db_upstream_schema.m` (integration, own protocol).
- Create `dev/db_downgrade_to_upstream.m`: backup, dry-run report, downgrade.

---

### Task 1: Memory-only derived lists (`db_derived_nodes`)

**Files:**
- Create: `toolbox/db/db_derived_nodes.m`
- Test: `dev/tests/test_db_derived_nodes.m`

**Interfaces:**
- Produces:
  - `sSubject = db_derived_nodes('AttachSubject', sSubject, SubjectsRoot)`. Rebuilds `sSubject.Surface(i).Manifold/Operator/Eigen/Wavelet` from the files in the subject's folder. The rebuild replaces the lists and never appends.
  - `sProtocolSubjects = db_derived_nodes('Attach', sProtocolSubjects, SubjectsRoot)` does the same for `DefaultSubject` and every `Subject` of one protocol.
  - `sProtocolSubjects = db_derived_nodes('Strip', sProtocolSubjects)` removes the derived fields from every surface of every protocol in the array.
  - `derivedNodes = db_derived_nodes('Scan', SubjectsRoot, SubjectDir)` returns a struct array `{Type, FileName, ParentSurface, ParentEigen, Entry}`.
  - `[sSurf, isFixed] = db_derived_nodes('NormalizeSurfaceArray', sSurf)` widens a Surface array to `db_template('Surface')`.
  - `fieldNames = db_derived_nodes('DerivedFields')` returns `{'Manifold','Operator','Eigen','Wavelet','Filterbank'}`.

- [ ] **Step 1: Write the failing tests**

`dev/tests/test_db_derived_nodes.m`:

```matlab
function tests = test_db_derived_nodes
% Tests for db_derived_nodes: derived surface lists rebuilt from disk, stripped before save.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    if ~brainstorm('status'); brainstorm nogui; end
    testCase.TestData.SubjectsRoot = tempname;
    mkdir(fullfile(testCase.TestData.SubjectsRoot, 'Subject01'));
end

function teardownOnce(testCase)
    rmdir(testCase.TestData.SubjectsRoot, 's');
end

function setup(testCase)
    % Fresh subject folder per test
    subjectDir = fullfile(testCase.TestData.SubjectsRoot, 'Subject01');
    delete(fullfile(subjectDir, '*.mat'));
    delete(fullfile(subjectDir, '._*'));
end

% ---- helpers ----
function sSubject = local_subject()
    sSubject = db_template('Subject');
    sSubject.FileName = 'Subject01/brainstormsubject.mat';
    sSurface = db_template('Surface');
    sSurface.FileName = 'Subject01/tess_cortex_pial_low.mat';
    sSurface.SurfaceType = 'Cortex';
    % Upstream-shaped surface: none of the derived fields
    sSubject.Surface = rmfield(sSurface, intersect(db_derived_nodes('DerivedFields'), fieldnames(sSurface)));
end

function local_save(testCase, fileName, contents)
    save(fullfile(testCase.TestData.SubjectsRoot, 'Subject01', fileName), '-struct', 'contents');
end

function local_save_all_kinds(testCase)
    local_save(testCase, 'manifold_261026_1200.mat', struct('Comment', 'Manifold', ...
        'ParentSurface', 'Subject01/tess_cortex_pial_low.mat', 'Provenance', struct('Gauge', 'Canonical')));
    local_save(testCase, 'operator_261026_1201.mat', struct('Comment', 'Dirac', 'Variant', 'Dirac', ...
        'ParentSurface', 'Subject01/tess_cortex_pial_low.mat', 'Provenance', struct('Tau', 0.5)));
    local_save(testCase, 'eigen_261026_1202.mat', struct('Comment', 'Dirac 1200', 'Variant', 'Dirac', 'nModes', 1200, ...
        'ParentSurface', 'Subject01/tess_cortex_pial_low.mat', 'OperatorFile', 'Subject01/operator_261026_1201.mat', ...
        'Provenance', struct('Tau', 0.5)));
    local_save(testCase, 'wavelet_261026_1203.mat', struct('Comment', 'Wavelet', ...
        'ParentEigen', 'Subject01/eigen_261026_1202.mat'));
end

% ---- tests ----
function testAttachAllKinds(testCase)
    local_save_all_kinds(testCase);
    sSubject = db_derived_nodes('AttachSubject', local_subject(), testCase.TestData.SubjectsRoot);
    sSurface = sSubject.Surface(1);
    verifyEqual(testCase, {sSurface.Manifold.FileName}, {'Subject01/manifold_261026_1200.mat'});
    verifyEqual(testCase, sSurface.Manifold.Gauge, 'Canonical');
    verifyEqual(testCase, sSurface.Operator.Variant, 'Dirac');
    verifyEqual(testCase, sSurface.Operator.Tau, 0.5);
    verifyEqual(testCase, sSurface.Eigen.nModes, 1200);
    verifyEqual(testCase, sSurface.Eigen.OperatorFile, 'Subject01/operator_261026_1201.mat');
    verifyEqual(testCase, sSurface.Wavelet.ParentEigen, 'Subject01/eigen_261026_1202.mat');
    % Entry field sets match the template exactly (homogeneous with db_add_* appends)
    templateSurface = db_template('Surface');
    verifyEqual(testCase, fieldnames(sSurface.Eigen), fieldnames(templateSurface.Eigen));
    verifyEqual(testCase, fieldnames(sSurface.Wavelet), fieldnames(templateSurface.Wavelet));
end

function testEigenFallsBackToK(testCase)
    local_save(testCase, 'eigen_261026_1300.mat', struct('Comment', 'LBO', 'Variant', 'LBO', 'K', 400, ...
        'ParentSurface', 'Subject01/tess_cortex_pial_low.mat'));
    sSubject = db_derived_nodes('AttachSubject', local_subject(), testCase.TestData.SubjectsRoot);
    verifyEqual(testCase, sSubject.Surface(1).Eigen.nModes, 400);
    verifyEmpty(testCase, sSubject.Surface(1).Eigen.Tau);
end

function testIgnoresAppleDouble(testCase)
    local_save_all_kinds(testCase);
    fid = fopen(fullfile(testCase.TestData.SubjectsRoot, 'Subject01', '._eigen_261026_1202.mat'), 'w');
    fwrite(fid, uint8([0 5 22 7 0 2 0 0])); fclose(fid);
    sSubject = db_derived_nodes('AttachSubject', local_subject(), testCase.TestData.SubjectsRoot);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
end

function testOrphanIsReportedNotAttached(testCase)
    local_save(testCase, 'eigen_261026_1400.mat', struct('Comment', 'Old', 'Variant', 'LBO', 'nModes', 400, ...
        'ParentSurface', 'Subject01/tess_cortex_renamed.mat'));
    local_save(testCase, 'wavelet_261026_1401.mat', struct('Comment', 'W', 'ParentEigen', 'Subject01/eigen_gone.mat'));
    [output, sSubject] = evalc('db_derived_nodes(''AttachSubject'', local_subject(), testCase.TestData.SubjectsRoot)');
    verifyEmpty(testCase, sSubject.Surface(1).Eigen);
    verifyEmpty(testCase, sSubject.Surface(1).Wavelet);
    verifySubstring(testCase, output, 'BST> Derived file not attached: Subject01/eigen_261026_1400.mat');
    verifySubstring(testCase, output, 'BST> Derived file not attached: Subject01/wavelet_261026_1401.mat');
end

function testScanReadsOnlySmallFields(testCase)
    local_save(testCase, 'operator_261026_1500.mat', struct('Comment', 'LBO', 'Variant', 'LBO', ...
        'ParentSurface', 'Subject01/tess_cortex_pial_low.mat', 'Operator', {{rand(1500)}}));
    derivedNodes = db_derived_nodes('Scan', testCase.TestData.SubjectsRoot, 'Subject01');
    verifyNumElements(testCase, derivedNodes, 1);
    verifyFalse(testCase, isfield(derivedNodes.Entry, 'Operator'));
    verifyEqual(testCase, derivedNodes.Entry.Variant, 'LBO');
end

function testAttachIsIdempotent(testCase)
    local_save_all_kinds(testCase);
    sSubject = db_derived_nodes('AttachSubject', local_subject(), testCase.TestData.SubjectsRoot);
    sSubject = db_derived_nodes('AttachSubject', sSubject, testCase.TestData.SubjectsRoot);
    verifyNumElements(testCase, sSubject.Surface(1).Manifold, 1);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
    verifyNumElements(testCase, sSubject.Surface(1).Wavelet, 1);
end

function testStripRemovesDerivedFields(testCase)
    local_save_all_kinds(testCase);
    sProtocolSubjects = db_template('ProtocolSubjects');
    sProtocolSubjects.Subject = db_derived_nodes('AttachSubject', local_subject(), testCase.TestData.SubjectsRoot);
    sProtocolSubjects = db_derived_nodes('Strip', sProtocolSubjects);
    remaining = intersect(db_derived_nodes('DerivedFields'), fieldnames(sProtocolSubjects.Subject(1).Surface));
    verifyEmpty(testCase, remaining);
    verifyEqual(testCase, sProtocolSubjects.Subject(1).Surface(1).FileName, 'Subject01/tess_cortex_pial_low.mat');
end

function testStripHandlesEmptyProtocol(testCase)
    sProtocolSubjects = [db_template('ProtocolSubjects'), db_template('ProtocolSubjects')];
    sProtocolSubjects(2).DefaultSubject = [];
    sStripped = db_derived_nodes('Strip', sProtocolSubjects);
    verifySize(testCase, sStripped, [1 2]);
end
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_derived_nodes.m'); disp(table(results))"`
Expected: every test FAILS with `Undefined function 'db_derived_nodes'`.

- [ ] **Step 3: Write the module**

`toolbox/db/db_derived_nodes.m`:

```matlab
function varargout = db_derived_nodes(varargin)
% DB_DERIVED_NODES: Derived surface nodes (manifold, operator, eigen, wavelet), kept out of the saved database.
%
% USAGE:  sSubject          = db_derived_nodes('AttachSubject', sSubject, SubjectsRoot)
%         sProtocolSubjects = db_derived_nodes('Attach', sProtocolSubjects, SubjectsRoot)
%         sProtocolSubjects = db_derived_nodes('Strip', sProtocolSubjects)
%         derivedNodes      = db_derived_nodes('Scan', SubjectsRoot, SubjectDir)
%         [sSurf, isFixed]  = db_derived_nodes('NormalizeSurfaceArray', sSurf)
%         fieldNames        = db_derived_nodes('DerivedFields')
%
% Derived files are listed as children of their surface (sSubject.Surface(i).Manifold, .Operator,
% .Eigen, .Wavelet). The lists live in memory only: 'Attach' rebuilds them from the files in the
% subject folder (each file names its parent: ParentSurface, or ParentEigen for wavelets), and
% db_save calls 'Strip' before writing protocol.mat and brainstorm.mat, so the saved database keeps
% upstream Brainstorm's format.

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
% =============================================================================@
%
% Authors: Diellor Basha, 2026

eval(macro_method);
end


%% ===== DERIVED FIELDS =====
% Surface fields that never go to disk ('Filterbank' only exists in databases migrated to 5.05)
function fieldNames = DerivedFields()
    fieldNames = {'Manifold', 'Operator', 'Eigen', 'Wavelet', 'Filterbank'};
end


%% ===== STRIP =====
function sProtocolSubjects = Strip(sProtocolSubjects)
    for iProtocol = 1:numel(sProtocolSubjects)
        if isfield(sProtocolSubjects(iProtocol), 'DefaultSubject') && ~isempty(sProtocolSubjects(iProtocol).DefaultSubject)
            sProtocolSubjects(iProtocol).DefaultSubject.Surface = StripSurface(sProtocolSubjects(iProtocol).DefaultSubject.Surface);
        end
        if isfield(sProtocolSubjects(iProtocol), 'Subject')
            for iSubject = 1:numel(sProtocolSubjects(iProtocol).Subject)
                sProtocolSubjects(iProtocol).Subject(iSubject).Surface = StripSurface(sProtocolSubjects(iProtocol).Subject(iSubject).Surface);
            end
        end
    end
end

function sSurf = StripSurface(sSurf)
    if ~isstruct(sSurf)
        return;
    end
    derivedPresent = intersect(DerivedFields(), fieldnames(sSurf), 'stable');
    if ~isempty(derivedPresent)
        sSurf = rmfield(sSurf, derivedPresent);
    end
end


%% ===== ATTACH (one protocol) =====
function sProtocolSubjects = Attach(sProtocolSubjects, SubjectsRoot)
    if isempty(sProtocolSubjects) || isempty(SubjectsRoot)
        return;
    end
    if isfield(sProtocolSubjects, 'DefaultSubject') && ~isempty(sProtocolSubjects.DefaultSubject)
        sProtocolSubjects.DefaultSubject = AttachSubject(sProtocolSubjects.DefaultSubject, SubjectsRoot);
    end
    if isfield(sProtocolSubjects, 'Subject') && ~isempty(sProtocolSubjects.Subject)
        subjectCells = arrayfun(@(sSubject) AttachSubject(sSubject, SubjectsRoot), sProtocolSubjects.Subject, 'UniformOutput', 0);
        sProtocolSubjects.Subject = reshape([subjectCells{:}], size(sProtocolSubjects.Subject));
    end
end


%% ===== ATTACH (one subject) =====
function sSubject = AttachSubject(sSubject, SubjectsRoot)
    if isempty(sSubject) || ~isfield(sSubject, 'Surface') || isempty(sSubject.Surface) || isempty(sSubject.FileName)
        return;
    end
    templateSurface = db_template('Surface');
    sSurf = NormalizeSurfaceArray(sSubject.Surface, templateSurface);
    % Rebuild, never append: start every list empty
    for iSurf = 1:numel(sSurf)
        sSurf(iSurf).Manifold = templateSurface.Manifold;
        sSurf(iSurf).Operator = templateSurface.Operator;
        sSurf(iSurf).Eigen    = templateSurface.Eigen;
        sSurf(iSurf).Wavelet  = templateSurface.Wavelet;
    end
    derivedNodes = Scan(SubjectsRoot, bst_fileparts(sSubject.FileName));
    surfaceFiles = {sSurf.FileName};
    % Manifold, operator and eigen files hang from their surface
    for iNode = 1:numel(derivedNodes)
        node = derivedNodes(iNode);
        if strcmp(node.Type, 'wavelet')
            continue;
        end
        iSurf = find(file_compare(node.ParentSurface, surfaceFiles), 1);
        if isempty(iSurf)
            disp(['BST> Derived file not attached: ' node.FileName]);
            continue;
        end
        switch node.Type
            case 'manifold',  sSurf(iSurf).Manifold(end+1) = node.Entry;
            case 'operator',  sSurf(iSurf).Operator(end+1) = node.Entry;
            case 'eigen',     sSurf(iSurf).Eigen(end+1)    = node.Entry;
        end
    end
    % Wavelet files hang from the surface that holds their parent eigen file
    for iNode = 1:numel(derivedNodes)
        node = derivedNodes(iNode);
        if ~strcmp(node.Type, 'wavelet')
            continue;
        end
        iSurfParent = [];
        for iSurf = 1:numel(sSurf)
            if ~isempty(sSurf(iSurf).Eigen) && any(file_compare(node.ParentEigen, {sSurf(iSurf).Eigen.FileName}))
                iSurfParent = iSurf;
                break;
            end
        end
        if isempty(iSurfParent)
            disp(['BST> Derived file not attached: ' node.FileName]);
            continue;
        end
        sSurf(iSurfParent).Wavelet(end+1) = node.Entry;
    end
    sSubject.Surface = sSurf;
end


%% ===== SCAN =====
% List the derived files of one subject folder with their list entries (small variables only)
function derivedNodes = Scan(SubjectsRoot, SubjectDir)
    derivedNodes = repmat(struct('Type', '', 'FileName', '', 'ParentSurface', '', 'ParentEigen', '', 'Entry', []), 1, 0);
    templateSurface = db_template('Surface');
    listing = dir(bst_fullfile(SubjectsRoot, SubjectDir, '*.mat'));
    % Sorted by name = by creation time (names end in _YYMMDD_HHMM)
    [~, iSorted] = sort({listing.name});
    listing = listing(iSorted);
    for iFile = 1:numel(listing)
        fileName = listing(iFile).name;
        % macOS AppleDouble companions on external drives
        if strncmp(fileName, '._', 2)
            continue;
        end
        fileType = file_gettype(fileName);
        if ~ismember(fileType, {'manifold', 'operator', 'eigen', 'wavelet'})
            continue;
        end
        relativeFile = bst_fullfile(SubjectDir, fileName);
        fileMeta = LoadMeta(bst_fullfile(SubjectsRoot, relativeFile), ...
            {'Comment', 'ParentSurface', 'ParentEigen', 'Variant', 'nModes', 'K', 'OperatorFile', 'Provenance'});
        node.Type          = fileType;
        node.FileName      = file_short(relativeFile);
        node.ParentSurface = MetaValue(fileMeta, 'ParentSurface', '');
        node.ParentEigen   = MetaValue(fileMeta, 'ParentEigen', '');
        provenance         = MetaValue(fileMeta, 'Provenance', struct());
        switch fileType
            case 'manifold'
                node.Entry = NewEntry(templateSurface.Manifold, 'FileName', node.FileName, ...
                    'Comment', MetaValue(fileMeta, 'Comment', 'Manifold'), 'Gauge', MetaValue(provenance, 'Gauge', ''));
            case 'operator'
                node.Entry = NewEntry(templateSurface.Operator, 'FileName', node.FileName, ...
                    'Comment', MetaValue(fileMeta, 'Comment', 'Operator'), 'Variant', MetaValue(fileMeta, 'Variant', ''), ...
                    'Tau', MetaValue(provenance, 'Tau', []));
            case 'eigen'
                nModes = MetaValue(fileMeta, 'nModes', MetaValue(fileMeta, 'K', []));
                node.Entry = NewEntry(templateSurface.Eigen, 'FileName', node.FileName, ...
                    'Comment', MetaValue(fileMeta, 'Comment', 'Eigen'), 'Variant', MetaValue(fileMeta, 'Variant', ''), ...
                    'nModes', nModes, 'Tau', MetaValue(provenance, 'Tau', []), ...
                    'OperatorFile', MetaValue(fileMeta, 'OperatorFile', ''));
            case 'wavelet'
                node.Entry = NewEntry(templateSurface.Wavelet, 'FileName', node.FileName, ...
                    'Comment', MetaValue(fileMeta, 'Comment', 'Wavelet'), 'ParentEigen', node.ParentEigen);
        end
        derivedNodes(end+1) = node; %#ok<AGROW>
    end
end

% Load only the named variables (missing ones are simply absent); unreadable file = empty struct
function fileMeta = LoadMeta(fileFull, variableNames)
    fileMeta = struct();
    try
        variablesInFile = who('-file', fileFull);
        variablesToLoad = intersect(variableNames, variablesInFile, 'stable');
        if ~isempty(variablesToLoad)
            fileMeta = load(fileFull, variablesToLoad{:});
        end
    catch
        disp(['BST> Derived file could not be read: ' fileFull]);
    end
end

function value = MetaValue(sInput, fieldName, defaultValue)
    if isstruct(sInput) && isfield(sInput, fieldName) && ~isempty(sInput.(fieldName))
        value = sInput.(fieldName);
    else
        value = defaultValue;
    end
end

% One list entry with exactly the template's fields, in template order
function entry = NewEntry(templateEntries, varargin)
    entryFields = fieldnames(templateEntries);
    entryValues = cell(numel(entryFields), 1);
    for iArg = 1:2:numel(varargin)
        entryValues{strcmp(entryFields, varargin{iArg})} = varargin{iArg + 1};
    end
    entry = cell2struct(entryValues, entryFields, 1);
end


%% ===== NORMALIZE SURFACE ARRAY =====
% Backfill any missing template fields on every element of a Surface struct array, keeping existing
% values and template field order, so that "Surface(i) = db_template('Surface')" assignments succeed.
% (Moved from db_update, where the fork's 5.04-5.06 migrations used it.)
function [sSurf, isFixed] = NormalizeSurfaceArray(sSurf, templateSurface)
    if (nargin < 2) || isempty(templateSurface)
        templateSurface = db_template('Surface');
    end
    isFixed = 0;
    if isempty(sSurf)
        return;
    end
    templateFields = fieldnames(templateSurface);
    if all(isfield(sSurf, templateFields))
        return;
    end
    isFixed = 1;
    originalSize = size(sSurf);
    % A struct array cannot gain fields one element at a time: rebuild it element by element
    surfaceCells = cell(1, numel(sSurf));
    for iSurf = 1:numel(sSurf)
        sOne = struct_copy_fields(sSurf(iSurf), templateSurface, 0);
        extraFields = setdiff(fieldnames(sOne), templateFields, 'stable');
        surfaceCells{iSurf} = orderfields(sOne, [templateFields; extraFields]);
    end
    sSurf = reshape([surfaceCells{:}], originalSize);
end
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_derived_nodes.m'); disp(table(results))"`
Expected: 9 passed, 0 failed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/db/db_derived_nodes.m dev/tests/test_db_derived_nodes.m
git commit -m "feat(db): db_derived_nodes — derived surface lists rebuilt from disk, stripped before save"
```

---

### Task 2: `db_save` writes upstream-format surfaces

**Files:**
- Modify: `toolbox/db/db_save.m:50` and `:87`
- Test: `dev/tests/test_db_upstream_schema.m` (created here; Tasks 3–4 add tests to it)

**Interfaces:**
- Consumes: `db_derived_nodes('Strip', sProtocolSubjects)` (Task 1).
- Produces: `protocol.mat` and `brainstorm.mat` without derived fields, while `GlobalData.DataBase.ProtocolSubjects` keeps them.

- [ ] **Step 1: Write the failing test**

`dev/tests/test_db_upstream_schema.m`:

```matlab
function tests = test_db_upstream_schema
% Integration tests: the fork saves an upstream-format database and rebuilds derived lists on load.
% Uses its own protocol 'TestUpstreamSchema', created and deleted here.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    if ~brainstorm('status'); brainstorm nogui; end
    gui_brainstorm('DeleteProtocol', 'TestUpstreamSchema');
    gui_brainstorm('CreateProtocol', 'TestUpstreamSchema', 0, 0);
    [~, iSubject] = db_add_subject('Subject01');
    % Smallest real surface: an icosphere (never a hand-built mesh)
    [sphereVertices, sphereFaces] = tess_sphere(162);
    sSurfaceMat = db_template('surfacemat');
    sSurfaceMat.Comment  = 'cortex_test';
    sSurfaceMat.Vertices = sphereVertices;
    sSurfaceMat.Faces    = sphereFaces;
    testCase.TestData.SurfaceFile = db_add(iSubject, sSurfaceMat);
    testCase.TestData.iSubject    = iSubject;
    % One derived file on disk, registered the way the fork does it today
    sEigenMat = struct('Comment', 'LBO test', 'Variant', 'LBO', 'nModes', 10, 'Phi', zeros(162, 10));
    db_add_eigen(iSubject, testCase.TestData.SurfaceFile, sEigenMat, 'LBO test');
end

function teardownOnce(~)
    gui_brainstorm('DeleteProtocol', 'TestUpstreamSchema');
end

function testSavedProtocolHasNoDerivedFields(testCase)
    db_save(1);
    ProtocolInfo = bst_get('ProtocolInfo');
    savedProtocol = load(bst_fullfile(ProtocolInfo.STUDIES, 'protocol.mat'), 'ProtocolSubjects');
    savedSurface  = savedProtocol.ProtocolSubjects.Subject(1).Surface;
    verifyEmpty(testCase, intersect(db_derived_nodes('DerivedFields'), fieldnames(savedSurface)));
    % Memory still lists the derived file
    sSubject = bst_get('Subject', testCase.TestData.iSubject);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
end

function testSavedSettingsHaveNoDerivedFields(testCase)
    db_save(1);
    savedSettings = load(bst_get('BrainstormDbFile'), 'ProtocolsListSubjects');
    for iProtocol = 1:numel(savedSettings.ProtocolsListSubjects)
        for iSubject = 1:numel(savedSettings.ProtocolsListSubjects(iProtocol).Subject)
            savedSurface = savedSettings.ProtocolsListSubjects(iProtocol).Subject(iSubject).Surface;
            if isstruct(savedSurface)
                verifyEmpty(testCase, intersect(db_derived_nodes('DerivedFields'), fieldnames(savedSurface)));
            end
        end
    end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_upstream_schema.m'); disp(table(results))"`
Expected: `testSavedProtocolHasNoDerivedFields` FAILS (the saved surface still has `Eigen` etc.). `testSavedSettingsHaveNoDerivedFields` may pass, because `db_save` empties loaded protocols in `brainstorm.mat`; it guards the unloaded case.

- [ ] **Step 3: Strip in `db_save`**

In `toolbox/db/db_save.m`, replace line 50:

```matlab
BstMat.ProtocolsListSubjects = GlobalData.DataBase.ProtocolSubjects;
```

with:

```matlab
% Derived surface lists (manifold/operator/eigen/wavelet) are rebuilt from disk on load: never saved
BstMat.ProtocolsListSubjects = db_derived_nodes('Strip', GlobalData.DataBase.ProtocolSubjects);
```

and replace line 87:

```matlab
    ProtocolMat.ProtocolSubjects  = GlobalData.DataBase.ProtocolSubjects(iProtocol);
```

with:

```matlab
    ProtocolMat.ProtocolSubjects  = db_derived_nodes('Strip', GlobalData.DataBase.ProtocolSubjects(iProtocol));
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_upstream_schema.m'); disp(table(results))"`
Expected: 2 passed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/db/db_save.m dev/tests/test_db_upstream_schema.m
git commit -m "fix(db): db_save writes upstream-format surfaces (derived lists stay in memory)"
```

---

### Task 3: Rebuild the derived lists on every load

**Files:**
- Modify: `toolbox/core/bst_startup.m:242`, `toolbox/db/db_load_protocol.m:100`, `toolbox/db/private/db_parse_subject.m` (end of function)
- Test: `dev/tests/test_db_upstream_schema.m` (add tests)

**Interfaces:**
- Consumes: `db_derived_nodes('Attach', sProtocolSubjects, SubjectsRoot)` and `db_derived_nodes('AttachSubject', sSubject, SubjectsRoot)` (Task 1).
- Produces: derived lists present in memory after startup, protocol load and any subject parse (reload, fix, import).

- [ ] **Step 1: Write the failing tests** (append to `dev/tests/test_db_upstream_schema.m`)

```matlab
function testReloadSubjectKeepsDerivedNodes(testCase)
    db_reload_subjects(testCase.TestData.iSubject);
    sSubject = bst_get('Subject', testCase.TestData.iSubject);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
    verifyEqual(testCase, sSubject.Surface(1).Eigen.Variant, 'LBO');
end

function testReloadProtocolKeepsDerivedNodes(testCase)
    db_save(1);
    db_reload_database('current');
    sSubject = bst_get('Subject', testCase.TestData.iSubject);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
end

function testLoadFromProtocolFileKeepsDerivedNodes(testCase)
    % Switch away and back: the protocol is unloaded, then read from protocol.mat
    db_save(1);
    iTestProtocol = bst_get('iProtocol');
    gui_brainstorm('CreateProtocol', 'TestUpstreamSchemaOther', 0, 0);
    gui_brainstorm('SetCurrentProtocol', iTestProtocol);
    sSubject = bst_get('Subject', testCase.TestData.iSubject);
    verifyNumElements(testCase, sSubject.Surface(1).Eigen, 1);
    gui_brainstorm('DeleteProtocol', 'TestUpstreamSchemaOther');
    gui_brainstorm('SetCurrentProtocol', bst_get('Protocol', 'TestUpstreamSchema'));
end
```

- [ ] **Step 2: Run them to verify they fail**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_upstream_schema.m'); disp(table(results))"`
Expected: `testReloadSubjectKeepsDerivedNodes` and `testReloadProtocolKeepsDerivedNodes` FAIL (0 Eigen entries: the reload drops them, the known gotcha). `testLoadFromProtocolFileKeepsDerivedNodes` FAILS after Task 2 (the saved file has no lists).

- [ ] **Step 3: Attach at the three load points**

At the end of `toolbox/db/private/db_parse_subject.m` (after `sSubject` is fully built, before the function returns), add:

```matlab
% Derived surface lists (manifold/operator/eigen/wavelet) are not in the database: rebuild from disk
sSubject = db_derived_nodes('AttachSubject', sSubject, subjectsDir);
```

(`subjectsDir` is the function's first argument, the protocol's SUBJECTS folder: `function sSubject = db_parse_subject( subjectsDir, subjectSubDir, sizeProgress )`. If the function has early `return`s after `sSubject` is built, attach before each of them.)

In `toolbox/db/db_load_protocol.m`, replace line 100:

```matlab
        GlobalData.DataBase.ProtocolSubjects(iProtocols(i))   = ProtocolMat.ProtocolSubjects;
```

with:

```matlab
        GlobalData.DataBase.ProtocolSubjects(iProtocols(i))   = db_derived_nodes('Attach', ProtocolMat.ProtocolSubjects, GlobalData.DataBase.ProtocolInfo(iProtocols(i)).SUBJECTS);
```

In `toolbox/core/bst_startup.m`, inside the same block, after `isDbLoaded = 1;` (line 254 on `development`, which follows the `isProtocolLoaded` assignment at 249–253), add:

```matlab
        % Derived surface lists are not in the database: rebuild them for the loaded protocols
        for iProtocolLoaded = find(GlobalData.DataBase.isProtocolLoaded(:)')
            GlobalData.DataBase.ProtocolSubjects(iProtocolLoaded) = db_derived_nodes('Attach', ...
                GlobalData.DataBase.ProtocolSubjects(iProtocolLoaded), GlobalData.DataBase.ProtocolInfo(iProtocolLoaded).SUBJECTS);
        end
```

(This runs before `db_update` at line ~300. On a 5.07 database `db_update` then refuses as upstream does. That is intended: Task 6's script downgrades first.)

- [ ] **Step 4: Run the tests to verify they pass**

Run: `cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests('dev/tests/test_db_upstream_schema.m'); disp(table(results))"`
Expected: 5 passed.

- [ ] **Step 5: Commit**

```bash
git add toolbox/db/private/db_parse_subject.m toolbox/db/db_load_protocol.m toolbox/core/bst_startup.m dev/tests/test_db_upstream_schema.m
git commit -m "fix(db): rebuild derived surface lists from disk on startup, protocol load and subject parse"
```

---

### Task 4: Back to upstream's database version (5.03)

**Files:**
- Modify: `toolbox/core/bst_startup.m:208`
- Modify: `toolbox/db/db_update.m` (delete blocks `< 5.04`, `< 5.05`, `< 5.06`, `< 5.07` and the helpers `NormalizeSurfaceArray`, `BackfillDerivedSpecs`, `i_backfill`, `i_fill_spec`, `i_loadmeta`, `i_provtau`)
- Modify: every other caller of `db_update`'s local `NormalizeSurfaceArray`, if any (`grep -rn NormalizeSurfaceArray toolbox`), to call `db_derived_nodes('NormalizeSurfaceArray', ...)`
- Test: `dev/tests/test_db_upstream_schema.m` (add test)

**Interfaces:**
- Consumes: Tasks 1–3 (lists come from disk, so no migration is needed).
- Produces: `GlobalData.DataBase.DbVersion == 5.03` after startup on an upstream database.

- [ ] **Step 1: Write the failing test** (append)

```matlab
function testDatabaseVersionIsUpstream(testCase)
    global GlobalData
    verifyEqual(testCase, GlobalData.DataBase.DbVersion, 5.03);
    db_save(1);
    ProtocolInfo = bst_get('ProtocolInfo');
    savedProtocol = load(bst_fullfile(ProtocolInfo.STUDIES, 'protocol.mat'), 'DbVersion');
    verifyEqual(testCase, savedProtocol.DbVersion, 5.03);
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: same command as Task 3.
Expected: `testDatabaseVersionIsUpstream` FAILS: actual 5.07.

- [ ] **Step 3: Change the version and delete the fork migrations**

In `toolbox/core/bst_startup.m:208`: `CurrentDbVersion = 5.07;` becomes `CurrentDbVersion = 5.03;`.

In `toolbox/db/db_update.m`, delete the four `if (CurrentDbVersion < 5.04)` … `< 5.07` blocks (lines 376–499 on `development`, including their header comments) and the local functions listed above (lines ~609–755). Keep upstream's blocks (5.02 at 331, 5.03 at 366) and everything else. Then run `git diff upstream/master -- toolbox/db/db_update.m` and check that only upstream's own content remains, or say what else differs.

- [ ] **Step 4: Run all database tests and the existing derived-node tests**

Run:
```bash
cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema && matlab -batch "addpath(pwd); results = runtests({'dev/tests/test_db_derived_nodes.m','dev/tests/test_db_upstream_schema.m','dev/tests/test_tess_manifold.m','dev/tests/test_tess_operators.m','dev/tests/test_tess_eigen.m','dev/tests/test_bst_eigen_resolve.m'}); disp(table(results(~[results.Passed])))" 2>&1 | tail -30
```
Expected: no failed rows. A test in the last four that fails *only* because it assumed a reload drops nodes may be updated to the new behaviour. Any other failure is a bug in this branch; fix it before committing.

- [ ] **Step 5: Commit**

```bash
git add toolbox/core/bst_startup.m toolbox/db/db_update.m dev/tests/test_db_upstream_schema.m
git commit -m "fix(db): back to upstream database version 5.03; drop the fork's 5.04-5.07 migrations"
```

---

### Task 5: Check the result in Brainstorm on a copy of a real protocol

**Files:** none changed (verification only; write findings to the commit message of Task 6 or a note).

- [ ] **Step 1:** Copy one real protocol that has derived nodes (e.g. `omega-tutorial-cortical-flow`) to `/Volumes/Pergammon/pergammon/tmp/upstream-db-schema/brainstorm_db/`. Use `rsync -a`, which only reads the source. Copy `brainstorm.mat` to a separate settings folder there.
- [ ] **Step 2:** Run Task 6's script on that copy (dry run, then for real), pointing it at the copy.
- [ ] **Step 3:** Start the fork (this branch) on the copy, with the settings-folder isolation used for the PET tutorial run of 2026-10-26: `java.lang.System.setProperty('user.home', <copy settings folder>)` before `brainstorm nogui`. Check that every derived node the 5.07 database listed is in the tree (`bst_get('Subject', i).Surface(j).Eigen` etc.), and that `tess_eigen` finds an existing eigen node rather than recomputing (`bst_get('EigenFileForSurface', ...)` returns it).
- [ ] **Step 4:** Start upstream Brainstorm (`pr/pet-pvc-5`, or a checkout of `upstream/master`) on the same copy. It must open without "created with a newer version", and show the subject and surfaces.
- [ ] **Step 5:** Run the GUI checks of `tree_callbacks` by hand if a display is available: right-click a derived node, then View, then Delete. After a reload, the deleted node must stay gone.

---

### Task 6: The downgrade script (written and tested on a copy; not run on the real database)

**Files:**
- Create: `dev/db_downgrade_to_upstream.m`

**Interfaces:**
- Consumes: `db_derived_nodes('Scan'|'AttachSubject'|'Strip')` (Task 1). Runs without starting Brainstorm's GUI or database (`brainstorm` on the path only, for `db_template`, `file_*`, `bst_fullfile`).
- Produces: `report = db_downgrade_to_upstream(DbDir, SettingsFile, 'DryRun', true|false)`. `report` is a struct array, one row per protocol: `Protocol`, `OldVersion`, `ListedBefore`, `RebuiltFromDisk`, `MissingAfter` (cell of files listed at 5.07 that disk does not rebuild), `Orphans` (cell of derived files with no parent), `BackupFolder`.

- [ ] **Step 1: Write the script**

```matlab
function report = db_downgrade_to_upstream(DbDir, SettingsFile, varargin)
% DB_DOWNGRADE_TO_UPSTREAM: Take a database migrated by the fork (DbVersion 5.04-5.07) back to upstream's 5.03.
%
% USAGE:  report = db_downgrade_to_upstream(DbDir, SettingsFile)                    % dry run: report only
%         report = db_downgrade_to_upstream(DbDir, SettingsFile, 'DryRun', false)  % back up, then rewrite
%
% For every protocol folder in DbDir (data/protocol.mat) and for the settings file (brainstorm.mat):
%   1. Compares the derived surface lists saved at 5.0x with the lists rebuilt from the subject folders.
%   2. (DryRun=false) Copies protocol.mat and brainstorm.mat to DbDir/../db_downgrade_backup_<timestamp>/,
%      removes the derived fields, sets DbVersion = 5.03, saves (v7).
% Derived .mat files are never touched. Run with Brainstorm closed.

    options = struct('DryRun', true);
    for iArg = 1:2:numel(varargin)
        options.(varargin{iArg}) = varargin{iArg + 1};
    end
    upstreamVersion = 5.03;
    backupFolder = '';
    if ~options.DryRun
        backupFolder = fullfile(fileparts(DbDir), ['db_downgrade_backup_' datestr(now, 'yymmdd_HHMMSS')]);
        mkdir(backupFolder);
        copyfile(SettingsFile, fullfile(backupFolder, 'brainstorm.mat'));
    end

    report = repmat(struct('Protocol', '', 'OldVersion', [], 'ListedBefore', 0, 'RebuiltFromDisk', 0, ...
        'MissingAfter', {{}}, 'Orphans', {{}}, 'BackupFolder', backupFolder), 1, 0);
    protocolFolders = dir(DbDir);
    protocolFolders = protocolFolders([protocolFolders.isdir] & ~strncmp({protocolFolders.name}, '.', 1));
    for iProtocol = 1:numel(protocolFolders)
        protocolName = protocolFolders(iProtocol).name;
        protocolFile = fullfile(DbDir, protocolName, 'data', 'protocol.mat');
        subjectsRoot = fullfile(DbDir, protocolName, 'anat');
        if ~exist(protocolFile, 'file')
            continue;
        end
        protocolMat = load(protocolFile);
        row = struct('Protocol', protocolName, 'OldVersion', protocolMat.DbVersion, 'ListedBefore', 0, ...
            'RebuiltFromDisk', 0, 'MissingAfter', {{}}, 'Orphans', {{}}, 'BackupFolder', backupFolder);
        % Saved lists vs lists rebuilt from disk
        listedBefore = ListDerivedFiles(protocolMat.ProtocolSubjects);
        rebuilt = protocolMat.ProtocolSubjects;
        orphanOutput = evalc('rebuilt = db_derived_nodes(''Attach'', rebuilt, subjectsRoot);');
        rebuiltFiles = ListDerivedFiles(rebuilt);
        row.ListedBefore    = numel(listedBefore);
        row.RebuiltFromDisk = numel(rebuiltFiles);
        row.MissingAfter    = setdiff(listedBefore, rebuiltFiles);
        row.Orphans         = regexp(orphanOutput, '(?<=Derived file not attached: )\S+', 'match');
        if ~options.DryRun && (protocolMat.DbVersion > upstreamVersion)
            mkdir(fullfile(backupFolder, protocolName));
            copyfile(protocolFile, fullfile(backupFolder, protocolName, 'protocol.mat'));
            protocolMat.ProtocolSubjects = db_derived_nodes('Strip', protocolMat.ProtocolSubjects);
            protocolMat.DbVersion = upstreamVersion;
            save(protocolFile, '-struct', 'protocolMat', '-v7');
        end
        report(end+1) = row; %#ok<AGROW>
        fprintf('%-32s %.2f  listed %3d  rebuilt %3d  missing %3d  orphans %3d\n', protocolName, row.OldVersion, ...
            row.ListedBefore, row.RebuiltFromDisk, numel(row.MissingAfter), numel(row.Orphans));
    end

    if ~options.DryRun
        settingsMat = load(SettingsFile);
        settingsMat.ProtocolsListSubjects = db_derived_nodes('Strip', settingsMat.ProtocolsListSubjects);
        settingsMat.DbVersion = upstreamVersion;
        save(SettingsFile, '-struct', 'settingsMat', '-v7');
        fprintf('Backup: %s\n', backupFolder);
    end
end

% Relative file names of every derived entry listed in one protocol's subjects
function derivedFiles = ListDerivedFiles(sProtocolSubjects)
    derivedFiles = {};
    allSubjects = sProtocolSubjects.Subject;
    if ~isempty(sProtocolSubjects.DefaultSubject)
        allSubjects = [{sProtocolSubjects.DefaultSubject}, num2cell(allSubjects)];
    else
        allSubjects = num2cell(allSubjects);
    end
    for iSubject = 1:numel(allSubjects)
        sSurf = allSubjects{iSubject}.Surface;
        for listName = {'Manifold', 'Operator', 'Eigen', 'Wavelet'}
            if isfield(sSurf, listName{1})
                for iSurf = 1:numel(sSurf)
                    if ~isempty(sSurf(iSurf).(listName{1}))
                        derivedFiles = [derivedFiles, {sSurf(iSurf).(listName{1}).FileName}]; %#ok<AGROW>
                    end
                end
            end
        end
    end
    derivedFiles = unique(derivedFiles);
end
```

- [ ] **Step 2: Dry run on the Task 5 copy**

Run: `matlab -batch "cd ~/workspace/research/code/brainstorm3@refactor-upstream-db-schema; brainstorm setpath; addpath('dev'); r = db_downgrade_to_upstream('/Volumes/Pergammon/pergammon/tmp/upstream-db-schema/brainstorm_db', '/Volumes/Pergammon/pergammon/tmp/upstream-db-schema/settings/.brainstorm/brainstorm.mat');"`
Expected: one line per protocol. `missing 0` for the copied protocol. Any `missing` or `orphans` is investigated before Step 3: open the file and check its `ParentSurface`.

- [ ] **Step 3: Real run on the copy, then repeat Task 5 Steps 3–4**

Run the same with `'DryRun', false`. Expected: the backup folder is printed, `protocol.mat` and `brainstorm.mat` show `DbVersion` 5.03, and Task 5 Steps 3–4 pass on the copy.

- [ ] **Step 4: Commit**

```bash
git add dev/db_downgrade_to_upstream.m
git commit -m "feat(dev): db_downgrade_to_upstream — back up, report, and take a fork-migrated database back to 5.03"
```

---

## After the plan (decisions for Diellor, not tasks)

1. **Merge** `refactor/upstream-db-schema` into `development` (PR on the fork).
2. **Other branches at 5.07** (`fix/preventad-import-inputstruct`, `feat/preventad-multimodal`, `fix/fibers-acpc-registration`, `feature/import-fibers-pipeline`, `feature/source-flow-processes`, `refactor/pet-metadata`): merge `development` into those still in use, and retire the rest. Any of them run on the database after the downgrade would raise it to 5.07 again.
3. **Real downgrade** (needs his OK): with Brainstorm closed everywhere, dry-run `db_downgrade_to_upstream` on `/Volumes/SpikeData-2/workspace/library/datasets/brainstorm_db` and `~/.brainstorm/brainstorm.mat`, show him the report, then run it with `'DryRun', false`.
4. **Then** load `TutorialPET` into the protocol list.
