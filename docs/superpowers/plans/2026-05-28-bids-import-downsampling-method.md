# BIDS Import — Cortex Downsampling Method Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let `process_import_bids` choose the cortex downsampling method (reducepatch vs icosphere) and default it to icosphere/ico5 (20484-vertex manifold cortex), threaded through to `import_anatomy_fs` for FreeSurfer subjects.

**Architecture:** All changes are in one file, `toolbox/process/functions/process_import_bids.m`. Two new GUI radios (`downsamplemethod`, `icolevel`) feed two pure, macro-dispatched helpers — `GetIcoVertexCount` (level→total vertices) and `ResolveAnatDownsample` (format+options→`nVertices`/`Method`/warning) — whose outputs are passed into the existing per-subject anatomy import. Icosphere is FreeSurfer-only; other formats fall back to reducepatch with a warning.

**Tech Stack:** MATLAB; Brainstorm process framework (`eval(macro_method)` dispatch); MATLAB MCP tools `evaluate_matlab_code` (run a test function by name), `check_matlab_code` (static analysis). Tests are script-style functions in `dev/tests/` that print `ALL TESTS PASSED`.

**Spec:** `docs/superpowers/specs/2026-05-28-bids-import-downsampling-method-design.md`

**Reference reading (do not modify):**
- `toolbox/io/import_anatomy_fs.m:1` — signature `(iSubject, FsDir, nVertices, isInteractive, sFid, isExtraMaps=0, isVolumeAtlas=1, isKeepMri=0, Method=[])`.
- `toolbox/io/import_anatomy_fs.m:141-156` — non-interactive icosphere path: `nVertHemi = round(nVertices/2)`, `tess_downsize` snaps to nearest ico count.
- `dev/ico-downsize.md` — ico levels: ico3/4/5/6 = 642/2562/10242/40962 per hemisphere (×2 = 1284/5124/20484/81924 total).

---

## Task 1: GUI options + `GetIcoVertexCount` helper

**Files:**
- Modify: `toolbox/process/functions/process_import_bids.m` (GetDescription ~lines 94-97; add subfunction)
- Test: `dev/tests/test_process_import_bids_options.m` (create)

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_process_import_bids_options.m`:

```matlab
function test_process_import_bids_options
% Verify the BIDS importer exposes the cortex downsampling method/level options
% (default icosphere/ico5) and the GetIcoVertexCount level->count helper.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

sProcess = process_import_bids('GetDescription');

% New options exist, old one preserved
assert(isfield(sProcess.options, 'downsamplemethod'), 'Missing option: downsamplemethod');
assert(isfield(sProcess.options, 'icolevel'),         'Missing option: icolevel');
assert(isfield(sProcess.options, 'nvertices'),        'nvertices option must still exist');

% Correct widget types
assert(strcmp(sProcess.options.downsamplemethod.Type, 'radio_linelabel'), 'downsamplemethod must be radio_linelabel.');
assert(strcmp(sProcess.options.icolevel.Type,         'radio_linelabel'), 'icolevel must be radio_linelabel.');

% Defaults: icosphere / ico5
assert(strcmp(sProcess.options.downsamplemethod.Value, 'icosphere'), 'Default downsamplemethod must be icosphere.');
assert(strcmp(sProcess.options.icolevel.Value,         'ico5'),      'Default icolevel must be ico5.');

% GetIcoVertexCount mappings (total cortex vertices)
assert(process_import_bids('GetIcoVertexCount', 'ico3') == 1284,  'ico3 -> 1284');
assert(process_import_bids('GetIcoVertexCount', 'ico4') == 5124,  'ico4 -> 5124');
assert(process_import_bids('GetIcoVertexCount', 'ico5') == 20484, 'ico5 -> 20484');
assert(process_import_bids('GetIcoVertexCount', 'ico6') == 81924, 'ico6 -> 81924');

% Unknown level errors
threw = false;
try, process_import_bids('GetIcoVertexCount', 'ico9'); catch, threw = true; end
assert(threw, 'Unknown ico level should error.');

fprintf('ALL TESTS PASSED: test_process_import_bids_options\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`): `test_process_import_bids_options`
Expected: FAIL — `downsamplemethod` field missing (and/or `GetIcoVertexCount` errors as an unknown subfunction).

- [ ] **Step 3: Add the two GUI radios + relabel `nvertices`**

In `process_import_bids.m` GetDescription, replace the existing "Option: Number of vertices" block (currently):

```matlab
    % Option: Number of vertices
    sProcess.options.nvertices.Comment = 'Number of vertices (cortex): ';
    sProcess.options.nvertices.Type    = 'value';
    sProcess.options.nvertices.Value   = {15000, '', 0};
```

with:

```matlab
    % Cortex downsampling method (icosphere is FreeSurfer-only; reducepatch is the legacy path)
    sProcess.options.downsamplemethod.Comment = {'Reducepatch', 'Icosphere', 'Cortex downsampling:'; ...
                                                 'reducepatch', 'icosphere', ''};
    sProcess.options.downsamplemethod.Type    = 'radio_linelabel';
    sProcess.options.downsamplemethod.Value   = 'icosphere';
    % Icosphere resolution level (used only when method = icosphere)
    sProcess.options.icolevel.Comment = {'ico3', 'ico4', 'ico5', 'ico6', 'Icosphere level:'; ...
                                         'ico3', 'ico4', 'ico5', 'ico6', ''};
    sProcess.options.icolevel.Type    = 'radio_linelabel';
    sProcess.options.icolevel.Value   = 'ico5';
    % Option: Number of vertices (reducepatch path only)
    sProcess.options.nvertices.Comment = 'Number of vertices (cortex, reducepatch): ';
    sProcess.options.nvertices.Type    = 'value';
    sProcess.options.nvertices.Value   = {15000, '', 0};
```

- [ ] **Step 4: Add the `GetIcoVertexCount` subfunction**

Add this subfunction to `process_import_bids.m` (place it after `GetDescription`/`FormatComment`, before `Run`):

```matlab
%% ===== GET ICOSPHERE VERTEX COUNT =====
% Total cortex vertex count (both hemispheres) for a FreeSurfer/MNE icosphere level.
% import_anatomy_fs takes round(nVertices/2) per hemisphere and tess_downsize snaps to
% the nearest ico count, so returning the exact total makes the chosen level explicit.
function n = GetIcoVertexCount(level) %#ok<DEFNU>
    switch lower(level)
        case 'ico3', n = 1284;
        case 'ico4', n = 5124;
        case 'ico5', n = 20484;
        case 'ico6', n = 81924;
        otherwise,   error('Unknown icosphere level: %s (expected ico3/ico4/ico5/ico6).', level);
    end
end
```

- [ ] **Step 5: Run the test to verify it passes**

Run (MATLAB MCP `evaluate_matlab_code`): `test_process_import_bids_options`
Expected: `ALL TESTS PASSED: test_process_import_bids_options`

- [ ] **Step 6: Static analysis**

Run (MATLAB MCP `check_matlab_code`) on `toolbox/process/functions/process_import_bids.m`.
Expected: no new errors (pre-existing style lint only).

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_import_bids.m dev/tests/test_process_import_bids_options.m
git commit -m "BIDS import: add cortex downsampling method/level options (default icosphere/ico5)"
```

---

## Task 2: `ResolveAnatDownsample` pure helper

**Files:**
- Modify: `toolbox/process/functions/process_import_bids.m` (add subfunction)
- Test: `dev/tests/test_import_bids_resolve_pure.m` (create)

This helper centralizes Component 5's decision so the cross-format branching is unit-testable without a full import. It depends on `GetIcoVertexCount` from Task 1.

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_import_bids_resolve_pure.m`:

```matlab
function test_import_bids_resolve_pure
% Verify ResolveAnatDownsample decides (nVertices, Method, warning) per anatomy format:
%   FreeSurfer + icosphere   -> ico total count + 'icosphere', no warning
%   FreeSurfer + reducepatch -> passthrough nVertices + 'reducepatch', no warning
%   non-FreeSurfer + icosphere -> reducepatch fallback + non-empty warning
%   non-FreeSurfer + reducepatch -> passthrough, no warning
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status')
    brainstorm nogui
end

% FreeSurfer + icosphere/ico5 -> 20484, icosphere, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'FreeSurfer', 'icosphere', 'ico5', 15000);
assert(nv == 20484 && strcmp(m, 'icosphere') && isempty(w), 'FreeSurfer icosphere ico5 must give 20484/icosphere/no-warning.');

% FreeSurfer + reducepatch -> passthrough 15000, reducepatch, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'FreeSurfer', 'reducepatch', 'ico5', 15000);
assert(nv == 15000 && strcmp(m, 'reducepatch') && isempty(w), 'FreeSurfer reducepatch must pass through 15000/reducepatch.');

% CAT12 + icosphere -> reducepatch fallback + warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'CAT12', 'icosphere', 'ico5', 15000);
assert(nv == 15000 && strcmp(m, 'reducepatch') && ~isempty(w), 'Non-FreeSurfer icosphere must fall back to reducepatch with a warning.');

% CAT12 + reducepatch -> passthrough, no warning
[nv, m, w] = process_import_bids('ResolveAnatDownsample', 'CAT12', 'reducepatch', 'ico5', 12000);
assert(nv == 12000 && strcmp(m, 'reducepatch') && isempty(w), 'Non-FreeSurfer reducepatch must pass through, no warning.');

fprintf('ALL TESTS PASSED: test_import_bids_resolve_pure\n');
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `evaluate_matlab_code`): `test_import_bids_resolve_pure`
Expected: FAIL — `ResolveAnatDownsample` is an unknown subfunction.

- [ ] **Step 3: Add the `ResolveAnatDownsample` subfunction**

Add to `process_import_bids.m` (place directly after `GetIcoVertexCount`):

```matlab
%% ===== RESOLVE ANATOMY DOWNSAMPLE =====
% Decide the (nVertices, Method) to pass to the anatomy importer for one subject, given the
% requested downsampling and the subject's anatomy format. Icosphere is FreeSurfer-only; for any
% other format we fall back to reducepatch and return a non-empty warning string.
% DEFERRED: extend icosphere to import_anatomy_cat/bs/bv/civet (they share the same per-hemisphere
% reducepatch pattern). Pure / no side effects -> unit-testable.
function [nVertArg, methodArg, warnMsg] = ResolveAnatDownsample(anatFormat, downsampleMethod, icoLevel, nVertices) %#ok<DEFNU>
    warnMsg = '';
    if strcmpi(downsampleMethod, 'icosphere')
        if strcmpi(anatFormat, 'FreeSurfer')
            nVertArg  = GetIcoVertexCount(icoLevel);
            methodArg = 'icosphere';
        else
            % Icosphere not implemented for this importer yet -> reducepatch fallback + warning
            nVertArg  = nVertices;
            methodArg = 'reducepatch';
            warnMsg   = sprintf(['Icosphere downsampling is currently FreeSurfer-only; ' ...
                'anatomy format "%s" imported with reducepatch at %d vertices.'], anatFormat, nVertices);
        end
    else
        nVertArg  = nVertices;
        methodArg = 'reducepatch';
    end
end
```

- [ ] **Step 4: Run the test to verify it passes**

Run (MATLAB MCP `evaluate_matlab_code`): `test_import_bids_resolve_pure`
Expected: `ALL TESTS PASSED: test_import_bids_resolve_pure`

- [ ] **Step 5: Static analysis**

Run (MATLAB MCP `check_matlab_code`) on `toolbox/process/functions/process_import_bids.m`.
Expected: no new errors (pre-existing style lint only).

- [ ] **Step 6: Commit**

```bash
git add toolbox/process/functions/process_import_bids.m dev/tests/test_import_bids_resolve_pure.m
git commit -m "BIDS import: add ResolveAnatDownsample helper (FreeSurfer-only icosphere, reducepatch fallback)"
```

---

## Task 3: Wire options into Run, defaults, interactive ask, and the import call site

**Files:**
- Modify: `toolbox/process/functions/process_import_bids.m` (Run ~line 142; `Def_OPTIONS` ~line 180; interactive ask ~lines 468-475; anatomy `switch` ~lines 484-497)

No new unit test (this is integration wiring; verified by static analysis, the unchanged pure tests, and the live check in Task 4).

- [ ] **Step 1: Read the two new options in `Run`**

In `Run`, immediately after the line `OPTIONS.RegisterMethod   = sProcess.options.anatregister.Value;`, add:

```matlab
    OPTIONS.DownsampleMethod = sProcess.options.downsamplemethod.Value;
    OPTIONS.IcoLevel         = sProcess.options.icolevel.Value;
```

- [ ] **Step 2: Add defaults to `Def_OPTIONS`**

In `ImportBidsDataset`, change the `Def_OPTIONS` struct (currently ending `'RegisterMethod', 'spm12');`) to:

```matlab
    Def_OPTIONS = struct(...
        'nVertices',        [], ...
        'isInteractive',    1, ...
        'ChannelAlign',     0, ...
        'SelectedSubjects', [], ...
        'isGroupSessions',  1, ...
        'MniMethod',        'maff8', ...  % {'maff8','segment','no'}
        'RegisterMethod',   'spm12', ...  % {'smp12','no'}
        'DownsampleMethod', 'icosphere', ...  % {'reducepatch','icosphere'}
        'IcoLevel',         'ico5');          % {'ico3','ico4','ico5','ico6'}
```

- [ ] **Step 3: Replace the interactive vertex-count ask with a method/level ask**

In `ImportBidsDataset`, replace the existing block (currently):

```matlab
            % Ask for number of vertices (so it is not asked multiple times)
            if isempty(OPTIONS.nVertices)
                OPTIONS.nVertices = java_dialog('input', 'Number of vertices on the cortex surface:', 'Import FreeSurfer folder', [], '15000');
                if isempty(OPTIONS.nVertices)
                    return;
                end
                OPTIONS.nVertices = str2double(OPTIONS.nVertices);
            end
```

with:

```matlab
            % Ask the cortex downsampling method/level once (so it is not asked per subject).
            % nVertices empty is the "not yet specified" sentinel for interactive callers.
            if isempty(OPTIONS.nVertices)
                iMethod = java_dialog('radio', 'Cortex downsampling method:', 'Import BIDS dataset', [], ...
                    {'<HTML><B>reducepatch</B>: Matlab decimation to ~N vertices', ...
                     '<HTML><B>icosphere</B>: FreeSurfer/MNE-style uniform grid (FreeSurfer only)'}, 2);
                if isempty(iMethod)
                    return;
                end
                if (iMethod == 2)
                    OPTIONS.DownsampleMethod = 'icosphere';
                    iLevel = java_dialog('radio', 'Icosphere resolution (total cortex vertices):', 'Import BIDS dataset', [], ...
                        {'ico3  -   1284 vertices', ...
                         'ico4  -   5124 vertices', ...
                         'ico5  -  20484 vertices', ...
                         'ico6  -  81924 vertices'}, 3);
                    if isempty(iLevel)
                        return;
                    end
                    icoVals = {'ico3', 'ico4', 'ico5', 'ico6'};
                    OPTIONS.IcoLevel  = icoVals{iLevel};
                    OPTIONS.nVertices = GetIcoVertexCount(OPTIONS.IcoLevel);   % satisfies the "asked once" sentinel
                else
                    OPTIONS.DownsampleMethod = 'reducepatch';
                    OPTIONS.nVertices = java_dialog('input', 'Number of vertices on the cortex surface:', 'Import BIDS dataset', [], '15000');
                    if isempty(OPTIONS.nVertices)
                        return;
                    end
                    OPTIONS.nVertices = str2double(OPTIONS.nVertices);
                end
            end
```

- [ ] **Step 4: Resolve downsampling at the import call site and thread `Method`**

In `ImportBidsDataset`, replace the existing anatomy `switch` block (currently):

```matlab
            % Import subject anatomy
            switch (SubjectAnatFormat{iSubj})
                case 'FreeSurfer'
                    errorMsg = import_anatomy_fs(iSubject, SubjectAnatDir{iSubj}, OPTIONS.nVertices, isInteractiveAnat, sMriFid, 0);
                case 'CAT12'
                    errorMsg = import_anatomy_cat(iSubject, SubjectAnatDir{iSubj}, OPTIONS.nVertices, isInteractiveAnat, sMriFid, 1, 2, 1);
                case 'BrainSuite'
                    errorMsg = import_anatomy_bs(iSubject, SubjectAnatDir{iSubj}, OPTIONS.nVertices, isInteractiveAnat, sMriFid);
                case 'BrainVISA'
                    errorMsg = import_anatomy_bv(iSubject, SubjectAnatDir{iSubj}, OPTIONS.nVertices, isInteractiveAnat, sMriFid);
                case 'CIVET'
                    errorMsg = import_anatomy_civet(iSubject, SubjectAnatDir{iSubj}, OPTIONS.nVertices, isInteractiveAnat, sMriFid, 0);
                otherwise
                    errorMsg = ['Invalid file format: ' SubjectAnatFormat{iSubj}];
            end
```

with:

```matlab
            % Resolve the cortex downsampling for this subject's anatomy format.
            % (Icosphere is FreeSurfer-only; other formats fall back to reducepatch + warning.)
            [nVertArg, methodArg, anatWarn] = ResolveAnatDownsample(SubjectAnatFormat{iSubj}, OPTIONS.DownsampleMethod, OPTIONS.IcoLevel, OPTIONS.nVertices);
            if ~isempty(anatWarn)
                anatWarn = [anatWarn ' (subject "' SubjectName{iSubj} '")'];
                Messages = [Messages, 10, anatWarn];
                disp(['BST> ' anatWarn]);
            end
            % Import subject anatomy
            switch (SubjectAnatFormat{iSubj})
                case 'FreeSurfer'
                    errorMsg = import_anatomy_fs(iSubject, SubjectAnatDir{iSubj}, nVertArg, isInteractiveAnat, sMriFid, 0, 1, 0, methodArg);
                case 'CAT12'
                    errorMsg = import_anatomy_cat(iSubject, SubjectAnatDir{iSubj}, nVertArg, isInteractiveAnat, sMriFid, 1, 2, 1);
                case 'BrainSuite'
                    errorMsg = import_anatomy_bs(iSubject, SubjectAnatDir{iSubj}, nVertArg, isInteractiveAnat, sMriFid);
                case 'BrainVISA'
                    errorMsg = import_anatomy_bv(iSubject, SubjectAnatDir{iSubj}, nVertArg, isInteractiveAnat, sMriFid);
                case 'CIVET'
                    errorMsg = import_anatomy_civet(iSubject, SubjectAnatDir{iSubj}, nVertArg, isInteractiveAnat, sMriFid, 0);
                otherwise
                    errorMsg = ['Invalid file format: ' SubjectAnatFormat{iSubj}];
            end
```

Note: the FreeSurfer call gains args 7/8/9 = `0, 1, 0, methodArg` (isExtraMaps=0, isVolumeAtlas=1, isKeepMri=0) — the verified current effective defaults, so reducepatch behavior is preserved exactly. For non-FreeSurfer formats `methodArg` is always `'reducepatch'` and `nVertArg == OPTIONS.nVertices`, so those calls behave exactly as before.

- [ ] **Step 5: Static analysis**

Run (MATLAB MCP `check_matlab_code`) on `toolbox/process/functions/process_import_bids.m`.
Expected: no new errors (pre-existing style lint only). In particular, confirm no "undefined function/variable" for `nVertArg`, `methodArg`, `anatWarn`, `SubjectName`.

- [ ] **Step 6: Re-run both pure tests (regression)**

Run (MATLAB MCP `evaluate_matlab_code`): `test_process_import_bids_options` then `test_import_bids_resolve_pure`
Expected: both print `ALL TESTS PASSED`.

- [ ] **Step 7: Commit**

```bash
git add toolbox/process/functions/process_import_bids.m
git commit -m "BIDS import: thread downsampling method/level into anatomy import (Run, defaults, interactive ask, call site)"
```

---

## Task 4: Live verification on OMEGA sub-0002 (default import → ico5 manifold cortex)

**Files:**
- Create: `dev/tests/test_import_bids_ico_live.m` (a manual live check; not part of the automated suite)

This is the integration proof that the importer's **default** reaches `tess_downsize('icosphere', 10242)`. It runs the real (modified) `process_import_bids` on the local OMEGA dataset limited to sub-0002, in a disposable protocol, then tears it down. It is heavyweight (links recordings) and depends on the local dataset; if the substrate is unavailable during execution, mark this task DONE_WITH_CONCERNS and defer to the OMEGA-rebuild increment that consumes this change.

- [ ] **Step 1: Write the live check script**

Create `dev/tests/test_import_bids_ico_live.m`:

```matlab
function results = test_import_bids_ico_live(varargin)
% TEST_IMPORT_BIDS_ICO_LIVE: Confirm process_import_bids defaults produce an ico5 manifold cortex.
% Runs the real importer on OMEGA sub-0002 (local) with DEFAULT options (no method set ->
% GetDescription default icosphere/ico5), checks the cortex, and deletes the disposable protocol.
%
% USAGE: results = test_import_bids_ico_live()
%        results = test_import_bids_ico_live('BidsDir', '/path/to/omega-tutorial')
Def.BidsDir      = '/Users/diellorbasha/workspace/library/datasets/omega-tutorial';
Def.ProtocolName = 'BidsIcoDefaultCheck';
Def.SubjTag      = 'sub-0002';
OPT = Def;
for i = 1:2:numel(varargin), OPT.(varargin{i}) = varargin{i+1}; end
assert(exist(OPT.BidsDir, 'dir') == 7, 'BIDS dataset folder not found: %s', OPT.BidsDir);

if ~brainstorm('status')
    brainstorm nogui
end

% Disposable, repeatable protocol
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);
gui_brainstorm('CreateProtocol', OPT.ProtocolName, 0, 0);
bst_report('Start');

% Run the real importer with DEFAULT downsampling (do NOT set downsamplemethod/icolevel)
bst_process('CallProcess', 'process_import_bids', [], [], ...
    'bidsdir',      {OPT.BidsDir, 'BIDS'}, ...
    'selectsubj',   OPT.SubjTag, ...
    'channelalign', 0);

% Find the imported subject's cortex and check it
ProtocolSubjects = bst_get('ProtocolSubjects');
iCortexSubj = [];
for iS = 1:numel(ProtocolSubjects.Subject)
    if ~isempty(ProtocolSubjects.Subject(iS).iCortex)
        iCortexSubj = iS; break;
    end
end
assert(~isempty(iCortexSubj), 'No subject with a cortex surface was imported.');
sSubject   = ProtocolSubjects.Subject(iCortexSubj);
CortexFile = sSubject.Surface(sSubject.iCortex).FileName;
TessMat    = in_tess_bst(CortexFile);

nVertices = size(TessMat.Vertices, 1);
nFaces    = size(TessMat.Faces, 1);
VertConn  = tess_vertconn(TessMat.Vertices, TessMat.Faces);
nComp     = max(conncomp(graph(VertConn)));
[~, ~, isManifold] = tess_manifold(TessMat.Vertices, TessMat.Faces);

results = struct('subject', sSubject.Name, 'cortexFile', CortexFile, ...
    'nVertices', nVertices, 'nFaces', nFaces, 'nComp', nComp, 'isManifold', logical(isManifold));
results.pass = (nVertices == 20484) && (nComp == 2) && isequal(logical(isManifold), true);

fprintf('Cortex: %s\n  vertices=%d (expect 20484), faces=%d, components=%d (expect 2), manifold=%d (expect 1)\n', ...
    sSubject.Name, nVertices, nFaces, nComp, isManifold);

% Tear down the disposable protocol
gui_brainstorm('DeleteProtocol', OPT.ProtocolName);

if results.pass
    fprintf('ALL TESTS PASSED: test_import_bids_ico_live\n');
else
    error('test_import_bids_ico_live FAILED: cortex did not match ico5/2-component/manifold.');
end
end
```

- [ ] **Step 2: Run the live check**

Run (MATLAB MCP `evaluate_matlab_code`): `test_import_bids_ico_live`
Expected: console shows `vertices=20484 ... components=2 ... manifold=1` and `ALL TESTS PASSED: test_import_bids_ico_live`.

- [ ] **Step 3: Static analysis**

Run (MATLAB MCP `check_matlab_code`) on `dev/tests/test_import_bids_ico_live.m`.
Expected: no errors (pre-existing style lint only).

- [ ] **Step 4: Commit**

```bash
git add dev/tests/test_import_bids_ico_live.m
git commit -m "BIDS import: add live check that default import yields ico5 manifold cortex on OMEGA sub-0002"
```

---

## Self-review notes (author)

- **Spec coverage:** Component 1 (options) → Task 1; Component 2 (`GetIcoVertexCount`) → Task 1; Component 3 (Run wiring) → Task 3 Step 1; Component 4 (Def_OPTIONS + interactive ask) → Task 3 Steps 2-3; Component 5 (call site + non-FS warning) → Task 2 (decision logic) + Task 3 Step 4. Error handling → Task 1 (level error), Task 2/3 (warning), retained nVertices validation untouched. Testing items 1-4 → Tasks 1,2,3,4.
- **Plan-level refinement vs spec:** the spec inlined Component 5's decision at the call site; this plan extracts it into the pure `ResolveAnatDownsample` helper for testability. Same behavior, more testable — consistent with the spec's pure-helper philosophy.
- **Type consistency:** `GetIcoVertexCount(level)->n`; `ResolveAnatDownsample(anatFormat, downsampleMethod, icoLevel, nVertices)->[nVertArg, methodArg, warnMsg]`; GUI options `downsamplemethod`/`icolevel`; OPTIONS fields `DownsampleMethod`/`IcoLevel`. Used identically everywhere.
