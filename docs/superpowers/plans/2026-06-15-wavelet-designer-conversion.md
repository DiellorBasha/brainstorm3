# Wavelet Designer Conversion Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Convert the existing (conflated) FilterDesigner into a correct, localized **Wavelet Designer**: rename filter→wavelet across the GUI/DB/tree, make the input localization-only, and specify the seed direction in the local `manifold_` frame (U,V,N) instead of the world frame.

**Architecture:** A mechanical rename pass (filter→wavelet) preserving behavior, followed by two behavioral changes — drop the filter-style source-map input, and replace the world az/el seed direction with a local-frame embedding `d = cosθ(cosφ·U + sinφ·V) + sinθ·N` derived from the surface's manifold node. The Dirac math (`bst_dirac_filter`, `bst_dirac_eigenmodes_filter`, `bst_filterbank_tiles`) is unchanged — only where the seed direction comes from changes.

**Tech Stack:** MATLAB (Brainstorm toolbox), Brainstorm DB layer (`db_template`/`db_update`/`bst_get`/`db_add_*`), Java-Swing GUI, the `nxr`-backed `tess_manifold` + `view_manifold('DeriveVertexFrame')` frame derivation, and the `bst-java` fork (tree node type). Tests run via the MATLAB MCP against the live dev protocol (Subject01 `cortex_20484V` has a Dirac eigen node) and with small synthetic structs for pure logic.

**Repos / branches:** `brainstorm3` on `feat/filterbank-designer`; `bst-java` fork on a paired branch. Never push/merge the fork upstream.

**Conventions:** Tests live in `dev/tests/`, are plain functions that `error()` on failure, run with the MATLAB MCP `run_matlab_file`. Brainstorm panel functions dispatch subfunctions via `eval(macro_method)`, so `panel_wavelet_designer('Fn', args)` calls local `Fn`. Commit after each task with the message in its final step.

**Scope note — rename ONLY these (leave everything else with "filter"/"wavelet" in the name alone):** `panel_filter_designer`, `view_filter_designer`, `db_add_filterbank`, the `filterbankmat` template, the `Surface.Filterbank` list, `bst_get('FilterbankFile')`, the `filterbank` file type, the `filterbank` tree node + its `tree_callbacks`/`node_create_subject` wiring, the `FilterDesigner*` appdata keys in `figure_3d.m`, and the 3 tests (`test_filterbank_schema`, `test_db_add_filterbank`, `test_filter_designer_session`). **Do NOT rename:** `bst_filterbank_tiles` (still the spectrum-tiling module), `bst_dirac_filter` (the synthesis core), the `toolbox/math/eigfilter/*` library, `bst_sensor_cwt`, or `process_evt_detect_burst_wavelet`.

---

## File Structure

**Rename (git mv):**
- `toolbox/gui/panel_filter_designer.m` → `toolbox/gui/panel_wavelet_designer.m`
- `toolbox/gui/view_filter_designer.m` → `toolbox/gui/view_wavelet_designer.m`
- `toolbox/db/db_add_filterbank.m` → `toolbox/db/db_add_wavelet.m`
- `dev/tests/test_filterbank_schema.m` → `dev/tests/test_wavelet_schema.m`
- `dev/tests/test_db_add_filterbank.m` → `dev/tests/test_db_add_wavelet.m`
- `dev/tests/test_filter_designer_session.m` → `dev/tests/test_wavelet_designer_session.m`

**Modify:**
- `toolbox/db/db_template.m` — `filterbankmat`→`waveletmat`; `Surface.Filterbank`→`Surface.Wavelet`.
- `toolbox/db/db_update.m`, `toolbox/core/bst_startup.m` — DB v5.06 migration adding `Surface.Wavelet`.
- `toolbox/core/bst_get.m` — `FilterbankFile`→`WaveletFile`.
- `toolbox/io/file_gettype.m`, `toolbox/io/file_fullpath.m` — `filterbank`→`wavelet` type.
- `toolbox/tree/node_create_subject.m`, `toolbox/tree/tree_callbacks.m` — wavelet node nesting/menus.
- `toolbox/gui/figure_3d.m` — `FilterDesigner*` appdata keys → `WaveletDesigner*`.
- `~/workspace/research/code/bst-java/.../BstNode.java` (+ renderer) — `filterbank`→`wavelet` node type.

**New behavioral edits (in the renamed GUI files):**
- `panel_wavelet_designer.m` — remove input radios + `SetSeedSource`; `SeedDirection` uses the local frame via a pure `i_embed_direction` helper; relabel sliders.
- `view_wavelet_designer.m` — find-or-create the manifold, derive the frame, pass it to the panel.

---

## Task 1: DB schema + migration rename (filterbank → wavelet)

**Files:**
- Modify: `toolbox/db/db_template.m`
- Modify: `toolbox/core/bst_startup.m:208`
- Modify: `toolbox/db/db_update.m`
- Rename+Modify: `dev/tests/test_filterbank_schema.m` → `dev/tests/test_wavelet_schema.m`

- [ ] **Step 1: Rename the schema test and update it to assert the wavelet schema**

```bash
git mv dev/tests/test_filterbank_schema.m dev/tests/test_wavelet_schema.m
```

Replace the body of `dev/tests/test_wavelet_schema.m` with:

```matlab
function test_wavelet_schema()
% Schema regression for the waveletmat template and the Surface.Wavelet list.
% Authors: Diellor Basha, 2026
    nFail = 0;
    w = db_template('waveletmat');
    need = {'Comment','ParentEigen','Variant','Tiles','Tiling','Provenance'};
    for f = need
        if ~isfield(w, f{1}); fprintf('MISSING waveletmat.%s\n', f{1}); nFail = nFail+1; end
    end
    s = db_template('surface');
    if ~isfield(s, 'Wavelet'); fprintf('MISSING surface.Wavelet\n'); nFail = nFail+1; end
    if isfield(s,'Wavelet')
        sub = fieldnames(s.Wavelet);
        for f = {'FileName','Comment','ParentEigen'}
            if ~ismember(f{1}, sub); fprintf('MISSING surface.Wavelet.%s\n', f{1}); nFail = nFail+1; end
        end
        if ~isempty(s.Wavelet); fprintf('surface.Wavelet must start EMPTY\n'); nFail = nFail+1; end
    end
    if isfield(s, 'Filterbank'); fprintf('surface.Filterbank should be GONE\n'); nFail = nFail+1; end
    fprintf('\n==== test_wavelet_schema: %d failed ====\n', nFail);
    if nFail > 0, error('test_wavelet_schema FAILED'); end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_wavelet_schema.m`
Expected: FAIL — `Unknown data template : waveletmat` (the template is still `filterbankmat`).

- [ ] **Step 3: Rename the template in db_template.m**

In `toolbox/db/db_template.m`, change the `case 'filterbankmat'` line to `case 'waveletmat'` (the struct body — `Comment/ParentEigen/Variant/Tiles/Tiling/Provenance` — stays identical):

```matlab
    case 'waveletmat'
        template = struct(...
              'Comment',     '', ...
              'ParentEigen', '', ...   % file_short of the eigen_ node this wavelet applies in
              'Variant',     '', ...   % operator variant inherited from eigen (e.g. 'Dirac')
              'Tiles',       [], ...   % 1xN struct array (Kernel,Params,Direction,Chirality,Axis)
              'Tiling',      [], ...   % struct(Wavelet, Opts) — single design + tiling options
              'Provenance',  []);
```

In the `case 'surface'` block, rename the `Filterbank` field to `Wavelet`:

```matlab
                          'Filterbank',  struct('FileName',{},'Comment',{},'ParentEigen',{}));
```
becomes
```matlab
                          'Wavelet',     struct('FileName',{},'Comment',{},'ParentEigen',{}));
```

- [ ] **Step 4: Bump the DB version and add the 5.06 migration**

In `toolbox/core/bst_startup.m`, change `CurrentDbVersion = 5.05;` to:

```matlab
CurrentDbVersion = 5.06;
```

In `toolbox/db/db_update.m`, immediately after the closing `end` of the `if (CurrentDbVersion < 5.05)` block, add (mirrors the 5.05 block; `NormalizeSurfaceArray` adds whatever the current template has, i.e. `Wavelet`):

```matlab
if (CurrentDbVersion < 5.06)
    isSurfFixed = 0;
    templateSurface = db_template('Surface');
    for iProt = 1:length(ProtocolsListSubjects)
        subjFields = fieldnames(ProtocolsListSubjects(iProt));
        for iField = 1:length(subjFields)
            subjField = subjFields{iField};
            for iSubj = 1:length(ProtocolsListSubjects(iProt).(subjField))
                sSurf = ProtocolsListSubjects(iProt).(subjField)(iSubj).Surface;
                if isempty(sSurf)
                    continue;
                end
                [sSurf, isFix] = NormalizeSurfaceArray(sSurf, templateSurface);
                if isFix
                    ProtocolsListSubjects(iProt).(subjField)(iSubj).Surface = sSurf;
                    isSurfFixed = 1;
                end
            end
        end
    end
    if isSurfFixed
        disp('BST> Database structure: Adding wavelet support to surfaces...');
        SaveProtocolSubjects();
        disp('BST> Database structure: Done.');
    end
end
```

- [ ] **Step 5: Run the schema test + apply the migration to the live protocol**

Run (MCP `evaluate_matlab_code`):

```matlab
rehash;
db_update(5.06);
sSubject = bst_get('Subject', 1);
assert(isfield(sSubject.Surface, 'Wavelet'), 'Wavelet field not backfilled');
fprintf('OK: Surface.Wavelet present\n');
```
Expected: prints `OK: Surface.Wavelet present` (and the "Adding wavelet support" line the first time).

Then run `dev/tests/test_wavelet_schema.m`.
Expected: PASS — `0 failed`.

- [ ] **Step 6: Commit**

```bash
git add toolbox/db/db_template.m toolbox/core/bst_startup.m toolbox/db/db_update.m dev/tests/test_wavelet_schema.m
git commit -m "refactor(db): rename filterbankmat->waveletmat, Surface.Filterbank->Wavelet (v5.06)"
```

---

## Task 2: File type + accessor + db_add rename

**Files:**
- Modify: `toolbox/io/file_gettype.m`, `toolbox/io/file_fullpath.m`
- Modify: `toolbox/core/bst_get.m`
- Rename+Modify: `toolbox/db/db_add_filterbank.m` → `toolbox/db/db_add_wavelet.m`
- Rename+Modify: `dev/tests/test_db_add_filterbank.m` → `dev/tests/test_db_add_wavelet.m`

- [ ] **Step 1: Rename the db_add round-trip test and update it**

```bash
git mv dev/tests/test_db_add_filterbank.m dev/tests/test_db_add_wavelet.m
```

Replace the body of `dev/tests/test_db_add_wavelet.m` with:

```matlab
function test_db_add_wavelet()
% Round-trip: save a wavelet under an eigen node, resolve it, reload it, delete it.
% Requires Brainstorm running with a Dirac eigen node (Subject01 surface 5).
% Authors: Diellor Basha, 2026
    nFail = 0;
    sSubject = bst_get('Subject', 1);
    iSurf = 5;
    assert(isfield(sSubject.Surface(iSurf),'Eigen') && ~isempty(sSubject.Surface(iSurf).Eigen), ...
        'No eigen node on Subject01 surface 5; compute one first.');
    EigenFile = sSubject.Surface(iSurf).Eigen(1).FileName;

    wavelet = struct('Kernel','mexhat','Params',struct('t',0.01),'Direction',[1 0 0], ...
                     'Chirality',0,'Axis',[0 0 1]);
    opts    = struct('N',4,'Spacing','geometric','LambdaRange',[1 256],'Chiralities',[]);
    w = db_template('waveletmat');
    w.ParentEigen = file_short(EigenFile);
    w.Variant     = 'Dirac';
    w.Tiles       = bst_filterbank_tiles(wavelet, opts);
    w.Tiling      = struct('Wavelet', wavelet, 'Opts', opts);

    iW = db_add_wavelet(1, EigenFile, w, 'TEST wavelet');
    nFail = nFail + chk('returns an index', ~isempty(iW));

    sSubject = bst_get('Subject', 1);
    ws = sSubject.Surface(iSurf).Wavelet;
    nFail = nFail + chk('appears in Surface.Wavelet', ~isempty(ws));
    nFail = nFail + chk('ParentEigen matches', any(strcmp({ws.ParentEigen}, file_short(EigenFile))));

    newFile = ws(end).FileName;
    [~,~,iSurf2,iW2] = bst_get('WaveletFile', newFile);
    nFail = nFail + chk('accessor resolves it', ~isempty(iW2) && iSurf2==iSurf);
    R = load(file_fullpath(newFile));
    nFail = nFail + chk('reloaded Tiles count', numel(R.Tiles)==numel(w.Tiles));
    nFail = nFail + chk('reloaded Tiling.Wavelet', strcmp(R.Tiling.Wavelet.Kernel,'mexhat'));

    file_delete(file_fullpath(newFile), 1);
    sSubject.Surface(iSurf).Wavelet(end) = [];
    ProtocolSubjects = bst_get('ProtocolSubjects');
    ProtocolSubjects.Subject(1) = sSubject;
    bst_set('ProtocolSubjects', ProtocolSubjects);
    db_save();

    fprintf('\n==== test_db_add_wavelet: %d failed ====\n', nFail);
    if nFail > 0, error('test_db_add_wavelet FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/tests/test_db_add_wavelet.m`
Expected: FAIL — `Undefined function 'db_add_wavelet'`.

- [ ] **Step 3: Register the `wavelet` file type**

In `toolbox/io/file_gettype.m`, change the filterbank branch:

```matlab
        elseif ~isempty(strfind(fileName, '_filterbank_'))
            fileType = 'filterbank';
```
to
```matlab
        elseif ~isempty(strfind(fileName, '_wavelet_'))
            fileType = 'wavelet';
```

In `toolbox/io/file_fullpath.m`, change `'filterbank'` to `'wavelet'` in the anatomy `case {...}` list (the line containing `'manifold', 'operator', 'eigen', 'filterbank'`):

```matlab
    case {'brainstormsubject', 'subject', 'subjectimage', 'anatomy', 'scalp', 'outerskull', 'innerskull', 'cortex', 'fibers', 'fem', 'other', 'tess', 'manifold', 'operator', 'eigen', 'wavelet'}
```

- [ ] **Step 4: Rename the bst_get accessor**

In `toolbox/core/bst_get.m`, find `case 'FilterbankFile'` and replace the whole block by swapping `Filterbank`→`Wavelet` and the local var `FilterbankFile`→`WaveletFile` throughout:

```matlab
    case 'WaveletFile'
        if isempty(GlobalData) || isempty(GlobalData.DataBase) || isempty(GlobalData.DataBase.iProtocol) || (GlobalData.DataBase.iProtocol == 0)
            return;
        end
        ProtocolSubjects = GlobalData.DataBase.ProtocolSubjects(GlobalData.DataBase.iProtocol);
        if isempty(ProtocolSubjects)
            return
        end
        if (nargin == 2)
            WaveletFile = varargin{2};
        else
            error('Invalid call to bst_get().');
        end
        WaveletFile = file_short(WaveletFile);
        if ~isempty(ProtocolSubjects.DefaultSubject)
            for iSurf = 1:length(ProtocolSubjects.DefaultSubject.Surface)
                if isfield(ProtocolSubjects.DefaultSubject.Surface(iSurf), 'Wavelet') && ...
                        ~isempty(ProtocolSubjects.DefaultSubject.Surface(iSurf).Wavelet)
                    iW = find(file_compare(WaveletFile, {ProtocolSubjects.DefaultSubject.Surface(iSurf).Wavelet.FileName}), 1);
                    if ~isempty(iW)
                        argout1 = ProtocolSubjects.DefaultSubject; argout2 = 0; argout3 = iSurf; argout4 = iW;
                        return
                    end
                end
            end
        end
        for iSubj = 1:length(ProtocolSubjects.Subject)
            for iSurf = 1:length(ProtocolSubjects.Subject(iSubj).Surface)
                if isfield(ProtocolSubjects.Subject(iSubj).Surface(iSurf), 'Wavelet') && ...
                        ~isempty(ProtocolSubjects.Subject(iSubj).Surface(iSurf).Wavelet)
                    iW = find(file_compare(WaveletFile, {ProtocolSubjects.Subject(iSubj).Surface(iSurf).Wavelet.FileName}), 1);
                    if ~isempty(iW)
                        argout1 = ProtocolSubjects.Subject(iSubj); argout2 = iSubj; argout3 = iSurf; argout4 = iW;
                        return
                    end
                end
            end
        end
```

- [ ] **Step 5: Rename db_add_filterbank → db_add_wavelet**

```bash
git mv toolbox/db/db_add_filterbank.m toolbox/db/db_add_wavelet.m
```

In `toolbox/db/db_add_wavelet.m`: rename the function to `db_add_wavelet`, the output var `iFilterbank`→`iWavelet`, the output filename prefix `'filterbank_'`→`'wavelet_'`, and the DB list field `.Filterbank`→`.Wavelet`. The substantive lines become:

```matlab
function iWavelet = db_add_wavelet(iSubject, ParentEigenFile, WaveletMat, Comment)
% DB_ADD_WAVELET: Save a wavelet_*.mat and register it as a child of an eigen node.
%   ...
    if (nargin < 4) || isempty(Comment); Comment = 'Wavelet'; end
    [~, iSubjectE, iSurface] = bst_get('EigenFile', ParentEigenFile);
    if isempty(iSurface)
        error('db_add_wavelet:noEigen', 'Parent eigen node not found: %s', ParentEigenFile);
    end
    iSubject = iSubjectE;
    ProtocolSubjects = bst_get('ProtocolSubjects');
    if (iSubject == 0); sSubject = ProtocolSubjects.DefaultSubject;
    else;               sSubject = ProtocolSubjects.Subject(iSubject); end
    ProtocolInfo = bst_get('ProtocolInfo');
    c = clock;
    strTime = sprintf('%02.0f%02.0f%02.0f_%02.0f%02.0f', c(1)-2000, c(2:5));
    OutputFile = ['wavelet_' strTime '.mat'];
    OutputFileFull = file_unique(bst_fullfile(ProtocolInfo.SUBJECTS, bst_fileparts(sSubject.FileName), OutputFile));
    WaveletMat.ParentEigen = file_short(ParentEigenFile);
    WaveletMat.Comment     = Comment;
    bst_save(OutputFileFull, WaveletMat, 'v7');
    templateSurface = db_template('Surface');
    tFields = fieldnames(templateSurface);
    if ~isempty(sSubject.Surface) && ~all(isfield(sSubject.Surface, tFields))
        sNew = cell(1, numel(sSubject.Surface));
        for k = 1:numel(sSubject.Surface)
            s = struct_copy_fields(sSubject.Surface(k), templateSurface, 0);
            extraFields = setdiff(fieldnames(s), tFields, 'stable');
            sNew{k} = orderfields(s, [tFields; extraFields]);
        end
        sSubject.Surface = reshape([sNew{:}], size(sSubject.Surface));
    end
    newEntry.FileName    = file_short(OutputFileFull);
    newEntry.Comment     = Comment;
    newEntry.ParentEigen = file_short(ParentEigenFile);
    sSubject.Surface(iSurface).Wavelet(end+1) = newEntry;
    iWavelet = numel(sSubject.Surface(iSurface).Wavelet);
    if (iSubject == 0); ProtocolSubjects.DefaultSubject = sSubject;
    else;               ProtocolSubjects.Subject(iSubject) = sSubject; end
    bst_set('ProtocolSubjects', ProtocolSubjects);
    panel_protocols('UpdateNode', 'Subject', iSubject);
    db_save();
end
```

- [ ] **Step 6: Run the round-trip test**

Run: `dev/tests/test_db_add_wavelet.m`
Expected: PASS — all `PASS`, `0 failed`, test cleans up its node.

- [ ] **Step 7: Commit**

```bash
git add toolbox/io/file_gettype.m toolbox/io/file_fullpath.m toolbox/core/bst_get.m toolbox/db/db_add_wavelet.m dev/tests/test_db_add_wavelet.m
git commit -m "refactor(db): wavelet file type + WaveletFile accessor + db_add_wavelet"
```

---

## Task 3: bst-java wavelet node type

**Files (bst-java fork):**
- Modify: `brainstorm/src/org/brainstorm/tree/BstNode.java` (or the file holding the `filterbank` type, if added earlier)

> The tree cell renderer already falls back to `ICON_OBJECT` for unregistered types (as eigen/operator/manifold do), so a brand-new `wavelet` type renders fine. If a `filterbank` branch was added to the fork earlier, rename it; otherwise this task only needs to ensure `CreateNode('wavelet', …)` is accepted.

- [ ] **Step 1: Create the paired fork branch**

```bash
cd ~/workspace/research/code/bst-java
git checkout -b feat/wavelet-node 2>/dev/null || git checkout feat/wavelet-node
```

- [ ] **Step 2: Find any existing 'filterbank' reference in the fork**

```bash
grep -rni "filterbank" ~/workspace/research/code/bst-java/brainstorm/src || echo "no filterbank in fork — default icon path is used"
```

- [ ] **Step 3: Rename/register the type**

If the grep found a `filterbank` branch (icon or tooltip), change the string `"filterbank"` → `"wavelet"`. If nothing was found, no Java change is required — `CreateNode('wavelet', …)` will render with the default `ICON_OBJECT`. (Do not invent new icon assets.)

- [ ] **Step 4: Rebuild the jar (only if Step 3 changed Java) and reload**

Follow the fork's existing build, deploy the jar where Brainstorm loads it, then restart Brainstorm (`brainstorm stop; brainstorm nogui` or full restart) so the JVM picks up the new class. If no Java changed, skip — just `db_reload_subjects(0)` later when the tree renders wavelet nodes.

- [ ] **Step 5: Commit (only if Java changed)**

```bash
cd ~/workspace/research/code/bst-java && git add -A && git commit -m "refactor(tree): wavelet node type (was filterbank)"
```

---

## Task 4: GUI + tree rename (FilterDesigner → WaveletDesigner)

**Files:**
- Rename: `toolbox/gui/panel_filter_designer.m` → `toolbox/gui/panel_wavelet_designer.m`
- Rename: `toolbox/gui/view_filter_designer.m` → `toolbox/gui/view_wavelet_designer.m`
- Modify: `toolbox/gui/figure_3d.m`, `toolbox/tree/node_create_subject.m`, `toolbox/tree/tree_callbacks.m`
- Rename+Modify: `dev/tests/test_filter_designer_session.m` → `dev/tests/test_wavelet_designer_session.m`

- [ ] **Step 1: git mv the GUI files and the session test**

```bash
git mv toolbox/gui/panel_filter_designer.m toolbox/gui/panel_wavelet_designer.m
git mv toolbox/gui/view_filter_designer.m  toolbox/gui/view_wavelet_designer.m
git mv dev/tests/test_filter_designer_session.m dev/tests/test_wavelet_designer_session.m
```

- [ ] **Step 2: Substitute identifiers in the two GUI files**

In BOTH `panel_wavelet_designer.m` and `view_wavelet_designer.m`, replace every occurrence:
- `panel_filter_designer` → `panel_wavelet_designer`
- `view_filter_designer` → `view_wavelet_designer`
- `'FilterDesigner'` → `'WaveletDesigner'` (panel name, the gui_show tab, all dispatch strings)
- `FilterDesignerState` → `WaveletDesignerState`
- `FilterDesignerPick` → `WaveletDesignerPick`
- `FilterDesignerTemp` → `WaveletDesignerTemp`
- `FilterDesignerStudy` → `WaveletDesignerStudy`
- `FilterDesignerDS` → `WaveletDesignerDS`
- `FilterDesignerClosing` → `WaveletDesignerClosing`
- `db_template('filterbankmat')` → `db_template('waveletmat')`
- `db_add_filterbank(` → `db_add_wavelet(`
- `bst_get('FilterbankFile'` → `bst_get('WaveletFile'`
- in `view_wavelet_designer.m`, the eigen-vs-node resolution: `bst_get('FilterbankFile', NodeFile)` → `bst_get('WaveletFile', NodeFile)`, and `loadBank.ParentEigen` usage unchanged.

Keep `bst_filterbank_tiles(` calls as-is (the module name is unchanged).

- [ ] **Step 3: Update the appdata keys in figure_3d.m**

In `toolbox/gui/figure_3d.m`, in the click-to-seed hook block, replace `FilterDesignerPick` → `WaveletDesignerPick`:

```matlab
        if isappdata(hFig, 'WaveletDesignerPick') && strcmpi(clickAction, 'rotate') ...
                && strcmpi(get(hFig, 'SelectionType'), 'normal')
            iVertex = panel_coordinates('SelectPoint', hFig, 0);
            if ~isempty(iVertex)
                pickFn = getappdata(hFig, 'WaveletDesignerPick');
                pickFn(iVertex);
            end
            return;
        end
```

- [ ] **Step 4: Update tree node creation (node_create_subject.m)**

In `toolbox/tree/node_create_subject.m`, the nested loop under each eigen node: rename `Filterbank`→`Wavelet` and `CreateNode('filterbank', …)`→`CreateNode('wavelet', …)`:

```matlab
                        if isfield(sSubject.Surface(iSurface), 'Wavelet')
                            eigName = char(sSubject.Surface(iSurface).Eigen(iE).FileName);
                            ws = sSubject.Surface(iSurface).Wavelet;
                            for iW = 1:numel(ws)
                                if file_compare(ws(iW).ParentEigen, eigName)
                                    [wCreated, wNode] = CreateNode('wavelet', ...
                                        char(ws(iW).Comment), char(ws(iW).FileName), ...
                                        iW, iSubject, iSearch);
                                    if wCreated
                                        chNode.add(wNode);
                                    end
                                end
                            end
                        end
```

- [ ] **Step 5: Update tree_callbacks.m menus + delete callback**

In `toolbox/tree/tree_callbacks.m`:
- The eigen popup item: `'Design filterbank...'` → `'Design wavelet...'`, and the callback `@(h,ev)bst_call(@view_filter_designer, filenameRelative)` → `@(h,ev)bst_call(@view_wavelet_designer, filenameRelative)`.
- The `case 'filterbank'` popup → `case 'wavelet'`, with "Edit / re-open" calling `view_wavelet_designer` and "Delete" calling `WaveletDelete_Callback`.
- Rename `FilterbankDelete_Callback` → `WaveletDelete_Callback`, swapping `bst_get('FilterbankFile', …)`→`bst_get('WaveletFile', …)` and `.Filterbank`→`.Wavelet`:

```matlab
function WaveletDelete_Callback(filenameRelative)
    if ~java_dialog('confirm', ['Delete wavelet?' 10 filenameRelative], 'Delete wavelet')
        return;
    end
    file_delete(file_fullpath(filenameRelative), 1);
    [~, iSubject, iSurface, iW] = bst_get('WaveletFile', filenameRelative);
    if ~isempty(iSubject)
        ProtocolSubjects = bst_get('ProtocolSubjects');
        if iSubject == 0
            ProtocolSubjects.DefaultSubject.Surface(iSurface).Wavelet(iW) = [];
        else
            ProtocolSubjects.Subject(iSubject).Surface(iSurface).Wavelet(iW) = [];
        end
        bst_set('ProtocolSubjects', ProtocolSubjects);
        db_save();
        panel_protocols('UpdateNode', 'Subject', iSubject);
    end
end
```
- In `EigenDelete_Callback`, the cascade block: `.Filterbank`→`.Wavelet` (delete child wavelets when the eigen node is deleted).

- [ ] **Step 6: Update the session test to the new names**

In `dev/tests/test_wavelet_designer_session.m`, replace `view_filter_designer`→`view_wavelet_designer`, `panel_filter_designer`→`panel_wavelet_designer`, `'FilterDesigner'`→`'WaveletDesigner'`, `.Filterbank`→`.Wavelet`, `bst_get('FilterbankFile'`→`bst_get('WaveletFile'`. The function name → `test_wavelet_designer_session`.

- [ ] **Step 7: Run the session test**

First reload the tree so wavelet nodes render: in MCP `evaluate_matlab_code` run `rehash; db_reload_subjects(0); panel_protocols('UpdateTree');`
Then run `dev/tests/test_wavelet_designer_session.m`.
Expected: PASS — `0 failed` (open → seed → save → nested `wavelet_` node → delete → teardown).

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/panel_wavelet_designer.m toolbox/gui/view_wavelet_designer.m toolbox/gui/figure_3d.m toolbox/tree/node_create_subject.m toolbox/tree/tree_callbacks.m dev/tests/test_wavelet_designer_session.m
git commit -m "refactor(gui): rename FilterDesigner -> WaveletDesigner (panel, view, tree, hook)"
```

---

## Task 5: Input becomes localization-only

**Files:**
- Modify: `toolbox/gui/panel_wavelet_designer.m`

- [ ] **Step 1: Remove the input radios from the Input section**

In `CreatePanel`, delete the two radio lines and their button group, replacing them with a one-line hint. The Section 1 header + input becomes:

```matlab
    jSec1 = gui_river([2 2], [2 8 3 6], '1. Input');
    gui_component('label', jSec1, 'br', '<HTML><I>Click a vertex on the cortex to place the wavelet.</I>');
```

Delete the lines that created `jInputDelta`, `jInputSource`, and `grpIn`.

- [ ] **Step 2: Drop the input handles from the ctrl struct and remove SetSeedSource**

In the `ctrl = struct(...)` literal, remove `'jInputDelta',jInputDelta, 'jInputSource',jInputSource,`. Delete the entire `SetSeedSource` subfunction (the filter-style source seeding path).

- [ ] **Step 3: Verify the panel still builds and seeds via click**

Run (MCP `evaluate_matlab_code`):

```matlab
rehash;
try, gui_hide('WaveletDesigner'); catch, end
hOld = findall(0,'Type','figure','Tag','3DViz');
for k=1:numel(hOld); try, bst_figures('DeleteFigure', hOld(k), []); catch, end; end
EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
hFig = view_wavelet_designer(EigenFile); drawnow;
panel_wavelet_designer('SetSeedVertex','WaveletDesigner',100); drawnow;
S = panel_wavelet_designer('GetState','WaveletDesigner');
fprintf('builds + seeds OK: %d (tiles=%d)\n', ~isempty(S.SeedCoeffs), numel(S.Tiles));
panel_wavelet_designer('OnCancel','WaveletDesigner'); drawnow;
```
Expected: `builds + seeds OK: 1 (tiles=1)`.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/panel_wavelet_designer.m
git commit -m "feat(wavelet): input is localization-only (drop filter-style source-map input)"
```

---

## Task 6: Seed direction in the local manifold frame

**Files:**
- Modify: `toolbox/gui/view_wavelet_designer.m` (find-or-create manifold, derive frame, pass to panel)
- Modify: `toolbox/gui/panel_wavelet_designer.m` (`i_embed_direction` helper, `SeedDirection`, slider labels/range)
- Test: `dev/tests/test_wavelet_direction.m`

- [ ] **Step 1: Write the failing embedding test**

Create `dev/tests/test_wavelet_direction.m`:

```matlab
function test_wavelet_direction()
% Pure test of the local-frame direction embedding d = cos(t)(cos(p)U+sin(p)V)+sin(t)N.
% Authors: Diellor Basha, 2026
    nFail = 0;
    U = [1 0 0]; V = [0 1 0]; N = [0 0 1];     % orthonormal frame
    % theta=+90 -> exactly the (outward) normal
    d = panel_wavelet_designer('EmbedDirection', 0, 90, U, V, N);
    nFail = nFail + chk('tilt +90 = normal', max(abs(d - N)) < 1e-12);
    % theta=0, phi=0 -> U
    d = panel_wavelet_designer('EmbedDirection', 0, 0, U, V, N);
    nFail = nFail + chk('tilt 0, phi 0 = U', max(abs(d - U)) < 1e-12);
    % theta=0, phi=90 -> V
    d = panel_wavelet_designer('EmbedDirection', 90, 0, U, V, N);
    nFail = nFail + chk('tilt 0, phi 90 = V', max(abs(d - V)) < 1e-12);
    % general: unit length, correct formula
    p = 37; t = 20;
    d = panel_wavelet_designer('EmbedDirection', p, t, U, V, N);
    expect = cosd(t)*(cosd(p)*U + sind(p)*V) + sind(t)*N; expect = expect/norm(expect);
    nFail = nFail + chk('general formula', max(abs(d - expect)) < 1e-12);
    nFail = nFail + chk('unit length', abs(norm(d)-1) < 1e-12);
    fprintf('\n==== test_wavelet_direction: %d failed ====\n', nFail);
    if nFail > 0, error('test_wavelet_direction FAILED'); end
end

function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run it to verify it fails**

Run: `dev/tests/test_wavelet_direction.m`
Expected: FAIL — `EmbedDirection` is not a dispatched subfunction yet.

- [ ] **Step 3: Add the pure embedding helper (dispatchable)**

In `panel_wavelet_designer.m`, add:

```matlab
%% ===== PURE: embed a local-frame direction into an ambient 3-vector =====
function d = EmbedDirection(phiDeg, thetaDeg, U, V, N) %#ok<DEFNU>
    d = cosd(thetaDeg) * (cosd(phiDeg)*U(:).' + sind(phiDeg)*V(:).') + sind(thetaDeg)*N(:).';
    nrm = norm(d); if nrm > 0; d = d / nrm; end
end
```

- [ ] **Step 4: Run the embedding test to verify it passes**

Run: `dev/tests/test_wavelet_direction.m`
Expected: PASS — all `PASS`, `0 failed`.

- [ ] **Step 5: Derive the manifold frame in the orchestrator and pass it to the panel**

In `view_wavelet_designer.m`, after resolving `EigenMat`/`SurfaceFile` and computing `nV`, add (before building the panel):

```matlab
    % --- local cortical frame (U,V,N) per vertex from the manifold node ---
    bst_progress('text', 'Loading cortical frame (manifold)...');
    ManifoldMat = tess_manifold(SurfaceFile);   % find-or-load-or-create
    Frame = view_manifold('DeriveVertexFrame', ManifoldMat.Embedded, ManifoldMat.Gauge, nV);
```

Change the panel construction to pass `Frame`:

```matlab
    bstPanel = panel_wavelet_designer('CreatePanel', EigenMat, EigenFile, hFig, ctxFn, Frame);
```

- [ ] **Step 6: Store the frame in panel state and rewrite SeedDirection**

In `panel_wavelet_designer.m` `CreatePanel`, add the `Frame` argument and store `Frame.U/V/N` in the session state:

```matlab
function bstPanelNew = CreatePanel(EigenMat, EigenFile, hFig, ctxFn, Frame) %#ok<DEFNU>
```
In the `S = struct(...)` initializer, add fields:
```matlab
               'FrameU', Frame.U, 'FrameV', Frame.V, 'FrameN', Frame.N, ...
```

Rewrite `SeedDirection` to use the local frame at the clicked vertex (φ = azimuth slider, θ = elevation slider):

```matlab
function d = SeedDirection(S, ctrl)
    d = [1 0 0];
    if ~S.isDirac || isempty(ctrl.jAz) || isempty(S.iVertex); return; end
    v = S.iVertex;
    if v > size(S.FrameU,1); return; end
    phi = double(ctrl.jAz.getValue());   % in-plane angle (deg)
    th  = double(ctrl.jEl.getValue());   % tilt toward normal (deg)
    d = EmbedDirection(phi, th, S.FrameU(v,:), S.FrameV(v,:), S.FrameN(v,:));
end
```

> Note: `SeedDirection` now needs the seed vertex. `SetSeedVertex` sets `S.iVertex` *before* calling `SeedDirection`; reorder so the vertex is stored first. In `SetSeedVertex`, set `S.iVertex = iVertex; SetState(ctrl,S);` then compute `d = SeedDirection(S, ctrl);` (re-fetch `S` if needed) before building the delta.

- [ ] **Step 7: Relabel the sliders + set the tilt default to outward normal**

In `CreatePanel`, change the two direction sliders so azimuth→in-plane and elevation→tilt (default 90):

```matlab
        gui_component('label', jSec1, 'br', '<HTML><I>Seed direction (local frame)</I>');
        [jAz, jAzVal] = i_labeled_slider(jSec1, '<HTML>In-plane angle: 0&deg;',  '0',  '360', 0, 360, 0);
        [jEl, jElVal] = i_labeled_slider(jSec1, '<HTML>Tilt to normal: 90&deg;', '-90','90', -90, 90, 90);
```

Update `OnDirSlider` label text to `'<HTML>In-plane angle: %d&deg;'` and `'<HTML>Tilt to normal: %d&deg;'`.

- [ ] **Step 8: Verify live — default seed is the outward normal**

Run (MCP `evaluate_matlab_code`):

```matlab
rehash;
try, gui_hide('WaveletDesigner'); catch, end
hOld = findall(0,'Type','figure','Tag','3DViz');
for k=1:numel(hOld); try, bst_figures('DeleteFigure', hOld(k), []); catch, end; end
EigenFile = bst_get('Subject',1).Surface(5).Eigen(1).FileName;
hFig = view_wavelet_designer(EigenFile); drawnow;
ctrl = bst_get('PanelControls','WaveletDesigner');
panel_wavelet_designer('SetSeedVertex','WaveletDesigner',100); drawnow;
S = panel_wavelet_designer('GetState','WaveletDesigner');
d = panel_wavelet_designer('SeedDirection', S, ctrl);
fprintf('default seed vs outward normal at v100: cos=%.4f (expect ~1)\n', dot(d, S.FrameN(100,:))/norm(S.FrameN(100,:)));
panel_wavelet_designer('OnCancel','WaveletDesigner'); drawnow;
```
Expected: `cos=~1.0` (default tilt 90° → the seed equals the outward normal at the clicked vertex).

- [ ] **Step 9: Run the direction test + the session test**

Run `dev/tests/test_wavelet_direction.m` and `dev/tests/test_wavelet_designer_session.m`.
Expected: both PASS — `0 failed`.

- [ ] **Step 10: Commit**

```bash
git add toolbox/gui/view_wavelet_designer.m toolbox/gui/panel_wavelet_designer.m dev/tests/test_wavelet_direction.m
git commit -m "feat(wavelet): seed direction in the local manifold frame (in-plane + tilt)"
```

---

## Self-Review

**Spec coverage:**
- Rename tool + DB artifact filter→wavelet → Tasks 1, 2, 3, 4. ✓
- Input localization-only (drop source-map input) → Task 5. ✓
- Seed direction in local manifold frame (two angles, default outward normal) → Task 6. ✓
- Keep kernel/scale/chirality/tiling/save unchanged → untouched by all tasks (only `SeedDirection` + input change in the panel). ✓
- New manifold dependency (find-or-create on open) → Task 6 Step 5. ✓
- DB v5.06 migration → Task 1. ✓
- bst-java node type → Task 3. ✓
- Tests renamed + direction/embedding test added → Tasks 1,2,4,6. ✓
- Gauge caveat (φ gauge-dependent) → documented in the panel direction label/help (Task 6 Step 7 uses "local frame"); the manifold's default gauge is used by `tess_manifold(SurfaceFile)`. ✓

**Placeholder scan:** No "TBD/TODO". Rename steps list exact old→new identifiers per file (a rename is a complete substitution). Behavioral steps (5,6) show full code. The bst-java task is conditional but explicit about both branches (found / not found). ✓

**Type consistency:** `waveletmat`/`Surface.Wavelet`/`WaveletFile`/`db_add_wavelet`/`wavelet` file type used consistently across Tasks 1–4 and the tests. `EmbedDirection(phiDeg, thetaDeg, U, V, N)` defined in Task 6 Step 3, used in Task 6 Steps 6, 8 and the test. `S.FrameU/FrameV/FrameN` set in CreatePanel (Step 6) and read in `SeedDirection` + the live check. `bst_filterbank_tiles(wavelet, opts)` signature unchanged and used identically in `test_db_add_wavelet`. The saved `Tiling = struct(Wavelet, Opts)` shape matches `test_db_add_wavelet`'s assertions. ✓

---

## Build order summary

1. Task 1 — DB schema + v5.06 migration (waveletmat, Surface.Wavelet)
2. Task 2 — file type + WaveletFile accessor + db_add_wavelet
3. Task 3 — bst-java wavelet node type (before the tree renders wavelet nodes)
4. Task 4 — GUI + tree rename (FilterDesigner → WaveletDesigner)
5. Task 5 — input localization-only
6. Task 6 — seed direction in the local manifold frame

Tasks 1–4 are the behavior-preserving rename (verified by the renamed schema/db_add/session tests). Tasks 5–6 are the two behavioral corrections (verified by the build/seed check, the pure embedding test, and the outward-normal live check).
