# Source-Vector Quiver Overlay — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a quiver (arrow) overlay to Brainstorm's existing "Display on cortex" 3D source view that shows the ambient 3D direction of an unconstrained source field and updates as the user time-steps.

**Architecture:** Two pure, headless-testable subfunctions (`ComputeSourceVectorGlyphs`, `SelectSourceVectorIdx`) compute glyph geometry; a find-or-create `PlotSourceVectors(hFig,iTess)` draws/updates a tagged `quiver3` object on the 3D axes, mirroring the existing `PlotGrid`. It is called from `UpdateSurfaceColor` (which already runs on every time-cursor change), so time-stepping is automatic. A `CheckBoxMenuItem` in `DisplayFigurePopup` toggles it via `SetShowSourceVectors`. No inverse changes.

**Tech Stack:** MATLAB R2023b, Brainstorm (fork at `/Users/diellorbasha/workspace/research/code/brainstorm3`), Java/Swing GUI via `gui_component`, MATLAB MCP for headless test runs. Brainstorm must be running (`brainstorm('status')==1`) for dispatch/tests.

**Conventions:** Subfunctions are dispatched via `eval(macro_method)` — callable as `figure_3d('FuncName',args)`. New display subfunctions get a `%#ok<DEFNU>` tag. Vertices/positions are in meters (SCS). Arrows are unit-normalized direction only; amplitude stays on the colormap. Display is the ambient pullback-bundle field `f*Tℝ³` — the Dirac operator is NOT used here.

**Spec:** `docs/superpowers/specs/2026-06-11-dirac-source-vector-quiver-viewer-design.md`

---

## File Structure

| File | Change | Responsibility |
|---|---|---|
| `toolbox/gui/figure_3d.m` | Modify | Add `ComputeSourceVectorGlyphs`, `SelectSourceVectorIdx`, `PlotSourceVectors`, `SetShowSourceVectors` subfunctions; one call in `UpdateSurfaceColor`; one `CheckBoxMenuItem` in `DisplayFigurePopup`. |
| `dev/tests/test_source_vector_glyphs.m` | Create | Headless regression test for the two pure glyph subfunctions. |
| `dev/tests/make_alpha_dirac_source.m` | Create | One-off script: create an unconstrained Dirac source results node for the alpha block (validation fixture) if none exists. |

All four new subfunctions live in `figure_3d.m` because they share the file's `eval(macro_method)` dispatch and operate on its figures; this matches how `PlotGrid`/`PlotSensors3D` are organized (one focused subfunction each, same file).

---

## Task 1: Pure glyph geometry — `ComputeSourceVectorGlyphs`

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (add subfunction near `PlotGrid`, ~line 2963)
- Test: `dev/tests/test_source_vector_glyphs.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_source_vector_glyphs.m`:

```matlab
function test_source_vector_glyphs()
% TEST_SOURCE_VECTOR_GLYPHS  Headless regression for the source-vector quiver helpers.
% Requires Brainstorm running (figure_3d on path). Tests pure subfunctions via dispatch.
% Authors: Diellor Basha, 2026
    nPass = 0; nFail = 0;

    % ===== ComputeSourceVectorGlyphs =====
    P   = [0 0 0; 1 0 0; 2 0 0];           % anchors
    Nrm = [0 0 1; 0 0 1; 0 0 1];           % +z normals
    V3  = [3 0 0;  0 0 0;  0 4 0];         % vtx1 +x(mag3), vtx2 zero, vtx3 +y(mag4)
    G = figure_3d('ComputeSourceVectorGlyphs', P, Nrm, V3, [], 0.5, 0.1);
    [nPass,nFail] = chk('glyph offset lifts base along +z', all(abs(G.Z - 0.1) < 1e-12), nPass,nFail);
    [nPass,nFail] = chk('glyph unit-normalized +x * scale', abs(G.U(1)-0.5)<1e-12 && abs(G.V(1))<1e-12 && abs(G.W(1))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('glyph zero vector -> zero arrow',   abs(G.U(2))<1e-12 && abs(G.V(2))<1e-12 && abs(G.W(2))<1e-12, nPass,nFail);
    [nPass,nFail] = chk('glyph unit-normalized +y * scale',  abs(G.V(3)-0.5)<1e-12 && abs(G.U(3))<1e-12, nPass,nFail);
    G2 = figure_3d('ComputeSourceVectorGlyphs', P, Nrm, V3, [1 3], 1, 0);
    [nPass,nFail] = chk('glyph idx subselect count==2', numel(G2.X)==2, nPass,nFail);

    fprintf('\n==== test_source_vector_glyphs: %d passed, %d failed ====\n', nPass, nFail);
    if nFail > 0, error('test_source_vector_glyphs: %d test(s) FAILED.', nFail); end
end

function [p,f] = chk(name, c, p, f)
    if c, fprintf('  PASS  %s\n', name); p=p+1; else, fprintf('  FAIL  %s\n', name); f=f+1; end
end
```

- [ ] **Step 2: Run the test to verify it fails**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_source_vector_glyphs.m`
Expected: FAIL — `Undefined function or variable 'ComputeSourceVectorGlyphs'` (dispatch error from `figure_3d`).

- [ ] **Step 3: Add the subfunction**

In `toolbox/gui/figure_3d.m`, immediately after the end of `PlotGrid` (the `end` at line 2963), insert:

```matlab
%% ===== COMPUTE SOURCE VECTOR GLYPHS =====
% Build quiver3 glyph geometry for an ambient (Cartesian) source vector field.
%   P         [nVert x 3] anchor positions (displayed surface vertices, meters)
%   Nrm       [nVert x 3] vertex normals (for lifting arrows off the surface)
%   V3        [nVert x 3] per-vertex source 3-vector (ambient SCS components)
%   idx       [k x 1] vertex indices to draw ([] = all)
%   ScaleLen  scalar arrow length in meters (unit direction is scaled by this)
%   OffsetLen scalar lift in meters along the normal (avoids depth-buffer occlusion)
% Returns struct G with column vectors X,Y,Z (arrow bases) and U,V,W (arrow vectors).
% Direction is unit-normalized with an eps-guard; below-eps vectors get zero arrows.
function G = ComputeSourceVectorGlyphs(P, Nrm, V3, idx, ScaleLen, OffsetLen) %#ok<DEFNU>
    if isempty(idx)
        idx = (1:size(P,1))';
    end
    idx = idx(:);
    Pi = P(idx,:);   Ni = Nrm(idx,:);   Vi = V3(idx,:);
    % Unit normals (guard against zero-length)
    nN = sqrt(sum(Ni.^2,2));   nN(nN < eps) = 1;
    Ni = Ni ./ nN;
    % Arrow bases lifted along the normal
    Base = Pi + OffsetLen * Ni;
    % Unit-normalized directions (eps-guard -> zero arrow)
    mag = sqrt(sum(Vi.^2,2));
    good = mag > eps;
    Dir = zeros(size(Vi));
    Dir(good,:) = Vi(good,:) ./ mag(good);
    Vec = ScaleLen * Dir;
    G = struct('X', Base(:,1), 'Y', Base(:,2), 'Z', Base(:,3), ...
               'U', Vec(:,1),  'V', Vec(:,2),  'W', Vec(:,3));
end
```

- [ ] **Step 4: Run the test to verify it passes (the `SelectSourceVectorIdx` cases are added in Task 2)**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_source_vector_glyphs.m`
Expected: PASS — `5 passed, 0 failed`.

- [ ] **Step 5: Lint the modified file**

Run (MATLAB MCP `check_matlab_code`): `toolbox/gui/figure_3d.m`
Expected: no new errors introduced by the added subfunction (pre-existing Brainstorm warnings may remain).

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/figure_3d.m dev/tests/test_source_vector_glyphs.m
git commit -m "feat(gui): add ComputeSourceVectorGlyphs helper for source-vector quiver"
```

---

## Task 2: Decimation set — `SelectSourceVectorIdx`

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (add subfunction after `ComputeSourceVectorGlyphs`)
- Test: `dev/tests/test_source_vector_glyphs.m` (extend)

- [ ] **Step 1: Extend the test with failing cases**

In `dev/tests/test_source_vector_glyphs.m`, insert these lines just before the `fprintf('\n==== ...` summary line:

```matlab
    % ===== SelectSourceVectorIdx =====
    mag = (1:100)';
    [nPass,nFail] = chk('select all when MaxArrows empty', numel(figure_3d('SelectSourceVectorIdx', mag, []))==100, nPass,nFail);
    idxDec = figure_3d('SelectSourceVectorIdx', mag, 10);
    [nPass,nFail] = chk('select decimates to <=MaxArrows', numel(idxDec) <= 10, nPass,nFail);
    [nPass,nFail] = chk('select decimation is a stable step', isequal(idxDec, (1:ceil(100/10):100)'), nPass,nFail);
    [nPass,nFail] = chk('select all when MaxArrows>=nVert', numel(figure_3d('SelectSourceVectorIdx', mag, 1000))==100, nPass,nFail);
```

- [ ] **Step 2: Run the test to verify the new cases fail**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_source_vector_glyphs.m`
Expected: FAIL — `Undefined function or variable 'SelectSourceVectorIdx'`.

- [ ] **Step 3: Add the subfunction**

In `toolbox/gui/figure_3d.m`, immediately after the end of `ComputeSourceVectorGlyphs`, insert:

```matlab
%% ===== SELECT SOURCE VECTOR INDICES =====
% Choose which vertices get an arrow. Default: all vertices. With MaxArrows set,
% decimate by a stable step so arrows keep fixed anchors across time frames.
%   mag       [nVert x 1] per-vertex magnitude (reserved for future threshold gating)
%   MaxArrows scalar cap ([] / non-finite / >= nVert -> all vertices)
function idx = SelectSourceVectorIdx(mag, MaxArrows) %#ok<DEFNU>
    nVert = numel(mag);
    if isempty(MaxArrows) || ~isfinite(MaxArrows) || (MaxArrows >= nVert)
        idx = (1:nVert)';
    else
        step = ceil(nVert / max(MaxArrows,1));
        idx  = (1:step:nVert)';
    end
end
```

- [ ] **Step 4: Run the test to verify all pass**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_source_vector_glyphs.m`
Expected: PASS — `9 passed, 0 failed`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/gui/figure_3d.m dev/tests/test_source_vector_glyphs.m
git commit -m "feat(gui): add SelectSourceVectorIdx decimation helper"
```

---

## Task 3: Draw/update overlay — `PlotSourceVectors`

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (add subfunction after `SelectSourceVectorIdx`)

No headless unit test (needs a live figure); validated in Task 6. This task adds the function and lints it.

- [ ] **Step 1: Add the subfunction**

In `toolbox/gui/figure_3d.m`, immediately after the end of `SelectSourceVectorIdx`, insert:

```matlab
%% ===== PLOT SOURCE VECTORS =====
% Find-or-create the ambient source-vector quiver overlay for surface iTess.
% Reads the un-oriented 3-vector at the current time and draws unit-normalized
% direction arrows over the amplitude colormap. Mirrors PlotGrid's find-or-create.
function hQuiver = PlotSourceVectors(hFig, iTess) %#ok<DEFNU>
    QuiverTag = 'SourceVectors';
    hAxes   = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    hQuiver = findobj(hAxes, 'Tag', QuiverTag);
    % Get surfaces
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (iTess > numel(TessInfo))
        return;
    end
    sTess = TessInfo(iTess);
    % Off / not a source surface -> remove any existing quiver and stop
    showVec = isfield(sTess,'ShowSourceVectors') && ~isempty(sTess.ShowSourceVectors) && sTess.ShowSourceVectors;
    isSrc   = ~isempty(sTess.DataSource) && strcmpi(sTess.DataSource.Type, 'Source') && ~isempty(sTess.DataSource.FileName);
    if ~showVec || ~isSrc
        if ~isempty(hQuiver), delete(hQuiver); end
        hQuiver = [];
        return;
    end
    % Resolve dataset / result / un-oriented 3-vector at the current time
    [tmp__, iFig, iDS] = bst_figures('GetFigure', hFig); %#ok<ASGLU>
    if isempty(iDS), if ~isempty(hQuiver), delete(hQuiver); end, hQuiver = []; return; end
    iResult = bst_memory('GetResultInDataSet', iDS, sTess.DataSource.FileName);
    if isempty(iResult), if ~isempty(hQuiver), delete(hQuiver); end, hQuiver = []; return; end
    [V3col, nComponents] = bst_memory('GetResultsValues', iDS, iResult, [], 'CurrentTimeIndex', 0);
    if (nComponents ~= 3)   % only unconstrained fields carry ambient vectors
        if ~isempty(hQuiver), delete(hQuiver); end
        hQuiver = [];
        return;
    end
    V3 = reshape(V3col, 3, [])';            % [nVert x 3] ambient components
    % Anchors + normals from the DISPLAYED patch (matches inflation/smoothing)
    hPatch = sTess.hPatch;
    P   = get(hPatch, 'Vertices');
    Nrm = get(hPatch, 'VertexNormals');
    if isempty(Nrm) || (size(Nrm,1) ~= size(P,1))
        Nrm = repmat([0 0 1], size(P,1), 1);
    end
    % Decimation ("threshold the quiver number"); default = entire field
    if isfield(sTess,'SourceVectorMaxArrows'), MaxArrows = sTess.SourceVectorMaxArrows; else, MaxArrows = []; end
    mag = sqrt(sum(V3.^2,2));
    idx = SelectSourceVectorIdx(mag, MaxArrows);
    % Arrow length + surface lift (meters)
    if isfield(sTess,'SourceVectorScale') && ~isempty(sTess.SourceVectorScale)
        ScaleLen = sTess.SourceVectorScale;
    else
        ScaleLen = 0.004;
    end
    OffsetLen = 0.0015;
    G = ComputeSourceVectorGlyphs(P, Nrm, V3, idx, ScaleLen, OffsetLen);
    % Create or update the quiver object
    if isempty(hQuiver)
        hQuiver = quiver3(G.X, G.Y, G.Z, G.U, G.V, G.W, 0, ...
            'Parent', hAxes, 'Color', [0 0 0], 'LineWidth', 1, ...
            'MaxHeadSize', 0.5, 'AutoScale', 'off', 'Tag', QuiverTag);
    else
        set(hQuiver, 'XData', G.X, 'YData', G.Y, 'ZData', G.Z, ...
                     'UData', G.U, 'VData', G.V, 'WData', G.W);
    end
end
```

- [ ] **Step 2: Lint the modified file**

Run (MATLAB MCP `check_matlab_code`): `toolbox/gui/figure_3d.m`
Expected: no new errors from `PlotSourceVectors`.

- [ ] **Step 3: Smoke-check it dispatches and no-ops without a source figure**

Run (MATLAB MCP `evaluate_matlab_code`):

```matlab
hF = figure('Tag','tmpNoAxes');
try
    r = figure_3d('PlotSourceVectors', hF, 1);   % no 'Surface' appdata -> early return
    fprintf('PlotSourceVectors no-op OK, returned isempty=%d\n', isempty(r));
catch e
    fprintf('ERROR: %s\n', e.message);
end
close(hF);
```
Expected: `PlotSourceVectors no-op OK, returned isempty=1` (no error).

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "feat(gui): add PlotSourceVectors quiver overlay (find-or-create)"
```

---

## Task 4: Auto-update on time step — hook into `UpdateSurfaceColor`

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (`UpdateSurfaceColor`, before its closing `end` at ~line 2881)

- [ ] **Step 1: Add the overlay-update call**

In `toolbox/gui/figure_3d.m`, in `UpdateSurfaceColor(hFig, iTess)`, immediately BEFORE the function's closing `end` (the `end` at line 2881, after the existing color-update logic), insert:

```matlab
    % Update the ambient source-vector quiver overlay (no-op unless enabled)
    PlotSourceVectors(hFig, iTess);
```

(`UpdateSurfaceColor` already returns early when `hPatch` is invalid, so this line only runs when there is a real patch. It is reached on every time-cursor change because `UpdateSurfaceColormap` calls `UpdateSurfaceColor` at `panel_surface.m:2241`.)

- [ ] **Step 2: Lint the modified file**

Run (MATLAB MCP `check_matlab_code`): `toolbox/gui/figure_3d.m`
Expected: no new errors.

- [ ] **Step 3: Verify existing color display is unbroken (regression)**

Run (MATLAB MCP `run_matlab_file`): `dev/tests/test_source_vector_glyphs.m`
Expected: PASS — `9 passed, 0 failed` (helpers unaffected). Full GUI redraw is checked in Task 6.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "feat(gui): refresh source-vector overlay on every surface-color update"
```

---

## Task 5: Toggle — `SetShowSourceVectors` + popup menu item

**Files:**
- Modify: `toolbox/gui/figure_3d.m` (add `SetShowSourceVectors` after `PlotSourceVectors`; add menu item in `DisplayFigurePopup`)

- [ ] **Step 1: Add the setter subfunction**

In `toolbox/gui/figure_3d.m`, immediately after the end of `PlotSourceVectors`, insert:

```matlab
%% ===== SET SHOW SOURCE VECTORS =====
% Toggle the ambient source-vector quiver overlay for surface iTess.
function SetShowSourceVectors(hFig, iTess, isShow) %#ok<DEFNU>
    TessInfo = getappdata(hFig, 'Surface');
    if isempty(TessInfo) || (iTess > numel(TessInfo))
        return;
    end
    TessInfo(iTess).ShowSourceVectors = isShow;
    setappdata(hFig, 'Surface', TessInfo);
    if isShow
        PlotSourceVectors(hFig, iTess);
    else
        hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
        delete(findobj(hAxes, 'Tag', 'SourceVectors'));
    end
end
```

- [ ] **Step 2: Add the checkbox menu item**

In `toolbox/gui/figure_3d.m`, in `DisplayFigurePopup`, locate the final display call near the end of the function: `gui_popup(jPopup, hFig);`. Immediately BEFORE that line, insert:

```matlab
    % ==== Source vectors (ambient quiver overlay) ====
    iSrcTess = find(arrayfun(@(t) ~isempty(t.DataSource) && strcmpi(t.DataSource.Type,'Source') && ~isempty(t.DataSource.FileName), TessInfo), 1);
    if ~isempty(iSrcTess)
        isShowVec = isfield(TessInfo(iSrcTess),'ShowSourceVectors') && ~isempty(TessInfo(iSrcTess).ShowSourceVectors) && TessInfo(iSrcTess).ShowSourceVectors;
        jItem = gui_component('CheckBoxMenuItem', jPopup, [], 'Show source vectors (quiver)', [], [], @(h,ev)SetShowSourceVectors(hFig, iSrcTess, ~isShowVec));
        jItem.setSelected(isShowVec);
    end
```

(`TessInfo` is already fetched at the top of `DisplayFigurePopup` via `getappdata(hFig,'Surface')`. The icon arg is `[]` to avoid depending on a specific `IconLoader` constant.)

- [ ] **Step 3: Lint the modified file**

Run (MATLAB MCP `check_matlab_code`): `toolbox/gui/figure_3d.m`
Expected: no new errors.

- [ ] **Step 4: Commit**

```bash
git add toolbox/gui/figure_3d.m
git commit -m "feat(gui): add 'Show source vectors (quiver)' toggle to 3D figure popup"
```

---

## Task 6: Validation fixture + manual end-to-end check

**Files:**
- Create: `dev/tests/make_alpha_dirac_source.m`

This creates an unconstrained (3-component) Dirac source results node for the alpha block so the overlay has a real field to draw, then guides a manual GUI check. The kernel comes from the existing `bst_inverse_dirac` (no inverse change); the results struct is templated from an existing results file in the study to guarantee field-completeness.

- [ ] **Step 1: Write the fixture script**

Create `dev/tests/make_alpha_dirac_source.m`:

```matlab
function ResultsFile = make_alpha_dirac_source()
% MAKE_ALPHA_DIRAC_SOURCE  Create an unconstrained Dirac source node for the
% alpha block (data_block001_band) so the quiver overlay has a vector field.
% Reuses an existing results_*.mat in the same study as a struct template and
% swaps in the Dirac kernel from bst_inverse_dirac. Returns the new file path.
% Authors: Diellor Basha, 2026
    cond = 'Subject01/S01_AEF_20131218_01_notch';
    sStudy = bst_get('StudyWithCondition', cond);
    % Locate the band-passed alpha data node
    iData = find(~cellfun(@isempty, regexp({sStudy.Data.FileName}, 'data_block001_band', 'once')), 1);
    if isempty(iData), error('alpha band data node not found in %s', cond); end
    DataFile = sStudy.Data(iData).FileName;
    % Inputs for the inverse
    HMfile  = sStudy.HeadModel(sStudy.iHeadModel).FileName;
    NCfile  = sStudy.NoiseCov(1).FileName;
    HMos    = in_bst_headmodel(HMfile, 0);
    Chan    = in_bst_channel(sStudy.Channel.FileName);
    G       = double(HMos.Gain);
    iMEG    = all(isfinite(G),2) & strcmpi({Chan.Channel.Type}', 'MEG');
    HMf     = HMos; HMf.Gain = G(iMEG,:);
    NC      = load(file_fullpath(NCfile));
    OPT = struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3, ...
                 'InverseMeasure','amplitude','nModes',400,'Tau',0.5);
    OPT.NoiseCovMat.NoiseCov = NC.NoiseCov(iMEG,iMEG);
    OPT.ChannelTypes = {Chan.Channel(iMEG).Type};
    R = bst_inverse_dirac(HMf, OPT);                 % R.ImagingKernel [3*nVert x nMEG]
    % Template from an existing source results file in the same study
    iSrc = find(arrayfun(@(r) ~isempty(strfind(lower(r.FileName),'results')), sStudy.Result), 1);
    if isempty(iSrc), error('no existing results_*.mat to template from in %s', cond); end
    Tmpl = load(file_fullpath(sStudy.Result(iSrc).FileName));
    % Build the unconstrained kernel results (nComponents=3)
    nMEG = sum(iMEG);
    Tmpl.ImagingKernel = R.ImagingKernel;
    Tmpl.ImageGridAmp  = [];
    Tmpl.nComponents   = 3;
    Tmpl.GoodChannel   = find(iMEG(:)');
    Tmpl.DataFile      = DataFile;
    Tmpl.HeadModelFile = HMfile;
    Tmpl.SurfaceFile   = HMos.SurfaceFile;
    Tmpl.Function      = 'dirac';
    Tmpl.Comment       = 'Dirac vector (amplitude) | alpha band';
    Tmpl.Time          = [];
    % Save + register
    ResultsFile = bst_process('GetNewFilename', fileparts(file_fullpath(DataFile)), 'results_dirac_vec');
    bst_save(ResultsFile, Tmpl, 'v6');
    db_add_data(sStudy.iStudy, ResultsFile, Tmpl);
    db_reload_studies(sStudy.iStudy);
    ResultsFile = file_short(ResultsFile);
    fprintf('Created unconstrained Dirac source node: %s\n', ResultsFile);
end
```

- [ ] **Step 2: Run the fixture and confirm it loads as unconstrained**

Run (MATLAB MCP `evaluate_matlab_code`):

```matlab
ResultsFile = make_alpha_dirac_source();
[iDS, iRes] = bst_memory('LoadResultsFile', ResultsFile, 0);
[V3, nComp] = bst_memory('GetResultsValues', iDS, iRes, [], 'CurrentTimeIndex', 0);
fprintf('nComponents=%d, V3 length=%d (expect 3*nVert)\n', nComp, numel(V3));
```
Expected: `nComponents=3` and `V3 length` = 3×(#cortex vertices).

- [ ] **Step 3: Manual GUI validation (record results in the commit message)**

Do this in the running Brainstorm GUI:
1. In the database tree, right-click the new "Dirac vector (amplitude) | alpha band" result → **Cortical activations → Display on cortex**. Confirm the colormap source map appears.
2. Right-click the 3D figure → confirm **"Show source vectors (quiver)"** is present; click it. Expected: black unit arrows appear over the cortex, lifted slightly off the surface (no z-fighting), with amplitude still read from the colormap.
3. Drag the time cursor across the burst (≈21.5–24 s in the 20–25 s block, peak ≈22.6 s). Expected: arrows update at each time step and the field direction evolves cycle-by-cycle over the right parieto-occipital generator.
4. Toggle the menu item off → arrows disappear; on → reappear. Confirm the existing transparency/threshold sliders still behave.
5. Note interactivity at full density (~20k arrows). If time-stepping is sluggish, set a cap once and re-check:
   ```matlab
   % find the source figure's TessInfo index, then:
   TessInfo = getappdata(hFig,'Surface');
   TessInfo(iSrcTess).SourceVectorMaxArrows = 800;
   setappdata(hFig,'Surface',TessInfo);
   panel_surface('UpdateSurfaceColormap', hFig);
   ```
   Expected: ≤800 arrows, smooth stepping. If a cap is needed, record the chosen default; a follow-up can wire it into the popup.

- [ ] **Step 4: Commit the fixture and a note of the validation outcome**

```bash
git add dev/tests/make_alpha_dirac_source.m
git commit -m "test(gui): alpha Dirac vector source fixture + manual quiver validation

Validated: overlay toggles, draws ambient arrows lifted off cortex, and
updates per time step across the alpha burst. [record density/perf note here]"
```

---

## Self-Review

**1. Spec coverage:**
- Quiver overlay in the existing Display-on-cortex figure → Tasks 3–5 (`PlotSourceVectors` on the same `Axes3D`, popup toggle). ✓
- Unit-normalized direction arrows, amplitude on colormap → Task 1 (`ComputeSourceVectorGlyphs` normalizes; colormap untouched). ✓
- ε-guard for near-zero vectors → Task 1 (`good = mag > eps`). ✓
- Slight normal offset, opaque cortex by default → Task 1/3 (`OffsetLen=0.0015`; no transparency change). ✓
- Entire field by default + threshold/decimate control → Task 2 (`SelectSourceVectorIdx` default all; `SourceVectorMaxArrows`). ✓
- Ambient Cartesian components straight to quiver → Task 1/3 (no projection; `reshape(V3col,3,[])'` → `U,V,W`). ✓
- Time-stepping via existing redraw hook → Task 4 (call in `UpdateSurfaceColor`). ✓
- Constrained/scalar source disables overlay → Task 3 (`nComponents ~= 3` deletes/no-ops). ✓
- Kernel and full-matrix results both supported → Task 3 (uses `GetResultsValues`, which handles both). ✓
- No inverse changes → confirmed; only `figure_3d.m` + test/fixture files touched. ✓
- Pullback-bundle / Dirac-out-of-scope framing → respected (renderer is pure ambient). ✓

**2. Placeholder scan:** No "TBD/TODO/handle edge cases". The single bracketed `[record ...]` in the Task 6 commit message is an instruction to paste a real measured result, not code. Acceptable.

**3. Type consistency:** `ComputeSourceVectorGlyphs` returns struct `G` with `X,Y,Z,U,V,W` — consumed identically in `PlotSourceVectors` (`quiver3(G.X,...,G.W)` and the `set(...)` update). `SelectSourceVectorIdx(mag, MaxArrows)` — same signature in test (Task 2) and caller (Task 3). `SetShowSourceVectors(hFig, iTess, isShow)` — same in setter (Task 5) and menu callback (Task 5). Tag `'SourceVectors'` is identical in `PlotSourceVectors` (create/find) and `SetShowSourceVectors` (delete). Consistent. ✓
