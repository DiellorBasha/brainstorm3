# Atom Geodesic Regions + Phase Filtering Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Let the user attach a geodesic cortical region to the active atom (localizing a time-only phase marker) and filter the atom list/cortex by oscillation phase.

**Architecture:** Three layers over the existing Dynamics stack. (1) Data model: one new per-occurrence `region` cell field on `atomgroup` plus a `bst_dynamics('AttachRegion', …)` verb that pads the occurrence arrays and writes a seed+region. (2) Capture UI: the user grows a disk with the existing Scout Area tool, then a Dynamics "Capture region → active atom" action snapshots the selected scout's vertices into the selected list occurrence; `view_dynamics('Redraw')` renders captured regions as translucent patches and tolerates partial localization. (3) Phase filter: a "Show phases" checkbox submenu drives `st.showPhase`, filtering both the flat window list and the cortex markers/regions.

**Tech Stack:** MATLAB, Brainstorm GUI (`gui_component`, Java Swing via `BstPanel`), `tess_scout_area` (heat-distance geodesic engine), `panel_scout` (Scout Area tool + `GetSelectedScouts`).

## Global Constraints

- The new `region` field is added to the `db_template('atomgroup')` template itself, so old tables normalize forward (Load uses `struct_copy_fields(template, …)`) and the regression `isequal(fieldnames(T.Groups), fieldnames(db_template('atomgroup')))` keeps holding.
- `region{i}` is a snapshot **copy** of vertex indices; atoms never reference a live Scout (avoids the bulk-scout-delete atlas hazard).
- Phase filtering is **non-destructive** display state — it never deletes atoms.
- The phase→index map is fixed: `peak=1, trough=2, rising=3, falling=4`; groups with an empty `phase` are always shown.
- Per-occurrence arrays stay aligned to `N = size(G.times,2)`; unlocalized occurrences hold `NaN` seeds / `{}` regions, and `line()`/Redraw skip them (no marker, no highlight).
- Do not start implementation on `development`; this work is on branch `feature/atom-regions-phasefilter` (already created, spec committed there).
- Tests run headless in MATLAB with Brainstorm running, via `test_dynamics_atoms`. Tests that need a 3D figure reuse the existing unconstrained-kernel fixture and SKIP together when it is absent.

---

### Task 1: Data model — `region` field + `AttachRegion` verb

**Files:**
- Modify: `toolbox/db/db_template.m` (the `case 'atomgroup'` struct, after the `'hemi'` field ~line 350)
- Modify: `toolbox/dynamics/bst_dynamics.m` (add an `AttachRegion` verb + 3 padding helpers)
- Test: `dev/test_dynamics_atoms.m` (extend T1's assertions; add an always-run `AttachRegion` check)

**Interfaces:**
- Consumes: `db_template('atomgroup')` (existing per-occurrence fields `vertices [1×N]`, `pos [N×3]`, `hemi [1×N]`, `strength [1×N]`, `charge [1×N]`, `times [1×N|2×N]`).
- Produces:
  - `atomgroup.region` — cell `{1×N}`, `region{i} = [v1 … vk]` cortex vertex indices; `{}` ⇒ point atom.
  - `G = bst_dynamics('AttachRegion', G, o, regionVerts, seed, pos, hemi)` — pads `vertices/pos/hemi/strength/charge` to length `N=size(G.times,2)` (NaN fill) and `region` to `N` (`[]` fill), then sets occurrence `o`: `region{o}=regionVerts(:)'`, `vertices(o)=seed`, `pos(o,:)=pos`, `hemi(o)=hemi`; sets `G.type='simple'`. Returns the updated group.

- [ ] **Step 1: Write the failing test**

In `dev/test_dynamics_atoms.m`, immediately after the T1 block (after the line `pass = pass && ok1;` near line 38) and **before** the `[linkFile, relData] = i_find_kernel();` line, insert an always-run model test for the new field + verb:

```matlab
    % ---------- T1b: region field default + AttachRegion padding/localization ----------
    hasRegion = isfield(db_template('atomgroup'), 'region');
    % a time-only marker: 3 occurrences, no space yet
    M = bst_dynamics('NewGroup', 'mk');
    M.type='simple';  M.times=[0.10 0.20 0.30];  M.phase='peak';
    M2 = bst_dynamics('AttachRegion', M, 2, [11 12 13 14], 12, [0.01 0.02 0.03], 1);
    padOK = isequal(size(M2.pos),[3 3]) && (numel(M2.vertices)==3) && (numel(M2.region)==3) ...
         && isequal(M2.region{2}, [11 12 13 14]) && (M2.vertices(2)==12) ...
         && isequal(M2.pos(2,:), [0.01 0.02 0.03]) && (M2.hemi(2)==1) ...
         && isnan(M2.vertices(1)) && isempty(M2.region{1}) && all(isnan(M2.pos(1,:)));
    % round-trip preserves the region cell + schema
    Tr = bst_dynamics('New','rt');  Tr = bst_dynamics('AddGroup', Tr, M2);
    fr = fullfile(bst_get('BrainstormTmpDir'), 'dyn_region_unit.mat');
    bst_dynamics('Save', fr, Tr);  Tr2 = bst_dynamics('Load', fr);
    rtOK = isequal(Tr2.Groups(1).region{2}, [11 12 13 14]) ...
        && isequal(fieldnames(Tr2.Groups), fieldnames(db_template('atomgroup')));
    ok1b = hasRegion && padOK && rtOK;
    fprintf('T1b region/AttachRegion: field=%d pad=%d roundtrip=%d => %s\n', hasRegion, padOK, rtOK, PF{ok1b+1});
    pass = pass && ok1b;
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB, Brainstorm running):
```matlab
test_dynamics_atoms
```
Expected: FAIL — `T1b region/AttachRegion: field=0 ...` and an error or `=> FAIL` because `region` is not a template field and `AttachRegion` is an undefined verb (`Unknown command 'AttachRegion'`).

- [ ] **Step 3: Add the `region` field to the template**

In `toolbox/db/db_template.m`, in the `case 'atomgroup'` struct, add the `region` field right after the `'hemi'` line. Wrap the empty cell as `{{}}` so `struct(...)` stores a `{}` field value (the standard empty-cell-in-struct idiom):

```matlab
            'hemi',        [], ...      % SPACE: [1 x N] 1=L 2=R
            'region',      {{}}, ...    % SPACE: cell {1 x N}, region{i}=[v...] geodesic-disk vertices for occurrence i ({}=point atom)
```

- [ ] **Step 4: Add the `AttachRegion` verb + padding helpers**

In `toolbox/dynamics/bst_dynamics.m`, add this function (and its 3 helpers) after the `Extrema` block (after `i_local_ext`, ~line 144):

```matlab
%% ===== ATTACH REGION (localize occurrence o with a geodesic region + seed) =====
% Pads the per-occurrence arrays to full length N = size(G.times,2), then writes occurrence o's
% seed (vertices/pos/hemi) and region. A time-only marker (empty vertices/pos) becomes partially
% localized: occurrence o gains a finite seed + region; the other occurrences stay NaN/[] (time only).
%   o           occurrence column to localize (1..N)
%   regionVerts [1 x k] cortex vertex indices of the geodesic disk (snapshot copy)
%   seed        scalar seed vertex index
%   pos         [1 x 3] seed position (SCS)
%   hemi        scalar 1=L 2=R
function G = AttachRegion(G, o, regionVerts, seed, pos, hemi) %#ok<DEFNU>
    N = size(G.times, 2);
    if (o < 1) || (o > N)
        error('bst_dynamics:AttachRegion', 'Occurrence %d out of range (N=%d).', o, N);
    end
    G.vertices = i_pad_row(G.vertices, N);
    G.hemi     = i_pad_row(G.hemi,     N);
    G.strength = i_pad_row(G.strength, N);
    G.charge   = i_pad_row(G.charge,   N);
    G.pos      = i_pad_pos(G.pos,      N);
    G.region   = i_pad_cell(G.region,  N);
    G.vertices(o) = double(seed);
    G.hemi(o)     = double(hemi);
    G.pos(o, :)   = double(pos(:)');
    G.region{o}   = double(regionVerts(:)');
    G.type = 'simple';
end

% Pad a [1 x m] numeric row to length N with NaN ([] -> all-NaN).
function v = i_pad_row(v, N)
    if isempty(v),          v = nan(1, N);
    elseif (numel(v) < N),  v(end+1:N) = NaN;  end
end
% Pad a [m x 3] position matrix to N rows with NaN ([] -> all-NaN).
function p = i_pad_pos(p, N)
    if isempty(p),           p = nan(N, 3);
    elseif (size(p,1) < N),  p(end+1:N, :) = NaN;  end
end
% Pad a {1 x m} cell to length N with [] ({} -> all-empty).
function c = i_pad_cell(c, N)
    if isempty(c),          c = cell(1, N);
    elseif (numel(c) < N),  c(end+1:N) = {[]};  end
end
```

- [ ] **Step 5: Run test to verify it passes**

Run:
```matlab
test_dynamics_atoms
```
Expected: `T1b region/AttachRegion: field=1 pad=1 roundtrip=1 => PASS`. The rest of the suite is unchanged (T2/T3 may still SKIP if no kernel link).

- [ ] **Step 6: Commit**

```bash
git add toolbox/db/db_template.m toolbox/dynamics/bst_dynamics.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): per-occurrence region field + bst_dynamics AttachRegion verb"
```

---

### Task 2: Capture region → active atom + region rendering

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (Capture menu item + Record-section button in `CreatePanel`; new `OnCaptureRegion`)
- Modify: `toolbox/gui/view_dynamics.m` (`Redraw`: region patches + partial localization)
- Test: `dev/test_dynamics_atoms.m` (add T6, kernel-gated)

**Interfaces:**
- Consumes: `bst_dynamics('AttachRegion', G, o, regionVerts, seed, pos, hemi)` (Task 1); `panel_scout('GetSelectedScouts')` → `[sSelScouts, iSelScouts, sSurf, iSurf]` where `sSelScouts(1).Vertices`/`.Seed` and `sSurf.FileName`; `file_compare(f1,f2)`; `tess_scout_area(SurfaceFile, Seed, Radius)` → `[Vertices, phi]`; the panel's `i_cs()`, `i_apply(st)`, `st.occMap` (`[g o _]` per right-list row), `st.hFig`, `st.file`, `getappdata(st.hFig,'DynamicsSurf')`.
- Produces: `panel_bst_dynamics('OnCaptureRegion')` writes the selected occurrence's region/seed and redraws+saves. `view_dynamics('Redraw', hFig, T)` now draws `AtomRegion<g>_<o>` patches and tolerates NaN (unlocalized) seed rows.

- [ ] **Step 1: Write the failing test**

In `dev/test_dynamics_atoms.m`, after the T5 block (after `pass = pass && ok5;` near line 133) and **before** the `% cleanup` block, add T6:

```matlab
    % T6: AttachRegion + Redraw -> region patch + seed marker on the cortex (capture render path)
    st = getappdata(0,'DynamicsTarget');  Tt = st.T;
    gPh = find(arrayfun(@(k) ~isempty(Tt.Groups(k).phase) && ~isempty(Tt.Groups(k).times), 1:Tt.nGroups), 1);
    SurfaceFile = Tt.SurfaceFile;
    SurfT = in_tess_bst(SurfaceFile, 0);
    seed = round(size(SurfT.Vertices,1)/3);
    regionVerts = tess_scout_area(SurfaceFile, seed, 0.008);     % 8 mm geodesic disk
    pos = SurfT.Vertices(seed,:);  hemi = 1 + (pos(2) < 0);
    G6  = bst_dynamics('AttachRegion', Tt.Groups(gPh), 1, regionVerts, seed, pos, hemi);
    Tt.Groups(gPh) = G6;  st.T = Tt;  setappdata(0, 'DynamicsTarget', st);
    view_dynamics('Redraw', hFig, st.T);  drawnow;
    nReg6 = numel(findobj(hFig, '-regexp', 'Tag', sprintf('^AtomRegion%d_', gPh)));
    nMk6  = numel(findobj(hFig, 'Tag', sprintf('AtomMarker%d', gPh)));
    fr6 = fullfile(bst_get('BrainstormTmpDir'), 'dyn_region_t6.mat');
    bst_dynamics('Save', fr6, st.T);  Tr6 = bst_dynamics('Load', fr6);
    rt6 = isequal(Tr6.Groups(gPh).region{1}, double(regionVerts(:)'));
    ok6 = ~isempty(regionVerts) && (nReg6>=1) && (nMk6==1) && isfinite(G6.vertices(1)) && rt6;
    fprintf('T6 capture-render: regionVerts=%d patch=%d marker=%d roundtrip=%d => %s\n', numel(regionVerts), nReg6, nMk6, rt6, PF{ok6+1});
    pass = pass && ok6;
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
test_dynamics_atoms
```
Expected: `T6 capture-render: ... patch=0 marker=0 ... => FAIL` — `Redraw` does not yet draw `AtomRegion*` patches, and a phase child with a NaN-padded seed currently errors or draws no `AtomMarker` (old `Redraw` indexes `Surf.VertNormals(vtx,:)` with the whole `vertices` vector).

- [ ] **Step 3: Rewrite `Redraw` to draw regions + tolerate partial localization**

In `toolbox/gui/view_dynamics.m`, replace the entire `Redraw` function (currently ~lines 128–168) with:

```matlab
%% ===== REDRAW MARKERS + REGIONS FROM A (possibly edited) TABLE =====
function Redraw(hFig, T)
    if isempty(hFig) || ~ishandle(hFig), return; end
    hAxes = findobj(hFig, '-depth', 1, 'Tag', 'Axes3D');
    if isempty(hAxes), return; end
    hAxes = hAxes(1);
    set(hAxes, 'NextPlot', 'add');   % low-level line()/patch() never resets the cortex
    % Cache the surface (normals + faces) per figure to avoid reloading on every edit
    Surf = getappdata(hFig, 'DynamicsSurf');
    if isempty(Surf)
        SurfaceFile = T.SurfaceFile;  if isempty(SurfaceFile), SurfaceFile = T.Groups(1).SurfaceFile; end
        Surf = in_tess_bst(SurfaceFile, 0);
        setappdata(hFig, 'DynamicsSurf', Surf);
    end
    hasNorm = isfield(Surf, 'VertNormals') && ~isempty(Surf.VertNormals);
    nVert   = size(Surf.Vertices, 1);
    % Clear old markers, regions, selection
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomMarker'));
    delete(findobj(hAxes, '-regexp', 'Tag', '^AtomRegion'));
    delete(findobj(hAxes, 'Tag', 'AtomSel'));
    % One marker set + region patches per spatial group
    GroupsPosOff = cell(1, numel(T.Groups));
    for g = 1:numel(T.Groups)
        G = T.Groups(g);
        if isempty(G.pos)
            GroupsPosOff{g} = zeros(0,3);   % temporal-only group (no localized occurrence)
            continue;
        end
        nOcc = size(G.pos, 1);
        col  = G.color;  if isempty(col), col = [1 0 1]; end
        % per-occurrence offset seed position; unlocalized (NaN) rows stay NaN -> no marker
        po = nan(nOcc, 3);
        for o = 1:nOcc
            p = G.pos(o, :);
            if any(~isfinite(p)), continue; end
            if hasNorm && (o <= numel(G.vertices)) && isfinite(G.vertices(o))
                nrm = Surf.VertNormals(G.vertices(o), :);
            else
                nrm = p ./ max(sqrt(sum(p.^2)), eps);
            end
            po(o, :) = p + 0.002 * nrm;     % ~2 mm out, clears the surface
        end
        GroupsPosOff{g} = po;
        line(po(:,1), po(:,2), po(:,3), 'Parent', hAxes, ...
            'Marker','o', 'MarkerFaceColor',col, 'MarkerEdgeColor',[0 0 0], ...
            'MarkerSize',6, 'LineStyle','none', 'Tag', sprintf('AtomMarker%d', g));
        % captured regions: a translucent patch of faces fully inside region{o}
        if isfield(G, 'region') && ~isempty(G.region)
            for o = 1:min(nOcc, numel(G.region))
                rv = G.region{o};
                if isempty(rv), continue; end
                inReg = false(nVert, 1);  inReg(rv) = true;
                fIn = all(inReg(Surf.Faces), 2);
                if ~any(fIn), continue; end
                patch('Faces', Surf.Faces(fIn,:), 'Vertices', Surf.Vertices, 'Parent', hAxes, ...
                    'FaceColor', col, 'FaceAlpha', 0.35, 'EdgeColor', 'none', ...
                    'Tag', sprintf('AtomRegion%d_%d', g, o));
            end
        end
    end
    setappdata(hFig, 'GroupsPosOff', GroupsPosOff);
    % Selection marker (hidden until an occurrence is picked)
    line(NaN, NaN, NaN, 'Parent', hAxes, ...
        'Marker','o', 'MarkerSize',13, 'MarkerEdgeColor',[1 1 0], ...
        'LineWidth',2, 'LineStyle','none', 'Visible','off', 'Tag','AtomSel');
end
```

- [ ] **Step 4: Add the Capture UI + `OnCaptureRegion` to the panel**

In `toolbox/gui/panel_bst_dynamics.m`, `CreatePanel`:

(a) Add a Capture item to the Atoms menu, right after the `'Record at cursor'` menu item (after line ~75):
```matlab
    gui_component('MenuItem', jMenuAtoms, [], 'Capture region -> active atom', IconLoader.ICON_SCOUT_NEW, [], @(h,e)bst_call(@OnCaptureRegion));
```

(b) Add a Capture button to the Record section, right after the `'Record at cursor'` button (after line ~135, the `jCtrl.add(jRec);` line should stay last):
```matlab
    gui_component('button', jRec, 'br hfill', 'Capture region -> active atom', [], 'Snapshot the selected Scout''s vertices into the selected atom (localizes a time-only marker)', @(h,e)bst_call(@OnCaptureRegion));
```

(c) Add the `OnCaptureRegion` function. Place it immediately after `OnRecord` (after line ~358):
```matlab
%% ===== CAPTURE: the selected Scout's geodesic region -> the active (selected) atom =====
% The active atom is the occurrence selected in the right-hand list. Snapshots the currently
% selected Scout's vertices into that occurrence (region + seed), localizing a time-only marker.
function OnCaptureRegion() %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~ishandle(st.hFig), return; end
    if isempty(st.occMap)
        java_dialog('warning', 'Select an atom in the list first.', 'Capture region');  return;
    end
    row = ctrl.jListOccur.getSelectedIndex() + 1;
    if (row < 1) || (row > size(st.occMap,1))
        java_dialog('warning', 'Select an atom in the list first.', 'Capture region');  return;
    end
    g = st.occMap(row,1);  o = st.occMap(row,2);
    if (size(st.T.Groups(g).times,1) ~= 1)
        java_dialog('warning', 'Select a single atom, not a time window.', 'Capture region');  return;
    end
    % the geodesic region = the currently selected Scout (grown with the Scout "Area" tool)
    [sScout, ~, sSurf] = panel_scout('GetSelectedScouts');
    if isempty(sScout) || isempty(sScout(1).Vertices)
        java_dialog('warning', 'Grow a region with the Scout "Area" tool first.', 'Capture region');  return;
    end
    sScout = sScout(1);
    if ~isempty(sSurf) && ~isempty(st.T.SurfaceFile) && ~file_compare(sSurf.FileName, st.T.SurfaceFile)
        java_dialog('warning', 'The selected region is on a different surface than the atoms.', 'Capture region');  return;
    end
    if ~isempty(sScout.Seed), seed = double(sScout.Seed(1)); else, seed = double(sScout.Vertices(1)); end
    Surf = getappdata(st.hFig, 'DynamicsSurf');
    if isempty(Surf), Surf = in_tess_bst(st.T.SurfaceFile, 0);  setappdata(st.hFig, 'DynamicsSurf', Surf); end
    pos  = Surf.Vertices(seed, :);
    hemi = 1 + (pos(2) < 0);                                       % SCS Y>0 = left
    st.T.Groups(g) = bst_dynamics('AttachRegion', st.T.Groups(g), o, sScout.Vertices, seed, pos, hemi);
    i_apply(st);                                                   % redraw markers/regions + rebuild tree
    if ~isempty(st.file), try, bst_dynamics('Save', st.file, st.T); catch, end; end %#ok<CTCH>
    bst_progress('text', sprintf('Captured %d-vertex region into "%s"', numel(sScout.Vertices), st.T.Groups(g).label));
end
```

- [ ] **Step 5: Run test to verify it passes**

Run:
```matlab
test_dynamics_atoms
```
Expected: `T6 capture-render: regionVerts=<k> patch=1 marker=1 roundtrip=1 => PASS`, and T1b plus the prior suite still pass.

- [ ] **Step 6: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/gui/view_dynamics.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): capture Scout geodesic region into the active atom + render region patches"
```

---

### Task 3: Show-phases filter (list + cortex)

**Files:**
- Modify: `toolbox/gui/panel_bst_dynamics.m` (Show-phases submenu + `jPhaseItems` control; `st.showPhase` state; `OnTogglePhase`; `i_phase_index`; filter in `i_window_atoms`; pass `showPhase` from `i_apply`)
- Modify: `toolbox/gui/view_dynamics.m` (`Redraw` gains an optional `showPhase` arg + its own `i_phase_index`)
- Test: `dev/test_dynamics_atoms.m` (add T7, kernel-gated)

**Interfaces:**
- Consumes: `st.showPhase` `[1×4]` logical (peak/trough/rising/falling); the `BstPanel` controls now include `jPhaseItems` (a `javax.swing.JCheckBoxMenuItem[4]`); `i_apply(st)`; `i_window_atoms(T, gBand, w, showPhase)`; `view_dynamics('Redraw', hFig, T, showPhase)`.
- Produces: `panel_bst_dynamics('OnTogglePhase', iPhase)` flips `st.showPhase(iPhase)` from the menu item state and refreshes. The window list and cortex omit phases toggled off.

- [ ] **Step 1: Write the failing test**

In `dev/test_dynamics_atoms.m`, after the T6 block (after `pass = pass && ok6;`) and before `% cleanup`, add T7:

```matlab
    % T7: Show-phases filter hides peak rows (list) + peak markers/regions (cortex); non-destructive
    ctrl = bst_get('PanelControls', 'Dynamics');
    st = getappdata(0,'DynamicsTarget');  Tt = st.T;  parents = {Tt.Groups.parent};
    gW = find(cellfun(@isempty,parents) & arrayfun(@(k) strcmp(Tt.Groups(k).bandName,'alpha') && size(Tt.Groups(k).times,1)==2, 1:Tt.nGroups), 1);
    gPeak = find(arrayfun(@(k) strcmpi(i_t7_str(Tt.Groups(k).phase),'peak') && strcmpi(i_t7_str(Tt.Groups(k).bandName),'alpha'), 1:Tt.nGroups), 1);
    % localize peak occurrence 1 so it actually has a cortex marker + region to hide
    SurfT = in_tess_bst(Tt.SurfaceFile, 0);  seed = round(size(SurfT.Vertices,1)/2);
    rvP = tess_scout_area(Tt.SurfaceFile, seed, 0.008);
    Tt.Groups(gPeak) = bst_dynamics('AttachRegion', Tt.Groups(gPeak), 1, rvP, seed, SurfT.Vertices(seed,:), 1+(SurfT.Vertices(seed,2)<0));
    st.T = Tt;  setappdata(0,'DynamicsTarget',st);
    view_dynamics('Redraw', hFig, st.T, st.showPhase);  drawnow;
    nPeakMkOn  = numel(findobj(hFig, 'Tag', sprintf('AtomMarker%d', gPeak)));
    nPeakRegOn = numel(findobj(hFig, '-regexp', 'Tag', sprintf('^AtomRegion%d_', gPeak)));
    nPeakRows  = i_t7_selwin_peakrows(ctrl, gW);       % select window 1, count 'peak' rows
    % toggle peak OFF
    ctrl.jPhaseItems(1).setSelected(false);
    panel_bst_dynamics('OnTogglePhase', 1);  drawnow;
    nPeakOff   = i_t7_selwin_peakrows(ctrl, gW);
    nPeakMkOff = numel(findobj(hFig, 'Tag', sprintf('AtomMarker%d', gPeak)));
    nPeakRegOff= numel(findobj(hFig, '-regexp', 'Tag', sprintf('^AtomRegion%d_', gPeak)));
    % toggle peak back ON -> restored
    ctrl.jPhaseItems(1).setSelected(true);
    panel_bst_dynamics('OnTogglePhase', 1);  drawnow;
    nPeakBack  = i_t7_selwin_peakrows(ctrl, gW);
    ok7 = (nPeakRows>0) && (nPeakMkOn==1) && (nPeakRegOn>=1) ...
       && (nPeakOff==0) && (nPeakMkOff==0) && (nPeakRegOff==0) && (nPeakBack==nPeakRows);
    fprintf('T7 phase-filter: rowsOn=%d mkOn=%d regOn=%d | rowsOff=%d mkOff=%d | back=%d => %s\n', ...
        nPeakRows, nPeakMkOn, nPeakRegOn, nPeakOff, nPeakMkOff, nPeakBack, PF{ok7+1});
    pass = pass && ok7;
```

Then add these two test helpers at the end of the file (after `i_find_kernel`, near line 174):

```matlab
%% ===== T7 HELPERS =====
function s = i_t7_str(x)
    if isempty(x), s = ''; else, s = char(x); end
end
% Select window w=1 of band group gW in the tree, return how many right-list rows contain 'peak'.
function n = i_t7_selwin_peakrows(ctrl, gW)
    st = getappdata(0, 'DynamicsTarget');
    iWinNode = find(arrayfun(@(k) strcmp(st.nodeInfo(k).kind,'window') && st.nodeInfo(k).g==gW && st.nodeInfo(k).w==1, 1:numel(st.nodeInfo)), 1);
    ctrl.jTree.setSelectionPath(javax.swing.tree.TreePath(st.nodeList{iWinNode}.getPath()));  drawnow;
    model = ctrl.jListOccur.getModel();  n = 0;
    for r = 0:(model.getSize()-1)
        if ~isempty(strfind(char(model.getElementAt(r)), 'peak')), n = n + 1; end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
test_dynamics_atoms
```
Expected: FAIL — `ctrl.jPhaseItems` does not exist (the submenu/control is not built yet) and `panel_bst_dynamics('OnTogglePhase', 1)` is an unknown verb, so the test errors or reports `=> FAIL`.

- [ ] **Step 3: Build the Show-phases submenu + control**

In `toolbox/gui/panel_bst_dynamics.m`, `CreatePanel`, add the submenu inside the Atoms menu. Insert it after the `'Set color'` item and before the first `jMenuAtoms.addSeparator();` (after line ~69):

```matlab
    jMenuAtoms.addSeparator();
    jMenuPhases = gui_component('Menu', jMenuAtoms, [], 'Show phases', IconLoader.ICON_EVT_TYPE, [], []);
    phaseNames  = {'peak','trough','rising','falling'};
    jPhaseItems = javaArray('javax.swing.JCheckBoxMenuItem', 4);
    for ip = 1:4
        jit = gui_component('checkboxmenuitem', jMenuPhases, [], phaseNames{ip}, [], [], @(h,e)bst_call(@()OnTogglePhase(ip)));
        jit.setSelected(true);
        jPhaseItems(ip) = jit;
    end
```

Then add `jPhaseItems` to the `BstPanel(...)` struct at the end of `CreatePanel` (extend the existing struct argument list, after `'jPeaks',jPeaks`):

```matlab
        'jSpaceParams',jSpaceParams, 'jSpacePot',jSpacePot, 'jSpaceStr',jSpaceStr, 'jPeaks',jPeaks, 'jPhaseItems',jPhaseItems));
```

- [ ] **Step 4: Add `showPhase` state, `OnTogglePhase`, and `i_phase_index`**

In `toolbox/gui/panel_bst_dynamics.m`, `SetTarget` (line ~415), add `showPhase` to the target struct so the state exists from the start:

```matlab
    setappdata(0, 'DynamicsTarget', struct('hFig',hFig, 'T',T, 'file',file, 'curGroup',0, ...
        'nodeList',{ {} }, 'nodeInfo',[], 'occMap',[], 'Lambda',[], 'showPhase',[1 1 1 1]));
```

Add `OnTogglePhase` and `i_phase_index` (place after `OnCaptureRegion`):

```matlab
%% ===== SHOW-PHASES FILTER (display-only; never deletes atoms) =====
function OnTogglePhase(ip) %#ok<DEFNU>
    [ctrl, st] = i_cs();
    if isempty(ctrl) || isempty(st) || ~isfield(ctrl,'jPhaseItems') || isempty(ctrl.jPhaseItems), return; end
    sp = i_field(st, 'showPhase', [1 1 1 1]);
    sp(ip) = ctrl.jPhaseItems(ip).isSelected();
    st.showPhase = sp;  setappdata(0, 'DynamicsTarget', st);
    i_apply(st);                                                  % rebuild list + redraw cortex
end

% Phase name -> filter index (peak=1 trough=2 rising=3 falling=4; 0 = not a phase group -> always shown).
function k = i_phase_index(ph)
    switch lower(i_str(ph))
        case 'peak',    k = 1;
        case 'trough',  k = 2;
        case 'rising',  k = 3;
        case 'falling', k = 4;
        otherwise,      k = 0;
    end
end
```

- [ ] **Step 5: Filter the window list + pass `showPhase` to Redraw**

(a) In `i_window_atoms` (line ~542), change the signature and skip filtered phases. Replace the function header and the child loop guard:

Change the header from `function [rows, occMap] = i_window_atoms(T, gBand, w)` to:
```matlab
function [rows, occMap] = i_window_atoms(T, gBand, w, showPhase)
    if (nargin < 4) || isempty(showPhase), showPhase = [1 1 1 1]; end
```
Inside the `for c = children(:)'` loop, immediately after `Gc = T.Groups(c);`, add:
```matlab
        pk = i_phase_index(Gc.phase);
        if (pk >= 1) && ~showPhase(pk), continue; end          % phase filtered out
```

(b) Update the call site in `TreeSel_Callback` (line ~500) to pass the filter:
```matlab
            [rows, occMap] = i_window_atoms(st.T, info.g, info.w, i_field(st,'showPhase',[1 1 1 1]));
```

(c) In `i_apply` (line ~649), pass `showPhase` into Redraw:
```matlab
function i_apply(st)
    setappdata(0, 'DynamicsTarget', st);
    if ~isempty(st.hFig) && ishandle(st.hFig)
        try, view_dynamics('Redraw', st.hFig, st.T, i_field(st,'showPhase',[1 1 1 1])); catch, end %#ok<CTCH>
    end
    BuildTree();
end
```

- [ ] **Step 6: Honor `showPhase` in `view_dynamics('Redraw')`**

In `toolbox/gui/view_dynamics.m`, add the optional arg to `Redraw` and skip filtered phase groups.

Change the header from `function Redraw(hFig, T)` to:
```matlab
function Redraw(hFig, T, showPhase)
    if (nargin < 3) || isempty(showPhase), showPhase = [1 1 1 1]; end
```
Inside the `for g = 1:numel(T.Groups)` loop, immediately after `G = T.Groups(g);`, add the phase-filter skip (before the `if isempty(G.pos)` check):
```matlab
        pk = i_phase_index(G.phase);
        if (pk >= 1) && ~showPhase(pk)
            GroupsPosOff{g} = zeros(0,3);   % phase filtered out: no marker, no region
            continue;
        end
```
Add a local `i_phase_index` helper at the end of `view_dynamics.m` (after `AtomsFromResult`/`Redraw`):
```matlab
% Phase name -> filter index (peak=1 trough=2 rising=3 falling=4; 0 = not a phase group).
function k = i_phase_index(ph)
    if isempty(ph), k = 0; return; end
    switch lower(char(ph))
        case 'peak',    k = 1;
        case 'trough',  k = 2;
        case 'rising',  k = 3;
        case 'falling', k = 4;
        otherwise,      k = 0;
    end
end
```

- [ ] **Step 7: Run test to verify it passes**

Run:
```matlab
test_dynamics_atoms
```
Expected: `T7 phase-filter: rowsOn=<n> mkOn=1 regOn=<m> | rowsOff=0 mkOff=0 | back=<n> => PASS`, and the full suite ends `==== SUITE: PASS ====` (8 checks: T1, T1b, T2, T3, T4, T5, T6, T7 — T2/T3/T6/T7 SKIP together when no kernel link).

- [ ] **Step 8: Commit**

```bash
git add toolbox/gui/panel_bst_dynamics.m toolbox/gui/view_dynamics.m dev/test_dynamics_atoms.m
git commit -m "feat(dynamics): Show-phases filter for the atom list + cortex markers/regions"
```

---

## Self-Review

**1. Spec coverage:**
- `region` field + back-compat → Task 1 (template + Load normalization via existing `struct_copy_fields`).
- `AttachRegion` verb (padding + localization) → Task 1.
- Capture reuse of Scout Area tool + `GetSelectedScouts` snapshot + guard rails (no occurrence / no scout / surface mismatch / window-not-atom) → Task 2.
- Capture exposed twice (Atoms menu + Record button) → Task 2.
- Region patch rendering + partial-localization tolerance → Task 2.
- Show-phases submenu (checkboxmenuitem) → Task 3.
- Filter applied in `i_window_atoms` (list) AND `Redraw` (cortex), non-destructive → Task 3.
- Tests T6 (capture render) + T7 (phase filter) → Tasks 2 & 3; model test T1b → Task 1.
All spec sections map to a task. No gaps.

**2. Placeholder scan:** No TBD/TODO; every code step shows complete code; every run step shows the exact command and expected line.

**3. Type consistency:** `AttachRegion(G,o,regionVerts,seed,pos,hemi)` signature is identical in Task 1 (definition), the T1b/T6/T7 tests, and `OnCaptureRegion`. `region` is a cell `{1×N}` everywhere. `showPhase` is `[1×4]` in `SetTarget`, `OnTogglePhase`, `i_window_atoms`, `i_apply`, and `Redraw`. `i_phase_index` returns `peak=1…falling=4, else 0` in both `panel_bst_dynamics` and `view_dynamics`. `Redraw` is called as `('Redraw', hFig, T)` (Task 2 test, default all-on) and `('Redraw', hFig, T, showPhase)` (Task 3) — the optional-arg default covers both. Tags `AtomMarker%d` / `AtomRegion%d_%d` match between `Redraw` and the test `findobj` queries.
