# Atom Localization Accessor (Phase 1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add `bst_atom` — a uniform `(center, extent, weighting)` Localization accessor that reads and writes any of the four atom axes (time, frequency, source, scale) over the existing `atomgroup` storage, so every later phase can treat the axes uniformly.

**Architecture:** A new I/O-free `toolbox/dynamics/bst_atom.m` (verb-dispatched like `bst_dynamics`) maps between the heterogeneous stored fields (`times` onset/offset, `band` `[fLo fHi]`, `vertices`+new `radius`, `scale` `[k1 k2]`) and a uniform Localization struct. Frequency and scale are group-level; time and source are per-occurrence. One new backward-compatible `radius` field gives the source axis its window parameter. No panel/detector/UI changes — those are Phase 2+.

**Tech Stack:** MATLAB, Brainstorm; `eval(macro_method)` verb dispatch; `db_template('atomgroup')`; tested headless via a `dev/test_bst_atom.m` suite with Brainstorm running.

## Global Constraints

- This is Phase 1 of the architecture analysis `docs/superpowers/specs/2026-06-24-atom-tensor-architecture-analysis.md`. **Non-breaking:** storage fields are unchanged except one additive `radius` field; existing `dynamics_*` tables and the `test_dynamics_atoms` suite must keep passing.
- The Localization primitive is `(center, extent, weighting)`; `weighting` defaults to `'hard'` (the `'soft'`/wavelet form is future — reserve the field, do not implement decay).
- Three per-axis states: `'unlocalized'` (axis not pinned), `'point'` (extent 0), `'window'` (extent > 0).
- Axis units and locality are fixed: time = seconds (per-occurrence), frequency = Hz (group-level), source = vertex id + 3-D pos, extent = geodesic radius in metres (per-occurrence), scale = eigenvalue (group-level).
- `center`/`extent` are numeric; `label` is an optional human-readable string (`bandName`/`scaleName` today; the atlas layer is future).
- New file is I/O-free (operates on an in-memory group struct), verb-dispatched via `eval(macro_method)`, lives in `toolbox/dynamics/` (auto-added to path by `brainstorm.m` genpath).
- Do not start implementation on `development`; branch first.
- Tests run headless in MATLAB with Brainstorm live (`brainstorm nogui`, TutorialAuditory). Do not `clear` or restart Brainstorm; `rehash` and re-run. The integration test reuses the unconstrained-kernel fixture and SKIPs with the suite when absent.

---

### Task 1: `radius` field + `bst_atom('Get')` read accessor

**Files:**
- Modify: `toolbox/db/db_template.m` (add `radius` to the `atomgroup` template, after `region` ~line 349)
- Create: `toolbox/dynamics/bst_atom.m`
- Create/Test: `dev/test_bst_atom.m`

**Interfaces:**
- Consumes: `db_template('atomgroup')` (existing fields `times [1|2 x N]`, `band [1x2]`, `bandName`, `scale [1x2]`, `scaleName`, `vertices [1xN]`, `pos [Nx3]`, `region {1xN}`).
- Produces:
  - `atomgroup.radius` — `[1 x N]` per-occurrence geodesic radius in metres (source extent); `[]`/`NaN` when unset.
  - `A = bst_atom('Axes')` — `1x4` struct array, fields `name` (`'time'|'freq'|'source'|'scale'`), `perOcc` (logical), `unit` (char).
  - `loc = bst_atom('NewLoc', axis)` — empty Localization: fields `axis,center,extent,weighting,label,state,pos` = `('',NaN,NaN,'hard','','unlocalized',[])`.
  - `loc = bst_atom('Get', G, axis, occ)` — Localization for that axis (occ optional, default 1; ignored for group-level axes). `state ∈ {'unlocalized','point','window'}`; `pos` set only for `'source'`.

- [ ] **Step 1: Write the failing test**

Create `dev/test_bst_atom.m`:

```matlab
function test_bst_atom()
% TEST_BST_ATOM: unit + round-trip tests for the (center,extent) localization accessor.
%
% USAGE:  test_bst_atom   % Brainstorm running
%
% Authors: Diellor Basha, 2026

    PF = {'FAIL','PASS'};  pass = true;

    % ---------- T1: Axes metadata + Get on a synthetic group (4 axes, 3 states) ----------
    A = bst_atom('Axes');
    axesOK = (numel(A)==4) && strcmp(A(1).name,'time') && strcmp(A(2).name,'freq') ...
          && strcmp(A(3).name,'source') && strcmp(A(4).name,'scale') ...
          && (A(1).perOcc && ~A(2).perOcc && A(3).perOcc && ~A(4).perOcc);

    G = bst_dynamics('NewGroup', 'g');
    G.times = [0.10 0.20; 0.30 0.40];      % extended (2xN): occ1 window [.10 .30]->c=.20 w=.10
    G.band  = [8 12];                       % freq window: c=10 w=2 (alpha)
    G.bandName = 'alpha';
    G.scale = [];                           % scale unlocalized
    G.vertices = [100 200];  G.pos = [1 2 3; 4 5 6];  G.radius = [0.005 0];  % occ1 window r=5mm, occ2 point

    lt = bst_atom('Get', G, 'time', 1);
    okT = strcmp(lt.state,'window') && abs(lt.center-0.20)<1e-9 && abs(lt.extent-0.10)<1e-9;
    lf = bst_atom('Get', G, 'freq');        % group-level -> occ ignored
    okF = strcmp(lf.state,'window') && (lf.center==10) && (lf.extent==2) && strcmp(lf.label,'alpha');
    ls1 = bst_atom('Get', G, 'source', 1);
    okS1 = strcmp(ls1.state,'window') && (ls1.center==100) && abs(ls1.extent-0.005)<1e-12 && isequal(ls1.pos,[1 2 3]);
    ls2 = bst_atom('Get', G, 'source', 2);
    okS2 = strcmp(ls2.state,'point') && (ls2.center==200) && (ls2.extent==0);
    lk = bst_atom('Get', G, 'scale');
    okK = strcmp(lk.state,'unlocalized') && ~isfinite(lk.center);
    % weighting default + radius field present
    okW = strcmp(lt.weighting,'hard') && isfield(db_template('atomgroup'),'radius');

    ok1 = axesOK && okT && okF && okS1 && okS2 && okK && okW;
    fprintf('T1 Get: axes=%d time=%d freq=%d src1=%d src2=%d scale=%d weight/field=%d => %s\n', ...
        axesOK, okT, okF, okS1, okS2, okK, okW, PF{ok1+1});
    pass = pass && ok1;

    fprintf('\n==== SUITE: %s ====\n', PF{pass+1});
end
```

- [ ] **Step 2: Run test to verify it fails**

Run (MATLAB, Brainstorm live):
```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev'); rehash; test_bst_atom
```
Expected: FAIL/error — `bst_atom` is undefined and `db_template('atomgroup')` has no `radius` field.

- [ ] **Step 3: Add the `radius` field to the template**

In `toolbox/db/db_template.m`, in `case 'atomgroup'`, add after the `region` line (~349):

```matlab
            'radius',      [], ...      % SPACE: [1 x N] per-occurrence geodesic radius [m] (source extent); []/NaN = unset
```

- [ ] **Step 4: Create `bst_atom.m` with `Axes`/`NewLoc`/`Get`**

Create `toolbox/dynamics/bst_atom.m`:

```matlab
function varargout = bst_atom( varargin )
% BST_ATOM: uniform (center, extent, weighting) localization accessor over an atom group.
%
% Every atom axis (time, frequency, source, scale) localizes the same way -- a center plus a
% window (extent), or a bare point. This accessor maps between an atomgroup's heterogeneous
% stored fields and a uniform Localization struct, so the panel and detectors can treat all
% four axes identically. I/O-free (operates on an in-memory group).
%
% Localization struct:
%   .axis      'time'|'freq'|'source'|'scale'
%   .center    numeric center (s | Hz | vertex id | eigenvalue)
%   .extent    numeric half-window (s | Hz | metres geodesic radius | eigenvalue); 0 = point
%   .weighting 'hard' (default) | 'soft'   (soft = wavelet decay; reserved, future)
%   .label     optional human-readable name ('' today; atlas layer is future)
%   .state     'unlocalized' | 'point' | 'window'
%   .pos       [1x3] seed position, source axis only (else [])
%
% USAGE:
%   A   = bst_atom('Axes')                 % axis metadata (name/perOcc/unit)
%   loc = bst_atom('NewLoc', axis)         % empty localization
%   loc = bst_atom('Get', G, axis, occ)    % read (occ default 1; ignored for group axes)
%   G   = bst_atom('Set', G, axis, occ, loc)  % write (Task 2)
%
% SEE ALSO: bst_dynamics, panel_bst_dynamics
%
% Authors: Diellor Basha, 2026

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm
%
% Copyright (c) University of Southern California & McGill University
% This software is distributed under the terms of the GNU General Public License
% as published by the Free Software Foundation. Further details on the GPLv3
% license can be found at http://www.gnu.org/copyleft/gpl.html.
%
% FOR RESEARCH PURPOSES ONLY. THE SOFTWARE IS PROVIDED "AS IS," AND THE
% UNIVERSITY OF SOUTHERN CALIFORNIA AND ITS COLLABORATORS DO NOT MAKE ANY
% WARRANTY, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO WARRANTIES OF
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE, NOR DO THEY ASSUME ANY
% LIABILITY OR RESPONSIBILITY FOR THE USE OF THIS SOFTWARE.
%
% For more information type "brainstorm license" at command prompt.
% =============================================================================@

eval(macro_method);
end


%% ===== AXES: canonical axis metadata =====
function A = Axes()
    A = struct('name',   {'time','freq','source','scale'}, ...
               'perOcc', {true,   false,  true,    false}, ...
               'unit',   {'s',    'Hz',   'vertex','eigenvalue'});
end


%% ===== NEW (empty) localization =====
function loc = NewLoc(axis)
    if (nargin < 1), axis = ''; end
    loc = struct('axis',axis, 'center',NaN, 'extent',NaN, 'weighting','hard', ...
                 'label','', 'state','unlocalized', 'pos',[]);
end


%% ===== GET: read one axis localization from the group =====
function loc = Get(G, axis, occ)
    if (nargin < 3) || isempty(occ), occ = 1; end
    loc = NewLoc(axis);
    switch axis
        case 'time'
            if isempty(G.times) || (occ > size(G.times,2)), return; end
            col = G.times(:, occ);
            if any(~isfinite(col)), return; end
            if (size(G.times,1) >= 2)
                loc.center = mean(col(1:2));  loc.extent = (col(2) - col(1)) / 2;
            else
                loc.center = col(1);          loc.extent = 0;
            end
            loc.state = i_state(loc.extent);
        case 'freq'
            if (numel(G.band) < 2), return; end
            loc.center = mean(G.band(1:2));  loc.extent = (G.band(2) - G.band(1)) / 2;
            loc.label  = G.bandName;         loc.state  = i_state(loc.extent);
        case 'source'
            if isempty(G.vertices) || (occ > numel(G.vertices)) || ~isfinite(G.vertices(occ)), return; end
            loc.center = double(G.vertices(occ));
            if ~isempty(G.pos) && (occ <= size(G.pos,1)), loc.pos = G.pos(occ, :); end
            hasR = isfield(G,'radius') && ~isempty(G.radius) && (occ <= numel(G.radius)) && isfinite(G.radius(occ));
            hasReg = ~isempty(G.region) && (occ <= numel(G.region)) && ~isempty(G.region{occ});
            if hasR
                loc.extent = G.radius(occ);  loc.state = i_state(loc.extent);
            elseif hasReg
                loc.extent = NaN;  loc.state = 'window';   % region materialized but radius unrecorded
            else
                loc.extent = 0;    loc.state = 'point';
            end
        case 'scale'
            if (numel(G.scale) < 2), return; end
            loc.center = mean(G.scale(1:2));  loc.extent = (G.scale(2) - G.scale(1)) / 2;
            loc.label  = G.scaleName;         loc.state  = i_state(loc.extent);
        otherwise
            error('bst_atom:Get', 'Unknown axis "%s".', axis);
    end
end


%% ===== state from extent =====
function s = i_state(extent)
    if ~isfinite(extent),   s = 'unlocalized';
    elseif (extent == 0),   s = 'point';
    else,                   s = 'window';
    end
end
```

- [ ] **Step 5: Run test to verify it passes**

Run:
```matlab
rehash; test_bst_atom
```
Expected: `T1 Get: axes=1 time=1 freq=1 src1=1 src2=1 scale=1 weight/field=1 => PASS` and `==== SUITE: PASS ====`.

- [ ] **Step 6: Verify the additive field did not break the atom suite**

Run:
```matlab
test_dynamics_atoms
```
Expected: `==== SUITE: PASS ====` (8/8). The new `radius` template field rides along via `struct_copy_fields`; the schema-equality checks (T1, T6) still hold because the field is in the template.

- [ ] **Step 7: Commit**

```bash
git add toolbox/db/db_template.m toolbox/dynamics/bst_atom.m dev/test_bst_atom.m
git commit -m "feat(atom): bst_atom Get accessor + radius field (uniform center/extent read)"
```

---

### Task 2: `bst_atom('Set')` write accessor + round-trip

**Files:**
- Modify: `toolbox/dynamics/bst_atom.m` (add `Set` + helpers)
- Modify: `dev/test_bst_atom.m` (add T2 round-trip)

**Interfaces:**
- Consumes: `bst_atom('Get'/'NewLoc'/'Axes')` (Task 1); the `radius` field (Task 1).
- Produces: `G = bst_atom('Set', G, axis, occ, loc)` — writes `loc.center`/`loc.extent` (and `loc.label`, and `loc.pos` for source) into the underlying field(s). For group-level axes (`freq`,`scale`) `occ` is ignored. Time promotes a simple group to extended when `extent > 0` (other occurrences become zero-width windows). Keeps `G.type` consistent with `size(times,1)`. Pads per-occurrence arrays (`times`,`vertices`,`pos`,`radius`) to `occ` with NaN when writing past the current length.

- [ ] **Step 1: Write the failing test**

In `dev/test_bst_atom.m`, before the final `fprintf('\n==== SUITE` line, add:

```matlab
    % ---------- T2: Set then Get round-trips on every axis; writes the right fields ----------
    G2 = bst_dynamics('NewGroup', 'g2');
    G2.times = [0.5];   % simple, 1 occurrence
    % time: set a window on occ 1 -> promotes to extended
    G2 = bst_atom('Set', G2, 'time', 1, i_loc('time', 0.40, 0.10));
    rtT = (size(G2.times,1)==2) && strcmp(G2.type,'extended');
    gt  = bst_atom('Get', G2, 'time', 1);
    rtT = rtT && abs(gt.center-0.40)<1e-9 && abs(gt.extent-0.10)<1e-9;
    % freq: set alpha band (group-level)
    G2 = bst_atom('Set', G2, 'freq', [], i_loc_lbl('freq', 10, 2, 'alpha'));
    gf  = bst_atom('Get', G2, 'freq');
    rtF = isequal(G2.band,[8 12]) && abs(gf.center-10)<1e-9 && abs(gf.extent-2)<1e-9 && strcmp(gf.label,'alpha');
    % source: set seed + radius on occ 1 (with pos)
    lcS = i_loc('source', 250, 0.006);  lcS.pos = [7 8 9];
    G2 = bst_atom('Set', G2, 'source', 1, lcS);
    gs  = bst_atom('Get', G2, 'source', 1);
    rtS = (G2.vertices(1)==250) && abs(G2.radius(1)-0.006)<1e-12 && isequal(G2.pos(1,:),[7 8 9]) ...
       && strcmp(gs.state,'window') && abs(gs.extent-0.006)<1e-12;
    % scale: set eigen-band (group-level)
    G2 = bst_atom('Set', G2, 'scale', [], i_loc_lbl('scale', 50, 10, 'gyrus'));
    gk  = bst_atom('Get', G2, 'scale');
    rtK = isequal(G2.scale,[40 60]) && abs(gk.center-50)<1e-9 && abs(gk.extent-10)<1e-9 && strcmp(gk.label,'gyrus');

    ok2 = rtT && rtF && rtS && rtK;
    fprintf('T2 Set round-trip: time=%d freq=%d source=%d scale=%d => %s\n', rtT, rtF, rtS, rtK, PF{ok2+1});
    pass = pass && ok2;
```

And add these test helpers at the end of `dev/test_bst_atom.m` (after the main function `end`):

```matlab
function loc = i_loc(axis, c, w)
    loc = bst_atom('NewLoc', axis);  loc.center = c;  loc.extent = w;
end
function loc = i_loc_lbl(axis, c, w, lbl)
    loc = i_loc(axis, c, w);  loc.label = lbl;
end
```

- [ ] **Step 2: Run test to verify it fails**

Run:
```matlab
rehash; test_bst_atom
```
Expected: FAIL/error — `bst_atom('Set', ...)` is an unknown verb (`Unknown command 'Set'`).

- [ ] **Step 3: Add `Set` + helpers to `bst_atom.m`**

In `toolbox/dynamics/bst_atom.m`, add after the `Get` block (before `i_state`):

```matlab
%% ===== SET: write one axis localization into the group =====
function G = Set(G, axis, occ, loc)
    if (nargin < 4), error('bst_atom:Set','Set requires (G, axis, occ, loc).'); end
    if isempty(occ), occ = 1; end
    c = loc.center;  w = loc.extent;  if ~isfinite(w), w = 0; end
    switch axis
        case 'time'
            G.times = i_pad_cols(G.times, occ);
            if (w > 0) && (size(G.times,1) < 2)
                G.times = [G.times; G.times];          % promote simple->extended (others zero-width)
            end
            if (size(G.times,1) >= 2)
                G.times(1, occ) = c - w;  G.times(2, occ) = c + w;
            else
                G.times(1, occ) = c;
            end
            G.type = i_type(G.times);
        case 'freq'
            G.band = [c - w, c + w];
            if ~isempty(loc.label), G.bandName = loc.label; end
        case 'source'
            G.vertices = i_pad_vec(G.vertices, occ);  G.vertices(occ) = c;
            G.radius   = i_pad_vec(G.radius,   occ);  G.radius(occ)   = w;
            if ~isempty(loc.pos)
                G.pos = i_pad_pos(G.pos, occ);  G.pos(occ, :) = loc.pos(:)';
            end
        case 'scale'
            G.scale = [c - w, c + w];
            if ~isempty(loc.label), G.scaleName = loc.label; end
        otherwise
            error('bst_atom:Set', 'Unknown axis "%s".', axis);
    end
end

% group type consistent with the times row count
function t = i_type(times)
    if (size(times,1) >= 2), t = 'extended'; else, t = 'simple'; end
end
% pad a [r x m] times matrix to >= n columns with NaN
function M = i_pad_cols(M, n)
    if isempty(M), M = nan(1, n); elseif (size(M,2) < n), M(:, end+1:n) = NaN; end
end
% pad a [1 x m] row vector to >= n with NaN
function v = i_pad_vec(v, n)
    if isempty(v), v = nan(1, n); elseif (numel(v) < n), v(end+1:n) = NaN; end
end
% pad a [m x 3] position matrix to >= n rows with NaN
function p = i_pad_pos(p, n)
    if isempty(p), p = nan(n, 3); elseif (size(p,1) < n), p(end+1:n, :) = NaN; end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run:
```matlab
rehash; test_bst_atom
```
Expected: `T2 Set round-trip: time=1 freq=1 source=1 scale=1 => PASS` and `==== SUITE: PASS ====`.

- [ ] **Step 5: Commit**

```bash
git add toolbox/dynamics/bst_atom.m dev/test_bst_atom.m
git commit -m "feat(atom): bst_atom Set accessor + Get/Set round-trip (4 axes)"
```

---

### Task 3: Integration — read a real Detect + Record table through `bst_atom`

**Files:**
- Modify: `dev/test_bst_atom.m` (add T3, kernel-gated)

**Interfaces:**
- Consumes: `bst_atom('Get')` (Task 1); `process_evt_refphase('Compute')`, `bst_dynamics`; the unconstrained-kernel fixture loader pattern from `dev/test_dynamics_atoms.m` (`i_find_kernel`).
- Produces: nothing new — proves the accessor reads what the real detectors/recorders produce (freq from the band, time from the phase-marker windows/points, source from a recorded occurrence).

- [ ] **Step 1: Write the failing test**

In `dev/test_bst_atom.m`, before the final SUITE line, add T3 plus a local kernel-finder (copied from the atom suite so this file is self-contained):

```matlab
    % ---------- T3: read a real refphase-detected band group through bst_atom ----------
    [linkFile, relData] = i_find_kernel_atom();
    if isempty(linkFile)
        fprintf('T3: SKIPPED (no unconstrained kernel link)\n');
        fprintf('\n==== SUITE: %s ====\n', PF{pass+1});  return;
    end
    % build a band-window group exactly as OnDetect does: refphase on alpha
    DataMat = in_bst_data(relData, 'F', 'Time');
    ChannelMat = in_bst_channel(bst_get('ChannelFileForStudy', relData));
    iMEG = channel_find(ChannelMat.Channel, 'MEG');
    if isstruct(DataMat.F)
        IO = db_template('ImportOptions');  IO.ImportMode='Time'; IO.UseCtfComp=1; IO.UseSsp=1;
        IO.EventsMode='ignore'; IO.DisplayMessages=0; IO.RemoveBaseline='no';
        [F, TimeVector] = in_fread(DataMat.F, ChannelMat, 1, [], iMEG, IO);
    else
        F = DataMat.F(iMEG,:);  TimeVector = DataMat.Time;
    end
    OPTIONS = process_evt_refphase('Compute');  OPTIONS.freqRange = [8 13];
    [evt, ~] = process_evt_refphase('Compute', F, TimeVector, OPTIONS);
    W = bst_dynamics('NewGroup', 'alpha (8-13 Hz)');
    W.times = evt;  W.band = [8 13];  W.bandName = 'alpha';
    % freq axis: center 10.5, extent 2.5, label alpha
    lf = bst_atom('Get', W, 'freq');
    okF = abs(lf.center-10.5)<1e-9 && abs(lf.extent-2.5)<1e-9 && strcmp(lf.label,'alpha') && strcmp(lf.state,'window');
    % time axis: occ 1 is an extended window with center = mean(onset,offset), extent>0
    lt = bst_atom('Get', W, 'time', 1);
    okT = strcmp(lt.state,'window') && (lt.extent>0) && abs(lt.center-mean(evt(:,1)))<1e-9;
    % source axis: a band-window group has no source -> unlocalized
    ls = bst_atom('Get', W, 'source', 1);
    okS = strcmp(ls.state,'unlocalized');
    ok3 = (size(evt,2)>0) && okF && okT && okS;
    fprintf('T3 real detect: nWin=%d freq=%d time=%d srcUnloc=%d => %s\n', size(evt,2), okF, okT, okS, PF{ok3+1});
    pass = pass && ok3;
```

Add the kernel-finder helper at the end of `dev/test_bst_atom.m`:

```matlab
function [linkFile, relData] = i_find_kernel_atom()
    linkFile = '';
    relData = 'Subject01/S01_AEF_20131218_01_notch/data_block001_02.mat';
    [sStudy, ~] = bst_get('DataFile', relData);
    if isempty(sStudy), return; end
    comments = {sStudy.Result.Comment};  fnames = {sStudy.Result.FileName};
    isMN = ~cellfun(@isempty, regexp(comments, 'MN: MEG\(Unconstr\)', 'once')) & ...
           ~cellfun(@isempty, regexp(fnames,   'KERNEL', 'once'));
    for j = find(isMN)
        try
            r = in_bst_results(fnames{j}, 0, 'nComponents','ImagingKernel');
            if (r.nComponents==3) && ~isempty(r.ImagingKernel)
                linkFile = ['link|' fnames{j} '|' relData];  return;
            end
        catch
        end
    end
end
```

- [ ] **Step 2: Run test to verify it fails**

First confirm it is a genuine assertion (not vacuous): temporarily it will run only if the kernel exists. Run:
```matlab
rehash; test_bst_atom
```
Expected: with the fixture present, T3 runs; before any code change it already has `bst_atom('Get')` from Tasks 1–2, so T3 should PASS immediately (it adds no new production code — it is an integration assertion). If the fixture is absent it prints `T3: SKIPPED`. Either way T1/T2 must still PASS.

(If T3 fails on `okT`/`okF`, that is a real accessor bug to fix in `bst_atom.m`, not a test artifact.)

- [ ] **Step 3: No production code needed — T3 is integration-only**

T3 exercises `bst_atom('Get')` against real `process_evt_refphase` output. If it passes, proceed. If it fails, fix the mapping in `bst_atom.m` `Get` (re-run Task 1/2 tests after any fix).

- [ ] **Step 4: Run the full accessor suite + the atom suite**

Run:
```matlab
rehash; test_bst_atom
test_dynamics_atoms
```
Expected: `test_bst_atom` ends `==== SUITE: PASS ====` (T1, T2, T3); `test_dynamics_atoms` stays `==== SUITE: PASS ====` (8/8).

- [ ] **Step 5: Commit**

```bash
git add dev/test_bst_atom.m
git commit -m "test(atom): integration — read real refphase detect output via bst_atom Get"
```

---

## Self-Review

**1. Spec coverage (Phase 1 scope of the analysis doc):**
- Localization `(center, extent, weighting)` primitive → `NewLoc` + `Get`/`Set` (Tasks 1–2); `weighting='hard'` default reserved (Global Constraints).
- Three-state axis (unlocalized/point/window) → `i_state` + tests for all three (T1 scale=unlocalized, src2=point, time/freq/src1=window).
- Per-occurrence vs group-level axes → `Axes().perOcc`; `Get`/`Set` ignore `occ` for freq/scale (T1 `lf=Get(G,'freq')`, T2 `Set(...,'freq',[],...)`).
- Source extent = geodesic radius → new `radius` field (Task 1) + source Get/Set (Tasks 1–2).
- Non-breaking → only additive `radius`; `test_dynamics_atoms` re-run guards it (Task 1 Step 6).
- Optional label (atlas-ready) → `loc.label` from `bandName`/`scaleName` (T1/T2/T3).
- Reads real detector output → T3 integration.
- Out of Phase 1 (correctly deferred to later plans): panel navigator, Navigate/Detect/Save separation, scale-axis detector, atlas layer, soft/wavelet weighting.

**2. Placeholder scan:** none — every code step is complete; every run step has an exact command + expected line.

**3. Type consistency:** the Localization struct fields (`axis,center,extent,weighting,label,state,pos`) are identical across `NewLoc`, `Get`, `Set`, and all test helpers. `bst_atom` verbs (`Axes`,`NewLoc`,`Get`,`Set`) match between definition and call sites. `radius` is `[1xN]` everywhere. `Axes().perOcc` order (time,freq,source,scale = true,false,true,false) matches the T1 assertion. `i_loc`/`i_loc_lbl` test helpers build the same struct `Get` returns.
