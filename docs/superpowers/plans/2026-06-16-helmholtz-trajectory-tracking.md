# Helmholtz Per-Hemisphere + Live Trajectory Tracking — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Make `view_helmholtz` count + gate singular points per hemisphere, and add a "Track trajectory" toggle that draws geodesic trajectories (accumulated over contiguous play) with per-track Karcher-mean centroids.

**Architecture:** Two pure, testable helpers (`bst_persistence_gate`, `bst_vortex_link_step`) plus GUI state in `view_helmholtz`. Cores gain a `.hemi` tag. The Track overlay accumulates as the time cursor steps forward by one: each frame's gated cores are linked to the previous frame's track heads (same chirality+hemisphere, greedy nearest), connected by mesh-edge Dijkstra geodesics; per-hemisphere nxr contexts give Karcher-mean centroids.

**Tech Stack:** MATLAB (Brainstorm); nxr-compute (`nxr.manifold.context`, `nxr.manifold.query.center`); MATLAB `graph`/`shortestpath`. Tests: plain functions with `chk(label,cond)` run via the MATLAB MCP.

**Commits:** user-managed on `development`; commit steps OPTIONAL/user-gated.

**Preconditions:** MATLAB up, Brainstorm running, `TutorialAuditory` protocol, `addpath(genpath('<root>/toolbox'))` + `dev/tests`.

---

### Task 1: Tag cores with `.hemi`

**Files:** Modify `toolbox/math/bst_dirac_helmholtz.m`; Test `dev/tests/test_dirac_helmholtz.m`.

- [ ] **Step 1: Failing test** — append before the final `fprintf` in `test_dirac_helmholtz.m`:

```matlab
    % --- (8) cores carry a hemisphere tag consistent with Op.vH ---
    Op8 = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    Ht8 = bst_dirac_helmholtz('Frame', Op8, J(:,1));
    nFail = nFail + chk('cores have hemi field', isfield(Ht8.Cores,'hemi'));
    okHemi = true;
    for c = Ht8.Cores
        okHemi = okHemi && ismember(c.iVertex, Op8.vH{c.hemi});
    end
    nFail = nFail + chk('hemi matches Op.vH membership', okHemi);
```

- [ ] **Step 2: Run, expect FAIL** (`Reference to non-existent field 'hemi'`). Run `dev/tests/test_dirac_helmholtz.m`.

- [ ] **Step 3: Implement** — in `bst_dirac_helmholtz.m`:

In `i_empty_cores`, add `'hemi'`:
```matlab
function s = i_empty_cores()
    s = struct('iVertex',{},'charge',{},'chirality',{},'omega',{},'persistence',{}, ...
               'isGlobal',{},'birth',{},'death',{},'pos',{},'hemi',{});
end
```
Change `i_make_core` to take and store `hemi`:
```matlab
function s = i_make_core(vg, om, chirality, persistence, isGlobal, birth, death, pos, hemi)
    if om ~= 0, ch = sign(om); else, ch = -chirality; end
    s = struct('iVertex',vg, 'charge',ch, 'chirality',chirality, 'omega',om, ...
               'persistence',persistence, 'isGlobal',logical(isGlobal), ...
               'birth',birth, 'death',death, 'pos',pos, 'hemi',hemi);
end
```
In `FindCoresOp`, pass `hh` (the loop index is the hemisphere):
```matlab
            cores(end+1) = i_make_core(vg, omega(vg), C.chirality(k), C.persistence(k), ...
                C.isGlobal(k), C.birth(k), C.death(k), ...
                i_subvertex(vloc, fl, nb, Vloc, Vn), hh); %#ok<AGROW>
```
In legacy `FindCores`, pass `1`:
```matlab
        cores(end+1) = i_make_core(v, omega(v), C.chirality(k), C.persistence(k), ...
            C.isGlobal(k), C.birth(k), C.death(k), nan(1,3), 1); %#ok<AGROW>
```

- [ ] **Step 4: Run, expect PASS** (`test_dirac_helmholtz: 0 failed`, now 25 checks). Also re-run `dev/tests/test_vortex_track.m` and `test_vortex_persistence.m` (should still pass; `process_vortex_track` ignores `.hemi`).

- [ ] **Step 5: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): tag cores with hemisphere index"`

---

### Task 2: Per-hemisphere persistence gate + readout

**Files:** Create `toolbox/math/bst_persistence_gate.m`; Modify `toolbox/gui/view_helmholtz.m`; Test `dev/tests/test_helmholtz_track.m`.

- [ ] **Step 1: Failing test** — create `dev/tests/test_helmholtz_track.m`:

```matlab
function test_helmholtz_track()
% Unit tests for the pure helpers behind the Helmholtz trajectory overlay.
% Author: Diellor Basha, 2026
    nFail = 0;
    mkc = @(v,ch,p,h,xyz) struct('iVertex',v,'charge',ch,'chirality',sign(ch), ...
                                 'persistence',p,'hemi',h,'pos',xyz);

    % --- bst_persistence_gate: per-hemisphere threshold ---
    % hemi1 max finite=5; hemi2 max finite=0.5. A global gate (max=5) at frac .5 -> thr 2.5
    % would drop hemi2's real core (0.5); per-hemi keeps it.
    mk = [ mkc(1, 1, inf, 1,[0 0 0]), mkc(2, 1, 5,   1,[0 0 0]), ...
           mkc(3,-1, inf, 2,[0 0 0]), mkc(4,-1, 0.5, 2,[0 0 0]), mkc(5,-1,0.05,2,[0 0 0]) ];
    g = bst_persistence_gate(mk, 0.5);
    kv = [g.iVertex];
    nFail = nFail + chk('globals kept (1,3)', all(ismember([1 3], kv)));
    nFail = nFail + chk('hemi1 core 2 kept (>=2.5)', ismember(2, kv));
    nFail = nFail + chk('hemi2 core 4 kept per-hemi (>=0.25)', ismember(4, kv));
    nFail = nFail + chk('hemi2 core 5 dropped (<0.25)', ~ismember(5, kv));

    % --- bst_vortex_link_step: chirality + hemi gated greedy match ---
    prev = [ mkc(10, 1, 1, 1,[0 0 0]), mkc(11,-1,1,1,[0.10 0 0]) ];
    cur  = [ mkc(20, 1, 1, 1,[0.001 0 0]), ...   % matches prev(1)
             mkc(21,-1, 1, 2,[0.10 0 0]), ...     % hemi mismatch -> no match to prev(2)
             mkc(22, 1, 1, 1,[0.20 0 0]) ];       % too far -> birth
    m = bst_vortex_link_step(prev, cur, 0.012);
    nFail = nFail + chk('prev1 -> cur1', m(1)==1);
    nFail = nFail + chk('prev2 unmatched (hemi)', m(2)==0);
    nFail = nFail + chk('chirality mismatch blocks', bst_vortex_link_step(mkc(1,1,1,1,[0 0 0]), mkc(2,-1,1,1,[0 0 0]), 0.012)==0);
    nFail = nFail + chk('within radius same chir+hemi matches', bst_vortex_link_step(mkc(1,1,1,1,[0 0 0]), mkc(2,1,1,1,[0.005 0 0]), 0.012)==1);

    fprintf('\n==== test_helmholtz_track: %d failed ====\n', nFail);
    if nFail > 0, error('test_helmholtz_track FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run, expect FAIL** (`Undefined function 'bst_persistence_gate'`). Run `dev/tests/test_helmholtz_track.m`.

- [ ] **Step 3: Implement `bst_persistence_gate.m`**:

```matlab
function mk = bst_persistence_gate(mk, frac)
% BST_PERSISTENCE_GATE  Per-hemisphere persistence gate for singular-point markers.
% Keeps a core if it is global (Inf persistence) or its persistence >= frac * (max
% FINITE persistence within its OWN hemisphere). frac<=0 keeps everything.
%   mk : struct array with fields .hemi and .persistence
% Author: Diellor Basha, 2026
    if isempty(mk) || frac <= 0, return; end
    hh   = [mk.hemi];
    keep = false(1, numel(mk));
    for h = unique(hh)
        idx = find(hh == h);
        pr  = [mk(idx).persistence];
        mxf = max([pr(isfinite(pr)), 0]);
        if mxf <= 0
            keep(idx) = true;
        else
            keep(idx) = isinf(pr) | (pr >= frac * mxf);
        end
    end
    mk = mk(keep);
end
```

(`bst_vortex_link_step` is implemented in Task 3; this test file references it but Task 2 only needs the gate checks to pass — run just the gate portion now, or implement both before running. To keep steps clean, implement `bst_vortex_link_step` here too:)

Create `toolbox/math/bst_vortex_link_step.m`:
```matlab
function match = bst_vortex_link_step(prevHeads, curCores, maxJump)
% BST_VORTEX_LINK_STEP  One-frame greedy match of track heads to current cores.
% match(h) = index into curCores matched to prevHeads(h) (0 = none). A pair is
% eligible only if same chirality SIGN and same hemisphere and within maxJump (m).
%   prevHeads, curCores : struct arrays with .pos (1x3) .charge .hemi
% Author: Diellor Basha, 2026
    nH = numel(prevHeads);  nC = numel(curCores);
    match = zeros(1, nH);
    if nH == 0 || nC == 0, return; end
    pairs = zeros(0,3);                       % [h, c, dist]
    for h = 1:nH
        for c = 1:nC
            if sign(prevHeads(h).charge) ~= sign(curCores(c).charge), continue; end
            if prevHeads(h).hemi ~= curCores(c).hemi, continue; end
            d = norm(prevHeads(h).pos - curCores(c).pos);
            if d <= maxJump, pairs(end+1,:) = [h c d]; end %#ok<AGROW>
        end
    end
    if isempty(pairs), return; end
    pairs = sortrows(pairs, 3);
    usedC = false(1, nC);
    for p = 1:size(pairs,1)
        h = pairs(p,1);  c = pairs(p,2);
        if match(h) ~= 0 || usedC(c), continue; end
        match(h) = c;  usedC(c) = true;
    end
end
```

- [ ] **Step 4: Run, expect PASS** (`test_helmholtz_track: 0 failed`).

- [ ] **Step 5: Wire the gate + per-hemi readout into `view_helmholtz.m`.**

Replace the gate block (currently lines ~126-135) with:
```matlab
    % --- component markers, pruned per hemisphere by the persistence gate;
    %     each hemisphere's global core is always kept ---
    mk = comp.Markers;
    if ~isempty(mk) && St.GateFrac > 0
        mk = bst_persistence_gate(mk, St.GateFrac);
    end
```
Replace the `i_readout` `'vortex'` and `'source'` cases with per-hemisphere lines:
```matlab
        case 'vortex'
            txt = i_count_str(mk, 'vortices', 'antivortices', Ht);
        case 'source'
            txt = i_count_str(mk, 'sources', 'sinks', []);
```
Add this helper at the bottom of `view_helmholtz.m` (before the final `end` of the file, alongside other locals):
```matlab
function txt = i_count_str(mk, posName, negName, Ht)
    if isempty(mk)
        txt = sprintf('0 %s, 0 %s', posName, negName);
    else
        hh = [mk.hemi];  parts = {};
        for h = unique(hh)
            m = mk(hh==h);
            np = sum([m.charge] > 0);  nn = sum([m.charge] < 0);
            tag = sprintf('H%d', h);
            parts{end+1} = sprintf('%s: %d(+), %d(-)', tag, np, nn); %#ok<AGROW>
        end
        txt = strjoin(parts, '  |  ');
        if ~isempty(Ht)
            pr = [mk.persistence];  tp = max(pr(isfinite(pr)));  if isempty(tp), tp = 0; end
            txt = sprintf('%s  (top persistence %.2g)', txt, tp);
        end
    end
end
```

- [ ] **Step 6: Lint** both new files + `view_helmholtz.m` (MCP `check_matlab_code`); idioms only.

- [ ] **Step 7: Verify readout/gate on a real frame** (MCP `evaluate_matlab_code`):
```matlab
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat'; [sS,~]=bst_get('DataFile',df);
link=['link|' sS.Result(1).FileName '|' df]; sR=in_bst_results(link,0); DMx=in_bst_data(sR.DataFile);
[~,iT]=min(abs(DMx.Time-22.6)); Jt=sR.ImagingKernel*DMx.F(sR.GoodChannel,iT);
Sf=in_tess_bst(sR.SurfaceFile,0); D=tess_operators(sR.SurfaceFile,'Dirac'); L=tess_operators(sR.SurfaceFile,'Laplace-Beltrami');
Op=bst_dirac_helmholtz('Prepare',D,L,Sf); Ht=bst_dirac_helmholtz('Frame',Op,Jt);
g=bst_persistence_gate(Ht.Cores,0.5); hh=[g.hemi];
fprintf('per-hemi kept: H1=%d, H2=%d (total %d of %d)\n', sum(hh==1),sum(hh==2),numel(g),numel(Ht.Cores));
```
Expected: nonzero kept in BOTH hemispheres (per-hemi gate does not starve the weaker side).

- [ ] **Step 8: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): per-hemisphere persistence gate + readout"`

---

### Task 3: Panel toggle + SetTrack + state scaffolding

**Files:** Modify `toolbox/gui/panel_helmholtz.m`, `toolbox/gui/view_helmholtz.m`.

- [ ] **Step 1: Add the checkbox to `panel_helmholtz.m`.** After the readout line (~line 52), add:
```matlab
    jTrack = gui_component('checkbox', jSec, 'br', 'Track trajectory');
    java_setcb(jTrack, 'ActionPerformedCallback', @(h,e) OnTrack(panelName));
```
Add `jTrack` to the controls struct (line ~57-59):
```matlab
    ctrl = struct('hFig',hFig, 'jVec',jVec, 'jMark',jMark, 'jReadout',jReadout, ...
                  'jSmoothOn',jSmoothOn, 'Lambda',Lambda, 'jThresh',jThresh, 'jTrack',jTrack);
```
Add the callback (near `OnMarkers`):
```matlab
function OnTrack(panelName) %#ok<DEFNU>
    ctrl = bst_get('PanelControls', panelName);
    if i_valid(ctrl); view_helmholtz('SetTrack', ctrl.hFig, ctrl.jTrack.isSelected()); end
end
```

- [ ] **Step 2: Register `SetTrack` in `view_helmholtz` dispatch.** In the top dispatch (lines ~23-30), add `'SetTrack'` to BOTH the outer command list and the handle-requiring inner list:
```matlab
    if (nargin >= 1) && ischar(SrcResultsFile) && any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','Close','UpdateFrame'}))
        if any(strcmp(SrcResultsFile, {'SetComponent','SetVectors','SetMarkers','SetSmoothing','SetGate','SetTrack','UpdateFrame'})) && ...
                (isempty(varargin) || isempty(varargin{1}) || ~all(ishandle(varargin{1})))
            return;
        end
        feval(SrcResultsFile, varargin{:});
        return;
    end
```

- [ ] **Step 3: Extend the state struct.** In the `St = struct(...)` constructor (lines ~70-74), add the track fields and stash geometry (Surf is already loaded at line ~51):
```matlab
    St = struct('Op',Op, 'srcDS',iDSf, 'srcResult',iResult, 'Component','Total', ...
                'ShowVectors',true, 'ShowMarkers',true, 'iTess',iTess, 'nV',nV, ...
                'EigenMat',EigenMat, 'Mass',{OpMat.Mass}, 'Lambda',Lambda, ...
                'Smooth',struct('on',false,'name','heat','params',struct()), 'GateFrac',0, ...
                'Cache',containers.Map('KeyType','double','ValueType','any'), ...
                'Track',false, 'LastIT',[], 'Tracks',[], 'Graph',[], 'Ctx',{cell(1,numel(Op.vH))}, ...
                'V2H',[], 'Vertices',Surf.Vertices, 'Faces',double(Surf.Faces), ...
                'TrackComp','', 'TrackSmooth',false);
```

- [ ] **Step 4: Add the `SetTrack` handler** (near `SetMarkers`):
```matlab
function SetTrack(hFig, isOn) %#ok<DEFNU>
    St = getappdata(hFig, 'HelmholtzState'); if isempty(St); return; end
    St.Track = logical(isOn);
    if ~St.Track
        hAx = findobj(hFig,'-depth',1,'Tag','Axes3D');
        if ~isempty(hAx)
            delete(findobj(hAx(1),'Tag','HelmholtzTrack'));
            delete(findobj(hAx(1),'Tag','HelmholtzCentroid'));
        end
        St.Tracks = []; St.LastIT = [];
    end
    setappdata(hFig, 'HelmholtzState', St);
    UpdateFrame(hFig);
end
```

- [ ] **Step 5: Verify panel + dispatch build** (MCP `evaluate_matlab_code`):
```matlab
bp = panel_helmholtz('CreatePanel', [], (1:50)'); disp(class(bp));
view_helmholtz('SetTrack', [], 1);   % no-op (no handle) -> must not error
fprintf('panel + SetTrack dispatch OK\n');
```
Expected: `BstPanel`, then `panel + SetTrack dispatch OK` (the no-handle SetTrack returns silently).

- [ ] **Step 6: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): Track-trajectory toggle scaffolding"`

---

### Task 4: Trajectory accumulation + geodesic drawing

**Files:** Modify `toolbox/gui/view_helmholtz.m`.

- [ ] **Step 1: Call the updater from `UpdateFrame`.** Immediately after the `i_readout(...)` line at the end of `UpdateFrame`, add:
```matlab
    if St.Track
        St = i_track_update(hFig, St, hAx, mk, iT);
        setappdata(hFig, 'HelmholtzState', St);
    end
```

- [ ] **Step 2: Implement the updater + helpers** (add as locals in `view_helmholtz.m`):
```matlab
function St = i_track_update(hFig, St, hAx, mk, iT) %#ok<INUSL>
    MAXJUMP = 0.012;                                  % m, max core jump per frame
    V = St.Vertices;
    if isempty(St.Graph),  St.Graph = i_build_graph(V, St.Faces);  end
    if isempty(St.V2H)
        St.V2H = zeros(size(V,1),1);
        for h = 1:numel(St.Op.vH), St.V2H(St.Op.vH{h}) = h; end
    end
    % reset on toggle-on, non-contiguous time, backward step, or component/smoothing change
    reset = isempty(St.LastIT) || (iT ~= St.LastIT + 1) ...
            || ~strcmp(St.TrackComp, St.Component) || (St.TrackSmooth ~= St.Smooth.on);
    if reset
        St.Tracks = i_seed_tracks(mk, St.V2H);
    else
        heads = i_heads(St.Tracks, V);
        match = bst_vortex_link_step(heads.s, i_cores_struct(mk, St.V2H, V), MAXJUMP);
        cur   = i_cores_struct(mk, St.V2H, V);
        usedC = false(1, numel(cur));
        for hi = 1:numel(match)
            ti = heads.idx(hi);
            if match(hi) > 0
                c = match(hi);  usedC(c) = true;
                p = shortestpath(St.Graph, St.Tracks(ti).coreVerts(end), cur(c).iVertex);
                if numel(p) >= 2, St.Tracks(ti).path = [St.Tracks(ti).path, p(2:end)]; end
                St.Tracks(ti).coreVerts(end+1) = cur(c).iVertex;
                St.Tracks(ti).persist = cur(c).persistence;
            else
                St.Tracks(ti).open = false;
            end
        end
        born = i_seed_tracks(mk(~usedC), St.V2H);
        St.Tracks = [St.Tracks, born];
    end
    % draw polylines
    delete(findobj(hAx,'Tag','HelmholtzTrack'));
    for t = St.Tracks
        if numel(t.path) < 2, continue; end
        col = [1 0 0]; if t.chirality < 0, col = [0 0 1]; end
        line('Parent',hAx,'XData',V(t.path,1),'YData',V(t.path,2),'ZData',V(t.path,3), ...
             'Color',col,'LineWidth',2,'Tag','HelmholtzTrack','Clipping','off');
    end
    St.LastIT = iT;  St.TrackComp = St.Component;  St.TrackSmooth = St.Smooth.on;
end

function G = i_build_graph(V, F)
    E = [F(:,[1 2]); F(:,[2 3]); F(:,[3 1])];
    E = unique(sort(E,2), 'rows');
    w = sqrt(sum((V(E(:,1),:) - V(E(:,2),:)).^2, 2));
    G = graph(E(:,1), E(:,2), w, size(V,1));
end

function T = i_seed_tracks(mk, V2H)
    T = i_empty_track();
    for k = 1:numel(mk)
        v = mk(k).iVertex;
        T(end+1) = struct('coreVerts',v, 'path',v, 'chirality',sign(mk(k).charge), ...
                          'hemi',V2H(v), 'open',true, 'persist',mk(k).persistence, ...
                          'centroid',[]); %#ok<AGROW>
    end
end
function t = i_empty_track()
    t = struct('coreVerts',{},'path',{},'chirality',{},'hemi',{},'open',{},'persist',{},'centroid',{});
end
function h = i_heads(Tracks, V) %#ok<INUSD>
    h.idx = find([Tracks.open]);
    s = struct('pos',{},'charge',{},'hemi',{});
    for i = h.idx
        v = Tracks(i).coreVerts(end);
        s(end+1) = struct('pos',V(v,:), 'charge',Tracks(i).chirality, 'hemi',Tracks(i).hemi); %#ok<AGROW>
    end
    h.s = s;
end
function cur = i_cores_struct(mk, V2H, V)
    cur = struct('pos',{},'charge',{},'hemi',{},'iVertex',{},'persistence',{});
    for k = 1:numel(mk)
        v = mk(k).iVertex;
        cur(end+1) = struct('pos',V(v,:), 'charge',mk(k).charge, 'hemi',V2H(v), ...
                            'iVertex',v, 'persistence',mk(k).persistence); %#ok<AGROW>
    end
end
```

- [ ] **Step 3: Lint** `view_helmholtz.m` (idioms only).

- [ ] **Step 4: Integration — accumulate over frames headlessly.** Drive `UpdateFrame` logic directly to confirm tracks grow and stay within a hemisphere (MCP `evaluate_matlab_code`):
```matlab
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat'; [sS,~]=bst_get('DataFile',df);
link=['link|' sS.Result(1).FileName '|' df]; sR=in_bst_results(link,0); DMx=in_bst_data(sR.DataFile);
Sf=in_tess_bst(sR.SurfaceFile,0); D=tess_operators(sR.SurfaceFile,'Dirac'); L=tess_operators(sR.SurfaceFile,'Laplace-Beltrami');
Op=bst_dirac_helmholtz('Prepare',D,L,Sf); win=find(DMx.Time>=22.55 & DMx.Time<=22.65);
G=[]; F=double(Sf.Faces); V=Sf.Vertices; V2H=zeros(size(V,1),1); for h=1:numel(Op.vH),V2H(Op.vH{h})=h;end
E=unique(sort([F(:,[1 2]);F(:,[2 3]);F(:,[3 1])],2),'rows'); Gr=graph(E(:,1),E(:,2),sqrt(sum((V(E(:,1),:)-V(E(:,2),:)).^2,2)),size(V,1));
heads=[]; tracks={};
prev=[];
for n=1:numel(win)
  Ht=bst_dirac_helmholtz('Frame',Op,sR.ImagingKernel*DMx.F(sR.GoodChannel,win(n)));
  mk=bst_persistence_gate(Ht.Cores,0.6);
  if ~isempty(prev)
     cur=arrayfun(@(c)struct('pos',V(c.iVertex,:),'charge',c.charge,'hemi',V2H(c.iVertex)),mk);
     m=bst_vortex_link_step(prev,cur,0.012);
     nLinks=sum(m>0);
     if n==2, fprintf('frame2 links=%d; sample geodesic len=%d\n', nLinks, numel(shortestpath(Gr, prev(find(m>0,1)).v, cur(m(find(m>0,1))).iVertex))); end
  end
  prev=arrayfun(@(c)struct('pos',V(c.iVertex,:),'charge',c.charge,'hemi',V2H(c.iVertex),'v',c.iVertex),mk);
end
fprintf('stepped %d frames OK\n', numel(win));
```
Expected: nonzero links on frame 2 and a geodesic length >= 2 (a real path); "stepped N frames OK".

- [ ] **Step 5: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): live geodesic trajectory accumulation"`

---

### Task 5: Karcher-mean centroids

**Files:** Modify `toolbox/gui/view_helmholtz.m`.

- [ ] **Step 1: Compute + draw centroids in `i_track_update`.** Before the final `St.LastIT = iT;` line, add:
```matlab
    % --- Karcher-mean centroids for the top-N longest-lived, strongest tracks ---
    delete(findobj(hAx,'Tag','HelmholtzCentroid'));
    openT = find([St.Tracks.open] & arrayfun(@(t)numel(t.coreVerts), St.Tracks) >= 3);
    if ~isempty(openT)
        [~, ord] = sort(arrayfun(@(i)St.Tracks(i).persist, openT), 'descend');
        openT = openT(ord(1:min(5, numel(ord))));            % cap to top-5
        for i = openT
            hh = St.Tracks(i).hemi;
            [ctx, loc] = i_hemi_context(St, hh);
            lv = loc(St.Tracks(i).coreVerts);  lv = lv(lv > 0);
            if numel(lv) < 3, continue; end
            try
                c = nxr.manifold.query.center(ctx, lv(:));
            catch, continue; end
            if isscalar(c), p = St.Vertices(St.Op.vH{hh}(c), :); else, p = c(:)'; end
            col = [1 0 0]; if St.Tracks(i).chirality < 0, col = [0 0 1]; end
            line('Parent',hAx,'XData',p(1),'YData',p(2),'ZData',p(3),'Marker','d', ...
                 'MarkerSize',12,'MarkerFaceColor',col,'MarkerEdgeColor','k','LineStyle','none', ...
                 'Tag','HelmholtzCentroid','Clipping','off');
        end
    end
```
Add the cached-context helper (local):
```matlab
function [ctx, loc] = i_hemi_context(St, hh)
% Cache one nxr context per hemisphere + a global->local vertex map (persisted on hFig
% via the returned St is NOT possible here; cache lives in appdata keyed by hemi).
    key = sprintf('HelmCtx%d', hh);
    ctx = getappdata(0, key);
    if isempty(ctx)
        vH = St.Op.vH{hh};  nVh = numel(vH);
        isV = false(size(St.Vertices,1),1); isV(vH) = true;
        fMask = all(isV(St.Faces), 2);
        mapV = zeros(size(St.Vertices,1),1); mapV(vH) = 1:nVh;
        Floc = mapV(St.Faces(fMask,:));  Vloc = St.Vertices(vH,:);
        ctx = nxr.manifold.context(Vloc, int32(Floc));
        setappdata(0, key, ctx);
    end
    loc = zeros(size(St.Vertices,1),1);
    loc(St.Op.vH{hh}) = 1:numel(St.Op.vH{hh});
end
```
(Note: contexts are cached in root appdata keyed by hemisphere to survive the per-call `St`; they're rebuilt once per MATLAB session. This is acceptable for a single loaded surface; if multiple surfaces are tracked in one session, clear with `rmappdata(0,'HelmCtx1')`.)

- [ ] **Step 2: Lint** `view_helmholtz.m`.

- [ ] **Step 3: Verify Karcher centroid on a real track** (MCP `evaluate_matlab_code`):
```matlab
% reuse Op, Sf, V, V2H from a fresh build; pick one hemisphere's vortex vertices over 4 frames
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat'; [sS,~]=bst_get('DataFile',df);
link=['link|' sS.Result(1).FileName '|' df]; sR=in_bst_results(link,0); DMx=in_bst_data(sR.DataFile);
Sf=in_tess_bst(sR.SurfaceFile,0); D=tess_operators(sR.SurfaceFile,'Dirac'); L=tess_operators(sR.SurfaceFile,'Laplace-Beltrami');
Op=bst_dirac_helmholtz('Prepare',D,L,Sf); vH=Op.vH{1};
isV=false(size(Sf.Vertices,1),1); isV(vH)=true; fM=all(isV(double(Sf.Faces)),2);
mapV=zeros(size(Sf.Vertices,1),1); mapV(vH)=1:numel(vH);
Floc=mapV(double(Sf.Faces(fM,:))); Vloc=Sf.Vertices(vH,:);
ctx=nxr.manifold.context(Vloc,int32(Floc));
lv=mapV([vH(100); vH(200); vH(300)]);
c=nxr.manifold.query.center(ctx, lv(:));
fprintf('Karcher center returned: %s (class %s)\n', mat2str(size(c)), class(c));
```
Expected: returns a scalar vertex index or a 3x1 coord without error (proves the per-hemisphere context + center call work).

- [ ] **Step 4: Commit (optional)** `git add -A && git commit -m "feat(helmholtz): per-trajectory Karcher-mean centroids"`

---

### Task 6: Full GUI integration + regression

**Files:** none (verification only)

- [ ] **Step 1: Run all suites** — `test_helmholtz_track.m`, `test_dirac_helmholtz.m`, `test_vortex_track.m`, `test_vortex_persistence.m`: all `0 failed`.

- [ ] **Step 2: Open the live view, enable Track, step frames, screenshot.** (MCP `evaluate_matlab_code`)
```matlab
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat'; [sS,~]=bst_get('DataFile',df);
link=['link|' sS.Result(1).FileName '|' df];
hFig = view_helmholtz(link);
view_helmholtz('SetComponent', hFig, 'Solen');
view_helmholtz('SetGate', hFig, 0.5);
view_helmholtz('SetTrack', hFig, 1);
% step the global time cursor across a few contiguous frames
panel_time('SetCurrentTime', 22.58);
for tt = 22.58:0.005:22.62, panel_time('SetCurrentTime', tt); end
png='/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks/helmholtz_track_view.png';
saveas(hFig, png); fprintf('saved %s\n', png);
```
Expected: the figure shows per-hemisphere singular points (persistence-sized), trajectory polylines accumulating within hemispheres, and diamond centroids; the panel readout splits H1/H2. Inspect the PNG.

(If `view_helmholtz` triggers a long `tess_eigen('Dirac')` build on first open, allow it to finish; it is cached afterward.)

- [ ] **Step 3: Final lint sweep** of the three modified/created source files.

- [ ] **Step 4: Commit (optional)** `git add -A && git commit -m "test(helmholtz): trajectory tracking integration"`

---

## Self-review notes

- **Spec coverage:** A (Task 1 hemi tag, Task 2 gate+readout) ✓; B (Task 3 toggle/scaffold, Task 4 accumulate+geodesic) ✓; C (Task 5 Karcher centroid, capped top-5, per-hemi context) ✓; reset semantics (Task 4 `reset` predicate) ✓; tests (Task 2 gate+link unit, Task 1 hemi, Task 6 integration) ✓.
- **Naming consistency:** `bst_persistence_gate`, `bst_vortex_link_step` used identically in tests and `i_track_update`; track struct fields `coreVerts/path/chirality/hemi/open/persist/centroid` consistent across `i_seed_tracks`/`i_heads`/`i_track_update`/centroid block; `St.V2H`, `St.Graph`, `St.Ctx`, `St.Vertices`, `St.Faces` defined in Task 3 and consumed in Tasks 4-5.
- **Placeholder scan:** none.
- **Risk:** Karcher per-frame solve cost — bounded by top-5/len>=3 cap; contexts cached in root appdata. Geodesic uses `shortestpath` on a cached `graph` (hemispheres are disconnected components, so a match never crosses hemispheres — `bst_vortex_link_step` already enforces same-hemi).
```
