# Vortex Tracking (Phase 2) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** A batch process that links per-frame Dirac vortex cores (and sources/sinks) into trajectories over a time window and saves them as a Brainstorm dipoles file, viewable in `view_dipoles`.

**Architecture:** A pure greedy chirality-consistent linker (`bst_vortex_track`) turns per-frame core lists into trajectories. `process_vortex_track` decomposes a Dirac source over a window (Phase 1 detection per frame), links vortices and sources, and assembles one dipoles file (each track = one dipole `.Index` group; `.Loc` = sub-vertex pos; `.Amplitude` = chirality·persistence·normal; `.Goodness` = persistence). Reuses `view_dipoles`/`panel_dipoles` — no new viewer.

**Tech Stack:** MATLAB (Brainstorm). Tests are plain functions with a `chk(label,cond)` failure counter, run via the MATLAB MCP.

**Commits:** user-managed on `development`; commit steps are OPTIONAL/user-gated.

**Preconditions:** MATLAB up, Brainstorm running, `TutorialAuditory` protocol, toolbox on path (`addpath(genpath('<root>/toolbox'))` + `dev/tests`). Phase 1 already merged.

---

### Task 1: `bst_vortex_track` — greedy chirality-consistent linker

**Files:**
- Create: `toolbox/math/bst_vortex_track.m`
- Test: `dev/tests/test_vortex_track.m`

- [ ] **Step 1: Write the failing test**

Create `dev/tests/test_vortex_track.m`:

```matlab
function test_vortex_track()
% Pure unit tests for bst_vortex_track on hand-built core sequences.
% Author: Diellor Basha, 2026
    nFail = 0;
    mk = @(v,ch,p,xyz) struct('iVertex',v,'chirality',ch,'persistence',p,'pos',xyz);

    % A) one core drifting in small steps over 5 frames -> one length-5 track
    A = cell(1,5);
    for t=1:5, A{t} = mk(t, 1, 1, [0.001*t 0 0]); end
    T = bst_vortex_track(A, 'MaxJump', 0.010);
    nFail = nFail + chk('A: single track',        numel(T)==1);
    nFail = nFail + chk('A: length 5',            numel(T(1).frames)==5);
    nFail = nFail + chk('A: birth1 death5',       T(1).birthFrame==1 && T(1).deathFrame==5);

    % B) two separated stationary cores -> two length-5 tracks
    B = cell(1,5);
    for t=1:5, B{t} = [mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0])]; end
    Tb = bst_vortex_track(B, 'MaxJump', 0.010);
    nFail = nFail + chk('B: two tracks',          numel(Tb)==2);
    nFail = nFail + chk('B: both length 5',       all(arrayfun(@(x)numel(x.frames)==5, Tb)));

    % C) birth + death: core A (frames 1-2), core B far away (frames 3-5)
    C = {mk(1,1,1,[0 0 0]), mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0]), mk(2,1,1,[0.05 0 0]), mk(2,1,1,[0.05 0 0])};
    Tc = bst_vortex_track(C, 'MaxJump', 0.010);
    nFail = nFail + chk('C: two tracks',          numel(Tc)==2);
    db = sort([Tc.deathFrame]); bb = sort([Tc.birthFrame]);
    nFail = nFail + chk('C: a death at 2',        any([Tc.deathFrame]==2));
    nFail = nFail + chk('C: a birth at 3',        any([Tc.birthFrame]==3));

    % D) chirality mismatch must not link
    D = {mk(1,1,1,[0 0 0]), mk(2,-1,1,[0.001 0 0])};
    Td = bst_vortex_track(D, 'MaxJump', 0.010);
    nFail = nFail + chk('D: chirality blocks link', numel(Td)==2);

    % E) jump beyond MaxJump must not link
    E = {mk(1,1,1,[0 0 0]), mk(2,1,1,[0.05 0 0])};
    Te = bst_vortex_track(E, 'MaxJump', 0.010);
    nFail = nFail + chk('E: long jump splits',    numel(Te)==2);

    fprintf('\n==== test_vortex_track: %d failed ====\n', nFail);
    if nFail > 0, error('test_vortex_track FAILED'); end
end
function n = chk(label, cond)
    if cond; fprintf('  PASS %s\n', label); n = 0; else; fprintf('  FAIL %s\n', label); n = 1; end
end
```

- [ ] **Step 2: Run, verify it fails**

Run (MCP `run_matlab_file`): `dev/tests/test_vortex_track.m` -> `Undefined function 'bst_vortex_track'`.

- [ ] **Step 3: Implement `bst_vortex_track.m`**

Create `toolbox/math/bst_vortex_track.m`:

```matlab
function Tracks = bst_vortex_track(coresPerFrame, varargin)
% BST_VORTEX_TRACK  Link per-frame cores into trajectories (greedy, chirality-consistent).
%
% USAGE: Tracks = bst_vortex_track(coresPerFrame, 'MinPersistence', p, 'MaxJump', d)
% INPUT:
%   coresPerFrame : {1 x nT} cell; coresPerFrame{t} = struct array of cores with
%                   fields .iVertex .chirality(+/-1) .persistence .pos (1x3, meters)
%   'MinPersistence' (default 0)  : drop cores below this persistence (Inf always kept)
%   'MaxJump' (default 0.010 m)   : max core displacement between consecutive frames
% OUTPUT: Tracks struct array (one per trajectory):
%   .frames [1xL] .iVertex [1xL] .pos [Lx3] .persistence [1xL]
%   .chirality (+/-1) .birthFrame .deathFrame
% Author: Diellor Basha, 2026

    MinPersistence = 0;  MaxJump = 0.010;
    for k = 1:2:numel(varargin)
        switch lower(varargin{k})
            case 'minpersistence', MinPersistence = varargin{k+1};
            case 'maxjump',        MaxJump        = varargin{k+1};
            otherwise, error('bst_vortex_track: unknown option %s', varargin{k});
        end
    end

    nT = numel(coresPerFrame);
    Tracks = i_empty_tracks();
    openIdx = [];                         % indices into Tracks still open
    for t = 1:nT
        cur = coresPerFrame{t};
        if ~isempty(cur)
            cur = cur([cur.persistence] >= MinPersistence);   % Inf passes
        end
        nc = numel(cur);
        matched = false(1, nc);

        if ~isempty(openIdx) && nc > 0
            H = numel(openIdx);
            pairs = zeros(0,3);            % [headSlot, curIdx, dist]
            for hs = 1:H
                tr = Tracks(openIdx(hs));
                hp = tr.pos(end,:);  hc = tr.chirality;
                for c = 1:nc
                    if cur(c).chirality ~= hc, continue; end
                    d = norm(cur(c).pos - hp);
                    if d <= MaxJump, pairs(end+1,:) = [hs, c, d]; end %#ok<AGROW>
                end
            end
            usedH = false(1,H);
            if ~isempty(pairs)
                pairs = sortrows(pairs, 3);
                for p = 1:size(pairs,1)
                    hs = pairs(p,1);  c = pairs(p,2);
                    if usedH(hs) || matched(c), continue; end
                    usedH(hs) = true;  matched(c) = true;
                    ti = openIdx(hs);
                    Tracks(ti) = i_extend(Tracks(ti), t, cur(c));
                end
            end
            newOpen = openIdx(usedH);
            for hs = find(~usedH), Tracks(openIdx(hs)).deathFrame = t-1; end
            openIdx = newOpen;
        elseif ~isempty(openIdx)          % nc == 0: close everything open
            for hs = 1:numel(openIdx), Tracks(openIdx(hs)).deathFrame = t-1; end
            openIdx = [];
        end

        for c = 1:nc                       % births
            if ~matched(c)
                Tracks(end+1) = i_new(t, cur(c)); %#ok<AGROW>
                openIdx(end+1) = numel(Tracks); %#ok<AGROW>
            end
        end
    end
    for hs = 1:numel(openIdx), Tracks(openIdx(hs)).deathFrame = nT; end
end

function s = i_empty_tracks()
    s = struct('frames',{},'iVertex',{},'pos',{},'persistence',{}, ...
               'chirality',{},'birthFrame',{},'deathFrame',{});
end
function tr = i_new(t, c)
    tr = struct('frames',t, 'iVertex',c.iVertex, 'pos',c.pos, 'persistence',c.persistence, ...
                'chirality',c.chirality, 'birthFrame',t, 'deathFrame',t);
end
function tr = i_extend(tr, t, c)
    tr.frames(end+1)      = t;
    tr.iVertex(end+1)     = c.iVertex;
    tr.pos(end+1,:)       = c.pos;
    tr.persistence(end+1) = c.persistence;
    tr.deathFrame         = t;
end
```

- [ ] **Step 4: Run, verify pass**

Run: `dev/tests/test_vortex_track.m` -> `==== test_vortex_track: 0 failed ====`.

- [ ] **Step 5: Lint** `toolbox/math/bst_vortex_track.m` (MCP `check_matlab_code`).

- [ ] **Step 6: Commit (optional, user-gated)**

```bash
git add toolbox/math/bst_vortex_track.m dev/tests/test_vortex_track.m
git commit -m "feat(vortex): greedy chirality-consistent core tracker"
```

---

### Task 2: `Decompose` also emits per-frame sources

**Files:**
- Modify: `toolbox/math/bst_dirac_helmholtz.m` (`Decompose`, ~lines 150–161)
- Test: `dev/tests/test_dirac_helmholtz.m` (extend)

- [ ] **Step 1: Write the failing test addition**

In `test_dirac_helmholtz.m`, just before the final `fprintf('\n==== ...')`, append:

```matlab
    % --- (7) Decompose emits per-frame Sources matching Frame ---
    H7 = bst_dirac_helmholtz(Dirac, LBO, Surf, J);
    nFail = nFail + chk('Decompose has Sources cell', isfield(H7,'Sources') && iscell(H7.Sources) && numel(H7.Sources)==2);
    Op7 = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    Ht7 = bst_dirac_helmholtz('Frame', Op7, J(:,1));
    nFail = nFail + chk('Decompose Sources col1 == Frame', numel(H7.Sources{1})==numel(Ht7.Sources));
```

- [ ] **Step 2: Run, verify it fails**

Run: `dev/tests/test_dirac_helmholtz.m` -> FAIL at "Decompose has Sources cell".

- [ ] **Step 3: Add `H.Sources` to `Decompose`**

In `Decompose`, add `H.Sources = cell(1, nT);` next to `H.Cores = cell(1, nT);`, and inside the loop add `H.Sources{t}=Ht.Sources;` next to `H.Cores{t}=Ht.Cores;`. Result:

```matlab
    H.Cores = cell(1, nT);
    H.Sources = cell(1, nT);
    for t = 1:nT
        Ht = Frame(Op, J(:,t));
        H.Curl(:,t)=Ht.Curl; H.Div(:,t)=Ht.Div; H.Psi(:,t)=Ht.Psi;
        H.Phi(:,t)=Ht.Phi;   H.Fmag(:,t)=Ht.Fmag;  H.Cores{t}=Ht.Cores;  H.Sources{t}=Ht.Sources;
    end
```

- [ ] **Step 4: Run, verify pass**

Run: `dev/tests/test_dirac_helmholtz.m` -> `0 failed` (now 23 checks).

- [ ] **Step 5: Commit (optional, user-gated)**

```bash
git add toolbox/math/bst_dirac_helmholtz.m dev/tests/test_dirac_helmholtz.m
git commit -m "feat(helmholtz): Decompose emits per-frame Sources for tracking"
```

---

### Task 3: `process_vortex_track` — results -> dipoles

**Files:**
- Create: `toolbox/process/functions/process_vortex_track.m`

- [ ] **Step 1: Implement the process**

Create `toolbox/process/functions/process_vortex_track.m`:

```matlab
function varargout = process_vortex_track( varargin )
% PROCESS_VORTEX_TRACK  Track Dirac vortex cores (and sources/sinks) over time.
% Detects persistence-ranked cores per frame (bst_dirac_helmholtz) and links them
% into trajectories (bst_vortex_track), saved as a dipoles file (view_dipoles).
% @=============================================================================
% Author: Diellor Basha, 2026
    eval(macro_method);
end

%% ===== DESCRIPTION =====
function sProcess = GetDescription() %#ok<DEFNU>
    sProcess.Comment     = 'Track vortex cores (Dirac)';
    sProcess.Category    = 'File';
    sProcess.SubGroup    = 'Sources';
    sProcess.Index       = 0;
    sProcess.InputTypes  = {'results'};
    sProcess.OutputTypes = {'dipoles'};
    sProcess.nInputs     = 1;
    sProcess.nMinFiles   = 1;
    sProcess.options.timewindow.Comment = 'Time window:';
    sProcess.options.timewindow.Type    = 'timewindow';
    sProcess.options.timewindow.Value   = [];
    sProcess.options.minpers.Comment = 'Min persistence (0 = keep all): ';
    sProcess.options.minpers.Type    = 'value';
    sProcess.options.minpers.Value   = {0, '', 6};
    sProcess.options.maxjump.Comment = 'Max core jump per frame: ';
    sProcess.options.maxjump.Type    = 'value';
    sProcess.options.maxjump.Value   = {10, 'mm', 1};
    sProcess.options.tracksrc.Comment = 'Also track sources/sinks (phi)';
    sProcess.options.tracksrc.Type    = 'checkbox';
    sProcess.options.tracksrc.Value   = 1;
end

function Comment = FormatComment(sProcess) %#ok<DEFNU>
    Comment = sProcess.Comment;
end

%% ===== RUN =====
function OutputFiles = Run(sProcess, sInput) %#ok<DEFNU>
    OutputFiles = {};
    MinPers = sProcess.options.minpers.Value{1};
    MaxJump = sProcess.options.maxjump.Value{1} / 1000;     % mm -> m
    TrackSrc = sProcess.options.tracksrc.Value;
    TimeWindow = sProcess.options.timewindow.Value{1};

    % --- source file (must be unconstrained / 3-component) ---
    sRes = in_bst_results(sInput.FileName, 0);
    if isempty(sRes.nComponents) || (sRes.nComponents ~= 3)
        bst_report('Error', sProcess, sInput, 'Vortex tracking requires an unconstrained (3-component) Dirac source.');
        return;
    end
    % data over the window
    if isempty(sRes.DataFile), DataMat.Time = sRes.Time; DataMat.F = [];
    else, DataMat = in_bst_data(sRes.DataFile); end
    if isempty(sRes.ImageGridAmp)
        F = DataMat.F(sRes.GoodChannel, :);
        J = sRes.ImagingKernel * F;
    else
        J = sRes.ImageGridAmp;
    end
    Time = sRes.Time; if isempty(Time) || numel(Time)~=size(J,2), Time = DataMat.Time; end
    if ~isempty(TimeWindow)
        sb = bst_closest(TimeWindow, Time); iW = sb(1):sb(2);
    else
        iW = 1:size(J,2);
    end
    J = J(:, iW); tW = Time(iW);

    % --- operators + per-frame decomposition ---
    SurfaceFile = sRes.SurfaceFile;
    Surf = in_tess_bst(SurfaceFile, 0);
    Dirac = i_op(SurfaceFile, 'Dirac');  LBO = i_op(SurfaceFile, 'Laplace-Beltrami');
    Op = bst_dirac_helmholtz('Prepare', Dirac, LBO, Surf);
    nT = numel(tW);
    coresV = cell(1,nT); coresS = cell(1,nT);
    bst_progress('start', 'Vortex tracking', 'Decomposing frames...', 1, nT);
    for t = 1:nT
        Ht = bst_dirac_helmholtz('Frame', Op, J(:,t));
        coresV{t} = Ht.Cores;  coresS{t} = Ht.Sources;
        bst_progress('inc', 1);
    end
    bst_progress('stop');

    % --- link ---
    Tv = bst_vortex_track(coresV, 'MinPersistence', MinPers, 'MaxJump', MaxJump);
    Ts = i_empty_tracks();
    if TrackSrc, Ts = bst_vortex_track(coresS, 'MinPersistence', MinPers, 'MaxJump', MaxJump); end
    if isempty(Tv) && isempty(Ts)
        bst_report('Warning', sProcess, sInput, 'No vortex cores found in the selected window.');
    end

    % --- assemble dipoles ---
    VN = Surf.VertNormals;
    DipolesMat = i_tracks_to_dipoles(Tv, Ts, tW, VN);
    DipolesMat.Comment  = sprintf('Vortex tracks (%d+%d)', numel(Tv), numel(Ts));
    DipolesMat.DataFile = sInput.FileName;
    DipolesMat = bst_history('add', DipolesMat, 'vortextrack', ['Generated from: ' sInput.FileName]);

    % --- save + db ---
    [sStudy, iStudy] = bst_get('AnyFile', sInput.FileName);
    ProtocolInfo = bst_get('ProtocolInfo');
    DipoleFile = file_unique(bst_fullfile(ProtocolInfo.STUDIES, bst_fileparts(sStudy.FileName), 'dipoles_vortextrack.mat'));
    bst_save(DipoleFile, DipolesMat);
    Bst = db_template('Dipoles');
    Bst.FileName = file_short(DipoleFile);
    Bst.Comment  = DipolesMat.Comment;
    Bst.DataFile = sInput.FileName;
    sStudy.Dipoles(end+1) = Bst;
    bst_set('Study', iStudy, sStudy);
    panel_protocols('UpdateNode', 'Study', iStudy);
    db_save();
    OutputFiles{1} = DipoleFile;
end

%% ===== helpers =====
function s = i_empty_tracks()
    s = struct('frames',{},'iVertex',{},'pos',{},'persistence',{}, ...
               'chirality',{},'birthFrame',{},'deathFrame',{});
end

function D = i_tracks_to_dipoles(Tv, Ts, tW, VN)
% Each track -> one dipole .Index group; Loc=pos, Amplitude=chirality*normPers*normal.
    allP = [arrayfun(@(x)max(x.persistence(isfinite(x.persistence))), Tv), ...
            arrayfun(@(x)max(x.persistence(isfinite(x.persistence))), Ts)];
    maxP = max([allP(~isnan(allP)), eps]);
    D = struct(); D.Time = tW(:)'; D.Dipole = repmat(i_empty_dip(), 0, 1);
    names = {};
    gi = 0;
    [D, gi, names] = i_add_set(D, gi, names, Tv, tW, VN, maxP, 'Vortex');
    [D, gi, names] = i_add_set(D, gi, names, Ts, tW, VN, maxP, 'Source');  %#ok<ASGLU>
    D.DipoleNames = names;
    if isempty(D.Dipole), D.Subset = []; else, D.Subset = unique([D.Dipole.Index]); end
end

function d = i_empty_dip()
    d = struct('Index',0,'Time',0,'Origin',[0 0 0],'Loc',[0;0;0],'Amplitude',[0;0;0], ...
               'Goodness',0,'Errors',0,'Noise',[],'SingleError',[],'ErrorMatrix',[], ...
               'ConfVol',[],'Probability',[],'NoiseEstimate',[],'Perform',0);
end

function [D, gi, names] = i_add_set(D, gi, names, T, tW, VN, maxP, kind)
    for k = 1:numel(T)
        gi = gi + 1;  tr = T(k);
        if strcmp(kind,'Vortex'); lbl = sprintf('Vortex%s #%d', i_sgn(tr.chirality), gi);
        else;                     lbl = sprintf('%s #%d', i_srcname(tr.chirality), gi); end
        names{gi} = lbl; %#ok<AGROW>
        for j = 1:numel(tr.frames)
            d = i_empty_dip();
            d.Index = gi;  d.Time = tW(tr.frames(j));
            d.Loc = tr.pos(j,:)';
            pn = tr.persistence(j); if ~isfinite(pn), pn = maxP; end
            d.Amplitude = (tr.chirality * (pn/maxP) * VN(tr.iVertex(j),:))';
            d.Goodness = min(pn, maxP)/maxP;  d.Perform = pn;
            D.Dipole(end+1) = d; %#ok<AGROW>
        end
    end
end

function s = i_sgn(ch),     if ch>=0, s='+'; else, s='-'; end, end
function s = i_srcname(ch), if ch>=0, s='Source'; else, s='Sink'; end, end

function Op = i_op(SurfaceFile, variant)
    [sSubject,~,iSurf] = bst_get('SurfaceFile', SurfaceFile);
    Op = [];
    if ~isempty(iSurf) && isfield(sSubject.Surface(iSurf),'Operator')
        for k = 1:numel(sSubject.Surface(iSurf).Operator)
            S = load(file_fullpath(sSubject.Surface(iSurf).Operator(k).FileName));
            if strcmpi(S.Variant, variant), Op = S; break; end
        end
    end
    if isempty(Op), tess_operators(SurfaceFile, variant); Op = i_op(SurfaceFile, variant); end
end
```

- [ ] **Step 2: Lint** `toolbox/process/functions/process_vortex_track.m` (MCP `check_matlab_code`); resolve beyond Brainstorm idioms.

- [ ] **Step 3: Smoke-check the process registers**

Run (MCP `evaluate_matlab_code`):
```matlab
sP = process_vortex_track('GetDescription');
fprintf('Comment=%s In=%s Out=%s opts=%s\n', sP.Comment, sP.InputTypes{1}, sP.OutputTypes{1}, strjoin(fieldnames(sP.options)',','));
```
Expected: prints the comment, `results`/`dipoles`, and the option names (timewindow, minpers, maxjump, tracksrc).

- [ ] **Step 4: Commit (optional, user-gated)**

```bash
git add toolbox/process/functions/process_vortex_track.m
git commit -m "feat(vortex): batch process tracking cores into a dipoles file"
```

---

### Task 4: Integration verification on the real S01 alpha window

**Files:** none (verification only)

- [ ] **Step 1: Run all three test suites**

Run: `test_vortex_track.m`, `test_vortex_persistence.m`, `test_dirac_helmholtz.m` -> all `0 failed`.

- [ ] **Step 2: End-to-end on a real Dirac source over an alpha sub-window**

Run (MCP `evaluate_matlab_code`) — builds an in-memory unconstrained Dirac results struct, runs the Run logic via a temporary saved results file, then tracks:
```matlab
df='Subject01/S01_AEF_20131218_01_notch/data_block001_band.mat';
[sS,~]=bst_get('DataFile',df); cm=in_bst_channel(sS.Channel(1).FileName); ty={cm.Channel.Type};
HMos=in_bst_headmodel([fileparts(df) '/headmodel_surf_os_meg.mat'],0);
G=double(HMos.Gain); iMEG=all(isfinite(G),2)&strcmpi(ty(:),'MEG');
NC=load(file_fullpath([fileparts(df) '/noisecov_full.mat'])); Cn=NC.NoiseCov(iMEG,iMEG); Cn=(Cn+Cn')/2;
HMf=HMos; HMf.Gain=G(iMEG,:);
OPT=struct('NoiseMethod','reg','NoiseReg',0.1,'SnrMethod','fixed','SnrFixed',3,'InverseMeasure','dspm2018');
OPT.NoiseCovMat.NoiseCov=Cn; OPT.ChannelTypes=ty(iMEG);
Rd=bst_inverse_dirac(HMf,OPT);
DM=in_bst_data(df); win=DM.Time>=22.4 & DM.Time<=22.8; tW=DM.Time(win);
J=Rd.ImagingKernel*double(DM.F(iMEG,win));
Surf=in_tess_bst(HMos.SurfaceFile,0);
D=tess_operators(HMos.SurfaceFile,'Dirac'); L=tess_operators(HMos.SurfaceFile,'Laplace-Beltrami');
Op=bst_dirac_helmholtz('Prepare',D,L,Surf);
nT=numel(tW); cv=cell(1,nT);
for t=1:nT, Ht=bst_dirac_helmholtz('Frame',Op,J(:,t)); cv{t}=Ht.Cores; end
Tv=bst_vortex_track(cv,'MaxJump',0.012);
L4=arrayfun(@(x)numel(x.frames),Tv);
fprintf('%d frames, %d vortex tracks; longest track %d frames; max-persistence track length %d\n', ...
  nT, numel(Tv), max(L4), L4(find(arrayfun(@(x)any(isinf(x.persistence)),Tv),1)));
```
Expected: several tracks, at least one spanning most of the window (the dominant superior-parietal vortex persists), longest track length close to `nT`.

- [ ] **Step 3: Full process + viewer (manual)**

If a saved unconstrained Dirac source file exists in the DB, run the process through the pipeline (or call `process_vortex_track('Run', sProcess, sInput)` with a constructed `sProcess`/`sInput`) and confirm a `dipoles_vortextrack*.mat` node appears and `view_dipoles` opens it with the time slider animating the cores. (If no saved Dirac source exists, note this and rely on Step 2 for algorithmic verification.)

- [ ] **Step 4: Final lint sweep** of the three new/modified source files.

- [ ] **Step 5: Commit (optional, user-gated)**

```bash
git add -A
git commit -m "test(vortex): tracking regression + integration checks"
```

---

## Self-review notes

- **Spec coverage:** linker (Task 1) ✓; Decompose Sources (Task 2) ✓; process results->dipoles with both vortices+sources, dipole encoding (Task 3) ✓; greedy chirality+gating birth/death (Task 1 algorithm) ✓; tests unit+regression+integration (Tasks 1,2,4) ✓; error handling for non-3-component + empty (Task 3 Run) ✓.
- **Naming consistency:** `bst_vortex_track` returns `.frames/.iVertex/.pos/.persistence/.chirality/.birthFrame/.deathFrame`, consumed identically in `i_add_set`. `i_empty_tracks` duplicated in the process (standalone file) intentionally — the process must run without depending on a private linker helper.
- **Dipole fields** match `process_dipole_scanning` (`Index/Time/Origin/Loc/Amplitude/Goodness/Perform/...`) so `view_dipoles`/`db_template('Dipoles')` accept the file.
- **Placeholder scan:** none.
- **Risk:** Step 3 needs a saved Dirac source in the DB; Step 2 gives algorithmic verification independent of that, so the plan is verifiable regardless.
```
