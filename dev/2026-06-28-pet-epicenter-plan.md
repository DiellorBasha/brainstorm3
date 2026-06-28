# Geometry-based Epicenter Detection Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build `pet_epicenter.m` — a per-subject detector of cortical amyloid/tau concentration foci (dominant epicenter + persistence-ranked secondaries) from a surface SUVR map — plus a cohort analysis script that validates the foci against known seeding regions.

**Architecture:** Discrete Morse-Smale on the cortical mesh. Heat-smooth the SUVR field (LBO `(M+t·K)⁻¹` solve via `tess_operators`), find local maxima over the mesh adjacency (`tess_vertconn`), grow gradient-ascent basins, and rank/filter foci by topological persistence (super-level merge tree). The covariant gradient field (`bst_gradient`) is an optional flow output. No external deps beyond Brainstorm.

**Tech Stack:** MATLAB (R2023b), Brainstorm toolbox. Reuses `tess_operators` (Laplace-Beltrami `Operator`/`Mass`/`GlobalVertices`), `tess_vertconn`, `in_tess_bst`, and (for cohort mapping) `tess_interp_tess2tess`. Testing = synthetic-field-on-real-surface validation run through the MATLAB MCP tools (`check_matlab_code` for lint, `evaluate_matlab_code` for behavior).

## Global Constraints
- MATLAB only; no toolbox beyond base MATLAB + Brainstorm (no Image Processing / Stats toolbox calls — use `convn`/`corrcoef`, not `imgaussfilt`/`corr`).
- Never call `clear`; edited `.m` files auto-reload via `rehash path`.
- Surface vertices are per-hemisphere in the operators: `LBO.GlobalVertices{hh}` maps `Operator{hh}`/`Mass{hh}` (per-hemi `[nVh×nVh]`) to full-surface vertex indices.
- Files live in `toolbox/anatomy/` (function) and `dev/benchmarks/` (cohort script). Commit messages end with the two trailer lines used across this project (Co-Authored-By + Claude-Session).
- Each task: write the validation snippet, run it (fail), implement, `check_matlab_code` (lint clean), run validation (pass), commit.

---

### Task 1: `pet_epicenter` scaffold — heat-smooth + local-maxima detection

**Files:**
- Create: `toolbox/anatomy/pet_epicenter.m`

**Interfaces:**
- Produces: `[foci, basinLabel, info] = pet_epicenter(SurfaceFile, suvrMap, Opts)`.
  - `SurfaceFile` : Brainstorm cortex surface file (string).
  - `suvrMap`     : `[nVert × 1]` per-vertex SUVR (NaN allowed for the medial wall).
  - `Opts`        : `.HeatT` (heat time, default `2e-5`), `.MinPersist` (default `[]`→0.15·range), `.nFociMax` (default 10).
  - `foci`        : struct array sorted by persistence desc, fields `.vertex .peak .persistence .basinArea`.
  - `basinLabel`  : `[nVert × 1]` basin index (0 = unassigned/NaN).
  - `info`        : struct `.smoothed [nVert×1] .nFoci`.
  - This task delivers `.smoothed` + the raw local-maxima list inside `info.maxVerts`; basins/persistence are added in Tasks 2–3.

- [ ] **Step 1: Write the failing validation snippet**

Run via `evaluate_matlab_code`:

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0002');
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
sW=in_tess_bst(wf); V=sW.Vertices; seed=5000;
f=exp(-sum((V-V(seed,:)).^2,2)/(2*0.012^2));        % ~12mm Gaussian bump at vertex 5000
[~,~,info]=pet_epicenter(wf, f, struct('HeatT',2e-5));
d=sqrt(sum((V(info.maxVerts,:)-V(seed,:)).^2,1))*1000;   % mm from each detected max to the seed
assert(min(d)<6, 'dominant max not within 6mm of the planted seed');
fprintf('PASS: nearest detected max %.1f mm from seed\n', min(d));
```

- [ ] **Step 2: Run it — expect FAIL** (`Unrecognized function 'pet_epicenter'`).

- [ ] **Step 3: Implement the scaffold + heat-smooth + maxima**

Create `toolbox/anatomy/pet_epicenter.m`:

```matlab
function [foci, basinLabel, info] = pet_epicenter(SurfaceFile, suvrMap, Opts)
% PET_EPICENTER: cortical concentration foci (epicenters) of a surface SUVR map via discrete
% Morse-Smale geometry. Heat-smooth (LBO) -> local maxima (mesh adjacency) -> gradient-ascent
% basins -> persistence-ranked dominant + secondary foci.
%
% USAGE: [foci, basinLabel, info] = pet_epicenter(SurfaceFile, suvrMap, Opts)
%
% Author: Diellor Basha, 2026
    if (nargin<3)||isempty(Opts), Opts=struct(); end
    Def=struct('HeatT',2e-5,'MinPersist',[],'nFociMax',10);
    fn=fieldnames(Def); for i=1:numel(fn), if ~isfield(Opts,fn{i}), Opts.(fn{i})=Def.(fn{i}); end; end

    sSurf=in_tess_bst(SurfaceFile); nV=size(sSurf.Vertices,1);
    f=double(suvrMap(:)); f(~isfinite(f))=NaN;

    % --- LBO heat smoothing (tangential), per hemisphere ---
    fs=local_heatsmooth(SurfaceFile, f, Opts.HeatT);

    % --- vertex adjacency + local maxima (f(v) strictly >= all finite neighbours) ---
    if isfield(sSurf,'VertConn') && ~isempty(sSurf.VertConn), A=sSurf.VertConn;
    else, A=tess_vertconn(sSurf.Vertices, sSurf.Faces); end
    A=A | A'; A=A - diag(diag(A));
    maxVerts=local_localmax(fs, A);

    foci=struct('vertex',{},'peak',{},'persistence',{},'basinArea',{});
    basinLabel=zeros(nV,1);
    info=struct('smoothed',fs,'nFoci',numel(maxVerts),'maxVerts',maxVerts,'A',{A});
end

% ===== LBO heat smoothing of a vertex field (per hemisphere), cached factorization =====
function g=local_heatsmooth(SurfaceFile, f, t)
    LBO=tess_operators(SurfaceFile,'Laplace-Beltrami');
    g=f;
    for hh=1:numel(LBO.Operator)
        gv=double(LBO.GlobalVertices{hh}(:)); M=LBO.Mass{hh}; K=LBO.Operator{hh};
        fh=f(gv); fh(~isfinite(fh))=0;
        g(gv)=(M + t*K) \ (M*fh);
    end
end

% ===== local maxima over the adjacency: f(v) >= every neighbour (ties resolved by index) =====
function mv=local_localmax(f, A)
    n=numel(f); mv=[];
    for v=1:n
        if ~isfinite(f(v)), continue; end
        nb=find(A(:,v)); nb=nb(isfinite(f(nb)));
        if isempty(nb) || all(f(v) > f(nb)) || (all(f(v) >= f(nb)) && v < min(nb(f(nb)==f(v))))
            mv(end+1)=v; %#ok<AGROW>
        end
    end
    mv=mv(:);
end
```

- [ ] **Step 4: Lint** with `check_matlab_code` on `toolbox/anatomy/pet_epicenter.m` — expect `code_issues: []`.

- [ ] **Step 5: Run the Step-1 snippet — expect PASS** ("nearest detected max < 6 mm from seed").

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/pet_epicenter.m
git commit -m "feat(pet): pet_epicenter scaffold - heat-smooth + local-maxima detection

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 2: Gradient-ascent basins (Morse-Smale segmentation)

**Files:**
- Modify: `toolbox/anatomy/pet_epicenter.m`

**Interfaces:**
- Consumes: `info.maxVerts`, `info.A` (adjacency), `info.smoothed` from Task 1.
- Produces: `basinLabel [nVert×1]` (each finite vertex → index into `maxVerts` of the maximum its steepest ascent reaches; 0 for NaN vertices), and a `local_basins(fs, A, maxVerts)` helper.

- [ ] **Step 1: Write the failing validation snippet**

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0002');
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
sW=in_tess_bst(wf); V=sW.Vertices; s1=3000; s2=15000;
f=exp(-sum((V-V(s1,:)).^2,2)/(2*0.012^2)) + 0.9*exp(-sum((V-V(s2,:)).^2,2)/(2*0.012^2));
[~,basinLabel,info]=pet_epicenter(wf, f, struct('HeatT',2e-5));
nB=numel(unique(basinLabel(basinLabel>0)));
assert(nB>=2, 'expected >=2 basins for two bumps');
% the two seeds must fall in different basins
assert(basinLabel(s1)~=basinLabel(s2), 'two seeds landed in the same basin');
fprintf('PASS: %d basins; seeds separated\n', nB);
```

- [ ] **Step 2: Run it — expect FAIL** (`basinLabel` all zeros → assertion `nB>=2`).

- [ ] **Step 3: Implement basins** — in `pet_epicenter`, replace the `basinLabel=zeros(nV,1);` line with a call, and add the helper.

Replace:
```matlab
    basinLabel=zeros(nV,1);
```
with:
```matlab
    basinLabel=local_basins(fs, A, maxVerts);
```

Add this helper (after `local_localmax`):
```matlab
% ===== steepest-ascent basins: each vertex flows to the max it climbs to =====
function lbl=local_basins(f, A, maxVerts)
    n=numel(f); nextUp=zeros(n,1);
    for v=1:n
        if ~isfinite(f(v)), continue; end
        nb=find(A(:,v)); nb=nb(isfinite(f(nb)));
        [mx,k]=max(f(nb));
        if isempty(nb) || mx<=f(v), nextUp(v)=v; else, nextUp(v)=nb(k); end
    end
    maxIdx=zeros(n,1); maxIdx(maxVerts)=1:numel(maxVerts);
    lbl=zeros(n,1);
    for v=1:n
        if ~isfinite(f(v)), continue; end
        path=v;
        while nextUp(path(end))~=path(end), path(end+1)=nextUp(path(end)); end %#ok<AGROW>
        peak=path(end); lbl(path)=maxIdx(peak);   % label whole ascent path (memoization-light)
    end
end
```

- [ ] **Step 4: Lint** — expect clean.

- [ ] **Step 5: Run the Step-1 snippet — expect PASS** (">=2 basins; seeds separated").

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/pet_epicenter.m
git commit -m "feat(pet): pet_epicenter gradient-ascent basins (Morse-Smale segmentation)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 3: Persistence ranking + filtering, and the `foci` output

**Files:**
- Modify: `toolbox/anatomy/pet_epicenter.m`

**Interfaces:**
- Consumes: `info.maxVerts`, `info.A`, `info.smoothed`, `basinLabel`.
- Produces: the populated `foci` struct array sorted by persistence desc (`.vertex .peak .persistence .basinArea`), filtered by `Opts.MinPersist`/`Opts.nFociMax`. Adds `local_persistence(f, A, maxVerts)`.

- [ ] **Step 1: Write the failing validation snippet**

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0002');
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
sW=in_tess_bst(wf); V=sW.Vertices; s1=3000; s2=15000;
f=1.0*exp(-sum((V-V(s1,:)).^2,2)/(2*0.012^2)) + 0.4*exp(-sum((V-V(s2,:)).^2,2)/(2*0.012^2));
f=f + 0.02*sin((1:size(V,1))'*12.9898);             % low-amplitude pseudo-noise
foci=pet_epicenter(wf, f, struct('HeatT',2e-5,'MinPersist',0.1));
assert(numel(foci)>=2, 'expected dominant + secondary focus');
d1=sqrt(sum((V(foci(1).vertex,:)-V(s1,:)).^2))*1000;  % dominant near the strong bump
assert(d1<6 && foci(1).persistence>=foci(2).persistence, 'dominant focus wrong or unsorted');
fprintf('PASS: dominant %.1fmm from strong seed, persist %.2f>=%.2f\n', d1, foci(1).persistence, foci(2).persistence);
```

- [ ] **Step 2: Run it — expect FAIL** (`foci` is empty `0×0` struct → `numel(foci)>=2`).

- [ ] **Step 3: Implement persistence + foci population** — replace the empty `foci=struct(...)` initializer block in `pet_epicenter` with:

```matlab
    % --- persistence of each maximum (super-level merge tree) ---
    pers=local_persistence(fs, A, maxVerts);
    if isempty(Opts.MinPersist)
        rng=max(fs(isfinite(fs)))-min(fs(isfinite(fs))); thr=0.15*rng;
    else, thr=Opts.MinPersist; end
    keep=pers>=thr; mk=maxVerts(keep); pk=pers(keep);
    [pk,ord]=sort(pk,'descend'); mk=mk(ord);
    if numel(mk)>Opts.nFociMax, mk=mk(1:Opts.nFociMax); pk=pk(1:Opts.nFociMax); end
    foci=struct('vertex',{},'peak',{},'persistence',{},'basinArea',{});
    for i=1:numel(mk)
        foci(i)=struct('vertex',mk(i),'peak',fs(mk(i)),'persistence',pk(i), ...
                       'basinArea',nnz(basinLabel==find(maxVerts==mk(i),1)));
    end
```

(Move this block to AFTER `basinLabel=local_basins(...)` so `basinLabel` exists; keep the `info=struct(...)` line last.)

Add the helper:
```matlab
% ===== persistence: peak height minus the level at which its basin merges into a higher one =====
function pers=local_persistence(f, A, maxVerts)
    n=numel(f); [~,order]=sort(f,'descend'); order=order(isfinite(f(order)));
    comp=zeros(n,1); peakOf=zeros(n,1);            % union-find component + its peak vertex
    isMax=false(n,1); isMax(maxVerts)=true; nextId=0;
    pers=inf(numel(maxVerts),1); maxId=zeros(n,1); maxId(maxVerts)=1:numel(maxVerts);
    for oi=1:numel(order)
        v=order(oi);
        nbUp=find(A(:,v)); nbUp=nbUp(comp(nbUp)>0);  % already-added (higher) neighbours
        roots=unique(arrayfun(@(x) local_find(comp,x), nbUp));
        if isempty(roots)
            nextId=nextId+1; comp(v)=v; peakOf(v)=v;  % new component born at a maximum
        else
            [~,hi]=max(f(arrayfun(@(r)peakOf(r),roots))); keepRoot=roots(hi);
            comp(v)=keepRoot;
            for r=roots(:)'
                if r==keepRoot, continue; end
                pv=peakOf(r);
                if isMax(pv), pers(maxId(pv))=f(pv)-f(v); end  % this peak dies at level f(v)
                comp(r)=keepRoot;                              % merge into the higher peak
            end
        end
    end
    pers(~isfinite(pers))=max(f(isfinite(f)))-min(f(isfinite(f)));  % the global max never dies
end
function r=local_find(comp,x)
    r=x; while comp(r)~=r && comp(r)~=0, r=comp(r); end
end
```

- [ ] **Step 4: Lint** — expect clean.

- [ ] **Step 5: Run the Step-1 snippet — expect PASS** (dominant near strong seed, persistence sorted).

- [ ] **Step 6: Commit**

```bash
git add toolbox/anatomy/pet_epicenter.m
git commit -m "feat(pet): pet_epicenter persistence ranking + ranked foci output

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 4: End-to-end validation on a real amyloid-positive subject

**Files:**
- Modify: `toolbox/anatomy/pet_epicenter.m` (docstring/Opts finalization only; no behavior change expected)

**Interfaces:**
- Consumes: the full `pet_epicenter` from Tasks 1–3.
- Produces: confidence that on a real amyloid+ subject the dominant focus is anatomically sensible (a Desikan cortical region with high SUVR), and that NaN medial-wall vertices are handled (no foci on the masked wall).

- [ ] **Step 1: Write the validation snippet** (real subject 0166, the strongest amyloid+ case; uses `pet_suvr` + the surface projection we already have)

```matlab
rehash path;
[sS,~]=bst_get('Subject','sub-MTL0166');
af=@(c) sS.Anatomy(find(strcmp({sS.Anatomy.Comment},c),1)).FileName;
wf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
pf=sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},'cortex_pial_low\.mat$','once')),1)).FileName;
sW=in_tess_bst(wf); sP=in_tess_bst(pf); sT1=in_mri_bst(af('MRI T1'));
sSuvr=pet_suvr(in_mri_bst(af('PET 18FNAV4694_mean_pvc')), in_mri_bst(af('ASEG')));
Vmid=0.5*(sW.Vertices+sP.Vertices); vox=cs_convert(sSuvr,'scs','voxel',Vmid);
cs=size(sSuvr.Cube); ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3);
suvr=nan(size(Vmid,1),1); suvr(ok)=interpn(double(sSuvr.Cube),vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN);
th=sqrt(sum((sP.Vertices-sW.Vertices).^2,2))*1000; suvr(th<1)=NaN;   % mask medial wall
foci=pet_epicenter(wf, suvr, struct('HeatT',2e-5));
assert(~isempty(foci), 'no foci detected on a strong amyloid+ subject');
assert(suvr(foci(1).vertex) > median(suvr,'omitnan'), 'dominant focus is not a high-SUVR vertex');
assert(th(foci(1).vertex) >= 1, 'dominant focus landed on the masked medial wall');
fprintf('PASS: %d foci; dominant SUVR=%.2f (cohort-median %.2f), persist=%.2f\n', ...
        numel(foci), suvr(foci(1).vertex), median(suvr,'omitnan'), foci(1).persistence);
```

- [ ] **Step 2: Run it — expect PASS** (foci found, dominant is a high-SUVR non-wall vertex). If it fails on NaN handling, fix `local_localmax`/`local_basins` to skip `~isfinite` (already coded) and re-run.

- [ ] **Step 3: Finalize the docstring** — ensure the header documents `Opts.HeatT/.MinPersist/.nFociMax`, the `foci`/`basinLabel`/`info` outputs, and `SEE ALSO: pet_suvr, tess_operators, tess_vertconn, bst_gradient`.

- [ ] **Step 4: Lint** — expect clean.

- [ ] **Step 5: Commit**

```bash
git add toolbox/anatomy/pet_epicenter.m
git commit -m "test(pet): validate pet_epicenter on amyloid+ subject + finalize docstring

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

### Task 5: Cohort analysis script — known-region validation, topology, density map

**Files:**
- Create: `dev/benchmarks/analyze_pet_epicenters.m`

**Interfaces:**
- Consumes: `pet_epicenter` (Tasks 1–4), `pet_suvr`, `tess_interp_tess2tess`.
- Produces: `R = analyze_pet_epicenters(subjects, tracer)` — per-subject foci + region; a cohort focus-density map (template space); the known-region validation fraction; the amyloid-vs-tau multifocality summary; a figure `dev/benchmarks/pet_epicenters_<tracer>.png`.

- [ ] **Step 1: Write the script**

Create `dev/benchmarks/analyze_pet_epicenters.m`:

```matlab
function R = analyze_pet_epicenters(subjects, tracer)
% ANALYZE_PET_EPICENTERS: cohort epicenter detection + validation vs known seeding regions.
% Per subject: surface SUVR -> pet_epicenter -> dominant focus region (Desikan); cohort focus
% density on the template (registration spheres); fraction of dominant foci in the expected
% regions (amyloid: precuneus/posteriorcingulate/isthmuscingulate/medialorbitofrontal; tau:
% entorhinal/inferiortemporal/fusiform); multifocality (number of persistent foci).
%
% Author: Diellor Basha, 2026
    if (nargin<2)||isempty(tracer), tracer='18FNAV4694'; end
    if strcmp(tracer,'18FNAV4694'), pet='PET 18FNAV4694_mean_pvc'; ...
        expect={'precuneus','posteriorcingulate','isthmuscingulate','medialorbitofrontal'};
    else, pet='PET 18Fflortaucipir_mean_pvc'; expect={'entorhinal','inferiortemporal','fusiform'}; end
    if (nargin<1)||isempty(subjects)
        ps=bst_get('ProtocolSubjects'); nm={ps.Subject.Name};
        subjects=nm(~cellfun('isempty',regexp(nm,'^sub-MTL\d+$','once')));
    end
    here=bst_fileparts(mfilename('fullpath'));
    sDef=bst_get('Subject',0);
    tWf=sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_white_low\.mat$','once')),1)).FileName;
    tW=in_tess_bst(tWf); dens=zeros(size(tW.Vertices,1),1);
    R=struct('subj',{},'tracer',{},'domRegion',{},'domPersist',{},'nFoci',{},'inExpected',{});
    for s=1:numel(subjects)
        subj=subjects{s}; [sS,~]=bst_get('Subject',subj); cmt={sS.Anatomy.Comment};
        af=@(c) sS.Anatomy(find(strcmp(cmt,c),1)).FileName;
        sf=@(rx) sS.Surface(find(~cellfun('isempty',regexp({sS.Surface.FileName},rx,'once')),1)).FileName;
        if ~any(strcmp(cmt,pet)), continue; end
        try
            wf=sf('cortex_white_low\.mat$'); pf=sf('cortex_pial_low\.mat$');
            sW=in_tess_bst(wf); sP=in_tess_bst(pf); sSuvr=pet_suvr(in_mri_bst(af(pet)),in_mri_bst(af('ASEG')));
            Vmid=0.5*(sW.Vertices+sP.Vertices); vox=cs_convert(sSuvr,'scs','voxel',Vmid);
            cs=size(sSuvr.Cube); ok=all(vox>=1,2)&vox(:,1)<=cs(1)&vox(:,2)<=cs(2)&vox(:,3)<=cs(3);
            suvr=nan(size(Vmid,1),1); suvr(ok)=interpn(double(sSuvr.Cube),vox(ok,1),vox(ok,2),vox(ok,3),'linear',NaN);
            th=sqrt(sum((sP.Vertices-sW.Vertices).^2,2))*1000; suvr(th<1)=NaN;
            foci=pet_epicenter(wf, suvr, struct('HeatT',2e-5));
            if isempty(foci), continue; end
            dv=foci(1).vertex; reg=local_region(sW, dv);
            inExp=any(cellfun(@(e) ~isempty(strfind(lower(reg),e)), expect));   %#ok<STREMP>
            R(end+1)=struct('subj',subj,'tracer',tracer,'domRegion',reg,'domPersist',foci(1).persistence,'nFoci',numel(foci),'inExpected',inExp); %#ok<AGROW>
            % accumulate cohort density: map each focus vertex to the template
            Wm=tess_interp_tess2tess(wf, tWf, 0, 0, 0);
            ind=zeros(size(sW.Vertices,1),1); ind([foci.vertex])=1; dens=dens+(Wm*ind > 0.5);
        catch ME, fprintf('  %s ERR %s\n',subj,ME.message); end
    end
    frac=mean([R.inExpected]); 
    fprintf('\n%s: %d subjects | dominant focus in expected region: %.0f%% | median #foci=%.1f\n', ...
        tracer, numel(R), 100*frac, median([R.nFoci]));
    tbl=tabulate({R.domRegion});  % region frequencies of the dominant focus
    [~,o]=sort(cell2mat(tbl(:,2)),'descend');
    fprintf('  top dominant-focus regions:\n'); for k=1:min(6,size(tbl,1)), fprintf('     %-26s %d\n', tbl{o(k),1}, tbl{o(k),2}); end
    % figure: cohort focus-density on the template pial (lateral + medial of LH)
    f=figure('Visible','off','Position',[40 40 900 420]);
    tP=in_tess_bst(sDef.Surface(find(~cellfun('isempty',regexp({sDef.Surface.FileName},'cortex_pial_low\.mat$','once')),1)).FileName);
    [rH,lH]=tess_hemisplit(tP);
    for s2=1:2, vw={[180 -10],[0 -10]}; idx=double(lH(:)); keep=all(ismember(tP.Faces,idx),2);
        rmp=zeros(size(tP.Vertices,1),1); rmp(idx)=1:numel(idx); Fh=rmp(tP.Faces(keep,:));
        subplot(1,2,s2); patch('Faces',Fh,'Vertices',tP.Vertices(idx,:),'FaceVertexCData',dens(idx),'FaceColor','interp','EdgeColor','none');
        clim([0 max(dens)]); colormap(hot); view(vw{s2}); axis equal off vis3d; camlight headlight; lighting gouraud;
        title(sprintf('%s focus density (LH %s)',tracer,{'lateral','medial'}{s2}),'Interpreter','none');
    end
    png=fullfile(here,sprintf('pet_epicenters_%s.png',tracer)); print(f,png,'-dpng','-r110'); close(f);
    fprintf('  figure -> %s\n', png); R=struct('rows',{R},'density',dens);
end

function reg=local_region(sW, v)
    A=sW.Atlas; ai=find(~cellfun('isempty',regexp({A.Name},'Desikan|aparc','once')),1); reg='unknown';
    for i=1:numel(A(ai).Scouts), if any(A(ai).Scouts(i).Vertices==v), reg=A(ai).Scouts(i).Label; return; end; end
end
```

- [ ] **Step 2: Lint** `dev/benchmarks/analyze_pet_epicenters.m` — expect clean.

- [ ] **Step 3: Run on the amyloid+ subgroup first** (where epicenters are meaningful):

```matlab
rehash path;
R=analyze_pet_epicenters({'sub-MTL0166','sub-MTL0079','sub-MTL0018','sub-MTL0054','sub-MTL0268','sub-MTL0284','sub-MTL0311'}, '18FNAV4694');
```
Expected: a printed "dominant focus in expected region: X%" line, a top-regions list dominated by precuneus/posterior-cingulate, and a saved figure. Sanity gate: ≥50% of amyloid+ dominant foci in the expected DMN/precuneus regions.

- [ ] **Step 4: Run the full cohort** (all 66) for amyloid and tau:

```matlab
Ra=analyze_pet_epicenters([], '18FNAV4694');
Rt=analyze_pet_epicenters([], '18Fflortaucipir');
```
Expected: two figures + the validation/topology printout. (Whole-cohort amyloid is mostly negative, so the *amyloid+ subgroup* result from Step 3 is the scientific anchor; the full-cohort run produces the density map + the topology comparison.)

- [ ] **Step 5: Commit** (script + the two figures)

```bash
git add dev/benchmarks/analyze_pet_epicenters.m dev/benchmarks/pet_epicenters_18FNAV4694.png dev/benchmarks/pet_epicenters_18Fflortaucipir.png
git commit -m "bench(pet): cohort epicenter analysis - known-region validation + density map

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>
Claude-Session: https://claude.ai/code/session_01BTbrKJ8SkHiQg9hX2gkaPi"
```

---

## Self-Review

**Spec coverage:** Method steps 1–5 (heat-smooth, gradient/maxima, basins, persistence) → Tasks 1–3; the covariant-gradient flow is noted optional and deferred (mesh-adjacency discrete Morse is the implemented realization — consistent with the design's "Morse-Smale / persistence" core). Group analysis & validation (location vs known regions, topology, density map) → Task 5. Robustness-vs-argmax and persistence-as-marker are not yet separate tasks — they're cheap follow-ups on `R` (the script returns per-subject `domPersist`/`nFoci`); add them after Task 5 if desired. Data scope (subject-space detection, template density via spheres) → Task 5.

**Placeholder scan:** all steps carry real MATLAB + real assertions; no TBD/TODO.

**Type consistency:** `foci` fields (`.vertex .peak .persistence .basinArea`) are defined in Task 3 and consumed unchanged in Tasks 4–5; `info.maxVerts`/`info.A` defined in Task 1 and consumed in Tasks 2–3; `pet_epicenter(SurfaceFile, suvrMap, Opts)` signature stable across all tasks.

**Known caveat to watch during execution:** `tess_interp_tess2tess` indexing for the density accumulation maps full-surface focus indicators; verify the `Wm` orientation (`[nTemplate × nSubject]`) on first run and transpose if the density comes out empty (same check used in `group_surface_suvr`).
