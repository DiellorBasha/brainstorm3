# Eigenmode Accuracy Benchmark Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a synthetic-on-real-cortex benchmark that shows the LBO-eigenmode inverse (`eig_mne/log`, `eig_dspm/log`) is competitive in localization accuracy with `wMNE`/`dSPM`/`sLORETA` across focal/patch/distributed regimes, an SNR sweep, and an eigenmode-count (K) sweep — producing CSV statistics and five MATLAB figures.

**Architecture:** Plain MATLAB scripts under `dev/benchmarks/` that orchestrate already-tested primitives (`bst_benchmark_sources`, `bst_benchmark_simulate`, `bst_benchmark_inverse`, `bst_benchmark_metrics`). The two real protocols (`TutorialAuditory`, `TutorialNeuromag`) supply real cortex + head model + noise covariance; a fixture step computes cortex eigenmodes (1000/hemisphere) and the composed eigenmode head model. The only library change is adding an optional `nModes` passthrough to `bst_benchmark_inverse` for the K-sweep.

**Tech Stack:** MATLAB R2023b, Brainstorm (running, `brainstorm nogui`), MATLAB MCP for execution. No Statistics Toolbox dependency (Wilcoxon signed-rank and box/strip plots are implemented locally).

---

## Background the engineer must know

- **Run everything through the MATLAB MCP** (`mcp__plugin_brainstorm-dev_MATLAB__run_matlab_file` for scripts, `evaluate_matlab_code` for one-liners). Brainstorm is already started in the session. The two protocols are already imported.
- **Repo root:** `/Users/diellorbasha/workspace/research/code/brainstorm3`. All benchmark files live in `dev/benchmarks/`. Tests live in `dev/tests/` and follow the existing convention: a parameterless `function test_name`, `addpath(repoRoot)`, `if ~brainstorm('status'); brainstorm nogui; end`, `assert(...)`, ending with `disp('ALL TESTS PASSED')`.
- **Primitive signatures (already implemented, do not change except where stated):**
  - `S = bst_benchmark_sources(Surface, regime, Name,Value...)` where `Surface=struct('Vertices',[nV×3],'VertConn',[nV×nV sparse])`, `regime∈{'focal','patch','distributed'}`. Returns `S.GT [nV×1]`, `S.Sources [nV×nTime]`, `S.SeedVertex`. Options: `'nTime'`, `'Seed'`, `'SeedVertex'`, `'Radius'` (patch), `'Sigma'` (distributed).
  - `Sim = bst_benchmark_simulate(L, Sources, NoiseCov, 'SNR',snrDb, 'Seed',k)` → `Sim.F [nGoodCh×nTime]`. `L` is the **good-channel** constrained leadfield `[nGoodCh×nV]`.
  - `Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)` → struct with fields `.wmne .dspm .sloreta .eig_mne_log .eig_dspm_log`, each `[nV×nTime]`. **This task adds an optional trailing `nModes`.**
  - `M = bst_benchmark_metrics(gt, est, GridLoc, seedVertex)` → `.LocError` (mm), `.Correlation`, `.NRMSE`, `.AUC`, `.SpatialDispersion` (mm). `gt,est` are `[nV×1]`; `GridLoc` is `[nV×3]` in **metres**.
  - `bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, 'Method',m,'Prior','log','SNR',snrLin,'nModes',K)` already honors `'nModes'` (caps to first K of the globally-sorted composed model).
- **K semantics:** `tess_eigenmodes` computes `nModes` **per connected component** (per hemisphere). `bst_eigenmode_leadfield` keeps the globally-lowest-eigenvalue modes in ascending order, so capping the inverse to first `Ktot` modes ≈ `Ktot/2` per hemisphere. The K-sweep therefore uses **total** counts `{600,1200,2000}` ≈ `{300,600,1000}`/hemisphere. Compute eigenmodes once at `nmodes=1000` (per hemi) so 2000 total are available.
- **Reducing time before scoring:** estimates are `[nV×nTime]`; score at the peak-source time sample `tStar = argmax_t ||Sources(:,t)||`. Use `gt=S.GT`, `estMap=Est.(method)(:,tStar)`, `seedVertex=S.SeedVertex`.

---

## File structure

| File | Responsibility |
|---|---|
| `toolbox/inverse/bst_benchmark_inverse.m` (modify) | Accept optional `nModes`, pass to the eigenmode inverse |
| `dev/benchmarks/bench_config.m` (create) | Return the benchmark config struct; `bench_config('smoke')` returns the fast preset |
| `dev/benchmarks/bench_fixtures.m` (create) | Ensure each anatomy has cortex eigenmodes + composed eigenmode head model (idempotent) |
| `dev/benchmarks/bench_run.m` (create) | Run the synthetic loop, write `synthetic.csv` |
| `dev/benchmarks/bench_stats.m` (create) | Aggregate medians/IQR + paired Wilcoxon, write `stats.csv` + `stats.md` |
| `dev/benchmarks/bench_figures.m` (create) | Produce the five figures (PNG + .fig) |
| `dev/benchmarks/benchmark_eigenmodes.m` (create) | Top-level driver: config→fixtures→run→stats→figures→`REPORT.md` |
| `dev/tests/test_benchmark_inverse_nmodes_e2e.m` (create) | Verify `nModes` passthrough changes mode count and stays finite |
| `dev/tests/test_bench_config_pure.m` (create) | Validate config fields and smoke preset |
| `dev/tests/test_bench_stats_pure.m` (create) | Validate aggregation + signed-rank on synthetic rows |
| `dev/tests/test_bench_figures_pure.m` (create) | Validate five PNGs are produced from a synthetic CSV |
| `dev/tests/test_benchmark_eigenmodes_smoke.m` (create) | End-to-end smoke on the `smoke` preset |

---

## Task 1: Add `nModes` passthrough to `bst_benchmark_inverse`

**Files:**
- Modify: `toolbox/inverse/bst_benchmark_inverse.m:1` (signature) and `:63-64` (eigenmode call)
- Test: `dev/tests/test_benchmark_inverse_nmodes_e2e.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_inverse_nmodes_e2e
% Verify the optional nModes arg caps the eigenmode reconstruction and stays finite.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Find a study with base surface HM + noise cov + cortex eigenmodes (created by bench_fixtures).
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
T = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.HeadModel) || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        if (~isfield(hm,'isEigenmode')||~hm.isEigenmode) && strcmpi(hm.HeadModelType,'surface'); iBase = ih; break; end
    end
    if isempty(iBase); continue; end
    hmB = in_bst_headmodel(s.HeadModel(iBase).FileName,0,'SurfaceFile');
    [~, isE] = in_tess_eigenmodes(hmB.SurfaceFile);
    if isE; T = struct('iStudy',iS,'iBase',iBase); break; end
end
if isempty(T); disp('SKIP: no study with base HM + eigenmodes (run bench_fixtures first).'); return; end

s = sStudies.Study(T.iStudy);
baseHmFile = s.HeadModel(T.iBase).FileName;
ncFile     = s.NoiseCov(1).FileName;
chFile     = bst_get('ChannelFileForStudy', s.FileName);

baseHM = in_bst_headmodel(baseHmFile, 1);   % [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess_bst(baseHM.SurfaceFile);
if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
    Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
end
S  = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn),'focal','nTime',5,'Seed',1);
NC = load(file_fullpath(ncFile)); C = NC.NoiseCov(goodMask,goodMask);
Sim = bst_benchmark_simulate(L, S.Sources, C, 'SNR', 6, 'Seed', 1);

Est600  = bst_benchmark_inverse(Sim.F, baseHmFile, ncFile, chFile, goodMask, 6, 600);
Est2000 = bst_benchmark_inverse(Sim.F, baseHmFile, ncFile, chFile, goodMask, 6, 2000);

assert(isfield(Est600,'eig_mne_log') && all(isfinite(Est600.eig_mne_log(:))), 'eig estimate must be finite at K=600.');
assert(isfield(Est2000,'eig_mne_log') && all(isfinite(Est2000.eig_mne_log(:))), 'eig estimate must be finite at K=2000.');
assert(norm(Est600.eig_mne_log(:) - Est2000.eig_mne_log(:)) > 0, 'different K must give different eig reconstruction.');
assert(all(isfinite(Est600.wmne(:))), 'standard methods must remain finite.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run via MCP `run_matlab_file` on `dev/tests/test_benchmark_inverse_nmodes_e2e.m`.
Expected: FAIL — `bst_benchmark_inverse` errors on the 7th argument (`Too many input arguments`), OR SKIP if fixtures not yet built. If SKIP, run Task 3 fixtures first on one anatomy, then re-run; it must then FAIL on the extra argument.

- [ ] **Step 3: Modify `bst_benchmark_inverse`**

Change the signature line 1 from:

```matlab
function Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)
```

to:

```matlab
function Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR, nModes)
% nModes (optional): cap eigenmode count for the K-sweep (default [] = all composed modes).
if nargin < 7; nModes = []; end
```

Change the eigenmode call (currently lines 63-64) from:

```matlab
        [Inv, errE] = bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, ...
            'Method', eigCfg{k,1}, 'Prior', 'log', 'SNR', snrLin);
```

to:

```matlab
        [Inv, errE] = bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, ...
            'Method', eigCfg{k,1}, 'Prior', 'log', 'SNR', snrLin, 'nModes', nModes);
```

(`bst_inverse_eigenmodes` already treats `nModes=[]` as "all", and `Phi(:,1:Inv.nModes)` on line 67 picks up the capped count automatically.)

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_benchmark_inverse_nmodes_e2e.m`. Expected: `ALL TESTS PASSED` (requires Task 3 fixtures on at least one anatomy).

- [ ] **Step 5: Commit**

```bash
git add toolbox/inverse/bst_benchmark_inverse.m dev/tests/test_benchmark_inverse_nmodes_e2e.m
git commit -m "feat(benchmark): nModes passthrough in bst_benchmark_inverse for K-sweep"
```

---

## Task 2: `bench_config`

**Files:**
- Create: `dev/benchmarks/bench_config.m`
- Test: `dev/tests/test_bench_config_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_bench_config_pure
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));

C = bench_config();
req = {'anatomies','methods','regimes','snr_db','k_total','nReplicates','seed','outDir','nModes_eig'};
for i=1:numel(req); assert(isfield(C,req{i}), 'config missing field %s', req{i}); end
assert(numel(C.anatomies)==2, 'expect 2 anatomies (Auditory, Neuromag).');
assert(isequal(sort(C.k_total), [600 1200 2000]), 'K-sweep must be {600,1200,2000} total.');
assert(C.nModes_eig==1000, 'eigenmodes computed at 1000 per hemisphere.');
assert(all(ismember({'focal','patch','distributed'}, C.regimes)), 'three regimes required.');

Csm = bench_config('smoke');
assert(numel(Csm.anatomies)==1, 'smoke uses 1 anatomy.');
assert(isequal(Csm.regimes, {'focal'}), 'smoke uses focal only.');
assert(numel(Csm.snr_db)==2 && Csm.nReplicates==2, 'smoke is tiny.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_bench_config_pure.m`. Expected: FAIL — `Undefined function 'bench_config'`.

- [ ] **Step 3: Create `dev/benchmarks/bench_config.m`**

```matlab
function C = bench_config(preset)
% BENCH_CONFIG: Configuration for the eigenmode accuracy benchmark.
% USAGE:  C = bench_config()         % full benchmark
%         C = bench_config('smoke')  % fast pipeline-validation preset
if nargin < 1 || isempty(preset); preset = 'full'; end

% Each anatomy: a protocol name + the subject whose cortex/headmodel to use.
audi = struct('label','auditory','protocol','TutorialAuditory','subject','Subject01');
neur = struct('label','neuromag','protocol','TutorialNeuromag','subject','Subject01');

C = struct();
C.methods      = {'wmne','dspm','sloreta','eig_mne_log','eig_dspm_log'};
C.eigMethods   = {'eig_mne_log','eig_dspm_log'};
C.stdMethods   = {'wmne','dspm','sloreta'};
C.regimes      = {'focal','patch','distributed'};
C.regimeOpts   = struct('focal',{{}}, 'patch',{{'Radius',2}}, 'distributed',{{'Sigma',0.01}});
C.snr_db       = [2 4 6 10 20];
C.k_total      = [600 1200 2000];     % total modes ~ {300,600,1000}/hemisphere
C.nModes_eig   = 1000;                % per hemisphere, at compute time
C.nReplicates  = 15;
C.nTime        = 20;
C.seed         = 20260602;            % master seed (date-derived, fixed for reproducibility)
C.plateauTol_mm= 1.0;                 % focal LocError improvement below this => plateau

switch lower(preset)
    case 'smoke'
        C.anatomies   = {audi};
        C.regimes     = {'focal'};
        C.regimeOpts  = struct('focal',{{}});
        C.snr_db      = [4 10];
        C.k_total     = [600 2000];
        C.nReplicates = 2;
        C.outDir      = fullfile(fileparts(mfilename('fullpath')), 'smoke_run');
    otherwise
        C.anatomies   = {audi, neur};
        C.outDir      = fullfile(fileparts(mfilename('fullpath')), 'eigenmode_accuracy_run');
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_bench_config_pure.m`. Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/bench_config.m dev/tests/test_bench_config_pure.m
git commit -m "feat(benchmark): bench_config with full + smoke presets"
```

---

## Task 3: `bench_fixtures` — ensure eigenmodes + composed leadfield per anatomy

**Files:**
- Create: `dev/benchmarks/bench_fixtures.m`
- Test: covered by re-running Task 1's `test_benchmark_inverse_nmodes_e2e.m` after fixtures (no separate pure test — this mutates the DB and is validated by the downstream e2e).

- [ ] **Step 1: Create `dev/benchmarks/bench_fixtures.m`**

```matlab
function info = bench_fixtures(anat, nModesPerHemi)
% BENCH_FIXTURES: Idempotently ensure an anatomy's protocol has cortex eigenmodes
% and a composed eigenmode head model, ready for bst_benchmark_inverse.
%
% USAGE:  info = bench_fixtures(anat, nModesPerHemi)
%   anat          : struct with .protocol, .subject (from bench_config)
%   nModesPerHemi : eigenmodes to compute per hemisphere (e.g. 1000)
% RETURNS info: .iStudy .baseHmFile .ncFile .chFile .eigHmFile

% Select the protocol
iProt = bst_get('Protocol', anat.protocol);
if isempty(iProt); error('bench_fixtures: protocol %s not found.', anat.protocol); end
gui_brainstorm('SetCurrentProtocol', iProt);

% Find a study with a base surface head model + noise covariance + a data/raw file
sStudies = bst_get('ProtocolStudies');
T = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.HeadModel) || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Data); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        if (~isfield(hm,'isEigenmode')||~hm.isEigenmode) && strcmpi(hm.HeadModelType,'surface'); iBase = ih; break; end
    end
    if ~isempty(iBase); T = struct('iStudy',iS,'iBase',iBase); break; end
end
if isempty(T); error('bench_fixtures: no study with base surface HM + noise cov + data in %s.', anat.protocol); end
s = sStudies.Study(T.iStudy);
baseHmFile = s.HeadModel(T.iBase).FileName;
ncFile     = s.NoiseCov(1).FileName;
chFile     = bst_get('ChannelFileForStudy', s.FileName);
dataFile   = s.Data(1).FileName;

% (1) Cortex eigenmodes: compute if absent (idempotent via 'overwrite'=0)
hmB = in_bst_headmodel(baseHmFile, 0, 'SurfaceFile');
[~, hasEig] = in_tess_eigenmodes(hmB.SurfaceFile);
if ~hasEig
    bst_process('CallProcess', 'process_eigenmodes', [], [], ...
        'subjectname', anat.subject, ...
        'surfacetype', 'cortex', ...
        'nmodes',      {nModesPerHemi, '', 0}, ...
        'masstype',    'barycentric', ...
        'removedc',    1, ...
        'repair',      0, ...
        'overwrite',   0);
end

% (2) Composed eigenmode head model: build if absent (search again after compute)
sStudies = bst_get('ProtocolStudies');
s = sStudies.Study(T.iStudy);
eigHmFile = '';
for ih = 1:numel(s.HeadModel)
    try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode; eigHmFile = s.HeadModel(ih).FileName; break; end
end
if isempty(eigHmFile)
    bst_process('CallProcess', 'process_eigenmode_leadfield', dataFile, [], ...
        'nmodes', {0, '', 0});   % 0 = all available (~2*nModesPerHemi total)
    sStudies = bst_get('ProtocolStudies');
    s = sStudies.Study(T.iStudy);
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        if isfield(hm,'isEigenmode') && hm.isEigenmode; eigHmFile = s.HeadModel(ih).FileName; break; end
    end
end
if isempty(eigHmFile); error('bench_fixtures: failed to create eigenmode head model in %s.', anat.protocol); end

info = struct('iStudy',T.iStudy,'baseHmFile',baseHmFile,'ncFile',ncFile, ...
              'chFile',chFile,'eigHmFile',eigHmFile,'dataFile',dataFile, ...
              'surfaceFile',hmB.SurfaceFile);
end
```

- [ ] **Step 2: Run fixtures on the Auditory anatomy (manual verification)**

Via MCP `evaluate_matlab_code`:

```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks');
C = bench_config();
info = bench_fixtures(C.anatomies{1}, C.nModes_eig);
disp(info);
```

Expected: prints a struct with non-empty `eigHmFile`; no error. (First run computes ~1000 modes/hemi — up to ~1–2 min.)

- [ ] **Step 3: Verify Task 1's e2e test now passes**

Run `dev/tests/test_benchmark_inverse_nmodes_e2e.m`. Expected: `ALL TESTS PASSED` (no longer SKIP).

- [ ] **Step 4: Commit**

```bash
git add dev/benchmarks/bench_fixtures.m
git commit -m "feat(benchmark): bench_fixtures ensures cortex eigenmodes + composed leadfield"
```

---

## Task 4: `bench_run` — synthetic loop → `synthetic.csv`

**Files:**
- Create: `dev/benchmarks/bench_run.m`
- Test: covered by the smoke e2e (Task 7); this task includes a manual smoke check.

- [ ] **Step 1: Create `dev/benchmarks/bench_run.m`**

```matlab
function csvPath = bench_run(C)
% BENCH_RUN: Execute the synthetic-on-real-cortex benchmark; write synthetic.csv.
% USAGE:  csvPath = bench_run(bench_config())
if ~exist(C.outDir,'dir'); mkdir(C.outDir); end
csvPath = fullfile(C.outDir, 'synthetic.csv');

rows = {};   % each row: {anatomy, regime, snr_db, replicate, method, K, locerror_mm, corr, nrmse, auc, dispersion_mm}
for ia = 1:numel(C.anatomies)
    anat = C.anatomies{ia};
    info = bench_fixtures(anat, C.nModes_eig);

    % Constrained good-channel leadfield + cortex geometry + noise cov
    baseHM = in_bst_headmodel(info.baseHmFile, 1);          % [nCh x nVert]
    goodMask = all(isfinite(double(baseHM.Gain)),2);
    L = double(baseHM.Gain(goodMask,:));
    Surf = in_tess_bst(info.surfaceFile);
    if ~isfield(Surf,'VertConn') || isempty(Surf.VertConn)
        Surf.VertConn = tess_vertconn(Surf.Vertices, Surf.Faces);
    end
    GridLoc = Surf.Vertices;                                 % metres
    Surface = struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn);
    NC = load(file_fullpath(info.ncFile)); Cnoise = NC.NoiseCov(goodMask,goodMask);

    for ir = 1:numel(C.regimes)
        regime = C.regimes{ir};
        ropts  = C.regimeOpts.(regime);
        for isnr = 1:numel(C.snr_db)
            snr = C.snr_db(isnr);
            for rep = 1:C.nReplicates
                seed = mod(C.seed + 1000*ia + 100*ir + 10*isnr + rep, 2^31-1);
                S = bst_benchmark_sources(Surface, regime, 'nTime', C.nTime, 'Seed', seed, ropts{:});
                Sim = bst_benchmark_simulate(L, S.Sources, Cnoise, 'SNR', snr, 'Seed', seed);
                [~, tStar] = max(sum(S.Sources.^2,1));
                gt = S.GT(:);

                for jk = 1:numel(C.k_total)
                    Ktot = C.k_total(jk);
                    try
                        Est = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, ...
                            info.chFile, goodMask, snr, Ktot);
                    catch ME
                        warning('bench_run: inverse failed (%s/%s/snr%d/rep%d/K%d): %s', ...
                            anat.label, regime, snr, rep, Ktot, ME.message);
                        continue;
                    end
                    % Eigenmode methods: one row per K. Standard methods: only on first K.
                    if jk == 1; theseMethods = C.methods; else; theseMethods = C.eigMethods; end
                    for im = 1:numel(theseMethods)
                        meth = theseMethods{im};
                        if ~isfield(Est, meth); continue; end
                        estMap = Est.(meth)(:, tStar);
                        if ~all(isfinite(estMap)); continue; end
                        M = bst_benchmark_metrics(gt, estMap, GridLoc, S.SeedVertex);
                        if ismember(meth, C.eigMethods); Kcol = Ktot; else; Kcol = NaN; end
                        rows(end+1,:) = {anat.label, regime, snr, rep, meth, Kcol, ...
                            M.LocError, M.Correlation, M.NRMSE, M.AUC, M.SpatialDispersion}; %#ok<AGROW>
                    end
                end
            end
        end
        fprintf('BENCH> %s / %s done (%d rows so far)\n', anat.label, regime, size(rows,1));
    end
end

% Write CSV
hdr = {'anatomy','regime','snr_db','replicate','method','K','locerror_mm', ...
       'correlation','nrmse','auc','spatial_dispersion_mm'};
fid = fopen(csvPath,'w');
fprintf(fid, '%s\n', strjoin(hdr, ','));
for i = 1:size(rows,1)
    fprintf(fid, '%s,%s,%g,%d,%s,%g,%g,%g,%g,%g,%g\n', rows{i,1}, rows{i,2}, rows{i,3}, ...
        rows{i,4}, rows{i,5}, rows{i,6}, rows{i,7}, rows{i,8}, rows{i,9}, rows{i,10}, rows{i,11});
end
fclose(fid);
fprintf('BENCH> wrote %d rows to %s\n', size(rows,1), csvPath);
end
```

- [ ] **Step 2: Manual smoke check**

Via MCP `evaluate_matlab_code`:

```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks');
C = bench_config('smoke');
p = bench_run(C);
T = readtable(p);
disp(head(T));
fprintf('rows=%d methods=%s\n', height(T), strjoin(unique(T.method)', ','));
```

Expected: a non-empty table; `methods` includes all five; `K` is `NaN` for standard rows and `{600,2000}` for eig rows; all `locerror_mm` finite and ≥ 0.

- [ ] **Step 3: Commit**

```bash
git add dev/benchmarks/bench_run.m
git commit -m "feat(benchmark): bench_run synthetic loop writing synthetic.csv"
```

---

## Task 5: `bench_stats` — medians/IQR + paired Wilcoxon

**Files:**
- Create: `dev/benchmarks/bench_stats.m`
- Test: `dev/tests/test_bench_stats_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_bench_stats_pure
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));

% Build a tiny synthetic CSV: two methods, focal, snr 10, 6 replicates.
tmp = fullfile(tempdir, sprintf('bench_stats_test_%d', feature('getpid')));
if ~exist(tmp,'dir'); mkdir(tmp); end
csvPath = fullfile(tmp,'synthetic.csv');
fid = fopen(csvPath,'w');
fprintf(fid,'anatomy,regime,snr_db,replicate,method,K,locerror_mm,correlation,nrmse,auc,spatial_dispersion_mm\n');
for r=1:6
    fprintf(fid,'auditory,focal,10,%d,wmne,NaN,%g,0.9,0.1,0.99,30\n', r, 10+r);     % 11..16
    fprintf(fid,'auditory,focal,10,%d,eig_mne_log,2000,%g,0.9,0.1,0.99,30\n', r, 12+r); % 13..18 (eig +2mm)
end
fclose(fid);

St = bench_stats(csvPath, tmp);
% Aggregated medians: wmne focal = median(11..16)=13.5; eig=median(13..18)=15.5
w = St.summary(strcmp(St.summary.method,'wmne') & strcmp(St.summary.regime,'focal'), :);
e = St.summary(strcmp(St.summary.method,'eig_mne_log') & strcmp(St.summary.regime,'focal'), :);
assert(abs(w.median_locerror_mm - 13.5) < 1e-9, 'wmne median LocError wrong.');
assert(abs(e.median_locerror_mm - 15.5) < 1e-9, 'eig median LocError wrong.');
% Paired Wilcoxon eig vs wmne: constant +2mm difference over 6 pairs -> p < 0.05
cmp = St.compare(strcmp(St.compare.eig_method,'eig_mne_log') & strcmp(St.compare.std_method,'wmne') ...
                 & strcmp(St.compare.regime,'focal'), :);
assert(~isempty(cmp), 'comparison row missing.');
assert(abs(cmp.median_diff_mm - 2) < 1e-9, 'median paired difference must be +2 mm.');
assert(cmp.p_value < 0.05, 'consistent +2mm shift over 6 pairs must be significant.');
assert(exist(fullfile(tmp,'stats.csv'),'file')==2, 'stats.csv must be written.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_bench_stats_pure.m`. Expected: FAIL — `Undefined function 'bench_stats'`.

- [ ] **Step 3: Create `dev/benchmarks/bench_stats.m`**

```matlab
function St = bench_stats(csvPath, outDir)
% BENCH_STATS: Aggregate benchmark rows -> median/IQR per (method,regime,snr,K);
% paired Wilcoxon signed-rank (each eig method vs each standard method, matched on
% anatomy/regime/snr/replicate at the eig method's largest K). Writes stats.csv + stats.md.
if nargin < 2 || isempty(outDir); outDir = fileparts(csvPath); end
T = readtable(csvPath);
T.method = string(T.method); T.regime = string(T.regime); T.anatomy = string(T.anatomy);

% ---- Summary: median + IQR of LocError per (method,regime) pooled over anatomy/snr/replicate
%      (eig methods summarized at their max K = plateau proxy) ----
eigMethods = ["eig_mne_log","eig_dspm_log"];
stdMethods = ["wmne","dspm","sloreta"];
regimes = unique(T.regime,'stable');
Kmax = max(T.K(~isnan(T.K)));

methods = [stdMethods, eigMethods];
S = table();
for ir = 1:numel(regimes)
    for im = 1:numel(methods)
        m = methods(im);
        if ismember(m, eigMethods)
            sel = T.method==m & T.regime==regimes(ir) & T.K==Kmax;
        else
            sel = T.method==m & T.regime==regimes(ir);
        end
        v = T.locerror_mm(sel);
        if isempty(v); continue; end
        row = table(string(regimes(ir)), m, median(v), iqr_local(v), mean(v), numel(v), ...
            'VariableNames', {'regime','method','median_locerror_mm','iqr_locerror_mm','mean_locerror_mm','n'});
        S = [S; row]; %#ok<AGROW>
    end
end

% ---- Paired comparisons: eig (at Kmax) vs std, matched keys ----
Cmp = table();
keyVars = {'anatomy','regime','snr_db','replicate'};
for ir = 1:numel(regimes)
    for ie = 1:numel(eigMethods)
        eTab = T(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==Kmax, :);
        for is = 1:numel(stdMethods)
            sTab = T(T.method==stdMethods(is) & T.regime==regimes(ir), :);
            [d, n] = paired_diffs(eTab, sTab, keyVars, 'locerror_mm');
            if n < 1; continue; end
            p = local_signrank(d);
            row = table(string(regimes(ir)), eigMethods(ie), stdMethods(is), median(d), n, p, ...
                'VariableNames', {'regime','eig_method','std_method','median_diff_mm','n_pairs','p_value'});
            Cmp = [Cmp; row]; %#ok<AGROW>
        end
    end
end

% ---- K-sweep curve table (focal emphasis): median LocError per (regime, eig method, K) ----
Ksweep = table();
Kvals = unique(T.K(~isnan(T.K)));
for ir = 1:numel(regimes)
    for ie = 1:numel(eigMethods)
        for ik = 1:numel(Kvals)
            v = T.locerror_mm(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==Kvals(ik));
            if isempty(v); continue; end
            Ksweep = [Ksweep; table(string(regimes(ir)), eigMethods(ie), Kvals(ik), median(v), ...
                'VariableNames',{'regime','eig_method','K','median_locerror_mm'})]; %#ok<AGROW>
        end
    end
end

St = struct('summary',S,'compare',Cmp,'ksweep',Ksweep);
writetable(S,   fullfile(outDir,'stats.csv'));
writetable(Cmp, fullfile(outDir,'stats_compare.csv'));
writetable(Ksweep, fullfile(outDir,'stats_ksweep.csv'));

% Markdown summary
fid = fopen(fullfile(outDir,'stats.md'),'w');
fprintf(fid,'# Eigenmode accuracy benchmark — statistics\n\n## Median LocError (mm) by method × regime\n\n');
fprintf(fid,'| regime | method | median | IQR | n |\n|---|---|---|---|---|\n');
for i=1:height(S)
    fprintf(fid,'| %s | %s | %.2f | %.2f | %d |\n', S.regime(i), S.method(i), ...
        S.median_locerror_mm(i), S.iqr_locerror_mm(i), S.n(i));
end
fprintf(fid,'\n## Paired Wilcoxon (eig vs standard), LocError\n\n');
fprintf(fid,'| regime | eig | vs | median diff (mm) | n | p |\n|---|---|---|---|---|---|\n');
for i=1:height(Cmp)
    fprintf(fid,'| %s | %s | %s | %+.2f | %d | %.4g |\n', Cmp.regime(i), Cmp.eig_method(i), ...
        Cmp.std_method(i), Cmp.median_diff_mm(i), Cmp.n_pairs(i), Cmp.p_value(i));
end
fclose(fid);
end

% ---- helpers ----
function q = iqr_local(v); q = quantile_local(v,0.75) - quantile_local(v,0.25); end

function q = quantile_local(v, p)
v = sort(v(:)); n = numel(v);
if n==1; q=v; return; end
h = (n-1)*p + 1; lo = floor(h); hi = min(lo+1,n);
q = v(lo) + (h-lo)*(v(hi)-v(lo));
end

function [d, n] = paired_diffs(eTab, sTab, keyVars, valVar)
% Match rows by keyVars; return eig - std differences.
ke = strcat(string(eTab.anatomy),'|',string(eTab.regime),'|',string(eTab.snr_db),'|',string(eTab.replicate));
ks = strcat(string(sTab.anatomy),'|',string(sTab.regime),'|',string(sTab.snr_db),'|',string(sTab.replicate));
[tf, loc] = ismember(ke, ks);
d = eTab.(valVar)(tf) - sTab.(valVar)(loc(tf));
n = numel(d);
end

function p = local_signrank(d)
% Two-sided Wilcoxon signed-rank p-value via normal approximation with tie correction.
d = d(d~=0); n = numel(d);
if n < 1; p = NaN; return; end
R = local_tiedrank(abs(d));
Wp = sum(R(d>0));
mu = n*(n+1)/4;
% tie correction term
[~,~,g] = unique(abs(d));
tieTerm = 0;
for k = 1:max(g); tk = sum(g==k); tieTerm = tieTerm + (tk^3 - tk); end
sigma = sqrt(n*(n+1)*(2*n+1)/24 - tieTerm/48);
if sigma == 0; p = 1; return; end
z = (Wp - mu) / sigma;
p = erfc(abs(z)/sqrt(2));   % two-sided
end

function r = local_tiedrank(x)
[xs, ord] = sort(x(:)); r = zeros(numel(x),1); i = 1; n = numel(xs);
while i <= n
    j = i; while j < n && xs(j+1)==xs(i); j = j+1; end
    r(ord(i:j)) = (i+j)/2; i = j+1;
end
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_bench_stats_pure.m`. Expected: `ALL TESTS PASSED`.

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/bench_stats.m dev/tests/test_bench_stats_pure.m
git commit -m "feat(benchmark): bench_stats medians/IQR + paired Wilcoxon (toolbox-free)"
```

---

## Task 6: `bench_figures` — five figures

**Files:**
- Create: `dev/benchmarks/bench_figures.m`
- Test: `dev/tests/test_bench_figures_pure.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_bench_figures_pure
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));

% Minimal synthetic CSV spanning the dimensions the figures need.
tmp = fullfile(tempdir, sprintf('bench_fig_test_%d', feature('getpid')));
if ~exist(tmp,'dir'); mkdir(tmp); end
csvPath = fullfile(tmp,'synthetic.csv');
fid = fopen(csvPath,'w');
fprintf(fid,'anatomy,regime,snr_db,replicate,method,K,locerror_mm,correlation,nrmse,auc,spatial_dispersion_mm\n');
regs = {'focal','patch','distributed'}; snrs=[4 10]; meths={'wmne','dspm','sloreta'}; Ks=[600 2000];
for ir=1:3; for is=1:2; for rep=1:4
    for im=1:3
        fprintf(fid,'auditory,%s,%d,%d,%s,NaN,%g,0.9,0.1,0.99,30\n', regs{ir}, snrs(is), rep, meths{im}, 12+rep+im);
    end
    for ik=1:2
        fprintf(fid,'auditory,%s,%d,%d,eig_mne_log,%d,%g,0.9,0.1,0.99,30\n', regs{ir}, snrs(is), rep, Ks(ik), 14+rep-ik);
        fprintf(fid,'auditory,%s,%d,%d,eig_dspm_log,%d,%g,0.9,0.1,0.99,30\n', regs{ir}, snrs(is), rep, Ks(ik), 15+rep-ik);
    end
end; end; end
fclose(fid);

figDir = fullfile(tmp,'figures');
files = bench_figures(csvPath, figDir, []);   % [] => skip the cortex-render figure (no protocol)
for i=1:numel(files)
    assert(exist(files{i},'file')==2, 'figure not written: %s', files{i});
end
assert(numel(files) >= 4, 'expect at least 4 data-driven figures without cortex render.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_bench_figures_pure.m`. Expected: FAIL — `Undefined function 'bench_figures'`.

- [ ] **Step 3: Create `dev/benchmarks/bench_figures.m`**

```matlab
function files = bench_figures(csvPath, figDir, renderCtx)
% BENCH_FIGURES: Produce the five benchmark figures (PNG + .fig).
% USAGE: files = bench_figures(csvPath, figDir, renderCtx)
%   renderCtx : [] to skip the cortex-render figure, or struct with fields
%               .Vertices .Faces .gt .estMaps (struct method->[nV x 1]) .titleStr
if nargin < 3; renderCtx = []; end
if ~exist(figDir,'dir'); mkdir(figDir); end
T = readtable(csvPath);
T.method = string(T.method); T.regime = string(T.regime);
files = {};

eigMethods = ["eig_mne_log","eig_dspm_log"];
stdMethods = ["wmne","dspm","sloreta"];
allMethods = [stdMethods, eigMethods];
Kmax = max(T.K(~isnan(T.K)));
col = lines(numel(allMethods));

% Helper: LocError vector for a method (eig at Kmax) optionally filtered by regime/snr
    function v = locOf(m, regime, snr)
        sel = T.method==m;
        if ismember(m, eigMethods); sel = sel & T.K==Kmax; end
        if nargin>=2 && ~isempty(regime); sel = sel & T.regime==regime; end
        if nargin>=3 && ~isempty(snr);    sel = sel & T.snr_db==snr; end
        v = T.locerror_mm(sel);
    end

% ----- F1: Distribution (strip + median/IQR) of LocError per method -----
f1 = figure('Color','w','Position',[100 100 760 460]); hold on;
for im=1:numel(allMethods)
    v = locOf(allMethods(im));
    x = im + (rand(numel(v),1)-0.5)*0.3;
    scatter(x, v, 14, col(im,:), 'filled', 'MarkerFaceAlpha',0.45);
    med = median(v); q1 = quantile(v,0.25); q3 = quantile(v,0.75);
    plot([im-0.32 im+0.32],[med med],'k-','LineWidth',2);
    plot([im im],[q1 q3],'k-','LineWidth',1);
end
set(gca,'XTick',1:numel(allMethods),'XTickLabel',cellstr(allMethods),'XTickLabelRotation',20);
ylabel('Localization error (mm)'); title('LocError distribution by method (eig at K_{max})'); grid on;
files{end+1} = save_fig(f1, figDir, 'f1_distribution');

% ----- F2: SNR sweep, one panel per regime -----
regimes = unique(T.regime,'stable'); snrs = unique(T.snr_db)';
f2 = figure('Color','w','Position',[100 100 1100 360]);
for ir=1:numel(regimes)
    subplot(1,numel(regimes),ir); hold on;
    for im=1:numel(allMethods)
        mu = arrayfun(@(s) mean(locOf(allMethods(im), regimes(ir), s)), snrs);
        sd = arrayfun(@(s) std(locOf(allMethods(im), regimes(ir), s)),  snrs);
        errorbar(snrs, mu, sd, '-o', 'Color', col(im,:), 'MarkerFaceColor', col(im,:), 'CapSize',4);
    end
    xlabel('SNR (dB)'); ylabel('LocError (mm)'); title(regimes(ir)); grid on;
    if ir==numel(regimes); legend(cellstr(allMethods),'Location','northeastoutside'); end
end
sgtitle('LocError vs SNR by regime');
files{end+1} = save_fig(f2, figDir, 'f2_snr_sweep');

% ----- F4: Per-regime grouped bars (median LocError, IQR whiskers) -----
f4 = figure('Color','w','Position',[100 100 900 460]); hold on;
nM = numel(allMethods); nR = numel(regimes); bw = 0.8/nM;
for im=1:nM
    meds = arrayfun(@(ir) median(locOf(allMethods(im), regimes(ir))), 1:nR);
    q1   = arrayfun(@(ir) quantile(locOf(allMethods(im), regimes(ir)),0.25), 1:nR);
    q3   = arrayfun(@(ir) quantile(locOf(allMethods(im), regimes(ir)),0.75), 1:nR);
    xpos = (1:nR) + (im-(nM+1)/2)*bw;
    bar(xpos, meds, bw*0.9, 'FaceColor', col(im,:), 'EdgeColor','none');
    errorbar(xpos, meds, meds-q1, q3-meds, 'k', 'LineStyle','none', 'CapSize',3);
end
set(gca,'XTick',1:nR,'XTickLabel',cellstr(regimes)); ylabel('Median LocError (mm)');
legend(cellstr(allMethods),'Location','northeastoutside'); title('LocError by regime × method'); grid on;
files{end+1} = save_fig(f4, figDir, 'f4_per_regime');

% ----- F5: K-sweep curve (focal emphasis), eig methods vs std reference lines -----
Kvals = unique(T.K(~isnan(T.K)))';
f5 = figure('Color','w','Position',[100 100 1100 360]);
for ir=1:numel(regimes)
    subplot(1,numel(regimes),ir); hold on;
    for ie=1:numel(eigMethods)
        med = arrayfun(@(k) median(T.locerror_mm(T.method==eigMethods(ie) & T.regime==regimes(ir) & T.K==k)), Kvals);
        plot(Kvals, med, '-o', 'LineWidth',1.5, 'DisplayName', char(eigMethods(ie)));
    end
    for is=1:numel(stdMethods)
        ref = median(locOf(stdMethods(is), regimes(ir)));
        yline(ref, '--', char(stdMethods(is)), 'LabelHorizontalAlignment','left');
    end
    xlabel('Total modes K'); ylabel('Median LocError (mm)'); title(regimes(ir)); grid on;
    if ir==1; legend('Location','northeast'); end
end
sgtitle('K-sweep: eigenmode LocError vs total modes (std = dashed reference)');
files{end+1} = save_fig(f5, figDir, 'f5_ksweep');

% ----- F3: Example reconstructions on cortex (only if renderCtx provided) -----
if ~isempty(renderCtx)
    methodsR = fieldnames(renderCtx.estMaps);
    f3 = figure('Color','w','Position',[100 100 1200 420]);
    npan = 1 + numel(methodsR);
    subplot(1,npan,1);
    render_map(renderCtx.Vertices, renderCtx.Faces, renderCtx.gt); title('ground truth');
    for i=1:numel(methodsR)
        subplot(1,npan,1+i);
        render_map(renderCtx.Vertices, renderCtx.Faces, renderCtx.estMaps.(methodsR{i}));
        title(strrep(methodsR{i},'_','\_'));
    end
    sgtitle(renderCtx.titleStr);
    files{end+1} = save_fig(f3, figDir, 'f3_cortex_reconstruction');
end
end

% ---- helpers ----
function p = save_fig(h, figDir, name)
p = fullfile(figDir, [name '.png']);
print(h, p, '-dpng', '-r150');
savefig(h, fullfile(figDir, [name '.fig']));
close(h);
end

function render_map(V, F, map)
map = abs(map(:)); if max(map)>0; map = map/max(map); end
patch('Vertices',V,'Faces',F,'FaceVertexCData',map,'FaceColor','interp', ...
      'EdgeColor','none','FaceLighting','gouraud'); 
axis equal off; view(0,90); camlight headlight; colormap(hot); caxis([0 1]);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_bench_figures_pure.m`. Expected: `ALL TESTS PASSED` (4 PNGs written; cortex render skipped with `renderCtx=[]`).

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/bench_figures.m dev/tests/test_bench_figures_pure.m
git commit -m "feat(benchmark): bench_figures (distribution, SNR sweep, per-regime, K-sweep, cortex render)"
```

---

## Task 7: `benchmark_eigenmodes` driver + smoke e2e

**Files:**
- Create: `dev/benchmarks/benchmark_eigenmodes.m`
- Test: `dev/tests/test_benchmark_eigenmodes_smoke.m`

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_eigenmodes_smoke
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
if ~brainstorm('status'); brainstorm nogui; end

out = benchmark_eigenmodes('smoke');
assert(exist(fullfile(out,'synthetic.csv'),'file')==2, 'synthetic.csv missing.');
assert(exist(fullfile(out,'stats.csv'),'file')==2, 'stats.csv missing.');
assert(exist(fullfile(out,'REPORT.md'),'file')==2, 'REPORT.md missing.');
png = dir(fullfile(out,'figures','*.png'));
assert(numel(png) >= 4, 'expect at least 4 figures (>=5 if cortex render ran).');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run test to verify it fails**

Run `dev/tests/test_benchmark_eigenmodes_smoke.m`. Expected: FAIL — `Undefined function 'benchmark_eigenmodes'`.

- [ ] **Step 3: Create `dev/benchmarks/benchmark_eigenmodes.m`**

```matlab
function outDir = benchmark_eigenmodes(preset)
% BENCHMARK_EIGENMODES: Top-level driver for the eigenmode accuracy benchmark.
% USAGE:  outDir = benchmark_eigenmodes()         % full
%         outDir = benchmark_eigenmodes('smoke')  % fast preset
if nargin < 1 || isempty(preset); preset = 'full'; end
C = bench_config(preset);
outDir = C.outDir;
if ~exist(outDir,'dir'); mkdir(outDir); end

% 1) Run synthetic benchmark
csvPath = bench_run(C);

% 2) Statistics
St = bench_stats(csvPath, outDir);

% 3) Build one cortex-render context from a representative focal case (anatomy 1)
renderCtx = [];
try
    renderCtx = build_render_ctx(C);
catch ME
    warning('benchmark_eigenmodes: cortex render skipped: %s', ME.message);
end

% 4) Figures
figDir = fullfile(outDir,'figures');
figFiles = bench_figures(csvPath, figDir, renderCtx);

% 5) REPORT.md
write_report(outDir, C, St, figFiles);
fprintf('BENCH> report written to %s\n', fullfile(outDir,'REPORT.md'));
end

function renderCtx = build_render_ctx(C)
anat = C.anatomies{1};
info = bench_fixtures(anat, C.nModes_eig);
baseHM = in_bst_headmodel(info.baseHmFile, 1);
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess_bst(info.surfaceFile);
if ~isfield(Surf,'VertConn')||isempty(Surf.VertConn); Surf.VertConn = tess_vertconn(Surf.Vertices,Surf.Faces); end
NC = load(file_fullpath(info.ncFile)); Cnoise = NC.NoiseCov(goodMask,goodMask);
S = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn), 'focal','nTime',C.nTime,'Seed',C.seed);
Sim = bst_benchmark_simulate(L, S.Sources, Cnoise, 'SNR', 10, 'Seed', C.seed);
Est = bst_benchmark_inverse(Sim.F, info.baseHmFile, info.ncFile, info.chFile, goodMask, 10, max(C.k_total));
[~,tStar] = max(sum(S.Sources.^2,1));
estMaps = struct();
for m = C.methods
    if isfield(Est, m{1}); estMaps.(m{1}) = Est.(m{1})(:,tStar); end
end
renderCtx = struct('Vertices',Surf.Vertices,'Faces',Surf.Faces,'gt',S.GT(:), ...
    'estMaps',estMaps,'titleStr',sprintf('%s — focal source, SNR 10 dB', anat.label));
end

function write_report(outDir, C, St, figFiles)
fid = fopen(fullfile(outDir,'REPORT.md'),'w');
fprintf(fid,'# Eigenmode Source-Mapping Accuracy Benchmark — Report\n\n');
fprintf(fid,'Anatomies: %s. Methods: %s. Regimes: %s. SNR (dB): %s. K (total): %s.\n\n', ...
    strjoin(cellfun(@(a)a.label,C.anatomies,'uni',0),', '), strjoin(C.methods,', '), ...
    strjoin(C.regimes,', '), mat2str(C.snr_db), mat2str(C.k_total));

% Plateau-K finding (focal, eig_mne_log)
fprintf(fid,'## K-sweep (focal, eig\\_mne\\_log)\n\n| K | median LocError (mm) |\n|---|---|\n');
ks = St.ksweep(St.ksweep.regime=="focal" & St.ksweep.eig_method=="eig_mne_log", :);
ks = sortrows(ks,'K');
plateauK = NaN;
for i=1:height(ks)
    fprintf(fid,'| %d | %.2f |\n', ks.K(i), ks.median_locerror_mm(i));
    if i>1 && isnan(plateauK) && (ks.median_locerror_mm(i-1)-ks.median_locerror_mm(i) <= C.plateauTol_mm)
        plateauK = ks.K(i-1);
    end
end
if isnan(plateauK) && ~isempty(ks); plateauK = ks.K(end); end
fprintf(fid,'\nPlateau-K (focal): **%d total modes** (improvement <= %.1f mm beyond this).\n\n', plateauK, C.plateauTol_mm);

% Competitiveness verdict from paired comparisons
fprintf(fid,'## Competitiveness (paired Wilcoxon, eig vs standard)\n\n');
fprintf(fid,'| regime | eig | vs | median diff (mm) | p |\n|---|---|---|---|---|\n');
for i=1:height(St.compare)
    c = St.compare(i,:);
    fprintf(fid,'| %s | %s | %s | %+.2f | %.4g |\n', c.regime, c.eig_method, c.std_method, c.median_diff_mm, c.p_value);
end
fprintf(fid,'\n_Eigenmode is "competitive" where median diff is small and p is not significant, or where the diff favours eig._\n\n');

% Figures
fprintf(fid,'## Figures\n\n');
for i=1:numel(figFiles)
    [~,nm,ext] = fileparts(figFiles{i});
    fprintf(fid,'![%s](figures/%s%s)\n\n', nm, nm, ext);
end
fclose(fid);
end
```

- [ ] **Step 4: Run test to verify it passes**

Run `dev/tests/test_benchmark_eigenmodes_smoke.m`. Expected: `ALL TESTS PASSED` (writes `dev/benchmarks/smoke_run/` with CSVs, ≥4 figures, REPORT.md).

- [ ] **Step 5: Commit**

```bash
git add dev/benchmarks/benchmark_eigenmodes.m dev/tests/test_benchmark_eigenmodes_smoke.m
git commit -m "feat(benchmark): benchmark_eigenmodes driver + smoke e2e + REPORT.md"
```

---

## Task 8: Full run + rename output to dated folder

**Files:**
- No new files; produces `dev/benchmarks/eigenmode_accuracy_run/` artifacts.

- [ ] **Step 1: Execute the full benchmark**

Via MCP `evaluate_matlab_code` (allow several minutes; the 1000-mode eigendecomposition per anatomy is the slow part on first run):

```matlab
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3');
addpath('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks');
out = benchmark_eigenmodes();   % full preset
fprintf('DONE: %s\n', out);
```

Expected: `DONE: .../dev/benchmarks/eigenmode_accuracy_run`; `REPORT.md`, `synthetic.csv`, `stats.csv`, `stats_compare.csv`, `stats_ksweep.csv`, and 5 PNGs present.

- [ ] **Step 2: Sanity-check the outputs**

```matlab
T = readtable('/Users/diellorbasha/workspace/research/code/brainstorm3/dev/benchmarks/eigenmode_accuracy_run/synthetic.csv');
fprintf('rows=%d anatomies=%s methods=%s regimes=%s\n', height(T), ...
    strjoin(unique(string(T.anatomy))',','), strjoin(unique(string(T.method))',','), strjoin(unique(string(T.regime))',','));
assert(all(isfinite(T.locerror_mm)), 'all LocError finite');
assert(numel(unique(T.anatomy))==2, 'both anatomies present');
```

Expected: 2 anatomies, 5 methods, 3 regimes; all finite.

- [ ] **Step 3: Review the five figures and REPORT.md**

Open `dev/benchmarks/eigenmode_accuracy_run/figures/*.png` and `REPORT.md`. Confirm: F5 shows whether focal LocError plateaus with K; the paired-Wilcoxon table populates; F3 cortex render shows ground truth vs methods.

- [ ] **Step 4: Commit the results**

```bash
git add dev/benchmarks/eigenmode_accuracy_run
git commit -m "results(benchmark): full eigenmode accuracy run (synthetic, stats, 5 figures, report)"
```

---

## Self-Review

**1. Spec coverage**

- Synthetic-on-real-cortex track → Tasks 3–4. ✔
- Methods {wmne,dspm,sloreta,eig_mne_log,eig_dspm_log} → `bench_config`, scored in `bench_run`. ✔
- Regimes focal/patch/distributed → `bench_config.regimes` + `regimeOpts`. ✔
- SNR sweep [2 4 6 10 20] → `bench_config.snr_db`, looped in `bench_run`. ✔
- K-sweep {600,1200,2000}≈{300,600,1000}/hemi → Task 1 passthrough + `bench_run` K loop + Task 2 config. ✔
- Metrics LocError/Corr/NRMSE/AUC → recorded by `bench_run` via `bst_benchmark_metrics`. ✔
- Stats: median/IQR + paired Wilcoxon → `bench_stats` (toolbox-free). ✔
- Figures 1–5 (distribution, SNR sweep, cortex reconstruction, per-regime, K-sweep) → `bench_figures`. ✔
- Outputs under dated/named folder + REPORT.md → `benchmark_eigenmodes` + Task 8. ✔
- Reproducibility (deterministic seeds) → `bench_run` seed derivation; `bench_config.seed`. ✔
- Error handling (try/catch, never abort) → `bench_run` inverse try/catch; `bench_figures` render guard. ✔
- Smoke preset → `bench_config('smoke')` + Task 7 e2e. ✔
- Plateau-K selection → `write_report` computes it from the K-sweep medians. ✔

**2. Placeholder scan:** No "TBD/TODO"; every code step is complete and runnable. The cortex-render path is explicitly optional (`renderCtx=[]`), not a placeholder.

**3. Type consistency:** CSV header columns are identical across `bench_run` (writer), `bench_stats` (`readtable`), `bench_figures` (`readtable`), and all test fixtures: `anatomy,regime,snr_db,replicate,method,K,locerror_mm,correlation,nrmse,auc,spatial_dispersion_mm`. `bst_benchmark_inverse`'s new 7th arg `nModes` matches the `bst_inverse_eigenmodes` `'nModes'` option. Method field names (`wmne,dspm,sloreta,eig_mne_log,eig_dspm_log`) match `bst_benchmark_inverse`'s output struct fields exactly.

**Out-of-scope (per spec):** phantom track, real evoked-data validation, prior-family comparison, vector eigenmodes/dispersion — none included. ✔
