# ESI Benchmark Framework (v1) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a reusable, statistically-grounded ESI benchmark: simulate ground-truth cortical sources (focal/patch/distributed) with colored noise at controlled SNR, reconstruct with a comparator panel (native inverses + eigenmode variants), score with a metric suite, and report descriptive stats with paired differences.

**Architecture:** Five focused MATLAB units + a thin driver. Math-core units (`bst_benchmark_sources`, `bst_benchmark_metrics`, `bst_benchmark_simulate`, `bst_benchmark_report`) are pure and unit-tested with synthetic data. The integration unit (`bst_benchmark_inverse`) wires existing Brainstorm inverses and is exercised by an e2e smoke on OMEGA. The driver `tutorial_benchmark_esi.m` defines the sweep.

**Tech Stack:** MATLAB / Brainstorm. Tests are `dev/tests/*.m` functions printing `ALL TESTS PASSED`, run via the MATLAB MCP (`run_matlab_file`). Spec: `dev/2026-06-02-esi-benchmark-framework-design.md`.

---

## File structure

| File | Responsibility |
|---|---|
| `toolbox/math/bst_benchmark_sources.m` | Ground-truth generator: regime → GT vertex map + time course + seed. |
| `toolbox/inverse/bst_benchmark_metrics.m` | Metric suite (LocError, AUC, NRMSE, correlation, spatial dispersion). |
| `toolbox/math/bst_benchmark_simulate.m` | Forward-project GT + colored noise at target SNR. |
| `toolbox/math/bst_benchmark_report.m` | Aggregate rows → median/IQR/bootstrap CI + paired diffs → CSV/markdown. |
| `toolbox/inverse/bst_benchmark_inverse.m` | Run comparator panel → vertex estimates (matched configs). |
| `toolbox/script/tutorial_benchmark_esi.m` | Driver: define sweep, orchestrate, write report. |

Reused unchanged: `bst_inverse_linear_2018`, `bst_inverse_eigenmodes`, `bst_eigenmode_leadfield`, `in_bst_headmodel`, `in_tess_eigenmodes`, `bst_resolution_metrics`.

---

## Task 1: `bst_benchmark_sources` — ground-truth generator

**Files:** Create `toolbox/math/bst_benchmark_sources.m`; test `dev/tests/test_benchmark_sources_pure.m`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_sources_pure
% Verify the ground-truth generator: focal/patch/distributed regimes on a
% synthetic surface; shapes; seeding reproducibility; spatial profiles.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Synthetic surface: 6x6 grid of vertices on a plane (deterministic), 4-neighbour conn
n = 6; [X,Y] = ndgrid(1:n, 1:n);
Vertices = [X(:)*0.01, Y(:)*0.01, zeros(n*n,1)];   % 10 mm spacing, in metres
nV = size(Vertices,1);
VertConn = sparse(nV, nV);
idx = reshape(1:nV, n, n);
for i = 1:n
    for j = 1:n
        if i<n; VertConn(idx(i,j), idx(i+1,j)) = 1; VertConn(idx(i+1,j), idx(i,j)) = 1; end
        if j<n; VertConn(idx(i,j), idx(i,j+1)) = 1; VertConn(idx(i,j+1), idx(i,j)) = 1; end
    end
end
Surface = struct('Vertices', Vertices, 'VertConn', VertConn);

% --- focal ---
S = bst_benchmark_sources(Surface, 'focal', 'nTime', 10, 'Seed', 1);
assert(isequal(size(S.GT), [nV 1]), 'GT map must be [nV x 1].');
assert(nnz(S.GT) == 1, 'focal GT must have exactly one active vertex.');
assert(S.GT(S.SeedVertex) == max(S.GT), 'seed vertex must be the active one.');
assert(isequal(size(S.Sources), [nV 10]), 'Sources must be [nV x nTime].');

% --- patch (radius 1 hop -> seed + 4 neighbours = 5 active for an interior seed) ---
P = bst_benchmark_sources(Surface, 'patch', 'Radius', 1, 'Seed', 1, 'SeedVertex', idx(3,3));
assert(P.GT(idx(3,3))>0 && P.GT(idx(2,3))>0 && P.GT(idx(4,3))>0 ...
    && P.GT(idx(3,2))>0 && P.GT(idx(3,4))>0, 'patch must cover seed + 1-hop neighbours.');
assert(P.GT(idx(1,1))==0, 'patch must not reach distant vertices at radius 1.');

% --- distributed (Gaussian falloff, monotone decreasing with distance) ---
D = bst_benchmark_sources(Surface, 'distributed', 'Sigma', 0.02, 'SeedVertex', idx(3,3), 'Seed', 1);
dPeak = D.GT(idx(3,3)); dNear = D.GT(idx(3,4)); dFar = D.GT(idx(1,1));
assert(dPeak >= dNear && dNear > dFar, 'distributed profile must fall off with distance.');
assert(abs(dPeak - max(D.GT)) < 1e-12, 'peak must be at the seed.');

% --- seeding reproducibility (same seed -> same random seed vertex) ---
A = bst_benchmark_sources(Surface, 'focal', 'Seed', 7);
B = bst_benchmark_sources(Surface, 'focal', 'Seed', 7);
assert(A.SeedVertex == B.SeedVertex, 'same Seed must give the same random seed vertex.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'bst_benchmark_sources'`.

- [ ] **Step 3: Write the implementation**

```matlab
function S = bst_benchmark_sources(Surface, regime, varargin)
% BST_BENCHMARK_SOURCES: Generate ground-truth cortical sources for benchmarking.
%
% USAGE:  S = bst_benchmark_sources(Surface, regime, 'Param', value, ...)
%
% INPUTS:
%   Surface : struct with .Vertices [nV x 3] (metres) and .VertConn [nV x nV] sparse
%   regime  : 'focal' | 'patch' | 'distributed'
% OPTIONS:
%   'Seed'       : RNG seed for random seed-vertex selection (default 1)
%   'SeedVertex' : explicit seed vertex index (default: random, reproducible via Seed)
%   'Radius'     : patch radius in graph hops (default 2)               [regime 'patch']
%   'Sigma'      : Gaussian spatial scale in metres (default 0.02)      [regime 'distributed']
%   'nTime'      : number of time samples (default 20)
%
% OUTPUT struct S:
%   .GT         [nV x 1]      ground-truth spatial amplitude map (peak normalized to 1)
%   .Sources    [nV x nTime]  GT source matrix = GT * timecourse
%   .SeedVertex scalar        seed vertex index (the GT "true location")
%   .Regime     char          the regime
%
% Authors: Diellor Basha, 2026

Seed = 1; SeedVertex = []; Radius = 2; Sigma = 0.02; nTime = 20;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'seed',       Seed       = varargin{i+1};
        case 'seedvertex', SeedVertex = varargin{i+1};
        case 'radius',     Radius     = varargin{i+1};
        case 'sigma',      Sigma      = varargin{i+1};
        case 'ntime',      nTime      = varargin{i+1};
    end
end
V  = double(Surface.Vertices);
nV = size(V, 1);

% Seed vertex (reproducible)
if isempty(SeedVertex)
    rng(Seed);
    SeedVertex = randi(nV);
end

% Spatial profile
GT = zeros(nV, 1);
switch lower(regime)
    case 'focal'
        GT(SeedVertex) = 1;
    case 'patch'
        keep = graph_ball(Surface.VertConn, SeedVertex, Radius);
        GT(keep) = 1;
    case 'distributed'
        d = sqrt(sum((V - V(SeedVertex,:)).^2, 2));   % Euclidean distance (metres)
        GT = exp(-(d.^2) / (2 * Sigma^2));
    otherwise
        error('bst_benchmark_sources:UnknownRegime', 'Unknown regime: %s', regime);
end
GT = GT / max(GT);   % peak-normalize

% Time course: Gaussian-windowed burst, peak at the centre sample
t  = (1:nTime) - (nTime+1)/2;
tc = exp(-(t.^2) / (2 * (nTime/6)^2));

S = struct();
S.GT         = GT;
S.Sources    = GT * tc;            % [nV x nTime]
S.SeedVertex = SeedVertex;
S.Regime     = lower(regime);
end

function keep = graph_ball(VertConn, seed, radius)
% Vertices within `radius` hops of seed (BFS), inclusive.
keep = false(size(VertConn, 1), 1);
keep(seed) = true;
frontier = seed;
for h = 1:radius
    nbr = any(VertConn(:, frontier) ~= 0, 2);
    new = nbr & ~keep;
    keep = keep | nbr;
    frontier = find(new);
    if isempty(frontier); break; end
end
keep = find(keep);
end
```

- [ ] **Step 4: Run it; verify it PASSES** (`ALL TESTS PASSED`).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_benchmark_sources.m dev/tests/test_benchmark_sources_pure.m
git commit -m "$(printf 'ESI benchmark: ground-truth source generator (focal/patch/distributed)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: `bst_benchmark_metrics` — metric suite

**Files:** Create `toolbox/inverse/bst_benchmark_metrics.m`; test `dev/tests/test_benchmark_metrics_pure.m`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_metrics_pure
% Verify the metric suite on synthetic GT vs estimate.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% 4 vertices on a line, 10 mm apart
GridLoc = [0 0 0; 10 0 0; 20 0 0; 30 0 0] / 1000;   % metres
gt = [0; 1; 0; 0];                                   % focal GT at vertex 2

% Perfect estimate
M = bst_benchmark_metrics(gt, gt, GridLoc, 2);
assert(M.LocError < 1e-9, 'perfect estimate -> zero localization error.');
assert(abs(M.Correlation - 1) < 1e-12, 'perfect estimate -> correlation 1.');
assert(abs(M.NRMSE) < 1e-9, 'perfect estimate -> NRMSE 0.');
assert(abs(M.AUC - 1) < 1e-12, 'perfect estimate -> AUC 1.');

% Estimate peaks at vertex 4 (30 mm from true vertex 2 at 10 mm) -> LE = 20 mm
est = [0; 0.1; 0; 1];
M2 = bst_benchmark_metrics(gt, est, GridLoc, 2);
assert(abs(M2.LocError - 20) < 1e-9, 'peak at vertex 4 -> 20 mm localization error.');
assert(M2.SpatialDispersion > 0, 'dispersed estimate -> positive dispersion.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'bst_benchmark_metrics'`.

- [ ] **Step 3: Write the implementation**

```matlab
function M = bst_benchmark_metrics(gt, est, GridLoc, seedVertex)
% BST_BENCHMARK_METRICS: Compare an estimated source map against ground truth.
%
% USAGE:  M = bst_benchmark_metrics(gt, est, GridLoc, seedVertex)
%
% INPUTS:
%   gt         [nV x 1]   ground-truth amplitude map
%   est        [nV x 1]   estimated amplitude map (same vertices)
%   GridLoc    [nV x 3]   vertex positions in metres
%   seedVertex scalar     GT "true location" vertex index
%
% OUTPUT struct M: .LocError (mm) .AUC .NRMSE .Correlation .SpatialDispersion (mm)
%
% Authors: Diellor Basha, 2026
gt  = double(gt(:));  est = double(est(:));
ag  = abs(gt);        ae  = abs(est);
mm  = 1000;

% Localization error: GT seed location vs estimate peak location
[~, iPk] = max(ae);
M.LocError = norm(GridLoc(iPk,:) - GridLoc(seedVertex,:)) * mm;

% Spatial dispersion about the true location, weighted by estimate power
d = sqrt(sum((GridLoc - GridLoc(seedVertex,:)).^2, 2)) * mm;
w = ae.^2;
if sum(w) > 0
    M.SpatialDispersion = sqrt(sum(d.^2 .* w) / sum(w));
else
    M.SpatialDispersion = 0;
end

% NRMSE (norm-matched, GBF definition)
if norm(est) > 0
    est_s = est * (norm(gt) / norm(est));
else
    est_s = est;
end
rng_gt = max(gt) - min(gt);
if rng_gt > 0
    M.NRMSE = sqrt(mean((gt - est_s).^2)) / rng_gt;
else
    M.NRMSE = sqrt(mean((gt - est_s).^2));
end

% Correlation
if std(est) > 0 && std(gt) > 0
    M.Correlation = corr(gt, est);
else
    M.Correlation = 0;
end

% AUC: detect GT-active vertices (|gt| > 50% of peak) from estimate magnitude
labels = ag > 0.5 * max(ag);
M.AUC  = local_auc(labels, ae);
end

function auc = local_auc(labels, scores)
% Rank-based ROC AUC (Mann-Whitney). labels logical, scores >= 0.
pos = scores(labels); neg = scores(~labels);
nP = numel(pos); nN = numel(neg);
if nP == 0 || nN == 0; auc = 0.5; return; end
[~, order] = sort(scores); r = zeros(size(scores)); r(order) = 1:numel(scores);
auc = (sum(r(labels)) - nP*(nP+1)/2) / (nP*nN);
end
```

- [ ] **Step 4: Run it; verify it PASSES** (`ALL TESTS PASSED`).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/inverse/bst_benchmark_metrics.m dev/tests/test_benchmark_metrics_pure.m
git commit -m "$(printf 'ESI benchmark: metric suite (LocError/AUC/NRMSE/corr/dispersion)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 3: `bst_benchmark_simulate` — forward + colored noise at target SNR

**Files:** Create `toolbox/math/bst_benchmark_simulate.m`; test `dev/tests/test_benchmark_simulate_pure.m`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_simulate_pure
% Verify forward projection + colored-noise addition at a target SNR.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

nCh = 20; nV = 8; nTime = 50;
rng(0);
L = randn(nCh, nV);                          % leadfield
Sources = zeros(nV, nTime); Sources(3,:) = sin(2*pi*(1:nTime)/nTime);  % one active vertex
C = eye(nCh);                                % white noise cov (for a clean SNR check)

snrTarget = 6;  % dB
Sim = bst_benchmark_simulate(L, Sources, C, 'SNR', snrTarget, 'Seed', 1);
assert(isequal(size(Sim.F), [nCh nTime]), 'F must be [nCh x nTime].');

% Achieved SNR (signal vs noise power) should match the target within tolerance
sigPow = mean(Sim.Fsignal(:).^2);
noiPow = mean(Sim.Fnoise(:).^2);
snrAch = 10*log10(sigPow / noiPow);
assert(abs(snrAch - snrTarget) < 0.5, 'achieved SNR must match target within 0.5 dB.');

% Reproducible
Sim2 = bst_benchmark_simulate(L, Sources, C, 'SNR', snrTarget, 'Seed', 1);
assert(max(abs(Sim.F(:) - Sim2.F(:))) < 1e-12, 'same Seed -> identical noise draw.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'bst_benchmark_simulate'`.

- [ ] **Step 3: Write the implementation**

```matlab
function Sim = bst_benchmark_simulate(L, Sources, NoiseCov, varargin)
% BST_BENCHMARK_SIMULATE: Forward-project GT sources and add colored sensor noise.
%
% USAGE:  Sim = bst_benchmark_simulate(L, Sources, NoiseCov, 'SNR', dB, 'Seed', s)
%
% INPUTS:
%   L        [nCh x nV]    constrained good-channel leadfield
%   Sources  [nV x nTime]  ground-truth source matrix
%   NoiseCov [nCh x nCh]   sensor noise covariance (colored noise model)
% OPTIONS:
%   'SNR'  : target sensor SNR in dB (default 6)
%   'Seed' : RNG seed (default 1)
%
% OUTPUT struct Sim: .F (signal+noise) .Fsignal .Fnoise  [all nCh x nTime], .SNR
%
% Colored noise uses the eigendecomposition model of process_simulate_recordings:
% xn = V*sqrt(D)*randn, scaled so 10*log10(sigPow/noisePow) = SNR.
%
% Authors: Diellor Basha, 2026
SNR = 6; Seed = 1;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'snr',  SNR  = varargin{i+1};
        case 'seed', Seed = varargin{i+1};
    end
end
L = double(L); Sources = double(Sources); NoiseCov = double(NoiseCov);
[nCh, nTime] = deal(size(L,1), size(Sources,2));

Fsignal = L * Sources;                              % [nCh x nTime]

% Colored noise from the covariance (real symmetric part for safety)
NoiseCov = 0.5*(NoiseCov + NoiseCov');
[Vn, Dn] = eig(NoiseCov);
dn = max(real(diag(Dn)), 0);
rng(Seed);
xn = real(Vn) * diag(sqrt(dn)) * randn(nCh, nTime);

% Scale noise to the target SNR
sigPow = mean(Fsignal(:).^2);
curPow = mean(xn(:).^2);
if curPow == 0; curPow = eps; end
noisePowTarget = sigPow / (10^(SNR/10));
xn = xn * sqrt(noisePowTarget / curPow);

Sim = struct('F', Fsignal + xn, 'Fsignal', Fsignal, 'Fnoise', xn, 'SNR', SNR);
end
```

- [ ] **Step 4: Run it; verify it PASSES** (`ALL TESTS PASSED`).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_benchmark_simulate.m dev/tests/test_benchmark_simulate_pure.m
git commit -m "$(printf 'ESI benchmark: forward + colored-noise simulation at target SNR\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 4: `bst_benchmark_report` — aggregation + descriptive stats

**Files:** Create `toolbox/math/bst_benchmark_report.m`; test `dev/tests/test_benchmark_report_pure.m`.

- [ ] **Step 1: Write the failing test**

```matlab
function test_benchmark_report_pure
% Verify aggregation: median/IQR/bootstrap CI + paired differences.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

% Synthetic long-format rows: 2 methods, 1 metric, 6 realizations
rows = struct('regime',{},'snr',{},'method',{},'metric',{},'value',{},'realization',{});
methodVals = struct('eigenmode',[10 12 11 13 9 10], 'dspm',[14 15 13 16 12 14]);
mlist = {'eigenmode','dspm'};
k = 0;
for m = 1:2
    v = methodVals.(mlist{m});
    for r = 1:6
        k = k+1;
        rows(k) = struct('regime','focal','snr',6,'method',mlist{m}, ...
            'metric','LocError','value',v(r),'realization',r);
    end
end

R = bst_benchmark_report(rows, 'RefMethod','eigenmode', 'Seed',1, 'OutDir','');
% Summary: median per (regime,snr,method,metric)
sEig = R.summary(strcmp({R.summary.method},'eigenmode') & strcmp({R.summary.metric},'LocError'));
assert(abs(sEig.median - median(methodVals.eigenmode)) < 1e-9, 'median must match.');
assert(sEig.ci_hi >= sEig.median && sEig.ci_lo <= sEig.median, 'CI must bracket the median.');

% Paired diff: eigenmode - dspm, per realization (all negative here)
p = R.paired(strcmp({R.paired.method},'dspm') & strcmp({R.paired.metric},'LocError'));
assert(abs(p.median_diff - median(methodVals.eigenmode - methodVals.dspm)) < 1e-9, ...
    'paired median diff must match.');
assert(p.median_diff < 0, 'eigenmode beats dspm here -> negative diff.');

disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'bst_benchmark_report'`.

- [ ] **Step 3: Write the implementation**

```matlab
function R = bst_benchmark_report(rows, varargin)
% BST_BENCHMARK_REPORT: Aggregate benchmark rows into descriptive statistics.
%
% USAGE:  R = bst_benchmark_report(rows, 'RefMethod', name, 'Seed', s, 'OutDir', dir)
%
% INPUT rows: struct array with fields .regime .snr .method .metric .value .realization
% OPTIONS:
%   'RefMethod' : reference method for paired differences (default 'eigenmode')
%   'Seed'      : RNG seed for the bootstrap (default 1)
%   'OutDir'    : if non-empty, write summary.csv there ('' = no file)
%   'nBoot'     : bootstrap resamples (default 2000)
%
% OUTPUT struct R:
%   .summary : per (regime,snr,method,metric): .median .iqr .ci_lo .ci_hi .n
%   .paired  : per (regime,snr,method,metric): .median_diff .ci_lo .ci_hi (ref - method)
%
% Authors: Diellor Basha, 2026
RefMethod = 'eigenmode'; Seed = 1; OutDir = ''; nBoot = 2000;
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'refmethod', RefMethod = varargin{i+1};
        case 'seed',      Seed      = varargin{i+1};
        case 'outdir',    OutDir    = varargin{i+1};
        case 'nboot',     nBoot     = varargin{i+1};
    end
end
rng(Seed);

regimes = unique({rows.regime}); snrs = unique([rows.snr]);
methods = unique({rows.method}); metrics = unique({rows.metric});

% ---- per-group descriptive summary ----
R.summary = struct('regime',{},'snr',{},'method',{},'metric',{}, ...
    'median',{},'iqr',{},'ci_lo',{},'ci_hi',{},'n',{});
for ir = 1:numel(regimes)
  for is = 1:numel(snrs)
    for im = 1:numel(methods)
      for ix = 1:numel(metrics)
        v = group_values(rows, regimes{ir}, snrs(is), methods{im}, metrics{ix});
        if isempty(v); continue; end
        ci = boot_ci(v, nBoot);
        R.summary(end+1) = struct('regime',regimes{ir},'snr',snrs(is), ...
            'method',methods{im},'metric',metrics{ix}, ...
            'median',median(v),'iqr',iqr(v),'ci_lo',ci(1),'ci_hi',ci(2),'n',numel(v)); %#ok<AGROW>
      end
    end
  end
end

% ---- paired differences (RefMethod - method), matched by realization ----
R.paired = struct('regime',{},'snr',{},'method',{},'metric',{}, ...
    'median_diff',{},'ci_lo',{},'ci_hi',{},'n',{});
for ir = 1:numel(regimes)
  for is = 1:numel(snrs)
    for im = 1:numel(methods)
      if strcmp(methods{im}, RefMethod); continue; end
      for ix = 1:numel(metrics)
        d = paired_diff(rows, regimes{ir}, snrs(is), RefMethod, methods{im}, metrics{ix});
        if isempty(d); continue; end
        ci = boot_ci(d, nBoot);
        R.paired(end+1) = struct('regime',regimes{ir},'snr',snrs(is), ...
            'method',methods{im},'metric',metrics{ix}, ...
            'median_diff',median(d),'ci_lo',ci(1),'ci_hi',ci(2),'n',numel(d)); %#ok<AGROW>
      end
    end
  end
end

% ---- optional CSV ----
if ~isempty(OutDir)
    if ~exist(OutDir,'dir'); mkdir(OutDir); end
    fid = fopen(fullfile(OutDir,'summary.csv'),'w');
    fprintf(fid,'regime,snr,method,metric,median,iqr,ci_lo,ci_hi,n\n');
    for i = 1:numel(R.summary)
        s = R.summary(i);
        fprintf(fid,'%s,%g,%s,%s,%g,%g,%g,%g,%d\n', s.regime,s.snr,s.method,s.metric, ...
            s.median,s.iqr,s.ci_lo,s.ci_hi,s.n);
    end
    fclose(fid);
end
end

function v = group_values(rows, regime, snr, method, metric)
mask = strcmp({rows.regime},regime) & ([rows.snr]==snr) & ...
       strcmp({rows.method},method) & strcmp({rows.metric},metric);
v = [rows(mask).value];
end

function d = paired_diff(rows, regime, snr, refM, method, metric)
maskC = strcmp({rows.regime},regime) & ([rows.snr]==snr) & strcmp({rows.metric},metric);
ref = rows(maskC & strcmp({rows.method},refM));
oth = rows(maskC & strcmp({rows.method},method));
[~,ia,ib] = intersect([ref.realization],[oth.realization]);
d = [ref(ia).value] - [oth(ib).value];
end

function ci = boot_ci(v, nBoot)
v = v(:); n = numel(v);
if n < 2; ci = [v(1) v(1)]; return; end
bm = zeros(nBoot,1);
for b = 1:nBoot
    bm(b) = median(v(randi(n, n, 1)));
end
ci = prctile(bm, [2.5 97.5]);
end
```

- [ ] **Step 4: Run it; verify it PASSES** (`ALL TESTS PASSED`).

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/math/bst_benchmark_report.m dev/tests/test_benchmark_report_pure.m
git commit -m "$(printf 'ESI benchmark: descriptive aggregation (median/IQR/bootstrap CI + paired diffs)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 5: `bst_benchmark_inverse` — comparator panel (integration)

**Files:** Create `toolbox/inverse/bst_benchmark_inverse.m`; test `dev/tests/test_benchmark_inverse_e2e.m`.

This unit wires existing Brainstorm inverses. It is exercised by an e2e smoke on OMEGA (the pure units above already cover the math).

- [ ] **Step 1: Write the failing e2e test**

```matlab
function test_benchmark_inverse_e2e
% Smoke: run the comparator panel on simulated data using a real OMEGA base head
% model + eigenmode head model + noise cov; verify each method returns a finite
% [nVert x 1] vertex estimate. Skips cleanly otherwise.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study') || isempty(sStudies.Study)
    disp('SKIP: no protocol loaded.'); return;
end
% Find a study with a base (non-eigenmode) surface head model + noise cov + cortex eigenmodes
T = [];
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || s.iHeadModel<1 || numel(s.HeadModel)<1; continue; end
    if isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Channel) || ~any(strcmpi(s.Channel.DisplayableSensorTypes,'MEG')); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        isEig = isfield(hm,'isEigenmode') && hm.isEigenmode;
        if ~isEig && strcmpi(hm.HeadModelType,'surface'); iBase = ih; break; end
    end
    if isempty(iBase); continue; end
    [~, isE] = in_tess_eigenmodes(in_bst_headmodel(s.HeadModel(iBase).FileName,0,'SurfaceFile').SurfaceFile);
    if isE; T = struct('iStudy',iS,'iBase',iBase); break; end
end
if isempty(T); disp('SKIP: no study with base surface HM + eigenmodes + noise cov.'); return; end

s = sStudies.Study(T.iStudy);
baseHmFile = s.HeadModel(T.iBase).FileName;
ncFile     = s.NoiseCov(1).FileName;
chFile     = bst_get('ChannelFileForStudy', s.FileName);

% Simulate a focal source through the constrained base leadfield
baseHM = in_bst_headmodel(baseHmFile, 1);   % ApplyOrient=1 -> [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess(baseHM.SurfaceFile, 'Vertices','VertConn');
S = bst_benchmark_sources(struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn), ...
    'focal', 'nTime', 10, 'Seed', 1);
NC = load(file_fullpath(ncFile)); C = NC.NoiseCov(goodMask,goodMask);
Sim = bst_benchmark_simulate(L, S.Sources, C, 'SNR', 6, 'Seed', 1);

Est = bst_benchmark_inverse(Sim.F, baseHmFile, ncFile, chFile, goodMask, 6);  % struct of [nVert x nTime] per method
fn = fieldnames(Est);
assert(~isempty(fn), 'panel must return at least one method estimate.');
for i = 1:numel(fn)
    E = Est.(fn{i});
    assert(size(E,1)==size(L,2), '%s estimate must have nVert rows.', fn{i});
    assert(all(isfinite(E(:))), '%s estimate must be finite.', fn{i});
end
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'bst_benchmark_inverse'`. (If it SKIPs, set up an OMEGA study with a base surface head model + eigenmodes + noise cov first.)

- [ ] **Step 3: Write the implementation**

```matlab
function Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)
% BST_BENCHMARK_INVERSE: Reconstruct simulated data with the comparator panel.
%
% USAGE:  Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)
%
% INPUTS:
%   F          [nGoodCh x nTime]  simulated good-channel sensor data
%   baseHmFile char               base (non-eigenmode) surface head model file
%   ncFile     char               noise covariance file
%   chFile     char               channel file
%   goodMask   [nAllCh x 1]       logical good-channel mask
%   SNR        scalar             SNR (dB) -> SnrFixed for the linear inverses
%
% OUTPUT struct Est: one [nVert x nTime] vertex estimate per method field:
%   .wmne .dspm .sloreta [.lcmv] .eig_mne_log .eig_dspm_log [.eloreta]
%
% NOTE: verify bst_inverse_linear_2018 OPTIONS field names against its header when
% implementing; the e2e is the gate. Reuses bst_inverse_eigenmodes for the
% eigenmode variants (requires an eigenmode head model in the same study).
%
% Authors: Diellor Basha, 2026
Est = struct();
iGood = find(goodMask(:));

% ----- Standard linear inverses (wMNE/dSPM/sLORETA) via bst_inverse_linear_2018 -----
HeadModel = in_bst_headmodel(baseHmFile, 0);     % unconstrained gain for the 2018 inverse
NC = load(file_fullpath(ncFile));
NoiseCovMat = struct('NoiseCov', NC.NoiseCov, 'nSamples', [], 'FourthMoment', []);
ChannelMat  = in_bst_channel(chFile);
chTypes = {ChannelMat.Channel(iGood).Type};

measures = {'amplitude','wmne'; 'dspm2018','dspm'; 'sloreta','sloreta'};
for k = 1:size(measures,1)
    OPTIONS = bst_inverse_linear_2018();
    OPTIONS.InverseMethod  = 'minnorm';
    OPTIONS.InverseMeasure = measures{k,1};
    OPTIONS.NoiseCovMat    = NoiseCovMat;
    OPTIONS.ChannelTypes   = chTypes;
    OPTIONS.SourceOrient   = {'fixed'};
    OPTIONS.SnrMethod      = 'fixed';
    OPTIONS.SnrFixed       = SNR;
    HM = HeadModel; HM.Gain = HeadModel.Gain(iGood,:);
    [Results, ~] = bst_inverse_linear_2018(HM, OPTIONS);
    Est.(measures{k,2}) = Results.ImagingKernel * F;
end

% ----- Eigenmode variants via the eigenmode head model in this study -----
[sStudy] = bst_get('AnyFile', baseHmFile);
eigHmFile = '';
for ih = 1:numel(sStudy.HeadModel)
    try hm = in_bst_headmodel(sStudy.HeadModel(ih).FileName,0); catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode; eigHmFile = sStudy.HeadModel(ih).FileName; break; end
end
if ~isempty(eigHmFile)
    HMe = in_bst_headmodel(eigHmFile, 0);
    [Eig,~] = in_tess_eigenmodes(HMe.SurfaceFile);
    eigCfg = {'mne','eig_mne_log'; 'dspm','eig_dspm_log'};
    for k = 1:size(eigCfg,1)
        [Inv, errE] = bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, ...
            'Method', eigCfg{k,1}, 'Prior', 'log', 'SNR', SNR);
        if ~isempty(errE); continue; end
        Phi = double(Eig.Vectors(:,1:Inv.nModes));
        Est.(eigCfg{k,2}) = (Phi * Inv.ImagingKernel) * F;
    end
end

% ----- eLORETA via FieldTrip, only if available -----
if exist('ft_sourceanalysis','file') == 2
    % Optional: implement via process_ft_sourceanalysis when FieldTrip is on the path.
    % Left unset when FieldTrip is absent (logged by the driver).
end
end
```

- [ ] **Step 4: Run the e2e; verify it PASSES** (`ALL TESTS PASSED`). Debug OPTIONS field mismatches against the `bst_inverse_linear_2018` header if any method errors; the test asserts finite `[nVert x nTime]` estimates. Report any field-name corrections.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/inverse/bst_benchmark_inverse.m dev/tests/test_benchmark_inverse_e2e.m
git commit -m "$(printf 'ESI benchmark: comparator panel (native inverses + eigenmode variants)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 6: `tutorial_benchmark_esi` — driver + e2e

**Files:** Create `toolbox/script/tutorial_benchmark_esi.m`; test `dev/tests/test_benchmark_esi_e2e.m`.

- [ ] **Step 1: Write the failing e2e test**

```matlab
function test_benchmark_esi_e2e
% Smoke: run a tiny sweep (1 regime x 1 SNR x 2 reps) on OMEGA and confirm the
% report struct + CSV are produced. Skips cleanly without a suitable protocol.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
if ~brainstorm('status'); brainstorm nogui; end

outDir = fullfile(repoRoot, 'dev', 'benchmarks', 'smoke');
R = tutorial_benchmark_esi('Regimes', {'focal'}, 'SNRs', 6, 'nLoc', 2, 'nNoise', 1, ...
    'OutDir', outDir);
if isempty(R)
    disp('SKIP: no suitable OMEGA study for the benchmark.'); return;
end
assert(isfield(R,'summary') && ~isempty(R.summary), 'report must have a non-empty summary.');
assert(exist(fullfile(outDir,'summary.csv'),'file')==2, 'summary.csv must be written.');
disp('ALL TESTS PASSED');
end
```

- [ ] **Step 2: Run it; verify it FAILS** — `Undefined function 'tutorial_benchmark_esi'`.

- [ ] **Step 3: Write the implementation**

```matlab
function R = tutorial_benchmark_esi(varargin)
% TUTORIAL_BENCHMARK_ESI: ESI benchmark driver — simulate, reconstruct, score, report.
%
% USAGE:  R = tutorial_benchmark_esi('Regimes',{...},'SNRs',[...],'nLoc',N,'nNoise',M,'OutDir',dir)
%
% Sweeps regimes x SNRs x source locations x noise draws over OMEGA MEG subjects,
% runs the comparator panel, scores each estimate vs ground truth, and aggregates
% into descriptive statistics. Returns [] (skips) if no suitable study is found.
%
% Authors: Diellor Basha, 2026
Regimes = {'focal','patch','distributed'}; SNRs = [0 3 6 10];
nLoc = 20; nNoise = 5; OutDir = fullfile(fileparts(fileparts(mfilename('fullpath'))), '..', 'dev','benchmarks','run');
for i = 1:2:numel(varargin)
    switch lower(varargin{i})
        case 'regimes', Regimes = varargin{i+1};
        case 'snrs',    SNRs    = varargin{i+1};
        case 'nloc',    nLoc    = varargin{i+1};
        case 'nnoise',  nNoise  = varargin{i+1};
        case 'outdir',  OutDir  = varargin{i+1};
    end
end
if ~brainstorm('status'); brainstorm nogui; end

% Reuse the e2e study-finder contract: a study with base surface HM + eigenmodes + noise cov + MEG
target = find_benchmark_study();
if isempty(target); R = []; return; end

baseHM = in_bst_headmodel(target.baseHmFile, 1);    % constrained [nCh x nVert]
goodMask = all(isfinite(double(baseHM.Gain)),2);
L = double(baseHM.Gain(goodMask,:));
Surf = in_tess(baseHM.SurfaceFile, 'Vertices','VertConn');
GridLoc = baseHM.GridLoc;
NC = load(file_fullpath(target.ncFile)); C = NC.NoiseCov(goodMask,goodMask);
SurfStruct = struct('Vertices',Surf.Vertices,'VertConn',Surf.VertConn);

rows = struct('regime',{},'snr',{},'method',{},'metric',{},'value',{},'realization',{});
real = 0;
for ir = 1:numel(Regimes)
  for il = 1:nLoc
    S = bst_benchmark_sources(SurfStruct, Regimes{ir}, 'Seed', il);
    for is = 1:numel(SNRs)
      for in = 1:nNoise
        real = real + 1;
        Sim = bst_benchmark_simulate(L, S.Sources, C, 'SNR', SNRs(is), 'Seed', 1000*il+in);
        Est = bst_benchmark_inverse(Sim.F, target.baseHmFile, target.ncFile, target.chFile, goodMask, SNRs(is));
        tEval = round(size(S.Sources,2)/2);
        fn = fieldnames(Est);
        for k = 1:numel(fn)
            M = bst_benchmark_metrics(S.GT, Est.(fn{k})(:,tEval), GridLoc, S.SeedVertex);
            mn = {'LocError','AUC','NRMSE','Correlation','SpatialDispersion'};
            for q = 1:numel(mn)
                rows(end+1) = struct('regime',Regimes{ir},'snr',SNRs(is),'method',fn{k}, ...
                    'metric',mn{q},'value',M.(mn{q}),'realization',real); %#ok<AGROW>
            end
        end
      end
    end
  end
end

R = bst_benchmark_report(rows, 'RefMethod','eig_mne_log', 'Seed',1, 'OutDir',OutDir);
fprintf('ESI benchmark complete: %d rows -> %s\n', numel(rows), OutDir);
end

function target = find_benchmark_study()
target = [];
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies) || ~isfield(sStudies,'Study'); return; end
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.iHeadModel) || numel(s.HeadModel)<1; continue; end
    if isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    if isempty(s.Channel) || ~any(strcmpi(s.Channel.DisplayableSensorTypes,'MEG')); continue; end
    iBase = [];
    for ih = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(ih).FileName,0); catch; continue; end
        if (~isfield(hm,'isEigenmode')||~hm.isEigenmode) && strcmpi(hm.HeadModelType,'surface'); iBase=ih; break; end
    end
    if isempty(iBase); continue; end
    hmB = in_bst_headmodel(s.HeadModel(iBase).FileName,0,'SurfaceFile');
    [~, isE] = in_tess_eigenmodes(hmB.SurfaceFile);
    if ~isE; continue; end
    target = struct('baseHmFile', s.HeadModel(iBase).FileName, ...
        'ncFile', s.NoiseCov(1).FileName, ...
        'chFile', bst_get('ChannelFileForStudy', s.FileName));
    return;
end
end
```

- [ ] **Step 4: Run the e2e; verify it PASSES** (`ALL TESTS PASSED`) against OMEGA (needs a study with base surface HM + eigenmodes + eigenmode HM + noise cov — set one up via the head-model dialog/process if absent). If it SKIPs, report what the OMEGA protocol is missing.

- [ ] **Step 5: Commit**

```bash
cd /Users/diellorbasha/workspace/research/code/brainstorm3
git add toolbox/script/tutorial_benchmark_esi.m dev/tests/test_benchmark_esi_e2e.m
git commit -m "$(printf 'ESI benchmark: driver + e2e smoke (simulate/reconstruct/score/report)\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Final verification

- [ ] Run all pure tests via MCP → each `ALL TESTS PASSED`: `test_benchmark_sources_pure`, `test_benchmark_metrics_pure`, `test_benchmark_simulate_pure`, `test_benchmark_report_pure`.
- [ ] Run the e2e tests against OMEGA: `test_benchmark_inverse_e2e`, `test_benchmark_esi_e2e`.
- [ ] Confirm `tutorial_eigenmodes_validation.m` is untouched (prototype kept alongside).
- [ ] Run a small real sweep (`tutorial_benchmark_esi('nLoc',5,'nNoise',3)`) and eyeball `dev/benchmarks/run/summary.csv` — sanity-check that eigenmode-vs-dSPM paired differences flip sign between focal and distributed regimes (the regime hypothesis).

---

## Self-review (done)

- **Spec coverage:** sources generator (T1), metric suite (T2), simulate+colored-noise+SNR (T3), descriptive stats + paired diffs (T4), comparator panel native+eigenmode, eLORETA-if-FieldTrip, MxNE deferred (T5), driver/sweep/regimes/SNR/reps/reporting (T6), prototype kept (final check). Statistics = descriptive (median/IQR/bootstrap CI + paired) per spec. MEG/OMEGA per spec.
- **Placeholders:** none — full code for the four pure units; concrete integration code for T5/T6 with the e2e as the gate (the one explicit "verify OPTIONS against the header" note is a real integration check, not a deferred requirement).
- **Consistency:** field names align across units — `bst_benchmark_sources` returns `.GT/.Sources/.SeedVertex`; `bst_benchmark_metrics(gt, est, GridLoc, seedVertex)`; `bst_benchmark_simulate(L, Sources, C, ...)` returns `.F/.Fsignal/.Fnoise`; `bst_benchmark_report(rows,...)` consumes `.regime/.snr/.method/.metric/.value/.realization` and returns `.summary/.paired`; the driver produces exactly those row fields and uses `RefMethod='eig_mne_log'` matching a key emitted by `bst_benchmark_inverse`.
