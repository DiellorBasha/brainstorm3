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
% Base MATLAB only (no Statistics Toolbox): percentiles via bst_prctile.
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
            'median',median(v),'iqr',local_iqr(v),'ci_lo',ci(1),'ci_hi',ci(2),'n',numel(v)); %#ok<AGROW>
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
ci = bst_prctile(bm, [2.5 97.5]);
end

function q = local_iqr(v)
q = bst_prctile(v, 75) - bst_prctile(v, 25);
end
