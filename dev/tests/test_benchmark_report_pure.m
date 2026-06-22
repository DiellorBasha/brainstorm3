function test_benchmark_report_pure
% Verify aggregation: median/IQR/bootstrap CI + paired differences.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot, fullfile(repoRoot,'dev','benchmarks'));
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
