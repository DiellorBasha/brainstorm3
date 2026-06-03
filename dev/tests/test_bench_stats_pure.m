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
