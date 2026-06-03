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
