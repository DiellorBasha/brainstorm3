function test_sa_figures
% sa_figures writes the expected PNG files from precomputed S (Component 1),
% P (crossing profiles), and optional C (Component 2). Uses tiny synthetic
% structs so it runs headless and fast (no DB needed).
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));

outDir = fullfile(tempdir, 'sa_fig_test');
if exist(outDir, 'dir'), rmdir(outDir, 's'); end

% Minimal synthetic Component-1 stats.
S = struct('dn_sulci_median',0.8,'dn_crown_median',0.3, ...
           'df_sulci_median',0.12,'df_crown_median',0.10,'singEnergyFrac',0.6, ...
           'dn',rand(100,1),'df',rand(100,1),'defined',true(100,1), ...
           'edgeSulcal',[true(50,1);false(50,1)]);
% Minimal synthetic crossing profiles.
al = linspace(0,0.01,12)';
P = struct('path',{(1:12)'}, 'arclen',{al}, ...
           'nAngle',{linspace(0,pi,12)'}, 'fAngle',{linspace(0,0.4,12)'}, ...
           's',{[ones(6,1); -ones(6,1)]});
% Minimal synthetic Component-2 stats.
C = struct('nPairs',40,'signFlipRate',0.7,'tPeak',0.1, ...
           'df',rand(40,1)*0.5,'dt',rand(40,1),'dx',rand(40,1)*2, ...
           'phaseDisc',struct('fiedler_median',0.2,'tang_median',0.6,'xyz_median',1.4));

files = sa_figures(S, P, C, outDir);
assert(~isempty(files), 'sa_figures returned no files.');
for k = 1:numel(files)
    assert(exist(files{k}, 'file') == 2, 'Expected figure file missing: %s', files{k});
end
rmdir(outDir, 's');
fprintf('PASSED: test_sa_figures (%d files).\n', numel(files));
end
