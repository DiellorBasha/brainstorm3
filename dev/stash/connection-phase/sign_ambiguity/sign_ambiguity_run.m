function out = sign_ambiguity_run(SurfaceFile, ResultsFile)
% SIGN_AMBIGUITY_RUN: Driver for the connection-eigenmode sulcal-wall continuity
% test. Runs Component 1 (data-free smoothness) on a cortex surface, and
% Component 2 (real-data phase continuity) if an unconstrained kernel is found.
% Writes figures + stats under dev/benchmarks/sign_ambiguity/sign_ambiguity_run/.
%
% USAGE:  out = sign_ambiguity_run()                       % auto-resolve
%         out = sign_ambiguity_run(SurfaceFile)            % given surface
%         out = sign_ambiguity_run(SurfaceFile, ResultsFile)
%
% Returns [] if no usable cortex is found (SKIP). Otherwise a struct with
% outDir, statsFile, figures, C1 (Component 1 stats), C2 (Component 2 stats or []).

if nargin < 1 || isempty(SurfaceFile), SurfaceFile = local_find_cortex(); end
if isempty(SurfaceFile), out = []; return; end
if nargin < 2, ResultsFile = local_find_kernel(SurfaceFile); end

thisDir = fileparts(mfilename('fullpath'));
outDir  = fullfile(thisDir, 'sign_ambiguity_run');
if ~exist(outDir, 'dir')
    [ok, msg] = mkdir(outDir);
    if ~ok, error('sign_ambiguity_run:mkdir', 'Cannot create output dir: %s', msg); end
end

% Shared frames + SulciMap.
F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

% Component 1 (always).
C1 = sa_smoothness(F, SulciMap);

% Component 2 (if a kernel is available).
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3, ...
              'Window',[0.05 0.15], 'MagRatio',0.5);
C2 = [];
if ~isempty(ResultsFile)
    C2 = sa_continuity(F, SulciMap, ResultsFile, opts);
end

% Crossing profiles (Component 1 illustration; uses a few facing-wall pairs).
% Overlay the constrained scalar from the peak current when Component 2 ran.
pairs = sa_sulcal_walls(F.Vtx, F.Nv, SulciMap, F.VertConn, opts);
if ~isempty(C2)
    P = sa_crossing_profile(F, pairs, 3, C2.J);
else
    P = sa_crossing_profile(F, pairs, 3);
end

% Figures + stats.
figures = sa_figures(C1, P, C2, outDir);
statsFile = fullfile(outDir, 'stats.md');
local_write_stats(statsFile, SurfaceFile, ResultsFile, C1, C2);
local_write_csv(fullfile(outDir, 'stats.csv'), C1, C2);

out = struct('outDir', outDir, 'statsFile', statsFile, 'figures', {figures}, ...
             'C1', C1, 'C2', C2, 'P', {P}, 'SurfaceFile', SurfaceFile, 'ResultsFile', ResultsFile);
fprintf('sign_ambiguity_run: wrote %s\n', outDir);
end

% ---- local helpers ----
function SurfaceFile = local_find_cortex()
SurfaceFile = '';
sSubjects = bst_get('ProtocolSubjects');
if isempty(sSubjects), return; end
allSubj = [sSubjects.Subject];
for iS = 1:numel(allSubj)
    surf = allSubj(iS).Surface;
    for iF = 1:numel(surf)
        if ~strcmpi(surf(iF).SurfaceType, 'Cortex'), continue; end
        try, T = load(file_fullpath(surf(iF).FileName), 'Vertices'); catch, continue; end
        if size(T.Vertices,1) == 20484, SurfaceFile = surf(iF).FileName; return; end
    end
end
if isempty(SurfaceFile)
    fprintf('sign_ambiguity_run: no 20484-vertex cortex found in current protocol -- skipping.\n');
end
end

function ResultsFile = local_find_kernel(SurfaceFile)
% Find an unconstrained (nComponents=3) link file paired with data.
ResultsFile = '';
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies), return; end
allStudy = [sStudies.Study];
for iS = 1:numel(allStudy)
    res = allStudy(iS).Result;
    for iR = 1:numel(res)
        fname = res(iR).FileName;
        if isempty(fname), continue; end
        % Only link files have an associated data file for LoadFull=1.
        if ~strncmp(fname, 'link|', 5), continue; end
        try, M = load(file_fullpath(file_resolve_link(fname)), 'nComponents'); catch, continue; end
        if isfield(M,'nComponents') && isequal(M.nComponents, 3)
            ResultsFile = fname; return;
        end
    end
end
end

function local_write_stats(statsFile, SurfaceFile, ResultsFile, C1, C2)
fid = fopen(statsFile, 'w');
if fid < 0, error('sign_ambiguity_run:io', 'Cannot write: %s', statsFile); end
c = onCleanup(@() fclose(fid));
fprintf(fid, '# Sign-ambiguity continuity test\n\n');
fprintf(fid, '- Surface: `%s`\n', SurfaceFile);
fprintf(fid, '- Kernel: `%s`\n\n', local_iff(isempty(ResultsFile),'(none)',ResultsFile));
fprintf(fid, '## Component 1 (data-free smoothness)\n\n');
fprintf(fid, '| metric | sulcal | crown |\n|---|---|---|\n');
fprintf(fid, '| normal var dn (median rad) | %.4f | %.4f |\n', C1.dn_sulci_median, C1.dn_crown_median);
fprintf(fid, '| Fiedler var df (median rad) | %.4f | %.4f |\n', C1.df_sulci_median, C1.df_crown_median);
fprintf(fid, '\n- singularity energy fraction of df^2: %.3f\n', C1.singEnergyFrac);
fprintf(fid, '- decoupling: dn sulci/crown ratio = %.2f, df sulci/crown ratio = %.2f\n\n', ...
    C1.dn_sulci_median/max(C1.dn_crown_median,eps), C1.df_sulci_median/max(C1.df_crown_median,eps));
if ~isempty(C2)
    fprintf(fid, '## Component 2 (real-data phase continuity)\n\n');
    fprintf(fid, '- peak time: %.0f ms, conditioned pairs: %d\n', 1000*C2.tPeak, C2.nPairs);
    fprintf(fid, '- constrained sign-flip rate: %.3f\n\n', C2.signFlipRate);
    fprintf(fid, '| frame | median phase disc (rad) | disc rate (>pi/2) |\n|---|---|---|\n');
    fprintf(fid, '| Fiedler | %.3f | %.3f |\n', C2.phaseDisc.fiedler_median, C2.phaseDiscRate.fiedler);
    fprintf(fid, '| tess_tangents | %.3f | %.3f |\n', C2.phaseDisc.tang_median, C2.phaseDiscRate.tang);
    fprintf(fid, '| global-xyz | %.3f | %.3f |\n', C2.phaseDisc.xyz_median, C2.phaseDiscRate.xyz);
end
end

function local_write_csv(csvFile, C1, C2)
fid = fopen(csvFile, 'w');
if fid < 0, error('sign_ambiguity_run:io', 'Cannot write: %s', csvFile); end
c = onCleanup(@() fclose(fid));
fprintf(fid, 'component,metric,value\n');
fprintf(fid, 'C1,dn_sulci_median,%.6f\n', C1.dn_sulci_median);
fprintf(fid, 'C1,dn_crown_median,%.6f\n', C1.dn_crown_median);
fprintf(fid, 'C1,df_sulci_median,%.6f\n', C1.df_sulci_median);
fprintf(fid, 'C1,df_crown_median,%.6f\n', C1.df_crown_median);
fprintf(fid, 'C1,singEnergyFrac,%.6f\n', C1.singEnergyFrac);
if ~isempty(C2)
    fprintf(fid, 'C2,nPairs,%d\n', C2.nPairs);
    fprintf(fid, 'C2,signFlipRate,%.6f\n', C2.signFlipRate);
    fprintf(fid, 'C2,phaseDisc_fiedler,%.6f\n', C2.phaseDisc.fiedler_median);
    fprintf(fid, 'C2,phaseDisc_tang,%.6f\n', C2.phaseDisc.tang_median);
    fprintf(fid, 'C2,phaseDisc_xyz,%.6f\n', C2.phaseDisc.xyz_median);
end
end

function s = local_iff(cond, a, b)
if cond, s = a; else, s = b; end
end
