function test_sa_continuity
% Component 2: end-to-end on the unconstrained kernel. Conditioned pair set is
% non-empty; constrained sign-flip rate exceeds the Fiedler phase-discontinuity
% RATE; Fiedler phase-disc rate < 0.6 (better than random on sulcal walls).
% Median ordering is NOT asserted (cross-wall geometry makes medians all near-random).
% SKIP if no cortex or no unconstrained kernel.
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(fileparts(thisDir));
addpath(repoRoot);
addpath(fullfile(repoRoot, 'dev', 'benchmarks', 'sign_ambiguity'));
if ~brainstorm('status'), brainstorm nogui; end
[isOk, errMsg] = bst_plugin('Install', 'nxr-compute');
assert(isOk, 'nxr-compute required: %s', errMsg);
bst_plugin('Load', 'nxr-compute');
SurfaceFile = find_cortex_20484V();
if isempty(SurfaceFile)
    fprintf('SKIP: no 20484-vertex cortex.\n'); return;
end
ResultsFile = find_unconstrained_kernel(SurfaceFile);
if isempty(ResultsFile)
    fprintf('SKIP: no unconstrained (nComponents=3) kernel on this surface.\n'); return;
end

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end
opts = struct('MaxDist',0.003, 'NormalDot',-0.7, 'Nring',3, ...
              'Window',[0.05 0.15], 'MagRatio',0.5);
C = sa_continuity(F, SulciMap, ResultsFile, opts);

assert(C.nPairs > 0, 'Conditioned pair set must be non-empty.');
assert(C.signFlipRate > C.phaseDiscRate.fiedler, ...
    'Constrained sign-flip rate (%.2f) should exceed Fiedler phase-disc rate (%.2f).', ...
    C.signFlipRate, C.phaseDiscRate.fiedler);
% The phase-discontinuity RATE (fraction of pairs with circular dist > pi/2)
% is what matters: Fiedler avoids the constrained flip artifact.
assert(C.phaseDiscRate.fiedler < C.signFlipRate, ...
    'Fiedler phase-disc rate (%.2f) should be below constrained sign-flip rate (%.2f).', ...
    C.phaseDiscRate.fiedler, C.signFlipRate);
assert(C.phaseDiscRate.fiedler < 0.6, ...
    'Fiedler phase-disc rate (%.2f) should be below 0.6 (better than random on sulcal walls).', ...
    C.phaseDiscRate.fiedler);
fprintf(['PASSED: test_sa_continuity (nPairs=%d, signFlip=%.2f, ', ...
    'fiedlerDiscRate=%.2f, phaseDisc fiedler/tang/xyz = %.2f/%.2f/%.2f).\n'], ...
    C.nPairs, C.signFlipRate, C.phaseDiscRate.fiedler, ...
    C.phaseDisc.fiedler_median, C.phaseDisc.tang_median, C.phaseDisc.xyz_median);
end

function SurfaceFile = find_cortex_20484V()
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
end

function ResultsFile = find_unconstrained_kernel(SurfaceFile)
% Search all studies for a link results file on this surface with nComponents==3.
% Link files are preferred (bare kernel files lack a data file for LoadFull=1).
ResultsFile = '';
sStudies = bst_get('ProtocolStudies');
if isempty(sStudies), return; end
allStudy = [sStudies.Study];
for iS = 1:numel(allStudy)
    res = allStudy(iS).Result;
    for iR = 1:numel(res)
        fname = res(iR).FileName;
        if isempty(fname), continue; end
        % Only consider link files (kernel+data pairs); bare kernels have no data.
        if ~strncmp(fname, 'link|', 5), continue; end
        try
            [kernelFile, ~] = file_resolve_link(fname);
            M = load(kernelFile, 'nComponents');
        catch, continue; end
        if isfield(M,'nComponents') && isequal(M.nComponents, 3)
            ResultsFile = fname; return;
        end
    end
end
end
