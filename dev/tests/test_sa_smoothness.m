function test_sa_smoothness
% Component 1 decoupling claim: on real cortex, normal angular variation is
% PARTITION-SENSITIVE (crown edges are more curved than sulcal walls), while
% the Fiedler covariant variation is LESS sensitive to the partition. SKIP if no cortex.
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

F = sa_frames(SurfaceFile);
TessMat = in_tess_bst(SurfaceFile);
if isfield(TessMat,'SulciMap') && ~isempty(TessMat.SulciMap)
    SulciMap = TessMat.SulciMap;
else
    SulciMap = tess_sulcimap(TessMat);
end

S = sa_smoothness(F, SulciMap);

% The normal's angular variation differs between partitions (crown edges are
% more curved — convex ridges — so variation is higher there than at sulcal
% walls). Verify this partition sensitivity is real.
assert(abs(S.dn_sulci_median - S.dn_crown_median) > 0.01, ...
    'Normal variation should differ between partitions (sulci=%.3f, crown=%.3f).', ...
    S.dn_sulci_median, S.dn_crown_median);
% Fiedler covariant variation must be LESS sensitive to the sulcal/crown
% partition than the normal field is -- this is the decoupling claim.
dnSens = abs(S.dn_sulci_median - S.dn_crown_median) / max(S.dn_crown_median, eps);
dfSens = abs(S.df_sulci_median - S.df_crown_median) / max(S.df_crown_median, eps);
assert(dfSens < dnSens, ...
    'Fiedler should be less partition-sensitive than normal (dfSens=%.3f, dnSens=%.3f).', ...
    dfSens, dnSens);
% Singularity-energy fraction is a valid fraction.
assert(S.singEnergyFrac >= 0 && S.singEnergyFrac <= 1, 'singEnergyFrac out of range.');
fprintf('PASSED: test_sa_smoothness (dnSens=%.3f, dfSens=%.3f, singFrac=%.2f).\n', ...
    dnSens, dfSens, S.singEnergyFrac);
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
