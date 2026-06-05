function test_sa_smoothness
% Component 1 decoupling claim: on real cortex, normal angular variation is
% ELEVATED on sulcal edges, while the Fiedler covariant variation is NOT
% (statistically flat across the sulcal/crown partition). SKIP if no cortex.
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

% Normal variation must be higher on sulcal edges than crown edges.
assert(S.dn_sulci_median > S.dn_crown_median, ...
    'Normal variation should be elevated at sulci (sulci=%.3f, crown=%.3f).', ...
    S.dn_sulci_median, S.dn_crown_median);
% Fiedler covariant variation must NOT be elevated at sulci: the sulcal median
% should not exceed the crown median by more than the normal field does. Use a
% loose ratio test (decoupling): df elevation ratio << dn elevation ratio.
dnRatio = S.dn_sulci_median / max(S.dn_crown_median, eps);
dfRatio = S.df_sulci_median / max(S.df_crown_median, eps);
assert(dfRatio < dnRatio, ...
    'Fiedler variation should be less folding-coupled than the normal (dfRatio=%.2f, dnRatio=%.2f).', ...
    dfRatio, dnRatio);
% Singularity-energy fraction is a valid fraction.
assert(S.singEnergyFrac >= 0 && S.singEnergyFrac <= 1, 'singEnergyFrac out of range.');
fprintf('PASSED: test_sa_smoothness (dnRatio=%.2f, dfRatio=%.2f, singFrac=%.2f).\n', ...
    dnRatio, dfRatio, S.singEnergyFrac);
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
