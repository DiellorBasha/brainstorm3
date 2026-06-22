function tutorial_eigenmodes_validation(ReportFile)
% TUTORIAL_EIGENMODES_VALIDATION: Validation harness for eigenmode source mapping.
%
% USAGE:  tutorial_eigenmodes_validation()
%         tutorial_eigenmodes_validation(ReportFile)
%
% DESCRIPTION:
%     Benchmarks the eigenmode (geometric basis function, GBF) source-mapping
%     method against Brainstorm's standard linear inverses (wMNE, dSPM,
%     sLORETA) on the currently loaded protocol (designed for the OMEGA
%     resting tutorial). Runs three validation levels and writes a markdown
%     report. Any level whose data is unavailable is SKIPPED with a logged
%     reason; the script never throws on missing data.
%
%     Level 1 (REQUIRED) - Resolution metrics
%       On one subject with a base surface head model + eigenmodes + noise
%       cov, builds the eigenmode MNE/log (GBF MAP) vertex kernel and the standard wMNE /
%       dSPM / sLORETA vertex kernels on the SAME base head model + noise cov
%       + good MEG channels, then scores each with bst_resolution_metrics
%       against the constrained base leadfield (apples-to-apples). Reports
%       median localization error, median spatial dispersion, and a depth-bias
%       slope (OverallAmplitude vs source depth).
%
%     Level 2 (BEST EFFORT) - Ground-truth simulation
%       Plants a known focal source (single cortical vertex) on the subject's
%       cortex, forward-projects through the base leadfield, adds sensor noise
%       at a target SNR, reconstructs with eigenmode MNE/log and standard dSPM,
%       and reports the distance from the reconstructed peak vertex to the
%       seed. Swept across 2 SNRs.
%
%     Level 3 (REQUIRED OMEGA part) - GBF vs dSPM on real data
%       On each available subject (up to 2) with imported recordings, computes
%       the eigenmode MNE/log and standard-dSPM |source| maps for the same data
%       window and reports their spatial correlation. Phantom localization is
%       attempted only if a phantom protocol is loaded (else SKIPPED).
%
% INPUT:
%     ReportFile : (optional) output markdown path. Default:
%                  <repo>/dev/tests/eigenmode-validation-results.md
%
% Authors: Diellor Basha, 2026

% ===== SETUP =====
thisDir  = fileparts(mfilename('fullpath'));
repoRoot = fileparts(thisDir);                       % toolbox/.. = repo root? No: toolbox/script
repoRoot = fileparts(repoRoot);                      % repo root
if nargin < 1 || isempty(ReportFile)
    ReportFile = fullfile(repoRoot, 'dev', 'tests', 'eigenmode-validation-results.md');
end
if ~brainstorm('status'); brainstorm nogui; end

L = {};                                              % report lines (cellstr)
addline = @(s) sap(s);                               %#ok<NASGU>  (use local below)

L{end+1} = '# Eigenmode source-mapping validation';
L{end+1} = '';
pInfo = bst_get('ProtocolInfo');
L{end+1} = sprintf('- Protocol: `%s`', pInfo.Comment);
L{end+1} = sprintf('- Date: %s', datestr(now, 'yyyy-mm-dd HH:MM'));
L{end+1} = sprintf('- Method under test: eigenmode (GBF) MNE, prior=log (GBF MAP estimate), full mode count K = nModes > nChannels');
L{end+1} = '';

summary = struct('L1', 'SKIP', 'L2', 'SKIP', 'L3omega', 'SKIP', 'L3phantom', 'SKIP');

% ===== LEVEL 1 =====
L{end+1} = '## Level 1 - Resolution metrics (REQUIRED)';
L{end+1} = '';
[l1lines, l1ok] = level1_resolution();
L = [L, l1lines];
if l1ok; summary.L1 = 'OK'; end
L{end+1} = '';

% ===== LEVEL 2 =====
L{end+1} = '## Level 2 - Ground-truth simulation (BEST EFFORT)';
L{end+1} = '';
[l2lines, l2ok] = level2_simulation();
L = [L, l2lines];
if l2ok; summary.L2 = 'OK'; end
L{end+1} = '';

% ===== LEVEL 3 =====
L{end+1} = '## Level 3 - OMEGA GBF-vs-dSPM (REQUIRED) + phantom (BEST EFFORT)';
L{end+1} = '';
[l3lines, l3omegaOk, l3phantomOk] = level3_omega_phantom();
L = [L, l3lines];
if l3omegaOk;   summary.L3omega   = 'OK'; end
if l3phantomOk; summary.L3phantom = 'OK'; end
L{end+1} = '';

% ===== WRITE REPORT =====
[outDir,~,~] = fileparts(ReportFile);
if ~exist(outDir, 'dir'); mkdir(outDir); end
fid = fopen(ReportFile, 'w');
if fid < 0; error('Cannot open report file: %s', ReportFile); end
fprintf(fid, '%s\n', L{:});
fclose(fid);

fprintf('\n==== eigenmode validation summary ====\n');
fprintf('L1(resolution)=%s | L2(simulation)=%s | L3(OMEGA GBF-vs-dSPM)=%s | L3(phantom)=%s\n', ...
    summary.L1, summary.L2, summary.L3omega, summary.L3phantom);
fprintf('Report: %s\n', ReportFile);
end

% =========================================================================
% Trivial no-op (kept so the inline addline handle resolves harmlessly).
function sap(~)
end

% =========================================================================
% LEVEL 1 helper
% =========================================================================
function [lines, ok] = level1_resolution()
lines = {}; ok = false;

% Find a study with a base (non-eigenmode) surface HM + eigenmodes + noise cov
[iStudy, iBaseHM] = find_base_study_with_eigenmodes();
if isempty(iStudy)
    lines{end+1} = 'SKIPPED: no study with a base surface head model + surface eigenmodes + noise cov found.';
    return;
end
sStudy = bst_get('Study', iStudy);
baseHMFile = sStudy.HeadModel(iBaseHM).FileName;
HM = in_bst_headmodel(baseHMFile, 0);
lines{end+1} = sprintf('- Study %d (`%s`), base HM `%s`.', iStudy, sStudy.Name, baseHMFile);

% Eigenmodes for this surface
[E, eok] = in_tess_eigenmodes(HM.SurfaceFile);
if ~eok
    lines{end+1} = 'SKIPPED: surface has no precomputed eigenmodes.';
    return;
end

% Noise cov + channel file + good MEG channels
NoiseCovFile = sStudy.NoiseCov(1).FileName;
ChannelFile  = sStudy.Channel.FileName;
ChMat = in_bst_channel(ChannelFile);
% Good channels: prefer a data ChannelFlag if present, else all-good
ChannelFlag = ones(numel(ChMat.Channel), 1);
if ~isempty(sStudy.Data)
    try
        D = in_bst_data(sStudy.Data(1).FileName, 'ChannelFlag');
        if ~isempty(D.ChannelFlag); ChannelFlag = D.ChannelFlag; end
    catch; end
end
iMeg = good_channel(ChMat.Channel, ChannelFlag, 'MEG');
if numel(iMeg) < 10
    lines{end+1} = 'SKIPPED: fewer than 10 good MEG channels.';
    return;
end
GoodMask = false(numel(ChMat.Channel), 1); GoodMask(iMeg) = true;
lines{end+1} = sprintf('- %d good MEG channels; %d eigenmodes available.', numel(iMeg), E.nModes);

% Constrained base leadfield Lc [nCh x nSrc] restricted to good MEG channels
LcFull = bst_gain_orient(double(HM.Gain), HM.GridOrient, getf(HM,'GridAtlas',[]));
Lc      = LcFull(iMeg, :);                            % [nMeg x nSrc]
GridLoc = HM.GridLoc;                                 % [nSrc x 3] meters
nSrc    = size(GridLoc, 1);

% Use the FULL available mode count. The eigenmode inverse + spectral prior are
% designed for the underdetermined K > nChannels regime (the prior regularizes
% the mode-space system); capping at nChannels defeats the spectral-prior regime.
K = E.nModes;

% --- Eigenmode-dSPM vertex kernel ---
CompHM = bst_eigenmode_leadfield(HM, E, 'nModes', K);
% bst_inverse_eigenmodes needs the composed HM as a file or struct? It loads a file.
% Save a temporary composed head model node, run, then remove it.
[InvE, errE, tmpCompFile, iStudyTmp] = run_eigenmode_inverse_tmp(CompHM, NoiseCovFile, ChannelFile, GoodMask, K, iStudy);
if isempty(InvE)
    lines{end+1} = sprintf('SKIPPED: eigenmode inverse failed (%s).', errE);
    return;
end
Phi = double(E.Vectors(:, 1:InvE.nModes));
KernE = Phi * InvE.ImagingKernel;                    % [nSrc x nMeg]

% --- Standard wMNE / dSPM / sLORETA vertex kernels ---
NC = load(file_fullpath(NoiseCovFile));
% Restrict the base HM to good MEG channels for an apples-to-apples comparison
HMstd = HM;
HMstd.Gain = double(HM.Gain(iMeg, :));
ChTypes = repmat({'MEG'}, 1, numel(iMeg));
NoiseCovMat = NC;
NoiseCovMat.NoiseCov = NC.NoiseCov(iMeg, iMeg);
if isfield(NC,'FourthMoment') && ~isempty(NC.FourthMoment)
    NoiseCovMat.FourthMoment = NC.FourthMoment(iMeg, iMeg);
end

measures = {'amplitude','dspm2018','sloreta'};
labels   = {'wMNE','dSPM','sLORETA'};
stdKern  = cell(1, numel(measures));
for m = 1:numel(measures)
    OPT = bst_inverse_linear_2018();
    OPT.InverseMethod  = 'minnorm';
    OPT.InverseMeasure = measures{m};
    OPT.DataTypes      = {'MEG'};
    OPT.ChannelTypes   = ChTypes;
    OPT.NoiseCovMat    = NoiseCovMat;
    OPT.SourceOrient   = {'fixed'};
    try
        R = bst_inverse_linear_2018(HMstd, OPT);
        stdKern{m} = R.ImagingKernel;                % [nSrc x nMeg] (fixed orient)
    catch ME
        stdKern{m} = [];
        lines{end+1} = sprintf('  - WARN: standard %s failed: %s', labels{m}, ME.message);
    end
end

% --- Depth proxy: distance from inner-skull / head centroid ---
depth = sqrt(sum((GridLoc - mean(GridLoc,1)).^2, 2)) * 1000;   % mm from cortex centroid

% --- Score all kernels ---
lines{end+1} = '';
lines{end+1} = '| Method | median LocError (mm) | median SpatialDisp (mm) | depth-bias slope (amp/mm) |';
lines{end+1} = '|--------|----------------------|--------------------------|----------------------------|';

allKern  = [{KernE}, stdKern];
allLabel = [{'eigenmode-MNE/log'}, labels];
for k = 1:numel(allKern)
    Kn = allKern{k};
    if isempty(Kn) || size(Kn,1) ~= nSrc
        lines{end+1} = sprintf('| %s | - | - | SKIPPED (kernel unavailable / shape %s) |', ...
            allLabel{k}, mat2str(size(Kn)));
        continue;
    end
    Mx = bst_resolution_metrics(Kn, Lc, GridLoc);
    medLoc  = median(Mx.LocError);
    medDisp = median(Mx.SpatialDispersion);
    % depth-bias slope = LSQ slope of (normalized) OverallAmplitude vs depth
    amp = Mx.OverallAmplitude;
    if max(amp) > 0; amp = amp / max(amp); end
    p = polyfit(depth, amp, 1);
    lines{end+1} = sprintf('| %s | %.2f | %.2f | %+.4e |', allLabel{k}, medLoc, medDisp, p(1));
    if k == 1; ok = true; end                        % eigenmode row produced => Level 1 ran
end
lines{end+1} = '';
lines{end+1} = '_Depth proxy = distance of each source from the cortex centroid (mm); slope = LSQ fit of peak-PSF amplitude (normalized) vs that distance. A larger-magnitude slope indicates stronger dependence of resolution amplitude on this depth proxy (i.e. more depth bias); slopes near 0 indicate depth-uniform resolution. (Cortex-centroid distance is a crude proxy; treat magnitudes comparatively across methods, not absolutely.)_';

% Cleanup temporary composed head model node
cleanup_tmp_headmodel(tmpCompFile, iStudyTmp);
end

% =========================================================================
% LEVEL 2 helper - self-contained ground-truth simulation through leadfield
% =========================================================================
function [lines, ok] = level2_simulation()
lines = {}; ok = false;

[iStudy, iBaseHM] = find_base_study_with_eigenmodes();
if isempty(iStudy)
    lines{end+1} = 'SKIPPED: no study with base HM + eigenmodes + noise cov found.';
    return;
end
sStudy = bst_get('Study', iStudy);
HM = in_bst_headmodel(sStudy.HeadModel(iBaseHM).FileName, 0);
[E, eok] = in_tess_eigenmodes(HM.SurfaceFile);
if ~eok; lines{end+1} = 'SKIPPED: no eigenmodes.'; return; end

ChannelFile = sStudy.Channel.FileName;
NoiseCovFile = sStudy.NoiseCov(1).FileName;
ChMat = in_bst_channel(ChannelFile);
iMeg = good_channel(ChMat.Channel, ones(numel(ChMat.Channel),1), 'MEG');
GoodMask = false(numel(ChMat.Channel),1); GoodMask(iMeg) = true;

LcFull = bst_gain_orient(double(HM.Gain), HM.GridOrient, getf(HM,'GridAtlas',[]));
Lc = LcFull(iMeg, :);                                 % [nMeg x nSrc]
GridLoc = HM.GridLoc;
nSrc = size(GridLoc,1);
K = E.nModes;   % full mode count (K > nChannels regime), see Level 1 note

% Build eigenmode (MNE/log) and standard-dSPM vertex kernels once (reused per SNR)
CompHM = bst_eigenmode_leadfield(HM, E, 'nModes', K);
[InvE, errE, tmpCompFile, iStudyTmp] = run_eigenmode_inverse_tmp(CompHM, NoiseCovFile, ChannelFile, GoodMask, K, iStudy);
if isempty(InvE)
    lines{end+1} = sprintf('SKIPPED: eigenmode inverse failed (%s).', errE);
    return;
end
Phi = double(E.Vectors(:, 1:InvE.nModes));
KernE = Phi * InvE.ImagingKernel;                    % [nSrc x nMeg]

NC = load(file_fullpath(NoiseCovFile));
HMstd = HM; HMstd.Gain = double(HM.Gain(iMeg,:));
NoiseCovMat = NC; NoiseCovMat.NoiseCov = NC.NoiseCov(iMeg, iMeg);
if isfield(NC,'FourthMoment') && ~isempty(NC.FourthMoment)
    NoiseCovMat.FourthMoment = NC.FourthMoment(iMeg, iMeg);
end
OPT = bst_inverse_linear_2018();
OPT.InverseMethod='minnorm'; OPT.InverseMeasure='dspm2018';
OPT.DataTypes={'MEG'}; OPT.ChannelTypes=repmat({'MEG'},1,numel(iMeg));
OPT.NoiseCovMat=NoiseCovMat; OPT.SourceOrient={'fixed'};
try
    Rstd = bst_inverse_linear_2018(HMstd, OPT);
    KernS = Rstd.ImagingKernel;
catch ME
    cleanup_tmp_headmodel(tmpCompFile, iStudyTmp);
    lines{end+1} = sprintf('SKIPPED: standard dSPM build failed (%s).', ME.message);
    return;
end

% Deterministic seed vertices spread across the cortex (averaged for stability,
% since a single seed + single noise draw is sensitive to the random topography).
rng(42);
nSeed = 12;
seedIdx = unique(round(linspace(round(0.04*nSrc), round(0.96*nSrc), nSeed)));

% Noise covariance for adding realistic colored sensor noise
Cn = NC.NoiseCov(iMeg, iMeg); Cn = 0.5*(Cn+Cn');
[Un,Sn] = svd(Cn); sn = max(diag(Sn), 0);
Wnoise = Un * diag(sqrt(sn));                         % colored-noise generator

lines{end+1} = sprintf('- Self-contained simulation: %d focal seeds spread across the cortex; localization error averaged (median) over seeds.', numel(seedIdx));
lines{end+1} = '';
lines{end+1} = '| SNR (amp) | eigenmode-MNE/log median LocError (mm) | standard dSPM median LocError (mm) |';
lines{end+1} = '|-----------|--------------------------------------|--------------------------------------|';

SNRs = [3, 10];
for snr = SNRs
    eErr = nan(numel(seedIdx),1);
    sErr = nan(numel(seedIdx),1);
    for q = 1:numel(seedIdx)
        iSeed = seedIdx(q);
        seedLoc = GridLoc(iSeed, :);
        b0 = Lc(:, iSeed);
        if norm(b0) == 0; continue; end
        noise = Wnoise * randn(numel(iMeg), 1);
        if norm(noise) == 0; continue; end
        % scale signal so amplitude SNR = ||signal||/||noise|| = snr
        a = snr * norm(noise) / norm(b0);
        b = a * b0 + noise;                          % [nMeg x 1] sensor measurement
        [~, iE] = max(abs(KernE * b));
        [~, iS] = max(abs(KernS * b));
        eErr(q) = norm(GridLoc(iE,:) - seedLoc) * 1000;
        sErr(q) = norm(GridLoc(iS,:) - seedLoc) * 1000;
    end
    if all(isnan(eErr))
        lines{end+1} = sprintf('| %d | SKIPPED (degenerate signal/noise) | - |', snr); %#ok<AGROW>
        continue;
    end
    lines{end+1} = sprintf('| %d | %.2f | %.2f |', snr, ...
        median(eErr,'omitnan'), median(sErr,'omitnan')); %#ok<AGROW>
    ok = true;
end
lines{end+1} = '';
lines{end+1} = '_Ground truth = leadfield forward of a single cortical vertex + colored sensor noise; error = distance from reconstructed peak to seed, median over seeds. (The GUI simulate processes were not driven headlessly; this self-contained forward simulation provides the same ground-truth comparison.)_';

cleanup_tmp_headmodel(tmpCompFile, iStudyTmp);
end

% =========================================================================
% LEVEL 3 helper
% =========================================================================
function [lines, omegaOk, phantomOk] = level3_omega_phantom()
lines = {}; omegaOk = false; phantomOk = false;

% --- OMEGA: GBF-dSPM vs standard-dSPM spatial correlation on real data ---
lines{end+1} = '### OMEGA: spatial correlation of |source| maps (eigenmode-MNE/log vs standard dSPM)';
lines{end+1} = '';
studies = find_studies_with_imported_data(2);
if isempty(studies)
    lines{end+1} = 'SKIPPED: no subject study with imported recordings + base HM + eigenmodes found.';
else
    lines{end+1} = '| Subject study | nVert | spatial corr (peak time) | spatial corr (time-avg) |';
    lines{end+1} = '|---------------|-------|--------------------------|--------------------------|';
    for ii = 1:numel(studies)
        st = studies(ii);
        [row, okSub] = level3_one_subject(st.iStudy, st.iBaseHM, st.iData);
        lines{end+1} = row; %#ok<AGROW>
        if okSub; omegaOk = true; end
    end
end
lines{end+1} = '';

% --- Phantom (best effort) ---
lines{end+1} = '### Phantom localization';
lines{end+1} = '';
if has_phantom_protocol()
    lines{end+1} = 'NOTE: a phantom-named protocol was detected but per-protocol switching is not driven by this harness; phantom localization not computed in this run.';
else
    lines{end+1} = 'SKIPPED: no phantom protocol loaded (no protocol with name containing "phantom").';
end
end

% --- one subject for Level 3 OMEGA ---
function [row, ok] = level3_one_subject(iStudy, iBaseHM, iData)
ok = false; row = '';
sStudy = bst_get('Study', iStudy);
HM = in_bst_headmodel(sStudy.HeadModel(iBaseHM).FileName, 0);
[E, eok] = in_tess_eigenmodes(HM.SurfaceFile);
if ~eok; row = sprintf('| %s | - | SKIPPED (no eigenmodes) | - |', sStudy.Name); return; end

ChannelFile = sStudy.Channel.FileName;
NoiseCovFile = sStudy.NoiseCov(1).FileName;
ChMat = in_bst_channel(ChannelFile);
D = in_bst_data(sStudy.Data(iData).FileName);
iMeg = good_channel(ChMat.Channel, D.ChannelFlag, 'MEG');
GoodMask = false(numel(ChMat.Channel),1); GoodMask(iMeg) = true;

% Data window: up to 200 time points to bound cost
F = double(D.F(iMeg, :));
nT = size(F,2);
maxT = min(nT, 200);
ti = round(linspace(1, nT, maxT));
F = F(:, ti);

K = E.nModes;   % full mode count (K > nChannels regime), see Level 1 note

% eigenmode (MNE/log) vertex map
CompHM = bst_eigenmode_leadfield(HM, E, 'nModes', K);
[InvE, errE, tmpCompFile, iStudyTmp] = run_eigenmode_inverse_tmp(CompHM, NoiseCovFile, ChannelFile, GoodMask, K, iStudy);
if isempty(InvE)
    row = sprintf('| %s | - | SKIPPED (eig inverse: %s) | - |', sStudy.Name, errE); return;
end
Phi = double(E.Vectors(:, 1:InvE.nModes));
KernE = Phi * InvE.ImagingKernel;                    % [nSrc x nMeg]
JE = KernE * F;                                      % [nSrc x nT]

% standard dSPM vertex map
NC = load(file_fullpath(NoiseCovFile));
HMstd = HM; HMstd.Gain = double(HM.Gain(iMeg,:));
NoiseCovMat = NC; NoiseCovMat.NoiseCov = NC.NoiseCov(iMeg,iMeg);
if isfield(NC,'FourthMoment') && ~isempty(NC.FourthMoment)
    NoiseCovMat.FourthMoment = NC.FourthMoment(iMeg,iMeg);
end
OPT = bst_inverse_linear_2018();
OPT.InverseMethod='minnorm'; OPT.InverseMeasure='dspm2018';
OPT.DataTypes={'MEG'}; OPT.ChannelTypes=repmat({'MEG'},1,numel(iMeg));
OPT.NoiseCovMat=NoiseCovMat; OPT.SourceOrient={'fixed'};
try
    Rstd = bst_inverse_linear_2018(HMstd, OPT);
    KernS = Rstd.ImagingKernel;
catch ME
    cleanup_tmp_headmodel(tmpCompFile, iStudyTmp);
    row = sprintf('| %s | - | SKIPPED (std dSPM: %s) | - |', sStudy.Name, ME.message); return;
end
JS = KernS * F;

% |source| maps
AE = abs(JE); AS = abs(JS);
nVert = size(AE,1);

% peak-time spatial correlation (time of max global field in standard dSPM)
gfp = sqrt(sum(AS.^2, 1));
[~, iPk] = max(gfp);
cPeak = corr_safe(AE(:,iPk), AS(:,iPk));
% time-averaged map correlation
cAvg = corr_safe(mean(AE,2), mean(AS,2));

row = sprintf('| %s | %d | %.3f | %.3f |', sStudy.Name, nVert, cPeak, cAvg);
ok = true;
cleanup_tmp_headmodel(tmpCompFile, iStudyTmp);
end

% =========================================================================
% Shared helpers
% =========================================================================
function [iStudy, iBaseHM] = find_base_study_with_eigenmodes()
% First study that has a base (non-eigenmode) surface HM, surface eigenmodes,
% and a noise covariance.
iStudy = []; iBaseHM = [];
sStudies = bst_get('ProtocolStudies');
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.HeadModel); continue; end
    if ~isfield(s,'NoiseCov') || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    for k = 1:numel(s.HeadModel)
        try
            hm = in_bst_headmodel(s.HeadModel(k).FileName, 0);
        catch; continue; end
        isEig = isfield(hm,'isEigenmode') && hm.isEigenmode==1;
        if isEig; continue; end
        if ~strcmpi(hm.HeadModelType,'surface'); continue; end
        if isempty(hm.GridLoc) || isempty(hm.GridOrient); continue; end
        [~, eok] = in_tess_eigenmodes(hm.SurfaceFile);
        if ~eok; continue; end
        iStudy = iS; iBaseHM = k; return;
    end
end
end

function studies = find_studies_with_imported_data(maxN)
% Studies with imported (non-raw) recordings + base surface HM + eigenmodes + noise cov.
studies = struct('iStudy',{},'iBaseHM',{},'iData',{});
sStudies = bst_get('ProtocolStudies');
for iS = 1:numel(sStudies.Study)
    s = sStudies.Study(iS);
    if isempty(s.HeadModel) || isempty(s.Data); continue; end
    if ~isfield(s,'NoiseCov') || isempty(s.NoiseCov) || isempty(s.NoiseCov(1).FileName); continue; end
    % base HM with eigenmodes
    iBaseHM = [];
    for k = 1:numel(s.HeadModel)
        try hm = in_bst_headmodel(s.HeadModel(k).FileName, 0); catch; continue; end
        if isfield(hm,'isEigenmode') && hm.isEigenmode==1; continue; end
        if ~strcmpi(hm.HeadModelType,'surface'); continue; end
        if isempty(hm.GridLoc); continue; end
        [~, eok] = in_tess_eigenmodes(hm.SurfaceFile);
        if eok; iBaseHM = k; break; end
    end
    if isempty(iBaseHM); continue; end
    % imported (recordings) data
    iData = [];
    for d = 1:numel(s.Data)
        if strcmpi(s.Data(d).DataType, 'recordings'); iData = d; break; end
    end
    if isempty(iData); continue; end
    studies(end+1) = struct('iStudy',iS,'iBaseHM',iBaseHM,'iData',iData); %#ok<AGROW>
    if numel(studies) >= maxN; break; end
end
end

function [Inv, errMsg, tmpFile, iStudyTmp] = run_eigenmode_inverse_tmp(CompHM, NoiseCovFile, ChannelFile, GoodMask, K, iStudy)
% Save the composed head model to a temporary file in the study folder and run
% the eigenmode inverse (file-based: bst_inverse_eigenmodes resolves the path
% via in_bst_headmodel, no DB registration required). Returns the temp file so
% the caller can delete it.
Inv = []; errMsg = ''; tmpFile = ''; iStudyTmp = iStudy;
try
    CompHM.Comment = sprintf('TMP eigenmode validation (%d modes)', K);
    sStudy = bst_get('Study', iStudy);
    studyPath = bst_fileparts(file_fullpath(sStudy.FileName));
    tmpFull = bst_process('GetNewFilename', studyPath, 'headmodel_eigenmode_tmpval');
    bst_save(tmpFull, CompHM, 'v7');
    tmpFile = file_short(tmpFull);
    % Headline eigenmode method = MNE with the log spectral prior (the GBF MAP
    % estimate). This is the configuration designed for the K > nChannels regime.
    [Inv, errMsg] = bst_inverse_eigenmodes(tmpFile, NoiseCovFile, ChannelFile, GoodMask, ...
        'Method','mne', 'Prior','log', 'SNR',3, 'nModes',K);
catch ME
    errMsg = ME.message;
end
end

function cleanup_tmp_headmodel(tmpFile, iStudy) %#ok<INUSD>
if isempty(tmpFile); return; end
try
    f = file_fullpath(tmpFile);
    if exist(f,'file'); file_delete(f, 1); end
catch; end
end

function tf = has_phantom_protocol()
tf = false;
try
    global GlobalData %#ok<GVMIS>
    if isfield(GlobalData,'DataBase') && isfield(GlobalData.DataBase,'ProtocolInfo')
        PI = GlobalData.DataBase.ProtocolInfo;
        for i = 1:numel(PI)
            if ~isempty(strfind(lower(PI(i).Comment), 'phantom')); tf = true; return; end %#ok<STREMP>
        end
    end
catch; end
end

function c = corr_safe(a, b)
a = a(:); b = b(:);
a = a - mean(a); b = b - mean(b);
da = norm(a); db = norm(b);
if da == 0 || db == 0; c = 0; else; c = (a'*b)/(da*db); end
end

function v = getf(s, f, d)
if isfield(s,f) && ~isempty(s.(f)); v = s.(f); else; v = d; end
end
