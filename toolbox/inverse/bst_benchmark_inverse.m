function Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)
% BST_BENCHMARK_INVERSE: Reconstruct simulated data with the comparator panel.
%
% USAGE:  Est = bst_benchmark_inverse(F, baseHmFile, ncFile, chFile, goodMask, SNR)
%
% INPUTS:
%   F          [nGoodCh x nTime]  simulated good-channel sensor data
%   baseHmFile char               base (non-eigenmode) surface head model file
%   ncFile     char               noise covariance file
%   chFile     char               channel file
%   goodMask   [nAllCh x 1]       logical good-channel mask
%   SNR        scalar             SNR (dB) -> SnrFixed for the linear inverses
%
% OUTPUT struct Est: one [nVert x nTime] vertex estimate per method field:
%   .wmne .dspm .sloreta .eig_mne_log .eig_dspm_log
%
% Authors: Diellor Basha, 2026
Est = struct();
iGood = find(goodMask(:));

% ----- Standard linear inverses via bst_inverse_linear_2018 -----
HeadModel = in_bst_headmodel(baseHmFile, 0);     % unconstrained gain [nCh x 3*nVert]
NC = load(file_fullpath(ncFile));
NoiseCovMat = struct('NoiseCov', NC.NoiseCov(iGood,iGood), ...
                     'nSamples', [], ...
                     'FourthMoment', []);
ChannelMat  = in_bst_channel(chFile);
chTypes = {ChannelMat.Channel(iGood).Type};

measures = {'amplitude','wmne'; 'dspm2018','dspm'; 'sloreta','sloreta'};
for k = 1:size(measures,1)
    OPTIONS = bst_inverse_linear_2018();
    OPTIONS.InverseMethod  = 'minnorm';
    OPTIONS.InverseMeasure = measures{k,1};
    OPTIONS.NoiseCovMat    = NoiseCovMat;
    OPTIONS.ChannelTypes   = chTypes;
    OPTIONS.SourceOrient   = {'fixed'};
    OPTIONS.SnrMethod      = 'fixed';
    OPTIONS.SnrFixed       = SNR;
    HM = HeadModel; HM.Gain = HeadModel.Gain(iGood,:);
    [Results, ~] = bst_inverse_linear_2018(HM, OPTIONS);
    Est.(measures{k,2}) = Results.ImagingKernel * F;
end

% ----- Eigenmode variants via the eigenmode head model in this study -----
sStudy = bst_get('AnyFile', baseHmFile);
eigHmFile = '';
for ih = 1:numel(sStudy.HeadModel)
    try hm = in_bst_headmodel(sStudy.HeadModel(ih).FileName,0); catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode; eigHmFile = sStudy.HeadModel(ih).FileName; break; end
end
if ~isempty(eigHmFile)
    HMe = in_bst_headmodel(eigHmFile, 0);
    [Eig,~] = in_tess_eigenmodes(HMe.SurfaceFile);
    eigCfg = {'mne','eig_mne_log'; 'dspm','eig_dspm_log'};
    for k = 1:size(eigCfg,1)
        [Inv, errE] = bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, ...
            'Method', eigCfg{k,1}, 'Prior', 'log', 'SNR', SNR);
        if ~isempty(errE); continue; end
        Phi = double(Eig.Vectors(:,1:Inv.nModes));
        Est.(eigCfg{k,2}) = (Phi * Inv.ImagingKernel) * F;
    end
end

% ----- eLORETA via FieldTrip, only if available -----
if exist('ft_sourceanalysis','file') == 2
    % Optional; left unset when FieldTrip is absent.
end
end
