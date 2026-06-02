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
%   SNR        scalar             SNR in dB (converted internally to a linear amplitude ratio for the inverse regularizers)
%
% OUTPUT struct Est: one [nVert x nTime] vertex estimate per method field:
%   .wmne .dspm .sloreta .eig_mne_log .eig_dspm_log
%
% Authors: Diellor Basha, 2026
Est = struct();
iGood = find(goodMask(:));
snrLin = 10^(SNR/20);   % dB power-SNR -> linear amplitude-SNR for the inverse regularizers

% ----- Standard linear inverses via bst_inverse_linear_2018 -----
HeadModel = in_bst_headmodel(baseHmFile, 0);     % unconstrained [nCh x 3*nVert]; 'fixed' orientation applied inside bst_inverse_linear_2018
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
    OPTIONS.SnrFixed       = snrLin;
    HM = HeadModel; HM.Gain = HeadModel.Gain(iGood,:);
    try
        [Results, ~] = bst_inverse_linear_2018(HM, OPTIONS);
        Est.(measures{k,2}) = Results.ImagingKernel * F;
    catch ME
        warning('bst_benchmark_inverse:%s failed: %s', measures{k,2}, ME.message);
    end
end

% ----- Eigenmode variants via the eigenmode head model in this study -----
sStudy = bst_get('AnyFile', baseHmFile);
eigHmFile = '';
for ih = 1:numel(sStudy.HeadModel)
    try hm = in_bst_headmodel(sStudy.HeadModel(ih).FileName,0); catch; continue; end
    if isfield(hm,'isEigenmode') && hm.isEigenmode; eigHmFile = sStudy.HeadModel(ih).FileName; break; end
end
if isempty(eigHmFile)
    warning('bst_benchmark_inverse: no eigenmode head model found in study; skipping eig_mne_log, eig_dspm_log.');
end
if ~isempty(eigHmFile)
    eigCfg = {'mne','eig_mne_log'; 'dspm','eig_dspm_log'};
    for k = 1:size(eigCfg,1)
        [Inv, errE] = bst_inverse_eigenmodes(eigHmFile, ncFile, chFile, goodMask, ...
            'Method', eigCfg{k,1}, 'Prior', 'log', 'SNR', snrLin);
        if ~isempty(errE); continue; end
        [EigK,~] = in_tess_eigenmodes(Inv.SurfaceFile);
        Phi = double(EigK.Vectors(:, 1:Inv.nModes));
        Est.(eigCfg{k,2}) = (Phi * Inv.ImagingKernel) * F;
    end
end

% eLORETA is available only via FieldTrip (process_ft_sourceanalysis); not wired in this v1 panel.
end
